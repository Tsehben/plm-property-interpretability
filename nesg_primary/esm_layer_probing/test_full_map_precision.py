"""Offline checks; synthetic tensors are not NESG precision-pilot results."""

import importlib.util
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import torch

from common import read_extraction, sha
from pilot_full_map_precision import SOURCE_SHA, compare, roundtrip, select_sample


class PrecisionPilotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source = Path(__file__).resolve().parents[2] / "interprot_ref/raw/sae_model.py"
        assert sha(source) == SOURCE_SHA
        spec = importlib.util.spec_from_file_location("reference_sae", source)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        torch.manual_seed(123)
        cls.sae = module.SparseAutoencoder(8, 32, k=4).eval()

    def test_label_blind_sampling_is_stable_and_covers_edges(self):
        data = read_extraction()  # Development extraction manifest only, no labels.
        selected = select_sample(data)
        self.assertLessEqual(len(selected), 11)
        self.assertEqual(selected.sequence.str.contains("X").sum(), 3)
        self.assertEqual(set(selected.split), {"train", "validation"})
        self.assertEqual(
            list(selected.record_id),
            list(select_sample(data.sample(frac=1, random_state=42)).record_id),
        )
        for split in ("train", "validation"):
            full = data[data.split == split]
            part = selected[selected.split == split]
            self.assertEqual(
                (part.length.min(), part.length.max()),
                (full.length.min(), full.length.max()),
            )

    def test_file_roundtrip_preserves_axes_metadata_and_exact_fp32(self):
        torch.manual_seed(4)
        values = torch.randn(33, 5, 1280)
        with TemporaryDirectory() as tmp:
            original, rounded = roundtrip(values, Path(tmp), {"record_id": "synthetic"})
            self.assertEqual(list(Path(tmp).iterdir()), [])
        self.assertTrue(torch.equal(original, values))
        self.assertTrue(torch.equal(rounded[7], values[7].half().float()))
        self.assertFalse(torch.equal(rounded, values))

    def test_overflow_fails(self):
        with TemporaryDirectory() as tmp, self.assertRaises(AssertionError):
            roundtrip(torch.full((33, 1, 1280), 70000.0), Path(tmp), {})

    def test_identical_inputs_have_identical_sae_support(self):
        x = torch.randn(20, 8)
        result = compare(self.sae, x, x.clone())
        self.assertEqual(result["mean_active_jaccard"], 1)
        self.assertEqual(result["latent_relative_l2"], 0)
        self.assertTrue(result["within_engineering_tolerances"])

    def test_large_change_is_not_called_compatible(self):
        x = torch.randn(20, 8)
        self.assertFalse(compare(self.sae, x, -x)["within_engineering_tolerances"])


if __name__ == "__main__":
    unittest.main()
