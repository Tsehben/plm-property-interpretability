"""Engineering tests; synthetic tensors are never reported as ESM results."""

import os
from types import SimpleNamespace
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch, MagicMock
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score, average_precision_score
from common import read_local, composition, token, write_json
from extract import (
    pool_hidden,
    tokenize,
    save_shard,
    verify_shard,
    batch_indices,
    plans,
)
from probe import bootstrap_weights, weighted_metrics
from hf_storage import private_api, commit_files, git_blob_sha, restore


class Stage1Tests(unittest.TestCase):
    def test_hf_resume_skips_identical_and_refuses_conflicts(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "extraction.json"
            path.write_text("{}")
            api = MagicMock()
            info = SimpleNamespace(sha="commit", private=True)
            api.repo_info.return_value = info
            remote = SimpleNamespace(
                path="esm_stage1/extraction.json", lfs=None, blob_id=git_blob_sha(path)
            )
            api.get_paths_info.return_value = [remote]
            with patch("hf_storage.private_api", return_value=(api, info)):
                receipt = commit_files("owner/private", {"extraction.json": path})
                self.assertEqual(receipt["uploaded_files"], 0)
                api.create_commit.assert_not_called()
                remote.blob_id = "different"
                with self.assertRaisesRegex(RuntimeError, "no overwrite"):
                    commit_files("owner/private", {"extraction.json": path})
                api.create_commit.assert_not_called()

    def test_hf_upload_is_atomic_and_failure_stops(self):
        with tempfile.TemporaryDirectory() as temp:
            paths = {
                name: Path(temp) / name
                for name in ["shard_0000.npz", "shard_0000.json"]
            }
            for path in paths.values():
                path.write_bytes(b"synthetic")
            api = MagicMock()
            info = SimpleNamespace(sha="parent", private=True)
            api.repo_info.return_value = info
            api.get_paths_info.return_value = []
            api.create_commit.return_value = SimpleNamespace(oid="new-commit")
            with patch("hf_storage.private_api", return_value=(api, info)):
                receipt = commit_files("owner/private", paths)
                self.assertEqual(receipt["uploaded_files"], 2)
                self.assertEqual(
                    len(api.create_commit.call_args.kwargs["operations"]), 2
                )
                self.assertEqual(
                    api.create_commit.call_args.kwargs["parent_commit"], "parent"
                )
                api.create_commit.side_effect = RuntimeError(
                    "credential-like-secret-error"
                )
                with self.assertRaises(RuntimeError) as error:
                    commit_files("owner/private", paths)
                self.assertNotIn("credential-like-secret-error", str(error.exception))

    def test_hf_restores_partial_run_without_completed_manifest(self):
        d = pd.DataFrame(
            {
                "record_id": ["one", "two"],
                "sha256": ["a", "b"],
                "split": ["train", "validation"],
                "length": [2, 3],
            }
        )
        with (
            tempfile.TemporaryDirectory() as source,
            tempfile.TemporaryDirectory() as dest,
        ):
            source = Path(source)
            save_shard(
                source, "shard_0000", d, np.zeros((2, 33, 1280), dtype=np.float32)
            )
            write_json(source / "extraction.json", {"shard_size": 2})
            api = MagicMock()
            api.list_repo_files.return_value = [
                "esm_stage1/" + p.name for p in source.iterdir()
            ]

            def fetch(repo, filename, **kwargs):
                self.assertEqual(kwargs["revision"], "frozen-commit")
                return str(source / Path(filename).name)

            with (
                patch(
                    "hf_storage.private_api",
                    return_value=(api, SimpleNamespace(sha="frozen-commit")),
                ),
                patch("hf_storage.token", return_value="synthetic"),
                patch("hf_storage.read_extraction", return_value=d),
                patch("hf_storage.verify_config", return_value={"shard_size": 2}),
                patch("hf_storage.safe_config"),
                patch("huggingface_hub.hf_hub_download", side_effect=fetch),
            ):
                restore("owner/private", dest)
            verify_shard(dest, "shard_0000", d)
            self.assertFalse((Path(dest) / "completed.json").exists())

    def test_installed_transformers_indexing_on_tiny_random_model(self):
        # API/indexing test only: random tiny weights are not 650M model results.
        from transformers import EsmConfig, EsmModel

        torch.manual_seed(0)
        config = EsmConfig(
            vocab_size=33,
            hidden_size=16,
            num_hidden_layers=2,
            num_attention_heads=2,
            intermediate_size=32,
            pad_token_id=1,
            mask_token_id=32,
            token_dropout=False,
            position_embedding_type="rotary",
        )
        model = EsmModel(config, add_pooling_layer=False).requires_grad_(False).eval()
        seen = {}
        handles = [
            block.register_forward_hook(
                lambda module, args, output, layer=i: seen.update(
                    {layer: output if isinstance(output, torch.Tensor) else output[0]}
                )
            )
            for i, block in enumerate(model.encoder.layer)
        ]
        with torch.inference_mode():
            output = model(
                input_ids=torch.tensor([[0, 4, 5, 2]]), output_hidden_states=True
            )
        self.assertEqual(len(output.hidden_states), 3)
        torch.testing.assert_close(output.hidden_states[1], seen[0])
        torch.testing.assert_close(
            output.hidden_states[2], model.encoder.emb_layer_norm_after(seen[1])
        )
        for handle in handles:
            handle.remove()

    def test_frozen_data_and_folds(self):
        d = read_local()
        self.assertEqual(len(d), 3052)
        self.assertEqual(d.groupby("audit_cluster").split.nunique().max(), 1)
        self.assertEqual(
            d[d.split == "train"]
            .groupby("audit_cluster")
            .inner_cv_fold.nunique()
            .max(),
            1,
        )

    def test_X_not_imputed(self):
        values = composition(["AX", "AC"])
        self.assertEqual(values[0, :20].sum(), 0.5)
        self.assertEqual(values[1, :20].sum(), 1.0)

    def test_pooling_excludes_special_and_padding(self):
        h = torch.tensor(
            [
                [[999.0], [2.0], [4.0], [999.0], [999.0]],
                [[999.0], [6.0], [999.0], [999.0], [999.0]],
            ]
        )
        actual = pool_hidden([h, h + 1], [2, 1]).numpy()
        np.testing.assert_array_equal(actual[:, :, 0], [[3, 4], [6, 7]])

    def test_tokenization_keeps_X(self):
        class Tokenizer:
            cls_token_id = 0
            eos_token_id = 2
            pad_token_id = 1
            unk_token_id = 3

            def convert_tokens_to_ids(self, tokens):
                return [{"A": 4, "X": 5}[c] for c in tokens]

            def __call__(self, seqs, **kwargs):
                n = max(map(len, seqs)) + 2
                ids = [
                    [0] + self.convert_tokens_to_ids(s) + [2] + [1] * (n - len(s) - 2)
                    for s in seqs
                ]
                return {
                    "input_ids": torch.tensor(ids),
                    "attention_mask": torch.tensor(
                        [[int(i != 1) for i in row] for row in ids]
                    ),
                }

        self.assertEqual(tokenize(Tokenizer(), ["AX", "A"])["input_ids"][0, 2], 5)

    def test_group_bootstrap_keeps_members(self):
        w = bootstrap_weights(np.array(["a", "a", "b", "c"]), repetitions=100)
        np.testing.assert_array_equal(w[:, 0], w[:, 1])
        self.assertTrue(np.all(w[:, [0, 2, 3]].sum(axis=1) == 3))

    def test_bootstrap_metrics_match_sklearn_with_ties(self):
        y = np.array([0, 1, 0, 1, 1, 0])
        scores = np.array([0.1, 0.1, 0.3, 0.5, 0.5, 0.8])
        weights = np.array([[1, 1, 1, 1, 1, 1], [2, 2, 0, 1, 1, 3], [0, 0, 0, 1, 1, 0]])
        auc, ap = weighted_metrics(y, scores, weights)
        for i in range(2):
            self.assertAlmostEqual(
                auc[i], roc_auc_score(y, scores, sample_weight=weights[i])
            )
            self.assertAlmostEqual(
                ap[i], average_precision_score(y, scores, sample_weight=weights[i])
            )
        self.assertTrue(np.isnan(auc[2]))

    def test_shards_reject_corruption_and_reordering(self):
        d = pd.DataFrame(
            {
                "record_id": ["one", "two"],
                "sha256": ["a", "b"],
                "split": ["train", "validation"],
            }
        )
        values = np.zeros((2, 33, 1280), dtype=np.float32)
        with tempfile.TemporaryDirectory() as temp:
            save_shard(temp, "shard_0000", d, values)
            with self.assertRaises(AssertionError):
                verify_shard(temp, "shard_0000", d.iloc[::-1])
            path = Path(temp) / "shard_0000.npz"
            path.write_bytes(path.read_bytes() + b"corruption")
            with self.assertRaises(AssertionError):
                verify_shard(temp, "shard_0000", d)

    def test_bucketing_unique_and_bounded(self):
        lengths = [10, 11, 30, 80]
        batches = list(batch_indices(lengths, 2, 100))
        self.assertEqual(sum(batches, []), list(range(4)))
        for b in batches:
            self.assertLessEqual((max(lengths[i] for i in b) + 2) * len(b), 100)
        d = pd.DataFrame({"record_id": ["z", "b", "a"], "length": [8, 3, 3]})
        np.testing.assert_array_equal(np.concatenate(plans(d, 2)), [2, 1, 0])

    def test_credentials_and_public_repo_refusal(self):
        with patch.dict(os.environ, {}, clear=True):
            with self.assertRaisesRegex(RuntimeError, "HF_TOKEN_PLM"):
                token()
        with (
            patch.dict(os.environ, {"HF_TOKEN_PLM": "synthetic-test-value"}),
            patch("huggingface_hub.HfApi") as api,
        ):
            api.return_value.repo_info.return_value.private = False
            with self.assertRaisesRegex(ValueError, "non-private"):
                private_api("test/repo")
            api.return_value.create_commit.assert_not_called()


if __name__ == "__main__":
    unittest.main(verbosity=2)
