"""Synthetic engineering checks only; no fabricated scientific results are retained."""

from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase, main
from unittest.mock import patch
import io
import json
import zipfile

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from common import write_json
from probe import fit_or_reuse, bootstrap_weights
from plot_train_validation import plot
from run_layerwise_auroc import make_bundle, unpack_bundle


class LayerPlotTests(TestCase):
    def data(self):
        n = 28
        return pd.DataFrame(
            dict(
                record_id=[f"synthetic{i}" for i in range(n)],
                sha256=[f"hash{i}" for i in range(n)],
                split=["train"] * 20 + ["validation"] * 8,
                label=np.arange(n) % 2,
                audit_cluster=np.arange(n),
                inner_cv_fold=[i // 4 for i in range(20)] + [-1] * 8,
            )
        )

    def test_real_small_fit_cache_and_training_scores(self):
        d = self.data()
        X = np.random.default_rng(4).normal(size=(len(d), 6))
        weights = bootstrap_weights(
            d.loc[d.split == "validation", "audit_cluster"], repetitions=20
        )
        with TemporaryDirectory() as tmp:
            out = Path(tmp)
            write_json(out / "probe_protocol.json", {"synthetic": True})
            metric = fit_or_reuse(X, d, "layer_01", out, weights, resume=True)
            fitted = joblib.load(out / "models/layer_01.joblib")
            self.assertAlmostEqual(
                metric["training_fit_auroc"],
                roc_auc_score(d.label[:20], fitted.decision_function(X[:20])),
            )
            self.assertEqual(
                len(pd.read_csv(out / "layer_01_training_oof.tsv", sep="\t")), 20
            )
            with patch("probe.fit_probe") as fit:
                self.assertEqual(
                    metric, fit_or_reuse(X, d, "layer_01", out, weights, resume=True)
                )
                fit.assert_not_called()
                (out / "layer_01_metrics.json").write_text("corrupt")
                with self.assertRaises(AssertionError):
                    fit_or_reuse(X, d, "layer_01", out, weights, resume=True)
                fit.assert_not_called()

    def test_published_bundle_omits_labels_and_detects_corruption(self):
        with TemporaryDirectory() as tmp, TemporaryDirectory() as dest:
            root = Path(tmp)
            self.data().to_csv(
                root / "layer_01_training_oof.tsv", sep="\t", index=False
            )
            write_json(root / "layer_01_complete.json", {"local_only": True})
            raw = make_bundle(root)
            unpack_bundle(raw, dest)
            table = pd.read_csv(Path(dest) / "layer_01_training_oof.tsv", sep="\t")
            self.assertTrue(
                {"label", "audit_cluster", "inner_cv_fold"}.isdisjoint(table.columns)
            )
            self.assertEqual(table.record_id.tolist(), self.data().record_id.tolist())
            with zipfile.ZipFile(io.BytesIO(raw)) as z:
                self.assertNotIn("layer_01_complete.json", z.namelist())
            bad = io.BytesIO()
            with zipfile.ZipFile(bad, "w") as z:
                z.writestr(
                    "bundle_checksums.json", json.dumps({"file.tsv": "incorrect"})
                )
                z.writestr("file.tsv", "broken")
            with self.assertRaises(AssertionError):
                unpack_bundle(bad.getvalue(), dest)

    def test_plot_requires_every_layer_and_preserves_both_train_estimands(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = [
                dict(
                    model="composition" if i == 0 else f"layer_{i:02}",
                    layer=i,
                    C=0.1,
                    training_n=2361,
                    validation_n=691,
                    cv_mean=0.65,
                    cv_sd=0.02,
                    training_oof_auroc=0.64,
                    training_fit_auroc=0.8,
                    auroc=0.63,
                    auroc_low=0.59,
                    auroc_high=0.67,
                    average_precision=0.75,
                )
                for i in range(34)
            ]
            pd.DataFrame(rows[:-1]).to_csv(root / "metrics.tsv", sep="\t", index=False)
            with self.assertRaises(AssertionError):
                plot(root)
            pd.DataFrame(rows).to_csv(root / "metrics.tsv", sep="\t", index=False)
            table = plot(root)
            self.assertEqual(table.layer.tolist(), list(range(1, 34)))
            self.assertEqual(len(list((root / "figures").iterdir())), 4)
            self.assertIn(
                "not unbiased nested-CV",
                (root / "TRAIN_VALIDATION_AUROC.md").read_text(),
            )


if __name__ == "__main__":
    main()
