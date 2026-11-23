"""Checks use development files and creation-time aggregates; never open test rows."""

import json
import unittest
from pathlib import Path
import nbformat
from analysis_helpers import read, seqhash, verify_frozen_snapshot, ReadGuard

ROOT = Path(__file__).resolve().parent
R = ROOT / "results"


class ValidationTests(unittest.TestCase):
    def test_notebook_executed_and_figures_embedded(self):
        nb = nbformat.read(ROOT / "NESG_Dataset_Validation.ipynb", as_version=4)
        nbformat.validate(nb)
        code = [c for c in nb.cells if c.cell_type == "code"]
        self.assertEqual(len(code), 13)
        self.assertTrue(all(c.execution_count is not None for c in code))
        outputs = [o for c in code for o in c.outputs]
        self.assertFalse(any(o.output_type == "error" for o in outputs))
        self.assertEqual(sum("image/png" in o.get("data", {}) for o in outputs), 10)

    def test_frozen_snapshots(self):
        lock = verify_frozen_snapshot(R)
        self.assertEqual(lock["status"], "LOCKED")
        self.assertAlmostEqual(lock["test_fraction"], 710 / 3762)
        self.assertIn(
            "reserved_test/final_cohort_and_split_manifest.tsv",
            lock["restricted_checksums"],
        )

    def test_assignment_integrity(self):
        d = read(R / "split_assignments.tsv")
        self.assertEqual(len(d), 3762)
        self.assertTrue(d.record_id.is_unique)
        self.assertTrue(d.sha256.is_unique)
        self.assertEqual(d.groupby("audit_cluster").split.nunique().max(), 1)
        self.assertEqual(
            d.split.value_counts().to_dict(),
            {"train": 2361, "reserved_test": 710, "validation": 691},
        )
        train = d[d.split == "train"]
        self.assertEqual(
            train.groupby("audit_cluster").inner_cv_fold.nunique().max(), 1
        )
        self.assertTrue(d.loc[d.split != "train", "inner_cv_fold"].eq(-1).all())

    def test_development_scores_sequences(self):
        d = read(R / "development_records.tsv.gz")
        self.assertEqual(len(d), 3052)
        self.assertEqual(set(d.split), {"train", "validation"})
        self.assertTrue((d.sequence.map(seqhash) == d.sha256).all())
        self.assertTrue(d.exp.isin([4, 5]).all())
        self.assertTrue(d.sol.isin([0, 4, 5]).all())
        self.assertTrue((d.label == (d.sol >= 4).astype(int)).all())

    def test_creation_time_detected_edge_check(self):
        x = json.loads((R / "integrity_checks.json").read_text())
        self.assertEqual(x["exact_cross_split_sequences"], 0)
        self.assertEqual(x["cross_split_groups"], 0)
        self.assertEqual(x["cross_split_detected_edges"], 0)
        self.assertEqual(x["all_union_edges_verified"], 447527)
        self.assertGreater(x["cohort_direct_edges"], 3762)

    def test_family_heldout_training_only(self):
        d = read(R / "training_family_heldout_feasibility_manifest.tsv")
        t = read(R / "training_manifest.tsv")
        self.assertEqual(set(d.record_id), set(t.record_id))
        domains = read(R / "development_domains.tsv.gz")
        x = d.merge(domains[["sha256", "pfam"]].drop_duplicates(), on="sha256")
        self.assertEqual(x.groupby("pfam").family_cv_fold.nunique().max(), 1)
        self.assertEqual(d.groupby("audit_cluster").family_cv_fold.nunique().max(), 1)

    def test_controlled_family_counts(self):
        d = read(R / "within_family_feasibility.tsv").set_index(
            ["annotation_policy", "training_groups_per_class"]
        )
        self.assertEqual(d.loc[("all_GA", 3), "validation_both_families"], 16)
        self.assertEqual(d.loc[("all_GA", 5), "validation_both_families"], 8)
        self.assertEqual(
            d.loc[("no_competing_overlap", 3), "validation_both_families"], 1
        )
        self.assertEqual(
            d.loc[("no_competing_overlap", 5), "validation_both_families"], 0
        )

    def test_locked_execution_access_log(self):
        accesses = json.loads((R / "notebook_project_file_access.json").read_text())
        self.assertGreater(len(accesses), 0)
        for entry in accesses:
            self.assertTrue(entry["path"].startswith("nesg_primary/"))
            self.assertNotIn("/reserved_test/", entry["path"])
            self.assertNotIn("/restricted_prelock/", entry["path"])

    def test_guard_rejects_paths_without_opening_them(self):
        guard = ReadGuard.__new__(ReadGuard)
        guard.repo = ROOT.parent
        guard.root = ROOT
        guard.locked = True
        guard.accesses = []
        for path in [
            R / "reserved_test/test_manifest.tsv",
            R / "restricted_prelock/raw_record_decisions.tsv.gz",
            ROOT.parent / "dataset_comparison/results/records.tsv.gz",
        ]:
            with self.assertRaises(PermissionError):
                guard.check("open", (str(path), "r", 0))


if __name__ == "__main__":
    unittest.main()
