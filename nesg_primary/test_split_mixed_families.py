"""Checks for the explicitly authorized held-out family-label report."""

import unittest
import json
from collections import defaultdict
import pandas as pd
from report_split_mixed_families import ROOT, REPO, RESULTS, OUT, summarize, digest


class SplitFamilyTests(unittest.TestCase):
    def test_repeated_domain_does_not_repeat_protein(self):
        manifest = pd.DataFrame(
            [["a", "s1", 1, "g1"], ["b", "s2", 0, "g2"]],
            columns=["record_id", "sha256", "label", "audit_cluster"],
        )
        domains = pd.DataFrame(
            [
                ["s1", "PF1", "Fam", "Biological family", False],
                ["s1", "PF1", "Fam", "Biological family", False],
                ["s2", "PF1", "Fam", "Biological family", False],
            ],
            columns=[
                "sha256",
                "pfam",
                "family_name",
                "description",
                "competing_overlap",
            ],
        )
        t, x = summarize(manifest, domains)
        self.assertEqual(len(t), 1)
        self.assertEqual(t.total_unique_proteins.iloc[0], 2)
        self.assertEqual(t.unique_soluble_proteins.iloc[0], 1)

    def test_unflagged_repeat_preserves_membership(self):
        manifest = pd.DataFrame(
            [["a", "s1", 1, "g"], ["b", "s2", 0, "g"]],
            columns=["record_id", "sha256", "label", "audit_cluster"],
        )
        domains = pd.DataFrame(
            [
                ["s1", "PF1", "Fam", "Family", True],
                ["s1", "PF1", "Fam", "Family", False],
                ["s2", "PF1", "Fam", "Family", False],
            ],
            columns=[
                "sha256",
                "pfam",
                "family_name",
                "description",
                "competing_overlap",
            ],
        )
        t, _ = summarize(manifest, domains, domains[~domains.competing_overlap])
        self.assertEqual(len(t), 1)
        self.assertTrue(t.has_competing_overlapping_assignments.iloc[0])
        self.assertEqual(t.total_unique_proteins.iloc[0], 2)
        self.assertEqual(t.homology_groups_shared_between_classes.iloc[0], 1)

    def test_complete_export_matches_independent_set_counts(self):
        domains = pd.read_csv(
            REPO / "dataset_comparison/results/domains_with_clans.tsv.gz",
            sep="\t",
            keep_default_na=False,
        )
        for split, path in [
            ("train", RESULTS / "training_manifest.tsv"),
            ("validation", RESULTS / "validation_manifest.tsv"),
            ("reserved_test", RESULTS / "reserved_test/test_manifest.tsv"),
        ]:
            manifest = pd.read_csv(path, sep="\t")
            lookup = {
                r.sha256: (r.record_id, r.label, r.audit_cluster)
                for r in manifest.itertuples()
            }
            for policy in ["all_GA", "no_competing_overlap"]:
                proteins = defaultdict(lambda: {0: set(), 1: set()})
                groups = defaultdict(lambda: {0: set(), 1: set()})
                subset = domains[domains.sha256.isin(lookup)]
                if policy == "no_competing_overlap":
                    subset = subset[~subset.competing_overlap]
                for r in subset.itertuples():
                    uid, label, group = lookup[r.sha256]
                    proteins[r.pfam][label].add(uid)
                    groups[r.pfam][label].add(group)
                expected = {
                    p: (
                        len(v[1]),
                        len(v[0]),
                        len(v[1] | v[0]),
                        len(groups[p][1]),
                        len(groups[p][0]),
                    )
                    for p, v in proteins.items()
                    if v[0] and v[1]
                }
                table = pd.read_csv(
                    OUT / f"{split}_{policy}_mixed_families.tsv", sep="\t"
                )
                actual = {
                    r.pfam_accession: (
                        r.unique_soluble_proteins,
                        r.unique_insoluble_proteins,
                        r.total_unique_proteins,
                        r.homology_groups_soluble,
                        r.homology_groups_insoluble,
                    )
                    for r in table.itertuples()
                }
                self.assertEqual(actual, expected)
                self.assertTrue(table.total_unique_proteins.is_monotonic_decreasing)
                csv = pd.read_csv(OUT / f"{split}_{policy}_mixed_families.csv")
                pd.testing.assert_frame_equal(csv, table)

    def test_inputs_and_outputs_unchanged_since_report(self):
        event = json.loads((OUT / "inspection_event.json").read_text())
        for v in event["inputs"].values():
            self.assertEqual(digest(REPO / v["path"]), v["sha256"])
        for p, h in event["outputs"].items():
            self.assertEqual(digest(ROOT / p), h)
        self.assertEqual(
            digest(ROOT / "SPLIT_MIXED_LABEL_FAMILIES.md"), event["report_sha256"]
        )


if __name__ == "__main__":
    unittest.main()
