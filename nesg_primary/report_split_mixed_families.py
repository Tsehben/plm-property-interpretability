"""User-authorized held-out-label inspection: per-split mixed Pfam family tables.

This script intentionally reads reserved-test labels. It does not modify cohorts,
annotations, groups, split assignments, model choices, or hypotheses. Do not use it
as a routine locked-test notebook step; each execution is an inspection.
"""

from pathlib import Path
from datetime import datetime, timezone
import hashlib
import itertools
import json
import pandas as pd

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
RESULTS = ROOT / "results"
OUT = RESULTS / "split_mixed_label_families"
FIELDS = [
    "pfam_accession",
    "pfam_name",
    "biological_family_name",
    "unique_soluble_proteins",
    "unique_insoluble_proteins",
    "total_unique_proteins",
    "homology_groups_soluble",
    "homology_groups_insoluble",
    "homology_groups_shared_between_classes",
    "has_competing_overlapping_assignments",
    "unique_proteins_with_competing_assignments",
]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def summarize(manifest, domains, retained=None):
    """One record per protein/family. Overlap flags always describe original GA hits."""
    if retained is None:
        retained = domains
    evidence = (
        domains.groupby(["sha256", "pfam"])
        .competing_overlap.any()
        .rename("original_competing_overlap")
        .reset_index()
    )
    meta = domains[["pfam", "family_name", "description"]].drop_duplicates()
    assert meta.pfam.is_unique
    memberships = retained[["sha256", "pfam"]].drop_duplicates()
    x = manifest[["record_id", "sha256", "label", "audit_cluster"]].merge(
        memberships, on="sha256", validate="one_to_many"
    )
    x = x.merge(evidence, on=["sha256", "pfam"], validate="one_to_one").merge(
        meta, on="pfam", validate="many_to_one"
    )
    assert not x.duplicated(["record_id", "pfam"]).any()
    rows = []
    for fam, g in x.groupby("pfam", sort=False):
        sol = g[g.label == 1]
        ins = g[g.label == 0]
        if sol.empty or ins.empty:
            continue
        gs = set(sol.audit_cluster)
        gi = set(ins.audit_cluster)
        rows.append(
            dict(
                pfam_accession=fam,
                pfam_name=g.family_name.iloc[0],
                biological_family_name=g.description.iloc[0],
                unique_soluble_proteins=sol.record_id.nunique(),
                unique_insoluble_proteins=ins.record_id.nunique(),
                total_unique_proteins=g.record_id.nunique(),
                homology_groups_soluble=len(gs),
                homology_groups_insoluble=len(gi),
                homology_groups_shared_between_classes=len(gs & gi),
                has_competing_overlapping_assignments=bool(
                    g.original_competing_overlap.any()
                ),
                unique_proteins_with_competing_assignments=g.loc[
                    g.original_competing_overlap, "record_id"
                ].nunique(),
            )
        )
    table = (
        pd.DataFrame(rows, columns=FIELDS)
        .sort_values(
            ["total_unique_proteins", "pfam_accession"], ascending=[False, True]
        )
        .reset_index(drop=True)
    )
    assert (
        table.total_unique_proteins
        == table.unique_soluble_proteins + table.unique_insoluble_proteins
    ).all()
    assert (
        (table[["unique_soluble_proteins", "unique_insoluble_proteins"]] >= 1)
        .all()
        .all()
    )
    assert (table.homology_groups_soluble <= table.unique_soluble_proteins).all()
    assert (table.homology_groups_insoluble <= table.unique_insoluble_proteins).all()
    return table, x


def markdown(table):
    columns = [
        "pfam_accession",
        "biological_family_name",
        "unique_soluble_proteins",
        "unique_insoluble_proteins",
        "total_unique_proteins",
        "homology_groups_soluble",
        "homology_groups_insoluble",
        "has_competing_overlapping_assignments",
        "unique_proteins_with_competing_assignments",
    ]
    headers = [
        "Pfam",
        "Biological family/domain name",
        "Soluble",
        "Insoluble",
        "Total",
        "Groups: soluble",
        "Groups: insoluble",
        "Competing overlap*",
        "Affected proteins*",
    ]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---", "---"] + ["---:"] * 5 + ["---", "---:"]) + " |",
    ]
    for r in table[columns].itertuples(index=False, name=None):
        vals = [
            ("Yes" if v else "No")
            if isinstance(v, bool)
            else str(v).replace("|", "/").replace("\n", " ")
            for v in r
        ]
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines) + "\n"


def main():
    OUT.mkdir(exist_ok=True)
    paths = {
        "train_manifest": RESULTS / "training_manifest.tsv",
        "validation_manifest": RESULTS / "validation_manifest.tsv",
        "reserved_test_manifest": RESULTS / "reserved_test/test_manifest.tsv",
        "domains": REPO / "dataset_comparison/results/domains_with_clans.tsv.gz",
        "test_lock": RESULTS / "test_lock.json",
    }
    before = {k: digest(p) for k, p in paths.items()}
    lock = json.loads(paths["test_lock"].read_text())
    assert (
        before["train_manifest"]
        == lock["rerun_safe_checksums"]["training_manifest.tsv"]
    )
    assert (
        before["validation_manifest"]
        == lock["rerun_safe_checksums"]["validation_manifest.tsv"]
    )
    assert (
        before["reserved_test_manifest"]
        == lock["restricted_checksums"]["reserved_test/test_manifest.tsv"]
    )
    pins = json.loads(
        (REPO / "dataset_comparison/results/output_checksums.json").read_text()
    )
    assert before["domains"] == pins["results/domains_with_clans.tsv.gz"]
    domains = pd.read_csv(paths["domains"], sep="\t", keep_default_na=False)
    assert domains.competing_overlap.isin([True, False]).all()
    tables = {}
    summary = []
    sensitivity_comparison = []
    manifests = {}
    for split, N, nsol, nins in [
        ("train", 2361, 1645, 716),
        ("validation", 691, 491, 200),
        ("reserved_test", 710, 478, 232),
    ]:
        m = pd.read_csv(paths[split + "_manifest"], sep="\t", keep_default_na=False)
        assert len(m) == N and m.label.value_counts().to_dict() == {1: nsol, 0: nins}
        assert set(m.split) == {split} and m.record_id.is_unique and m.sha256.is_unique
        assert (
            m.audit_cluster.ne("").all()
            and m.exp.isin([4, 5]).all()
            and m.sol.isin([0, 4, 5]).all()
        )
        assert (m.label == (m.sol >= 4).astype(int)).all()
        manifests[split] = m
        d = domains[domains.sha256.isin(m.sha256)].copy()
        assert (d.sequence_length == d.sha256.map(m.set_index("sha256").length)).all()
        assert (
            (d.sequence_score >= d.sequence_ga) & (d.domain_score >= d.domain_ga)
        ).all()
        assert (
            (d.ali_from >= 1)
            & (d.ali_from <= d.ali_to)
            & (d.ali_to <= d.sequence_length)
        ).all()
        # Verify existing overlap flags from recorded coordinates; no domain search is repeated.
        recalculated = set()
        for _, g in d[d.clan != ""].groupby(["sha256", "clan"]):
            for i, j in itertools.combinations(g.index, 2):
                a, b = d.loc[i], d.loc[j]
                if a.pfam != b.pfam and max(a.ali_from, b.ali_from) <= min(
                    a.ali_to, b.ali_to
                ):
                    recalculated.update([i, j])
        assert set(d.index[d.competing_overlap]) == recalculated
        for policy, retained in [
            ("all_GA", d),
            ("no_competing_overlap", d[~d.competing_overlap]),
        ]:
            table, membership = summarize(m, d, retained)
            tables[(split, policy)] = table
            stem = f"{split}_{policy}_mixed_families"
            table.to_csv(OUT / (stem + ".csv"), index=False)
            table.to_csv(OUT / (stem + ".tsv"), index=False, sep="\t")
            eligible = membership[membership.pfam.isin(table.pfam_accession)]
            summary.append(
                dict(
                    split=split,
                    annotation_policy=policy,
                    split_unique_proteins=len(m),
                    split_soluble=nsol,
                    split_insoluble=nins,
                    qualifying_families=len(table),
                    unique_proteins_in_qualifying_families=eligible.record_id.nunique(),
                    soluble_in_qualifying_families=eligible.loc[
                        eligible.label == 1, "record_id"
                    ].nunique(),
                    insoluble_in_qualifying_families=eligible.loc[
                        eligible.label == 0, "record_id"
                    ].nunique(),
                    qualifying_families_with_original_competing_hits=int(
                        table.has_competing_overlapping_assignments.sum()
                    ),
                )
            )
        primary = tables[(split, "all_GA")]
        sens = tables[(split, "no_competing_overlap")]
        assert set(sens.pfam_accession) <= set(primary.pfam_accession)
        compared = primary[
            [
                "pfam_accession",
                "biological_family_name",
                "unique_soluble_proteins",
                "unique_insoluble_proteins",
                "total_unique_proteins",
            ]
        ].merge(
            sens[
                [
                    "pfam_accession",
                    "unique_soluble_proteins",
                    "unique_insoluble_proteins",
                    "total_unique_proteins",
                ]
            ],
            on="pfam_accession",
            how="left",
            suffixes=("_all_GA", "_sensitivity_mixed_only"),
        )
        compared["split"] = split
        compared["remains_mixed_after_overlap_removal"] = compared.pfam_accession.isin(
            sens.pfam_accession
        )
        # Exact sensitivity counts, including families no longer mixed: do not misreport missing rows as no hits.
        sensitivity_members = (
            d[~d.competing_overlap][["sha256", "pfam"]]
            .drop_duplicates()
            .merge(m[["sha256", "label"]], on="sha256")
        )
        for lab, col in [
            (1, "sensitivity_soluble_including_nonmixed"),
            (0, "sensitivity_insoluble_including_nonmixed"),
        ]:
            counts = (
                sensitivity_members[sensitivity_members.label == lab]
                .groupby("pfam")
                .sha256.nunique()
            )
            compared[col] = compared.pfam_accession.map(counts).fillna(0).astype(int)
        assert (
            compared.sensitivity_soluble_including_nonmixed
            <= compared.unique_soluble_proteins_all_GA
        ).all()
        assert (
            compared.sensitivity_insoluble_including_nonmixed
            <= compared.unique_insoluble_proteins_all_GA
        ).all()
        sensitivity_comparison.append(compared)
    for a, b in itertools.combinations(manifests.values(), 2):
        assert not set(a.sha256) & set(b.sha256)
        assert not set(a.audit_cluster) & set(b.audit_cluster)
    for k, p in paths.items():
        assert digest(p) == before[k], f"Input modified: {k}"
    summary = pd.DataFrame(summary)
    summary.to_csv(OUT / "split_summary.tsv", sep="\t", index=False)
    summary.to_csv(OUT / "split_summary.csv", index=False)
    pd.concat(sensitivity_comparison, ignore_index=True).to_csv(
        OUT / "overlap_sensitivity_comparison.tsv", sep="\t", index=False
    )
    lines = [
        "# Split-specific mixed-label Pfam families",
        "",
        "**Held-out-label inspection disclosure.** The user explicitly requested this analysis of reserved-test labels. Producing and viewing the test-label-by-family tables constitutes inspection of held-out labels. The test can no longer be described as uninspected with respect to family–label composition. This is not model-performance evaluation, and these results must not be used to revise the cohort, split assignment, model selection or experimental hypotheses. No such changes were made. The UESolDS independent test was not accessed.",
        "",
        "The original `results/test_lock.json` is preserved as a historical freeze record; its unchanged `LOCKED` field does not negate this subsequent authorized inspection. Access restrictions for ordinary notebook reruns remain in place. This standalone inspection is recorded in `TEST_LABEL_EXPOSURE.md` and its machine-readable inspection event.",
        "",
        "## Definitions and validation",
        "",
        "- Each split is analyzed independently. A qualifying family has ≥1 soluble and ≥1 insoluble protein **in that split**; no cross-split support rule, ≥3/≥5 threshold or model result is used. Training is evaluated independently by the same definition.",
        "- Soluble is label 1; insoluble is label 0. Protein counts deduplicate stable record ID and Pfam accession. Each split has unique sequence hashes, so these also equal unique sequences. Repeat domains/assignments do not inflate counts.",
        "- Family names are the Pfam 37.0 biological descriptions. Full CSV/TSV files also contain the short Pfam name. These are domain families, not necessarily identical full-protein architectures.",
        "- Homology counts use the existing audit MMseqs2 component IDs. A group containing both labels counts once in each class; the CSV/TSV tables additionally report shared groups. Families, clans and homology groups remain distinct.",
        "- *Competing overlap* means a hit for this family overlaps a different-family hit from the same nonempty Pfam clan on the same protein (1-based inclusive aligned coordinates). Existing flags were checked directly against coordinates.",
        "- The primary tables retain **all GA-passing assignments**. The separate sensitivity removes both competing hits, then deduplicates memberships and independently reassesses mixed status. An unflagged repeat of the same family can preserve membership. No proteins or labels are removed from the cohort.",
        "- *Affected proteins* counts unique proteins in that table with an original competing hit for the reported family. Sensitivity-table flags retain the original evidence among surviving memberships; they do not imply that flagged hits were retained.",
        "- Tables are complete, sorted by total unique proteins descending, then accession. Summing rows can count multidomain proteins more than once; split-summary protein totals use unions.",
        "- Input SHA256 matches the frozen manifests and pinned domain artifact. Split sizes/labels, unique IDs/hashes, coordinate bounds, sequence lengths, both GA thresholds, overlap flags and absence of group/hash overlap between every pair of splits were asserted. Input hashes are unchanged after execution.",
        "",
        "## Counts",
        "",
        "| Split | Split proteins (soluble / insoluble) | All-GA mixed families | Mixed families after overlap removal |",
        "| --- | ---: | ---: | ---: |",
    ]
    for split in ["train", "validation", "reserved_test"]:
        row = summary[
            (summary.split == split) & (summary.annotation_policy == "all_GA")
        ].iloc[0]
        lines.append(
            f"| {split} | {row.split_unique_proteins} ({row.split_soluble} / {row.split_insoluble}) | {len(tables[(split, 'all_GA')])} | {len(tables[(split, 'no_competing_overlap')])} |"
        )
    for policy, title in [
        ("all_GA", "Primary: all GA-passing assignments"),
        ("no_competing_overlap", "Sensitivity: remove competing overlapping hits"),
    ]:
        lines += ["", f"## {title}", ""]
        for split in ["train", "validation", "reserved_test"]:
            table = tables[(split, policy)]
            stem = f"{split}_{policy}_mixed_families"
            prefix = "results/split_mixed_label_families/"
            lines += [
                f"### {split}: {len(table)} qualifying families",
                "",
                f"[Complete CSV]({prefix}{stem}.csv) · [Complete TSV]({prefix}{stem}.tsv)",
                "",
                markdown(table),
            ]
    lines += [
        "",
        "## Reproduce and limits",
        "",
        "`python nesg_primary/report_split_mixed_families.py` uses existing frozen manifests and domain annotations only; it does not rerun Pfam, cluster proteins, train models or change hypotheses. **Running this script is itself an authorized-label-inspection operation**, not a normal locked-test analysis step.",
        "",
        "Counts are descriptive. A mixed family with one example of each class is included as requested; this does not establish statistical power, assay matching, architectural identity or mechanistic evidence. The separate `overlap_sensitivity_comparison.tsv` includes exact surviving counts even when a family ceases to be mixed.",
    ]
    report = ROOT / "SPLIT_MIXED_LABEL_FAMILIES.md"
    report.write_text("\n".join(lines) + "\n")
    event = dict(
        timestamp_UTC=datetime.now(timezone.utc).isoformat(),
        authorization="Explicit user request for validation and reserved-test mixed-label Pfam tables, subsequently extended to training",
        exposure="Reserved-test labels inspected and reported by Pfam family; no longer uninspected for family-label composition",
        allowed_use="Descriptive reporting only; no cohort/split/model-selection/hypothesis revisions",
        inputs={
            k: {"path": str(p.relative_to(REPO)), "sha256": before[k]}
            for k, p in paths.items()
        },
        unchanged_inputs=True,
        UESolDS_independent_test_accessed=False,
        annotations_or_clusters_rerun=False,
        training_support_used=False,
        outputs={
            str(p.relative_to(ROOT)): digest(p)
            for p in sorted(OUT.glob("*"))
            if p.suffix in [".csv", ".tsv"]
        },
        report_sha256=digest(report),
        script_sha256=digest(Path(__file__)),
        summary=summary.to_dict("records"),
        pandas_version=pd.__version__,
    )
    (OUT / "inspection_event.json").write_text(json.dumps(event, indent=2) + "\n")
    with (ROOT / "TEST_LABEL_INSPECTIONS.jsonl").open("a") as f:
        f.write(
            json.dumps(
                {
                    "timestamp_UTC": event["timestamp_UTC"],
                    "report": "SPLIT_MIXED_LABEL_FAMILIES.md",
                    "event": "results/split_mixed_label_families/inspection_event.json",
                    "event_sha256": digest(OUT / "inspection_event.json"),
                    "exposure": event["exposure"],
                }
            )
            + "\n"
        )
    (ROOT / "TEST_LABEL_EXPOSURE.md").write_text(
        "# Subsequent reserved-test label inspection\n\n"
        + event["timestamp_UTC"]
        + " — The user explicitly requested and authorized a split-specific mixed-label Pfam report including the reserved test. Test labels have now been inspected by family. Producing and viewing this report is held-out-label inspection, even though no model performance was evaluated.\n\nSee [the complete report](SPLIT_MIXED_LABEL_FAMILIES.md) and `TEST_LABEL_INSPECTIONS.jsonl`. The original lock file is historical and remains unchanged; it must not be taken as evidence that no later label inspection occurred. Existing cohort, split and annotation files are unchanged. These results must not be used to revise dataset membership, model selection or experimental hypotheses. The ordinary notebook lock remains operational; further test uses require their own declared scope.\n"
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
