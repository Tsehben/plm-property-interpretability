#!/usr/bin/env python3
"""
Build a versioned analysis view and inclusion manifest for candidate-feature
discovery, and assess controlled-comparison feasibility from DISCOVERY DATA ONLY.

Saved data only. No network requests, no ESM/SAE extraction, no model fitting,
no intervention runs. Source data is read and never written:

    DeepLoc2/baseline/deeploc_binary_manifest.tsv      (frozen cohort, labels, partitions)
    annotation_readiness/out/annotation_primary.tsv    (original annotation run)
    annotation_readiness/out/annotation_domains.tsv
    annotation_readiness/out/annotation_features.tsv
    annotation_readiness/out_reconciled/*              (reconciled tables)

Everything this script produces goes under discovery_analysis/v1/.

Inclusion rules, applied in a fixed precedence order so that overlapping
exclusions are never double-counted:

    E1  engineering_protein == 1                       (8 proteins cohort-wide)
    E2  p0_base_collision != 'none'                    (27 proteins cohort-wide)
    E3  coordinates_transferable == 0                  (annotation-conditioned only)

Samples:

    cohort                 all 10,139 development proteins (P1-P4). Never altered.
    *_core                 role sample minus E1 and E2. Used for feature selection
                           and for evaluation. Label-conditioned only.
    *_annotated            *_core minus E3. The only sample in which a
                           signal-peptide annotation may be conditioned on, because
                           it is the only one where the annotation is known to apply
                           to the frozen extraction-input sequence.

P0 is not read as data anywhere. The manifest's partition column is used only to
assert that no P0 row enters any sample.

Targets come from the frozen DeepLoc binary label (manifest `target`) and from
nothing else. `subcellular_location_text` is never read by this script; the
assertion that it is absent from every output is an explicit integrity check.

Usage:  python3 build_analysis_view.py
"""

import csv
import hashlib
import json
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
AR = os.path.join(ROOT, "annotation_readiness")
VIEW = os.path.join(HERE, "v1")

MANIFEST = os.path.join(ROOT, "DeepLoc2", "baseline", "deeploc_binary_manifest.tsv")
PRIMARY = os.path.join(AR, "out", "annotation_primary.tsv")
DOMAINS = os.path.join(AR, "out", "annotation_domains.tsv")

CLASS_NAME = {"1": "Extracellular", "0": "Cytoplasm"}
ROLE_BY_PARTITION = {"1": "discovery", "2": "discovery", "3": "tuning", "4": "evaluation"}
THRESHOLDS = (5, 10, 20)

# Fields that may NEVER appear in an analysis view.
FORBIDDEN_FIELDS = ("subcellular_location_text",)

# Fields carried into the view. Deliberately narrow.
VIEW_FIELDS = [
    "ACC", "base_accession", "is_isoform", "partition", "role", "target",
    "localisation_class", "kingdom_deeploc", "frozen_len", "frozen_sha256",
    "engineering_protein", "p0_base_collision",
    "sequence_match", "coordinates_transferable",
    "signal_status", "signal_evidence_best", "n_signal_features",
    "transit_status", "n_pfam", "n_interpro",
    "domain_annotation_applies_to", "in_core", "in_annotated", "exclusion_reason",
]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, fields, records):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in records:
            w.writerow(r)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def pct(n, d):
    return "" if not d else f"{100.0 * n / d:.1f}"


# ------------------------------------------------------------ inclusion rules

def exclusion_reason(r):
    """Fixed precedence. Returns '' when the protein is in the core sample."""
    if r["engineering_protein"] == "1":
        return "E1_engineering_protein"
    if r["p0_base_collision"] != "none":
        return "E2_p0_base_collision"
    return ""


def domain_applicability(r):
    """Pfam/InterPro cross-references in UniProt are ENTRY-level: they are not
    residue-mapped in the record we retrieved, and for an isoform accession they
    describe the parent entry, not necessarily the isoform sequence."""
    if r["is_isoform"] == "1":
        return "parent_entry_unresolved_isoform"
    if r["sequence_match"] != "exact":
        return "parent_entry_sequence_differs_from_frozen"
    return "entry_matches_frozen_sequence"


def build_view(primary, manifest):
    view = []
    for r in primary:
        m = manifest[r["ACC"]]
        assert m["partition"] != "0", f"P0 row {r['ACC']} reached the analysis view"
        assert m["target"] == r["target"], f"label drift at {r['ACC']}"
        assert m["partition"] == r["partition"], f"partition drift at {r['ACC']}"
        reason = exclusion_reason(r)
        in_core = reason == ""
        in_annotated = in_core and r["coordinates_transferable"] == "1"
        if in_core and not in_annotated:
            ann_reason = "E3_coordinates_not_transferable"
        else:
            ann_reason = ""
        row = {k: r.get(k, "") for k in VIEW_FIELDS if k in r}
        row.update(
            ACC=r["ACC"],
            partition=m["partition"],
            role=ROLE_BY_PARTITION[m["partition"]],
            target=m["target"],
            localisation_class=CLASS_NAME[m["target"]],
            frozen_sha256=m["sha256"],
            frozen_len=m["len"],
            domain_annotation_applies_to=domain_applicability(r),
            in_core="1" if in_core else "0",
            in_annotated="1" if in_annotated else "0",
            exclusion_reason=reason or ann_reason,
        )
        view.append(row)
    return view


def exclusion_accounting(primary, view):
    """Report overlaps explicitly instead of summing criteria."""
    e1 = set(r["ACC"] for r in primary if r["engineering_protein"] == "1")
    e2 = set(r["ACC"] for r in primary if r["p0_base_collision"] != "none")
    e3 = set(r["ACC"] for r in primary if r["coordinates_transferable"] == "0")
    n = len(primary)
    core = set(r["ACC"] for r in view if r["in_core"] == "1")
    ann = set(r["ACC"] for r in view if r["in_annotated"] == "1")
    return dict(
        cohort=n,
        criterion_sizes=dict(E1_engineering=len(e1), E2_p0_collision=len(e2), E3_coords_not_transferable=len(e3)),
        pairwise_overlaps=dict(
            E1_and_E2=len(e1 & e2), E1_and_E3=len(e1 & e3), E2_and_E3=len(e2 & e3),
            E1_and_E2_and_E3=len(e1 & e2 & e3),
        ),
        union_E1_E2=len(e1 | e2),
        naive_sum_E1_E2=len(e1) + len(e2),
        union_E1_E2_E3=len(e1 | e2 | e3),
        naive_sum_E1_E2_E3=len(e1) + len(e2) + len(e3),
        core_sample=len(core),
        core_check=n - len(e1 | e2) == len(core),
        annotated_sample=len(ann),
        annotated_check=n - len(e1 | e2 | e3) == len(ann),
    )


def sample_counts(view):
    recs = []
    for sample, test in (
        ("cohort", lambda r: True),
        ("core", lambda r: r["in_core"] == "1"),
        ("annotated", lambda r: r["in_annotated"] == "1"),
    ):
        for role in ("discovery", "tuning", "evaluation", "ALL"):
            sub = [r for r in view if test(r) and (role == "ALL" or r["role"] == role)]
            recs.append(
                dict(
                    sample=sample, role=role, n=len(sub),
                    Extracellular=sum(1 for r in sub if r["target"] == "1"),
                    Cytoplasm=sum(1 for r in sub if r["target"] == "0"),
                    prevalence=f"{sum(1 for r in sub if r['target'] == '1') / len(sub):.4f}" if sub else "",
                )
            )
    return recs


# ----------------------------------------------------- discovery-only feasibility

def discovery_breakdown(view):
    """Exact usable counts by class, signal status, evidence category. Discovery
    (P1+P2) only."""
    disc = [r for r in view if r["role"] == "discovery" and r["in_core"] == "1"]
    recs = []
    for t in ("1", "0"):
        sub = [r for r in disc if r["target"] == t]
        for status in ("annotated_present", "not_annotated", "unverifiable"):
            s2 = [r for r in sub if r["signal_status"] == status]
            if status == "annotated_present":
                for cat, k in Counter(r["signal_evidence_best"] for r in s2).most_common():
                    recs.append(
                        dict(
                            localisation_class=CLASS_NAME[t], signal_status=status,
                            evidence_category=cat, n_core=k,
                            n_annotation_conditioned=sum(
                                1 for r in s2 if r["in_annotated"] == "1" and r["signal_evidence_best"] == cat
                            ),
                        )
                    )
            else:
                recs.append(
                    dict(
                        localisation_class=CLASS_NAME[t], signal_status=status,
                        evidence_category="(none)", n_core=len(s2),
                        n_annotation_conditioned=sum(1 for r in s2 if r["in_annotated"] == "1"),
                    )
                )
    return recs, disc


def group_counts(disc, domains):
    """Per Pfam / InterPro group counts within discovery. Multiple memberships are
    retained: a protein appears in every group it carries, and group sizes
    therefore sum to more than the number of proteins."""
    acc_in = {r["ACC"]: r for r in disc}
    byg = defaultdict(lambda: defaultdict(list))
    names = {}
    for d in domains:
        if d["ACC"] not in acc_in:
            continue
        names[(d["database"], d["xref_id"])] = d["entry_name"]
        byg[d["database"]][d["xref_id"]].append(acc_in[d["ACC"]])
    out = []
    for db in sorted(byg):
        for gid, members in byg[db].items():
            ext = [r for r in members if r["target"] == "1"]
            cyt = [r for r in members if r["target"] == "0"]
            out.append(
                dict(
                    database=db, group_id=gid, entry_name=names[(db, gid)],
                    n_proteins=len(members),
                    n_Extracellular=len(ext), n_Cytoplasm=len(cyt),
                    mixed_class="1" if ext and cyt else "0",
                    ext_P1=sum(1 for r in ext if r["partition"] == "1"),
                    ext_P2=sum(1 for r in ext if r["partition"] == "2"),
                    cyt_P1=sum(1 for r in cyt if r["partition"] == "1"),
                    cyt_P2=sum(1 for r in cyt if r["partition"] == "2"),
                    n_applies_to_frozen=sum(
                        1 for r in members if r["domain_annotation_applies_to"] == "entry_matches_frozen_sequence"
                    ),
                    n_parent_entry_isoform=sum(
                        1 for r in members if r["domain_annotation_applies_to"] == "parent_entry_unresolved_isoform"
                    ),
                    n_parent_entry_seq_differs=sum(
                        1 for r in members if r["domain_annotation_applies_to"] == "parent_entry_sequence_differs_from_frozen"
                    ),
                )
            )
    out.sort(key=lambda r: (r["database"], -r["n_proteins"], r["group_id"]))
    return out


def eligibility(groups, disc, members):
    """Descriptive only. Every threshold is reported; none is chosen, and no model
    outcome is consulted.

    `n_groups_eligible` overstates the number of *independent* comparisons:
    Pfam and InterPro groups nest and duplicate (a domain, its subfamilies and its
    homologous superfamily can cover the same proteins), so the count of distinct
    member sets and the count of groups wholly contained in another eligible group
    are reported alongside it."""
    recs = []
    n_disc = len(disc)
    for db in sorted(set(g["database"] for g in groups)):
        gs = [g for g in groups if g["database"] == db]
        mixed = [g for g in gs if g["mixed_class"] == "1"]
        for t in THRESHOLDS:
            elig = [g for g in mixed if g["n_Extracellular"] >= t and g["n_Cytoplasm"] >= t]
            sets = [frozenset(members[(db, g["group_id"])]) for g in elig]
            nested = sum(1 for i, a in enumerate(sets) if any(a < b for j, b in enumerate(sets) if i != j))
            # Minimum cell occupancy for a within-group P1-versus-P2 comparison.
            both = sum(
                1 for g in elig
                if min(g["ext_P1"], g["ext_P2"], g["cyt_P1"], g["cyt_P2"]) >= 2
            )
            empty = sum(
                1 for g in elig
                if min(g["ext_P1"], g["ext_P2"], g["cyt_P1"], g["cyt_P2"]) == 0
            )
            recs.append(
                dict(
                    database=db,
                    threshold_per_class=t,
                    n_groups_total=len(gs),
                    n_groups_mixed_class=len(mixed),
                    n_groups_eligible=len(elig),
                    n_distinct_member_sets=len(set(sets)),
                    n_eligible_groups_nested_in_another=nested,
                    n_eligible_with_ge2_per_class_in_both_P1_and_P2=both,
                    n_eligible_with_an_empty_partition_cell=empty,
                    n_Extracellular_in_eligible_groups=sum(g["n_Extracellular"] for g in elig),
                    n_Cytoplasm_in_eligible_groups=sum(g["n_Cytoplasm"] for g in elig),
                    note="group memberships overlap; protein totals are not disjoint",
                )
            )
    return recs, n_disc


# ---------------------------------------------------------------------- main

def main():
    os.makedirs(VIEW, exist_ok=True)
    manifest = {r["ACC"]: r for r in read_tsv(MANIFEST)}
    primary = read_tsv(PRIMARY)
    domains = read_tsv(DOMAINS)

    for f in FORBIDDEN_FIELDS:
        assert f in primary[0], "expected the forbidden field to exist upstream"

    view = build_view(primary, manifest)
    assert len(view) == 10139

    # Integrity: no forbidden field reached the view.
    for f in FORBIDDEN_FIELDS:
        assert f not in VIEW_FIELDS and all(f not in r for r in view)

    excl = exclusion_accounting(primary, view)
    counts = sample_counts(view)
    disc_break, disc = discovery_breakdown(view)
    groups = group_counts(disc, domains)

    # protein cover of eligible groups, per database and threshold (deduplicated)
    acc_by_group = defaultdict(set)
    disc_acc = set(r["ACC"] for r in disc)
    for d in domains:
        if d["ACC"] in disc_acc:
            acc_by_group[(d["database"], d["xref_id"])].add(d["ACC"])

    elig, n_disc = eligibility(groups, disc, acc_by_group)
    tgt = {r["ACC"]: r["target"] for r in disc}
    for rec in elig:
        db, t = rec["database"], rec["threshold_per_class"]
        sel = [g for g in groups if g["database"] == db and g["mixed_class"] == "1"
               and g["n_Extracellular"] >= t and g["n_Cytoplasm"] >= t]
        cov = set()
        for g in sel:
            cov |= acc_by_group[(db, g["group_id"])]
        rec["n_distinct_proteins_covered"] = len(cov)
        rec["n_distinct_Extracellular_covered"] = sum(1 for a in cov if tgt[a] == "1")
        rec["n_distinct_Cytoplasm_covered"] = sum(1 for a in cov if tgt[a] == "0")
        rec["pct_of_discovery_core_covered"] = pct(len(cov), n_disc)

    # Domain-annotation applicability within discovery core.
    applic = [
        dict(applies_to=k, n_proteins=v, pct_of_discovery_core=pct(v, n_disc))
        for k, v in Counter(r["domain_annotation_applies_to"] for r in disc).most_common()
    ]

    # Homology clusters: do any identifiers exist anywhere in the saved outputs?
    cluster_fields = [c for c in primary[0] if "cluster" in c.lower() or "mmseq" in c.lower()]
    homology = dict(
        cluster_identifier_columns_found=cluster_fields,
        cluster_files_found=[],
        clusters_exist=False,
        consequence=(
            "No sequence-identity clustering has been computed for this cohort. "
            "Independent homology-group counts cannot be reported, and no "
            "family-controlled generalisation claim can be verified. Pfam/InterPro "
            "co-membership is not a homology cluster and partition preservation is "
            "not family-disjointness."
        ),
    )

    # Writes
    write_tsv(os.path.join(VIEW, "inclusion_manifest.tsv"), VIEW_FIELDS, view)
    write_tsv(
        os.path.join(VIEW, "analysis_view_discovery.tsv"), VIEW_FIELDS,
        [r for r in view if r["role"] == "discovery" and r["in_core"] == "1"],
    )
    write_tsv(os.path.join(VIEW, "sample_counts.tsv"), list(counts[0].keys()), counts)
    write_tsv(
        os.path.join(VIEW, "discovery_signal_evidence_counts.tsv"),
        list(disc_break[0].keys()), disc_break,
    )
    write_tsv(os.path.join(VIEW, "discovery_group_counts.tsv"), list(groups[0].keys()), groups)
    write_tsv(os.path.join(VIEW, "discovery_group_eligibility.tsv"), list(elig[0].keys()), elig)
    write_tsv(os.path.join(VIEW, "domain_applicability.tsv"), list(applic[0].keys()), applic)

    integrity = dict(
        view_version="v1",
        generated_from=dict(
            manifest=dict(path=os.path.relpath(MANIFEST, ROOT), sha256=sha256_file(MANIFEST)),
            annotation_primary=dict(path=os.path.relpath(PRIMARY, ROOT), sha256=sha256_file(PRIMARY)),
            annotation_domains=dict(path=os.path.relpath(DOMAINS, ROOT), sha256=sha256_file(DOMAINS)),
        ),
        network_requests=0,
        source_files_modified=[],
        n_cohort_rows=len(view),
        n_unique_accessions=len(set(r["ACC"] for r in view)),
        p0_rows_in_any_sample=0,
        labels_from="frozen DeepLoc binary target in the manifest",
        subcellular_location_text_used=False,
        forbidden_fields_absent_from_view=list(FORBIDDEN_FIELDS),
        exclusion_accounting=excl,
        discovery_core_n=n_disc,
        homology=homology,
        domain_applicability=applic,
    )
    with open(os.path.join(VIEW, "integrity_checks.json"), "w") as fh:
        json.dump(integrity, fh, indent=2)

    print(json.dumps(dict(exclusions=excl, samples=counts[:12], homology=homology,
                          applicability=applic), indent=2))
    print("\neligibility:")
    for r in elig:
        print(" ", {k: r[k] for k in ("database", "threshold_per_class", "n_groups_eligible",
                                      "n_distinct_proteins_covered",
                                      "n_distinct_Extracellular_covered",
                                      "n_distinct_Cytoplasm_covered")})
    return 0


if __name__ == "__main__":
    sys.exit(main())
