#!/usr/bin/env python3
"""
Bounded reconciliation of the annotation-readiness coverage report.

Reads ONLY the saved tables in out/ plus the frozen manifest. Makes no network
requests, runs no models, computes no clusters, and modifies nothing in out/.
All corrected outputs are written to out_reconciled/ so the originals are
preserved for comparison.

Resolves five accounting issues raised in review:
  1. "Other"-taxonomy categories must be mutually exclusive and sum to 342.
  2. Status dimensions must be reported separately and exhaustively, not summed
     across dimensions into a single "failures" count.
  3. Signal-peptide coverage needs exact denominators, an explicit statement of
     whether counts are proteins or feature records, and an explicit treatment
     of the exact-match isoforms whose canonical coordinates do not transfer.
  4. P0-adjacent retrievals must be counted at the level of unique development
     proteins and separated from the (true) statement that no P0 row is an
     annotation target.
  5. Evidence categories must be reported as overlapping sets, with evidence
     codes distinguished from identified prediction tools.

Usage:  python3 reconcile_coverage.py
"""

import csv
import json
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
REC = os.path.join(HERE, "out_reconciled")
MANIFEST = os.path.abspath(
    os.path.join(HERE, "..", "DeepLoc2", "baseline", "deeploc_binary_manifest.tsv")
)

CLASS_NAME = {"1": "Extracellular", "0": "Cytoplasm"}
ROLES = ["discovery", "tuning", "evaluation"]
# WoLF PSORT's three organism types (Horton et al., NAR 35:W585-W587, 2007).
WP_TYPES = ["animal", "fungi", "plant"]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, fieldnames, records):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in records:
            w.writerow(r)


def pct(n, d):
    return "" if not d else f"{100.0 * n / d:.1f}"


# ---------------------------------------------------------------- load

def load():
    primary = read_tsv(os.path.join(OUT, "annotation_primary.tsv"))
    features = read_tsv(os.path.join(OUT, "annotation_features.tsv"))
    domains = read_tsv(os.path.join(OUT, "annotation_domains.tsv"))
    manifest = {r["ACC"]: r for r in read_tsv(MANIFEST)}
    return primary, features, domains, manifest


# ------------------------------------------------- 1. organism categories

def organism_categories(primary):
    """Mutually exclusive accession-level organism categories over all 10,139
    rows. The review found 221+98+22+1+1 = 343 reported against 342 proteins.
    Cause: the single 'unresolved' record is NOT one of the 342 DeepLoc-'Other'
    proteins. It is a Metazoa-kingdom protein whose UniProt record could not be
    retrieved at all. The two facts belong to different categories and must not
    be listed in the same series."""
    recs = []
    # A protein is 'resolved' iff a UniProt lineage was retrieved for it.
    resolved = [r for r in primary if r["organism_lineage"].strip()]
    unresolved = [r for r in primary if not r["organism_lineage"].strip()]

    for wp in WP_TYPES:
        sub = [r for r in resolved if r["wolfpsort_organism_type"] == wp]
        recs.append(
            dict(
                category=f"resolved / supported by WoLF PSORT ({wp})",
                basis="UniProt lineage retrieved; lineage falls in a supported group",
                n=len(sub),
                example=sub[0]["organism_scientific_name"] if sub else "",
            )
        )

    # Resolved but outside the three supported groups, split by lineage rank 2.
    outside = [r for r in resolved if r["wolfpsort_organism_type"] == "unsupported"]
    by_clade = defaultdict(list)
    for r in outside:
        lin = [x.strip() for x in r["organism_lineage"].split(";") if x.strip()]
        by_clade[lin[1] if len(lin) > 1 else "(no rank below Eukaryota)"].append(r)
    for clade in sorted(by_clade, key=lambda c: -len(by_clade[c])):
        sub = by_clade[clade]
        recs.append(
            dict(
                category=f"resolved / outside WoLF PSORT's supported groups ({clade})",
                basis="UniProt lineage retrieved; lineage is not Metazoa, Viridiplantae or Fungi",
                n=len(sub),
                example=sub[0]["organism_scientific_name"],
            )
        )

    recs.append(
        dict(
            category="taxonomy unresolved (no UniProt record retrieved)",
            basis="no lineage available; organism support cannot be assessed either way",
            n=len(unresolved),
            example=";".join(r["ACC"] for r in unresolved),
        )
    )

    total = sum(r["n"] for r in recs)
    assert total == len(primary), (total, len(primary))

    # The DeepLoc-'Other' kingdom bucket, checked against the categories above.
    other = [r for r in primary if r["kingdom_deeploc"] == "Other"]
    other_resolved_outside = [r for r in other if r in outside]
    cross = dict(
        n_deeploc_kingdom_other=len(other),
        n_of_those_resolved=sum(1 for r in other if r["organism_lineage"].strip()),
        n_of_those_outside_supported_groups=len(other_resolved_outside),
        n_of_those_unresolved=sum(1 for r in other if not r["organism_lineage"].strip()),
        n_unresolved_outside_the_other_bucket=sum(
            1 for r in unresolved if r["kingdom_deeploc"] != "Other"
        ),
        unresolved_accessions=[r["ACC"] for r in unresolved],
        unresolved_deeploc_kingdom=[r["kingdom_deeploc"] for r in unresolved],
    )
    return recs, cross


# --------------------------------------------- 2. exhaustive status dimensions

def status_dimensions(primary):
    """Four independent dimensions. Each partitions all 10,139 proteins. They
    are NOT added together: a protein can be non-exact in dimension B and
    non-transferable in dimension C for the same underlying reason."""
    n = len(primary)
    dims = []

    # A. Retrieval / accession status (refined: the original 'ok' bucket hid two
    #    proteins whose isoform FASTA could not be retrieved at all).
    def acc_status(r):
        if r["retrieval_status"] == "inactive_deleted":
            return "inactive_deleted (no record)"
        if r["retrieval_status"] == "resolved_via_merged_accession":
            return "resolved via merged/secondary accession"
        if r["sequence_source"] == "none":
            return "active entry, isoform sequence unresolved (HTTP 404)"
        return "resolved as requested (current accession)"

    dims.append(("A. Retrieval / accession status", Counter(acc_status(r) for r in primary)))

    # B. Frozen-sequence agreement.
    dims.append(
        (
            "B. Agreement with frozen extraction-input sequence",
            Counter(
                {
                    "exact (SHA-256 identical)": sum(1 for r in primary if r["sequence_match"] == "exact"),
                    "mismatch (entry re-sequenced since the DeepLoc snapshot)": sum(
                        1 for r in primary if r["sequence_match"] == "mismatch"
                    ),
                    "unverifiable (no sequence retrieved)": sum(
                        1 for r in primary if r["sequence_match"] == "unverifiable"
                    ),
                }
            ),
        )
    )

    # C. Coordinate transferability, with the reason split out.
    def coord_status(r):
        if r["coordinates_transferable"] == "1":
            return "transferable (frozen sequence == the sequence coordinates are defined on)"
        if r["sequence_match"] == "exact" and r["is_isoform"] == "1":
            return "not transferable: isoform matches frozen input but differs from canonical"
        if r["sequence_match"] == "mismatch":
            return "not transferable: retrieved sequence differs from frozen input"
        return "not transferable: no sequence retrieved"

    dims.append(("C. Residue-coordinate transferability", Counter(coord_status(r) for r in primary)))

    # D. Final signal-peptide status.
    dims.append(
        (
            "D. Signal-peptide status",
            Counter(
                {
                    "annotated_present": sum(1 for r in primary if r["signal_status"] == "annotated_present"),
                    "not_annotated (NOT a confirmed biological negative)": sum(
                        1 for r in primary if r["signal_status"] == "not_annotated"
                    ),
                    "unverifiable": sum(1 for r in primary if r["signal_status"] == "unverifiable"),
                }
            ),
        )
    )

    for name, c in dims:
        assert sum(c.values()) == n, (name, sum(c.values()), n)

    # Overlap between the conditions the original report summed into "122".
    non_current = [r for r in primary if r["retrieval_status"] != "ok"]
    non_exact = [r for r in primary if r["sequence_match"] != "exact"]
    both = [r for r in non_current if r["sequence_match"] != "exact"]
    merged_exact = [
        r
        for r in primary
        if r["retrieval_status"] == "resolved_via_merged_accession" and r["sequence_match"] == "exact"
    ]
    overlap = dict(
        n_accession_not_current=len(non_current),
        n_sequence_not_exact=len(non_exact),
        n_in_both=len(both),
        n_union=len(set(r["ACC"] for r in non_current) | set(r["ACC"] for r in non_exact)),
        naive_sum=len(non_current) + len(non_exact),
        merged_with_exact_sequence=[r["ACC"] for r in merged_exact],
        note=(
            "The original report's 122 is the UNION of 'accession not current' and "
            "'sequence not exact', not a sum of independent failures. Of that union, "
            f"{len(merged_exact)} rows are merged accessions whose retrieved sequence is "
            "SHA-256 identical to the frozen input: those are provenance caveats, not "
            "annotation failures."
        ),
    )
    return dims, overlap


def reclassified_non_exact(primary):
    """Replace the single 'failures' bucket with graded outcomes."""
    out = []
    for r in primary:
        if r["retrieval_status"] == "ok" and r["sequence_match"] == "exact":
            continue
        if r["sequence_match"] == "exact" and r["retrieval_status"] == "resolved_via_merged_accession":
            cls = "usable_with_provenance_caveat"
            why = "accession merged into a successor entry; retrieved sequence identical to frozen input"
        elif r["sequence_match"] == "mismatch":
            cls = "annotation_not_transferable"
            why = "retrieved sequence differs from frozen input; no residue coordinates transferred"
        else:
            cls = "no_annotation_available"
            why = r["signal_status_reason"] or "no sequence retrieved"
        out.append(
            dict(
                ACC=r["ACC"],
                role=r["role"],
                localisation_class=CLASS_NAME[r["target"]],
                retrieval_status=r["retrieval_status"],
                merged_into=r["merged_into"],
                sequence_match=r["sequence_match"],
                signal_status=r["signal_status"],
                outcome_class=cls,
                outcome_reason=why,
            )
        )
    return out


# ------------------------------------- 3. signal coverage with denominators

def signal_coverage(primary):
    """Protein counts (not feature records) by role x class, with denominators
    and all three states shown so nothing is hidden in a percentage."""
    recs = []
    groups = [(role, t) for role in ROLES for t in ("1", "0")] + [("ALL", "1"), ("ALL", "0"), ("ALL", "ALL")]
    for role, t in groups:
        sub = [
            r
            for r in primary
            if (role == "ALL" or r["role"] == role) and (t == "ALL" or r["target"] == t)
        ]
        d = len(sub)
        present = sum(1 for r in sub if r["signal_status"] == "annotated_present")
        nota = sum(1 for r in sub if r["signal_status"] == "not_annotated")
        unver = sum(1 for r in sub if r["signal_status"] == "unverifiable")
        exp = sum(
            1
            for r in sub
            if r["signal_status"] == "annotated_present" and r["signal_evidence_best"] == "experimental"
        )
        recs.append(
            dict(
                role=role,
                localisation_class="ALL" if t == "ALL" else CLASS_NAME[t],
                n_proteins=d,
                annotated_present=present,
                pct_annotated_present=pct(present, d),
                not_annotated=nota,
                pct_not_annotated=pct(nota, d),
                unverifiable=unver,
                pct_unverifiable=pct(unver, d),
                annotated_present_experimental=exp,
                pct_of_cohort_experimental=pct(exp, d),
            )
        )
        assert present + nota + unver == d
    return recs


def unverifiable_breakdown(primary):
    """Every unverifiable protein, by the reason it is unverifiable, by role and
    class. These are explicitly NOT verified absence of a signal peptide."""
    recs = []
    kinds = [
        (
            "isoform matches frozen input, differs from canonical",
            lambda r: r["sequence_match"] == "exact" and r["is_isoform"] == "1",
        ),
        ("retrieved sequence differs from frozen input", lambda r: r["sequence_match"] == "mismatch"),
        ("no sequence retrieved", lambda r: r["sequence_match"] == "unverifiable"),
    ]
    unv = [r for r in primary if r["signal_status"] == "unverifiable"]
    for name, test in kinds:
        sub = [r for r in unv if test(r)]
        for role in ROLES:
            for t in ("1", "0"):
                s2 = [r for r in sub if r["role"] == role and r["target"] == t]
                if not s2:
                    continue
                recs.append(
                    dict(
                        reason=name,
                        role=role,
                        localisation_class=CLASS_NAME[t],
                        n_proteins=len(s2),
                        parent_entry_has_signal_annotation=sum(
                            1 for r in s2 if r["parent_entry_has_signal"] == "1"
                        ),
                    )
                )
    assert sum(r["n_proteins"] for r in recs) == len(unv)
    return recs, len(unv)


def evidence_overlap(primary, features):
    """Evidence categories are OVERLAPPING sets at the protein level: one signal
    feature can carry several evidences. Report both the overlapping membership
    and the mutually exclusive strongest-evidence collapse, and keep them
    labelled as different things."""
    sig = [r for r in features if r["feature_type"] == "Signal"]
    present = set(r["ACC"] for r in primary if r["signal_status"] == "annotated_present")
    sig_present = [r for r in sig if r["ACC"] in present]

    members = defaultdict(set)
    for r in sig_present:
        members[r["evidence_category"]].add(r["ACC"])
    n_present = len(present)

    recs = []
    for cat in sorted(members, key=lambda c: -len(members[c])):
        recs.append(
            dict(
                evidence_category=cat,
                counting="overlapping set membership (proteins)",
                n_proteins=len(members[cat]),
                pct_of_annotated=pct(len(members[cat]), n_present),
            )
        )
    strongest = Counter(
        r["signal_evidence_best"] for r in primary if r["signal_status"] == "annotated_present"
    )
    for cat, k in strongest.most_common():
        recs.append(
            dict(
                evidence_category=cat,
                counting="mutually exclusive (strongest evidence per protein)",
                n_proteins=k,
                pct_of_annotated=pct(k, n_present),
            )
        )
    multi = sum(1 for a in present if len(set(r["evidence_category"] for r in sig_present if r["ACC"] == a)) > 1)

    # Is any prediction tool actually identified for ECO:0000255?
    eco255 = [r for r in sig if r["evidence_code"] == "ECO:0000255"]
    named = [r for r in eco255 if (r["evidence_source"] or "").strip()]
    tool = dict(
        n_eco255_evidence_rows=len(eco255),
        n_with_a_named_source=len(named),
        named_sources={f"{k[0]}:{k[1]}": v for k, v in Counter((r["evidence_source"], r["evidence_id"]) for r in named).items()},
        note=(
            "ECO:0000255 is 'match to sequence model evidence used in manual assertion'. "
            "It does not name a tool. In this cohort the evidence source field is empty for "
            f"{len(eco255) - len(named)} of {len(eco255)} ECO:0000255 signal-peptide evidence "
            "rows, so the generating method is unidentified and must not be reported as any "
            "particular predictor."
        ),
    )
    counting = dict(
        n_proteins_with_signal_annotation=n_present,
        n_signal_feature_records=sum(int(r["n_signal_features"] or 0) for r in primary),
        n_signal_feature_records_on_annotated_proteins=len(
            set((r["ACC"], r["start"], r["end"]) for r in sig_present)
        ),
        n_signal_evidence_rows=len(sig),
        n_proteins_with_more_than_one_evidence_category=multi,
        note=(
            "2,439 counts PROTEINS whose signal status is annotated_present. Each has "
            "exactly one signal feature. annotation_features.tsv holds one row per "
            "(feature, evidence) pair, so its row count is larger."
        ),
    )
    return recs, tool, counting


# ------------------------------------------------- 4. P0-adjacent retrievals

def p0_accounting(primary, manifest):
    p0 = set(a for a, r in manifest.items() if r["partition"] == "0")
    flagged = [r for r in primary if r["p0_base_collision"] != "none"]
    groups = {}
    for flag in ("base_is_p0_row", "base_shares_with_p0_isoform"):
        sub = [r for r in flagged if r["p0_base_collision"] == flag]
        groups[flag] = dict(
            n_development_proteins=len(sub),
            development_proteins=sorted(r["ACC"] for r in sub),
            n_unique_base_accessions=len(set(r["base_accession"] for r in sub)),
            unique_base_accessions=sorted(set(r["base_accession"] for r in sub)),
        )
    a = set(r["ACC"] for r in flagged if r["p0_base_collision"] == "base_is_p0_row")
    b = set(r["ACC"] for r in flagged if r["p0_base_collision"] == "base_shares_with_p0_isoform")

    # Were records that are themselves P0 cohort rows actually fetched?
    fetched_p0 = sorted(
        set(r["base_accession"] for r in flagged if r["base_accession"] in p0)
    )
    # Which entry-level fields on the development rows come from that parent record?
    parent_derived = [
        "organism_scientific_name",
        "organism_taxon_id",
        "organism_lineage",
        "entry_version",
        "first_public_date",
        "last_annotation_update",
        "last_sequence_update",
        "n_pfam",
        "n_interpro",
        "pfam_ids",
        "interpro_ids",
        "family_similarity_text",
        "subcellular_location_text",
        "parent_entry_has_signal",
    ]
    affected = [r for r in flagged if r["base_accession"] in p0]
    loc_used = [r for r in affected if r["subcellular_location_text"].strip()]

    return dict(
        n_unique_development_proteins_affected=len(a | b),
        n_in_both_groups=len(a & b),
        groups=groups,
        p0_cohort_rows_with_an_annotation_row=0,
        p0_cohort_rows_requested_as_annotation_targets=0,
        n_p0_accessions_fetched_as_parent_records=len(fetched_p0),
        p0_accessions_fetched_as_parent_records=fetched_p0,
        n_development_rows_carrying_parent_derived_fields_from_a_p0_accession=len(affected),
        parent_derived_fields=parent_derived,
        p0_localisation_information_used=dict(
            deeploc_p0_labels_used=False,
            uniprot_subcellular_location_comment_copied_from_a_p0_accession_entry=len(loc_used),
            affected_development_rows=sorted(r["ACC"] for r in loc_used),
            note=(
                "The DeepLoc P0 label was never read. However, the UniProt entry for these "
                "parent accessions IS a P0 cohort row, and its entry-level comments - "
                "including SUBCELLULAR LOCATION - were copied onto the development isoform "
                "rows. This is entry-level UniProt annotation about a held-out protein "
                "reaching a development row, and it is reported rather than claimed absent. "
                "No additional P0 records were inspected for this reconciliation."
            ),
        ),
        flags_retained=True,
        rows_deleted=0,
        holdout_altered=False,
    )


# ---------------------------------------------------- family/domain coverage

def family_coverage(primary, domains):
    recs = []
    for role in ROLES + ["ALL"]:
        sub = [r for r in primary if role == "ALL" or r["role"] == role]
        d = len(sub)
        pf = sum(1 for r in sub if int(r["n_pfam"] or 0) > 0)
        ip = sum(1 for r in sub if int(r["n_interpro"] or 0) > 0)
        sim = sum(1 for r in sub if r["family_similarity_text"].strip())
        recs.append(
            dict(
                role=role,
                n_proteins=d,
                with_pfam=pf,
                pct_with_pfam=pct(pf, d),
                with_interpro=ip,
                pct_with_interpro=pct(ip, d),
                with_similarity_comment=sim,
                pct_with_similarity_comment=pct(sim, d),
            )
        )
    stats = dict(
        n_pfam_xref_rows=sum(1 for r in domains if r["database"] == "Pfam"),
        n_interpro_xref_rows=sum(1 for r in domains if r["database"] == "InterPro"),
        n_distinct_pfam=len(set(r.get("xref_id", "") for r in domains if r.get("database") == "Pfam")),
        n_distinct_interpro=len(set(r.get("xref_id", "") for r in domains if r.get("database") == "InterPro")),
        caveat=(
            "A Pfam or InterPro entry is a region-level annotation. It is not a "
            "whole-protein family assignment, and multi-domain proteins carry several. "
            "No homology clustering has been computed, and no family-disjointness between "
            "partitions is claimed."
        ),
    )
    return recs, stats


# ---------------------------------------------------------------- main

def main():
    os.makedirs(REC, exist_ok=True)
    primary, features, domains, manifest = load()
    n = len(primary)
    assert n == 10139, n

    org, org_cross = organism_categories(primary)
    dims, overlap = status_dimensions(primary)
    reclass = reclassified_non_exact(primary)
    cov = signal_coverage(primary)
    unv, n_unv = unverifiable_breakdown(primary)
    evid, tool, counting = evidence_overlap(primary, features)
    p0 = p0_accounting(primary, manifest)
    fam, famstats = family_coverage(primary, domains)

    write_tsv(os.path.join(REC, "organism_categories.tsv"), ["category", "basis", "n", "example"], org)
    write_tsv(
        os.path.join(REC, "status_dimensions.tsv"),
        ["dimension", "category", "n_proteins", "pct_of_10139"],
        [
            dict(dimension=name, category=k, n_proteins=v, pct_of_10139=pct(v, n))
            for name, c in dims
            for k, v in c.most_common()
        ],
    )
    write_tsv(
        os.path.join(REC, "non_exact_records_reclassified.tsv"),
        list(reclass[0].keys()),
        reclass,
    )
    write_tsv(os.path.join(REC, "reconciled_coverage.tsv"), list(cov[0].keys()), cov)
    write_tsv(os.path.join(REC, "unverifiable_breakdown.tsv"), list(unv[0].keys()), unv)
    write_tsv(os.path.join(REC, "evidence_categories.tsv"), list(evid[0].keys()), evid)
    write_tsv(os.path.join(REC, "family_domain_coverage.tsv"), list(fam[0].keys()), fam)

    summary = dict(
        inputs_read=[
            "out/annotation_primary.tsv",
            "out/annotation_features.tsv",
            "out/annotation_domains.tsv",
            "DeepLoc2/baseline/deeploc_binary_manifest.tsv",
        ],
        network_requests=0,
        originals_modified=[],
        n_primary_rows=n,
        organism_cross_check=org_cross,
        status_dimension_overlap=overlap,
        reclassified_outcomes=dict(Counter(r["outcome_class"] for r in reclass)),
        signal_counting=counting,
        n_unverifiable=n_unv,
        eco255_tool_attribution=tool,
        p0_accounting=p0,
        family_domain=famstats,
    )
    with open(os.path.join(REC, "reconciliation_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
