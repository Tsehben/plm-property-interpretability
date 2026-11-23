#!/usr/bin/env python3
"""
Audit of the reference/training data shipped inside the WoLF PSORT 0.2 package.

Package-side only. This builds an index of what the distribution contains and
what identifier/sequence keys it exposes. It does NOT read the project cohort,
does not touch P0, and performs no overlap comparison; the point is to establish
exactly which overlap checks are technically possible later.
"""

import hashlib
import json
import os
import re
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "src_download", "wolfpsort_repo")
OUT = os.path.join(HERE, "reference_audit")
ORGS = ("animal", "fungi", "plant")


def read_psort_fasta(path):
    recs, ident, cls, seq = [], None, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if ident is not None:
                    recs.append((ident, cls, "".join(seq)))
                p = line[1:].split()
                ident, cls, seq = (p[0] if p else ""), (p[1] if len(p) > 1 else ""), []
            else:
                seq.append(re.sub(r"\s", "", line))
    if ident is not None:
        recs.append((ident, cls, "".join(seq)))
    return recs


def read_all_seq(path):
    recs, ident, cls, seq = [], None, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if ident is not None:
                    recs.append((ident, cls, "".join(seq)))
                p = line[1:].split("\t")
                head = p[0].split()
                ident = head[0]
                cls = p[1].strip() if len(p) > 1 else (head[1] if len(head) > 1 else "")
                seq = []
            else:
                seq.append(re.sub(r"\s", "", line))
    if ident is not None:
        recs.append((ident, cls, "".join(seq)))
    return recs


def main():
    os.makedirs(OUT, exist_ok=True)
    report = dict(package="WoLF PSORT Command Line 0.2", contents={}, keys={}, feasibility=[])

    # ------------------------------------------------- kNN reference sequences
    seq_index, per_org = {}, {}
    for org in ORGS:
        recs = read_psort_fasta(os.path.join(REPO, "data", f"{org}.psort"))
        cls = Counter(c for _, c, _ in recs)
        lens = sorted(len(s) for _, _, s in recs)
        per_org[org] = dict(
            file=f"data/{org}.psort", kind="FASTA: identifier + class label + full sequence",
            n_records=len(recs), n_distinct_identifiers=len({i for i, _, _ in recs}),
            n_distinct_sequences=len({s for _, _, s in recs}),
            n_classes=len(cls), n_single_site_classes=sum(1 for c in cls if "_" not in c),
            n_dual_site_classes=sum(1 for c in cls if "_" in c),
            n_extr=cls.get("extr", 0), n_cyto=cls.get("cyto", 0),
            median_length=lens[len(lens) // 2], min_length=lens[0], max_length=lens[-1],
            identifier_example=recs[0][0],
            identifier_namespace="UniProt entry name (ENTRY_SPECIES), not an accession")
        for i, c, s in recs:
            seq_index.setdefault(hashlib.sha256(s.encode()).hexdigest(),
                                 []).append(dict(organism=org, id=i, cls=c, length=len(s)))
    report["contents"]["knn_reference_sequences"] = per_org

    # ---------------------------------------------------- kNN feature vectors
    feats = {}
    for org in ORGS:
        p = os.path.join(REPO, "data", f"{org}.wolff")
        with open(p) as fh:
            head = fh.readline().split()
            rows = [ln.split() for ln in fh if ln.strip()]
        ids = {r[0] for r in rows}
        psort_ids = {i for i, _, _ in read_psort_fasta(os.path.join(REPO, "data", f"{org}.psort"))}
        feats[org] = dict(
            file=f"data/{org}.wolff",
            kind="transformed feature matrix: identifier + class + numeric features, NO sequence",
            n_rows=len(rows), n_columns=len(head),
            n_feature_columns=len(head) - 2, feature_names=head[2:],
            n_identifiers_with_a_sequence_record=len(ids & psort_ids),
            n_identifiers_without_a_sequence_record=len(ids - psort_ids))
    report["contents"]["knn_feature_vectors"] = feats

    # ---------------------------------------------------- class similarity / weights
    sim = {}
    for org in ORGS:
        with open(os.path.join(REPO, "data", f"{org}.wolfu")) as fh:
            classes = fh.readline().split()
            rows = [ln.split() for ln in fh if ln.strip()]
        with open(os.path.join(REPO, "data", f"{org}.wolfw")) as fh:
            txt = fh.read()
        k = re.search(r"k=(\d+)", txt)
        sim[org] = dict(similarity_file=f"data/{org}.wolfu", n_classes=len(classes),
                        classes=classes, matrix_rows=len(rows),
                        weight_file=f"data/{org}.wolfw",
                        k_for_knn=int(k.group(1)) if k else None)
    report["contents"]["class_similarity_and_feature_weights"] = sim

    # --------------------------------------------- classical PSORT training set
    al = read_all_seq(os.path.join(REPO, "bin", "psortModifiedForWolfFiles", "all.seq"))
    acc = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")
    report["contents"]["classical_psort_training_set"] = dict(
        file="bin/psortModifiedForWolfFiles/all.seq",
        kind="FASTA: identifier + 3-letter PSORT site + description + full sequence",
        n_records=len(al),
        n_identifiers_matching_uniprot_accession_syntax=sum(1 for i, _, _ in al if acc.match(i)),
        identifier_example=al[0][0],
        identifier_namespace="UniProt accession",
        classes=sorted({c for _, c, _ in al}),
        used_by="the PSORT-classical dbm tables (runWolfPsortHtmlTables); "
                "runWolfPsortSummary does not consult it at prediction time")

    # ------------------------------------------------------------ key summary
    dup = {h: v for h, v in seq_index.items() if len(v) > 1}
    report["keys"] = dict(
        union_distinct_reference_identifiers=len(
            {(r["id"]) for v in seq_index.values() for r in v}),
        union_distinct_reference_sequences=len(seq_index),
        n_sequences_shared_across_organism_models=sum(
            1 for v in dup.values() if len({r["organism"] for r in v}) > 1),
        n_sequences_duplicated_within_a_model=sum(
            1 for v in dup.values() if len({r["organism"] for r in v}) == 1),
        reference_sequence_sha256_index="reference_audit/reference_sequence_sha256.tsv")

    with open(os.path.join(OUT, "reference_sequence_sha256.tsv"), "w") as fh:
        fh.write("sha256\torganism\tidentifier\tclass\tlength\n")
        for h in sorted(seq_index):
            for r in seq_index[h]:
                fh.write(f"{h}\t{r['organism']}\t{r['id']}\t{r['cls']}\t{r['length']}\n")

    report["feasibility"] = [
        dict(check="exact sequence overlap (cohort frozen sequence vs reference sequence)",
             possible=True,
             how="SHA-256 of the whitespace-stripped, upper-cased residue string on both "
                 "sides; the package supplies full sequences for all three models",
             coverage=f"{len(seq_index)} distinct reference sequences"),
        dict(check="near-identity / homology overlap (e.g. >=90% identity)",
             possible=True,
             how="requires an aligner (BLAST/DIAMOND/MMseqs2), none of which is in this "
                 "package or this environment; the sequences needed to run one ARE present",
             coverage="not installed here; would be a separate, declared step"),
        dict(check="identifier overlap by UniProt accession",
             possible=False,
             how="the kNN reference files are keyed by UniProt ENTRY NAME, not accession; "
                 "mapping entry names to accessions needs the UniProt ID-mapping service, "
                 "i.e. a network query. Only bin/.../all.seq is accession-keyed, and that "
                 "file feeds the classical-PSORT tables, not the kNN prediction.",
             coverage="blocked without a network lookup"),
        dict(check="overlap against the kNN FEATURE matrix as such",
             possible=False,
             how="data/*.wolff contains only transformed features; it carries no sequence. "
                 "Overlap can be established only through the matching identifier in "
                 "data/*.psort, never from the feature vector itself.",
             coverage="not a sequence-level check"),
        dict(check="label agreement between the reference set and the cohort",
             possible=True,
             how="only for proteins matched by exact sequence first; the class vocabularies "
                 "differ (WoLF 'extr'/'cyto' + dual classes vs the frozen DeepLoc binary "
                 "labels) and would need an explicit, pre-registered mapping",
             coverage="conditional on the exact-sequence match above"),
    ]
    report["unresolved_limitations"] = [
        "No provenance or date is recorded in the package for how the reference sets were "
        "assembled, so their relationship to the DeepLoc 2.0 cohort cannot be established "
        "from the distribution alone.",
        "The reference identifiers are UniProt entry names, which are not stable keys: an "
        "entry name can be reassigned between releases, so even a successful name-based "
        "match would not prove record identity.",
        "A sequence that is absent from the reference set may still have close homologues in "
        "it. Exact-hash overlap is a lower bound on dependence, not a measure of independence.",
        "data/animal.psort has 9 duplicated identifiers and data/animal.wolff carries 5 "
        "identifiers with no sequence record, so the two files are not in exact bijection "
        "for the animal model.",
    ]

    json.dump(report, open(os.path.join(OUT, "reference_audit.json"), "w"), indent=2)
    print(json.dumps({k: report[k] for k in ("keys", "feasibility")}, indent=2))
    for org in ORGS:
        print(org, per_org[org]["n_records"], "seqs |", feats[org]["n_rows"], "feature rows |",
              "k =", sim[org]["k_for_knn"], "|", sim[org]["n_classes"], "classes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
