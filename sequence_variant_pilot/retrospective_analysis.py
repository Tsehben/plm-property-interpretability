#!/usr/bin/env python3
"""
RETROSPECTIVE, EXPLORATORY analyses of the completed P3 sequence-variant pilot.

Saved outputs only. No model is run, no variant is generated, no experimental
setting is changed, P0 is never read. Nothing in v1/ is modified: everything
written here goes to v1_retrospective/, to keep these analyses separate from the
pre-registered ones.

NOT pre-registered. The frozen protocol specified no uncertainty quantification
and no tie decomposition for this pilot; both were decided after the results were
seen. Intervals below therefore carry no pre-registered error control and are
descriptive of this sample.

(a) Uncertainty for the four mean paired differences, by protein-level resampling
    stratified by localisation class, with the same BCa machinery and the same
    pre-registered degenerate fallback as the P4 evaluation.

(b) Tie decomposition: for every paired cell, is a zero difference the result of
    the two arms generating the SAME sequence (unchanged generation), of DIFFERENT
    sequences receiving the SAME predictor score (evaluator insensitivity), or
    neither.
"""

import glob
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
V1 = os.path.join(HERE, "v1")
OUT = os.path.join(HERE, "v1_retrospective")
sys.path.insert(0, os.path.join(ROOT, "candidate_discovery"))
sys.path.insert(0, os.path.join(ROOT, "evaluation_p4"))
import run_discovery as RD  # noqa: E402
from analyse_eval import bca  # noqa: E402

CONTRASTS = [("B", "A"), ("C", "A"), ("B", "D"), ("C", "E")]
LABEL = {"A": "A unmodified", "B": "B 1446 x2", "C": "C 1446 x0",
         "D": "D random +u", "E": "E random -u"}
SEED = 20261122
N_BOOT = 10000


def main():
    os.makedirs(OUT, exist_ok=True)
    frozen = json.load(open(os.path.join(V1, "frozen_protocol.json")))
    nrep = frozen["config"]["n_replicates"]
    scores = {r["query_id"]: r for r in RD.read_tsv(os.path.join(V1, "wolfpsort_scores.tsv"))}
    pp = RD.read_tsv(os.path.join(V1, "per_protein_contrasts.tsv"))

    per, cls = {}, {}
    for fp in sorted(glob.glob(os.path.join(V1, "per_protein", "*.json"))):
        d = json.load(open(fp))
        per[d["protein"]] = d
        cls[d["protein"]] = "Extracellular" if d["target"] == 1 else "Cytoplasm"
    accs = sorted(per)

    # ------------------------------------------------------------- (a) intervals
    rng = np.random.default_rng(SEED)
    rows = []
    for hi, lo in CONTRASTS:
        sub = [r for r in pp if r["contrast"] == f"{hi}-{lo}"]
        vals = [float(r["mean_paired_difference"]) for r in sub]
        labs = [r["localisation_class"] for r in sub]
        theta, ivs, boots = bca(vals, np.mean, rng, N_BOOT, [0.95], labs)
        iv = ivs[0]
        v = np.asarray(vals)
        rows.append(dict(
            contrast=f"{hi}-{lo}", description=f"{LABEL[hi]} minus {LABEL[lo]}",
            n_proteins=len(v),
            n_Extracellular=sum(1 for l in labs if l == "Extracellular"),
            n_Cytoplasm=sum(1 for l in labs if l == "Cytoplasm"),
            estimate_mean=f"{theta:.17g}",
            median_of_paired_differences=f"{np.median(v):.17g}",
            sd=f"{v.std(ddof=1):.17g}",
            ci95_lo=f"{iv['lo']:.17g}", ci95_hi=f"{iv['hi']:.17g}",
            ci95_method=iv["method"], ci95_excludes_zero=iv["excludes_zero"],
            ci95_fallback_reason=iv["fallback_reason"] or "",
            bootstrap_sd=f"{iv['boot_sd']:.17g}",
            n_proteins_positive=int((v > 0).sum()), n_proteins_negative=int((v < 0).sum()),
            n_proteins_exact_zero=int((v == 0).sum()),
            status="RETROSPECTIVE, EXPLORATORY: not pre-registered, no error control"))
    RD.write_tsv(os.path.join(OUT, "retrospective_intervals.tsv"),
                 list(rows[0].keys()), rows)

    # --------------------------------------------------------- (b) tie decomposition
    seq = {}
    marg = {}
    for a in accs:
        marg[(a, "orig")] = float(scores[f"{a}__orig"]["margin"])
        for v in per[a]["variants"]:
            seq[(a, v["arm"], v["replicate"])] = v["sequence"]
            marg[(a, v["arm"], v["replicate"])] = float(scores[v["variant_id"]]["margin"])

    tie_rows, cell_rows = [], []
    for hi, lo in CONTRASTS:
        c = Counter()
        for a in accs:
            for k in range(nrep):
                sh, sl = seq[(a, hi, k)], seq[(a, lo, k)]
                mh, ml = marg[(a, hi, k)], marg[(a, lo, k)]
                same_seq = sh == sl
                same_score = mh == ml
                if same_seq:
                    cat = ("tie: identical generated sequences (unchanged generation)"
                           if same_score else "IMPOSSIBLE: same sequence, different score")
                elif same_score:
                    cat = "tie: distinct sequences, identical predictor score (evaluator insensitivity)"
                else:
                    cat = "not a tie: distinct sequences, different predictor scores"
                c[cat] += 1
                cell_rows.append(dict(ACC=a, localisation_class=cls[a],
                                      contrast=f"{hi}-{lo}", replicate=k,
                                      identical_sequences=int(same_seq),
                                      identical_scores=int(same_score),
                                      hamming_between_arms=sum(
                                          1 for x, y in zip(sh, sl) if x != y),
                                      margin_hi=f"{mh:.17g}", margin_lo=f"{ml:.17g}",
                                      difference=f"{mh - ml:.17g}", category=cat))
        tot = sum(c.values())
        n_same_seq = c["tie: identical generated sequences (unchanged generation)"]
        n_eval = c["tie: distinct sequences, identical predictor score (evaluator insensitivity)"]
        n_diff = c["not a tie: distinct sequences, different predictor scores"]
        ties = n_same_seq + n_eval
        tie_rows.append(dict(
            contrast=f"{hi}-{lo}", n_paired_cells=tot,
            n_tied_cells=ties, pct_tied=f"{100*ties/tot:.1f}",
            n_tie_identical_sequences=n_same_seq,
            pct_of_ties_from_identical_sequences=f"{100*n_same_seq/ties:.1f}" if ties else "",
            n_tie_distinct_sequences_same_score=n_eval,
            pct_of_ties_from_evaluator_insensitivity=f"{100*n_eval/ties:.1f}" if ties else "",
            n_not_tied=n_diff,
            n_impossible=c["IMPOSSIBLE: same sequence, different score"]))
    RD.write_tsv(os.path.join(OUT, "tie_decomposition.tsv"), list(tie_rows[0].keys()), tie_rows)
    RD.write_tsv(os.path.join(OUT, "tie_cells.tsv"), list(cell_rows[0].keys()), cell_rows)

    # ------------------------------------------------- evaluator resolution overall
    variants = [(a, v) for a in accs for v in per[a]["variants"]]
    seqs = [v["sequence"] for _, v in variants]
    margins = [float(scores[v["variant_id"]]["margin"]) for _, v in variants]
    by_prot_seq = defaultdict(set)
    by_prot_marg = defaultdict(set)
    ident_orig = 0
    for a, v in variants:
        by_prot_seq[a].add(v["sequence"])
        by_prot_marg[a].add(float(scores[v["variant_id"]]["margin"]))
        ident_orig += int(v["identical_to_original"])
    # distinct-sequence pairs within a protein that share a margin
    coll = 0
    coll_tot = 0
    for a in accs:
        m = defaultdict(set)
        for v in per[a]["variants"]:
            m[float(scores[v["variant_id"]]["margin"])].add(v["sequence"])
        for mv, ss in m.items():
            n = len(ss)
            coll += n * (n - 1) // 2
        n_all = len(by_prot_seq[a])
        coll_tot += n_all * (n_all - 1) // 2
    resolution = dict(
        n_variants=len(variants),
        n_distinct_variant_sequences=len(set(seqs)),
        n_variants_identical_to_their_original=ident_orig,
        mean_distinct_sequences_per_protein=f"{np.mean([len(by_prot_seq[a]) for a in accs]):.2f}",
        mean_distinct_margins_per_protein=f"{np.mean([len(by_prot_marg[a]) for a in accs]):.2f}",
        n_distinct_margin_values_overall=len(set(margins)),
        distinct_sequence_pairs_within_protein=coll_tot,
        distinct_sequence_pairs_sharing_a_margin=coll,
        pct_distinct_sequence_pairs_sharing_a_margin=f"{100*coll/coll_tot:.1f}"
        if coll_tot else "",
        note="20 distinct sequences per protein at most (4 replicates x 5 arms); the "
             "predictor resolves far fewer distinct margins than there are distinct "
             "sequences, which is evaluator insensitivity rather than unchanged generation",
    )

    summary = dict(
        status="RETROSPECTIVE AND EXPLORATORY. Not pre-registered. Decided after the "
               "results were seen. No error control, no confirmatory reading.",
        resampling_assumptions=dict(
            unit="protein (n=40); replicates averaged within protein first",
            scheme="nonparametric resampling of proteins with replacement, stratified by "
                   "localisation class (20/20), 10,000 resamples, seed 20261122",
            interval="BCa with jackknife acceleration, with the same pre-registered "
                     "degenerate fallback to the percentile interval as the P4 evaluation",
            assumes="that the 40 proteins are exchangeable within class and independent. "
                    "They are not known to be: no homology clustering exists for this "
                    "cohort, and 5 of the 40 are members of the predictor's own reference "
                    "set. Unmeasured homology dependence can make these intervals too "
                    "narrow; the direction and size of that error are unknown.",
            not_units="replicates, steps and arms are not independent sample units",
            zero_policy="an interval containing zero is NOT evidence of equivalence. It "
                        "means this sample does not resolve the sign at this width. No "
                        "equivalence bound was pre-specified and none is claimed.",
        ),
        intervals=rows, tie_decomposition=tie_rows, evaluator_resolution=resolution,
        output_hashes={fn: RD.sha256_file(os.path.join(OUT, fn)) for fn in (
            "retrospective_intervals.tsv", "tie_decomposition.tsv", "tie_cells.tsv")},
    )
    with open(os.path.join(OUT, "retrospective_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(dict(intervals=rows, ties=tie_rows, resolution=resolution), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
