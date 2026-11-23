#!/usr/bin/env python3
"""
Addendum to the frozen P4 evaluation. Saved outputs only: reads
v1/per_protein/*.json and v1/primary_results.tsv; writes v1_addendum/.
No model, no new sampling, no P0, nothing in v1/ is modified.

Purpose: report the MEAN CANDIDATE KL and the MEAN MATCHED-CONTROL KL on the
EXACT paired observations behind each primary contrast, using the identical
aggregation as the reported difference:

  1. (protein, position) cell -> average the AVAILABLE matched controls' KL;
     keep the candidate's KL for the same cell.
  2. average candidate KLs, control means, and their differences within protein
     over that protein's paired positions.
  3. average across proteins.

Because the mean is linear and the same cells enter all three series,
mean_candidate_KL - mean_control_KL reproduces the reported paired difference
exactly; this is asserted, not assumed.

These candidate means are NOT the all-sampled candidate means. The all-sampled
series is a different population (every active candidate cell, including those
with no feasible control) and the two must not be compared.
"""

import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
V1 = os.path.join(HERE, "v1")
OUT = os.path.join(HERE, "v1_addendum")
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "candidate_discovery"))
import run_discovery as RD  # noqa: E402

CANDS = [2564, 3900]
MULTS = [0.0, 2.0]
KL = "kl_aa20_renormalised_nats"


def f(v):
    return float(v)


def read_tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, ln.rstrip("\n").split("\t"))) for ln in fh]


def main():
    os.makedirs(OUT, exist_ok=True)
    meta, cls = {}, {}
    for fp in sorted(glob.glob(os.path.join(V1, "per_protein", "*.json"))):
        d = json.load(open(fp))
        meta[d["protein"]] = d
        cls[d["protein"]] = "Extracellular" if d["target"] == 1 else "Cytoplasm"
    assert len(meta) == 400, len(meta)
    reported = {(int(r["candidate"]), float(r["multiplier"])): r
                for r in read_tsv(os.path.join(V1, "primary_results.tsv"))}

    levels, perprot, cells = [], [], []
    for c in CANDS:
        for m in MULTS:
            pc, pk, pd_, labs = [], [], [], []
            pos_count = defaultdict(int)
            pos_count_cells = defaultdict(int)
            npos_total = 0
            for a in sorted(meta):
                cvs, kvs, dvs = [], [], []
                for r in meta[a]["rows"]:
                    if (r["arm"] != "candidate" or int(r["candidate"]) != c
                            or float(r["multiplier"]) != m or r["status"] != "edited"):
                        continue
                    cs = [q for q in meta[a]["rows"]
                          if q["arm"] == "control" and int(q["candidate"]) == c
                          and float(q["multiplier"]) == m
                          and q["residue_position"] == r["residue_position"]
                          and q["status"] == "edited"]
                    if not cs:
                        continue        # MISSING, never imputed as zero
                    cv = f(r[KL])
                    kv = float(np.mean([f(q[KL]) for q in cs]))
                    cvs.append(cv)
                    kvs.append(kv)
                    dvs.append(cv - kv)
                    pos_count_cells[r["residue_position"]] += 1
                    cells.append(dict(
                        ACC=a, localisation_class=cls[a], candidate=c, multiplier=m,
                        residue_position=r["residue_position"],
                        n_controls_available=len(cs),
                        candidate_kl_nats=f"{cv:.17g}",
                        mean_control_kl_nats=f"{kv:.17g}",
                        difference_nats=f"{cv - kv:.17g}"))
                if dvs:
                    pc.append(float(np.mean(cvs)))
                    pk.append(float(np.mean(kvs)))
                    pd_.append(float(np.mean(dvs)))
                    labs.append(cls[a])
                    npos_total += len(dvs)
                    pos_count[a] = len(dvs)
            pc, pk, pd_ = np.array(pc), np.array(pk), np.array(pd_)
            rep = reported[(c, m)]
            recon = float(pc.mean() - pk.mean())
            assert len(pc) == int(rep["n_proteins"]), (len(pc), rep["n_proteins"])
            assert npos_total == int(rep["n_paired_positions"])
            assert abs(float(pd_.mean()) - float(rep["estimate"])) < 1e-18
            assert abs(recon - float(rep["estimate"])) < 1e-18
            levels.append(dict(
                candidate=c, multiplier=m,
                condition="ablation (x0)" if m == 0.0 else "amplification (x2)",
                population="paired observations only (candidate cells with >=1 feasible "
                           "matched control); NOT the all-sampled candidate population",
                n_proteins=len(pc), n_paired_positions=npos_total,
                n_Extracellular=sum(1 for l in labs if l == "Extracellular"),
                n_Cytoplasm=sum(1 for l in labs if l == "Cytoplasm"),
                mean_candidate_kl_nats=f"{pc.mean():.17g}",
                median_candidate_kl_nats=f"{np.median(pc):.17g}",
                mean_control_kl_nats=f"{pk.mean():.17g}",
                median_control_kl_nats=f"{np.median(pk):.17g}",
                mean_paired_difference_nats=f"{pd_.mean():.17g}",
                reported_estimate_nats=rep["estimate"],
                difference_reconstructs_reported_estimate=bool(
                    abs(recon - float(rep["estimate"])) < 1e-18),
                ratio_candidate_over_control=f"{pc.mean()/pk.mean():.6f}"
                if pk.mean() else "",
            ))

    # ------------------------------------- class and position coverage, all contrasts
    cov = []
    for c in CANDS:
        for m in MULTS:
            sub = [x for x in cells if x["candidate"] == c and x["multiplier"] == m]
            prots = {x["ACC"] for x in sub}
            # denominator: candidate-active cells (paired or not)
            active = 0
            for a in meta:
                for r in meta[a]["rows"]:
                    if (r["arm"] == "candidate" and int(r["candidate"]) == c
                            and float(r["multiplier"]) == m and r["status"] == "edited"):
                        active += 1
            row = dict(candidate=c, multiplier=m,
                       n_candidate_active_cells=active, n_paired_cells=len(sub),
                       pct_of_active_cells_paired=f"{100.0*len(sub)/max(active,1):.1f}",
                       n_proteins_paired=len(prots),
                       n_proteins_paired_Extracellular=len(
                           {x["ACC"] for x in sub if x["localisation_class"] == "Extracellular"}),
                       n_proteins_paired_Cytoplasm=len(
                           {x["ACC"] for x in sub if x["localisation_class"] == "Cytoplasm"}))
            for p in (10, 20, 30):
                row[f"n_paired_cells_pos{p}"] = sum(1 for x in sub if x["residue_position"] == p)
                row[f"n_paired_cells_pos{p}_Extracellular"] = sum(
                    1 for x in sub if x["residue_position"] == p
                    and x["localisation_class"] == "Extracellular")
            cov.append(row)

    RD.write_tsv(os.path.join(OUT, "paired_levels.tsv"), list(levels[0].keys()), levels)
    RD.write_tsv(os.path.join(OUT, "paired_coverage_detail.tsv"), list(cov[0].keys()), cov)
    RD.write_tsv(os.path.join(OUT, "paired_cells.tsv"), list(cells[0].keys()), cells)
    summary = dict(
        source="evaluation_p4/v1 (unmodified)",
        aggregation="controls averaged within cell -> candidate/control/difference averaged "
                    "within protein over paired positions -> averaged across proteins",
        caution="Candidate means here are restricted to paired cells. Do not compare them "
                "with the all-sampled candidate means in v1/all_sampled_candidate_outcomes.tsv.",
        levels=levels, coverage=cov,
        output_hashes={fn: RD.sha256_file(os.path.join(OUT, fn)) for fn in
                       ("paired_levels.tsv", "paired_coverage_detail.tsv", "paired_cells.tsv")})
    with open(os.path.join(OUT, "addendum_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(dict(levels=levels, coverage=cov), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
