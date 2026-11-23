#!/usr/bin/env python3
"""
Bounded discovery-only redundancy and confounding analysis of the 20 existing
candidate features.

Saved data only: the pooled per-protein summaries written by run_discovery.py,
the saved annotation view, and the SAE checkpoint's decoder weights. No new ESM
extraction, no feature screening, no interventions, no P3/P4, no network.

The hashed candidate list is read and re-verified, never rewritten. All outputs
go to candidate_discovery/v2_redundancy/; v1 is untouched.

Usage:  python3 redundancy_analysis.py
"""

import csv
import json
import os
import sys
import time
from collections import Counter

import numpy as np
from scipy.stats import rankdata, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
V1 = os.path.join(HERE, "v1")
OUT = os.path.join(HERE, "v2_redundancy")
sys.path.insert(0, HERE)
import run_discovery as RD  # noqa: E402

INCLUSION = os.path.join(ROOT, "discovery_analysis", "v1", "inclusion_manifest.tsv")

MIN_PER_CLASS = 10          # below this a cell is reported as sparse, not scored
CLASS_NAME = {1: "Extracellular", 0: "Cytoplasm"}


def auroc(x, y):
    """Signed AUROC with average ranks. >0.5 means higher x associates with y==1.
    Returns None when either class is empty."""
    n1 = int((y == 1).sum())
    n0 = int((y == 0).sum())
    if n1 == 0 or n0 == 0:
        return None
    r = rankdata(x)
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def cell(x, y, label, extra=None):
    n1 = int((y == 1).sum())
    n0 = int((y == 0).sum())
    d = dict(stratum=label, n=int(len(y)), n_Extracellular=n1, n_Cytoplasm=n0)
    if extra:
        d.update(extra)
    if n1 == 0 or n0 == 0:
        d.update(auroc="", status="unsuitable: single class present")
    elif n1 < MIN_PER_CLASS or n0 < MIN_PER_CLASS:
        d.update(auroc="", status=f"unsuitable: sparse (<{MIN_PER_CLASS} in a class)")
    else:
        d.update(auroc=f"{auroc(x, y):.6f}", status="scored")
    return d


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)

    sel_path = os.path.join(V1, "selected_feature_ids.json")
    sel_hash_before = RD.sha256_file(sel_path)
    sel = json.load(open(sel_path))["selected_feature_ids"]
    rank_of = {f: i + 1 for i, f in enumerate(sel)}

    rows, _ = RD.load_sample()
    P = RD.load_pooled(rows)
    X = P["mean_nterm"][:, sel]                      # (N, 20)
    A = (P["nact_nterm"][:, sel] > 0)                # per-protein activity
    y = np.array([r["target"] for r in rows])
    part = np.array([r["partition"] for r in rows])
    king = np.array([r["kingdom"] for r in rows])
    olen = np.array([r["original_len"] for r in rows], dtype=float)
    N = len(rows)

    incl = {r["ACC"]: r for r in RD.read_tsv(INCLUSION)}
    sig = np.array([incl[r["ACC"]]["signal_status"] for r in rows])
    evid = np.array([incl[r["ACC"]]["signal_evidence_best"] for r in rows])
    in_ann = np.array([incl[r["ACC"]]["in_annotated"] == "1" for r in rows])

    # ------------------------------------------------------------ 2. redundancy
    def pairwise_spearman(M):
        rho, _ = spearmanr(M)
        return np.atleast_2d(rho)

    rho_all = pairwise_spearman(X)
    rho_ext = pairwise_spearman(X[y == 1])
    rho_cyt = pairwise_spearman(X[y == 0])

    from safetensors.torch import load_file
    W = load_file(RD.SAE_FILE)["w_dec"].numpy()      # (4096, 1280); row f = feature f
    D = W[sel]                                       # (20, 1280)
    norms = np.linalg.norm(D, axis=1)
    Dn = D / norms[:, None]
    cos = Dn @ Dn.T

    # Activation overlap: Jaccard over the sets of DISCOVERY PROTEINS in which the
    # feature is active at least once within the first 30 positions.
    inter = A.T.astype(np.int64) @ A.astype(np.int64)
    sizes = A.sum(axis=0).astype(np.int64)
    union = sizes[:, None] + sizes[None, :] - inter
    jac = np.where(union > 0, inter / np.maximum(union, 1), 0.0)

    pairs = []
    for i in range(len(sel)):
        for j in range(i + 1, len(sel)):
            pairs.append(
                dict(
                    feature_a=sel[i], feature_b=sel[j],
                    rank_a=rank_of[sel[i]], rank_b=rank_of[sel[j]],
                    spearman_all=f"{rho_all[i, j]:.6f}",
                    spearman_Extracellular=f"{rho_ext[i, j]:.6f}",
                    spearman_Cytoplasm=f"{rho_cyt[i, j]:.6f}",
                    decoder_cosine=f"{cos[i, j]:.6f}",
                    jaccard_active_proteins=f"{jac[i, j]:.6f}",
                    n_active_a=int(sizes[i]), n_active_b=int(sizes[j]),
                    n_active_both=int(inter[i, j]),
                )
            )
    RD.write_tsv(os.path.join(OUT, "candidate_pairwise.tsv"), list(pairs[0].keys()), pairs)

    iu = np.triu_indices(len(sel), 1)
    redundancy = dict(
        n_pairs=len(pairs),
        definition_spearman="Spearman correlation of per-protein mean_nterm across the "
                            "stated protein population (behaviour)",
        definition_decoder_cosine="cosine between SAE decoder rows w_dec[f] (intervention "
                                  "direction in residual space); rows are not unit-norm in "
                                  "this checkpoint, norms reported separately",
        definition_jaccard="|{proteins where A active in the first 30}  INTERSECT  {same for B}| / "
                           "|union of those sets|, over the 5,190 discovery proteins",
        spearman_all=dict(median=float(np.median(rho_all[iu])), min=float(rho_all[iu].min()),
                          max=float(rho_all[iu].max()),
                          n_pairs_above_0p7=int((rho_all[iu] > 0.7).sum()),
                          n_pairs_above_0p9=int((rho_all[iu] > 0.9).sum())),
        spearman_within_Extracellular=dict(median=float(np.median(rho_ext[iu])),
                                           min=float(rho_ext[iu].min()),
                                           max=float(rho_ext[iu].max()),
                                           n_pairs_above_0p7=int((rho_ext[iu] > 0.7).sum())),
        spearman_within_Cytoplasm=dict(median=float(np.median(rho_cyt[iu])),
                                       min=float(rho_cyt[iu].min()),
                                       max=float(rho_cyt[iu].max()),
                                       n_pairs_above_0p7=int((rho_cyt[iu] > 0.7).sum())),
        decoder_cosine=dict(median=float(np.median(cos[iu])), min=float(cos[iu].min()),
                            max=float(cos[iu].max()),
                            n_pairs_above_0p5=int((cos[iu] > 0.5).sum()),
                            n_pairs_above_0p3=int((cos[iu] > 0.3).sum())),
        jaccard=dict(median=float(np.median(jac[iu])), min=float(jac[iu].min()),
                     max=float(jac[iu].max())),
        decoder_row_norms={str(f): float(n) for f, n in zip(sel, norms)},
    )

    # ------------------------------------------------- 3. associations under controls
    strat = []
    quant = np.quantile(olen, [0.25, 0.5, 0.75])     # discovery-only boundaries
    bins = [(0, quant[0]), (quant[0], quant[1]), (quant[1], quant[2]), (quant[2], np.inf)]
    for k, f in enumerate(sel):
        x = X[:, k]
        for p in ("1", "2"):
            mp = part == p
            strat.append(dict(feature_id=f, partition=p, control="none (full partition)",
                              **cell(x[mp], y[mp], "all")))
            m = mp & (olen <= RD.MAX_RESIDUES)
            strat.append(dict(feature_id=f, partition=p, control="original length <= 1022",
                              **cell(x[m], y[m], "not truncated")))
            for kg in sorted(set(king)):
                m = mp & (king == kg)
                strat.append(dict(feature_id=f, partition=p, control="kingdom",
                                  **cell(x[m], y[m], kg)))
            for bi, (lo, hi) in enumerate(bins):
                m = mp & (olen > lo) & (olen <= hi)
                lbl = f"Q{bi+1}: {lo:.0f} < len <= {hi:.0f}" if np.isfinite(hi) else f"Q{bi+1}: len > {lo:.0f}"
                strat.append(dict(feature_id=f, partition=p, control="original-length quartile",
                                  **cell(x[m], y[m], lbl)))
    RD.write_tsv(os.path.join(OUT, "stratified_auroc.tsv"), list(strat[0].keys()), strat)

    strat_meta = dict(
        length_quartile_boundaries_from="the 5,190 discovery proteins only",
        q25=float(quant[0]), q50=float(quant[1]), q75=float(quant[2]),
        min_per_class_to_score=MIN_PER_CLASS,
        pooling="none. Subgroup AUROCs are reported per cell and are NOT averaged; any "
                "pooled figure would require a stated weighting and population, and none "
                "is computed here.",
        caveat="Stratification conditions on one measured variable at a time. It is not "
               "complete confounding removal: residual within-stratum variation, unmeasured "
               "confounders and homology dependence all remain.",
    )

    # ------------------------------------------- 4. signal-annotation counterexamples
    base = in_ann
    cells = [
        ("Extracellular, annotated_present", base & (y == 1) & (sig == "annotated_present")),
        ("Cytoplasm, annotated_present", base & (y == 0) & (sig == "annotated_present")),
        ("Extracellular, not_annotated", base & (y == 1) & (sig == "not_annotated")),
        ("Cytoplasm, not_annotated", base & (y == 0) & (sig == "not_annotated")),
        ("Extracellular, annotated_present, EXPERIMENTAL evidence",
         base & (y == 1) & (sig == "annotated_present") & (evid == "experimental")),
        ("Cytoplasm, annotated_present, EXPERIMENTAL evidence",
         base & (y == 0) & (sig == "annotated_present") & (evid == "experimental")),
    ]
    four = []
    for k, f in enumerate(sel):
        x = X[:, k]
        for label, m in cells:
            v = x[m]
            four.append(
                dict(
                    feature_id=f, cell=label, n=int(m.sum()),
                    prevalence_active=f"{A[m, k].mean():.6f}" if m.sum() else "",
                    mean=f"{v.mean():.6f}" if v.size else "",
                    median=f"{np.median(v):.6f}" if v.size else "",
                    p25=f"{np.percentile(v, 25):.6f}" if v.size else "",
                    p75=f"{np.percentile(v, 75):.6f}" if v.size else "",
                    max=f"{v.max():.6f}" if v.size else "",
                    frac_zero=f"{(v == 0).mean():.6f}" if v.size else "",
                )
            )
    RD.write_tsv(os.path.join(OUT, "signal_annotation_cells.tsv"), list(four[0].keys()), four)
    cell_ns = {label: int(m.sum()) for label, m in cells}

    # ---------------------------------------------- 5. interval-exception bookkeeping
    enr = {r["feature_id"]: r for r in RD.read_tsv(os.path.join(V1, "signal_interval_enrichment.tsv"))}
    exc = []
    for f in sel:
        r = enr[str(f)]
        diff = float(r["mean_per_protein_difference"])
        ing = int(r["n_proteins_inside_greater"]); outg = int(r["n_proteins_outside_greater"])
        exc.append(
            dict(
                feature_id=f, rank=rank_of[f],
                mean_inside=r["mean_inside_interval"], mean_outside=r["mean_outside_interval"],
                mean_per_protein_difference=r["mean_per_protein_difference"],
                median_per_protein_difference=r["median_per_protein_difference"],
                n_proteins=r["n_proteins"],
                n_inside_greater=ing, n_outside_greater=outg,
                frac_inside_greater=f"{ing / int(r['n_proteins']):.4f}",
                n_zero_inside=r["n_proteins_zero_inside"],
                enriched_by_mean="1" if diff > 0 else "0",
                enriched_by_majority="1" if ing > outg else "0",
            )
        )
    RD.write_tsv(os.path.join(OUT, "interval_exceptions.tsv"), list(exc[0].keys()), exc)

    by_mean = [e["feature_id"] for e in exc if e["enriched_by_mean"] == "0"]
    by_maj = [e["feature_id"] for e in exc if e["enriched_by_majority"] == "0"]
    weakest = sorted(exc, key=lambda e: float(e["frac_inside_greater"]))[:3]

    summary = dict(
        selected_feature_ids=sel,
        selected_list_sha256=sel_hash_before,
        candidate_list_unchanged=RD.sha256_file(sel_path) == sel_hash_before,
        n_proteins=N, n_Extracellular=int((y == 1).sum()), n_Cytoplasm=int((y == 0).sum()),
        redundancy=redundancy,
        stratification=strat_meta,
        signal_cell_counts=cell_ns,
        interval_exceptions=dict(
            not_enriched_by_mean_difference=by_mean,
            not_enriched_by_majority_of_proteins=by_maj,
            three_weakest_by_fraction_inside_greater=[
                dict(feature_id=e["feature_id"], frac_inside_greater=e["frac_inside_greater"],
                     n_zero_inside=e["n_zero_inside"]) for e in weakest
            ],
            comparison_region=(
                "FULL VISIBLE SEQUENCE. 'outside' is every visible residue not in "
                "[start, end], so the two regions differ in length, in position and in "
                "sequence context. This is not a matched-position window comparison, and "
                "the saved artefacts cannot produce one: per-residue activations were "
                "streamed and not retained."
            ),
        ),
        inputs=dict(
            pooled_dir="candidate_discovery/v1/pooled",
            selected_list=sel_hash_before,
            sae_file_sha256=RD.SAE_SHA,
            inclusion_manifest=RD.sha256_file(INCLUSION),
        ),
        network_requests=0,
        new_esm_extraction=False,
        seconds=time.time() - t0,
    )
    summary["output_hashes"] = {
        f: RD.sha256_file(os.path.join(OUT, f))
        for f in ("candidate_pairwise.tsv", "stratified_auroc.tsv",
                  "signal_annotation_cells.tsv", "interval_exceptions.tsv")
    }
    assert summary["candidate_list_unchanged"]
    with open(os.path.join(OUT, "redundancy_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps({k: v for k, v in summary.items() if k != "redundancy"}, indent=2))
    print(json.dumps(redundancy, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
