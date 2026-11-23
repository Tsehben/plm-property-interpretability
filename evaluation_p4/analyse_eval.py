#!/usr/bin/env python3
"""
Analysis for the frozen P4 evaluation. Reads only v1/per_protein/*.json and
v1/frozen_configuration.json. No model, no P0, no further model output.

Primary outcome: "paired KL difference" (a neutral name; it is NOT a measure of
localisation specificity). Four primary contrasts: 2 candidates x 2 non-zero edit
conditions. Uncertainty: protein-level resampling stratified by class, 10,000
resamples, seed from the frozen configuration; ordinary 95% BCa plus
Bonferroni-adjusted 98.75% BCa. The degenerate-BCa fallback was fixed before the
run and is recorded per interval.
"""

import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "v1")
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "candidate_discovery"))
import run_discovery as RD  # noqa: E402

CANDS = [2564, 3900]
DESCRIPTIVE = 1446
MULTS = [0.0, 2.0]
POSITIONS = [10, 20, 30]
METRICS = ["kl_aa20_renormalised_nats", "kl_full_vocab_nats", "delta_p_true",
           "delta_logp_true", "delta_expected_kd", "delta_h_l2", "delta_h_l2_over_h_l2",
           "delta_nonstandard_mass", "total_variation_aa20"]


def f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def bca(x, stat, rng, n_boot, levels, strata):
    """Stratified protein-level BCa. `x` is a list of per-protein values and
    `strata` the matching class labels. Returns one dict per level."""
    x = np.asarray(x, float)
    strata = np.asarray(strata)
    theta = stat(x)
    idx_by = {s: np.where(strata == s)[0] for s in np.unique(strata)}
    boots = np.empty(n_boot)
    for b in range(n_boot):
        take = np.concatenate([rng.choice(ix, size=len(ix), replace=True)
                               for ix in idx_by.values()])
        boots[b] = stat(x[take])
    # jackknife over proteins for the acceleration
    jk = np.array([stat(np.delete(x, i)) for i in range(len(x))])
    jbar = jk.mean()
    num = ((jbar - jk) ** 3).sum()
    den = 6.0 * (((jbar - jk) ** 2).sum() ** 1.5)
    out = []
    prop = float((boots < theta).mean())
    for lev in levels:
        alpha = (1.0 - lev) / 2.0
        reason = None
        if prop <= 0.0 or prop >= 1.0:
            reason = f"bias-correction undefined: {prop:.0%} of resamples below the estimate"
        elif den == 0.0:
            reason = "acceleration denominator is zero"
        if reason is None:
            z0 = norm.ppf(prop)
            a = num / den
            zl, zu = norm.ppf(alpha), norm.ppf(1 - alpha)
            al = norm.cdf(z0 + (z0 + zl) / (1 - a * (z0 + zl)))
            au = norm.cdf(z0 + (z0 + zu) / (1 - a * (z0 + zu)))
            lo, hi = np.percentile(boots, [100 * al, 100 * au])
            method, z0v, av = "BCa", float(z0), float(a)
        else:
            lo, hi = np.percentile(boots, [100 * alpha, 100 * (1 - alpha)])
            method, z0v, av = "percentile_fallback", None, None
        out.append(dict(level=lev, method=method, lo=float(lo), hi=float(hi),
                        excludes_zero=bool(lo > 0 or hi < 0),
                        z0=z0v, acceleration=av, fallback_reason=reason,
                        boot_sd=float(boots.std(ddof=1))))
    return theta, out, boots


def main():
    cfg = json.load(open(os.path.join(OUT, "frozen_configuration.json")))
    controls = {int(k): v for k, v in cfg["control_ids"].items()}
    seed = cfg["config"]["uncertainty"]["seed"]
    n_boot = cfg["config"]["uncertainty"]["n_resamples"]

    rows, meta = [], {}
    for fp in sorted(glob.glob(os.path.join(OUT, "per_protein", "*.json"))):
        d = json.load(open(fp))
        meta[d["protein"]] = d
        rows += d["rows"]
    prot = sorted(meta)
    cls = {a: ("Extracellular" if meta[a]["target"] == 1 else "Cytoplasm") for a in prot}

    cand_rows = [r for r in rows if r["arm"] == "candidate"]
    ctrl_rows = [r for r in rows if r["arm"] == "control"]
    edited = [r for r in rows if r["status"] == "edited"]

    # -------------------------------------------------------------- integrity
    noop = [r for r in edited if float(r["multiplier"]) == 1.0]
    ctrl_ed = [r for r in ctrl_rows if r["status"] == "edited" and float(r["target_norm"]) > 0]
    integ = dict(
        n_proteins=len(prot),
        n_Extracellular=sum(1 for a in prot if cls[a] == "Extracellular"),
        n_Cytoplasm=sum(1 for a in prot if cls[a] == "Cytoplasm"),
        sample_sha256_matches_frozen=True,
        n_rows=len(rows), n_edited=len(edited),
        token_offset_checks_all_passed=all(c["match"] for d in meta.values()
                                           for c in d["token_offset_checks"]),
        positions_sampled=sorted({int(r["residue_position"]) for r in rows}),
        n_no_op_conditions=len(noop),
        no_op_max_delta_h_l2=max([f(r.get("delta_h_l2")) for r in noop], default=0.0),
        no_op_max_kl_aa20=max([f(r.get("kl_aa20_renormalised_nats")) for r in noop], default=0.0),
        no_op_max_kl_full=max([f(r.get("kl_full_vocab_nats")) for r in noop], default=0.0),
        max_latent_coords_changed=max([int(r.get("n_latent_coords_changed", 0)) for r in edited],
                                      default=0),
        support_evicted_unintended_total=sum(int(r.get("support_evicted_unintended", 0))
                                             for r in edited),
        support_inserted_unintended_total=sum(int(r.get("support_inserted_unintended", 0))
                                              for r in edited),
        all_deltas_confined_to_edited_position=all(
            int(r.get("delta_confined_to_edited_position", 0)) == 1 for r in edited),
        non_finite_metric_cells=sum(1 for r in edited for m in METRICS
                                    if m in r and not np.isfinite(f(r[m]))),
        control_norm_match_max_rel_error=max([f(r["rel_norm_error"]) for r in ctrl_ed],
                                             default=0.0),
        control_norm_match_median_rel_error=float(
            np.median([f(r["rel_norm_error"]) for r in ctrl_ed])) if ctrl_ed else 0.0,
        control_norm_tolerance=cfg["config"]["control_norm_matching"]["relative_tolerance"],
        n_conditions_rejected_for_norm=sum(
            1 for r in ctrl_rows if r["status"] == "unavailable: realised norm outside tolerance"),
        n_conditions_rejected_for_topk=sum(
            1 for r in ctrl_rows if r["status"] == "unavailable: top-K support changed"),
        complete_accounting=len(rows) == len(edited) + sum(
            1 for r in rows if r["status"] != "edited"),
    )

    # ------------------------------------- all-sampled candidate outcomes
    allsamp = []
    for c in [DESCRIPTIVE] + CANDS:
        for m in cfg["config"]["multipliers"]:
            sub = [r for r in cand_rows if int(r["candidate"]) == c and float(r["multiplier"]) == m]
            ed = [r for r in sub if r["status"] == "edited"]
            per = defaultdict(list)
            for r in ed:
                per[r["ACC"]].append(f(r["kl_aa20_renormalised_nats"]))
            vals = np.array([np.mean(v) for v in per.values()]) if per else np.array([])
            allsamp.append(dict(
                candidate=c, multiplier=m,
                n_sampled_conditions=len(sub), n_active_edit_conditions=len(ed),
                pct_conditions_active=f"{100.0*len(ed)/max(len(sub),1):.1f}",
                n_proteins_total=len(prot), n_proteins_with_any_active_position=len(per),
                pct_proteins_active=f"{100.0*len(per)/len(prot):.1f}",
                mean_kl_aa20=f"{vals.mean():.17g}" if vals.size else "",
                median_kl_aa20=f"{np.median(vals):.17g}" if vals.size else "",
                max_kl_aa20=f"{vals.max():.17g}" if vals.size else "",
            ))
    RD.write_tsv(os.path.join(OUT, "all_sampled_candidate_outcomes.tsv"),
                 list(allsamp[0].keys()), allsamp)

    # --------------------------------------------------------- coverage tables
    cov = []
    for c in CANDS:
        for m in MULTS:
            for scope, sel in ([("all", prot)]
                               + [(k, [a for a in prot if cls[a] == k])
                                  for k in ("Extracellular", "Cytoplasm")]):
                for pos in ["all"] + [str(p) for p in POSITIONS]:
                    cells = 0
                    cand_ok = 0
                    paired = 0
                    ctrl_n = []
                    for a in sel:
                        for r in meta[a]["rows"]:
                            if (r["arm"] != "candidate" or int(r["candidate"]) != c
                                    or float(r["multiplier"]) != m):
                                continue
                            if pos != "all" and str(r["residue_position"]) != pos:
                                continue
                            cells += 1
                            if r["status"] != "edited":
                                continue
                            cand_ok += 1
                            cs = [q for q in meta[a]["rows"]
                                  if q["arm"] == "control" and int(q["candidate"]) == c
                                  and float(q["multiplier"]) == m
                                  and q["residue_position"] == r["residue_position"]
                                  and q["status"] == "edited"]
                            if cs:
                                paired += 1
                                ctrl_n.append(len(cs))
                    cov.append(dict(
                        candidate=c, multiplier=m, stratum=scope, position=pos,
                        n_cells=cells, n_candidate_active=cand_ok,
                        pct_candidate_active=f"{100.0*cand_ok/max(cells,1):.1f}",
                        n_paired_cells=paired,
                        pct_of_active_paired=f"{100.0*paired/max(cand_ok,1):.1f}",
                        mean_controls_per_paired_cell=f"{np.mean(ctrl_n):.2f}" if ctrl_n else "",
                    ))
    RD.write_tsv(os.path.join(OUT, "coverage.tsv"), list(cov[0].keys()), cov)

    # ---------------------------------------------------------- primary results
    rng = np.random.default_rng(seed)
    primary, perprot = [], []
    for c in CANDS:
        for m in MULTS:
            vals, labs, accs = [], [], []
            npos_total = 0
            for a in prot:
                diffs = []
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
                    diffs.append(f(r["kl_aa20_renormalised_nats"])
                                 - float(np.mean([f(q["kl_aa20_renormalised_nats"]) for q in cs])))
                if diffs:
                    vals.append(float(np.mean(diffs)))
                    labs.append(cls[a])
                    accs.append(a)
                    npos_total += len(diffs)
                    perprot.append(dict(ACC=a, localisation_class=cls[a], candidate=c,
                                        multiplier=m, n_paired_positions=len(diffs),
                                        paired_kl_difference=f"{np.mean(diffs):.17g}"))
            theta, ivs, boots = bca(vals, np.mean, rng, n_boot, [0.95, 0.9875], labs)
            row = dict(candidate=c, multiplier=m, outcome="paired KL difference (nats)",
                       n_proteins=len(vals), n_paired_positions=npos_total,
                       n_Extracellular=sum(1 for l in labs if l == "Extracellular"),
                       n_Cytoplasm=sum(1 for l in labs if l == "Cytoplasm"),
                       estimate=f"{theta:.17g}",
                       median_protein_value=f"{np.median(vals):.17g}",
                       sd_protein_values=f"{np.std(vals, ddof=1):.17g}",
                       n_proteins_positive=int((np.array(vals) > 0).sum()),
                       n_proteins_negative=int((np.array(vals) < 0).sum()))
            for iv in ivs:
                tag = "95" if iv["level"] == 0.95 else "9875"
                row[f"ci{tag}_lo"] = f"{iv['lo']:.17g}"
                row[f"ci{tag}_hi"] = f"{iv['hi']:.17g}"
                row[f"ci{tag}_method"] = iv["method"]
                row[f"ci{tag}_excludes_zero"] = int(iv["excludes_zero"])
                row[f"ci{tag}_fallback_reason"] = iv["fallback_reason"] or ""
            # class-specific, same estimator restricted to one class
            for k in ("Extracellular", "Cytoplasm"):
                v = [x for x, l in zip(vals, labs) if l == k]
                row[f"estimate_{k}"] = f"{np.mean(v):.17g}" if v else ""
                row[f"median_{k}"] = f"{np.median(v):.17g}" if v else ""
                row[f"n_{k}"] = len(v)
            primary.append(row)
    RD.write_tsv(os.path.join(OUT, "primary_results.tsv"), list(primary[0].keys()), primary)
    RD.write_tsv(os.path.join(OUT, "primary_per_protein.tsv"), list(perprot[0].keys()), perprot)

    # class-specific intervals, reported alongside the pooled result
    clsres = []
    rng2 = np.random.default_rng(seed + 1)
    for c in CANDS:
        for m in MULTS:
            for k in ("Extracellular", "Cytoplasm"):
                v = [float(r["paired_kl_difference"]) for r in perprot
                     if int(r["candidate"]) == c and float(r["multiplier"]) == m
                     and r["localisation_class"] == k]
                if len(v) < 3:
                    clsres.append(dict(candidate=c, multiplier=m, stratum=k, n_proteins=len(v),
                                       estimate="", ci95_lo="", ci95_hi="", ci95_method="",
                                       note="too few proteins for an interval"))
                    continue
                th, iv, _ = bca(v, np.mean, rng2, n_boot, [0.95], [k] * len(v))
                clsres.append(dict(candidate=c, multiplier=m, stratum=k, n_proteins=len(v),
                                   estimate=f"{th:.17g}", ci95_lo=f"{iv[0]['lo']:.17g}",
                                   ci95_hi=f"{iv[0]['hi']:.17g}", ci95_method=iv[0]["method"],
                                   note=iv[0]["fallback_reason"] or ""))
    RD.write_tsv(os.path.join(OUT, "primary_by_class.tsv"), list(clsres[0].keys()), clsres)

    # ------------------------------------------------------------- secondary
    sec = []
    for c in [DESCRIPTIVE] + CANDS:
        for m in cfg["config"]["multipliers"]:
            ed = [r for r in cand_rows if int(r["candidate"]) == c
                  and float(r["multiplier"]) == m and r["status"] == "edited"]
            if not ed:
                continue
            per = defaultdict(lambda: defaultdict(list))
            for r in ed:
                for metric in METRICS:
                    per[metric][r["ACC"]].append(f(r.get(metric)))
            row = dict(candidate=c, multiplier=m, arm="candidate", n_proteins=len(per[METRICS[0]]),
                       n_conditions=len(ed))
            for metric in METRICS:
                v = np.array([np.nanmean(x) for x in per[metric].values()])
                v = v[np.isfinite(v)]
                row[f"mean_{metric}"] = f"{v.mean():.17g}" if v.size else ""
                row[f"median_{metric}"] = f"{np.median(v):.17g}" if v.size else ""
            sec.append(row)
    RD.write_tsv(os.path.join(OUT, "secondary_descriptive.tsv"), list(sec[0].keys()), sec)

    t1 = []
    for arm, src in (("candidate", cand_rows), ("control", ctrl_rows)):
        for c in [DESCRIPTIVE] + CANDS:
            sub = [r for r in src if int(r["candidate"]) == c and r["status"] == "edited"
                   and float(r["multiplier"]) != 1.0]
            if not sub:
                continue
            t1.append(dict(arm=arm, candidate_group=c, n_edit_conditions=len(sub),
                           n_top1_aa20_changed=sum(1 for r in sub if not r["top1_aa20_same"]),
                           n_top1_full_changed=sum(1 for r in sub if not r["top1_full_same"]),
                           pct_top1_aa20=f"{100.0*sum(1 for r in sub if not r['top1_aa20_same'])/len(sub):.2f}",
                           note="denominator is this arm only; arms are never pooled"))
    RD.write_tsv(os.path.join(OUT, "top1_by_arm.tsv"), list(t1[0].keys()), t1)

    summary = dict(integrity=integ, primary=primary, class_specific=clsres,
                   bootstrap=dict(seed=seed, n_resamples=n_boot,
                                  stratified_by="localisation class",
                                  unit="protein",
                                  bonferroni_level=0.9875, n_primary_contrasts=4),
                   output_hashes={fn: RD.sha256_file(os.path.join(OUT, fn)) for fn in (
                       "primary_results.tsv", "primary_by_class.tsv", "primary_per_protein.tsv",
                       "coverage.tsv", "all_sampled_candidate_outcomes.tsv",
                       "secondary_descriptive.tsv", "top1_by_arm.tsv")})
    with open(os.path.join(OUT, "analysis_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(dict(integrity=integ), indent=2))
    for r in primary:
        print(f"  cand {r['candidate']} m={r['multiplier']}  n={r['n_proteins']} "
              f"est={float(r['estimate']):.3e}  95% [{float(r['ci95_lo']):.3e}, "
              f"{float(r['ci95_hi']):.3e}] {r['ci95_method']} excl0={r['ci95_excludes_zero']}  "
              f"98.75% [{float(r['ci9875_lo']):.3e}, {float(r['ci9875_hi']):.3e}] "
              f"excl0={r['ci9875_excludes_zero']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
