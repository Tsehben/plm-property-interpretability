#!/usr/bin/env python3
"""
Aggregation for the P3 intervention pilot. Reads only v1/per_protein/*.json and
v1/frozen_configuration.json. Writes to v1/. No model is loaded and no P4 or P0
data is touched.

Aggregation rule, applied everywhere: positions are averaged WITHIN a protein
first, then proteins are the unit across summaries. Positions and control
features are never counted as independent proteins.

Usage:  python3 analyse_pilot.py
"""

import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "v1")
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "candidate_discovery"))
import run_discovery as RD  # noqa: E402

METRICS = [
    "kl_aa20_renormalised_nats", "kl_full_vocab_nats", "total_variation_aa20",
    "delta_p_true", "delta_logp_true", "delta_expected_kd",
    "delta_nonstandard_mass", "delta_h_l2", "delta_h_l2_over_h_l2",
]
NO_OP = 1.0


def f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def load():
    rows = []
    meta = []
    for p in sorted(glob.glob(os.path.join(OUT, "per_protein", "*.json"))):
        d = json.load(open(p))
        meta.append(dict(ACC=d["protein"], target=d["target"], original_len=d["original_len"],
                         positions=d["positions"], duplicates=d["duplicate_or_invalid_positions"],
                         token_checks_ok=all(c["match"] for c in d["token_offset_checks"])))
        rows += d["rows"]
    return rows, meta


def stats(a):
    a = np.asarray([x for x in a if np.isfinite(x)], dtype=float)
    if a.size == 0:
        return dict(n_proteins=0)
    return dict(n_proteins=int(a.size), mean=float(a.mean()),
                sd=float(a.std(ddof=1)) if a.size > 1 else 0.0,
                min=float(a.min()), p25=float(np.percentile(a, 25)),
                median=float(np.median(a)), p75=float(np.percentile(a, 75)),
                max=float(a.max()))


def per_protein_means(rows, keyfn, metric, filt=None):
    """Average the metric over positions within each protein, then return one
    value per protein."""
    acc = defaultdict(list)
    for r in rows:
        if filt and not filt(r):
            continue
        acc[(keyfn(r), r["ACC"])].append(f(r.get(metric)))
    out = defaultdict(list)
    for (k, a), vals in acc.items():
        v = [x for x in vals if np.isfinite(x)]
        if v:
            out[k].append(float(np.mean(v)))
    return out


def main():
    rows, meta = load()
    cfg = json.load(open(os.path.join(OUT, "frozen_configuration.json")))
    prot_target = {m["ACC"]: m["target"] for m in meta}
    cls = lambda r: "Extracellular" if prot_target[r["ACC"]] == 1 else "Cytoplasm"

    edited = [r for r in rows if r["status"] == "edited"]
    cand_rows = [r for r in rows if r["arm"] == "candidate"]
    ctrl_rows = [r for r in rows if r["arm"] == "control"]

    # ------------------------------------------------- integrity: no-op and edit
    noop = [r for r in edited if float(r["multiplier"]) == NO_OP]
    nonzero_noop = [r for r in noop if f(r.get("delta_h_l2")) != 0.0]
    integ = dict(
        n_rows=len(rows), n_edited=len(edited),
        n_no_op_conditions=len(noop),
        no_op_delta_h_all_exactly_zero=len(nonzero_noop) == 0,
        no_op_max_delta_h_l2=max([f(r.get("delta_h_l2")) for r in noop], default=0.0),
        no_op_max_kl_aa20=max([f(r.get("kl_aa20_renormalised_nats")) for r in noop], default=0.0),
        no_op_max_kl_full=max([f(r.get("kl_full_vocab_nats")) for r in noop], default=0.0),
        n_latent_coords_changed_max=max([int(r.get("n_latent_coords_changed", 0)) for r in edited], default=0),
        support_evicted_unintended_total=sum(int(r.get("support_evicted_unintended", 0)) for r in edited),
        support_inserted_unintended_total=sum(int(r.get("support_inserted_unintended", 0)) for r in edited),
        delta_confined_to_edited_position_all=all(
            int(r.get("delta_confined_to_edited_position", 0)) == 1 for r in edited),
        n_delta_not_confined=sum(1 for r in edited
                                 if int(r.get("delta_confined_to_edited_position", 0)) != 1),
        non_finite_metric_cells=sum(
            1 for r in edited for m in METRICS
            if m in r and not np.isfinite(f(r[m]))),
        token_offset_checks_all_passed=all(m["token_checks_ok"] for m in meta),
        proteins_with_duplicate_positions=[m["ACC"] for m in meta if m["duplicates"]],
    )

    # ------------------------------------------------------------ availability
    avail = []
    for arm, rs in (("candidate", cand_rows), ("control", ctrl_rows)):
        byst = defaultdict(int)
        for r in rs:
            byst[r["status"]] += 1
        for st, n in sorted(byst.items(), key=lambda t: -t[1]):
            avail.append(dict(arm=arm, status=st, n_conditions=n,
                              pct=f"{100.0 * n / len(rs):.1f}"))
    RD.write_tsv(os.path.join(OUT, "condition_availability.tsv"),
                 ["arm", "status", "n_conditions", "pct"], avail)

    # per candidate: active-edit denominators
    dens = []
    for c in cfg["config"]["candidates"]:
        cr = [r for r in cand_rows if int(r["candidate"]) == c and float(r["multiplier"]) != NO_OP]
        ed = [r for r in cr if r["status"] == "edited"]
        prot_active = len(set(r["ACC"] for r in ed))
        cells = set((r["ACC"], r["residue_position"]) for r in cr)
        cells_act = set((r["ACC"], r["residue_position"]) for r in ed)
        dens.append(dict(
            candidate=c, n_sampled_conditions=len(cr), n_active_edit_conditions=len(ed),
            pct_active=f"{100.0 * len(ed) / max(len(cr), 1):.1f}",
            n_protein_position_cells=len(cells), n_cells_with_active_candidate=len(cells_act),
            pct_cells_active=f"{100.0 * len(cells_act) / max(len(cells), 1):.1f}",
            n_proteins_total=len(meta), n_proteins_with_any_active_cell=prot_active,
        ))
    RD.write_tsv(os.path.join(OUT, "active_edit_denominators.tsv"), list(dens[0].keys()), dens)

    # ----------------------------------------------------------- dose response
    dose = []
    for c in cfg["config"]["candidates"]:
        for m in cfg["config"]["multipliers"]:
            filt = lambda r, c=c, m=m: (r["arm"] == "candidate" and int(r["candidate"]) == c
                                        and float(r["multiplier"]) == m and r["status"] == "edited")
            for metric in METRICS:
                d = per_protein_means(rows, lambda r: "all", metric, filt)
                s = stats(d.get("all", []))
                dose.append(dict(candidate=c, multiplier=m, metric=metric, stratum="all", **s))
            for metric in ("kl_aa20_renormalised_nats", "delta_p_true", "delta_expected_kd"):
                d = per_protein_means(rows, cls, metric, filt)
                for k in ("Extracellular", "Cytoplasm"):
                    dose.append(dict(candidate=c, multiplier=m, metric=metric, stratum=k,
                                     **stats(d.get(k, []))))
    RD.write_tsv(os.path.join(OUT, "dose_response.tsv"), list(dose[0].keys()), dose)

    # ------------------------------------------------------------ by position
    bypos = []
    for c in cfg["config"]["candidates"]:
        for m in cfg["config"]["multipliers"]:
            if m == NO_OP:
                continue
            for metric in ("kl_aa20_renormalised_nats", "delta_expected_kd"):
                filt = lambda r, c=c, m=m: (r["arm"] == "candidate" and int(r["candidate"]) == c
                                            and float(r["multiplier"]) == m and r["status"] == "edited")
                acc = defaultdict(list)
                for r in rows:
                    if filt(r):
                        acc[str(r["residue_position"]) if int(r["residue_position"]) in (10, 20, 30)
                            else "floor(L/2)"].append(f(r.get(metric)))
                for k, v in sorted(acc.items()):
                    bypos.append(dict(candidate=c, multiplier=m, metric=metric,
                                      mask_position=k, n_conditions=len(v), **stats(v)))
    RD.write_tsv(os.path.join(OUT, "by_mask_position.tsv"), list(bypos[0].keys()), bypos)

    # --------------------------------------- paired candidate-vs-control contrast
    # Pair only where BOTH the candidate condition and that control condition are
    # 'edited' at the same (protein, position, multiplier). Controls are averaged
    # within protein first (over positions and over control features), so the unit
    # is the protein, never the control feature.
    paired = []
    ckey = {}
    for r in cand_rows:
        if r["status"] == "edited":
            ckey[(r["ACC"], r["residue_position"], int(r["candidate"]), float(r["multiplier"]))] = r
    for c in cfg["config"]["candidates"]:
        ctrls = cfg["controls"][str(c)]
        if not ctrls:
            paired.append(dict(candidate=c, multiplier="", metric="", n_paired_proteins=0,
                               n_paired_conditions=0, mean_candidate="", mean_control="",
                               mean_paired_difference="", median_paired_difference="",
                               n_proteins_candidate_greater="", note="no matched control exists"))
            continue
        for m in cfg["config"]["multipliers"]:
            if m == NO_OP:
                continue
            pairs = defaultdict(lambda: defaultdict(list))
            n_cond = 0
            for r in ctrl_rows:
                if r["status"] != "edited" or int(r["candidate"]) != c or float(r["multiplier"]) != m:
                    continue
                k = (r["ACC"], r["residue_position"], c, m)
                cr = ckey.get(k)
                if cr is None:
                    continue
                n_cond += 1
                for metric in METRICS:
                    pairs[metric][r["ACC"]].append((f(cr.get(metric)), f(r.get(metric))))
            for metric in METRICS:
                per = []
                for a, vs in pairs[metric].items():
                    cand_v = np.mean([x for x, _ in vs if np.isfinite(x)])
                    ctrl_v = np.mean([y for _, y in vs if np.isfinite(y)])
                    if np.isfinite(cand_v) and np.isfinite(ctrl_v):
                        per.append((cand_v, ctrl_v))
                if not per:
                    continue
                cv = np.array([x for x, _ in per]); kv = np.array([y for _, y in per])
                d = cv - kv
                paired.append(dict(
                    candidate=c, multiplier=m, metric=metric,
                    n_paired_proteins=len(per), n_paired_conditions=n_cond,
                    mean_candidate=float(cv.mean()), mean_control=float(kv.mean()),
                    mean_paired_difference=float(d.mean()),
                    median_paired_difference=float(np.median(d)),
                    n_proteins_candidate_greater=int((d > 0).sum()),
                    note="paired within protein; controls averaged over positions and control "
                         "features before pairing",
                ))
    RD.write_tsv(os.path.join(OUT, "candidate_vs_control_paired.tsv"), list(paired[0].keys()), paired)

    # ------------------------------------------------------- realised norm match
    nm = []
    for c in cfg["config"]["candidates"]:
        for m in cfg["config"]["multipliers"]:
            if m == NO_OP:
                continue
            cd, kd_ = [], []
            for r in ctrl_rows:
                if r["status"] != "edited" or int(r["candidate"]) != c or float(r["multiplier"]) != m:
                    continue
                cr = ckey.get((r["ACC"], r["residue_position"], c, m))
                if cr is None:
                    continue
                cd.append(f(cr.get("delta_h_l2"))); kd_.append(f(r.get("delta_h_l2")))
            if not cd:
                continue
            cd, kd_ = np.array(cd), np.array(kd_)
            rel = np.abs(cd - kd_) / np.maximum(cd, 1e-12)
            nm.append(dict(candidate=c, multiplier=m, n_pairs=len(cd),
                           mean_candidate_delta_h_l2=float(cd.mean()),
                           mean_control_delta_h_l2=float(kd_.mean()),
                           median_rel_norm_difference=float(np.median(rel)),
                           p90_rel_norm_difference=float(np.percentile(rel, 90)),
                           max_rel_norm_difference=float(rel.max()),
                           note="latent delta matched exactly in magnitude and sign; the "
                                "residual-space norm is NOT matched and is reported"))
    RD.write_tsv(os.path.join(OUT, "realised_norm_match.tsv"), list(nm[0].keys()), nm)

    # ------------------------------------------------------------ per protein
    pp = []
    for mm in meta:
        for c in cfg["config"]["candidates"]:
            for m in cfg["config"]["multipliers"]:
                sub = [r for r in cand_rows if r["ACC"] == mm["ACC"] and int(r["candidate"]) == c
                       and float(r["multiplier"]) == m]
                ed = [r for r in sub if r["status"] == "edited"]
                row = dict(ACC=mm["ACC"],
                           localisation_class="Extracellular" if mm["target"] == 1 else "Cytoplasm",
                           original_len=mm["original_len"], candidate=c, multiplier=m,
                           n_positions_sampled=len(sub), n_positions_active=len(ed))
                for metric in ("kl_aa20_renormalised_nats", "delta_p_true",
                               "delta_expected_kd", "delta_h_l2_over_h_l2"):
                    v = [f(r.get(metric)) for r in ed]
                    v = [x for x in v if np.isfinite(x)]
                    row[f"mean_{metric}"] = f"{np.mean(v):.17g}" if v else ""
                pp.append(row)
    RD.write_tsv(os.path.join(OUT, "per_protein_summary.tsv"), list(pp[0].keys()), pp)

    summary = dict(integrity=integ,
                   n_proteins=len(meta),
                   n_Extracellular=sum(1 for m in meta if m["target"] == 1),
                   n_Cytoplasm=sum(1 for m in meta if m["target"] == 0),
                   active_edit_denominators=dens,
                   control_shortfall=cfg["control_shortfall"],
                   output_hashes={fn: RD.sha256_file(os.path.join(OUT, fn)) for fn in (
                       "condition_availability.tsv", "active_edit_denominators.tsv",
                       "dose_response.tsv", "by_mask_position.tsv",
                       "candidate_vs_control_paired.tsv", "realised_norm_match.tsv",
                       "per_protein_summary.tsv")})
    with open(os.path.join(OUT, "analysis_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(dict(integrity=integ, denominators=dens), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
