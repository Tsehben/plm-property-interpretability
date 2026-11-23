#!/usr/bin/env python3
"""
Paired contrasts for the norm-matched correction. Reads v1/per_protein/*.json
(candidate arm, reused unchanged) and v2_normmatched/per_protein/*.json (control
arm, rebuilt). Writes to v2_normmatched/. No model, no P4, no P0.

Aggregation, applied exactly in this order and reported at every level:
  1. within a (protein, position, candidate, multiplier) cell, average the
     AVAILABLE control features' effects;
  2. form candidate minus that control average at that position;
  3. average those paired differences within each protein;
  4. summarise across proteins.

A position with no available control contributes nothing and is counted.
"""

import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
V1 = os.path.join(HERE, "v1")
V2 = os.path.join(HERE, "v2_normmatched")
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "candidate_discovery"))
import run_discovery as RD  # noqa: E402

METRICS = ["kl_aa20_renormalised_nats", "kl_full_vocab_nats", "total_variation_aa20",
           "delta_p_true", "delta_logp_true", "delta_expected_kd", "delta_h_l2"]
CANDS = [2564, 3900]
MULTS = [0.0, 0.5, 1.5, 2.0]


def f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def main():
    cfg = json.load(open(os.path.join(V1, "frozen_configuration.json")))
    controls = {int(k): v for k, v in cfg["controls"].items()}

    cand = {}
    targets = {}
    for fp in sorted(glob.glob(os.path.join(V1, "per_protein", "*.json"))):
        d = json.load(open(fp))
        targets[d["protein"]] = d["target"]
        for r in d["rows"]:
            if r["arm"] == "candidate" and r["status"] == "edited":
                cand[(r["ACC"], r["residue_position"], int(r["candidate"]), float(r["multiplier"]))] = r

    ctrl = defaultdict(list)
    status = defaultdict(lambda: defaultdict(int))
    normerr = defaultdict(list)
    for fp in sorted(glob.glob(os.path.join(V2, "per_protein", "*.json"))):
        d = json.load(open(fp))
        for r in d["rows"]:
            status[(int(r["candidate"]), float(r["multiplier"]))][r["status"]] += 1
            if r["status"] == "edited":
                ctrl[(r["ACC"], r["residue_position"], int(r["candidate"]),
                      float(r["multiplier"]))].append(r)
                if f(r.get("target_norm")) > 0:
                    normerr[(int(r["candidate"]), float(r["multiplier"]))].append(
                        f(r["rel_norm_error"]))

    # ------------------------------------------------- availability, v1 vs v2
    v1_avail = defaultdict(int)
    for fp in sorted(glob.glob(os.path.join(V1, "per_protein", "*.json"))):
        for r in json.load(open(fp))["rows"]:
            if r["arm"] == "control" and r["status"] == "edited":
                v1_avail[(int(r["candidate"]), float(r["multiplier"]))] += 1
    av = []
    for c in CANDS:
        for m in MULTS:
            st = status[(c, m)]
            n_edit = st.get("edited", 0)
            av.append(dict(
                candidate=c, multiplier=m,
                n_control_conditions_total=sum(st.values()),
                n_available_v2=n_edit, n_available_v1_latent_matched=v1_avail[(c, m)],
                lost_vs_v1=v1_avail[(c, m)] - n_edit,
                unavailable_control_inactive=st.get("unavailable: control inactive at this position", 0),
                unavailable_negative_latent=st.get(
                    "unavailable: required negative latent change exceeds the control's latent value", 0),
                unavailable_topk=st.get("unavailable: top-K support changed", 0),
                unavailable_norm_tolerance=st.get("unavailable: realised norm outside tolerance", 0),
                candidate_inactive=st.get("candidate_inactive_no_edit", 0),
                max_rel_norm_error=f"{max(normerr[(c, m)], default=0.0):.3e}",
                median_rel_norm_error=f"{np.median(normerr[(c, m)]) if normerr[(c, m)] else 0.0:.3e}",
            ))
    RD.write_tsv(os.path.join(V2, "control_availability_v2.tsv"), list(av[0].keys()), av)

    # ---------------------------------------------------------- paired contrasts
    out, perprot = [], []
    for c in CANDS:
        for m in MULTS:
            # step 1-2: per position, average available controls, then difference
            bypos = {}
            cells_cand = [k for k in cand if k[2] == c and k[3] == m]
            n_cells_no_ctrl = 0
            for k in cells_cand:
                cs = ctrl.get(k, [])
                if not cs:
                    n_cells_no_ctrl += 1
                    continue
                d = {}
                for metric in METRICS:
                    cv = f(cand[k].get(metric))
                    kv = np.mean([f(r.get(metric)) for r in cs])
                    d[metric] = (cv, kv, cv - kv)
                d["_n_controls"] = len(cs)
                bypos[k] = d
            # step 3: average paired differences within protein
            prot = defaultdict(list)
            for k, d in bypos.items():
                prot[k[0]].append(d)
            for metric in METRICS:
                vals, cvs, kvs = [], [], []
                for acc, ds in prot.items():
                    vals.append(np.mean([x[metric][2] for x in ds]))
                    cvs.append(np.mean([x[metric][0] for x in ds]))
                    kvs.append(np.mean([x[metric][1] for x in ds]))
                if not vals:
                    continue
                v = np.array(vals)
                out.append(dict(
                    candidate=c, multiplier=m, metric=metric,
                    n_proteins=len(v), n_paired_positions=len(bypos),
                    n_candidate_cells=len(cells_cand),
                    n_candidate_cells_without_any_control=n_cells_no_ctrl,
                    mean_controls_per_position=f"{np.mean([d['_n_controls'] for d in bypos.values()]):.2f}"
                    if bypos else "",
                    mean_candidate=f"{np.mean(cvs):.17g}", mean_control=f"{np.mean(kvs):.17g}",
                    mean_paired_difference=f"{v.mean():.17g}",
                    median_paired_difference=f"{np.median(v):.17g}",
                    sd_paired_difference=f"{v.std(ddof=1) if len(v) > 1 else 0.0:.17g}",
                    min_paired_difference=f"{v.min():.17g}",
                    p25=f"{np.percentile(v, 25):.17g}", p75=f"{np.percentile(v, 75):.17g}",
                    max_paired_difference=f"{v.max():.17g}",
                    n_proteins_positive=int((v > 0).sum()),
                    n_proteins_negative=int((v < 0).sum()),
                    n_proteins_zero=int((v == 0).sum()),
                ))
                if metric == "kl_aa20_renormalised_nats":
                    for acc, ds in prot.items():
                        perprot.append(dict(
                            ACC=acc,
                            localisation_class="Extracellular" if targets[acc] == 1 else "Cytoplasm",
                            candidate=c, multiplier=m, metric=metric,
                            n_paired_positions=len(ds),
                            mean_candidate=f"{np.mean([x[metric][0] for x in ds]):.17g}",
                            mean_control=f"{np.mean([x[metric][1] for x in ds]):.17g}",
                            paired_difference=f"{np.mean([x[metric][2] for x in ds]):.17g}",
                        ))
    RD.write_tsv(os.path.join(V2, "paired_contrasts_v2.tsv"), list(out[0].keys()), out)
    RD.write_tsv(os.path.join(V2, "paired_per_protein_v2.tsv"), list(perprot[0].keys()), perprot)

    # --------------------------------------------- top-1 denominators, by arm
    t1 = []
    for arm, src in (("candidate", cand.values()), ("control", [r for v in ctrl.values() for r in v])):
        rs = [r for r in src if float(r["multiplier"]) != 1.0]
        for c in CANDS + ([1446] if arm == "candidate" else []):
            sub = [r for r in rs if int(r["candidate"]) == c]
            if not sub:
                continue
            t1.append(dict(
                arm=arm, candidate=c, n_edit_conditions=len(sub),
                n_top1_aa20_changed=sum(1 for r in sub if not r["top1_aa20_same"]),
                n_top1_full_changed=sum(1 for r in sub if not r["top1_full_same"]),
                pct_top1_aa20_changed=f"{100.0*sum(1 for r in sub if not r['top1_aa20_same'])/len(sub):.2f}",
                note="denominator is edit conditions of THIS arm only; candidate and control "
                     "edits are never pooled",
            ))
    RD.write_tsv(os.path.join(V2, "top1_denominators_by_arm.tsv"), list(t1[0].keys()), t1)

    summary = dict(
        corrected="control arm rebuilt with residual-space norm matching",
        candidate_arm="reused unchanged from v1",
        feature_1446="descriptive only, no eligible matched control",
        availability=av, top1=t1,
        output_hashes={fn: RD.sha256_file(os.path.join(V2, fn)) for fn in (
            "control_availability_v2.tsv", "paired_contrasts_v2.tsv",
            "paired_per_protein_v2.tsv", "top1_denominators_by_arm.tsv")},
    )
    with open(os.path.join(V2, "analysis_summary_v2.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(dict(availability=av, top1=t1), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
