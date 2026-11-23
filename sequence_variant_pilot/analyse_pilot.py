#!/usr/bin/env python3
"""
Analysis for the bounded P3 sequence-variant pilot.

Aggregation, applied in this order everywhere:
  1. replicate level: D = margin(variant) - margin(original) for one
     (protein, arm, replicate), valid only if BOTH the original and that variant
     produced a prediction;
  2. contrasts are formed WITHIN a replicate (the arms share that replicate's
     uniforms), then averaged over replicates within protein;
  3. protein level: average the valid replicates;
  4. across proteins: summarise, protein as the sampling unit.

Replicates, steps and arms are not independent sample units.

Every "typical protein" statement uses the MEDIAN OF THE PAIRED DIFFERENCES, never
a difference between separately computed medians. Both are printed side by side so
the distinction is checkable.

Exploratory P3 pilot evidence. No efficacy pass threshold, no confirmatory claim.
"""

import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "v1")
PER = os.path.join(OUT, "per_protein")
sys.path.insert(0, os.path.join(ROOT, "candidate_discovery"))
import run_discovery as RD  # noqa: E402

ARMS = ["A", "B", "C", "D", "E"]
ARM_LABEL = {"A": "A unmodified", "B": "B 1446 x2", "C": "C 1446 x0",
             "D": "D random +u (matched to B)", "E": "E random -u (matched to C)"}
CONTRASTS = [("B", "A"), ("C", "A"), ("B", "D"), ("C", "E")]


def f(v):
    return float(v)


def fields_of(rows):
    """Ordered union of keys, so a row with extra keys is never silently truncated."""
    out = []
    for r in rows:
        for k in r:
            if k not in out:
                out.append(k)
    return out


def dist(vals, prefix=""):
    v = np.asarray(vals, dtype=float)
    if v.size == 0:
        return {}
    return {
        f"{prefix}n": int(v.size),
        f"{prefix}mean": f"{v.mean():.17g}",
        f"{prefix}median": f"{np.median(v):.17g}",
        f"{prefix}sd": f"{v.std(ddof=1):.17g}" if v.size > 1 else "0",
        f"{prefix}p25": f"{np.percentile(v, 25):.17g}",
        f"{prefix}p75": f"{np.percentile(v, 75):.17g}",
        f"{prefix}min": f"{v.min():.17g}", f"{prefix}max": f"{v.max():.17g}",
        f"{prefix}n_positive": int((v > 0).sum()),
        f"{prefix}n_negative": int((v < 0).sum()),
        f"{prefix}n_exact_zero": int((v == 0).sum()),
    }


def main():
    frozen = json.load(open(os.path.join(OUT, "frozen_protocol.json")))
    scores = {r["query_id"]: r for r in RD.read_tsv(os.path.join(OUT, "wolfpsort_scores.tsv"))}
    overlap = {r["ACC"]: r for r in RD.read_tsv(os.path.join(OUT, "reference_overlap.tsv"))}

    per, cls = {}, {}
    for fp in sorted(glob.glob(os.path.join(PER, "*.json"))):
        d = json.load(open(fp))
        per[d["protein"]] = d
        cls[d["protein"]] = "Extracellular" if d["target"] == 1 else "Cytoplasm"
    accs = sorted(per)

    # ------------------------------------------------------- replicate-level D
    repl, missing = {}, []
    orig_margin = {}
    for a in accs:
        o = scores.get(f"{a}__orig")
        if o is None or o["status"] != "predicted":
            missing.append(dict(query_id=f"{a}__orig", reason=o["status"] if o else "absent"))
            continue
        orig_margin[a] = f(o["margin"])
        for v in per[a]["variants"]:
            s = scores.get(v["variant_id"])
            if s is None or s["status"] != "predicted":
                missing.append(dict(query_id=v["variant_id"],
                                    reason=s["status"] if s else "absent"))
                continue
            repl[(a, v["arm"], v["replicate"])] = dict(
                margin_variant=f(s["margin"]), margin_original=orig_margin[a],
                D=f(s["margin"]) - orig_margin[a],
                extr_variant=f(s["extr_score"]), cyto_variant=f(s["cyto_score"]),
                top_is_compound=int(s["top_is_compound"]),
                n_mutations=int(v["n_mutations"]),
                identity=f(v["sequence_identity"]),
                identical=int(v["identical_to_original"]),
                dkd_whole=f(v["delta_kd_whole"]), dkd_nterm=f(v["delta_kd_nterm30"]))

    # -------------------------------------------------- per-arm protein-level D
    arm_rows, perprot = [], []
    prot_D = defaultdict(dict)
    for arm in ARMS:
        pv, labs = [], []
        nrep = 0
        for a in accs:
            r = [repl[(a, arm, k)]["D"] for k in range(frozen["config"]["n_replicates"])
                 if (a, arm, k) in repl]
            if not r:
                continue
            nrep += len(r)
            m = float(np.mean(r))
            prot_D[a][arm] = m
            pv.append(m)
            labs.append(cls[a])
            perprot.append(dict(ACC=a, localisation_class=cls[a], arm=arm,
                                n_valid_replicates=len(r),
                                mean_D=f"{m:.17g}", median_D_within_protein=f"{np.median(r):.17g}",
                                exact_reference_match=overlap[a]["exact_reference_match"]))
        row = dict(arm=arm, label=ARM_LABEL[arm], n_valid_replicate_pairs=nrep,
                   **dist(pv))
        for c in ("Extracellular", "Cytoplasm"):
            sub = [x for x, l in zip(pv, labs) if l == c]
            row.update(dist(sub, prefix=f"{c}_"))
        arm_rows.append(row)
    RD.write_tsv(os.path.join(OUT, "arm_level_D.tsv"), fields_of(arm_rows), arm_rows)
    RD.write_tsv(os.path.join(OUT, "per_protein_D.tsv"), fields_of(perprot), perprot)

    # ---------------------------------------------------- replicate-paired contrasts
    con_rows, con_pp = [], []
    for hi, lo in CONTRASTS:
        pv, labs, ovl = [], [], []
        npairs = 0
        sep_med_hi, sep_med_lo = [], []
        for a in accs:
            d = []
            for k in range(frozen["config"]["n_replicates"]):
                x, y = repl.get((a, hi, k)), repl.get((a, lo, k))
                if x is None or y is None:
                    continue
                d.append(x["D"] - y["D"])
            if not d:
                continue
            npairs += len(d)
            m = float(np.mean(d))
            pv.append(m)
            labs.append(cls[a])
            ovl.append(int(overlap[a]["exact_reference_match"]))
            sep_med_hi.append(prot_D[a][hi])
            sep_med_lo.append(prot_D[a][lo])
            con_pp.append(dict(ACC=a, localisation_class=cls[a], contrast=f"{hi}-{lo}",
                               n_valid_replicate_pairs=len(d),
                               mean_paired_difference=f"{m:.17g}",
                               exact_reference_match=overlap[a]["exact_reference_match"]))
        row = dict(contrast=f"{hi}-{lo}",
                   description=f"{ARM_LABEL[hi]}  minus  {ARM_LABEL[lo]}",
                   n_valid_replicate_pairs=npairs, **dist(pv))
        row["median_OF_paired_differences"] = row.get("median", "")
        row["difference_OF_separate_medians_do_not_use"] = (
            f"{np.median(sep_med_hi) - np.median(sep_med_lo):.17g}" if pv else "")
        for c in ("Extracellular", "Cytoplasm"):
            row.update(dist([x for x, l in zip(pv, labs) if l == c], prefix=f"{c}_"))
        for tag, want in (("ref_match", 1), ("no_ref_match", 0)):
            row.update(dist([x for x, o in zip(pv, ovl) if o == want], prefix=f"{tag}_"))
        con_rows.append(row)
    RD.write_tsv(os.path.join(OUT, "contrasts.tsv"), fields_of(con_rows), con_rows)
    RD.write_tsv(os.path.join(OUT, "per_protein_contrasts.tsv"), fields_of(con_pp), con_pp)

    # ------------------------------------------------- sequence-level descriptives
    seq_rows = []
    for arm in ARMS:
        nm, idn, dkw, dkn, ident = [], [], [], [], []
        for a in accs:
            r = [repl[(a, arm, k)] for k in range(frozen["config"]["n_replicates"])
                 if (a, arm, k) in repl]
            if not r:
                continue
            nm.append(np.mean([x["n_mutations"] for x in r]))
            idn.append(np.mean([x["identity"] for x in r]))
            dkw.append(np.mean([x["dkd_whole"] for x in r]))
            dkn.append(np.mean([x["dkd_nterm"] for x in r]))
            ident.append(np.mean([x["identical"] for x in r]))
        seq_rows.append(dict(
            arm=arm, label=ARM_LABEL[arm],
            mean_n_mutations=f"{np.mean(nm):.6f}", median_n_mutations=f"{np.median(nm):.6f}",
            mean_sequence_identity=f"{np.mean(idn):.6f}",
            pct_variants_identical_to_original=f"{100*np.mean(ident):.2f}",
            mean_delta_kd_whole=f"{np.mean(dkw):.17g}",
            median_delta_kd_whole=f"{np.median(dkw):.17g}",
            mean_delta_kd_nterm30=f"{np.mean(dkn):.17g}",
            median_delta_kd_nterm30=f"{np.median(dkn):.17g}"))
    RD.write_tsv(os.path.join(OUT, "sequence_descriptives.tsv"), fields_of(seq_rows), seq_rows)

    # --------------------------------------- feature activity and perturbation norms
    act = defaultdict(list)
    norms = defaultdict(list)
    hyp = defaultdict(list)
    nsteps = defaultdict(int)
    for a in accs:
        for s in per[a]["steps"]:
            key = s["arm"]
            nsteps[key] += 1
            act[key].append(f(s["candidate_latent_clean"]) > 0)
            norms[key].append(f(s["delta_h_l2"]))
            if "hypothetical_candidate_norm_on_this_input" in s:
                hyp[key].append(f(s["hypothetical_candidate_norm_on_this_input"]))
    act_rows = []
    for arm in ARMS:
        nz = [x for x in norms[arm] if x > 0]
        act_rows.append(dict(
            arm=arm, label=ARM_LABEL[arm], n_steps=nsteps[arm],
            n_steps_candidate_active=int(sum(act[arm])),
            pct_steps_candidate_active=f"{100*np.mean(act[arm]):.2f}",
            n_steps_nonzero_perturbation=len(nz),
            mean_perturbation_norm_over_nonzero=f"{np.mean(nz):.17g}" if nz else "0",
            median_perturbation_norm_over_nonzero=f"{np.median(nz):.17g}" if nz else "0",
            max_perturbation_norm=f"{max(norms[arm]):.17g}"))
    # mean active steps per trajectory
    traj = defaultdict(list)
    for a in accs:
        c = defaultdict(int)
        for s in per[a]["steps"]:
            if f(s["candidate_latent_clean"]) > 0:
                c[(s["arm"], s["replicate"])] += 1
        for arm in ARMS:
            for k in range(frozen["config"]["n_replicates"]):
                traj[arm].append(c[(arm, k)])
    for row in act_rows:
        row["mean_active_steps_per_trajectory_of_10"] = f"{np.mean(traj[row['arm']]):.4f}"
    RD.write_tsv(os.path.join(OUT, "activity_and_norms.tsv"), fields_of(act_rows), act_rows)

    # ------------------------------------------------------------------ integrity
    conf = json.load(open(os.path.join(OUT, "generation_report.json")))
    wolf = json.load(open(os.path.join(OUT, "wolfpsort_summary.json")))
    noop = [s for a in accs for s in per[a]["steps"] if s["arm"] == "A"]
    allsteps = [s for a in accs for s in per[a]["steps"]]
    mut_ok = True
    for a in accs:
        for v in per[a]["variants"]:
            if v["mutated_positions"]:
                mut_ok &= all(int(p) in frozen["config"]["positions"]
                              for p in v["mutated_positions"].split(";"))
            mut_ok &= (v["length"] == per[a]["original_len"])
    nfinite = sum(1 for s in allsteps
                  if not (np.isfinite(f(s["delta_h_l2"])) and np.isfinite(f(s["prob_sampled"]))))
    # D/E realised vs hypothetical norms
    rel = [f(s["rel_norm_error"]) for s in allsteps if "rel_norm_error" in s]
    cos = [f(s["cosine_with_candidate_delta"]) for s in allsteps
           if "cosine_with_candidate_delta" in s and f(s["delta_h_l2"]) > 0]
    integ = dict(
        forward_passes=conf["forward_passes"],
        forward_passes_expected=conf["expected_passes"],
        accounting_complete=conf["forward_passes"] == conf["expected_passes"],
        n_steps_recorded=len(allsteps),
        n_steps_expected=frozen["sample"]["n"] * frozen["config"]["n_replicates"]
        * len(frozen["config"]["arms"]) * len(frozen["config"]["positions"]),
        arm_A_noop_max_delta_h=f"{max(f(s['delta_h_l2']) for s in noop):.17g}",
        arm_A_noop_max_latent_coords_changed=max(int(s["n_latent_coords_changed"])
                                                 for s in noop),
        max_unintended_topk_changes=max(int(s["support_evicted_unintended"])
                                        + int(s["support_inserted_unintended"])
                                        for s in allsteps),
        n_steps_delta_not_confined=sum(1 for s in allsteps
                                       if int(s["delta_confined_to_edited_position"]) != 1),
        n_latent_coords_changed_in_random_arms=max(
            int(s["n_latent_coords_changed"]) for s in allsteps if s["arm"] in ("D", "E")),
        n_non_finite=nfinite,
        all_lengths_constant_and_mutations_in_allowed_positions=bool(mut_ok),
        n_sequences_scored=wolf["n_queries"], n_missing_predictions=wolf["n_missing"],
        n_unparsable=wolf["n_unparsable"],
        control_norm_rel_error_max=f"{max(rel):.3e}" if rel else "n/a",
        random_direction_cosine_with_candidate_delta=dict(
            n=len(cos), mean=f"{np.mean(cos):.6f}" if cos else "n/a",
            max_abs=f"{max(abs(x) for x in cos):.6f}" if cos else "n/a"),
        n_valid_replicate_cells=len(repl),
        n_valid_replicate_cells_expected=frozen["sample"]["n"]
        * frozen["config"]["n_replicates"] * len(frozen["config"]["arms"]),
        missing_pairs=missing,
        generation_seconds=round(conf["seconds"], 1),
        wolfpsort_seconds=wolf["wolfpsort_seconds"],
    )

    # ------------------------------------------------------------ overlap detail
    ov_rows = [r for r in RD.read_tsv(os.path.join(OUT, "reference_overlap.tsv"))
               if r["exact_reference_match"] == "1"]

    out = dict(integrity=integ, arms=arm_rows, contrasts=con_rows,
               sequence_descriptives=seq_rows, activity_and_norms=act_rows,
               reference_overlap_matches=ov_rows,
               n_exact_reference_matches=len(ov_rows),
               omitted_class_verdict=wolf["omitted_class_verdict"],
               output_hashes={fn: RD.sha256_file(os.path.join(OUT, fn)) for fn in (
                   "arm_level_D.tsv", "contrasts.tsv", "per_protein_D.tsv",
                   "per_protein_contrasts.tsv", "sequence_descriptives.tsv",
                   "activity_and_norms.tsv", "wolfpsort_scores.tsv")})
    with open(os.path.join(OUT, "analysis_summary.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps(dict(integrity=integ, contrasts=con_rows, arms=arm_rows,
                          activity=act_rows, sequence=seq_rows), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
