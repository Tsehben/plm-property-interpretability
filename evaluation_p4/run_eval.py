#!/usr/bin/env python3
"""
Frozen P4 intervention evaluation — run once.

Everything is fixed before any P4 model output: configuration, code version,
control IDs and the 400-protein sample are hashed first. No P3 tuning, no new
candidates, no altered strengths, no sequence generation, no external predictor,
nothing committed. P0 remains reserved and is never read.

Fixed experiment
    positions   one-based residues 10, 20, 30, masked separately
    multipliers 0 (ablation), 1 (no-op), 2 (amplification)
    candidates  2564 and 3900 with their existing matched control IDs;
                1446 descriptive only (no eligible matched control exists)
    edit        validated direct decoder-delta, hook on esm.encoder.layer[7]
    controls    realised residual-space norm matching, relative tolerance 1e-5

Usage:  python3 run_eval.py
"""

import glob
import hashlib
import json
import os
import platform
import sys
import time
from collections import Counter

import numpy as np
import torch
from transformers import AutoTokenizer, EsmForMaskedLM

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "v1")
PER = os.path.join(OUT, "per_protein")
sys.path.insert(0, os.path.join(ROOT, "interprot_ref", "raw"))
sys.path.insert(0, os.path.join(ROOT, "candidate_discovery"))
sys.path.insert(0, os.path.join(ROOT, "intervention_pilot"))
from sae_model import SparseAutoencoder  # noqa: E402
import run_discovery as RD  # noqa: E402
import run_pilot as RP  # noqa: E402
from rerun_norm_matched_controls import Capture, NORM_TOL, TOL_JUSTIFICATION  # noqa: E402

PILOT_CFG = os.path.join(ROOT, "intervention_pilot", "v1", "frozen_configuration.json")

N_PER_CLASS = 200
POSITIONS = [10, 20, 30]
MULTIPLIERS = [0.0, 1.0, 2.0]
CONTRAST_MULTIPLIERS = [0.0, 2.0]
CANDIDATES_CONTROLLED = [2564, 3900]
CANDIDATE_DESCRIPTIVE = 1446
BOOTSTRAP_SEED = 20261122
BOOTSTRAP_N = 10000

CONFIG = dict(
    purpose="frozen P4 evaluation, run once at the frozen sample size",
    partition="P4 only; P0 reserved and never read",
    p4_history=(
        "P4 was used in the Check-in 1 analyses (probe fitting and reporting), so it is NOT "
        "a previously untouched final test set. It was excluded from the current candidate "
        "discovery (P1+P2) and from the intervention tuning pilot (P3). P0 remains reserved "
        "from fitting, selection, tuning and scoring."
    ),
    n_per_class=N_PER_CLASS,
    selection_rule=(
        "within each class, sort eligible P4 proteins by frozen sequence SHA-256 hex "
        "ascending and take the first 200. No random seed. No filtering by candidate "
        "activation, signal annotation or model outcome."
    ),
    eligibility="role == evaluation and in_core == 1 in the prepared inclusion manifest "
                "(engineering proteins and p0_base_collision records already excluded); "
                "original length 60-1022; standard 20 amino acids only",
    positions=POSITIONS,
    multipliers=MULTIPLIERS,
    candidates_controlled=CANDIDATES_CONTROLLED,
    candidate_descriptive=CANDIDATE_DESCRIPTIVE,
    edit="validated direct decoder-delta at the masked position only",
    control_norm_matching=dict(
        construction="delta_a_control = sign * T / (||w_dec[control]||_2 * std[p]), "
                     "T = candidate's realised ||delta_h|| on the same masked input and token",
        relative_tolerance=NORM_TOL,
        justification=TOL_JUSTIFICATION,
        prohibited="no clipping, no latent insertion, no silent control replacement, "
                   "no widening of the matching rule",
    ),
    primary_outcome=dict(
        name="paired KL difference",
        definition=(
            "for each usable protein/position/condition: average the AVAILABLE control KLs, "
            "subtract from the candidate KL, then average those paired position differences "
            "within protein; summarise across proteins"
        ),
        kl="KL(clean || edited), AA-20 renormalised, in nats",
        note="a neutral name. This is not a measure of localisation specificity.",
    ),
    primary_contrasts="4 = 2 candidates x 2 non-zero edit conditions (multipliers 0 and 2)",
    uncertainty=dict(
        method="protein-level resampling stratified by localisation class",
        n_resamples=BOOTSTRAP_N, seed=BOOTSTRAP_SEED,
        intervals="ordinary 95% BCa, plus Bonferroni-adjusted 98.75% BCa for the 4 primary "
                  "contrasts",
        degenerate_fallback=(
            "FIXED BEFORE RUNNING: if the BCa bias-correction z0 is non-finite (every "
            "resample on one side of the observed statistic) or the acceleration denominator "
            "is zero, fall back to the ordinary percentile interval at the same level and "
            "record interval_method = 'percentile_fallback' with the reason. A failed "
            "interval is never dropped and never left blank."
        ),
        units="proteins. Positions and control features are NOT independent sample units.",
    ),
    availability_policy=(
        "an inactive candidate yields a genuine no-edit candidate outcome; an unavailable "
        "matched-control comparison is MISSING and is never imputed as zero"
    ),
    design_note=(
        "the 200/200 design is deliberately balanced and does not represent natural cohort "
        "prevalence. Class-specific results are reported alongside the pooled "
        "sampled-population result. No post-result weighting is applied."
    ),
)


def config_hash():
    return hashlib.sha256(json.dumps(CONFIG, sort_keys=True, default=str).encode()).hexdigest()


def code_version():
    h = {}
    for p in [os.path.abspath(__file__),
              os.path.join(ROOT, "intervention_pilot", "run_pilot.py"),
              os.path.join(ROOT, "intervention_pilot", "rerun_norm_matched_controls.py"),
              os.path.join(ROOT, "candidate_discovery", "run_discovery.py"),
              os.path.join(ROOT, "interprot_ref", "raw", "sae_model.py")]:
        h[os.path.relpath(p, ROOT)] = RD.sha256_file(p)
    return h


def build_sample():
    incl = {r["ACC"]: r for r in RD.read_tsv(RP.INCLUSION)}
    manifest = {r["ACC"]: r for r in RD.read_tsv(RP.MANIFEST)}
    import csv
    seqs = {}
    with open(RP.RAW_CSV, newline="") as fh:
        for r in csv.DictReader(fh):
            seqs[r["ACC"]] = r["Sequence"]
    pool, excl = [], Counter()
    for acc, r in incl.items():
        if r["role"] != "evaluation":
            excl["not_P4"] += 1
            continue
        if r["in_core"] != "1":
            excl[f"excluded_by_view:{r['exclusion_reason']}"] += 1
            continue
        m = manifest[acc]
        assert m["partition"] == "4", acc
        s = seqs[acc]
        assert hashlib.sha256(s.encode()).hexdigest() == m["sha256"], acc
        if not (RP.LEN_MIN <= len(s) <= RP.LEN_MAX):
            excl["length_outside_60_1022"] += 1
            continue
        if set(s) - RP.STANDARD_AA:
            excl["nonstandard_residues"] += 1
            continue
        pool.append(dict(ACC=acc, seq=s, target=int(m["target"]), original_len=len(s),
                         kingdom=m["kingdom"], sha256=m["sha256"]))
    chosen = []
    quota = {}
    for t in (1, 0):
        cls = sorted([p for p in pool if p["target"] == t], key=lambda p: p["sha256"])
        quota[str(t)] = dict(available=len(cls), requested=N_PER_CLASS)
        chosen += cls[:N_PER_CLASS]
    for t, q in quota.items():
        assert q["available"] >= N_PER_CLASS, (
            f"QUOTA CANNOT BE MET for class {t}: {q}. Stopping rather than changing the rule.")
    assert len(chosen) == 2 * N_PER_CLASS
    return chosen, dict(pool_size=len(pool), quota=quota, exclusions=dict(excl))


def main():
    t0 = time.time()
    os.makedirs(PER, exist_ok=True)
    pilot = json.load(open(PILOT_CFG))
    controls = {int(k): v for k, v in pilot["controls"].items()}
    assert controls[1446] == [], "1446 must have no controls"

    sample, pool_info = build_sample()
    sample_hash = hashlib.sha256(
        "\n".join(f"{p['ACC']}\t{p['sha256']}\t{p['target']}" for p in sample).encode()
    ).hexdigest()

    frozen = dict(
        config=CONFIG, config_hash=config_hash(), code_version=code_version(),
        control_ids={str(k): v for k, v in controls.items()},
        control_ids_source="intervention_pilot/v1/frozen_configuration.json (unchanged)",
        control_ids_sha256=hashlib.sha256(
            json.dumps({str(k): v for k, v in controls.items()}, sort_keys=True).encode()).hexdigest(),
        sample=dict(n=len(sample), sample_sha256=sample_hash,
                    n_Extracellular=sum(1 for p in sample if p["target"] == 1),
                    n_Cytoplasm=sum(1 for p in sample if p["target"] == 0),
                    accessions=[p["ACC"] for p in sample],
                    sequence_sha256={p["ACC"]: p["sha256"] for p in sample},
                    **pool_info),
        inputs=dict(inclusion_manifest=RD.sha256_file(RP.INCLUSION),
                    deeploc_manifest=RD.sha256_file(RP.MANIFEST),
                    selected_feature_ids=RD.sha256_file(RP.SELECTED),
                    pilot_frozen_configuration=RD.sha256_file(PILOT_CFG),
                    sae_file=RD.SAE_SHA),
        environment=dict(python=platform.python_version(), torch=torch.__version__,
                         transformers=__import__("transformers").__version__,
                         numpy=np.__version__, device=str(RP.device_of()),
                         platform=platform.platform()),
        frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    )
    fc = os.path.join(OUT, "frozen_configuration.json")
    if os.path.exists(fc):
        old = json.load(open(fc))
        assert old["sample"]["sample_sha256"] == sample_hash, "sample drifted after freezing"
        assert old["config_hash"] == config_hash(), "configuration drifted after freezing"
        print("frozen configuration already present and matches", flush=True)
    else:
        with open(fc, "w") as fh:
            json.dump(frozen, fh, indent=2)
    RD.write_tsv(os.path.join(OUT, "eval_sample.tsv"),
                 ["ACC", "target", "localisation_class", "original_len", "kingdom", "sha256"],
                 [dict(ACC=p["ACC"], target=p["target"],
                       localisation_class="Extracellular" if p["target"] else "Cytoplasm",
                       original_len=p["original_len"], kingdom=p["kingdom"], sha256=p["sha256"])
                  for p in sample])
    print(f"FROZEN  config={config_hash()[:16]}  sample={sample_hash[:16]}  n={len(sample)} "
          f"({frozen['sample']['n_Extracellular']}/{frozen['sample']['n_Cytoplasm']})", flush=True)

    # ------------------------------------------------------------------ models
    dev = RP.device_of()
    tok = AutoTokenizer.from_pretrained(RD.CHECKPOINT, revision=RD.CHECKPOINT_REVISION)
    mlm = EsmForMaskedLM.from_pretrained(RD.CHECKPOINT, revision=RD.CHECKPOINT_REVISION,
                                         torch_dtype=torch.float32)
    mlm.eval().to(dev)
    for p in mlm.parameters():
        p.requires_grad_(False)
    from safetensors.torch import load_file
    assert RD.sha256_file(RD.SAE_FILE) == RD.SAE_SHA
    sae = SparseAutoencoder(RD.HIDDEN, RD.SAE_DIM)
    sae.load_state_dict(load_file(RD.SAE_FILE))
    sae.eval()
    dec_norm = sae.w_dec.detach().norm(dim=1)
    aa_ids = [tok.convert_tokens_to_ids(a) for a in sorted(RP.STANDARD_AA)]
    kd_vec = torch.tensor([RP.KD[a] for a in sorted(RP.STANDARD_AA)], dtype=torch.float64)
    mask_id = tok.mask_token_id

    n_pass = 0
    for pi, prot in enumerate(sample):
        fp = os.path.join(PER, f"{prot['ACC']}.json")
        if os.path.exists(fp):
            continue
        s, L = prot["seq"], prot["original_len"]
        enc = tok([s], return_tensors="pt", add_special_tokens=True)
        ids0, am = enc["input_ids"], enc["attention_mask"]
        assert int(am.sum()) == L + 2
        rows, checks = [], []
        for p1 in POSITIONS:
            tpos = p1
            tokname = tok.convert_ids_to_tokens(int(ids0[0, tpos]))
            checks.append(dict(residue_position_1based=p1, token_index=tpos,
                               sequence_residue=s[p1 - 1], token_at_index=tokname,
                               match=bool(tokname == s[p1 - 1])))
            assert tokname == s[p1 - 1]
            true_id = int(ids0[0, tpos])
            ids = ids0.clone()
            ids[0, tpos] = mask_id
            idsd, amd = ids.to(dev), am.to(dev)
            cap = Capture(sae, tpos)
            l_clean = RP.forward_with(mlm, idsd, amd, cap)
            n_pass += 1
            a_pos, std_pos = cap.acts_at_pos, cap.std_at_pos

            for cand in [CANDIDATE_DESCRIPTIVE] + CANDIDATES_CONTROLLED:
                a_c = float(a_pos[cand])
                base = dict(ACC=prot["ACC"], target=prot["target"], original_len=L,
                            residue_position=p1, token_index=tpos, masked_residue=s[p1 - 1],
                            candidate=cand, candidate_latent_clean=f"{a_c:.17g}")
                for m in MULTIPLIERS:
                    b = dict(base, multiplier=m)
                    if a_c == 0.0:
                        rows.append(dict(b, arm="candidate", feature=cand,
                                         status="inactive_no_edit"))
                        for ctrl in controls[cand]:
                            rows.append(dict(b, arm="control", feature=ctrl,
                                             status="candidate_inactive_no_edit"))
                        continue
                    hk = RP.Edit(sae, feature=cand, pos=tpos, new_value=a_c * m)
                    l_ed = RP.forward_with(mlm, idsd, amd, hk)
                    n_pass += 1
                    rows.append(dict(
                        b, arm="candidate", feature=cand, status="edited",
                        **{k: (f"{v:.17g}" if isinstance(v, float) else v)
                           for k, v in hk.info.items() if k != "nonzero_delta_positions"},
                        n_nonzero_delta_positions=len(hk.info["nonzero_delta_positions"]),
                        delta_confined_to_edited_position=int(
                            hk.info["nonzero_delta_positions"] in ([], [tpos])),
                        **{k: (f"{v:.17g}" if isinstance(v, float) else v)
                           for k, v in RP.metrics(l_clean, l_ed, tpos, aa_ids,
                                                  true_id, kd_vec).items()}))
                    T = hk.info["delta_h_l2"]
                    sgn = np.sign(hk.info["latent_delta"])
                    for ctrl in controls[cand]:
                        a_k = float(a_pos[ctrl])
                        req = 0.0 if T == 0.0 else sgn * T / (float(dec_norm[ctrl]) * std_pos)
                        tgt = a_k + req
                        if a_k == 0.0 and T != 0.0:
                            st = "unavailable: control inactive at this position"
                        elif tgt < 0:
                            st = ("unavailable: required negative latent change exceeds the "
                                  "control's latent value")
                        else:
                            st = "edited"
                        if st != "edited":
                            rows.append(dict(b, arm="control", feature=ctrl, status=st,
                                             control_latent_clean=f"{a_k:.17g}",
                                             target_norm=f"{T:.17g}",
                                             required_latent_delta=f"{req:.17g}"))
                            continue
                        hc = RP.Edit(sae, feature=ctrl, pos=tpos, new_value=tgt)
                        l_c = RP.forward_with(mlm, idsd, amd, hc)
                        n_pass += 1
                        realised = hc.info["delta_h_l2"]
                        rel = abs(realised - T) / T if T > 0 else 0.0
                        unint = (hc.info["support_evicted_unintended"]
                                 + hc.info["support_inserted_unintended"])
                        common = dict(b, arm="control", feature=ctrl,
                                      control_latent_clean=f"{a_k:.17g}",
                                      target_norm=f"{T:.17g}",
                                      required_latent_delta=f"{req:.17g}",
                                      realised_norm=f"{realised:.17g}",
                                      rel_norm_error=f"{rel:.17g}")
                        if unint:
                            rows.append(dict(common, status="unavailable: top-K support changed"))
                            continue
                        if rel > NORM_TOL:
                            rows.append(dict(common,
                                             status="unavailable: realised norm outside tolerance"))
                            continue
                        rows.append(dict(
                            common, status="edited",
                            **{k: (f"{v:.17g}" if isinstance(v, float) else v)
                               for k, v in hc.info.items() if k != "nonzero_delta_positions"},
                            n_nonzero_delta_positions=len(hc.info["nonzero_delta_positions"]),
                            delta_confined_to_edited_position=int(
                                hc.info["nonzero_delta_positions"] in ([], [tpos])),
                            **{k: (f"{v:.17g}" if isinstance(v, float) else v)
                               for k, v in RP.metrics(l_clean, l_c, tpos, aa_ids,
                                                      true_id, kd_vec).items()}))
        with open(fp, "w") as fh:
            json.dump(dict(protein=prot["ACC"], target=prot["target"], original_len=L,
                           positions=POSITIONS, token_offset_checks=checks, rows=rows), fh)
        if (pi + 1) % 10 == 0 or pi == 0:
            el = time.time() - t0
            print(f"  [{pi+1}/{len(sample)}] {prot['ACC']} passes={n_pass} {el:.0f}s "
                  f"({n_pass/max(el,1e-9):.1f} pass/s)", flush=True)

    with open(os.path.join(OUT, "run_report.json"), "w") as fh:
        json.dump(dict(forward_passes=n_pass, seconds=time.time() - t0,
                       n_proteins=len(sample)), fh, indent=2)
    print(f"DONE  {n_pass} passes  {time.time()-t0:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
