#!/usr/bin/env python3
"""
Correction to the P3 pilot: match the realised residual-space perturbation norm
||delta_h||, not the latent delta.

v1 is preserved untouched. Corrected outputs go to v2_normmatched/.

CONSTRUCTION. With the validated direct decoder-delta and a single edited latent
coordinate f at token position p,

    delta_h[p] = (topK(a_edited) - topK(a))[p] @ w_dec * std[p]
               = delta_a * w_dec[f] * std[p]          (support otherwise unchanged)

so  ||delta_h[p]||_2 = |delta_a| * ||w_dec[f]||_2 * std[p].

Given a candidate's REALISED target norm T measured on the same masked input and
the same token, the control's required latent change is

    delta_a_control = sign(delta_a_candidate) * T / ( ||w_dec[control]||_2 * std[p] )

evaluated with the same per-token std the SAE itself computes. The realised norm
is then recomputed from the actual edited forward pass and checked against T.

NUMERICAL MATCHING CRITERION, stated before the rerun:

    | ||delta_h_control|| - T | / T  <=  1e-5

Justification, from fp32 precision and the validated implementation only. The
intervention validation measured the direct decoder-delta against a float64
reference: median relative L2 error 3.367e-07, maximum 9.829e-07, minimum
significant digits 6.01 against fp32's ~7.2. A tolerance of 1e-5 sits about an
order of magnitude above the worst observed fp32 error for this operator, so it
admits genuine rounding while rejecting any real construction error. It is
derived from the arithmetic, not from any observed effect size.

CONSTRAINTS. Sign is preserved. Latents stay non-negative. Top-K support must be
unchanged apart from the edited coordinate. A control that cannot reach T under
those constraints is marked UNAVAILABLE: nothing is clipped, the matching rule is
not widened, and the operator is not switched.

Usage:  python3 rerun_norm_matched_controls.py
"""

import glob
import hashlib
import json
import os
import sys
import time
from collections import defaultdict

import numpy as np
import torch
from transformers import AutoTokenizer, EsmForMaskedLM

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
V1 = os.path.join(HERE, "v1")
OUT = os.path.join(HERE, "v2_normmatched")
PER = os.path.join(OUT, "per_protein")
sys.path.insert(0, os.path.join(ROOT, "interprot_ref", "raw"))
sys.path.insert(0, os.path.join(ROOT, "candidate_discovery"))
sys.path.insert(0, HERE)
from sae_model import SparseAutoencoder  # noqa: E402
import run_discovery as RD  # noqa: E402
import run_pilot as RP  # noqa: E402

NORM_TOL = 1e-5
TOL_JUSTIFICATION = (
    "Derived from fp32 precision and the validated edit implementation, not from any "
    "observed effect. The intervention validation measured the direct decoder-delta against "
    "a float64 reference at median relative L2 error 3.367e-07 and maximum 9.829e-07, with a "
    "minimum of 6.01 significant digits against fp32's ~7.2. 1e-5 is about an order of "
    "magnitude above that worst case: it admits fp32 rounding and rejects construction error."
)


def main():
    t0 = time.time()
    os.makedirs(PER, exist_ok=True)
    cfg = json.load(open(os.path.join(V1, "frozen_configuration.json")))

    # ---- the experiment is fixed: same proteins, masks, candidates, multipliers, controls
    reuse_ok = dict(
        config_hash_v1=cfg["config_hash"],
        inclusion_manifest=RD.sha256_file(RP.INCLUSION) == cfg["inputs"]["inclusion_manifest"],
        eligible_feature_ranking=RD.sha256_file(RP.RANKING) == cfg["inputs"]["eligible_feature_ranking"],
        selected_feature_ids=RD.sha256_file(RP.SELECTED) == cfg["inputs"]["selected_feature_ids"],
        deeploc_manifest=RD.sha256_file(RP.MANIFEST) == cfg["inputs"]["deeploc_manifest"],
        sae_file=RD.SAE_SHA == cfg["inputs"]["sae_file"],
    )
    assert all(v for k, v in reuse_ok.items() if k != "config_hash_v1"), reuse_ok
    controls = {int(k): v for k, v in cfg["controls"].items()}
    sample, _ = RP.build_sample()
    assert [p["ACC"] for p in sample] == cfg["sample"]["accessions"], "sample drifted"
    for p in sample:
        assert p["sha256"] == cfg["sample"]["sequence_sha256"][p["ACC"]]

    # ---- candidate outputs are REUSED unchanged (same inputs, same configuration)
    v1_rows = {}
    v1_meta = {}
    for fp in sorted(glob.glob(os.path.join(V1, "per_protein", "*.json"))):
        d = json.load(open(fp))
        v1_meta[d["protein"]] = d
        for r in d["rows"]:
            if r["arm"] == "candidate":
                v1_rows[(r["ACC"], r["residue_position"], int(r["candidate"]),
                         float(r["multiplier"]))] = r

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
    dec_norm = sae.w_dec.detach().norm(dim=1)          # ||w_dec[f]||_2 per feature

    aa_ids = [tok.convert_tokens_to_ids(a) for a in sorted(RP.STANDARD_AA)]
    kd_vec = torch.tensor([RP.KD[a] for a in sorted(RP.STANDARD_AA)], dtype=torch.float64)
    mask_id = tok.mask_token_id

    n_pass = 0
    for pi, prot in enumerate(sample):
        fp = os.path.join(PER, f"{prot['ACC']}.json")
        if os.path.exists(fp):
            print(f"  {prot['ACC']}: present, skipping", flush=True)
            continue
        s = prot["seq"]
        L = len(s)
        positions = v1_meta[prot["ACC"]]["positions"]
        enc = tok([s], return_tensors="pt", add_special_tokens=True)
        ids0, am = enc["input_ids"], enc["attention_mask"]
        rows = []
        for p1 in positions:
            tpos = p1
            assert tok.convert_ids_to_tokens(int(ids0[0, tpos])) == s[p1 - 1]
            true_id = int(ids0[0, tpos])
            ids = ids0.clone()
            ids[0, tpos] = mask_id
            idsd, amd = ids.to(dev), am.to(dev)

            cap = Capture(sae, tpos)
            l_clean = RP.forward_with(mlm, idsd, amd, cap)
            n_pass += 1
            a_pos, std_pos = cap.acts_at_pos, cap.std_at_pos

            for cand in (2564, 3900):
                for m in cfg["config"]["multipliers"]:
                    key = (prot["ACC"], p1, cand, float(m))
                    cr = v1_rows.get(key)
                    base = dict(ACC=prot["ACC"], target=prot["target"],
                                localisation_class="Extracellular" if prot["target"] else "Cytoplasm",
                                original_len=L, residue_position=p1, token_index=tpos,
                                candidate=cand, multiplier=m)
                    if cr is None or cr["status"] != "edited":
                        for ctrl in controls[cand]:
                            rows.append(dict(base, arm="control", feature=ctrl,
                                             status="candidate_inactive_no_edit"))
                        continue
                    T = float(cr["delta_h_l2"])
                    sgn = np.sign(float(cr["latent_delta"]))
                    for ctrl in controls[cand]:
                        a_k = float(a_pos[ctrl])
                        if T == 0.0:                        # the no-op
                            req = 0.0
                        else:
                            req = sgn * T / (float(dec_norm[ctrl]) * std_pos)
                        tgt = a_k + req
                        if a_k == 0.0 and T != 0.0:
                            st = "unavailable: control inactive at this position"
                        elif tgt < 0:
                            st = ("unavailable: required negative latent change exceeds the "
                                  "control's latent value")
                        else:
                            st = "edited"
                        if st != "edited":
                            rows.append(dict(base, arm="control", feature=ctrl, status=st,
                                             control_latent_clean=f"{a_k:.17g}",
                                             target_norm=f"{T:.17g}",
                                             required_latent_delta=f"{req:.17g}"))
                            continue
                        hk = RP.Edit(sae, feature=ctrl, pos=tpos, new_value=tgt)
                        l_c = RP.forward_with(mlm, idsd, amd, hk)
                        n_pass += 1
                        realised = hk.info["delta_h_l2"]
                        rel = abs(realised - T) / T if T > 0 else 0.0
                        unintended = (hk.info["support_evicted_unintended"]
                                      + hk.info["support_inserted_unintended"])
                        if unintended:
                            rows.append(dict(base, arm="control", feature=ctrl,
                                             status="unavailable: top-K support changed",
                                             control_latent_clean=f"{a_k:.17g}",
                                             target_norm=f"{T:.17g}",
                                             required_latent_delta=f"{req:.17g}",
                                             realised_norm=f"{realised:.17g}",
                                             rel_norm_error=f"{rel:.17g}"))
                            continue
                        if rel > NORM_TOL:
                            rows.append(dict(base, arm="control", feature=ctrl,
                                             status="unavailable: realised norm outside tolerance",
                                             control_latent_clean=f"{a_k:.17g}",
                                             target_norm=f"{T:.17g}",
                                             required_latent_delta=f"{req:.17g}",
                                             realised_norm=f"{realised:.17g}",
                                             rel_norm_error=f"{rel:.17g}"))
                            continue
                        rows.append(dict(
                            base, arm="control", feature=ctrl, status="edited",
                            control_latent_clean=f"{a_k:.17g}",
                            target_norm=f"{T:.17g}", required_latent_delta=f"{req:.17g}",
                            realised_norm=f"{realised:.17g}", rel_norm_error=f"{rel:.17g}",
                            **{k: (f"{v:.17g}" if isinstance(v, float) else v)
                               for k, v in hk.info.items() if k != "nonzero_delta_positions"},
                            n_nonzero_delta_positions=len(hk.info["nonzero_delta_positions"]),
                            delta_confined_to_edited_position=int(
                                hk.info["nonzero_delta_positions"] in ([], [tpos])),
                            **{k: (f"{v:.17g}" if isinstance(v, float) else v)
                               for k, v in RP.metrics(l_clean, l_c, tpos, aa_ids,
                                                      true_id, kd_vec).items()}))
        with open(fp, "w") as fh:
            json.dump(dict(protein=prot["ACC"], target=prot["target"], original_len=L,
                           positions=positions, rows=rows), fh)
        print(f"  [{pi+1}/{len(sample)}] {prot['ACC']} rows={len(rows)} passes={n_pass} "
              f"{time.time()-t0:.0f}s", flush=True)

    with open(os.path.join(OUT, "rerun_report.json"), "w") as fh:
        json.dump(dict(
            norm_tolerance=NORM_TOL, tolerance_justification=TOL_JUSTIFICATION,
            construction="delta_a_control = sign * T / (||w_dec[control]|| * std[p])",
            candidates_rerun=[2564, 3900],
            candidate_outputs="reused unchanged from v1 (inputs and configuration verified "
                              "identical by SHA-256)",
            feature_1446="descriptive only; no eligible matched control exists",
            reuse_checks=reuse_ok, forward_passes=n_pass, seconds=time.time() - t0,
            decoder_row_norms={str(f): float(dec_norm[f]) for f in
                               [2564, 3900] + controls[2564] + controls[3900]},
        ), fh, indent=2)
    print(f"done: {n_pass} passes, {time.time()-t0:.0f}s")
    return 0


class Capture:
    """Capture latents and the SAE's per-token std at one position."""

    def __init__(self, sae, pos):
        self.sae, self.pos = sae, pos
        self.acts_at_pos = None
        self.std_at_pos = None

    def __call__(self, module, inp, out):
        h = out[0] if isinstance(out, tuple) else out
        x = h[0].float().cpu()
        self.acts_at_pos = self.sae.get_acts(x)[self.pos].clone()
        self.std_at_pos = float(self.sae.encode(x)[2][self.pos])
        return out


if __name__ == "__main__":
    sys.exit(main())
