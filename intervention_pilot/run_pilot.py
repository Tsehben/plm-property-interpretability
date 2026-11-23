#!/usr/bin/env python3
"""
Bounded P3 intervention pilot — exploratory tuning, not evaluation.

Fixed candidates 1446, 2564, 3900. Validated direct decoder-delta edit:
    delta_h = ((topK(a_edited) - topK(a)) @ w_dec) * std ,  h' = h + delta_h
applied by forward hook on esm.encoder.layer[7], whose output is hidden_states[8].

P4 outcomes and P0 are never read. No sequence generation, no external predictor,
no feature screening, nothing committed. Existing results and the hashed candidate
list are read-only.

Stages:
  1  freeze the 24-protein P3 sample and the configuration, before any outcome
  2  select controls from P1+P2 discovery statistics only
  3  run the positional design and record full-precision outcomes
  4  aggregate

Usage:  python3 run_pilot.py [--stage {sample,run,all}] [--limit N]
"""

import argparse
import csv
import hashlib
import json
import os
import platform
import sys
import time
from collections import Counter, defaultdict

import numpy as np
import torch
from transformers import AutoTokenizer, EsmForMaskedLM

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "v1")
PER = os.path.join(OUT, "per_protein")
REF = os.path.join(ROOT, "interprot_ref", "raw")
sys.path.insert(0, REF)
sys.path.insert(0, os.path.join(ROOT, "candidate_discovery"))
from sae_model import SparseAutoencoder  # noqa: E402
import run_discovery as RD  # noqa: E402

INCLUSION = os.path.join(ROOT, "discovery_analysis", "v1", "inclusion_manifest.tsv")
RANKING = os.path.join(ROOT, "candidate_discovery", "v1", "eligible_feature_ranking.tsv")
SELECTED = os.path.join(ROOT, "candidate_discovery", "v1", "selected_feature_ids.json")
MANIFEST = os.path.join(ROOT, "DeepLoc2", "baseline", "deeploc_binary_manifest.tsv")
RAW_CSV = os.path.join(ROOT, "DeepLoc2", "raw", "Swissprot_Train_Validation_dataset.csv")

CANDIDATES = [1446, 2564, 3900]
MULTIPLIERS = [0.0, 0.5, 1.0, 1.5, 2.0]        # 1.0 is the no-op
NO_OP = 1.0
N_PER_CLASS = 12
LEN_MIN, LEN_MAX = 60, 1022
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
MASK_POSITIONS_FIXED = [10, 20, 30]            # one-based residue coordinates
N_CONTROLS = 5
MATCH_TOL = 0.20                               # prespecified +/-20% relative, both criteria
HOOK_LAYER = 7                                 # output == hidden_states[8]

# Kyte & Doolittle hydropathy index. Source: Kyte J, Doolittle RF, "A simple method
# for displaying the hydropathic character of a protein", J Mol Biol 157:105-132
# (1982), Table 1. Values are the published hydropathy indices, unscaled.
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

CONFIG = dict(
    purpose="exploratory P3 tuning pilot; not an evaluation and not a held-out result",
    partition="P3 (tuning) only; P4 and P0 never read",
    candidates=CANDIDATES,
    candidate_list_source="candidate_discovery/v1/selected_feature_ids.json (unchanged)",
    n_per_class=N_PER_CLASS,
    length_filter=f"{LEN_MIN} <= original_len <= {LEN_MAX}",
    residue_filter="standard 20 amino acids only (no B, U, Z, O, X)",
    exclusions="engineering proteins and p0_base_collision records, inherited from "
               "discovery_analysis/v1/inclusion_manifest.tsv (in_core == 1)",
    selection_rule="deterministic: within each class, sort by the frozen sequence SHA-256 "
                   "hex string ascending and take the first 12. No random seed is used and "
                   "no candidate activation or annotation value enters the ordering.",
    mask_positions="one-based residue positions 10, 20, 30 and floor(L/2); duplicates "
                   "collapsed and reported",
    edit="direct decoder-delta at the masked position only",
    multipliers=MULTIPLIERS,
    no_op_multiplier=NO_OP,
    inactive_policy="if the candidate's latent is 0 at the masked position, record an "
                    "inactive/no-edit condition; never insert the feature and never move "
                    "to a different position",
    control_rule=f"up to {N_CONTROLS} distinct control features per candidate from P1+P2 "
                 f"discovery statistics only, excluding all 20 selected candidates; "
                 f"activation frequency and mean active magnitude matched to +/-{MATCH_TOL:.0%} "
                 f"relative; NO AUROC restriction",
    control_edit="same masked input; latent delta matched in magnitude and sign to the "
                 "candidate's latent delta; non-negativity of latents respected; "
                 "unintended top-K support changes checked and reported",
    scope="local masked-token sensitivity at the edited position. This pilot does not test "
          "distant propagation, downstream layers beyond the forward pass, or biological "
          "relocation.",
    hydrophobicity="expected Kyte-Doolittle hydropathy over the AA-20 renormalised "
                   "distribution at the masked position; a biochemical readout, not a "
                   "localisation prediction. No directional hypothesis is imposed.",
)


def config_hash():
    return hashlib.sha256(json.dumps(CONFIG, sort_keys=True, default=str).encode()).hexdigest()


def device_of():
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


# --------------------------------------------------------------- stage 1: sample

def build_sample():
    incl = {r["ACC"]: r for r in RD.read_tsv(INCLUSION)}
    manifest = {r["ACC"]: r for r in RD.read_tsv(MANIFEST)}
    seqs = {}
    with open(RAW_CSV, newline="") as fh:
        for r in csv.DictReader(fh):
            seqs[r["ACC"]] = r["Sequence"]

    pool, excl = [], Counter()
    for acc, r in incl.items():
        if r["role"] != "tuning":
            excl["not_P3"] += 1
            continue
        if r["in_core"] != "1":
            excl[f"excluded_by_view:{r['exclusion_reason']}"] += 1
            continue
        m = manifest[acc]
        assert m["partition"] == "3", acc
        s = seqs[acc]
        assert hashlib.sha256(s.encode()).hexdigest() == m["sha256"], acc
        L = len(s)
        if not (LEN_MIN <= L <= LEN_MAX):
            excl["length_outside_60_1022"] += 1
            continue
        if set(s) - STANDARD_AA:
            excl["nonstandard_residues"] += 1
            continue
        pool.append(dict(ACC=acc, seq=s, target=int(m["target"]), original_len=L,
                         kingdom=m["kingdom"], sha256=m["sha256"]))

    chosen, shortfall = [], {}
    for t in (1, 0):
        cls = sorted([p for p in pool if p["target"] == t], key=lambda p: p["sha256"])
        take = cls[:N_PER_CLASS]
        chosen += take
        if len(take) < N_PER_CLASS:
            shortfall[str(t)] = dict(requested=N_PER_CLASS, available=len(cls), taken=len(take))
    return chosen, dict(
        pool_size=len(pool),
        pool_Extracellular=sum(1 for p in pool if p["target"] == 1),
        pool_Cytoplasm=sum(1 for p in pool if p["target"] == 0),
        exclusions=dict(excl),
        quota_shortfalls=shortfall or "none",
    )


def mask_positions(L):
    pos = MASK_POSITIONS_FIXED + [L // 2]
    uniq, dup = [], []
    for p in pos:
        if p in uniq:
            dup.append(p)
        elif 1 <= p <= L:
            uniq.append(p)
        else:
            dup.append(p)
    return uniq, dup


# ------------------------------------------------------------- stage 2: controls

def select_controls():
    sel_all = set(json.load(open(SELECTED))["selected_feature_ids"])
    rank = RD.read_tsv(RANKING)
    stats = {int(r["feature_id"]): dict(
        prevalence=float(r["prevalence_nterm"]),
        mam=float(r["mean_active_magnitude_nterm"]),
    ) for r in rank}
    out, avail = {}, []
    for c in CANDIDATES:
        cs = stats[c]
        pool = []
        for f, s in stats.items():
            if f in sel_all:
                continue
            if cs["prevalence"] == 0 or cs["mam"] == 0:
                continue
            dp = abs(s["prevalence"] - cs["prevalence"]) / cs["prevalence"]
            dm = abs(s["mam"] - cs["mam"]) / cs["mam"]
            if dp <= MATCH_TOL and dm <= MATCH_TOL:
                pool.append((dp + dm, f, dp, dm, s))
        pool.sort(key=lambda t: (t[0], t[1]))      # distance, then feature id
        take = pool[:N_CONTROLS]
        out[c] = [t[1] for t in take]
        for d, f, dp, dm, s in take:
            avail.append(dict(candidate=c, control=f, match_distance=f"{d:.6f}",
                              rel_diff_prevalence=f"{dp:.6f}", rel_diff_mean_active_magnitude=f"{dm:.6f}",
                              control_prevalence=f"{s['prevalence']:.6f}",
                              control_mean_active_magnitude=f"{s['mam']:.6f}",
                              candidate_prevalence=f"{cs['prevalence']:.6f}",
                              candidate_mean_active_magnitude=f"{cs['mam']:.6f}"))
        if len(take) < N_CONTROLS:
            avail.append(dict(candidate=c, control="", match_distance="",
                              rel_diff_prevalence="", rel_diff_mean_active_magnitude="",
                              control_prevalence="", control_mean_active_magnitude="",
                              candidate_prevalence="", candidate_mean_active_magnitude="",
                              ))
    shortfall = {str(c): dict(requested=N_CONTROLS, found=len(v)) for c, v in out.items()
                 if len(v) < N_CONTROLS}
    return out, avail, shortfall


# ------------------------------------------------------------------ stage 3: run

class Edit:
    """Forward hook implementing the validated direct decoder-delta edit at one
    position and one latent coordinate, with integrity bookkeeping."""

    def __init__(self, sae, feature=None, pos=None, new_value=None, capture_only=False):
        self.sae, self.feature, self.pos = sae, feature, pos
        self.new_value, self.capture_only = new_value, capture_only
        self.info = {}
        self.acts_at_pos = None

    def __call__(self, module, inp, out):
        is_tuple = isinstance(out, tuple)
        h = out[0] if is_tuple else out
        x = h[0].float().cpu()
        sae = self.sae
        a = sae.get_acts(x)
        self.acts_at_pos = a[self.pos].clone() if self.pos is not None else None
        if self.capture_only:
            self.info = dict(h_l2_norm=float(x.norm()))
            return out
        a_ed = a.clone()
        before = float(a[self.pos, self.feature])
        a_ed[self.pos, self.feature] = self.new_value
        delta_latent = self.new_value - before

        diff = (a_ed != a)
        sup_a = sae.topK_activation(a, sae.k) > 0
        sup_e = sae.topK_activation(a_ed, sae.k) > 0
        other = torch.ones(sae.d_hidden, dtype=torch.bool)
        other[self.feature] = False

        delta = ((sae.topK_activation(a_ed, sae.k) - sae.topK_activation(a, sae.k))
                 @ sae.w_dec) * sae.encode(x)[2]
        self.info = dict(
            latent_before=before, latent_after=float(self.new_value),
            latent_delta=delta_latent,
            n_latent_coords_changed=int(diff.sum()),
            support_evicted_unintended=int((sup_a & ~sup_e & other).sum()),
            support_inserted_unintended=int((~sup_a & sup_e & other).sum()),
            edited_feature_left_support=bool((sup_a & ~sup_e)[self.pos, self.feature]),
            delta_h_l2=float(delta.norm()),
            delta_h_max_abs=float(delta.abs().max()),
            h_l2_norm=float(x.norm()),
            delta_h_l2_over_h_l2=float(delta.norm() / x.norm()),
            nonzero_delta_positions=torch.unique(
                (delta.abs().sum(-1) > 0).nonzero().flatten()).tolist(),
        )
        new_h = h.clone()
        new_h[0] = (x + delta).to(h.dtype).to(h.device)
        return (new_h,) + out[1:] if is_tuple else new_h


def forward_with(mlm, ids, am, hook):
    handle = mlm.esm.encoder.layer[HOOK_LAYER].register_forward_hook(hook)
    try:
        with torch.inference_mode():
            out = mlm(input_ids=ids, attention_mask=am)
    finally:
        handle.remove()
    return out.logits.detach().clone()


def metrics(l_clean, l_edit, pos, aa_ids, true_id, kd_vec):
    la = l_clean[0, pos, :].cpu().double()
    lb = l_edit[0, pos, :].cpu().double()
    pa, pb = torch.softmax(la, -1), torch.softmax(lb, -1)
    qa, qb = torch.softmax(la[aa_ids], -1), torch.softmax(lb[aa_ids], -1)
    j = aa_ids.index(true_id) if true_id in aa_ids else None
    out = dict(
        kl_full_vocab_nats=float((pa * (pa.log() - pb.log())).sum()),
        kl_aa20_renormalised_nats=float((qa * (qa.log() - qb.log())).sum()),
        total_variation_aa20=float(0.5 * (qa - qb).abs().sum()),
        max_abs_logit_diff=float((la - lb).abs().max()),
        top1_full_same=int(la.argmax()) == int(lb.argmax()),
        top1_aa20_same=int(qa.argmax()) == int(qb.argmax()),
        top1_aa20_clean=int(aa_ids[int(qa.argmax())]),
        top1_aa20_edited=int(aa_ids[int(qb.argmax())]),
        nonstandard_mass_clean=float(1.0 - pa[aa_ids].sum()),
        nonstandard_mass_edited=float(1.0 - pb[aa_ids].sum()),
        expected_kd_clean=float((qa * kd_vec).sum()),
        expected_kd_edited=float((qb * kd_vec).sum()),
    )
    out["delta_expected_kd"] = out["expected_kd_edited"] - out["expected_kd_clean"]
    out["delta_nonstandard_mass"] = out["nonstandard_mass_edited"] - out["nonstandard_mass_clean"]
    if j is not None:
        out.update(
            p_true_clean=float(qa[j]), p_true_edited=float(qb[j]),
            delta_p_true=float(qb[j] - qa[j]),
            delta_logp_true=float(qb[j].log() - qa[j].log()),
            p_true_full_clean=float(pa[true_id]), p_true_full_edited=float(pb[true_id]),
        )
    return out


def run(sample, controls, dev, limit=None):
    os.makedirs(PER, exist_ok=True)
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

    aa_ids = [tok.convert_tokens_to_ids(a) for a in sorted(STANDARD_AA)]
    assert all(i is not None and i != tok.unk_token_id for i in aa_ids)
    kd_vec = torch.tensor([KD[a] for a in sorted(STANDARD_AA)], dtype=torch.float64)
    mask_id = tok.mask_token_id

    if limit:
        sample = sample[:limit]
    t0 = time.time()
    n_pass = 0
    for pi, prot in enumerate(sample):
        fp = os.path.join(PER, f"{prot['ACC']}.json")
        if os.path.exists(fp):
            print(f"  {prot['ACC']}: present, skipping (resume)", flush=True)
            continue
        s = prot["seq"]
        L = len(s)
        positions, dups = mask_positions(L)
        enc = tok([s], return_tensors="pt", add_special_tokens=True)
        ids0, am = enc["input_ids"], enc["attention_mask"]
        assert int(am.sum()) == L + 2
        rows = []
        tokchecks = []
        for p1 in positions:                       # one-based residue coordinate
            tpos = p1                              # <cls> at 0 => residue i at token i
            orig_tok = tok.convert_ids_to_tokens(int(ids0[0, tpos]))
            tokchecks.append(dict(residue_position_1based=p1, token_index=tpos,
                                  sequence_residue=s[p1 - 1], token_at_index=orig_tok,
                                  match=bool(orig_tok == s[p1 - 1])))
            assert orig_tok == s[p1 - 1], (p1, orig_tok, s[p1 - 1])
            true_id = int(ids0[0, tpos])
            ids = ids0.clone()
            ids[0, tpos] = mask_id
            idsd, amd = ids.to(dev), am.to(dev)

            cap = Edit(sae, pos=tpos, capture_only=True)
            l_clean = forward_with(mlm, idsd, amd, cap)
            n_pass += 1
            a_pos = cap.acts_at_pos                 # latents on the ACTUAL masked input

            for cand in CANDIDATES:
                a_c = float(a_pos[cand])
                for m in MULTIPLIERS:
                    delta_latent = a_c * m - a_c
                    base = dict(ACC=prot["ACC"], target=prot["target"],
                                localisation_class=RD.CLASS_NAME[prot["target"]]
                                if hasattr(RD, "CLASS_NAME") else ("Extracellular" if prot["target"] else "Cytoplasm"),
                                original_len=L, residue_position=p1, token_index=tpos,
                                masked_residue=s[p1 - 1], candidate=cand, multiplier=m,
                                candidate_latent_clean=f"{a_c:.17g}")
                    if a_c == 0.0:
                        rows.append(dict(base, arm="candidate", feature=cand,
                                         status="inactive_no_edit",
                                         target_latent_delta=f"{0.0:.17g}"))
                        for ctrl in controls[cand]:
                            rows.append(dict(base, arm="control", feature=ctrl,
                                             status="candidate_inactive_no_edit",
                                             target_latent_delta=f"{0.0:.17g}"))
                        continue
                    hk = Edit(sae, feature=cand, pos=tpos, new_value=a_c * m)
                    l_ed = forward_with(mlm, idsd, amd, hk)
                    n_pass += 1
                    rows.append(dict(base, arm="candidate", feature=cand, status="edited",
                                     target_latent_delta=f"{delta_latent:.17g}",
                                     **{k: (f"{v:.17g}" if isinstance(v, float) else v)
                                        for k, v in hk.info.items()
                                        if k != "nonzero_delta_positions"},
                                     n_nonzero_delta_positions=len(hk.info["nonzero_delta_positions"]),
                                     delta_confined_to_edited_position=int(
                                         hk.info["nonzero_delta_positions"] in ([], [tpos])),
                                     **{k: (f"{v:.17g}" if isinstance(v, float) else v)
                                        for k, v in metrics(l_clean, l_ed, tpos, aa_ids,
                                                            true_id, kd_vec).items()}))
                    for ctrl in controls[cand]:
                        a_k = float(a_pos[ctrl])
                        tgt = a_k + delta_latent
                        if a_k == 0.0:
                            st = "unavailable: control inactive at this position"
                        elif tgt < 0:
                            st = ("unavailable: matched negative delta would make the latent "
                                  "negative")
                        else:
                            st = "edited"
                        if st != "edited":
                            rows.append(dict(base, arm="control", feature=ctrl, status=st,
                                             control_latent_clean=f"{a_k:.17g}",
                                             target_latent_delta=f"{delta_latent:.17g}"))
                            continue
                        hkc = Edit(sae, feature=ctrl, pos=tpos, new_value=tgt)
                        l_c = forward_with(mlm, idsd, amd, hkc)
                        n_pass += 1
                        rows.append(dict(base, arm="control", feature=ctrl, status="edited",
                                         control_latent_clean=f"{a_k:.17g}",
                                         target_latent_delta=f"{delta_latent:.17g}",
                                         **{k: (f"{v:.17g}" if isinstance(v, float) else v)
                                            for k, v in hkc.info.items()
                                            if k != "nonzero_delta_positions"},
                                         n_nonzero_delta_positions=len(hkc.info["nonzero_delta_positions"]),
                                         delta_confined_to_edited_position=int(
                                             hkc.info["nonzero_delta_positions"] in ([], [tpos])),
                                         **{k: (f"{v:.17g}" if isinstance(v, float) else v)
                                            for k, v in metrics(l_clean, l_c, tpos, aa_ids,
                                                                true_id, kd_vec).items()}))
        with open(fp, "w") as fh:
            json.dump(dict(protein=prot["ACC"], target=prot["target"], original_len=L,
                           positions=positions, duplicate_or_invalid_positions=dups,
                           token_offset_checks=tokchecks, rows=rows), fh)
        el = time.time() - t0
        print(f"  [{pi+1}/{len(sample)}] {prot['ACC']} L={L} pos={positions} "
              f"rows={len(rows)} passes={n_pass} {el:.0f}s", flush=True)
    return dict(forward_passes=n_pass, seconds=time.time() - t0)


# ---------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["sample", "run", "all"], default="all")
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    sample, pool_info = build_sample()
    controls, avail, shortfall = select_controls()

    RD.write_tsv(os.path.join(OUT, "pilot_sample.tsv"),
                 ["ACC", "target", "localisation_class", "original_len", "kingdom", "sha256"],
                 [dict(ACC=p["ACC"], target=p["target"],
                       localisation_class="Extracellular" if p["target"] else "Cytoplasm",
                       original_len=p["original_len"], kingdom=p["kingdom"],
                       sha256=p["sha256"]) for p in sample])
    RD.write_tsv(os.path.join(OUT, "control_selection.tsv"), list(avail[0].keys()), avail)

    frozen = dict(
        config=CONFIG, config_hash=config_hash(),
        sample=dict(
            n=len(sample),
            n_Extracellular=sum(1 for p in sample if p["target"] == 1),
            n_Cytoplasm=sum(1 for p in sample if p["target"] == 0),
            accessions=[p["ACC"] for p in sample],
            sequence_sha256={p["ACC"]: p["sha256"] for p in sample},
            lengths={p["ACC"]: p["original_len"] for p in sample},
            **pool_info,
        ),
        controls=controls,
        control_shortfall=shortfall or "none",
        inputs=dict(
            inclusion_manifest=RD.sha256_file(INCLUSION),
            eligible_feature_ranking=RD.sha256_file(RANKING),
            selected_feature_ids=RD.sha256_file(SELECTED),
            deeploc_manifest=RD.sha256_file(MANIFEST),
            sae_file=RD.SAE_SHA,
        ),
        environment=dict(python=platform.python_version(), torch=torch.__version__,
                         transformers=__import__("transformers").__version__,
                         numpy=np.__version__, device=str(device_of()),
                         platform=platform.platform()),
    )
    with open(os.path.join(OUT, "frozen_configuration.json"), "w") as fh:
        json.dump(frozen, fh, indent=2)
    print(json.dumps(dict(config_hash=config_hash(), n=len(sample),
                          pool=pool_info, controls=controls,
                          control_shortfall=shortfall or "none"), indent=2), flush=True)
    if args.stage == "sample":
        return 0

    dev = device_of()
    rep = run(sample, controls, dev, limit=args.limit)
    with open(os.path.join(OUT, "run_report.json"), "w") as fh:
        json.dump(rep, fh, indent=2)
    print(json.dumps(rep, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
