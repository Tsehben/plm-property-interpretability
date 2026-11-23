#!/usr/bin/env python3
"""
Bounded exploratory candidate discovery on the prepared localisation discovery
sample (P1+P2 core, n = 5,190).

Stage 1  freeze the configuration and assert the sample against the prepared
         inclusion manifest and the frozen cohort manifest;
Stage 2  extract layer-8 ESM-2 residue activations, encode with the InterProt
         layer-8 SAE, and pool per protein (streamed; no dense per-residue array
         is ever stored);
Stage 3  rank eligible features by |AUROC - 0.5| on the primary summary against
         the frozen DeepLoc binary label, select up to 20 distinct feature ids,
         write and hash the selected list.

Nothing in this script touches P0, P3 or P4. No interventions, no sequence
generation, no external predictor, no layer sweep. Existing data and results are
read-only; every artefact is written under candidate_discovery/v1/.

Usage:  python3 run_discovery.py [--limit N] [--stage {extract,rank,all}]
"""

import argparse
import csv
import hashlib
import json
import os
import platform
import sys
import time
from collections import Counter

import numpy as np
import torch
from transformers import AutoTokenizer, EsmModel

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "v1")
POOLED = os.path.join(OUT, "pooled")

REF = os.path.join(ROOT, "interprot_ref", "raw")
sys.path.insert(0, REF)
from sae_model import SparseAutoencoder  # noqa: E402

# ----------------------------------------------------------------- frozen config
CHECKPOINT = "facebook/esm2_t33_650M_UR50D"
CHECKPOINT_REVISION = "08e4846e537177426273712802403f7ba8261b6c"
SAE_FILE = os.path.join(REF, "esm2_plm1280_l8_sae4096.safetensors")
SAE_SHA = "1172503b90c950651e787a5625d91ff1c702a04ff3d851e8e95b700201e56ebd"
SAE_REPO = "liambai/InterProt-ESM2-SAEs"
SAE_REVISION = "3ded2e51641df84d6b89f25430874e3a8eb42c24"

LAYER = 8
MAX_RESIDUES = 1022            # unchanged screening policy
HIDDEN = 1280
SAE_DIM = 4096
TOPK = 128
NTERM = 30                     # primary pooling window, first 30 amino-acid positions
ELIGIBILITY_MIN_PREVALENCE = 0.05
N_CANDIDATES = 20
SHARD = 256
BATCH = 4

VIEW = os.path.join(ROOT, "discovery_analysis", "v1", "analysis_view_discovery.tsv")
MANIFEST = os.path.join(ROOT, "DeepLoc2", "baseline", "deeploc_binary_manifest.tsv")
RAW_CSV = os.path.join(ROOT, "DeepLoc2", "raw", "Swissprot_Train_Validation_dataset.csv")

CONFIG = dict(
    sample="discovery core P1+P2 (engineering and p0_base_collision excluded)",
    expected_n=5190,
    checkpoint=CHECKPOINT,
    checkpoint_revision=CHECKPOINT_REVISION,
    layer=LAYER,
    layer_indexing="HF hidden_states[8] = output of encoder.layer[7]",
    dtype="float32",
    mode="eval, inference_mode, frozen",
    sae_repo=SAE_REPO,
    sae_revision=SAE_REVISION,
    sae_file=os.path.basename(SAE_FILE),
    sae_sha256=SAE_SHA,
    sae_d_model=HIDDEN,
    sae_d_hidden=SAE_DIM,
    sae_k=TOPK,
    max_residues=MAX_RESIDUES,
    truncation_policy="first 1022 residues, unchanged from the existing screening policy",
    token_policy="<cls> + L residue tokens + <eos>; pooling uses hidden[1:1+L] only; "
                 "special tokens and padding excluded",
    primary_summary="mean_nterm: mean SAE activation over the first 30 amino-acid positions",
    secondary_summary="mean_visible: mean SAE activation over all visible amino-acid positions",
    eligibility="feature active (>0) at least once within the first 30 positions in at "
                "least 5% of discovery proteins; the SAE's training-time dead list is NOT used",
    ranking_statistic="abs(AUROC - 0.5) of the primary summary against the frozen DeepLoc "
                      "binary label; ties broken by ascending feature id",
    n_candidates=N_CANDIDATES,
    label_source="frozen DeepLoc binary target in deeploc_binary_manifest.tsv",
    reserved_partitions_excluded=["0", "3", "4"],
)


def config_hash():
    return hashlib.sha256(json.dumps(CONFIG, sort_keys=True).encode()).hexdigest()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, fields, records):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in records:
            w.writerow(r)


def device_of():
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


# -------------------------------------------------------------- stage 1: sample

def load_sample():
    """Assert exact accessions, labels, partitions and sequence hashes."""
    view = read_tsv(VIEW)
    manifest = {r["ACC"]: r for r in read_tsv(MANIFEST)}

    # frozen extraction-input sequences, straight from the raw DeepLoc file
    seqs = {}
    with open(RAW_CSV, newline="") as fh:
        for r in csv.DictReader(fh):
            seqs[r["ACC"]] = r["Sequence"]

    checks = dict(n_rows=len(view), expected_n=CONFIG["expected_n"])
    assert len(view) == CONFIG["expected_n"], (len(view), CONFIG["expected_n"])
    assert len(set(r["ACC"] for r in view)) == len(view), "duplicate accession in the view"

    rows = []
    bad_hash, missing_seq, bad_label, bad_partition = [], [], [], []
    for r in view:
        acc = r["ACC"]
        m = manifest.get(acc)
        if m is None:
            missing_seq.append(acc)
            continue
        if m["partition"] not in ("1", "2"):
            bad_partition.append(acc)
        if m["target"] != r["target"]:
            bad_label.append(acc)
        s = seqs.get(acc)
        if s is None:
            missing_seq.append(acc)
            continue
        if hashlib.sha256(s.encode()).hexdigest() != m["sha256"]:
            bad_hash.append(acc)
            continue
        rows.append(
            dict(
                ACC=acc,
                seq=s,
                target=int(m["target"]),
                partition=m["partition"],
                kingdom=m["kingdom"],
                original_len=len(s),
                visible_len=min(len(s), MAX_RESIDUES),
                truncated=int(len(s) > MAX_RESIDUES),
            )
        )
    checks.update(
        sequence_hash_mismatches=bad_hash,
        missing_sequences=missing_seq,
        label_mismatches=bad_label,
        partition_violations=bad_partition,
        n_usable=len(rows),
    )
    assert not bad_hash and not missing_seq and not bad_label and not bad_partition, checks
    assert len(rows) == CONFIG["expected_n"]
    # No reserved partition may appear.
    assert set(r["partition"] for r in rows) == {"1", "2"}
    return rows, checks


# ------------------------------------------------------- stage 2: extract + pool

def esm_residue_acts(model, tok, seqs, device):
    enc = tok(seqs, padding=True, return_tensors="pt")
    ids, am = enc["input_ids"].to(device), enc["attention_mask"].to(device)
    assert (ids == tok.unk_token_id).sum().item() == 0, "unexpected <unk> token"
    for r, s in enumerate(seqs):
        assert int(am[r].sum()) == len(s) + 2, f"token/residue mismatch len {len(s)}"
    with torch.inference_mode():
        h = model(input_ids=ids, attention_mask=am, output_hidden_states=True).hidden_states[LAYER]
    return [h[r, 1 : 1 + len(s), :].float().cpu() for r, s in enumerate(seqs)], h


def reference_agreement(sae, model, tok, rows, device, n=4):
    """The authors' own path (SAE applied to the full (L+2) tensor, then sliced)
    must agree with ours (slice first, then SAE), to the last bit we can check."""
    sub = rows[:n]
    seqs = [r["seq"][:MAX_RESIDUES] for r in sub]
    acts, hfull = esm_residue_acts(model, tok, seqs, device)
    worst = 0.0
    for j, s in enumerate(seqs):
        full = hfull[j, : len(s) + 2, :].float().cpu()
        ref = sae.get_acts(full).numpy()[1 : 1 + len(s)]
        mine = sae.get_acts(acts[j]).numpy()
        worst = max(worst, float(np.abs(ref - mine).max()))
    return dict(n_sequences=len(sub), max_abs_difference=worst,
                accessions=[r["ACC"] for r in sub])


def extract(rows, device, limit=None):
    """Streams. Per protein we keep four 4096-vectors and nothing residue-level."""
    os.makedirs(POOLED, exist_ok=True)
    if limit:
        rows = rows[:limit]

    tok = AutoTokenizer.from_pretrained(CHECKPOINT, revision=CHECKPOINT_REVISION)
    model = EsmModel.from_pretrained(
        CHECKPOINT, revision=CHECKPOINT_REVISION, add_pooling_layer=False,
        torch_dtype=torch.float32,
    )
    model.eval().to(device)
    for p in model.parameters():
        p.requires_grad_(False)

    from safetensors.torch import load_file
    assert sha256_file(SAE_FILE) == SAE_SHA, "SAE checkpoint hash mismatch"
    sae = SparseAutoencoder(HIDDEN, SAE_DIM)
    sae.load_state_dict(load_file(SAE_FILE))
    sae.eval()
    assert sae.k == TOPK

    ref = reference_agreement(sae, model, tok, rows, device)
    print(f"reference SAE agreement: max|diff| = {ref['max_abs_difference']:.3e}", flush=True)

    # length-sorted batching; keys carried through
    order = sorted(range(len(rows)), key=lambda i: -rows[i]["visible_len"])
    shards = [order[i : i + SHARD] for i in range(0, len(order), SHARD)]

    t0 = time.time()
    n_done = n_res = 0
    nonfinite = 0
    for si, shard in enumerate(shards):
        fp = os.path.join(POOLED, f"pooled_shard{si:04d}.npz")
        if os.path.exists(fp):
            z = np.load(fp, allow_pickle=True)
            n_done += len(z["keys"])
            print(f"  shard {si:04d}: present, resuming past it", flush=True)
            continue
        keys = []
        mean_n = np.zeros((len(shard), SAE_DIM), np.float32)
        mean_v = np.zeros((len(shard), SAE_DIM), np.float32)
        nact_n = np.zeros((len(shard), SAE_DIM), np.uint16)
        nact_v = np.zeros((len(shard), SAE_DIM), np.uint16)
        w = 0
        for bstart in range(0, len(shard), BATCH):
            idxs = shard[bstart : bstart + BATCH]
            seqs = [rows[i]["seq"][:MAX_RESIDUES] for i in idxs]
            acts, _ = esm_residue_acts(model, tok, seqs, device)
            for j, i in enumerate(idxs):
                L = len(seqs[j])
                assert acts[j].shape == (L, HIDDEN)
                lat = sae.get_acts(acts[j]).numpy()        # (L, 4096) float32
                assert lat.shape == (L, SAE_DIM)
                if not np.isfinite(lat).all():
                    nonfinite += 1
                pos = lat > 0
                nt = min(NTERM, L)
                keys.append(rows[i]["ACC"])
                mean_n[w] = lat[:nt].mean(axis=0)
                mean_v[w] = lat.mean(axis=0)
                nact_n[w] = pos[:nt].sum(axis=0).astype(np.uint16)
                nact_v[w] = np.minimum(pos.sum(axis=0), 65535).astype(np.uint16)
                w += 1
                n_done += 1
                n_res += L
                del lat, pos
        np.savez_compressed(
            fp, keys=np.array(keys), mean_nterm=mean_n[:w], mean_visible=mean_v[:w],
            nact_nterm=nact_n[:w], nact_visible=nact_v[:w],
        )
        el = time.time() - t0
        print(f"  shard {si:04d}: {w} proteins  |  {n_done}/{len(rows)}  "
              f"{el:.0f}s  {n_done / max(el, 1e-9):.2f} seq/s", flush=True)

    return dict(
        reference_agreement=ref, n_proteins=n_done, n_visible_residues=n_res,
        nonfinite_proteins=nonfinite, seconds=time.time() - t0,
    )


def load_pooled(rows):
    import glob
    ks, mn, mv, an, av = [], [], [], [], []
    for f in sorted(glob.glob(os.path.join(POOLED, "pooled_shard*.npz"))):
        z = np.load(f, allow_pickle=True)
        ks.append(z["keys"]); mn.append(z["mean_nterm"]); mv.append(z["mean_visible"])
        an.append(z["nact_nterm"]); av.append(z["nact_visible"])
    keys = np.concatenate(ks)
    order = {str(k): i for i, k in enumerate(keys)}
    assert len(order) == len(keys), "duplicate key across shards"
    want = [r["ACC"] for r in rows]
    assert set(order) == set(want), "key alignment failure between pooled shards and the sample"
    idx = np.array([order[a] for a in want])
    return dict(
        keys=np.array(want),
        mean_nterm=np.concatenate(mn)[idx],
        mean_visible=np.concatenate(mv)[idx],
        nact_nterm=np.concatenate(an)[idx],
        nact_visible=np.concatenate(av)[idx],
    )


# ------------------------------------------------------------- stage 3: ranking

def auroc_columns(X, y):
    """Rank-based AUROC per column with proper average ranks for ties.
    AUROC > 0.5 means higher values associate with y == 1 (Extracellular)."""
    from scipy.stats import rankdata
    n1 = int(y.sum())
    n0 = int(len(y) - n1)
    R = rankdata(X, axis=0)                      # average ranks, ties handled
    s1 = R[y == 1].sum(axis=0)
    return (s1 - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def rank_features(P, rows):
    y = np.array([r["target"] for r in rows], dtype=np.int64)
    part = np.array([r["partition"] for r in rows])
    X = P["mean_nterm"]
    N = len(rows)

    active = P["nact_nterm"] > 0                 # per protein, any active position in window
    prevalence = active.mean(axis=0)
    eligible = prevalence >= ELIGIBILITY_MIN_PREVALENCE

    prev_p1 = active[part == "1"].mean(axis=0)
    prev_p2 = active[part == "2"].mean(axis=0)

    # mean activation magnitude when active, for later control matching
    with np.errstate(invalid="ignore"):
        mam = np.where(active.sum(axis=0) > 0,
                       (X * active).sum(axis=0) / np.maximum(active.sum(axis=0), 1),
                       0.0)
    mean_positions_active = P["nact_nterm"].astype(np.float64).mean(axis=0)
    prevalence_visible = (P["nact_visible"] > 0).mean(axis=0)

    auroc_all = auroc_columns(X, y)
    auroc_p1 = auroc_columns(X[part == "1"], y[part == "1"])
    auroc_p2 = auroc_columns(X[part == "2"], y[part == "2"])
    auroc_sec = auroc_columns(P["mean_visible"], y)

    feat = np.arange(SAE_DIM)
    stat = np.abs(auroc_all - 0.5)
    # deterministic: descending statistic, ascending feature id
    order = sorted(feat[eligible], key=lambda f: (-stat[f], f))
    rank_all = {f: i + 1 for i, f in enumerate(order)}
    o1 = sorted(feat[eligible], key=lambda f: (-abs(auroc_p1[f] - 0.5), f))
    o2 = sorted(feat[eligible], key=lambda f: (-abs(auroc_p2[f] - 0.5), f))
    rank_p1 = {f: i + 1 for i, f in enumerate(o1)}
    rank_p2 = {f: i + 1 for i, f in enumerate(o2)}

    records = []
    for f in order:
        records.append(
            dict(
                rank=rank_all[f], feature_id=int(f),
                auroc=f"{auroc_all[f]:.17g}",
                abs_auroc_minus_half=f"{stat[f]:.17g}",
                direction="higher_in_Extracellular" if auroc_all[f] > 0.5 else "higher_in_Cytoplasm",
                auroc_P1=f"{auroc_p1[f]:.17g}", auroc_P2=f"{auroc_p2[f]:.17g}",
                rank_P1=rank_p1[f], rank_P2=rank_p2[f],
                auroc_secondary_mean_visible=f"{auroc_sec[f]:.17g}",
                prevalence_nterm=f"{prevalence[f]:.17g}",
                prevalence_nterm_P1=f"{prev_p1[f]:.17g}",
                prevalence_nterm_P2=f"{prev_p2[f]:.17g}",
                prevalence_visible=f"{prevalence_visible[f]:.17g}",
                mean_active_magnitude_nterm=f"{mam[f]:.17g}",
                mean_positions_active_nterm=f"{mean_positions_active[f]:.17g}",
                mean_summary_all=f"{X[:, f].mean():.17g}",
                mean_summary_Extracellular=f"{X[y == 1, f].mean():.17g}",
                mean_summary_Cytoplasm=f"{X[y == 0, f].mean():.17g}",
            )
        )

    selected = [int(f) for f in order[:N_CANDIDATES]]
    assert len(set(selected)) == len(selected), "selected feature ids are not distinct"

    # stability descriptives over the eligible set
    from scipy.stats import spearmanr
    el = feat[eligible]
    sp_signed = spearmanr(auroc_p1[el], auroc_p2[el]).statistic
    sp_abs = spearmanr(np.abs(auroc_p1[el] - 0.5), np.abs(auroc_p2[el] - 0.5)).statistic
    top20_p1 = set(o1[:N_CANDIDATES])
    top20_p2 = set(o2[:N_CANDIDATES])
    stability = dict(
        definition_top20_overlap=(
            "|top-20 by |AUROC-0.5| computed on P1 alone  INTERSECT  top-20 computed on P2 "
            "alone|, both restricted to the eligible set defined on P1+P2 together"
        ),
        top20_P1_vs_P2_overlap=len(top20_p1 & top20_p2),
        selected_in_P1_top20=len(set(selected) & top20_p1),
        selected_in_P2_top20=len(set(selected) & top20_p2),
        selected_in_both_P1_and_P2_top20=len(set(selected) & top20_p1 & top20_p2),
        definition_rank_correlation=(
            "Spearman rank correlation across all eligible features between the P1 and P2 "
            "values of the statistic named; computed on the feature-level statistics, not on "
            "proteins"
        ),
        spearman_signed_auroc_P1_P2=float(sp_signed),
        spearman_abs_auroc_minus_half_P1_P2=float(sp_abs),
        note=(
            "Exploratory discovery statistics. No pass quota is applied, the ranking is not "
            "altered if agreement is weak, and these are not held-out performance figures. "
            "No significance claim is made: the statistic was screened over thousands of "
            "features without multiplicity control."
        ),
    )

    summary = dict(
        n_proteins=N,
        n_Extracellular=int(y.sum()),
        n_Cytoplasm=int((1 - y).sum()),
        n_features=SAE_DIM,
        n_eligible=int(eligible.sum()),
        n_ineligible=int((~eligible).sum()),
        eligibility_rule=CONFIG["eligibility"],
        eligible_prevalence_min=float(prevalence[eligible].min()) if eligible.any() else None,
        n_features_never_active_in_window=int((prevalence == 0).sum()),
        n_features_never_active_anywhere=int((prevalence_visible == 0).sum()),
        training_time_dead_list_used=False,
    )
    return records, selected, stability, summary, dict(
        auroc_all=auroc_all, auroc_p1=auroc_p1, auroc_p2=auroc_p2,
        prevalence=prevalence, eligible=eligible,
    )


# ---------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--stage", choices=["extract", "rank", "all"], default="all")
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    t_all = time.time()
    dev = device_of()
    rows, sample_checks = load_sample()
    print(f"sample: {len(rows)} proteins  device={dev}  config={config_hash()[:16]}", flush=True)

    ext = None
    if args.stage in ("extract", "all"):
        ext = extract(rows, dev, limit=args.limit)
        with open(os.path.join(OUT, "extraction_report.json"), "w") as fh:
            json.dump(ext, fh, indent=2)

    if args.stage == "extract":
        return 0

    P = load_pooled(rows)
    records, selected, stability, summary, _ = rank_features(P, rows)

    write_tsv(os.path.join(OUT, "eligible_feature_ranking.tsv"), list(records[0].keys()), records)
    cand = [r for r in records if int(r["feature_id"]) in selected]
    write_tsv(os.path.join(OUT, "candidate_features.tsv"), list(cand[0].keys()), cand)
    with open(os.path.join(OUT, "selected_feature_ids.json"), "w") as fh:
        json.dump(dict(selected_feature_ids=selected, n=len(selected),
                       config_hash=config_hash()), fh, indent=2)

    sel_hash = sha256_file(os.path.join(OUT, "selected_feature_ids.json"))
    rank_hash = sha256_file(os.path.join(OUT, "eligible_feature_ranking.tsv"))

    prov = dict(
        config=CONFIG,
        config_hash=config_hash(),
        sample_checks=sample_checks,
        extraction=ext,
        ranking_summary=summary,
        stability=stability,
        selected_feature_ids=selected,
        output_hashes={
            "selected_feature_ids.json": sel_hash,
            "eligible_feature_ranking.tsv": rank_hash,
        },
        inputs={
            "analysis_view_discovery.tsv": sha256_file(VIEW),
            "deeploc_binary_manifest.tsv": sha256_file(MANIFEST),
            "esm2_plm1280_l8_sae4096.safetensors": SAE_SHA,
        },
        environment=dict(
            python=platform.python_version(), torch=torch.__version__,
            transformers=__import__("transformers").__version__,
            numpy=np.__version__, scipy=__import__("scipy").__version__,
            device=str(dev), platform=platform.platform(), machine=platform.machine(),
        ),
        total_seconds=time.time() - t_all,
    )
    with open(os.path.join(OUT, "discovery_provenance.json"), "w") as fh:
        json.dump(prov, fh, indent=2)

    print(json.dumps(dict(summary=summary, stability=stability, selected=selected,
                          selected_list_sha256=sel_hash), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
