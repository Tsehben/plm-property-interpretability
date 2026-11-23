#!/usr/bin/env python3
"""
Stage 4: describe the already-selected, already-hashed candidate features.

Runs only after run_discovery.py has written and hashed selected_feature_ids.json.
It re-reads that file, records its SHA-256, and never alters it. Nothing here can
change the candidate list; that is enforced by comparing the hash at the start and
at the end.

Descriptions produced, all on DISCOVERY data only:
  1. activation distributions by localisation class;
  2. annotated_present versus not_annotated signal-peptide status, with
     not_annotated carried as its own state and never as confirmed absence;
  3. the experimentally evidenced subset separately;
  4. residue-level activation inside transferable signal-peptide intervals against
     the remaining visible sequence, aggregated per protein (streamed recompute;
     no dense per-residue array is stored);
  5. kingdom and original-length distributions associated with activation;
  6. counterexamples in every direction.

Usage:  python3 characterise_candidates.py
"""

import csv
import hashlib
import json
import os
import sys
import time

import numpy as np
import torch
from transformers import AutoTokenizer, EsmModel

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "v1")
POOLED = os.path.join(OUT, "pooled")
REF = os.path.join(ROOT, "interprot_ref", "raw")
sys.path.insert(0, REF)
from sae_model import SparseAutoencoder  # noqa: E402

import run_discovery as RD  # noqa: E402  (frozen config, sample loader, helpers)

FEATURES_TSV = os.path.join(ROOT, "annotation_readiness", "out", "annotation_features.tsv")
INCLUSION = os.path.join(ROOT, "discovery_analysis", "v1", "inclusion_manifest.tsv")

NOT_ANNOTATED_NOTE = (
    "not_annotated means a matching UniProt record was retrieved and carries no "
    "signal-peptide annotation. It is NOT a confirmed biological negative and is "
    "never pooled with unverifiable or treated as absence."
)


def q(a):
    a = np.asarray(a, dtype=np.float64)
    if a.size == 0:
        return dict(n=0)
    return dict(
        n=int(a.size),
        mean=float(a.mean()),
        sd=float(a.std(ddof=1)) if a.size > 1 else 0.0,
        p25=float(np.percentile(a, 25)),
        median=float(np.median(a)),
        p75=float(np.percentile(a, 75)),
        p95=float(np.percentile(a, 95)),
        max=float(a.max()),
        frac_zero=float((a == 0).mean()),
    )


def main():
    t0 = time.time()
    sel_path = os.path.join(OUT, "selected_feature_ids.json")
    sel_hash_before = RD.sha256_file(sel_path)
    sel = json.load(open(sel_path))["selected_feature_ids"]
    print(f"candidates: {sel}\nselected list sha256: {sel_hash_before}", flush=True)

    rows, _ = RD.load_sample()
    P = RD.load_pooled(rows)
    acc_index = {r["ACC"]: i for i, r in enumerate(rows)}
    y = np.array([r["target"] for r in rows])
    kingdom = np.array([r["kingdom"] for r in rows])
    olen = np.array([r["original_len"] for r in rows], dtype=np.float64)
    trunc = np.array([r["truncated"] for r in rows])

    incl = {r["ACC"]: r for r in RD.read_tsv(INCLUSION)}
    sig_status = np.array([incl[r["ACC"]]["signal_status"] for r in rows])
    sig_evid = np.array([incl[r["ACC"]]["signal_evidence_best"] for r in rows])
    in_ann = np.array([incl[r["ACC"]]["in_annotated"] == "1" for r in rows])

    X = P["mean_nterm"]
    A = P["nact_nterm"] > 0

    # ---------------------------------------------------------------- 1 + 2 + 3
    by_class, by_signal = [], []
    for f in sel:
        x = X[:, f]
        for name, m in (("Extracellular", y == 1), ("Cytoplasm", y == 0)):
            by_class.append(dict(feature_id=f, group=name, prevalence_nterm=float(A[m, f].mean()), **q(x[m])))
        # annotation-conditioned discovery sample only
        base = in_ann
        for name, m in (
            ("annotated_present", base & (sig_status == "annotated_present")),
            ("not_annotated", base & (sig_status == "not_annotated")),
            ("annotated_present_experimental_evidence_only",
             base & (sig_status == "annotated_present") & (sig_evid == "experimental")),
            ("annotated_present_computational_evidence_only",
             base & (sig_status == "annotated_present") & (sig_evid == "sequence_analysis_predictor")),
            ("Extracellular_and_not_annotated", base & (y == 1) & (sig_status == "not_annotated")),
            ("Cytoplasm_and_annotated_present", base & (y == 0) & (sig_status == "annotated_present")),
        ):
            by_signal.append(dict(feature_id=f, group=name, prevalence_nterm=float(A[m, f].mean()) if m.sum() else 0.0, **q(x[m])))

    # ------------------------------------------------------------------- 5
    by_kingdom, by_length = [], []
    from scipy.stats import spearmanr
    for f in sel:
        x = X[:, f]
        for k in sorted(set(kingdom)):
            m = kingdom == k
            by_kingdom.append(dict(feature_id=f, kingdom=k, prevalence_nterm=float(A[m, f].mean()), **q(x[m])))
        by_length.append(
            dict(
                feature_id=f,
                spearman_activation_vs_original_length=float(spearmanr(x, olen).statistic),
                mean_original_length_active=float(olen[A[:, f]].mean()) if A[:, f].any() else float("nan"),
                mean_original_length_inactive=float(olen[~A[:, f]].mean()) if (~A[:, f]).any() else float("nan"),
                n_truncated_proteins=int(trunc.sum()),
                mean_activation_truncated=float(x[trunc == 1].mean()) if (trunc == 1).any() else float("nan"),
                mean_activation_not_truncated=float(x[trunc == 0].mean()),
            )
        )

    # ------------------------------------------------------------------- 6
    counterex = []
    for f in sel:
        x = X[:, f]
        def top(mask, k=5, reverse=False):
            idx = np.where(mask)[0]
            if idx.size == 0:
                return []
            o = idx[np.argsort(x[idx])]
            o = o[:k] if reverse else o[::-1][:k]
            return [(rows[i]["ACC"], float(x[i])) for i in o]
        counterex.append(
            dict(
                feature_id=f,
                highest_activating_Cytoplasm=json.dumps(top(y == 0)),
                lowest_activating_Extracellular=json.dumps(top(y == 1, reverse=True)),
                highest_activating_not_annotated=json.dumps(top(in_ann & (sig_status == "not_annotated"))),
                zero_activation_annotated_present=int(
                    ((x == 0) & in_ann & (sig_status == "annotated_present")).sum()
                ),
                n_Extracellular_with_zero_activation=int(((x == 0) & (y == 1)).sum()),
                n_Cytoplasm_above_Extracellular_median=int(
                    ((y == 0) & (x > np.median(x[y == 1]))).sum()
                ),
            )
        )

    # ------------------------------------------------------------------- 4
    enr = signal_interval_enrichment(rows, acc_index, sel, incl)

    RD.write_tsv(os.path.join(OUT, "activation_by_class.tsv"), list(by_class[0].keys()), by_class)
    RD.write_tsv(os.path.join(OUT, "activation_by_signal_status.tsv"), list(by_signal[0].keys()), by_signal)
    RD.write_tsv(os.path.join(OUT, "activation_by_kingdom.tsv"), list(by_kingdom[0].keys()), by_kingdom)
    RD.write_tsv(os.path.join(OUT, "activation_vs_length.tsv"), list(by_length[0].keys()), by_length)
    RD.write_tsv(os.path.join(OUT, "counterexamples.tsv"), list(counterex[0].keys()), counterex)
    RD.write_tsv(os.path.join(OUT, "signal_interval_enrichment.tsv"), list(enr["rows"][0].keys()), enr["rows"])

    sel_hash_after = RD.sha256_file(sel_path)
    assert sel_hash_after == sel_hash_before, "the candidate list changed during characterisation"

    report = dict(
        selected_feature_ids=sel,
        selected_list_sha256=sel_hash_before,
        candidate_list_unchanged=True,
        not_annotated_note=NOT_ANNOTATED_NOTE,
        interpretation_limits=[
            "Annotation agreement is descriptive co-occurrence. It is not causal evidence "
            "and does not show that a feature encodes localisation or a signal peptide.",
            "Untyped InterPro entries are not treated as families anywhere in this report.",
            "No significance claim is made; the ranking screened thousands of features "
            "without multiplicity control.",
        ],
        signal_interval_scope=enr["scope"],
        output_hashes={
            f: RD.sha256_file(os.path.join(OUT, f))
            for f in (
                "activation_by_class.tsv", "activation_by_signal_status.tsv",
                "activation_by_kingdom.tsv", "activation_vs_length.tsv",
                "counterexamples.tsv", "signal_interval_enrichment.tsv",
            )
        },
        seconds=time.time() - t0,
    )
    with open(os.path.join(OUT, "characterisation_report.json"), "w") as fh:
        json.dump(report, fh, indent=2)
    print(json.dumps(dict(scope=enr["scope"], hashes=report["output_hashes"]), indent=2))
    return 0


def signal_interval_enrichment(rows, acc_index, sel, incl):
    """Streamed recompute. For each discovery protein that has a signal-peptide
    annotation whose coordinates are verified to apply to the frozen sequence,
    compute, per candidate feature, the mean activation inside the annotated
    interval and over the remaining visible residues. Only 2 x 20 scalars per
    protein are retained."""
    intervals = {}
    for r in RD.read_tsv(FEATURES_TSV):
        if r["feature_type"] != "Signal" or r["coordinates_apply_to_frozen_sequence"] != "1":
            continue
        if r["start_modifier"] != "EXACT" or r["end_modifier"] != "EXACT":
            continue
        intervals.setdefault(r["ACC"], (int(r["start"]), int(r["end"])))

    use = []
    skipped_truncation = 0
    for r in rows:
        iv = intervals.get(r["ACC"])
        if iv is None:
            continue
        if incl[r["ACC"]]["in_annotated"] != "1":
            continue
        s, e = iv
        if e > r["visible_len"] or s < 1 or e < s:
            skipped_truncation += 1
            continue
        if r["visible_len"] - (e - s + 1) < 10:      # need a usable comparison region
            skipped_truncation += 1
            continue
        use.append((r, s, e))

    scope = dict(
        n_proteins_used=len(use),
        n_skipped_interval_outside_visible_or_too_short=skipped_truncation,
        n_Extracellular=sum(1 for r, _, _ in use if r["target"] == 1),
        n_Cytoplasm=sum(1 for r, _, _ in use if r["target"] == 0),
        rule=("signal-peptide features with EXACT start and end whose coordinates are "
              "verified to apply to the frozen extraction-input sequence, in the "
              "annotation-conditioned discovery sample"),
        comparison="mean activation over residues [start, end] versus over all other visible residues",
    )
    print(f"interval enrichment over {len(use)} proteins", flush=True)

    dev = RD.device_of()
    tok = AutoTokenizer.from_pretrained(RD.CHECKPOINT, revision=RD.CHECKPOINT_REVISION)
    model = EsmModel.from_pretrained(
        RD.CHECKPOINT, revision=RD.CHECKPOINT_REVISION, add_pooling_layer=False,
        torch_dtype=torch.float32,
    )
    model.eval().to(dev)
    for p in model.parameters():
        p.requires_grad_(False)
    from safetensors.torch import load_file
    sae = SparseAutoencoder(RD.HIDDEN, RD.SAE_DIM)
    sae.load_state_dict(load_file(RD.SAE_FILE))
    sae.eval()

    selarr = np.array(sel)
    inside = np.zeros((len(use), len(sel)), np.float64)
    outside = np.zeros((len(use), len(sel)), np.float64)
    use.sort(key=lambda t: -t[0]["visible_len"])
    t0 = time.time()
    for b in range(0, len(use), RD.BATCH):
        chunk = use[b : b + RD.BATCH]
        seqs = [r["seq"][: RD.MAX_RESIDUES] for r, _, _ in chunk]
        acts, _ = RD.esm_residue_acts(model, tok, seqs, dev)
        for j, (r, s, e) in enumerate(chunk):
            lat = sae.get_acts(acts[j]).numpy()[:, selarr]       # (L, 20)
            m = np.zeros(lat.shape[0], bool)
            m[s - 1 : e] = True
            inside[b + j] = lat[m].mean(axis=0)
            outside[b + j] = lat[~m].mean(axis=0)
            del lat
        if b % (RD.BATCH * 50) == 0:
            print(f"  {b}/{len(use)}  {time.time() - t0:.0f}s", flush=True)

    out = []
    for k, f in enumerate(sel):
        d = inside[:, k] - outside[:, k]
        out.append(
            dict(
                feature_id=int(f),
                n_proteins=len(use),
                mean_inside_interval=float(inside[:, k].mean()),
                mean_outside_interval=float(outside[:, k].mean()),
                mean_per_protein_difference=float(d.mean()),
                median_per_protein_difference=float(np.median(d)),
                n_proteins_inside_greater=int((d > 0).sum()),
                n_proteins_outside_greater=int((d < 0).sum()),
                n_proteins_equal=int((d == 0).sum()),
                n_proteins_zero_inside=int((inside[:, k] == 0).sum()),
                n_proteins_zero_outside=int((outside[:, k] == 0).sum()),
                note="descriptive co-occurrence; no causal or encoding claim",
            )
        )
    return dict(rows=out, scope=scope)


if __name__ == "__main__":
    sys.exit(main())
