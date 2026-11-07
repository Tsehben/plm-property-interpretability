#!/usr/bin/env python3
"""
Layer-8 InterProt SAE compatibility pilot, shared by both datasets.

SAE  liambai/InterProt-ESM2-SAEs @ revision 3ded2e51641df84d6b89f25430874e3a8eb42c24
     file esm2_plm1280_l8_sae4096.safetensors
     sha256 1172503b90c950651e787a5625d91ff1c702a04ff3d851e8e95b700201e56ebd
     class SparseAutoencoder(d_model=1280, d_hidden=4096, k=128)   [constructor defaults,
     exactly as the authors' README instantiates it:  SparseAutoencoder(ESM_DIM, SAE_DIM)]

TENSOR CORRESPONDENCE.  The authors' reference inference code loads
  EsmModel.from_pretrained('facebook/esm2_t33_650M_UR50D'), runs it with
  output_hidden_states=True and reads outputs.hidden_states[LAYER].
That is the identical checkpoint, identical class and identical indexing convention we
already use, so the tensor we extract is the tensor this SAE was trained on -- established
by matching the implementation path, not merely the layer number and the dimension.

PREPROCESSING (authors', reproduced exactly via their own class):
  get_acts(x) = topK( LN(x) - b_pre ) @ w_enc + b_enc ,  k = 128, ReLU on kept values
  LN normalises over dim=-1 (per token, across the 1280 features) with eps=1e-5 and
  torch's default *unbiased* std.  Because LN is per-token, applying the SAE to a subset
  of tokens is token-identical to applying it to the full (L+2) sequence and slicing
  afterwards -- this script verifies that numerically rather than assuming it.
"""
import csv, hashlib, json, os, sys, time, warnings
import numpy as np
import torch
from transformers import AutoTokenizer, EsmModel

warnings.filterwarnings('ignore', message='.*encountered in matmul.*',
                        category=RuntimeWarning)

# Repo root: the parent of esm2_common/.  Override with MIT6393_ROOT if the
# checkout lives elsewhere.
ROOT = os.environ.get('MIT6393_ROOT',
                      os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REF  = f'{ROOT}/interprot_ref/raw'
sys.path.insert(0, REF)                       # reference class, read in full before use
from sae_model import SparseAutoencoder       # noqa: E402

CHECKPOINT = 'facebook/esm2_t33_650M_UR50D'
SAE_FILE   = f'{REF}/esm2_plm1280_l8_sae4096.safetensors'
SAE_SHA    = '1172503b90c950651e787a5625d91ff1c702a04ff3d851e8e95b700201e56ebd'
LAYER, MAX_RESIDUES, HIDDEN, SAE_DIM, TOPK = 8, 1022, 1280, 4096, 128
NTERM = 30

DATASETS = {
    'uesolds': dict(dir=f'{ROOT}/Violet969-PLM_Sol-3589a95', emb_tag='uesolds', nterm=False),
    'deeploc': dict(dir=f'{ROOT}/DeepLoc2',                  emb_tag='deeploc', nterm=True),
}


def full_counts():
    import csv as _csv
    out = {}
    with open(f'{ROOT}/Violet969-PLM_Sol-3589a95/dev_cohort/dev_cohort_manifest.tsv',
              newline='') as f:
        rs = list(_csv.DictReader(f, delimiter='\t'))
    out['uesolds'] = dict(proteins=len(rs),
                          residues=sum(min(int(r['seq_len']), MAX_RESIDUES) for r in rs))
    with open(f'{ROOT}/DeepLoc2/baseline/deeploc_binary_manifest.tsv', newline='') as f:
        rs = [r for r in _csv.DictReader(f, delimiter='\t')
              if r['split_role'] == 'development_P1_4']
    out['deeploc'] = dict(proteins=len(rs),
                          residues=sum(min(int(r['len']), MAX_RESIDUES) for r in rs))
    return out


FULL_COUNTS = None


def device_of():
    if torch.cuda.is_available(): return torch.device('cuda')
    if torch.backends.mps.is_available(): return torch.device('mps')
    return torch.device('cpu')


def load_sample(d):
    with open(f'{d}/sae_pilot/sae_pilot_sample.tsv', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def load_saved_emb(d, tag):
    """Existing pooled whole-input means from the completed fp16 extraction."""
    import glob
    ks, ws = [], []
    for f in sorted(glob.glob(f'{d}/esm2_probe/emb/{tag}_shard*.npz')):
        z = np.load(f, allow_pickle=True)
        ks.append(z['keys']); ws.append(z['whole'])
    K = np.concatenate(ks); W = np.concatenate(ws)
    return {str(k): W[i] for i, k in enumerate(K)}


def esm_residue_acts(model, tok, seqs, device, dtype):
    """Return list of (L, 1280) float32 numpy arrays: layer-8 residue activations."""
    enc = tok(seqs, padding=True, return_tensors='pt')
    ids, am = enc['input_ids'].to(device), enc['attention_mask'].to(device)
    assert (ids == tok.unk_token_id).sum().item() == 0, 'unexpected <unk> token'
    for r, s in enumerate(seqs):
        assert int(am[r].sum()) == len(s) + 2, f'token/residue mismatch len {len(s)}'
    with torch.inference_mode():
        h = model(input_ids=ids, attention_mask=am,
                  output_hidden_states=True).hidden_states[LAYER]
    out = []
    for r, s in enumerate(seqs):
        L = len(s)
        out.append(h[r, 1:1 + L, :].float().cpu().numpy())   # drop <cls>, <eos>, padding
    return out, h


def main():
    global FULL_COUNTS
    FULL_COUNTS = full_counts()
    t_all = time.time()
    dev = device_of()
    assert hashlib.sha256(open(SAE_FILE, 'rb').read()).hexdigest() == SAE_SHA, \
        'SAE checkpoint hash mismatch'

    print(f'device={dev}  torch={torch.__version__}')
    tok = AutoTokenizer.from_pretrained(CHECKPOINT)

    # ---- SAE (reference class, authors' weights, fp32 on CPU) ---------------
    from safetensors.torch import load_file
    sae = SparseAutoencoder(HIDDEN, SAE_DIM)
    sae.load_state_dict(load_file(SAE_FILE))
    sae.eval()
    assert sae.k == TOPK, sae.k
    dead = (sae.stats_last_nonzero > sae.dead_steps_threshold).sum().item()
    print(f'SAE loaded: d_model={sae.d_model} d_hidden={sae.d_hidden} k={sae.k}  '
          f'features marked dead at end of training: {dead}/{SAE_DIM}')

    report = {'sae': {'repo': 'liambai/InterProt-ESM2-SAEs',
                      'revision': '3ded2e51641df84d6b89f25430874e3a8eb42c24',
                      'file': 'esm2_plm1280_l8_sae4096.safetensors', 'sha256': SAE_SHA,
                      'd_model': HIDDEN, 'd_hidden': SAE_DIM, 'k': TOPK,
                      'dead_features_at_train_end': dead},
              'datasets': {}}

    for dtype_name in ['float32']:
        model = EsmModel.from_pretrained(CHECKPOINT, add_pooling_layer=False,
                                         torch_dtype=torch.float32)
        model.eval().to(dev)
        for p in model.parameters(): p.requires_grad_(False)
        break

    # fp16 copy, to test whether top-k survives the dtype we used at scale
    model16 = EsmModel.from_pretrained(CHECKPOINT, add_pooling_layer=False,
                                       torch_dtype=torch.float16)
    model16.eval().to(dev)
    for p in model16.parameters(): p.requires_grad_(False)

    for name, cfg in DATASETS.items():
        d, tag = cfg['dir'], cfg['emb_tag']
        rows = load_sample(d)
        saved = load_saved_emb(d, tag)
        outdir = f'{d}/sae_pilot'
        print(f'\n{"="*72}\n{name.upper()}  n={len(rows)}\n{"="*72}')

        R = {}   # per-protein results
        cos_mean, maxrel, absdiff, tiny = [], [], [], []
        pooled_w16, pooled_n16 = [], []
        act_counts, recon_rel, recon_base = [], [], []
        seen_feats = np.zeros(SAE_DIM, dtype=bool)
        nonfinite = 0
        topk_agree16 = []
        pooled_w, pooled_n, keys = [], [], []
        t0 = time.time(); n_res = 0

        B = 4
        for i in range(0, len(rows), B):
            batch = rows[i:i + B]
            seqs = [r['sequence'][:MAX_RESIDUES] for r in batch]
            acts, hfull = esm_residue_acts(model, tok, seqs, dev, torch.float32)
            acts16, _ = esm_residue_acts(model16, tok, seqs, dev, torch.float16)

            for j, r in enumerate(batch):
                key, A = r['key'], acts[j]
                L = len(seqs[j])
                assert A.shape == (L, HIDDEN), A.shape
                assert L == int(r['visible_len']), (L, r['visible_len'])
                assert hashlib.sha256(r['sequence'].encode()).hexdigest() == r['sha256']
                n_res += L

                # --- (a) raw residue mean must reproduce the saved embedding ---
                m = A.mean(axis=0)
                sv = saved[key]
                cos_mean.append(float(m @ sv / (np.linalg.norm(m) * np.linalg.norm(sv))))
                # relative to the VECTOR NORM: a per-component ratio explodes wherever a
                # component of the saved mean is itself ~0, which says nothing about agreement.
                maxrel.append(float(np.abs(m - sv).max() / np.linalg.norm(sv)))
                absdiff.append(float(np.abs(m - sv).max()))
                tiny.append(float(np.abs(sv).min()))

                # --- (b) SAE on individual residue vectors, authors' preprocessing ---
                x = torch.from_numpy(A)
                lat = sae.get_acts(x)                      # (L, 4096)
                latn = lat.numpy()
                if not np.isfinite(latn).all(): nonfinite += 1
                act_counts.append((latn > 0).sum(axis=1))
                seen_feats |= (latn > 0).any(axis=0)

                # --- (c) reconstruction error vs a mean-vector baseline ---------
                rec = sae.forward_val(x).numpy()
                num = float(((rec - A) ** 2).sum())
                recon_rel.append(num / float((A ** 2).sum()))
                mu = A.mean(axis=0, keepdims=True)
                recon_base.append(num / float(((A - mu) ** 2).sum()))

                # --- (d) does fp16 ESM change which features fire? --------------
                lat16 = sae.get_acts(torch.from_numpy(acts16[j])).numpy()
                s1 = (latn > 0); s2 = (lat16 > 0)
                topk_agree16.append(float((s1 & s2).sum() / max(s1.sum(), 1)))

                # --- (e) pooled SAE features over visible residues --------------
                keys.append(key)
                pooled_w.append(latn.mean(axis=0))
                pooled_w16.append(lat16.mean(axis=0))
                if cfg['nterm']:
                    pooled_n.append(latn[:min(NTERM, L)].mean(axis=0))
                    pooled_n16.append(lat16[:min(NTERM, L)].mean(axis=0))

            # --- (f) reference-equivalence check on the first batch ------------
            if i == 0:
                ref_ok, ref_max = [], 0.0
                for j, s in enumerate(seqs):
                    # authors' exact path: full (L+2) tensor straight from hidden_states
                    full = hfull[j, :len(s) + 2, :].float().cpu()
                    ref = sae.get_acts(full).numpy()[1:1 + len(s)]
                    mine = sae.get_acts(torch.from_numpy(acts[j])).numpy()
                    ref_max = max(ref_max, float(np.abs(ref - mine).max()))
                    ref_ok.append(bool((ref > 0).sum() == (mine > 0).sum()
                                       and ((ref > 0) == (mine > 0)).all()))
                print(f'  reference check (batch of {len(seqs)}): '
                      f'max |ref-ours| = {ref_max:.3e}, '
                      f'active-set identical on {sum(ref_ok)}/{len(ref_ok)}')

        rt = time.time() - t0
        ac = np.concatenate(act_counts)
        Pw = np.stack(pooled_w)
        Pn = np.stack(pooled_n) if cfg['nterm'] else None

        # alignment verification by id and hash
        assert keys == [r['key'] for r in rows], 'pooled order != sample order'
        assert len(set(keys)) == len(keys), 'duplicate keys'

        np.savez_compressed(f'{outdir}/sae_pilot_pooled.npz',
                            keys=np.array(keys), whole=Pw.astype(np.float32),
                            **({'nterm': Pn.astype(np.float32)} if Pn is not None else {}))

        Pw16 = np.stack(pooled_w16)
        cos_pool = [float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))
                    for a, b in zip(Pw, Pw16)]
        dens = float((Pw > 0).mean())
        sz_dense = SAE_DIM * 4 / 1e6                      # MB per protein per vector, fp32

        full = FULL_COUNTS[name]
        est_s = full['residues'] / (n_res / rt)
        est_mb = full['proteins'] * sz_dense * (2 if cfg['nterm'] else 1)

        res = dict(
            n_proteins=len(rows), n_residues=int(n_res), runtime_s=round(rt, 1),
            residues_per_s=round(n_res / rt, 1),
            mean_cosine_vs_saved=float(np.mean(cos_mean)),
            min_cosine_vs_saved=float(np.min(cos_mean)),
            max_absdiff_over_norm_vs_saved=float(np.max(maxrel)),
            max_abs_diff_vs_saved=float(np.max(absdiff)),
            smallest_saved_component_abs=float(np.min(tiny)),
            nonfinite_proteins=nonfinite,
            active_per_residue_mean=float(ac.mean()),
            active_per_residue_min=int(ac.min()), active_per_residue_max=int(ac.max()),
            recon_rel_err_vs_zero=float(np.mean(recon_rel)),
            recon_rel_err_vs_meanvec=float(np.mean(recon_base)),
            frac_latents_ever_active=float(seen_feats.mean()),
            n_latents_ever_active=int(seen_feats.sum()),
            fp16_topk_overlap_mean=float(np.mean(topk_agree16)),
            fp16_topk_overlap_min=float(np.min(topk_agree16)),
            pooled_cosine_fp32_vs_fp16_mean=float(np.mean(cos_pool)),
            pooled_cosine_fp32_vs_fp16_min=float(np.min(cos_pool)),
            pooled_density_frac_nonzero=dens,
            full_run_proteins=full['proteins'], full_run_residues=full['residues'],
            full_run_est_runtime_min=round(est_s / 60, 1),
            full_run_est_storage_MB=round(est_mb, 1),
            pooled_shape_whole=list(Pw.shape),
            pooled_shape_nterm=(list(Pn.shape) if Pn is not None else None),
        )
        report['datasets'][name] = res
        for k, v in res.items(): print(f'  {k:28s} {v}')

    with open(f'{ROOT}/esm2_common/sae_pilot_report.json', 'w') as f:
        json.dump(report, f, indent=2)
    print(f'\ntotal {time.time()-t_all:.1f}s -> esm2_common/sae_pilot_report.json')


if __name__ == '__main__':
    main()
