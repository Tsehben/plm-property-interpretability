#!/usr/bin/env python3
"""Pooled ESM-2-650M representations for ALL encoder layers 1-33, one forward pass.

Identical to esm2_common/esm_extract.py in checkpoint, layer indexing, tokenisation,
truncation, pooling and compute dtype. The ONLY differences:
  (a) every layer 1..33 is pooled, not just layer 8;
  (b) pooling happens ON DEVICE so only pooled vectors (not residue-level activations)
      cross the device boundary or reach disk.
Pooling arithmetic is unchanged: hidden[r, 1:1+L, :].float(), mean over residue axis;
<cls>, <eos> and padding excluded. nterm pools the first min(30, L) residues.

Output: one .npz per shard holding keys + whole[33] + nterm[33], float32.

Usage: extract_all_layers.py --pairs P --out D [--limit N] [--batch-tokens T] [--tag TAG]
"""
import os, sys, json, time, argparse
import numpy as np
import torch
from transformers import AutoTokenizer, EsmModel

CHECKPOINT, MAX_RESIDUES, HIDDEN, NTERM = 'facebook/esm2_t33_650M_UR50D', 1022, 1280, 30
LAYERS = list(range(1, 34))     # encoder blocks 1..33; hidden_states[0] excluded
# Pooled means are COMPUTED in float32 on device and STORED as float16: a measured
# 3.58 GB free-disk constraint (float32 would need 3.43 GB of 3.58 GB available).
# Applied identically to every layer; see layer_sweep/v1/resource_check.json.
STORE_DT = 'float16'

def pick_device(req):
    if req != 'auto': return torch.device(req)
    if torch.cuda.is_available(): return torch.device('cuda')
    if torch.backends.mps.is_available(): return torch.device('mps')
    return torch.device('cpu')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pairs', required=True); ap.add_argument('--out', required=True)
    ap.add_argument('--tag', default='sweep'); ap.add_argument('--limit', type=int)
    ap.add_argument('--batch-tokens', type=int, default=12288)
    ap.add_argument('--shard-size', type=int, default=64)
    ap.add_argument('--device', default='auto'); ap.add_argument('--dtype', default='float16')
    a = ap.parse_args()
    dev = pick_device(a.device)
    os.makedirs(a.out, exist_ok=True)

    pairs = []
    with open(a.pairs) as f:
        next(f)
        for line in f:
            k, s = line.rstrip('\n').split('\t'); pairs.append((k, s))
    if a.limit: pairs = pairs[:a.limit]

    tok = AutoTokenizer.from_pretrained(CHECKPOINT)
    tdt = {'float32': torch.float32, 'float16': torch.float16,
           'bfloat16': torch.bfloat16}[a.dtype]
    model = EsmModel.from_pretrained(CHECKPOINT, add_pooling_layer=False, torch_dtype=tdt)
    model.eval().to(dev)
    for p in model.parameters(): p.requires_grad_(False)
    print(f'device={dev} checkpoint={CHECKPOINT} layers={LAYERS[0]}..{LAYERS[-1]} '
          f'pairs={len(pairs):,} batch_tokens={a.batch_tokens}', flush=True)

    # identical length-sorted batching to esm_extract.py
    order = sorted(range(len(pairs)), key=lambda i: -min(len(pairs[i][1]), MAX_RESIDUES))
    batches, cur, cur_max = [], [], 0
    for i in order:
        L = min(len(pairs[i][1]), MAX_RESIDUES) + 2
        m = max(cur_max, L)
        if cur and m * (len(cur) + 1) > a.batch_tokens: batches.append(cur); cur, cur_max = [i], L
        else: cur.append(i); cur_max = m
    if cur: batches.append(cur)
    shards = [batches[i:i + a.shard_size] for i in range(0, len(batches), a.shard_size)]

    truncated, t0, done, res_done = [], time.time(), 0, 0
    unk = tok.unk_token_id
    for si, shard in enumerate(shards):
        fp = os.path.join(a.out, f'{a.tag}_shard{si:04d}.npz')
        if os.path.exists(fp):
            print(f'  shard {si:04d}: present, resume-skip', flush=True); continue
        keys, W, N = [], [], []
        for b in shard:
            seqs = [pairs[i][1][:MAX_RESIDUES] for i in b]
            for i in b:
                if len(pairs[i][1]) > MAX_RESIDUES: truncated.append(pairs[i][0])
            enc = tok(seqs, return_tensors='pt', padding=True, add_special_tokens=True)
            ids = enc['input_ids']
            assert (ids == unk).sum().item() == 0, 'unexpected <unk> token'
            for r, s in enumerate(seqs):
                assert int(enc['attention_mask'][r].sum()) == len(s) + 2, 'token/residue mismatch'
            with torch.inference_mode():
                hs = model(input_ids=ids.to(dev),
                           attention_mask=enc['attention_mask'].to(dev),
                           output_hidden_states=True).hidden_states
            nb = len(seqs)
            w = np.empty((nb, len(LAYERS), HIDDEN), np.float32)
            n = np.empty((nb, len(LAYERS), HIDDEN), np.float32)
            for li, layer in enumerate(LAYERS):
                h = hs[layer].float()                       # (B,T,1280) on device
                for r, s in enumerate(seqs):
                    L = len(s)
                    res = h[r, 1:1 + L, :]                  # drop <cls>, <eos>, padding
                    w[r, li] = res.mean(dim=0).cpu().numpy()
                    n[r, li] = res[:min(NTERM, L)].mean(dim=0).cpu().numpy()
                del h
            del hs
            for r in range(nb):
                keys.append(pairs[b[r]][0]); W.append(w[r]); N.append(n[r])
                done += 1; res_done += len(seqs[r])
        np.savez(fp, keys=np.array(keys), whole=np.asarray(W, np.float32).astype(STORE_DT),
                 nterm=np.asarray(N, np.float32).astype(STORE_DT), layers=np.array(LAYERS))
        el = time.time() - t0
        print(f'  shard {si:04d}/{len(shards)-1}: {done:,} seqs {el:,.0f}s '
              f'{done/el:.2f} seq/s  {os.path.getsize(fp)/1e6:.1f}MB', flush=True)

    el = time.time() - t0
    meta = {'checkpoint': CHECKPOINT, 'implementation': 'transformers.EsmModel',
            'transformers_version': __import__('transformers').__version__,
            'torch_version': torch.__version__, 'device': str(dev),
            'layers_extracted': LAYERS,
            'layer_indexing_note': 'HF hidden_states[i] = output of transformer block i '
                                   '= fair-esm representations[i]; hidden_states[0] is the '
                                   'pre-block embedding output and is NOT in the sweep',
            'frozen': True, 'eval_mode': True, 'input': 'sequence only',
            'compute_dtype': a.dtype, 'pooling_dtype': 'float32',
            'stored_dtype': STORE_DT,
            'stored_dtype_reason': 'measured 3.58 GB free-disk constraint; float32 needed 3.43 GB. Applied identically to all 33 layers. Measured effect: min cosine 0.9999996 per layer, see resource_check.json',
            'pooling_device': 'on-device; only pooled vectors leave the accelerator',
            'max_residues': MAX_RESIDUES, 'hidden_dim': HIDDEN,
            'pooling': 'mean over residue tokens only; <cls>, <eos> and padding excluded',
            'nterm_positions': NTERM, 'n_sequences': done,
            'n_truncated': len(set(truncated)), 'truncated_keys': sorted(set(truncated)),
            'elapsed_s': round(el, 1), 'seq_per_s': round(done / el, 3) if el else None}
    json.dump(meta, open(os.path.join(a.out, f'{a.tag}_extraction_metadata.json'), 'w'), indent=1)
    print(f'\ndone: {done:,} seqs in {el:,.0f}s ({done/el:.2f} seq/s)', flush=True)

if __name__ == '__main__':
    main()
