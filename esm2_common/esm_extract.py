#!/usr/bin/env python3
"""
Frozen ESM-2-650M representation extraction, shared by both datasets.

CHECKPOINT / IMPLEMENTATION / LAYER INDEXING  (recorded into extraction_metadata.json)
  checkpoint      facebook/esm2_t33_650M_UR50D   (ESM-2, 33 layers, 1280-dim)
  implementation  HuggingFace transformers, EsmModel
  mode            model.eval(), torch.inference_mode() - frozen, no gradients, no fine-tuning
  LAYER           hidden_states[8]
                  HF hidden_states[0] is the embedding output BEFORE transformer block 0;
                  hidden_states[i] is the OUTPUT of transformer block i. So hidden_states[8]
                  is the output of the 8th transformer block, which is the same tensor
                  fair-esm returns as representations[8] (repr_layers=[8]).
                  This is the layer studied by the supplied InterProt paper and the layer
                  its released SAE was trained on.
  input           sequence alone. No kingdom, no header text, no annotation of any kind.

TOKENS / POOLING
  The ESM-2 vocabulary has dedicated ids for the ambiguous residues B, U, Z, O and X, so
  every residue maps to exactly one token and nothing becomes <unk> (verified per batch).
  Tokenization is <cls> + L residue tokens + <eos>. Pooling uses ONLY the L residue
  positions: <cls>, <eos> and all padding are excluded. Each sequence's residue block is
  sliced as hidden[1 : 1+L] and the per-sequence token count is asserted against L.

TRUNCATION
  Sequences longer than MAX_RESIDUES (1022) are truncated to their first 1022 residues for
  this first pass. The affected keys are recorded.

RESUME
  Pooled vectors are written incrementally as shards (one .npz per shard). A completed
  shard is never recomputed, so an interrupted job resumes by re-running the same command.
"""
import os, sys, json, time, argparse
import numpy as np
import torch
from transformers import AutoTokenizer, EsmModel

CHECKPOINT   = 'facebook/esm2_t33_650M_UR50D'
LAYER        = 8
MAX_RESIDUES = 1022
HIDDEN       = 1280

def pick_device(requested):
    if requested != 'auto': return torch.device(requested)
    if torch.cuda.is_available(): return torch.device('cuda')
    if torch.backends.mps.is_available(): return torch.device('mps')
    return torch.device('cpu')

def load_pairs(path):
    """TSV: key <TAB> sequence (with header)."""
    out = []
    with open(path) as f:
        next(f)
        for line in f:
            k, s = line.rstrip('\n').split('\t')
            out.append((k, s))
    return out

def extract(pairs, out_dir, device, nterm=0, batch_tokens=16384, shard_size=512,
            layer=LAYER, limit=None, tag='run', dtype='float32'):
    os.makedirs(out_dir, exist_ok=True)
    if limit: pairs = pairs[:limit]

    tok = AutoTokenizer.from_pretrained(CHECKPOINT)
    torch_dtype = {'float32': torch.float32, 'float16': torch.float16,
                   'bfloat16': torch.bfloat16}[dtype]
    model = EsmModel.from_pretrained(CHECKPOINT, add_pooling_layer=False,
                                     torch_dtype=torch_dtype)
    model.eval().to(device)
    for p in model.parameters(): p.requires_grad_(False)       # frozen

    # length-sorted batching for efficiency; original keys are carried through
    order = sorted(range(len(pairs)), key=lambda i: -min(len(pairs[i][1]), MAX_RESIDUES))
    batches, cur, cur_max = [], [], 0
    for i in order:
        L = min(len(pairs[i][1]), MAX_RESIDUES) + 2
        m = max(cur_max, L)
        if cur and m * (len(cur) + 1) > batch_tokens:
            batches.append(cur); cur, cur_max = [i], L
        else:
            cur.append(i); cur_max = m
    if cur: batches.append(cur)

    shards = [batches[i:i + shard_size] for i in range(0, len(batches), shard_size)]
    truncated, t0, done_seqs, done_res = [], time.time(), 0, 0
    unk_id = tok.unk_token_id

    for si, shard in enumerate(shards):
        fp = os.path.join(out_dir, f'{tag}_shard{si:04d}.npz')
        if os.path.exists(fp):
            print(f'  shard {si:04d}: already present, skipping (resume)', flush=True)
            continue
        keys, whole, nt = [], [], []
        for b in shard:
            seqs = [pairs[i][1][:MAX_RESIDUES] for i in b]
            for i in b:
                if len(pairs[i][1]) > MAX_RESIDUES: truncated.append(pairs[i][0])
            enc = tok(seqs, return_tensors='pt', padding=True, add_special_tokens=True)
            ids = enc['input_ids']
            assert (ids == unk_id).sum().item() == 0, 'unexpected <unk> token'
            # exactly one residue token per residue, plus <cls> and <eos>
            for r, s in enumerate(seqs):
                assert int(enc['attention_mask'][r].sum()) == len(s) + 2, \
                    f'token/residue mismatch for {pairs[b[r]][0]}'
            with torch.inference_mode():
                h = model(input_ids=ids.to(device),
                          attention_mask=enc['attention_mask'].to(device),
                          output_hidden_states=True).hidden_states[layer]
            h = h.float().cpu().numpy()
            for r, s in enumerate(seqs):
                L = len(s)
                res = h[r, 1:1 + L, :]                 # drop <cls>, <eos> and all padding
                assert res.shape == (L, HIDDEN), f'bad slice for {pairs[b[r]][0]}'
                keys.append(pairs[b[r]][0]); whole.append(res.mean(axis=0))
                if nterm: nt.append(res[:min(nterm, L), :].mean(axis=0))
                done_seqs += 1; done_res += L
        d = {'keys': np.array(keys), 'whole': np.asarray(whole, np.float32)}
        if nterm: d['nterm'] = np.asarray(nt, np.float32)
        np.savez_compressed(fp, **d)
        el = time.time() - t0
        print(f'  shard {si:04d}/{len(shards)-1}: {done_seqs:,} seqs  {el:,.0f}s  '
              f'{done_seqs/el:.1f} seq/s  {done_res/el:,.0f} res/s', flush=True)
    return truncated, time.time() - t0, done_seqs

def gather(out_dir, tag, nterm=False):
    ks, ws, ns = [], [], []
    for f in sorted(os.listdir(out_dir)):
        if not (f.startswith(tag + '_shard') and f.endswith('.npz')): continue
        z = np.load(os.path.join(out_dir, f), allow_pickle=True)
        ks.append(z['keys']); ws.append(z['whole'])
        if nterm: ns.append(z['nterm'])
    K = np.concatenate(ks); W = np.concatenate(ws)
    N = np.concatenate(ns) if nterm else None
    return K, W, N

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pairs', required=True); ap.add_argument('--out', required=True)
    ap.add_argument('--tag', default='run'); ap.add_argument('--nterm', type=int, default=0)
    ap.add_argument('--device', default='auto'); ap.add_argument('--limit', type=int)
    ap.add_argument('--batch-tokens', type=int, default=16384)
    ap.add_argument('--shard-size', type=int, default=512)
    ap.add_argument('--dtype', default='float32',
                    choices=['float32', 'float16', 'bfloat16'])
    a = ap.parse_args()
    dev = pick_device(a.device)
    pairs = load_pairs(a.pairs)
    print(f'device={dev}  checkpoint={CHECKPOINT}  layer={LAYER}  '
          f'pairs={len(pairs):,}' + (f'  LIMIT={a.limit}' if a.limit else ''))
    trunc, el, n = extract(pairs, a.out, dev, nterm=a.nterm, batch_tokens=a.batch_tokens,
                           shard_size=a.shard_size, limit=a.limit, tag=a.tag,
                           dtype=a.dtype)
    meta = {'checkpoint': CHECKPOINT, 'implementation': 'transformers.EsmModel',
            'transformers_version': __import__('transformers').__version__,
            'torch_version': torch.__version__, 'device': str(dev),
            'layer_index': LAYER,
            'layer_indexing_note': 'HF hidden_states[8] = output of transformer block 8 '
                                   '= fair-esm representations[8]',
            'frozen': True, 'eval_mode': True, 'input': 'sequence only',
            'compute_dtype': a.dtype, 'stored_dtype': 'float32',
            'max_residues': MAX_RESIDUES, 'hidden_dim': HIDDEN,
            'pooling': 'mean over residue tokens only; <cls>, <eos> and padding excluded',
            'nterm_positions': a.nterm,
            'n_sequences': n, 'n_truncated': len(trunc),
            'truncated_keys': sorted(set(trunc)),
            'elapsed_s': round(el, 1), 'seq_per_s': round(n / el, 2) if el else None}
    with open(os.path.join(a.out, f'{a.tag}_extraction_metadata.json'), 'w') as f:
        json.dump(meta, f, indent=1)
    print(f'\ndone: {n:,} sequences in {el:,.0f}s ({n/el:.1f} seq/s), '
          f'truncated={len(trunc)}')

if __name__ == '__main__':
    main()
