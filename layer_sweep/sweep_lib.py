#!/usr/bin/env python3
"""Shared extraction/store library for the layer sweep (v2).

Used by BOTH run_sweep_v2.py and test_checkpoint_resume.py so the resume test exercises
exactly the code path the real run uses.
"""
import os, sys, csv, json, time, hashlib, subprocess, collections
import numpy as np

CHECKPOINT_ID = 'facebook/esm2_t33_650M_UR50D'
MAX_RES, HID, NTERM = 1022, 1280, 30
LAYERS = list(range(1, 34))
STORE_DTYPE = 'float16'
COMPUTE_DTYPE = 'float16'
POOLING_DTYPE = 'float32'
POOLING_DESC = 'mean over residue tokens only; <cls>, <eos> and padding excluded'
DISK_MIN_BYTES = 500 * 1024 * 1024          # hard stop threshold


def file_sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()


def extraction_code_sha():
    """Hash of the extraction code itself, so a changed implementation refuses to resume."""
    return file_sha(os.path.abspath(__file__))


def resources(path='.'):
    st = os.statvfs(path)
    free = st.f_bavail * st.f_frsize
    sw = subprocess.run(['sysctl', '-n', 'vm.swapusage'], capture_output=True, text=True).stdout
    swf = 0.0
    for tok in sw.replace('=', ' ').split():
        pass
    try:
        swf = float(sw.split('free =')[1].split('M')[0]) * 1024 * 1024
    except Exception:
        swf = -1.0
    mp = subprocess.run(['memory_pressure'], capture_output=True, text=True).stdout
    try:
        mfree = float(mp.split('System-wide memory free percentage:')[1].split('%')[0])
    except Exception:
        mfree = -1.0
    return {'disk_free_bytes': free, 'swap_free_bytes': swf, 'mem_free_pct': mfree}


def resource_verdict(r):
    """(ok, reason). Hard stop on disk; true memory exhaustion also stops."""
    if r['disk_free_bytes'] < DISK_MIN_BYTES:
        return False, (f"disk free {r['disk_free_bytes']/1e6:.0f} MB below the "
                       f"{DISK_MIN_BYTES/1e6:.0f} MB floor")
    if 0 <= r['mem_free_pct'] < 2.0 and 0 <= r['swap_free_bytes'] < 50e6:
        return False, (f"memory exhaustion: {r['mem_free_pct']:.1f}% free and "
                       f"{r['swap_free_bytes']/1e6:.0f} MB swap free")
    return True, ''


def load_cohort(repo, subset=None):
    man = [r for r in csv.DictReader(open(os.path.join(
        repo, 'discovery_analysis/v1/inclusion_manifest.tsv')), delimiter='\t')
        if r['in_core'] == '1']
    seqs = {}
    with open(os.path.join(repo, 'DeepLoc2/esm2_probe/deeploc_pairs.tsv')) as f:
        next(f)
        for line in f:
            k, s = line.rstrip('\n').split('\t'); seqs[k] = s
    for r in man:
        assert hashlib.sha256(seqs[r['ACC']].encode()).hexdigest() == r['frozen_sha256'], r['ACC']
    if subset:
        man = man[:subset]
    return man, seqs


def build_batches(acc, seqs, batch_tokens):
    order = sorted(range(len(acc)), key=lambda i: -min(len(seqs[acc[i]]), MAX_RES))
    batches, cur, cmax = [], [], 0
    for i in order:
        L = min(len(seqs[acc[i]]), MAX_RES) + 2
        m = max(cmax, L)
        if cur and m * (len(cur) + 1) > batch_tokens:
            batches.append(cur); cur, cmax = [i], L
        else:
            cur.append(i); cmax = m
    if cur: batches.append(cur)
    flat = [i for b in batches for i in b]
    assert len(flat) == len(set(flat)) == len(acc), 'batching must cover the population exactly'
    return batches


def compat_key(acc, batches, batch_tokens):
    """Everything that must match for a store to be reusable."""
    return {
        'model_checkpoint': CHECKPOINT_ID,
        'extraction_code_sha256': extraction_code_sha(),
        'layers': LAYERS,
        'hidden_dim': HID,
        'store_dtype': STORE_DTYPE,
        'compute_dtype': COMPUTE_DTYPE,
        'pooling_dtype': POOLING_DTYPE,
        'pooling': POOLING_DESC,
        'nterm_positions': NTERM,
        'truncation_max_residues': MAX_RES,
        'n_rows': len(acc),
        'row_order_sha256': hashlib.sha256('\n'.join(acc).encode()).hexdigest(),
        'batch_tokens': batch_tokens,
        'n_batches': len(batches),
        'batch_map_sha256': hashlib.sha256(
            json.dumps([list(map(int, b)) for b in batches]).encode()).hexdigest(),
        'sequence_to_row_sha256': hashlib.sha256(
            json.dumps({a: i for i, a in enumerate(acc)}, sort_keys=True).encode()).hexdigest(),
    }


def store_paths(store_dir):
    return {l: os.path.join(store_dir, f'layer_{l:02d}.{STORE_DTYPE}') for l in LAYERS}


def preallocate(store_dir, n_rows, log=print):
    """Create full-size files up front so a disk shortfall fails immediately, not at 80%."""
    os.makedirs(store_dir, exist_ok=True)
    nbytes = n_rows * 2 * HID * np.dtype(STORE_DTYPE).itemsize
    paths = store_paths(store_dir)
    need = sum(nbytes for p in paths.values() if not os.path.exists(p))
    r = resources(store_dir)
    if need and r['disk_free_bytes'] < need + DISK_MIN_BYTES:
        return False, (f"preallocation needs {need/1e9:.2f} GB plus a "
                       f"{DISK_MIN_BYTES/1e9:.2f} GB floor; only "
                       f"{r['disk_free_bytes']/1e9:.2f} GB free")
    zero = np.zeros(1 << 20, np.uint8).tobytes()
    for l, p in paths.items():
        if os.path.exists(p) and os.path.getsize(p) == nbytes: continue
        with open(p, 'wb') as f:
            w = 0
            while w < nbytes:
                c = min(len(zero), nbytes - w); f.write(zero[:c]); w += c
            f.flush(); os.fsync(f.fileno())
    log(f'  preallocated {len(paths)} layer files x {nbytes/1e6:.1f} MB = '
        f'{len(paths)*nbytes/1e9:.3f} GB')
    return True, ''


def open_store(store_dir, n_rows, mode='r+'):
    return {l: np.memmap(p, dtype=STORE_DTYPE, mode=mode, shape=(n_rows, 2 * HID))
            for l, p in store_paths(store_dir).items()}


def sync_store(store_dir, mm):
    """flush() then fsync() every layer file, BEFORE the checkpoint is advanced."""
    for l in LAYERS: mm[l].flush()
    for p in store_paths(store_dir).values():
        fd = os.open(p, os.O_RDWR)
        try: os.fsync(fd)
        finally: os.close(fd)


def commit_checkpoint(store_dir, key, completed):
    """Atomic: write temp, fsync, rename. Callers MUST sync_store() first."""
    ck = os.path.join(store_dir, 'checkpoint.json')
    tmp = ck + '.tmp'
    with open(tmp, 'w') as f:
        json.dump({'compat': key, 'completed_batches': sorted(int(x) for x in completed)}, f)
        f.flush(); os.fsync(f.fileno())
    os.replace(tmp, ck)
    fd = os.open(store_dir, os.O_RDONLY)
    try: os.fsync(fd)
    finally: os.close(fd)


def read_checkpoint(store_dir, key):
    """(completed_set, status). Reuse ONLY on an exact compatibility match."""
    ck = os.path.join(store_dir, 'checkpoint.json')
    if not os.path.exists(ck): return set(), 'no checkpoint: starting fresh'
    c = json.load(open(ck))
    old = c.get('compat', {})
    diff = [k for k in key if old.get(k) != key[k]] + [k for k in old if k not in key]
    if diff:
        return None, 'REFUSED: compatibility mismatch in ' + ', '.join(sorted(diff))
    return set(c.get('completed_batches', [])), 'verified: resuming'


def layer_hashes(mm):
    out = {}
    for l in LAYERS:
        a = np.asarray(mm[l])
        out[str(l)] = {'whole': hashlib.sha256(a[:, :HID].tobytes()).hexdigest(),
                       'nterm': hashlib.sha256(a[:, HID:].tobytes()).hexdigest(),
                       'both': hashlib.sha256(a.tobytes()).hexdigest()}
        del a
    return out


def extract_into_store(acc, seqs, batches, store_dir, mm, key, completed,
                       max_batches=None, flush_every=10, log=print, device=None):
    """Extract pooled vectors into the disk-backed store, checkpointing atomically.

    Returns (status, info). status is 'complete', 'stopped_max_batches' or
    'stopped_resources'. On a resource stop the store is synced and the checkpoint
    committed first, so no work is lost and nothing is deleted.
    """
    import torch
    from transformers import AutoTokenizer, EsmModel
    todo = [bi for bi in range(len(batches)) if bi not in completed]
    if max_batches is not None: todo = todo[:max_batches]
    if not todo:
        return 'complete', {'new_sequences': 0, 'elapsed_s': 0.0, 'truncated': 0}

    dev = device or torch.device('mps' if torch.backends.mps.is_available() else 'cpu')
    tok = AutoTokenizer.from_pretrained(CHECKPOINT_ID)
    tdt = {'float16': torch.float16, 'float32': torch.float32}[COMPUTE_DTYPE]
    model = EsmModel.from_pretrained(CHECKPOINT_ID, add_pooling_layer=False, torch_dtype=tdt)
    model.eval().to(dev)
    for p in model.parameters(): p.requires_grad_(False)
    log(f'  device={dev} batches_to_do={len(todo)} (store has {len(completed)}/{len(batches)})')

    t0, nseq, ntr, unk = time.time(), 0, 0, tok.unk_token_id
    status, info = 'complete', {}
    try:
        for n, bi in enumerate(todo, 1):
            b = batches[bi]
            ss = [seqs[acc[i]][:MAX_RES] for i in b]
            ntr += sum(1 for i in b if len(seqs[acc[i]]) > MAX_RES)
            enc = tok(ss, return_tensors='pt', padding=True, add_special_tokens=True)
            assert (enc['input_ids'] == unk).sum().item() == 0, 'unexpected <unk>'
            for r, s in enumerate(ss):
                assert int(enc['attention_mask'][r].sum()) == len(s) + 2, 'token/residue mismatch'
            with torch.inference_mode():
                hs = model(input_ids=enc['input_ids'].to(dev),
                           attention_mask=enc['attention_mask'].to(dev),
                           output_hidden_states=True).hidden_states
            for layer in LAYERS:
                h = hs[layer].float()
                for r, s in enumerate(ss):
                    L = len(s); res = h[r, 1:1 + L, :]
                    mm[layer][b[r], :HID] = res.mean(dim=0).cpu().numpy().astype(STORE_DTYPE)
                    mm[layer][b[r], HID:] = res[:min(NTERM, L)].mean(dim=0).cpu().numpy().astype(STORE_DTYPE)
                del h
            del hs
            completed.add(bi); nseq += len(b)

            if n % flush_every == 0 or n == len(todo):
                sync_store(store_dir, mm)                 # flush + fsync DATA ...
                commit_checkpoint(store_dir, key, completed)   # ... THEN mark complete
                r = resources(store_dir)
                ok, why = resource_verdict(r)
                if not ok:
                    status = 'stopped_resources'
                    info = {'reason': why, 'resources': r}
                    log(f'  RESOURCE STOP: {why}')
                    break
            if n % 50 == 0 or n == len(todo):
                el = time.time() - t0
                r = resources(store_dir)
                log(f'  {n}/{len(todo)} batches  {nseq:,} seqs  {el:,.0f}s  {nseq/el:.2f} seq/s  '
                    f'store={len(completed)}/{len(batches)}  disk={r["disk_free_bytes"]/1e9:.2f}GB '
                    f'memfree={r["mem_free_pct"]:.0f}%')
    finally:
        sync_store(store_dir, mm)
        commit_checkpoint(store_dir, key, completed)
        del model
        import gc; gc.collect()
        if dev.type == 'mps':
            try: torch.mps.empty_cache(); torch.mps.synchronize()
            except Exception: pass
        gc.collect()
        log('  released ESM model and device allocations')

    if status == 'complete' and max_batches is not None and len(completed) < len(batches):
        status = 'stopped_max_batches'
    info.update({'new_sequences': nseq, 'elapsed_s': round(time.time() - t0, 1),
                 'truncated': ntr, 'completed_batches': len(completed),
                 'total_batches': len(batches)})
    return status, info
