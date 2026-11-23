#!/usr/bin/env python3
"""Verify checkpoint/resume on a small deterministic subset BEFORE the full run.

A: uninterrupted extraction of the subset.
B: same subset extracted in three interrupted stages, resuming from the checkpoint.
PASS requires B to equal A byte-for-byte on all 33 layers and both poolings.
Also checks that the compatibility key refuses an incompatible store.

Usage: test_checkpoint_resume.py <repo_root> <scratch_dir> [--subset N] [--batch-tokens T]
"""
import os, sys, json, shutil, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import sweep_lib as S

REPO, SCRATCH = sys.argv[1], sys.argv[2]
def arg(f, d): return type(d)(sys.argv[sys.argv.index(f) + 1]) if f in sys.argv else d
SUB, BT = arg('--subset', 48), arg('--batch-tokens', 6144)
log = lambda *a: print(*a, flush=True)

man, seqs = S.load_cohort(REPO, subset=SUB)
acc = [r['ACC'] for r in man]
batches = S.build_batches(acc, seqs, BT)
key = S.compat_key(acc, batches, BT)
log(f'subset={len(acc)} sequences  batches={len(batches)}  batch_tokens={BT}')
log(f'compat key fields: {sorted(key)}')

def fresh(d):
    if os.path.exists(d): shutil.rmtree(d)
    os.makedirs(d, exist_ok=True)

def run(d, stages):
    fresh(d)
    ok, why = S.preallocate(d, len(acc), log=log)
    assert ok, why
    comp, st = S.read_checkpoint(d, key)
    assert comp == set(), st
    for k, mx in enumerate(stages, 1):
        mm = S.open_store(d, len(acc))
        comp, st = S.read_checkpoint(d, key)
        assert comp is not None, st
        log(f'  stage {k}: {st}; resuming with {len(comp)} batches done, max_batches={mx}')
        status, info = S.extract_into_store(acc, seqs, batches, d, mm, key, comp,
                                            max_batches=mx, flush_every=2, log=log)
        log(f'  stage {k} -> {status} {info}')
        del mm                                    # close memmaps, as a restart would
    mm = S.open_store(d, len(acc), mode='r')
    h = S.layer_hashes(mm)
    comp, _ = S.read_checkpoint(d, key)
    del mm
    return h, comp

log('\n=== RUN A: uninterrupted ===')
hA, cA = run(os.path.join(SCRATCH, 'storeA'), [None])
log('\n=== RUN B: interrupted at 2 batches, then 3, then to completion ===')
hB, cB = run(os.path.join(SCRATCH, 'storeB'), [2, 3, None])

log('\n=== COMPARISON ===')
assert cA == cB == set(range(len(batches))), 'both runs must complete every batch'
mismatch = [(l, p) for l in hA for p in ('whole', 'nterm', 'both') if hA[l][p] != hB[l][p]]
for l in sorted(hA, key=int)[:3] + sorted(hA, key=int)[-2:]:
    log(f'  layer {l:>2}: A.both={hA[l]["both"][:16]}  B.both={hB[l]["both"][:16]}  '
        f'{"MATCH" if hA[l]["both"]==hB[l]["both"] else "MISMATCH"}')
log(f'  ... {len(hA)} layers checked, {len(mismatch)} mismatches')

log('\n=== COMPATIBILITY REFUSAL CHECKS ===')
refus = []
for field, bad in [('model_checkpoint', 'facebook/esm2_t12_35M_UR50D'),
                   ('extraction_code_sha256', '0' * 64),
                   ('store_dtype', 'float32'),
                   ('pooling', 'mean over all tokens including <cls>'),
                   ('truncation_max_residues', 512),
                   ('nterm_positions', 15),
                   ('batch_map_sha256', '1' * 64),
                   ('sequence_to_row_sha256', '2' * 64),
                   ('row_order_sha256', '3' * 64)]:
    k2 = dict(key); k2[field] = bad
    c, st = S.read_checkpoint(os.path.join(SCRATCH, 'storeA'), k2)
    refus.append((field, c is None, st))
    log(f'  mutate {field:<28} -> {"REFUSED" if c is None else "ACCEPTED (BAD)"}')
all_ref = all(r[1] for r in refus)

ok = (not mismatch) and all_ref and cA == cB
json.dump({'subset_n': len(acc), 'batches': len(batches), 'batch_tokens': BT,
           'runA_stages': ['uninterrupted'], 'runB_stages': [2, 3, 'to completion'],
           'layers_compared': len(hA), 'hash_mismatches': len(mismatch),
           'byte_identical': not mismatch,
           'compat_refusal_fields': [r[0] for r in refus],
           'all_incompatible_stores_refused': all_ref,
           'PASS': bool(ok),
           'layer_hashes_runA': hA},
          open(os.path.join(SCRATCH, 'resume_test_report.json'), 'w'), indent=1)
log(f'\nRESULT: {"PASS" if ok else "FAIL"}  '
    f'(byte-identical={not mismatch}, refusals_all_ok={all_ref})')
sys.exit(0 if ok else 1)
