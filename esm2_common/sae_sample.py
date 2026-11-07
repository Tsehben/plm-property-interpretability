#!/usr/bin/env python3
"""Sample 64 pilot proteins per dataset for the layer-8 InterProt SAE pilot.

Reproducible: seed 42, stratified by label x length quartile (2 x 4 = 8 strata,
8 records each) so both labels and a range of lengths are covered.

Reserved and never read: UESolDS test_dataset.fasta, DeepLoc partition 0,
hpa_testset.csv.  Eligible pools are restricted to UESolDS *training* records
and DeepLoc *development* (P1-4) records only.
"""
import csv, hashlib, json, os, sys
import numpy as np

# Repo root: the parent of esm2_common/.  Override with MIT6393_ROOT if the
# checkout lives elsewhere.
ROOT = os.environ.get('MIT6393_ROOT',
                      os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SEED, N_PILOT, MAX_RESIDUES = 42, 64, 1022


def read_tsv(p):
    with open(p, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def read_pairs(p):
    d = {}
    with open(p, newline='') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            d[r['key']] = r['sequence']
    return d


def stratified(records, keyfn_label, keyfn_len, n=N_PILOT):
    """8 strata = 2 labels x 4 length quartiles, n/8 records each."""
    lens = np.array([keyfn_len(r) for r in records])
    qs = np.quantile(lens, [.25, .5, .75])
    per = n // 8
    rng = np.random.default_rng(SEED)
    out = []
    for lab in ['0', '1']:
        for q in range(4):
            lo = -np.inf if q == 0 else qs[q - 1]
            hi = np.inf if q == 3 else qs[q]
            pool = [r for r in records
                    if keyfn_label(r) == lab
                    and (keyfn_len(r) > lo or q == 0) and keyfn_len(r) <= hi]
            pool.sort(key=lambda r: r['_key'])          # deterministic order
            assert len(pool) >= per, f'stratum {lab}/{q} too small: {len(pool)}'
            idx = rng.choice(len(pool), size=per, replace=False)
            for i in sorted(idx):
                r = dict(pool[i]); r['_label'] = lab; r['_quartile'] = q
                out.append(r)
    out.sort(key=lambda r: r['_key'])
    return out, qs


def emit(rows, seqs, out_path, dataset):
    with open(out_path, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['key', 'label', 'seq_len', 'visible_len', 'truncated',
                    'length_quartile', 'sha256', 'sequence'])
        for r in rows:
            s = seqs[r['_key']]
            assert hashlib.sha256(s.encode()).hexdigest() == r['_sha'], \
                f'sequence hash mismatch for {r["_key"]}'
            vis = min(len(s), MAX_RESIDUES)
            w.writerow([r['_key'], r['_label'], len(s), vis,
                        int(len(s) > MAX_RESIDUES), r['_quartile'], r['_sha'], s])
    print(f'{dataset}: wrote {len(rows)} -> {out_path}')
    L = [len(seqs[r['_key']]) for r in rows]
    tr = sum(1 for x in L if x > MAX_RESIDUES)
    lab = {l: sum(1 for r in rows if r['_label'] == l) for l in '01'}
    print(f'  label 0={lab["0"]}  label 1={lab["1"]}   '
          f'len min={min(L)} median={int(np.median(L))} max={max(L)}  truncated={tr}')


def main():
    # ---------------- UESolDS: training records of the frozen cohort --------
    man = read_tsv(f'{ROOT}/Violet969-PLM_Sol-3589a95/dev_cohort/dev_cohort_manifest.tsv')
    seqs = read_pairs(f'{ROOT}/Violet969-PLM_Sol-3589a95/esm2_probe/uesolds_pairs.tsv')
    recs = [dict(r, _key=r['record_uid'], _sha=r['sha256'])
            for r in man if r['published_split'] == 'train']
    assert len(recs) == 10000, len(recs)
    rows, qs = stratified(recs, lambda r: r['label'], lambda r: int(r['seq_len']))
    print(f'UESolDS length quartiles: {qs.tolist()}')
    emit(rows, seqs, f'{ROOT}/Violet969-PLM_Sol-3589a95/sae_pilot/sae_pilot_sample.tsv', 'UESolDS')

    # ---------------- DeepLoc: development partitions 1-4 only --------------
    man = read_tsv(f'{ROOT}/DeepLoc2/baseline/deeploc_binary_manifest.tsv')
    seqs = read_pairs(f'{ROOT}/DeepLoc2/esm2_probe/deeploc_pairs.tsv')
    recs = [dict(r, _key=r['ACC'], _sha=r['sha256'])
            for r in man if r['split_role'] == 'development_P1_4']
    assert len(recs) == 10139, len(recs)
    assert all(r['partition'] != '0' for r in recs), 'partition 0 leaked into pool'
    rows, qs = stratified(recs, lambda r: r['target'], lambda r: int(r['len']))
    print(f'DeepLoc length quartiles: {qs.tolist()}')
    emit(rows, seqs, f'{ROOT}/DeepLoc2/sae_pilot/sae_pilot_sample.tsv', 'DeepLoc')


if __name__ == '__main__':
    main()
