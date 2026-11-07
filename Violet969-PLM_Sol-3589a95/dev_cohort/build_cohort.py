#!/usr/bin/env python3
"""
Fixed UESolDS development cohort for the ESM-2 comparison.

  TRAIN  : exactly 10,000 sampled without replacement from the 70,030 eligible
           training records (PDB_230301_28958 excluded, per the audit).
  VAL    : all 4,000 published validation proteins.
  TEST   : NOT TOUCHED - reserved for final evaluation.

Stratification (reusing the baseline's source-group definitions verbatim):
  - SG_target records      -> stratum = (center, label)       e.g. ('MCSG','1')
  - every other ID class   -> stratum = (id_class, label)     e.g. ('PDB_serial','1')
Allocation is proportional with largest remainders; records are sorted by ID
before sampling; seed 42. Records are eligible regardless of annotation coverage.

Label convention (verified in the audit): raw '1' = soluble, '0' = insoluble.
  utils/general.py:10          SOLUBILITY = ['0','1']
  datasets/transforms.py:72-87 SolubilityToInt -> SOLUBILITY.index()
"""
import os, sys, json, random, hashlib, collections, warnings
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.metrics import roc_auc_score

warnings.filterwarnings('ignore', message='.*encountered in matmul.*',
                        category=RuntimeWarning)

BASE = sys.argv[1]
OUT  = sys.argv[2]
BL   = os.path.join(BASE, 'baseline_source_vs_sequence')
sys.path.insert(0, BL)
import baseline as B          # reuse id_class / source_center / parse_fasta / load / seq_features

N_TRAIN   = 10000
SEED      = 42
C_FIXED   = 1.0
ESM_LIMIT = 1022              # ESM-2 context 1024 tokens incl. <cls> and <eos>
EVAL_CENTERS = ['EFI', 'JCSG', 'MCSG', 'SECSG', 'SGPP']   # the five eligible in the baseline

def stratum(r):
    """SG_target -> center x label ; everything else -> id_class x label."""
    return (r['center'], r['y']) if r['cls'] == 'SG_target' else (r['cls'], r['y'])

def source_group(r):
    """One flat source-group label per record, for the manifest."""
    return f"SG_target:{r['center']}" if r['cls'] == 'SG_target' else r['cls']

def largest_remainder(sizes, total):
    """Proportional allocation with largest remainders, capped at stratum size.
    Any quota a stratum cannot absorb is redistributed to strata with headroom."""
    N = sum(sizes.values())
    exact = {k: v / N * total for k, v in sizes.items()}
    alloc = {k: min(int(v), sizes[k]) for k, v in exact.items()}
    order = sorted(exact, key=lambda k: (-(exact[k] - int(exact[k])), k))
    while sum(alloc.values()) < total:
        progressed = False
        for k in order:
            if sum(alloc.values()) >= total: break
            if alloc[k] < sizes[k]:
                alloc[k] += 1; progressed = True
        if not progressed:
            raise RuntimeError('cannot reach target: not enough records')
    return alloc

def sha(s): return hashlib.sha256(s.encode()).hexdigest()

def main():
    tr, n_drop = B.load(os.path.join(BASE, 'embedding_dataset', 'Train_dataset.fasta'),
                        B.EXCLUDE_FROM_TRAIN)
    va, _      = B.load(os.path.join(BASE, 'embedding_dataset', 'validation_dataset.fasta'))
    print(f"eligible train pool = {len(tr)} (excluded {n_drop}: {sorted(B.EXCLUDE_FROM_TRAIN)})")
    print(f"published validation = {len(va)}")
    assert len(tr) == 70030, f'expected 70030 eligible training records, got {len(tr)}'

    # The published Train FASTA reuses some identifiers across DIFFERENT records
    # (169 IDs over 374 records; distinct sequences = alternative constructs of the
    # same structural-genomics target; no label conflicts). So original_id is NOT a
    # primary key. Every record gets its 0-based position in the published file, and
    # record_uid = "<split>:<index>" is the manifest key.
    for i, r in enumerate(tr): r['idx'] = i
    for i, r in enumerate(va): r['idx'] = i
    _c = collections.Counter(r['id'] for r in tr)
    _rep = {k: v for k, v in _c.items() if v > 1}
    print(f"train: unique IDs={len(_c)} | IDs reused across records={len(_rep)}"
          f" | extra records={sum(v-1 for v in _rep.values())}")
    print(f"validation: unique IDs={len({r['id'] for r in va})} of {len(va)}")

    # ---------------- 1. stratified sample ----------------
    strata = collections.defaultdict(list)
    for r in tr: strata[stratum(r)].append(r)
    sizes = {k: len(v) for k, v in strata.items()}
    alloc = largest_remainder(sizes, N_TRAIN)
    rng = random.Random(SEED)
    sample = []
    for k in sorted(strata):                       # deterministic stratum order
        pool = sorted(strata[k], key=lambda r: r['id'])   # sort BEFORE sampling
        n = alloc[k]
        if n: sample += rng.sample(pool, n)
    assert len(sample) == N_TRAIN, len(sample)
    assert len({r['idx'] for r in sample}) == N_TRAIN, 'duplicate records drawn'
    n_shared_id = N_TRAIN - len({r['id'] for r in sample})
    print(f"sampled {N_TRAIN} distinct records; of these, {n_shared_id} share an "
          f"original_id with another sampled record")
    sample.sort(key=lambda r: (r['id'], r['idx']))

    print(f"\n### sampling allocation ({len(strata)} strata -> {N_TRAIN})")
    print(f"{'stratum':<26}{'pool':>8}{'pool%':>8}{'alloc':>7}{'alloc%':>8}")
    for k in sorted(strata, key=lambda k: -sizes[k]):
        print(f"{k[0]+' / label '+str(k[1]):<26}{sizes[k]:>8}{sizes[k]/len(tr)*100:>7.2f}%"
              f"{alloc[k]:>7}{alloc[k]/N_TRAIN*100:>7.2f}%")
    print(f"{'TOTAL':<26}{len(tr):>8}{100.0:>7.2f}%{sum(alloc.values()):>7}{100.0:>7.2f}%")

    # ---------------- 2. manifest ----------------
    rows = []
    for split, recs in [('train', sample), ('validation', va)]:
        for r in recs:
            rows.append({'record_uid': f"{split}:{r['idx']}",
                         'fasta_record_index': r['idx'],
                         'original_id': r['id'], 'published_split': split,
                         'label': r['y'], 'source_group': source_group(r),
                         'id_class': r['cls'], 'center': r['center'],
                         'seq_len': r['len'], 'sha256': sha(r['seq'])})
    cols = ['record_uid','original_id','published_split','fasta_record_index','label',
            'source_group','id_class','center','seq_len','sha256']
    os.makedirs(OUT, exist_ok=True)
    man = os.path.join(OUT, 'dev_cohort_manifest.tsv')
    with open(man, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows: f.write('\t'.join(str(r[c]) for c in cols) + '\n')

    # hash uniqueness checks
    print("\n### sequence-hash uniqueness")
    byspl = collections.defaultdict(list)
    for r in rows: byspl[r['published_split']].append(r['sha256'])
    for s, h in byspl.items():
        dup = len(h) - len(set(h))
        print(f"  within {s:<11} n={len(h):<6} unique={len(set(h)):<6} duplicates={dup}")
    inter = set(byspl['train']) & set(byspl['validation'])
    print(f"  across train<->validation shared hashes = {len(inter)}")
    if inter:
        bad = [r for r in rows if r['sha256'] in inter]
        for r in sorted(bad, key=lambda r: r['sha256'])[:10]: print('    ', r)
    allh = [r['sha256'] for r in rows]
    print(f"  cohort-wide n={len(allh)} unique={len(set(allh))} duplicates={len(allh)-len(set(allh))}")

    # counts by split x label
    print("\n### counts by split and label")
    print(f"{'split':<13}{'n':>7}{'insoluble(0)':>14}{'soluble(1)':>12}{'%soluble':>10}")
    for s in ['train','validation']:
        v = [r for r in rows if r['published_split'] == s]
        n1 = sum(r['label'] for r in v); n0 = len(v) - n1
        print(f"{s:<13}{len(v):>7}{n0:>14}{n1:>12}{n1/len(v)*100:>9.1f}%")
    n1 = sum(r['label'] for r in rows); n0 = len(rows) - n1
    print(f"{'COHORT':<13}{len(rows):>7}{n0:>14}{n1:>12}{n1/len(rows)*100:>9.1f}%")

    # length / ESM truncation planning
    print(f"\n### sequences longer than {ESM_LIMIT} residues (ESM-2 truncation planning)")
    print(f"{'split':<13}{'n':>7}{'>1022':>8}{'%':>8}{'max_len':>9}{'median':>8}")
    for s in ['train','validation']:
        v = [r for r in rows if r['published_split'] == s]
        L = [r['seq_len'] for r in v]; over = sum(1 for x in L if x > ESM_LIMIT)
        print(f"{s:<13}{len(v):>7}{over:>8}{over/len(v)*100:>7.2f}%"
              f"{max(L):>9}{int(np.median(L)):>8}")
    L = [r['seq_len'] for r in rows]; over = sum(1 for x in L if x > ESM_LIMIT)
    print(f"{'COHORT':<13}{len(rows):>7}{over:>8}{over/len(rows)*100:>7.2f}%{max(L):>9}{int(np.median(L)):>8}")
    byl = collections.Counter(r['label'] for r in rows if r['seq_len'] > ESM_LIMIT)
    print(f"  of those: insoluble={byl[0]} soluble={byl[1]}")
    print(f"  total residues after truncation at {ESM_LIMIT}: {sum(min(x,ESM_LIMIT) for x in L):,}")

    # ---------------- 3. refit A/B/C on the 10,000 ----------------
    ytr = np.array([r['y'] for r in sample]); yva = np.array([r['y'] for r in va])
    enc = OneHotEncoder(handle_unknown='ignore', sparse_output=False)
    Ctr = np.array([[r['cls'], r['center']] for r in sample], dtype=object)
    Cva = np.array([[r['cls'], r['center']] for r in va], dtype=object)
    enc.fit(Ctr)                                           # fit on the 10,000 only
    Atr, Ava = enc.transform(Ctr), enc.transform(Cva)
    unseen = int((Ava.sum(axis=1) < 2).sum())
    print(f"\nsource one-hot dims={Atr.shape[1]}  validation rows with >=1 unseen category: {unseen}")

    sc = StandardScaler().fit(B.seq_features(sample))      # fit on the 10,000 only
    Btr, Bva = sc.transform(B.seq_features(sample)), sc.transform(B.seq_features(va))
    Dtr, Dva = np.hstack([Atr, Btr]), np.hstack([Ava, Bva])
    for nm, M in [('Atr',Atr),('Ava',Ava),('Btr',Btr),('Bva',Bva),('Dtr',Dtr),('Dva',Dva)]:
        assert np.isfinite(M).all(), f'non-finite values in {nm}'

    models, scores = {}, {}
    for name, Xtr, Xva in [('A_source_only', Atr, Ava),
                           ('B_sequence_only', Btr, Bva),
                           ('C_source_plus_sequence', Dtr, Dva)]:
        m = LogisticRegression(penalty='l2', C=C_FIXED, solver='lbfgs',
                               max_iter=5000, random_state=SEED).fit(Xtr, ytr)
        models[name], scores[name] = m, m.predict_proba(Xva)[:, 1]
        print(f"fitted {name}: {Xtr.shape[1]} features, n_iter={m.n_iter_[0]}")

    def auroc(y, s):
        return None if len(np.unique(y)) < 2 else float(roc_auc_score(y, s))

    results = []
    def ev(mask, gtype, label):
        y = yva[mask]; n1 = int(y.sum()); n0 = int(len(y) - n1)
        row = {'group_type': gtype, 'group': label, 'n': int(len(y)),
               'n_insoluble_0': n0, 'n_soluble_1': n1,
               'auroc_defined': not (n0 == 0 or n1 == 0)}
        for nm in models:
            a = auroc(y, scores[nm][mask])
            row[nm] = None if a is None else round(a, 4)
        results.append(row); return row

    ev(np.ones(len(va), bool), 'overall', 'ALL_validation')
    ev(np.array([r['cls'] == 'SG_target' for r in va]), 'SG_target_only', 'SG_target')
    for c in EVAL_CENTERS:
        ev(np.array([r['cls'] == 'SG_target' and r['center'] == c for r in va]), 'SG_center', c)

    # unweighted mean across the five centers - each center counts equally
    cen = [r for r in results if r['group_type'] == 'SG_center']
    mean_row = {'group_type': 'macro_mean', 'group': f'unweighted_mean_of_{len(cen)}_centers',
                'n': sum(r['n'] for r in cen),
                'n_insoluble_0': sum(r['n_insoluble_0'] for r in cen),
                'n_soluble_1': sum(r['n_soluble_1'] for r in cen), 'auroc_defined': True}
    for nm in models:
        vals = [r[nm] for r in cen if r[nm] is not None]
        mean_row[nm] = round(float(np.mean(vals)), 4) if vals else None
    results.append(mean_row)

    cols2 = ['group_type','group','n','n_insoluble_0','n_soluble_1',
             'A_source_only','B_sequence_only','C_source_plus_sequence','auroc_defined']
    res = os.path.join(OUT, 'dev_cohort_validation_auroc.tsv')
    with open(res, 'w') as f:
        f.write('\t'.join(cols2) + '\n')
        for r in results:
            f.write('\t'.join('undefined_single_class' if r.get(c) is None and c.startswith(('A_','B_','C_'))
                              else str(r.get(c)) for c in cols2) + '\n')

    print("\n" + "=" * 100)
    print("VALIDATION AUROC - models refit on the 10,000-protein development cohort")
    print("=" * 100)
    print(f"{'group':<34}{'n':>6}{'insol':>7}{'sol':>6}{'A_source':>11}{'B_seq':>9}{'C_both':>9}")
    print("-" * 100)
    for r in results:
        f_ = lambda v: 'undefined' if v is None else f'{v:.4f}'
        pre = '  ' if r['group_type'] == 'SG_center' else ''
        print(f"{pre+r['group']:<34}{r['n']:>6}{r['n_insoluble_0']:>7}{r['n_soluble_1']:>6}"
              f"{f_(r['A_source_only']):>11}{f_(r['B_sequence_only']):>9}{f_(r['C_source_plus_sequence']):>9}")
    print("-" * 100)

    json.dump({'results': results, 'allocation': {f'{k[0]}|{k[1]}': v for k, v in alloc.items()},
               'seed': SEED, 'C': C_FIXED, 'n_train': N_TRAIN, 'n_val': len(va),
               'eligible_pool': len(tr), 'excluded': sorted(B.EXCLUDE_FROM_TRAIN),
               'esm_limit': ESM_LIMIT, 'unseen_category_val_rows': unseen,
               'train_val_shared_hashes': len(inter)},
              open(os.path.join(OUT, 'dev_cohort_results.json'), 'w'), indent=1)
    with open(os.path.join(OUT, 'dev_cohort_train_10000.fasta'), 'w') as f:
        for r in sample: f.write(f">{r['id']} A-{r['y']}\n{r['seq']}\n")
    print(f"\nwrote {man}\nwrote {res}\nwrote {os.path.join(OUT,'dev_cohort_results.json')}")
    print(f"wrote {os.path.join(OUT,'dev_cohort_train_10000.fasta')}")

if __name__ == '__main__':
    main()
