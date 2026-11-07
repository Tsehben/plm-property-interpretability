#!/usr/bin/env python3
"""
DeepLoc 2.0: frozen ESM-2-650M layer-8 probes under the EXISTING four-fold protocol.

Development = published partitions 1-4; each held out once, trained on the other three.
Partition 0 and hpa_testset.csv stay reserved. L2 logistic regression, fixed C=1,
scalers fit INSIDE each training fold. Sequence alone is the model input.

Models
  matched_C            MATCHED baseline C: AA fractions + log(length) over the SAME
                       visible residues (first 1022) + AA fractions over positions 1-30
  fulllen_C            original full-length baseline C, for continuity
  a_ESM_whole          mean of hidden_states[8] over all visible residue tokens
  b_ESM_whole_nterm    a + mean of hidden_states[8] over positions 1-30 (CONTEXTUAL:
                       those positions attend to the whole visible sequence)
  c_matchedC_plus_ESM  matched_C + both ESM vectors
"""
import os, sys, csv, json, collections, warnings
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score
warnings.filterwarnings('ignore', message='.*encountered in matmul.*', category=RuntimeWarning)

D, COMMON, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, COMMON)
from esm_extract import gather, MAX_RESIDUES

STD_AA = 'ACDEFGHIKLMNPQRSTVWY'
C_FIXED, SEED, NTERM = 1.0, 42, 30
FOLDS = ['1', '2', '3', '4']

def fracs(s):
    n = len(s)
    if n == 0: return [0.0] * 20
    c = collections.Counter(s)
    return [c.get(a, 0) / n for a in STD_AA]

man = [r for r in csv.DictReader(open(os.path.join(D, 'baseline',
       'deeploc_binary_manifest.tsv')), delimiter='\t') if r['split_role'] == 'development_P1_4']
seqs = {}
with open(os.path.join(D, 'esm2_probe', 'deeploc_pairs.tsv')) as f:
    next(f)
    for line in f:
        k, s = line.rstrip('\n').split('\t'); seqs[k] = s

K, W, N = gather(os.path.join(D, 'esm2_probe', 'emb'), 'deeploc', nterm=True)
pos = {k: i for i, k in enumerate(K)}
print(f'manifest dev={len(man):,}  embeddings={len(K):,}  dim={W.shape[1]}')
missing = [r['ACC'] for r in man if r['ACC'] not in pos]
assert not missing, f'{len(missing)} manifest records have no embedding: {missing[:5]}'
print('alignment: every development record has an embedding; 0 missing')

rec = []
for r in man:
    s = seqs[r['ACC']]
    rec.append({'ACC': r['ACC'], 'y': int(r['target']), 'p': r['partition'],
                'single': r['has_additional_localization'] == '0', 'seq': s,
                'w': W[pos[r['ACC']]], 'n': N[pos[r['ACC']]],
                'trunc': len(s) > MAX_RESIDUES})
y = np.array([r['y'] for r in rec]); part = np.array([r['p'] for r in rec])
print(f'development={len(rec):,}  truncated at {MAX_RESIDUES}: {sum(r["trunc"] for r in rec):,}')

def f_matched(R):
    return np.array([fracs(r['seq'][:MAX_RESIDUES])
                     + [np.log(min(len(r['seq']), MAX_RESIDUES))]
                     + fracs(r['seq'][:NTERM]) for r in R])
def f_full(R):
    return np.array([fracs(r['seq']) + [np.log(len(r['seq']))]
                     + fracs(r['seq'][:NTERM]) for r in R])
def f_w(R):  return np.stack([r['w'] for r in R])
def f_wn(R): return np.hstack([np.stack([r['w'] for r in R]), np.stack([r['n'] for r in R])])
def f_c(R):  return np.hstack([f_matched(R), f_w(R), np.stack([r['n'] for r in R])])

MODELS = {'matched_C': f_matched, 'fulllen_C': f_full, 'a_ESM_whole': f_w,
          'b_ESM_whole_nterm': f_wn, 'c_matchedC_plus_ESM': f_c}
oof = {m: np.full(len(rec), np.nan) for m in MODELS}

for p in FOLDS:
    te = np.where(part == p)[0]; tr = np.where(part != p)[0]
    Rtr, Rte = [rec[i] for i in tr], [rec[i] for i in te]
    for m, fn in MODELS.items():
        Xtr, Xte = fn(Rtr), fn(Rte)
        sc = StandardScaler().fit(Xtr)                 # fold-local
        Xtr, Xte = sc.transform(Xtr), sc.transform(Xte)
        assert np.isfinite(Xtr).all() and np.isfinite(Xte).all()
        clf = LogisticRegression(penalty='l2', C=C_FIXED, solver='lbfgs',
                                 max_iter=5000, random_state=SEED).fit(Xtr, y[tr])
        oof[m][te] = clf.predict_proba(Xte)[:, 1]
    print(f'  fold P{p}: train={len(tr):,} test={len(te):,} pos={int(y[te].sum()):,}')
for m in MODELS: assert not np.isnan(oof[m]).any()

with open(os.path.join(OUT, 'deeploc_esm_oof_predictions.tsv'), 'w') as f:
    f.write('ACC\tpartition\ttarget\thas_additional_localization\t' + '\t'.join(MODELS) + '\n')
    for i, r in enumerate(rec):
        f.write(f"{r['ACC']}\t{r['p']}\t{r['y']}\t{0 if r['single'] else 1}\t"
                + '\t'.join(f'{oof[m][i]:.6f}' for m in MODELS) + '\n')

results = []
def block(title, mask):
    yy = y[mask]
    print('\n' + '=' * 112); print(title); print('=' * 112)
    print(f'  n={len(yy):,}  Cytoplasm(0)={int(len(yy)-yy.sum()):,}  '
          f'Extracellular(1)={int(yy.sum()):,}  prevalence={yy.mean():.4f}')
    print(f"\n{'model':<24}" + ''.join(f'{"P"+p:>9}' for p in FOLDS) + f'{"MEAN":>9}{"|":>3}'
          + ''.join(f'{"P"+p:>9}' for p in FOLDS) + f'{"MEAN":>9}')
    print(f"{'':<24}{'------- AUROC -------':^39}{'|':>3}{'--- average precision ---':^39}")
    print('-' * 112)
    store = {}
    for m in MODELS:
        au, ap = [], []
        for p in FOLDS:
            fm = mask & (part == p)
            au.append(roc_auc_score(y[fm], oof[m][fm]))
            ap.append(average_precision_score(y[fm], oof[m][fm]))
        store[m] = (au, ap)
        print(f'{m:<24}' + ''.join(f'{v:>9.4f}' for v in au) + f'{np.mean(au):>9.4f}{"|":>3}'
              + ''.join(f'{v:>9.4f}' for v in ap) + f'{np.mean(ap):>9.4f}')
        results.append({'subset': title, 'model': m, 'n': int(len(yy)),
                        'n_pos': int(yy.sum()), 'n_neg': int(len(yy) - yy.sum()),
                        'prevalence': round(float(yy.mean()), 4),
                        **{f'auroc_P{p}': round(float(v), 4) for p, v in zip(FOLDS, au)},
                        'auroc_mean': round(float(np.mean(au)), 4),
                        **{f'ap_P{p}': round(float(v), 4) for p, v in zip(FOLDS, ap)},
                        'ap_mean': round(float(np.mean(ap)), 4)})
    print('-' * 112)
    for m in ['a_ESM_whole', 'b_ESM_whole_nterm', 'c_matchedC_plus_ESM']:
        dau = [store[m][0][i] - store['matched_C'][0][i] for i in range(4)]
        dap = [store[m][1][i] - store['matched_C'][1][i] for i in range(4)]
        lbl = f'delta {m[:1]} - matched_C'
        print(f'{lbl:<24}' + ''.join(f'{v:>+9.4f}' for v in dau)
              + f'{np.mean(dau):>+9.4f}{"|":>3}' + ''.join(f'{v:>+9.4f}' for v in dap)
              + f'{np.mean(dap):>+9.4f}')
        results.append({'subset': title, 'model': f'delta_{m}_minus_matched_C',
                        'n': int(len(yy)), 'n_pos': int(yy.sum()),
                        'n_neg': int(len(yy) - yy.sum()),
                        'prevalence': round(float(yy.mean()), 4),
                        **{f'auroc_P{p}': round(float(v), 4) for p, v in zip(FOLDS, dau)},
                        'auroc_mean': round(float(np.mean(dau)), 4),
                        **{f'ap_P{p}': round(float(v), 4) for p, v in zip(FOLDS, dap)},
                        'ap_mean': round(float(np.mean(dap)), 4)})
    print(f"   folds where c beats matched_C: AUROC "
          f"{sum(1 for i in range(4) if store['c_matchedC_plus_ESM'][0][i] > store['matched_C'][0][i])}/4,"
          f" AP {sum(1 for i in range(4) if store['c_matchedC_plus_ESM'][1][i] > store['matched_C'][1][i])}/4")

block('ALL development proteins (P1-4)', np.ones(len(rec), bool))
block('SENSITIVITY: no additional localization annotation',
      np.array([r['single'] for r in rec]))

cols = (['subset', 'model', 'n', 'n_neg', 'n_pos', 'prevalence']
        + [f'auroc_P{p}' for p in FOLDS] + ['auroc_mean']
        + [f'ap_P{p}' for p in FOLDS] + ['ap_mean'])
with open(os.path.join(OUT, 'deeploc_esm_results.tsv'), 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in results: f.write('\t'.join(str(r[c]) for c in cols) + '\n')
json.dump({'results': results, 'C': C_FIXED, 'layer': 8, 'n_dev': len(rec),
           'n_truncated': int(sum(r['trunc'] for r in rec))},
          open(os.path.join(OUT, 'deeploc_esm_results.json'), 'w'), indent=1)
print(f"\nwrote {os.path.join(OUT,'deeploc_esm_results.tsv')}")
print(f"wrote {os.path.join(OUT,'deeploc_esm_oof_predictions.tsv')}")
