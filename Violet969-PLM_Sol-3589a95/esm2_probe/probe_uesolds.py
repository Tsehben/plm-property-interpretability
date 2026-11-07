#!/usr/bin/env python3
"""
UESolDS: frozen ESM-2-650M layer-8 probe vs composition baselines.

Fixed cohort: 10,000 training + 4,000 validation records (manifest keyed on record_uid).
Test split stays reserved. L2 logistic regression, fixed C=1, scaling fit on TRAINING ONLY.

Models
  ESM_whole            mean over residue tokens of hidden_states[8], whole visible input
  comp_len_matched     20 AA fractions + log(length) over the SAME visible residues
                       (first 1022) - the fair comparison
  comp_len_fulllength  20 AA fractions + log(length) over the complete sequence
                       - the original baseline, for continuity
"""
import os, sys, csv, json, collections, warnings
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
warnings.filterwarnings('ignore', message='.*encountered in matmul.*', category=RuntimeWarning)

U, COMMON, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, COMMON)
from esm_extract import gather, MAX_RESIDUES

STD_AA = 'ACDEFGHIKLMNPQRSTVWY'
C_FIXED, SEED = 1.0, 42
CENTERS = ['EFI', 'JCSG', 'MCSG', 'SECSG', 'SGPP']   # the five eligible from the baseline

def fracs(s):
    n = len(s)
    c = collections.Counter(s)
    return [c.get(a, 0) / n for a in STD_AA] + [np.log(n)]

man = list(csv.DictReader(open(os.path.join(U, 'dev_cohort', 'dev_cohort_manifest.tsv')),
                          delimiter='\t'))
seqs = {}
with open(os.path.join(U, 'esm2_probe', 'uesolds_pairs.tsv')) as f:
    next(f)
    for line in f:
        k, s = line.rstrip('\n').split('\t'); seqs[k] = s

K, W, _ = gather(os.path.join(U, 'esm2_probe', 'emb'), 'uesolds')
pos = {k: i for i, k in enumerate(K)}
print(f'manifest={len(man):,}  embeddings={len(K):,}  dim={W.shape[1]}')
missing = [r['record_uid'] for r in man if r['record_uid'] not in pos]
assert not missing, f'{len(missing)} manifest records have no embedding: {missing[:5]}'
print('alignment: every manifest record has an embedding; 0 missing')

rows = []
for r in man:
    s = seqs[r['record_uid']]
    rows.append({'uid': r['record_uid'], 'split': r['published_split'],
                 'y': int(r['label']), 'center': r['center'], 'cls': r['id_class'],
                 'seq': s, 'emb': W[pos[r['record_uid']]],
                 'trunc': len(s) > MAX_RESIDUES})
tr = [r for r in rows if r['split'] == 'train']
va = [r for r in rows if r['split'] == 'validation']
print(f'train={len(tr):,}  validation={len(va):,}  '
      f'truncated: train={sum(r["trunc"] for r in tr)} val={sum(r["trunc"] for r in va)}')

FEATS = {
    'ESM_whole':           lambda R: np.stack([r['emb'] for r in R]),
    'comp_len_matched':    lambda R: np.array([fracs(r['seq'][:MAX_RESIDUES]) for r in R]),
    'comp_len_fulllength': lambda R: np.array([fracs(r['seq']) for r in R]),
}
ytr = np.array([r['y'] for r in tr]); yva = np.array([r['y'] for r in va])
scores = {}
for name, fn in FEATS.items():
    Xtr, Xva = fn(tr), fn(va)
    sc = StandardScaler().fit(Xtr)                      # TRAINING ONLY
    Xtr, Xva = sc.transform(Xtr), sc.transform(Xva)
    assert np.isfinite(Xtr).all() and np.isfinite(Xva).all()
    m = LogisticRegression(penalty='l2', C=C_FIXED, solver='lbfgs',
                           max_iter=5000, random_state=SEED).fit(Xtr, ytr)
    scores[name] = m.predict_proba(Xva)[:, 1]
    print(f'fitted {name}: {Xtr.shape[1]} features, n_iter={m.n_iter_[0]}')

results = []
print('\n' + '=' * 104)
print('UESolDS VALIDATION AUROC - frozen ESM-2-650M layer 8 vs composition')
print('=' * 104)
hdr = f"{'model':<22}{'pooled':>9}" + ''.join(f'{c:>9}' for c in CENTERS) + f'{"MEAN5":>9}'
print(hdr); print('-' * len(hdr))
cmask = {c: np.array([r['cls'] == 'SG_target' and r['center'] == c for r in va])
         for c in CENTERS}
for name in FEATS:
    pooled = roc_auc_score(yva, scores[name])
    per = [roc_auc_score(yva[cmask[c]], scores[name][cmask[c]]) for c in CENTERS]
    print(f'{name:<22}{pooled:>9.4f}' + ''.join(f'{v:>9.4f}' for v in per)
          + f'{np.mean(per):>9.4f}')
    results.append({'model': name, 'pooled_auroc': round(float(pooled), 4),
                    **{f'auroc_{c}': round(float(v), 4) for c, v in zip(CENTERS, per)},
                    'unweighted_mean_5_centers': round(float(np.mean(per)), 4)})
print('-' * len(hdr))
cnt = []
for c in CENTERS:
    yy = yva[cmask[c]]
    cnt.append(f'{int(yy.sum())}/{int(len(yy) - yy.sum())}')
print(f"{'sol/insol per center':<22}{'':>9}" + ''.join(f'{x:>9}' for x in cnt))

with open(os.path.join(OUT, 'uesolds_esm_results.tsv'), 'w') as f:
    cols = ['model', 'pooled_auroc'] + [f'auroc_{c}' for c in CENTERS] + ['unweighted_mean_5_centers']
    f.write('\t'.join(cols) + '\n')
    for r in results: f.write('\t'.join(str(r[c]) for c in cols) + '\n')
with open(os.path.join(OUT, 'uesolds_esm_val_predictions.tsv'), 'w') as f:
    f.write('record_uid\tcenter\tid_class\tlabel\t' + '\t'.join(FEATS) + '\n')
    for i, r in enumerate(va):
        f.write(f"{r['uid']}\t{r['center']}\t{r['cls']}\t{r['y']}\t"
                + '\t'.join(f'{scores[n][i]:.6f}' for n in FEATS) + '\n')
json.dump({'results': results, 'n_train': len(tr), 'n_val': len(va), 'C': C_FIXED,
           'layer': 8, 'centers': CENTERS,
           'truncated_train': int(sum(r['trunc'] for r in tr)),
           'truncated_val': int(sum(r['trunc'] for r in va))},
          open(os.path.join(OUT, 'uesolds_esm_results.json'), 'w'), indent=1)
print(f"\nwrote {os.path.join(OUT,'uesolds_esm_results.tsv')}")
print(f"wrote {os.path.join(OUT,'uesolds_esm_val_predictions.tsv')}")
