#!/usr/bin/env python3
"""STAGE 1 of the layer sweep: fit on P1+P2, compare layers on P3. NO P4 PREDICTION.

Protocol frozen in layer_sweep/v1/frozen_protocol.json before this was run.
Usage: run_sweep_p3.py <repo_root> <v_dir>
"""
import os, sys, csv, json, time, hashlib, warnings, collections
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.exceptions import ConvergenceWarning

REPO, V = sys.argv[1], sys.argv[2]
STD_AA, MAX_RES, NTERM = 'ACDEFGHIKLMNPQRSTVWY', 1022, 30
C_FIXED, SEED, MAXIT = 1.0, 42, 5000
LAYERS = list(range(1, 34))

def fracs(s):
    n = len(s)
    if n == 0: return [0.0] * 20
    c = collections.Counter(s)
    return [c.get(a, 0) / n for a in STD_AA]

man = [r for r in csv.DictReader(open(os.path.join(REPO,
       'discovery_analysis/v1/inclusion_manifest.tsv')), delimiter='\t') if r['in_core'] == '1']
assert len(man) == 10104, len(man)
seqs = {}
with open(os.path.join(REPO, 'DeepLoc2/esm2_probe/deeploc_pairs.tsv')) as f:
    next(f)
    for line in f:
        k, s = line.rstrip('\n').split('\t'); seqs[k] = s
for r in man:
    assert hashlib.sha256(seqs[r['ACC']].encode()).hexdigest() == r['frozen_sha256']

ks, Ws, Ns = [], [], []
for fn in sorted(os.listdir(os.path.join(V, 'emb'))):
    if not (fn.startswith('sweep_shard') and fn.endswith('.npz')): continue
    z = np.load(os.path.join(V, 'emb', fn), allow_pickle=True)
    ks.append(z['keys']); Ws.append(z['whole']); Ns.append(z['nterm'])
    assert list(z['layers']) == LAYERS
K = np.concatenate(ks); W = np.concatenate(Ws); N = np.concatenate(Ns)
pos = {k: i for i, k in enumerate(K)}
print(f'embeddings: {W.shape} {W.dtype}  manifest core={len(man)}', flush=True)
assert all(r['ACC'] in pos for r in man), 'missing embedding for a core protein'

idx = np.array([pos[r['ACC']] for r in man])
y = np.array([int(r['target']) for r in man])
part = np.array([r['partition'] for r in man])
fit_m = np.isin(part, ['1', '2']); p3_m = part == '3'
print(f'fit P1+P2={fit_m.sum()}  compare P3={p3_m.sum()}  '
      f'(pos {y[fit_m].sum()}/{y[p3_m].sum()})', flush=True)

COMP = np.array([fracs(seqs[r['ACC']][:MAX_RES])
                 + [np.log(min(len(seqs[r['ACC']]), MAX_RES))]
                 + fracs(seqs[r['ACC']][:NTERM]) for r in man], dtype=np.float64)
print(f'comp features: {COMP.shape}', flush=True)

def fit_eval(Xtr, ytr, Xte, yte):
    sc = StandardScaler().fit(Xtr)
    a, b = sc.transform(Xtr), sc.transform(Xte)
    assert np.isfinite(a).all() and np.isfinite(b).all()
    with warnings.catch_warnings(record=True) as cw:
        warnings.simplefilter('always', ConvergenceWarning)
        clf = LogisticRegression(penalty='l2', C=C_FIXED, solver='lbfgs',
                                 max_iter=MAXIT, random_state=SEED).fit(a, ytr)
        conv = not any(issubclass(x.category, ConvergenceWarning) for x in cw)
    s = clf.predict_proba(b)[:, 1]
    return (float(roc_auc_score(yte, s)), float(average_precision_score(yte, s)),
            int(clf.n_iter_[0]), bool(conv), s)

rows, preds, t0 = [], {}, time.time()
# comp does not depend on the layer: fit once, reuse
au, ap, ni, cv, s_comp = fit_eval(COMP[fit_m], y[fit_m], COMP[p3_m], y[p3_m])
comp_row = dict(auroc=au, ap=ap, n_iter=ni, converged=cv)
print(f'comp (layer-independent): AUROC={au:.4f} AP={ap:.4f} n_iter={ni} conv={cv}', flush=True)

for li, layer in enumerate(LAYERS):
    E = np.hstack([W[idx, li, :].astype(np.float64), N[idx, li, :].astype(np.float64)])
    sets = {'comp': COMP, 'esm': E, 'comp_plus_esm': np.hstack([COMP, E])}
    for name, X in sets.items():
        if name == 'comp':
            r = dict(comp_row); s = s_comp
        else:
            a_, p_, ni_, cv_, s = fit_eval(X[fit_m], y[fit_m], X[p3_m], y[p3_m])
            r = dict(auroc=a_, ap=p_, n_iter=ni_, converged=cv_)
        rows.append({'layer': layer, 'feature_set': name, 'n_features': X.shape[1],
                     'n_fit': int(fit_m.sum()), 'n_P3': int(p3_m.sum()),
                     'n_P3_extracellular': int(y[p3_m].sum()),
                     'n_P3_cytoplasm': int((1 - y[p3_m]).sum()),
                     'auroc_P3': round(r['auroc'], 6), 'ap_P3': round(r['ap'], 6),
                     'n_iter': r['n_iter'], 'converged': r['converged']})
        preds[f'{name}_L{layer}'] = s
    del E, sets
    print(f"  layer {layer:>2}: comp={rows[-3]['auroc_P3']:.4f} "
          f"esm={rows[-2]['auroc_P3']:.4f} comp+esm={rows[-1]['auroc_P3']:.4f} "
          f"[{time.time()-t0:,.0f}s]", flush=True)

cols = ['layer','feature_set','n_features','n_fit','n_P3','n_P3_extracellular',
        'n_P3_cytoplasm','auroc_P3','ap_P3','n_iter','converged']
with open(os.path.join(V, 'p3_layer_curves.tsv'), 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in rows: f.write('\t'.join(str(r[c]) for c in cols) + '\n')

# prespecified selection: highest comp_plus_esm AUROC on P3, lower layer index breaks ties
cpe = [r for r in rows if r['feature_set'] == 'comp_plus_esm']
best = max(r['auroc_P3'] for r in cpe)
sel = min(r['layer'] for r in cpe if r['auroc_P3'] == best)
nonconv = [(r['layer'], r['feature_set']) for r in rows if not r['converged']]
sel_json = {'criterion': 'highest comp_plus_esm AUROC on P3; lower layer index breaks exact ties',
            'selected_layer': sel, 'selected_P3_auroc': best,
            'layer_8_P3_auroc': [r['auroc_P3'] for r in cpe if r['layer'] == 8][0],
            'frozen_before_any_P4_prediction': True,
            'n_nonconverged_fits': len(nonconv), 'nonconverged': nonconv,
            'max_n_iter_observed': max(r['n_iter'] for r in rows),
            'elapsed_s': round(time.time() - t0, 1)}
json.dump(sel_json, open(os.path.join(V, 'layer_selection.json'), 'w'), indent=1)
np.savez_compressed(os.path.join(V, 'p3_predictions.npz'),
                    acc=np.array([r['ACC'] for r in man])[p3_m],
                    y=y[p3_m], **{k: v for k, v in preds.items()})
print('\n' + json.dumps(sel_json, indent=1))
print(f"\nwrote p3_layer_curves.tsv, layer_selection.json, p3_predictions.npz")
