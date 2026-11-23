#!/usr/bin/env python3
"""STAGE 2: evaluate the SELECTED layer and layer 8 on P4.

Runs only after layer_selection.json exists; the selected layer is read from it and is
never re-chosen here. Usage: run_sweep_p4.py <repo_root> <v_dir>
"""
import os, sys, csv, json, hashlib, warnings, collections
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.exceptions import ConvergenceWarning

REPO, V = sys.argv[1], sys.argv[2]
STD_AA, MAX_RES, NTERM = 'ACDEFGHIKLMNPQRSTVWY', 1022, 30
C_FIXED, SEED, MAXIT, B = 1.0, 42, 5000, 10000
LAYERS = list(range(1, 34))
TOXIN = 'Q2VBP4'

sel_j = json.load(open(os.path.join(V, 'layer_selection.json')))
SELECTED = int(sel_j['selected_layer'])
EVAL_LAYERS = sorted({SELECTED, 8})
print(f'selected layer (frozen on P3) = {SELECTED}; evaluating layers {EVAL_LAYERS} on P4')

def fracs(s):
    n = len(s)
    if n == 0: return [0.0] * 20
    c = collections.Counter(s)
    return [c.get(a, 0) / n for a in STD_AA]

man = [r for r in csv.DictReader(open(os.path.join(REPO,
       'discovery_analysis/v1/inclusion_manifest.tsv')), delimiter='\t') if r['in_core'] == '1']
clu = {r['ACC']: r['cluster_id30'] for r in csv.DictReader(
    open(os.path.join(REPO, 'homology_audit/v1/cluster_membership.tsv')), delimiter='\t')}
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
K = np.concatenate(ks); W = np.concatenate(Ws); N = np.concatenate(Ns)
pos = {k: i for i, k in enumerate(K)}
idx = np.array([pos[r['ACC']] for r in man])
acc = np.array([r['ACC'] for r in man])
y = np.array([int(r['target']) for r in man])
part = np.array([r['partition'] for r in man])
fit_m, p4_m = np.isin(part, ['1', '2']), part == '4'
COMP = np.array([fracs(seqs[a][:MAX_RES]) + [np.log(min(len(seqs[a]), MAX_RES))]
                 + fracs(seqs[a][:NTERM]) for a in acc], dtype=np.float64)

fit_clusters = {clu[a] for a in acc[fit_m]}
p4_acc, p4_y = acc[p4_m], y[p4_m]
p4_clu = np.array([clu[a] for a in p4_acc])
overlap = np.array([c in fit_clusters for c in p4_clu])
print(f'P4 n={p4_m.sum()}  extracellular={p4_y.sum()}  clusters={len(set(p4_clu))}  '
      f'sharing a cluster with P1+P2: {overlap.sum()}')

def fit_predict(X):
    sc = StandardScaler().fit(X[fit_m])
    with warnings.catch_warnings(record=True) as cw:
        warnings.simplefilter('always', ConvergenceWarning)
        clf = LogisticRegression(penalty='l2', C=C_FIXED, solver='lbfgs',
                                 max_iter=MAXIT, random_state=SEED).fit(sc.transform(X[fit_m]), y[fit_m])
        conv = not any(issubclass(x.category, ConvergenceWarning) for x in cw)
    return clf.predict_proba(sc.transform(X[p4_m]))[:, 1], int(clf.n_iter_[0]), bool(conv)

scores, meta = {}, {}
for L in EVAL_LAYERS:
    li = LAYERS.index(L)
    E = np.hstack([W[idx, li, :].astype(np.float64), N[idx, li, :].astype(np.float64)])
    for nm, X in (('comp', COMP), ('esm', E), ('comp_plus_esm', np.hstack([COMP, E]))):
        s, ni, cv = fit_predict(X)
        scores[(L, nm)] = s; meta[(L, nm)] = (X.shape[1], ni, cv)
    del E

def cw_weights(mask):
    c = collections.Counter(p4_clu[mask])
    return np.array([1.0 / c[x] for x in p4_clu[mask]])

def metrics(mask, s, w=None):
    yy, ss = p4_y[mask], s[mask]
    if len(set(yy)) < 2: return None, None
    return (float(roc_auc_score(yy, ss, sample_weight=w)),
            float(average_precision_score(yy, ss, sample_weight=w)))

SUBSETS = {
    'P4_all': np.ones(p4_m.sum(), bool),
    'P4_excl_fitting_set_cluster_overlap': ~overlap,
    'P4_excl_toxin_cluster_Q2VBP4': p4_clu != TOXIN,
}
rows = []
for sname, smask in SUBSETS.items():
    cc = collections.Counter(p4_clu[smask])
    big, bign = cc.most_common(1)[0]
    for L in EVAL_LAYERS:
        for nm in ('comp', 'esm', 'comp_plus_esm'):
            au, ap = metrics(smask, scores[(L, nm)])
            wau, wap = metrics(smask, scores[(L, nm)], cw_weights(smask))
            nf, ni, cv = meta[(L, nm)]
            rows.append({'subset': sname, 'layer': L, 'feature_set': nm, 'n_features': nf,
                         'n': int(smask.sum()),
                         'n_extracellular': int(p4_y[smask].sum()),
                         'n_cytoplasm': int((1 - p4_y[smask]).sum()),
                         'n_clusters_id30': len(cc), 'largest_cluster': big,
                         'largest_cluster_n': bign,
                         'largest_cluster_pct': round(100 * bign / smask.sum(), 2),
                         'auroc_protein_weighted': round(au, 6),
                         'ap_protein_weighted': round(ap, 6),
                         'auroc_equal_cluster_weight': round(wau, 6),
                         'ap_equal_cluster_weight': round(wap, 6),
                         'n_iter': ni, 'converged': cv})
cols = list(rows[0])
with open(os.path.join(V, 'p4_results.tsv'), 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in rows: f.write('\t'.join(str(r[c]) for c in cols) + '\n')

def paired_cluster_boot(mask, sA, sB, b=B, seed=SEED):
    """Same cluster resample for both models; pairing preserved. Returns diffs + undefined count."""
    yy = p4_y[mask]; a, bb = sA[mask], sB[mask]; g = p4_clu[mask]
    gs = sorted(set(g)); gidx = {x: np.where(g == x)[0] for x in gs}
    rng = np.random.default_rng(seed)
    dau, dap, undef = [], [], 0
    for _ in range(b):
        drawn = rng.choice(len(gs), size=len(gs), replace=True)
        sel = np.concatenate([gidx[gs[d]] for d in drawn])
        if len(set(yy[sel])) < 2: undef += 1; continue
        dau.append(roc_auc_score(yy[sel], a[sel]) - roc_auc_score(yy[sel], bb[sel]))
        dap.append(average_precision_score(yy[sel], a[sel]) - average_precision_score(yy[sel], bb[sel]))
    return np.array(dau), np.array(dap), undef

COMPARISONS = []
for sname, smask in SUBSETS.items():
    COMPARISONS.append((sname, f'layer{SELECTED}_vs_layer8_comp_plus_esm',
                        (SELECTED, 'comp_plus_esm'), (8, 'comp_plus_esm'), smask))
    for L in EVAL_LAYERS:
        COMPARISONS.append((sname, f'layer{L}_comp_plus_esm_vs_comp',
                            (L, 'comp_plus_esm'), (L, 'comp'), smask))
crows = []
for sname, label, A, Bk, smask in COMPARISONS:
    if A == Bk: continue
    dau, dap, undef = paired_cluster_boot(smask, scores[A], scores[Bk])
    oau = metrics(smask, scores[A])[0] - metrics(smask, scores[Bk])[0]
    oap = metrics(smask, scores[A])[1] - metrics(smask, scores[Bk])[1]
    crows.append({'subset': sname, 'comparison': label,
                  'n': int(smask.sum()), 'n_clusters': len(set(p4_clu[smask])),
                  'delta_auroc': round(oau, 6),
                  'auroc_ci95_lo': round(float(np.percentile(dau, 2.5)), 6),
                  'auroc_ci95_hi': round(float(np.percentile(dau, 97.5)), 6),
                  'auroc_ci95_excludes_zero': int(np.percentile(dau, 2.5) > 0 or np.percentile(dau, 97.5) < 0),
                  'delta_ap': round(oap, 6),
                  'ap_ci95_lo': round(float(np.percentile(dap, 2.5)), 6),
                  'ap_ci95_hi': round(float(np.percentile(dap, 97.5)), 6),
                  'ap_ci95_excludes_zero': int(np.percentile(dap, 2.5) > 0 or np.percentile(dap, 97.5) < 0),
                  'n_resamples': B, 'n_resamples_undefined_excluded': undef})
    print(f"  {sname:<38} {label:<38} dAUROC={oau:+.4f} "
          f"[{crows[-1]['auroc_ci95_lo']:+.4f},{crows[-1]['auroc_ci95_hi']:+.4f}] undef={undef}", flush=True)
ccols = list(crows[0])
with open(os.path.join(V, 'p4_cluster_bootstrap.tsv'), 'w') as f:
    f.write('\t'.join(ccols) + '\n')
    for r in crows: f.write('\t'.join(str(r[c]) for c in ccols) + '\n')

np.savez_compressed(os.path.join(V, 'p4_predictions.npz'), acc=p4_acc, y=p4_y,
                    cluster=p4_clu, overlaps_fitting_set=overlap,
                    **{f'{nm}_L{L}': scores[(L, nm)] for (L, nm) in scores})
json.dump({'selected_layer': SELECTED, 'evaluated_layers': EVAL_LAYERS,
           'results': rows, 'comparisons': crows,
           'P4_status': 'development evaluation; P4 is not a fresh final test set',
           'n_P4': int(p4_m.sum()), 'n_P4_clusters': len(set(p4_clu)),
           'n_P4_overlapping_fitting_clusters': int(overlap.sum())},
          open(os.path.join(V, 'p4_results.json'), 'w'), indent=1)
print('\nwrote p4_results.tsv, p4_cluster_bootstrap.tsv, p4_results.json, p4_predictions.npz')
