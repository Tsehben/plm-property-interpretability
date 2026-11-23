#!/usr/bin/env python3
"""Layer sweep v2: checkpointed, disk-backed, resource-guarded. Protocol UNCHANGED.

Scientific settings are exactly layer_sweep/v1/frozen_protocol.json: cohort (eligible core,
10,104), layers 1-33, feature definitions, dtypes, classifier settings, training-only
preprocessing, and the P3 selection rule with lower-index tie-break frozen before any P4
prediction. v2 changes only HOW the pooled vectors are stored and when memory is released.

Resume/integrity behaviour is verified by layer_sweep/test_checkpoint_resume.py.

Usage: run_sweep_v2.py <repo_root> <v_dir> [--batch-tokens N] [--flush-every K]
"""
import os, sys, csv, json, time, hashlib, platform, warnings, collections, gc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, sklearn, transformers, torch
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.exceptions import ConvergenceWarning
import sweep_lib as S

REPO, V = sys.argv[1], sys.argv[2]
def arg(f, d): return type(d)(sys.argv[sys.argv.index(f) + 1]) if f in sys.argv else d
BT, FLUSH = arg('--batch-tokens', 6144), arg('--flush-every', 10)
C_FIXED, SEED, MAXIT, B, TOXIN = 1.0, 42, 5000, 10000, 'Q2VBP4'
STD_AA = 'ACDEFGHIKLMNPQRSTVWY'
STORE = os.path.join(V, 'emb_store')
log = lambda *a: print(*a, flush=True)
ENV = {'python': sys.version.split()[0], 'platform': platform.platform(),
       'machine': platform.machine(), 'numpy': np.__version__, 'torch': torch.__version__,
       'transformers': transformers.__version__, 'sklearn': sklearn.__version__}

def stop_report(stage, reason, extra):
    p = os.path.join(V, 'resource_stop_report.json')
    json.dump({'stage': stage, 'reason': reason, 'resources': S.resources(V),
               'nothing_deleted': True, 'store_retained': True,
               'resume_with': f'python3 layer_sweep/run_sweep_v2.py . {V} '
                              f'--batch-tokens {BT}', 'environment': ENV, **extra},
              open(p, 'w'), indent=1)
    log(f'\n*** STOPPED ({stage}): {reason} ***\nreport: {p}\n'
        'store and checkpoint retained; nothing deleted.')
    sys.exit(4)

# ---------------------------------------------------------------- cohort (frozen)
man, seqs = S.load_cohort(REPO)
assert len(man) == 10104, len(man)
clu = {r['ACC']: r['cluster_id30'] for r in csv.DictReader(
    open(os.path.join(REPO, 'homology_audit/v1/cluster_membership.tsv')), delimiter='\t')}
acc = [r['ACC'] for r in man]
y = np.array([int(r['target']) for r in man])
part = np.array([r['partition'] for r in man])
fit_m, p3_m, p4_m = np.isin(part, ['1', '2']), part == '3', part == '4'
batches = S.build_batches(acc, seqs, BT)
key = S.compat_key(acc, batches, BT)
log(f'core={len(acc)} fit={fit_m.sum()} P3={p3_m.sum()} P4={p4_m.sum()} batches={len(batches)}')
log(f'compat: model={key["model_checkpoint"]} code={key["extraction_code_sha256"][:12]} '
    f'batchmap={key["batch_map_sha256"][:12]} rows={key["row_order_sha256"][:12]}')

# ---------------------------------------------------------------- store + extract
r0 = S.resources(V); log(f'resources: disk={r0["disk_free_bytes"]/1e9:.2f}GB '
                         f'memfree={r0["mem_free_pct"]:.0f}% swapfree={r0["swap_free_bytes"]/1e6:.0f}MB')
os.makedirs(STORE, exist_ok=True)
ok, why = S.preallocate(STORE, len(acc), log=log)
if not ok: stop_report('preallocation', why, {})
completed, status = S.read_checkpoint(STORE, key)
if completed is None: stop_report('checkpoint', status, {})
log(f'checkpoint: {status} ({len(completed)}/{len(batches)} batches)')
mm = S.open_store(STORE, len(acc))
st, info = S.extract_into_store(acc, seqs, batches, STORE, mm, key, completed,
                                flush_every=FLUSH, log=log)
log(f'extraction: {st} {info}')
if st == 'stopped_resources': stop_report('extraction', info.get('reason', ''), {'info': info})
assert len(completed) == len(batches), f'store incomplete {len(completed)}/{len(batches)}'

prov = {'checkpoint': S.CHECKPOINT_ID, 'implementation': 'transformers.EsmModel', **ENV,
        'device': 'mps', 'layers_extracted': S.LAYERS,
        'layer_indexing_note': 'HF hidden_states[i] = output of transformer block i '
                               '= fair-esm representations[i]; hidden_states[0] not swept',
        'frozen': True, 'eval_mode': True, 'input': 'sequence only',
        'compute_dtype': S.COMPUTE_DTYPE, 'pooling_dtype': S.POOLING_DTYPE,
        'stored_dtype': S.STORE_DTYPE, 'pooling': S.POOLING_DESC,
        'pooling_device': 'on device; residue-level activations never retained',
        'storage': 'disk-backed memmap, one file per layer, cols [0:1280]=whole '
                   '[1280:2560]=first-30; incremental with fsync-before-commit checkpointing',
        'max_residues': S.MAX_RES, 'hidden_dim': S.HID, 'nterm_positions': S.NTERM,
        'batch_tokens': BT, 'n_sequences': len(acc), 'cohort': 'eligible core, in_core==1',
        'n_truncated': info.get('truncated'), 'extract_elapsed_s': info.get('elapsed_s'),
        'compat_key': key, 'per_layer_pooled_sha256': S.layer_hashes(mm)}
json.dump(prov, open(os.path.join(V, 'pooled_embedding_provenance.json'), 'w'), indent=1)
log('wrote pooled_embedding_provenance.json')

# model and device allocations were released inside extract_into_store()
gc.collect()
try: torch.mps.empty_cache()
except Exception: pass
log(f'pre-fitting resources: {S.resources(V)}')

# ---------------------------------------------------------------- probes, one layer at a time
def fracs(s):
    n = len(s)
    if n == 0: return [0.0] * 20
    c = collections.Counter(s)
    return [c.get(a, 0) / n for a in STD_AA]
COMP = np.array([fracs(seqs[a][:S.MAX_RES]) + [np.log(min(len(seqs[a]), S.MAX_RES))]
                 + fracs(seqs[a][:S.NTERM]) for a in acc], dtype=np.float64)

def fit_pred(X, te_m):
    sc = StandardScaler().fit(X[fit_m])
    a, b = sc.transform(X[fit_m]), sc.transform(X[te_m])
    assert np.isfinite(a).all() and np.isfinite(b).all()
    with warnings.catch_warnings(record=True) as cwl:
        warnings.simplefilter('always', ConvergenceWarning)
        clf = LogisticRegression(penalty='l2', C=C_FIXED, solver='lbfgs',
                                 max_iter=MAXIT, random_state=SEED).fit(a, y[fit_m])
        cv = not any(issubclass(x.category, ConvergenceWarning) for x in cwl)
    return (clf.predict_proba(b)[:, 1], int(clf.n_iter_[0]), bool(cv),
            sc.mean_.copy(), sc.scale_.copy(), clf.coef_[0].copy(), float(clf.intercept_[0]))

rows, p3_pred, params, t1 = [], {}, {}, time.time()
s_comp, ni_c, cv_c, m_c, sd_c, co_c, ic_c = fit_pred(COMP, p3_m)
params['P3_comp'] = (m_c, sd_c, co_c, ic_c)
log(f'comp (layer-independent) P3 AUROC={roc_auc_score(y[p3_m], s_comp):.4f} n_iter={ni_c}')
for layer in S.LAYERS:
    rr = S.resources(V); okr, whyr = S.resource_verdict(rr)
    if not okr: stop_report('fitting', whyr, {'layer_reached': layer, 'rows_so_far': len(rows)})
    E = np.asarray(mm[layer]).astype(np.float64)
    for nm, X in (('comp', COMP), ('esm', E), ('comp_plus_esm', np.hstack([COMP, E]))):
        if nm == 'comp':
            s, ni, cv = s_comp, ni_c, cv_c
        else:
            s, ni, cv, mu, sdv, co, ic = fit_pred(X, p3_m)
            params[f'P3_{nm}_L{layer}'] = (mu, sdv, co, ic)
        rows.append({'layer': layer, 'feature_set': nm, 'n_features': X.shape[1],
                     'n_fit': int(fit_m.sum()), 'n_P3': int(p3_m.sum()),
                     'n_P3_extracellular': int(y[p3_m].sum()),
                     'n_P3_cytoplasm': int((1 - y[p3_m]).sum()),
                     'auroc_P3': round(float(roc_auc_score(y[p3_m], s)), 6),
                     'ap_P3': round(float(average_precision_score(y[p3_m], s)), 6),
                     'n_iter': ni, 'converged': cv})
        p3_pred[f'{nm}_L{layer}'] = s.astype(np.float32)
        if nm == 'comp_plus_esm': del X
    del E; gc.collect()
    log(f"  L{layer:>2}: comp={rows[-3]['auroc_P3']:.4f} esm={rows[-2]['auroc_P3']:.4f} "
        f"comp+esm={rows[-1]['auroc_P3']:.4f}  [{time.time()-t1:,.0f}s]")
cols = ['layer','feature_set','n_features','n_fit','n_P3','n_P3_extracellular','n_P3_cytoplasm',
        'auroc_P3','ap_P3','n_iter','converged']
with open(os.path.join(V, 'p3_layer_curves.tsv'), 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in rows: f.write('\t'.join(str(r[c]) for c in cols) + '\n')
np.savez_compressed(os.path.join(V, 'p3_predictions.npz'),
                    acc=np.array(acc)[p3_m], y=y[p3_m], **p3_pred)

# ---------------- FREEZE before any P4 feature is built ----------------
cpe = [r for r in rows if r['feature_set'] == 'comp_plus_esm']
best = max(r['auroc_P3'] for r in cpe)
SELECTED = min(r['layer'] for r in cpe if r['auroc_P3'] == best)
nonconv = [(r['layer'], r['feature_set']) for r in rows if not r['converged']]
sel = {'criterion': 'highest comp_plus_esm AUROC on P3; lower layer index breaks exact ties',
       'selected_layer': SELECTED, 'selected_P3_auroc': best,
       'layer_8_P3_auroc': [r['auroc_P3'] for r in cpe if r['layer'] == 8][0],
       'frozen_before_any_P4_prediction': True, 'n_nonconverged_fits': len(nonconv),
       'nonconverged': nonconv, 'max_n_iter_observed': max(r['n_iter'] for r in rows),
       'p3_elapsed_s': round(time.time() - t1, 1)}
with open(os.path.join(V, 'layer_selection.json'), 'w') as f:
    json.dump(sel, f, indent=1); f.flush(); os.fsync(f.fileno())
log('\n*** SELECTION FROZEN *** ' + json.dumps(sel))

# ---------------------------------------------------------------- P4
EVAL = sorted({SELECTED, 8})
p4_acc = np.array(acc)[p4_m]; p4_y = y[p4_m]
p4_clu = np.array([clu[a] for a in p4_acc])
fit_clusters = {clu[a] for a in np.array(acc)[fit_m]}
overlap = np.array([c in fit_clusters for c in p4_clu])
log(f'P4 n={p4_m.sum()} clusters={len(set(p4_clu))} overlapping fit clusters={overlap.sum()}')
scores, meta = {}, {}
for L in EVAL:
    E = np.asarray(mm[L]).astype(np.float64)
    for nm, X in (('comp', COMP), ('esm', E), ('comp_plus_esm', np.hstack([COMP, E]))):
        s, ni, cv, mu, sdv, co, ic = fit_pred(X, p4_m)
        scores[(L, nm)] = s; meta[(L, nm)] = (X.shape[1], ni, cv)
        params[f'P4_{nm}_L{L}'] = (mu, sdv, co, ic)
    del E; gc.collect()
np.savez_compressed(os.path.join(V, 'probe_parameters.npz'),
                    artifact_class=np.array(['original_run_object']),
                    accession_order=np.array(acc), note=np.array([
                        'fitted during the sweep itself; NOT a reconstruction']),
                    **{f'scaler_mean__{k}': params[k][0] for k in params},
                    **{f'scaler_scale__{k}': params[k][1] for k in params},
                    **{f'coef__{k}': params[k][2] for k in params},
                    **{f'intercept__{k}': np.array([params[k][3]]) for k in params})
log(f'wrote probe_parameters.npz: {len(params)} probes (original-run objects)')

SUB = {'P4_all': np.ones(p4_m.sum(), bool),
       'P4_excl_fitting_set_cluster_overlap': ~overlap,
       'P4_excl_toxin_cluster_Q2VBP4': p4_clu != TOXIN}
def wts(m):
    c = collections.Counter(p4_clu[m]); return np.array([1.0 / c[x] for x in p4_clu[m]])
def met(m, s, w=None):
    if len(set(p4_y[m])) < 2: return None, None
    return (float(roc_auc_score(p4_y[m], s[m], sample_weight=w)),
            float(average_precision_score(p4_y[m], s[m], sample_weight=w)))
prs = []
for sn, sm in SUB.items():
    cc = collections.Counter(p4_clu[sm]); big, bign = cc.most_common(1)[0]
    for L in EVAL:
        for nm in ('comp', 'esm', 'comp_plus_esm'):
            au, ap = met(sm, scores[(L, nm)]); wau, wap = met(sm, scores[(L, nm)], wts(sm))
            nf, ni, cv = meta[(L, nm)]
            prs.append({'subset': sn, 'layer': L, 'feature_set': nm, 'n_features': nf,
                        'n': int(sm.sum()), 'n_extracellular': int(p4_y[sm].sum()),
                        'n_cytoplasm': int((1 - p4_y[sm]).sum()), 'n_clusters_id30': len(cc),
                        'largest_cluster': big, 'largest_cluster_n': bign,
                        'largest_cluster_pct': round(100 * bign / sm.sum(), 2),
                        'auroc_protein_weighted': round(au, 6), 'ap_protein_weighted': round(ap, 6),
                        'auroc_equal_cluster_weight': round(wau, 6),
                        'ap_equal_cluster_weight': round(wap, 6), 'n_iter': ni, 'converged': cv})
pc = list(prs[0])
with open(os.path.join(V, 'p4_results.tsv'), 'w') as f:
    f.write('\t'.join(pc) + '\n')
    for r in prs: f.write('\t'.join(str(r[c]) for c in pc) + '\n')

def pboot(m, sA, sB, b=B, seed=SEED):
    yy = p4_y[m]; a, bb = sA[m], sB[m]; g = p4_clu[m]
    gs = sorted(set(g)); gi = {x: np.where(g == x)[0] for x in gs}
    rng = np.random.default_rng(seed); dau, dap, und = [], [], 0
    for _ in range(b):
        d = rng.choice(len(gs), size=len(gs), replace=True)
        s = np.concatenate([gi[gs[x]] for x in d])
        if len(set(yy[s])) < 2: und += 1; continue
        dau.append(roc_auc_score(yy[s], a[s]) - roc_auc_score(yy[s], bb[s]))
        dap.append(average_precision_score(yy[s], a[s]) - average_precision_score(yy[s], bb[s]))
    return np.array(dau), np.array(dap), und
cmp_rows = []
for sn, sm in SUB.items():
    pairs = ([((SELECTED, 'comp_plus_esm'), (8, 'comp_plus_esm'),
               f'layer{SELECTED}_vs_layer8_comp_plus_esm')] if SELECTED != 8 else []) + \
            [((L, 'comp_plus_esm'), (L, 'comp'), f'layer{L}_comp_plus_esm_vs_comp') for L in EVAL] + \
            [((L, 'esm'), (L, 'comp'), f'layer{L}_esm_vs_comp') for L in EVAL]
    for A, Bk, lab in pairs:
        dau, dap, und = pboot(sm, scores[A], scores[Bk])
        oau = met(sm, scores[A])[0] - met(sm, scores[Bk])[0]
        oap = met(sm, scores[A])[1] - met(sm, scores[Bk])[1]
        cmp_rows.append({'subset': sn, 'comparison': lab, 'n': int(sm.sum()),
                         'n_clusters': len(set(p4_clu[sm])), 'delta_auroc': round(oau, 6),
                         'auroc_ci95_lo': round(float(np.percentile(dau, 2.5)), 6),
                         'auroc_ci95_hi': round(float(np.percentile(dau, 97.5)), 6),
                         'auroc_ci95_excludes_zero': int(np.percentile(dau, 2.5) > 0 or np.percentile(dau, 97.5) < 0),
                         'delta_ap': round(oap, 6),
                         'ap_ci95_lo': round(float(np.percentile(dap, 2.5)), 6),
                         'ap_ci95_hi': round(float(np.percentile(dap, 97.5)), 6),
                         'ap_ci95_excludes_zero': int(np.percentile(dap, 2.5) > 0 or np.percentile(dap, 97.5) < 0),
                         'n_resamples': B, 'n_resamples_undefined_excluded': und})
        log(f"  {sn:<38}{lab:<42} dAUROC={oau:+.4f} "
            f"[{cmp_rows[-1]['auroc_ci95_lo']:+.4f},{cmp_rows[-1]['auroc_ci95_hi']:+.4f}] undef={und}")
ccols = list(cmp_rows[0])
with open(os.path.join(V, 'p4_cluster_bootstrap.tsv'), 'w') as f:
    f.write('\t'.join(ccols) + '\n')
    for r in cmp_rows: f.write('\t'.join(str(r[c]) for c in ccols) + '\n')
np.savez_compressed(os.path.join(V, 'p4_predictions.npz'), acc=p4_acc, y=p4_y,
                    cluster=p4_clu, overlaps_fitting_set=overlap,
                    **{f'{nm}_L{L}': scores[(L, nm)].astype(np.float32) for (L, nm) in scores})
json.dump({'selected_layer': SELECTED, 'evaluated_layers': EVAL, 'results': prs,
           'comparisons': cmp_rows, 'n_P4': int(p4_m.sum()), 'environment': ENV,
           'n_P4_clusters': len(set(p4_clu)),
           'n_P4_overlapping_fitting_clusters': int(overlap.sum()),
           'P4_status': 'development evaluation; P4 is NOT a fresh final test set'},
          open(os.path.join(V, 'p4_results.json'), 'w'), indent=1)
log('\nALL DONE')
