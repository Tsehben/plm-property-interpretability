#!/usr/bin/env python3
"""Layer sweep, single process: extract -> P3 compare -> FREEZE selection -> P4 evaluate.

WHY ONE PROCESS. The frozen protocol intended pooled vectors for all 33 layers to be
streamed to disk (1.71 GB as float16). A MEASURED constraint prevents that: only 1.34 GB of
disk is free and the host's swap (16.4 GB of a 17.4 GB swap file) cannot shrink, so the free
space is not recoverable. Rather than reduce the cohort or the layer range, pooled vectors
are held in RAM (1.71 GB of 25.8 GB) for the duration, and what is PERSISTED is:
  - full provenance for all 33 layers, including a SHA-256 of each layer's pooled array
  - the pooled vectors for the layers actually used downstream (selected layer and layer 8)
  - every prediction and metric table
Residue-level activations are never retained: pooling happens on device, per batch.

FREEZE ORDERING is enforced by program order and is auditable: layer_selection.json is
written and fsync'd BEFORE any P4 feature matrix is built. Nothing after that point can
change the selected layer.

Usage: run_full_sweep.py <repo_root> <v_dir> [--batch-tokens N]
"""
import os, sys, csv, json, time, hashlib, warnings, collections
import numpy as np, torch
from transformers import AutoTokenizer, EsmModel
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.exceptions import ConvergenceWarning

REPO, V = sys.argv[1], sys.argv[2]
BATCH_TOKENS = int(sys.argv[sys.argv.index('--batch-tokens') + 1]) if '--batch-tokens' in sys.argv else 6144
CHECKPOINT, MAX_RES, HID, NTERM = 'facebook/esm2_t33_650M_UR50D', 1022, 1280, 30
LAYERS = list(range(1, 34))
STD_AA, C_FIXED, SEED, MAXIT, B = 'ACDEFGHIKLMNPQRSTVWY', 1.0, 42, 5000, 10000
TOXIN = 'Q2VBP4'
os.makedirs(V, exist_ok=True)
log = lambda *a: print(*a, flush=True)

# ---------------------------------------------------------------- cohort
man = [r for r in csv.DictReader(open(os.path.join(REPO,
       'discovery_analysis/v1/inclusion_manifest.tsv')), delimiter='\t') if r['in_core'] == '1']
assert len(man) == 10104, len(man)
clu = {r['ACC']: r['cluster_id30'] for r in csv.DictReader(
    open(os.path.join(REPO, 'homology_audit/v1/cluster_membership.tsv')), delimiter='\t')}
seqs = {}
with open(os.path.join(REPO, 'DeepLoc2/esm2_probe/deeploc_pairs.tsv')) as f:
    next(f)
    for line in f:
        k, s = line.rstrip('\n').split('\t'); seqs[k] = s
for r in man:
    assert hashlib.sha256(seqs[r['ACC']].encode()).hexdigest() == r['frozen_sha256'], r['ACC']
acc = [r['ACC'] for r in man]
y = np.array([int(r['target']) for r in man])
part = np.array([r['partition'] for r in man])
fit_m, p3_m, p4_m = np.isin(part, ['1', '2']), part == '3', part == '4'
log(f'core={len(man)}  fit P1+P2={fit_m.sum()}  P3={p3_m.sum()}  P4={p4_m.sum()}  '
    f'(extracellular {y[fit_m].sum()}/{y[p3_m].sum()}/{y[p4_m].sum()})')

# ---------------------------------------------------------------- extraction
dev = torch.device('mps' if torch.backends.mps.is_available() else 'cpu')
tok = AutoTokenizer.from_pretrained(CHECKPOINT)
model = EsmModel.from_pretrained(CHECKPOINT, add_pooling_layer=False, torch_dtype=torch.float16)
model.eval().to(dev)
for p in model.parameters(): p.requires_grad_(False)

W = np.zeros((len(acc), len(LAYERS), HID), np.float16)
N = np.zeros((len(acc), len(LAYERS), HID), np.float16)
order = sorted(range(len(acc)), key=lambda i: -min(len(seqs[acc[i]]), MAX_RES))
batches, cur, cmax = [], [], 0
for i in order:
    L = min(len(seqs[acc[i]]), MAX_RES) + 2
    m = max(cmax, L)
    if cur and m * (len(cur) + 1) > BATCH_TOKENS: batches.append(cur); cur, cmax = [i], L
    else: cur.append(i); cmax = m
if cur: batches.append(cur)
log(f'device={dev} batch_tokens={BATCH_TOKENS} batches={len(batches)}')

t0, done, ntr, unk = time.time(), 0, 0, tok.unk_token_id
for bi, b in enumerate(batches):
    ss = [seqs[acc[i]][:MAX_RES] for i in b]
    ntr += sum(1 for i in b if len(seqs[acc[i]]) > MAX_RES)
    enc = tok(ss, return_tensors='pt', padding=True, add_special_tokens=True)
    assert (enc['input_ids'] == unk).sum().item() == 0
    for r, s in enumerate(ss):
        assert int(enc['attention_mask'][r].sum()) == len(s) + 2
    with torch.inference_mode():
        hs = model(input_ids=enc['input_ids'].to(dev),
                   attention_mask=enc['attention_mask'].to(dev),
                   output_hidden_states=True).hidden_states
    for li, layer in enumerate(LAYERS):
        h = hs[layer].float()
        for r, s in enumerate(ss):
            L = len(s); res = h[r, 1:1 + L, :]
            W[b[r], li] = res.mean(dim=0).cpu().numpy().astype(np.float16)
            N[b[r], li] = res[:min(NTERM, L)].mean(dim=0).cpu().numpy().astype(np.float16)
        del h
    del hs
    done += len(b)
    if bi % 50 == 0 or bi == len(batches) - 1:
        el = time.time() - t0
        log(f'  batch {bi+1}/{len(batches)}  {done:,}/{len(acc):,} seqs  {el:,.0f}s  {done/el:.2f} seq/s')
extract_s = time.time() - t0
del model
log(f'extraction done: {done:,} seqs, {extract_s:,.0f}s, truncated={ntr}')

prov = {'checkpoint': CHECKPOINT, 'implementation': 'transformers.EsmModel',
        'transformers_version': __import__('transformers').__version__,
        'torch_version': torch.__version__, 'device': str(dev),
        'layers_extracted': LAYERS,
        'layer_indexing_note': 'HF hidden_states[i] = output of transformer block i '
                               '= fair-esm representations[i]; hidden_states[0] not swept',
        'frozen': True, 'eval_mode': True, 'input': 'sequence only',
        'compute_dtype': 'float16', 'pooling_dtype': 'float32', 'stored_dtype': 'float16',
        'pooling': 'mean over residue tokens only; <cls>, <eos> and padding excluded',
        'pooling_device': 'on device; residue-level activations never retained',
        'max_residues': MAX_RES, 'hidden_dim': HID, 'nterm_positions': NTERM,
        'n_sequences': done, 'n_truncated': ntr, 'batch_tokens': BATCH_TOKENS,
        'elapsed_s': round(extract_s, 1), 'seq_per_s': round(done / extract_s, 3),
        'cohort': 'eligible core, in_core==1', 'n_cohort': len(acc),
        'per_layer_pooled_sha256': {
            str(l): {'whole': hashlib.sha256(W[:, i].tobytes()).hexdigest(),
                     'nterm': hashlib.sha256(N[:, i].tobytes()).hexdigest()}
            for i, l in enumerate(LAYERS)},
        'accession_order_sha256': hashlib.sha256('\n'.join(acc).encode()).hexdigest()}
json.dump(prov, open(os.path.join(V, 'pooled_embedding_provenance.json'), 'w'), indent=1)
log('wrote pooled_embedding_provenance.json')

# ---------------------------------------------------------------- probes
def fracs(s):
    n = len(s)
    if n == 0: return [0.0] * 20
    c = collections.Counter(s)
    return [c.get(a, 0) / n for a in STD_AA]
COMP = np.array([fracs(seqs[a][:MAX_RES]) + [np.log(min(len(seqs[a]), MAX_RES))]
                 + fracs(seqs[a][:NTERM]) for a in acc], dtype=np.float64)

def fit_pred(X, te_m):
    sc = StandardScaler().fit(X[fit_m])
    a, b = sc.transform(X[fit_m]), sc.transform(X[te_m])
    assert np.isfinite(a).all() and np.isfinite(b).all()
    with warnings.catch_warnings(record=True) as cwl:
        warnings.simplefilter('always', ConvergenceWarning)
        clf = LogisticRegression(penalty='l2', C=C_FIXED, solver='lbfgs',
                                 max_iter=MAXIT, random_state=SEED).fit(a, y[fit_m])
        conv = not any(issubclass(x.category, ConvergenceWarning) for x in cwl)
    return clf.predict_proba(b)[:, 1], int(clf.n_iter_[0]), bool(conv)

def emb(li): return np.hstack([W[:, li, :].astype(np.float64), N[:, li, :].astype(np.float64)])

rows, p3_pred, t1 = [], {}, time.time()
s_comp, ni_c, cv_c = fit_pred(COMP, p3_m)
log(f'comp (layer-independent) P3 AUROC={roc_auc_score(y[p3_m], s_comp):.4f} n_iter={ni_c} conv={cv_c}')
for li, layer in enumerate(LAYERS):
    E = emb(li)
    for nm, X in (('comp', COMP), ('esm', E), ('comp_plus_esm', np.hstack([COMP, E]))):
        if nm == 'comp': s, ni, cv = s_comp, ni_c, cv_c
        else: s, ni, cv = fit_pred(X, p3_m)
        rows.append({'layer': layer, 'feature_set': nm, 'n_features': X.shape[1],
                     'n_fit': int(fit_m.sum()), 'n_P3': int(p3_m.sum()),
                     'n_P3_extracellular': int(y[p3_m].sum()),
                     'n_P3_cytoplasm': int((1 - y[p3_m]).sum()),
                     'auroc_P3': round(float(roc_auc_score(y[p3_m], s)), 6),
                     'ap_P3': round(float(average_precision_score(y[p3_m], s)), 6),
                     'n_iter': ni, 'converged': cv})
        p3_pred[f'{nm}_L{layer}'] = s.astype(np.float32)
    del E
    log(f"  L{layer:>2}: comp={rows[-3]['auroc_P3']:.4f} esm={rows[-2]['auroc_P3']:.4f} "
        f"comp+esm={rows[-1]['auroc_P3']:.4f}  [{time.time()-t1:,.0f}s]")
cols = ['layer','feature_set','n_features','n_fit','n_P3','n_P3_extracellular','n_P3_cytoplasm',
        'auroc_P3','ap_P3','n_iter','converged']
with open(os.path.join(V, 'p3_layer_curves.tsv'), 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in rows: f.write('\t'.join(str(r[c]) for c in cols) + '\n')
np.savez_compressed(os.path.join(V, 'p3_predictions.npz'),
                    acc=np.array(acc)[p3_m], y=y[p3_m], **p3_pred)

# ---------------- FREEZE: selection written before ANY P4 feature is built ----------------
cpe = [r for r in rows if r['feature_set'] == 'comp_plus_esm']
best = max(r['auroc_P3'] for r in cpe)
SELECTED = min(r['layer'] for r in cpe if r['auroc_P3'] == best)
nonconv = [(r['layer'], r['feature_set']) for r in rows if not r['converged']]
sel = {'criterion': 'highest comp_plus_esm AUROC on P3; lower layer index breaks exact ties',
       'selected_layer': SELECTED, 'selected_P3_auroc': best,
       'layer_8_P3_auroc': [r['auroc_P3'] for r in cpe if r['layer'] == 8][0],
       'frozen_before_any_P4_prediction': True,
       'n_nonconverged_fits': len(nonconv), 'nonconverged': nonconv,
       'max_n_iter_observed': max(r['n_iter'] for r in rows),
       'p3_elapsed_s': round(time.time() - t1, 1)}
fp = os.path.join(V, 'layer_selection.json')
with open(fp, 'w') as f:
    json.dump(sel, f, indent=1); f.flush(); os.fsync(f.fileno())
log('\n*** SELECTION FROZEN *** ' + json.dumps(sel))

# ---------------------------------------------------------------- P4
EVAL = sorted({SELECTED, 8})
p4_acc = np.array(acc)[p4_m]; p4_y = y[p4_m]
p4_clu = np.array([clu[a] for a in p4_acc])
fit_clusters = {clu[a] for a in np.array(acc)[fit_m]}
overlap = np.array([c in fit_clusters for c in p4_clu])
log(f'P4 n={p4_m.sum()} clusters={len(set(p4_clu))} sharing cluster with P1+P2={overlap.sum()}')
scores, meta = {}, {}
for L in EVAL:
    li = LAYERS.index(L); E = emb(li)
    for nm, X in (('comp', COMP), ('esm', E), ('comp_plus_esm', np.hstack([COMP, E]))):
        s, ni, cv = fit_pred(X, p4_m); scores[(L, nm)] = s; meta[(L, nm)] = (X.shape[1], ni, cv)
    del E
# save pooled vectors ONLY for the evaluated layers (disk constraint)
np.savez_compressed(os.path.join(V, 'pooled_evaluated_layers.npz'),
                    acc=np.array(acc), layers=np.array(EVAL),
                    **{f'whole_L{L}': W[:, LAYERS.index(L)] for L in EVAL},
                    **{f'nterm_L{L}': N[:, LAYERS.index(L)] for L in EVAL})
del W, N

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
           'comparisons': cmp_rows, 'n_P4': int(p4_m.sum()),
           'n_P4_clusters': len(set(p4_clu)),
           'n_P4_overlapping_fitting_clusters': int(overlap.sum()),
           'P4_status': 'development evaluation; P4 is NOT a fresh final test set'},
          open(os.path.join(V, 'p4_results.json'), 'w'), indent=1)
log('\nALL DONE')
