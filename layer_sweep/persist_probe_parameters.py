#!/usr/bin/env python3
"""RECONSTRUCT and persist probe parameters for all 33 layers x 3 feature sets.

=============================== WHAT THESE ARTIFACTS ARE ===============================
Everything this script writes is a RECONSTRUCTED artifact, NOT an original-run object.
The estimators and scalers fitted during run_full_sweep.py existed only in that process's
memory and were never serialised; they are gone. This script re-extracts the pooled
embeddings and refits under the same frozen settings to produce equivalent objects. They
are labelled `reconstructed` in every output. A reconstructed coefficient vector is
evidence about the original fit only to the degree that the agreement checks below pass.
========================================================================================

ADDITIVE ONLY. Does not change the frozen selection or evaluation protocol, recomputes no
metric, and modifies no file written by run_full_sweep.py.

AGREEMENT GATE - criteria declared here, before running, and checked once:
  G1  Every one of the 33 layers must reproduce BOTH pooled arrays byte-for-byte
      (SHA-256 equal to pooled_embedding_provenance.json). Checked BEFORE any fitting.
  G2  Every reconstructed prediction must agree with the prediction run_full_sweep.py
      already saved, to PRED_TOL = 1e-10 absolute.
If G1 or G2 fails: the discrepancy is written to reconstruction_discrepancy.json, the
reconstruction STOPS, and NO reconstructed parameter file is written. The original sweep,
its selection and its evaluation are left exactly as they are. THE STEP IS NOT RETRIED:
this script makes a single attempt by design, because repeating a stochastic hardware
reduction until it agrees would manufacture the agreement it claims to test.

Usage: persist_probe_parameters.py <repo_root> <v_dir> [--batch-tokens N]
"""
import os, sys, csv, json, time, hashlib, platform, warnings, collections
import numpy as np, torch, sklearn, transformers
from transformers import AutoTokenizer, EsmModel
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.exceptions import ConvergenceWarning

REPO, V = sys.argv[1], sys.argv[2]
BT = int(sys.argv[sys.argv.index('--batch-tokens') + 1]) if '--batch-tokens' in sys.argv else 6144
CHECKPOINT, MAX_RES, HID, NTERM = 'facebook/esm2_t33_650M_UR50D', 1022, 1280, 30
LAYERS = list(range(1, 34))
STD_AA, C_FIXED, SEED, MAXIT = 'ACDEFGHIKLMNPQRSTVWY', 1.0, 42, 5000
PRED_TOL = 1e-10
log = lambda *a: print(*a, flush=True)

ENV = {'python': sys.version.split()[0], 'platform': platform.platform(),
       'machine': platform.machine(), 'numpy': np.__version__, 'torch': torch.__version__,
       'transformers': transformers.__version__, 'sklearn': sklearn.__version__}
FIT_SETTINGS = {'estimator': 'sklearn.linear_model.LogisticRegression',
                'penalty': 'l2', 'C': C_FIXED, 'solver': 'lbfgs', 'max_iter': MAXIT,
                'random_state': SEED, 'class_weight': None, 'fit_intercept': True,
                'scaler': 'sklearn.preprocessing.StandardScaler (with_mean=True, with_std=True)',
                'scaler_fitted_on': 'P1+P2 rows only',
                'feature_dtype': 'float64', 'pooled_storage_dtype': 'float16',
                'extraction_batch_tokens': BT, 'extraction_compute_dtype': 'float16',
                'pooling_dtype': 'float32'}

def halt(reason, payload):
    payload = {'status': 'STOPPED', 'gate_failed': reason, 'artifacts_written': False,
               'original_sweep_modified': False, 'retried': False,
               'retry_policy': 'single attempt by design; not rerun to obtain agreement',
               'environment': ENV, 'fit_settings': FIT_SETTINGS, **payload}
    json.dump(payload, open(os.path.join(V, 'reconstruction_discrepancy.json'), 'w'), indent=1)
    log(f'\n*** RECONSTRUCTION STOPPED: {reason} ***')
    log('discrepancy preserved in reconstruction_discrepancy.json; no parameters written; '
        'original sweep, selection and evaluation untouched; NOT retried.')
    sys.exit(2)

prov = json.load(open(os.path.join(V, 'pooled_embedding_provenance.json')))
sel = json.load(open(os.path.join(V, 'layer_selection.json')))
SELECTED = int(sel['selected_layer']); EVAL = sorted({SELECTED, 8})

man = [r for r in csv.DictReader(open(os.path.join(REPO,
       'discovery_analysis/v1/inclusion_manifest.tsv')), delimiter='\t') if r['in_core'] == '1']
seqs = {}
with open(os.path.join(REPO, 'DeepLoc2/esm2_probe/deeploc_pairs.tsv')) as f:
    next(f)
    for line in f:
        k, s = line.rstrip('\n').split('\t'); seqs[k] = s
acc = [r['ACC'] for r in man]
if hashlib.sha256('\n'.join(acc).encode()).hexdigest() != prov['accession_order_sha256']:
    halt('accession_order_sha256 mismatch',
         {'detail': 'the cohort or its ordering differs from the original sweep'})
y = np.array([int(r['target']) for r in man])
part = np.array([r['partition'] for r in man])
fit_m, p3_m, p4_m = np.isin(part, ['1', '2']), part == '3', part == '4'

# ---------------------------------------------------------------- re-extract (one attempt)
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
    if cur and m * (len(cur) + 1) > BT: batches.append(cur); cur, cmax = [i], L
    else: cur.append(i); cmax = m
if cur: batches.append(cur)
t0, done, unk = time.time(), 0, tok.unk_token_id
for bi, b in enumerate(batches):
    ss = [seqs[acc[i]][:MAX_RES] for i in b]
    enc = tok(ss, return_tensors='pt', padding=True, add_special_tokens=True)
    assert (enc['input_ids'] == unk).sum().item() == 0
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
    if bi % 100 == 0 or bi == len(batches) - 1:
        log(f'  batch {bi+1}/{len(batches)} {done:,} seqs {time.time()-t0:,.0f}s')
del model
extract_s = round(time.time() - t0, 1)

# ---------------------------------------------------------------- GATE G1 (before fitting)
hm, bad = {}, []
for i, l in enumerate(LAYERS):
    gw = hashlib.sha256(W[:, i].tobytes()).hexdigest()
    gn = hashlib.sha256(N[:, i].tobytes()).hexdigest()
    ew, en = prov['per_layer_pooled_sha256'][str(l)]['whole'], prov['per_layer_pooled_sha256'][str(l)]['nterm']
    hm[str(l)] = {'whole_match': gw == ew, 'nterm_match': gn == en}
    if gw != ew: bad.append({'layer': l, 'pooling': 'whole', 'original': ew, 'reconstructed': gw})
    if gn != en: bad.append({'layer': l, 'pooling': 'nterm', 'original': en, 'reconstructed': gn})
n_ok = sum(1 for v in hm.values() if v['whole_match'] and v['nterm_match'])
log(f'GATE G1: pooled arrays byte-identical for {n_ok}/33 layers')
if bad:
    halt('G1_pooled_embedding_hash_mismatch', {
        'n_layers_matching': n_ok, 'n_layers_total': 33,
        'n_mismatched_arrays': len(bad), 'mismatches': bad[:66],
        'per_layer': hm, 'reextraction_elapsed_s': extract_s,
        'interpretation': 'the re-extraction did not reproduce the original pooled values '
                          'byte-for-byte, so any refitted estimator would be fitted on '
                          'different inputs and could not be presented as a reconstruction '
                          'of the original fit. No parameters were written.'})

# ---------------------------------------------------------------- refit (inputs verified)
def fracs(s):
    n = len(s)
    if n == 0: return [0.0] * 20
    c = collections.Counter(s)
    return [c.get(a, 0) / n for a in STD_AA]
COMP = np.array([fracs(seqs[a][:MAX_RES]) + [np.log(min(len(seqs[a]), MAX_RES))]
                 + fracs(seqs[a][:NTERM]) for a in acc], dtype=np.float64)
p3saved = np.load(os.path.join(V, 'p3_predictions.npz'))
p4saved = np.load(os.path.join(V, 'p4_predictions.npz'))

store, diag, worst = {}, [], 0.0
for li, layer in enumerate(LAYERS):
    E = np.hstack([W[:, li, :].astype(np.float64), N[:, li, :].astype(np.float64)])
    for nm, X in (('comp', COMP), ('esm', E), ('comp_plus_esm', np.hstack([COMP, E]))):
        sc = StandardScaler().fit(X[fit_m])
        with warnings.catch_warnings(record=True) as cwl:
            warnings.simplefilter('always', ConvergenceWarning)
            clf = LogisticRegression(penalty='l2', C=C_FIXED, solver='lbfgs',
                                     max_iter=MAXIT, random_state=SEED).fit(sc.transform(X[fit_m]), y[fit_m])
            cv = not any(issubclass(x.category, ConvergenceWarning) for x in cwl)
        key = f'{nm}_L{layer}'
        d3 = float(np.abs(clf.predict_proba(sc.transform(X[p3_m]))[:, 1] - p3saved[key]).max())
        d4 = float(np.abs(clf.predict_proba(sc.transform(X[p4_m]))[:, 1] - p4saved[key]).max()) \
             if layer in EVAL else None
        worst = max(worst, d3, d4 or 0.0)
        if d3 > PRED_TOL or (d4 is not None and d4 > PRED_TOL):
            halt('G2_prediction_mismatch', {
                'first_failing_fit': {'layer': layer, 'feature_set': nm},
                'pred_tolerance': PRED_TOL,
                'max_abs_pred_diff_P3': d3, 'max_abs_pred_diff_P4': d4,
                'n_fits_checked_before_failure': len(diag),
                'diagnostics_up_to_failure': diag,
                'G1_passed': True, 'reextraction_elapsed_s': extract_s,
                'interpretation': 'pooled inputs were byte-identical, so a prediction '
                                  'disagreement above tolerance indicates the refit did not '
                                  'reproduce the original estimator. No parameters were '
                                  'written and nothing was retried.'})
        store[key] = (sc.mean_, sc.scale_, clf.coef_[0], float(clf.intercept_[0]))
        diag.append({'layer': layer, 'feature_set': nm, 'n_features': int(X.shape[1]),
                     'n_iter': int(clf.n_iter_[0]), 'converged': bool(cv),
                     'coef_l2_norm': float(np.linalg.norm(clf.coef_[0])),
                     'intercept': float(clf.intercept_[0]),
                     'pooled_hash_reproduced': True,
                     'max_abs_pred_diff_vs_saved_P3': d3,
                     'max_abs_pred_diff_vs_saved_P4': d4})
    del E
    log(f'  L{layer:>2} reconstructed: max pred diff {max(x["max_abs_pred_diff_vs_saved_P3"] for x in diag[-3:]):.3e}')
del W, N

np.savez_compressed(os.path.join(V, 'probe_parameters_all_layers_RECONSTRUCTED.npz'),
                    artifact_class=np.array(['reconstructed_not_original_run_object']),
                    accession_order=np.array(acc), layers=np.array(LAYERS),
                    **{f'coef__{k}': store[k][2] for k in store},
                    **{f'intercept__{k}': np.array([store[k][3]]) for k in store},
                    **{f'scaler_mean__{k}': store[k][0] for k in store},
                    **{f'scaler_scale__{k}': store[k][1] for k in store})
cols = ['layer','feature_set','n_features','n_iter','converged','coef_l2_norm','intercept',
        'pooled_hash_reproduced','max_abs_pred_diff_vs_saved_P3','max_abs_pred_diff_vs_saved_P4']
with open(os.path.join(V, 'probe_convergence_diagnostics_RECONSTRUCTED.tsv'), 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for d in diag: f.write('\t'.join('' if d.get(c) is None else str(d.get(c)) for c in cols) + '\n')
json.dump({'status': 'COMPLETED',
           'artifact_class': 'RECONSTRUCTED',
           'artifact_class_note': 'These estimators and scalers are RECONSTRUCTED artifacts, '
             'not the original-run objects. The original fitted objects lived only in the '
             'run_full_sweep.py process memory and were never serialised. These were refitted '
             'under the same frozen settings on pooled inputs verified byte-identical to the '
             'original, and every reconstructed prediction was checked against the predictions '
             'the original run saved.',
           'retried': False,
           'retry_policy': 'single attempt by design; not rerun to obtain agreement',
           'environment': ENV, 'fit_settings': FIT_SETTINGS,
           'gates': {'G1_pooled_hash_all_33_layers': True,
                     'G2_prediction_tolerance': PRED_TOL,
                     'G2_worst_abs_pred_diff_observed': worst},
           'n_fits': len(diag), 'selected_layer': SELECTED, 'evaluated_layers': EVAL,
           'n_nonconverged': sum(1 for d in diag if not d['converged']),
           'nonconverged': [(d['layer'], d['feature_set']) for d in diag if not d['converged']],
           'reextraction_elapsed_s': extract_s, 'diagnostics': diag},
          open(os.path.join(V, 'probe_parameters_RECONSTRUCTED_report.json'), 'w'), indent=1)
log(f'\nRECONSTRUCTION COMPLETE: {len(diag)} fits, worst prediction diff {worst:.3e} '
    f'(tolerance {PRED_TOL:.0e})')
log('wrote probe_parameters_all_layers_RECONSTRUCTED.npz, '
    'probe_convergence_diagnostics_RECONSTRUCTED.tsv, probe_parameters_RECONSTRUCTED_report.json')
