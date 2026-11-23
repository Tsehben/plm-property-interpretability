#!/usr/bin/env python3
"""Bounded localisation-probe ablation pilot. Protocol frozen in frozen_protocol.json.

Readout: change in the FROZEN layer-26 ESM-only probe's extracellular log-odds.
Edit: validated direct decoder-delta at hidden_states[8] (= output of encoder.layer[7]).
P4 and P0 are never read.

Usage: run_probe_ablation.py <repo_root> <out_dir>
"""
import csv, json, os, sys, hashlib, collections
import numpy as np, torch
from transformers import AutoTokenizer, EsmModel
from sklearn.metrics import roc_auc_score

REPO, OUT = sys.argv[1], sys.argv[2]
sys.path.insert(0, os.path.join(REPO, 'interprot_ref/raw'))
from sae_model import SparseAutoencoder
from safetensors.torch import load_file

CKPT, HIDDEN, SAE_DIM, K = 'facebook/esm2_t33_650M_UR50D', 1280, 4096, 128
HOOK_LAYER, PROBE_LAYER, NTERM, MAX_RES = 7, 26, 30, 1022
FEATURE, POSITIONS = 1446, [10, 20, 30]
log = lambda *a: print(*a, flush=True)

# ---------------- inputs ----------------
sample = list(csv.DictReader(open(os.path.join(
    REPO, 'intervention_pilot/v1/pilot_sample.tsv')), delimiter='\t'))
assert len(sample) == 24, len(sample)
seqs = {}
with open(os.path.join(REPO, 'DeepLoc2/esm2_probe/deeploc_pairs.tsv')) as f:
    next(f)
    for line in f:
        k, s = line.rstrip('\n').split('\t'); seqs[k] = s
man = {r['ACC']: r for r in csv.DictReader(open(os.path.join(
    REPO, 'discovery_analysis/v1/inclusion_manifest.tsv')), delimiter='\t') if r['in_core'] == '1'}
for r in sample:
    assert man[r['ACC']]['partition'] == '3', 'pilot must be P3 only'
    assert hashlib.sha256(seqs[r['ACC']].encode()).hexdigest() == r['sha256']

P = np.load(os.path.join(REPO, 'layer_sweep/v1/probe_parameters.npz'), allow_pickle=True)
coef = P['coef__P3_esm_L26'].astype(np.float64)
intercept = float(P['intercept__P3_esm_L26'][0])
mu = P['scaler_mean__P3_esm_L26'].astype(np.float64)
sd = P['scaler_scale__P3_esm_L26'].astype(np.float64)
assert coef.shape == (2 * HIDDEN,) and mu.shape == sd.shape == (2 * HIDDEN,)
# the probe fitted for P3 and for P4 are the same fit (both on P1+P2); confirm
assert np.array_equal(coef, P['coef__P4_esm_L26'].astype(np.float64)), 'P3/P4 probe differ'

zp3 = np.load(os.path.join(REPO, 'layer_sweep/v1/p3_predictions.npz'), allow_pickle=True)
saved = {str(a): float(p) for a, p in zip(zp3['acc'], zp3['esm_L26'])}

dev = torch.device('mps' if torch.backends.mps.is_available() else 'cpu')
tok = AutoTokenizer.from_pretrained(CKPT)
model = EsmModel.from_pretrained(CKPT, add_pooling_layer=False, torch_dtype=torch.float32)
model.eval().to(dev)
for p_ in model.parameters(): p_.requires_grad_(False)
sae = SparseAutoencoder(HIDDEN, SAE_DIM)
sae.load_state_dict(load_file(os.path.join(REPO, 'interprot_ref/raw/esm2_plm1280_l8_sae4096.safetensors')))
sae.eval()
log(f'device={dev}  probe=layer{PROBE_LAYER} ESM-only  feature={FEATURE}  positions={POSITIONS}')

# ---------------- edit hook (validated direct decoder-delta) ----------------
class Edit:
    def __init__(self, mode, pos=None, vec=None):
        self.mode, self.pos, self.vec = mode, pos, vec
        self.info = {}
    def __call__(self, module, inp, out):
        is_tuple = isinstance(out, tuple)
        h = out[0] if is_tuple else out
        x = h[0].float().cpu()
        if self.mode == 'raw':
            a = sae.get_acts(x)
            self.info = {'latent': float(a[self.pos, FEATURE])}
            return out
        if self.mode == 'random':
            delta = torch.zeros_like(x)
            delta[self.pos] = self.vec
            self.info = {'delta_h_l2': float(delta.norm()),
                         'n_positions_changed': int((delta.abs().sum(-1) > 0).sum())}
        else:                                   # 'noop' or 'ablate'
            a = sae.get_acts(x)
            before = float(a[self.pos, FEATURE])
            a_ed = a.clone()
            a_ed[self.pos, FEATURE] = before if self.mode == 'noop' else 0.0
            delta = ((sae.topK_activation(a_ed, sae.k) - sae.topK_activation(a, sae.k))
                     @ sae.w_dec) * sae.encode(x)[2]
            self.info = {'latent_before': before,
                         'latent_after': float(a_ed[self.pos, FEATURE]),
                         'delta_h_l2': float(delta.norm()),
                         'n_latent_coords_changed': int((a_ed != a).sum()),
                         'n_positions_changed': int((delta.abs().sum(-1) > 0).sum())}
        new_h = h.clone()
        new_h[0] = (x + delta).to(h.dtype).to(h.device)
        return (new_h,) + out[1:] if is_tuple else new_h

def forward_logodds(ids, am, L, hook=None):
    handle = model.encoder.layer[HOOK_LAYER].register_forward_hook(hook) if hook else None
    try:
        with torch.inference_mode():
            hs = model(input_ids=ids, attention_mask=am, output_hidden_states=True).hidden_states
        h = hs[PROBE_LAYER][0].float().cpu().numpy()
    finally:
        if handle: handle.remove()
    res = h[1:1 + L]
    feat = np.concatenate([res.mean(axis=0), res[:min(NTERM, L)].mean(axis=0)]).astype(np.float64)
    z = (feat - mu) / sd
    return float(z @ coef + intercept), feat

def seed_for(acc, pos):
    return int.from_bytes(hashlib.sha256(f'{acc}:{pos}'.encode()).digest()[:8], 'big')

# ---- verify the hook site IS hidden_states[8] for THIS module class, not assumed ----
_probe_seq = seqs[sample[0]['ACC']][:64]
_e = tok([_probe_seq], return_tensors='pt', padding=True, add_special_tokens=True)
_cap = {}
def _grab(module, inp, out):
    _cap['h'] = (out[0] if isinstance(out, tuple) else out).detach().float().cpu().clone()
_h = model.encoder.layer[HOOK_LAYER].register_forward_hook(_grab)
with torch.inference_mode():
    _hs = model(input_ids=_e['input_ids'].to(dev), attention_mask=_e['attention_mask'].to(dev),
                output_hidden_states=True).hidden_states
_h.remove()
_ref = _hs[8].detach().float().cpu()
_hook_ok = bool(torch.equal(_cap['h'], _ref))
_maxdiff = float((_cap['h'] - _ref).abs().max())
log(f'hook-site check: encoder.layer[{HOOK_LAYER}] output == hidden_states[8] -> '
    f'{_hook_ok} (max abs diff {_maxdiff:.3e})')
assert _hook_ok, 'hook site is not hidden_states[8] for EsmModel'
_ndown = len(model.encoder.layer) - (HOOK_LAYER + 1)
log(f'  downstream blocks after the hook: {_ndown}; probe reads hidden_states[{PROBE_LAYER}]')

rows, per_protein, checks = [], [], []
for i, r in enumerate(sample, 1):
    acc = r['ACC']; s = seqs[acc][:MAX_RES]; L = len(s)
    enc = tok([s], return_tensors='pt', padding=True, add_special_tokens=True)
    ids, am = enc['input_ids'].to(dev), enc['attention_mask'].to(dev)
    assert int(am.sum()) == L + 2
    lo_clean, _ = forward_logodds(ids, am, L)
    p_clean = 1 / (1 + np.exp(-lo_clean))
    p_saved = saved[acc]
    lo_saved = np.log(p_saved / (1 - p_saved)) if 0 < p_saved < 1 else np.nan
    checks.append({'ACC': acc, 'p_clean': p_clean, 'p_saved': p_saved,
                   'abs_dp': abs(p_clean - p_saved),
                   'lo_clean': lo_clean, 'lo_saved': lo_saved,
                   'abs_dlo': abs(lo_clean - lo_saved) if np.isfinite(lo_saved) else np.nan})
    d_abl, d_rnd, acts = [], [], []
    for pos in POSITIONS:
        tpos = pos                       # token index: <cls> at 0, residue j at index j
        assert 1 <= pos <= L, f'{acc}: position {pos} outside length {L}'
        hr = Edit('raw', tpos); forward_logodds(ids, am, L, hr)
        latent = hr.info['latent']; active = latent > 0
        hn = Edit('noop', tpos); lo_noop, _ = forward_logodds(ids, am, L, hn)
        ha = Edit('ablate', tpos); lo_abl, _ = forward_logodds(ids, am, L, ha)
        nrm = ha.info['delta_h_l2']
        g = torch.Generator().manual_seed(seed_for(acc, pos))
        v = torch.randn(HIDDEN, generator=g, dtype=torch.float32)
        v = v / v.norm() * nrm if nrm > 0 else torch.zeros(HIDDEN)
        hx = Edit('random', tpos, v); lo_rnd, _ = forward_logodds(ids, am, L, hx)
        rows.append({'ACC': acc, 'target': r['target'],
                     'localisation_class': r['localisation_class'], 'position': pos,
                     'latent_1446': round(latent, 6), 'active_1446': int(active),
                     'logodds_clean': round(lo_clean, 6),
                     'logodds_noop': round(lo_noop, 6),
                     'logodds_ablate': round(lo_abl, 6),
                     'logodds_random': round(lo_rnd, 6),
                     'd_noop': round(lo_noop - lo_clean, 9),
                     'd_ablate': round(lo_abl - lo_clean, 6),
                     'd_random': round(lo_rnd - lo_clean, 6),
                     'delta_h_l2_ablate': round(nrm, 6),
                     'delta_h_l2_random': round(hx.info['delta_h_l2'], 6),
                     'noop_delta_h_l2': round(hn.info['delta_h_l2'], 12),
                     'noop_latent_coords_changed': hn.info['n_latent_coords_changed'],
                     'ablate_positions_changed': ha.info['n_positions_changed']})
        d_abl.append(lo_abl - lo_clean); d_rnd.append(lo_rnd - lo_clean); acts.append(int(active))
    per_protein.append({'ACC': acc, 'target': int(r['target']),
                        'localisation_class': r['localisation_class'],
                        'logodds_clean': round(lo_clean, 6),
                        'mean_d_ablate': round(float(np.mean(d_abl)), 6),
                        'mean_d_random': round(float(np.mean(d_rnd)), 6),
                        'paired_candidate_minus_control': round(float(np.mean(d_abl) - np.mean(d_rnd)), 6),
                        'n_positions_active_1446': int(sum(acts)),
                        'mean_logodds_ablate': round(lo_clean + float(np.mean(d_abl)), 6)})
    log(f'  [{i:>2}/24] {acc}  clean_lo={lo_clean:+.4f}  active={sum(acts)}/3  '
        f'mean_d_abl={np.mean(d_abl):+.5f}  mean_d_rnd={np.mean(d_rnd):+.5f}')

# ---------------- verification ----------------
dp = np.array([c['abs_dp'] for c in checks]); dl = np.array([c['abs_dlo'] for c in checks])
noop_max = max(abs(r['d_noop']) for r in rows)
noop_norm_max = max(r['noop_delta_h_l2'] for r in rows)
ver = {'clean_vs_saved': {'n': len(checks),
        'max_abs_prob_diff': float(dp.max()), 'median_abs_prob_diff': float(np.median(dp)),
        'max_abs_logodds_diff': float(np.nanmax(dl)), 'median_abs_logodds_diff': float(np.nanmedian(dl)),
        'tolerance_prob_max': 0.02, 'tolerance_logodds_median': 0.10, 'tolerance_logodds_max': 0.50,
        'passes_prob': bool(dp.max() <= 0.02),
        'passes_logodds_median': bool(np.nanmedian(dl) <= 0.10),
        'passes_logodds_max': bool(np.nanmax(dl) <= 0.50)},
       'no_op_vs_clean': {'max_abs_logodds_diff': noop_max,
        'max_delta_h_l2': noop_norm_max, 'tolerance': 1e-6,
        'passes': bool(noop_max < 1e-6 and noop_norm_max == 0.0)}}

y = np.array([p['target'] for p in per_protein])
auc_clean = float(roc_auc_score(y, [p['logodds_clean'] for p in per_protein]))
auc_abl = float(roc_auc_score(y, [p['mean_logodds_ablate'] for p in per_protein]))
d = np.array([p['paired_candidate_minus_control'] for p in per_protein])
da = np.array([p['mean_d_ablate'] for p in per_protein])
dr = np.array([p['mean_d_random'] for p in per_protein])
act = collections.Counter()
for r in rows: act[r['position']] += r['active_1446']

summary = {
 'label': 'exploratory P3 tuning pilot; NOT an evaluation and NOT a held-out result',
 'n_proteins': 24, 'positions': POSITIONS, 'feature': FEATURE,
 'P4_accessed': False, 'P0_accessed': False, 'masking_used': False,
 'hook_site_verified_equals_hidden_states_8': True,
 'n_downstream_blocks_after_hook': len(model.encoder.layer) - (HOOK_LAYER + 1),
 'verification': ver,
 'activity_of_1446': {'by_position': {str(k): v for k, v in sorted(act.items())},
   'total_active_cases': int(sum(act.values())), 'total_cases': len(rows),
   'proteins_with_any_active': int(sum(1 for p in per_protein if p['n_positions_active_1446'] > 0)),
   'proteins_with_all_three_active': int(sum(1 for p in per_protein if p['n_positions_active_1446'] == 3))},
 'primary': {'readout': 'change in extracellular log-odds, edited minus clean, averaged over the 3 positions within each protein',
   'mean_d_ablate': float(da.mean()), 'median_d_ablate': float(np.median(da)),
   'mean_d_random': float(dr.mean()), 'median_d_random': float(np.median(dr)),
   'paired_candidate_minus_control_mean': float(d.mean()),
   'paired_candidate_minus_control_median': float(np.median(d)),
   'paired_sd': float(d.std(ddof=1)), 'paired_min': float(d.min()), 'paired_max': float(d.max()),
   'n_positive': int((d > 0).sum()), 'n_negative': int((d < 0).sum()), 'n_zero': int((d == 0).sum())},
 'descriptive_auroc': {'clean': round(auc_clean, 6), 'ablated': round(auc_abl, 6),
   'delta': round(auc_abl - auc_clean, 6), 'n': 24,
   'note': 'descriptive only; 24 proteins, a manipulation check on our own probe'},
}
json.dump({'summary': summary, 'per_protein': per_protein, 'clean_checks': checks},
          open(os.path.join(OUT, 'pilot_results.json'), 'w'), indent=1)
cols = list(rows[0])
with open(os.path.join(OUT, 'pilot_conditions.tsv'), 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in rows: f.write('\t'.join(str(r[c]) for c in cols) + '\n')
pcols = list(per_protein[0])
with open(os.path.join(OUT, 'pilot_per_protein.tsv'), 'w') as f:
    f.write('\t'.join(pcols) + '\n')
    for r in per_protein: f.write('\t'.join(str(r[c]) for c in pcols) + '\n')
log('\n' + json.dumps(summary, indent=1))
