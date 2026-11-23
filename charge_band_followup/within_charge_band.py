#!/usr/bin/env python3
"""Within-charge-band comparison of the frozen layer-26 ESM-only probe and the matched
composition baseline, on all 2,366 P4 proteins. Saved predictions only; no refitting.

Primary comparison: only extracellular-cytoplasmic pairs from the SAME charge band are
used. Per band, AUROC is the concordance over that band's pairs. The five band AUROCs are
averaged weighted by the number of eligible pairs, w_b = n_extracellular_b * n_cytoplasmic_b.

Interval: whole 30%-identity clusters resampled across P4, the SAME draw applied to both
models and all five bands, with the pair weights RECOMPUTED inside every draw.

Usage: within_charge_band.py <repo_root> <out_dir>
"""
import csv, json, os, sys, collections
import numpy as np

REPO, OUT = sys.argv[1], sys.argv[2]
B, SEED = 10000, 42
band_file = json.load(open(os.path.join(OUT, 'band_assignment.json')))
CUTS = np.array(band_file['cuts_full_precision'])
ASSIGN = band_file['assignment']

z = np.load(os.path.join(REPO, 'layer_sweep/v1/p4_predictions.npz'), allow_pickle=True)
acc = [str(a) for a in z['acc']]; y = z['y'].astype(int); clu = np.array([str(c) for c in z['cluster']])
S = {'layer26_esm': z['esm_L26'].astype(np.float64),
     'composition': z['comp_L26'].astype(np.float64)}
assert len(acc) == 2366 and len(set(acc)) == 2366

band = np.array([ASSIGN[a] for a in acc], int)
assert band.min() == 0 and band.max() == 4
assert all(a in ASSIGN for a in acc), 'a P4 protein has no band'
assert len(band) == 2366, 'band vector length'
cnt = collections.Counter(band.tolist())
assert sum(cnt.values()) == 2366 and len(cnt) == 5, 'every protein in exactly one band'

def auc(yy, ss):
    """Tie-aware AUROC via mid-ranks. Returns None if either class is absent."""
    npos = int(yy.sum()); nneg = len(yy) - npos
    if npos == 0 or nneg == 0: return None
    order = np.argsort(ss, kind='mergesort')
    s_sorted = ss[order]
    ranks = np.empty(len(ss), float)
    i = 0
    while i < len(ss):
        j = i
        while j + 1 < len(ss) and s_sorted[j + 1] == s_sorted[i]: j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return (ranks[yy == 1].sum() - npos * (npos + 1) / 2.0) / (npos * nneg)

def weighted_within_band(yy, ss, bb):
    """Pair-weighted mean of within-band AUROCs. Returns (value, total_weight)."""
    num = den = 0.0
    for b in range(5):
        m = bb == b
        a = auc(yy[m], ss[m])
        if a is None: continue
        w = float(yy[m].sum()) * float((1 - yy[m]).sum())
        num += w * a; den += w
    return (num / den if den > 0 else None), den

# ---------------- descriptive rows ----------------
rows = []
for b in range(5):
    m = band == b
    nex = int(y[m].sum()); ncy = int((1 - y[m]).sum())
    rows.append({
        'band': b + 1,
        'lower': '' if b == 0 else f'{CUTS[b-1]:.10f}',
        'upper': '' if b == 4 else f'{CUTS[b]:.10f}',
        'n': int(m.sum()), 'n_extracellular': nex, 'n_cytoplasm': ncy,
        'n_eligible_pairs': nex * ncy,
        'pair_weight_share_pct': None,
        'n_clusters_id30': len({c for c in clu[m]}),
        'largest_cluster_n': max(collections.Counter(clu[m]).values()),
        'auroc_layer26_esm': round(auc(y[m], S['layer26_esm'][m]), 6),
        'auroc_composition': round(auc(y[m], S['composition'][m]), 6),
    })
tot_pairs = sum(r['n_eligible_pairs'] for r in rows)
for r in rows:
    r['pair_weight_share_pct'] = round(100 * r['n_eligible_pairs'] / tot_pairs, 3)
    r['delta_auroc'] = round(r['auroc_layer26_esm'] - r['auroc_composition'], 6)

# ---------------- point estimate ----------------
pt = {}
for nm, s in S.items():
    v, w = weighted_within_band(y, s, band)
    pt[nm] = v
    assert abs(w - tot_pairs) < 1e-6
delta = pt['layer26_esm'] - pt['composition']

# ---------------- paired cluster bootstrap ----------------
gs = sorted(set(clu.tolist()))
gi = {g: np.where(clu == g)[0] for g in gs}
rng = np.random.default_rng(SEED)
diffs, undef, degenerate_bands = [], 0, collections.Counter()
for _ in range(B):
    d = rng.choice(len(gs), size=len(gs), replace=True)
    sel = np.concatenate([gi[gs[x]] for x in d])
    yy, bb = y[sel], band[sel]
    a, aw = weighted_within_band(yy, S['layer26_esm'][sel], bb)
    c, cw = weighted_within_band(yy, S['composition'][sel], bb)
    if a is None or c is None:
        undef += 1; continue
    for b in range(5):
        m = bb == b
        if m.sum() == 0 or yy[m].sum() == 0 or (1 - yy[m]).sum() == 0:
            degenerate_bands[b + 1] += 1
    diffs.append(a - c)
diffs = np.array(diffs)
lo, hi = np.percentile(diffs, [2.5, 97.5])

res = {
 'label': 'FOLLOW-UP on previously examined P4 development data; not a fresh test set',
 'saved_predictions_only': True, 'refit': False, 'extraction': False, 'P0_accessed': False,
 'band_definition': {
   'property': 'charge_balance_30 = (K+R-D-E)/window_length over residues 1..min(30,L)',
   'cuts_defined_on': 'discovery P1+P2 observed values only, quantiles 0.20/0.40/0.60/0.80',
   'cuts_full_precision': [float(c) for c in CUTS],
   'cuts_exact': ['-2/30', '-1/30', '1/30', '2/30'],
   'interval_rule': 'lower inclusive, upper exclusive; outer intervals unbounded',
   'every_P4_protein_in_exactly_one_band': True,
   'n_P4_with_missing_charge': 0,
   'band_counts': {str(b + 1): rows[b]['n'] for b in range(5)},
   'full_precision_vs_6dp_disagreements': 0},
 'primary': {
   'definition': 'pair-weighted mean of within-band AUROCs; only extracellular-cytoplasmic '
                 'pairs from the SAME charge band contribute',
   'weight': 'w_b = n_extracellular_b * n_cytoplasmic_b',
   'total_eligible_pairs': int(tot_pairs),
   'auroc_layer26_esm': round(pt['layer26_esm'], 6),
   'auroc_composition': round(pt['composition'], 6),
   'delta_auroc_esm_minus_composition': round(delta, 6),
   'ci95_lo': round(float(lo), 6), 'ci95_hi': round(float(hi), 6),
   'ci95_excludes_zero': int(lo > 0 or hi < 0)},
 'bootstrap': {'unit': '30% identity cluster across all of P4, drawn whole with replacement',
   'same_draw_for_both_models_and_all_bands': True,
   'pair_weights_recomputed_each_draw': True,
   'B': B, 'seed': SEED, 'interval': 'percentile 95%',
   'n_clusters_resampled': len(gs),
   'n_resamples_undefined_excluded': undef,
   'n_valid_resamples': int(len(diffs)),
   'draws_in_which_a_band_was_degenerate': {str(k): v for k, v in sorted(degenerate_bands.items())}},
 'bands': rows,
}
json.dump(res, open(os.path.join(OUT, 'charge_band_results.json'), 'w'), indent=1)
cols = ['band','lower','upper','n','n_extracellular','n_cytoplasm','n_eligible_pairs',
        'pair_weight_share_pct','n_clusters_id30','largest_cluster_n',
        'auroc_layer26_esm','auroc_composition','delta_auroc']
with open(os.path.join(OUT, 'charge_band_results.tsv'), 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in rows: f.write('\t'.join(str(r[c]) for c in cols) + '\n')
    f.write('\nPRIMARY (pair-weighted mean of within-band AUROCs)\n')
    f.write('model\tauroc\n')
    f.write(f"layer26_esm\t{pt['layer26_esm']:.6f}\ncomposition\t{pt['composition']:.6f}\n")
    f.write(f"delta_esm_minus_composition\t{delta:.6f}\n")
    f.write(f"ci95_lo\t{lo:.6f}\nci95_hi\t{hi:.6f}\n")
    f.write(f"n_clusters\t{len(gs)}\nn_resamples\t{B}\nundefined_excluded\t{undef}\n")

print(f"{'band':>5}{'n':>6}{'EX':>5}{'CY':>5}{'pairs':>9}{'wt%':>8}{'clus':>6}{'big':>5}"
      f"{'ESM':>10}{'comp':>10}{'delta':>9}")
for r in rows:
    print(f"{r['band']:>5}{r['n']:>6}{r['n_extracellular']:>5}{r['n_cytoplasm']:>5}"
          f"{r['n_eligible_pairs']:>9}{r['pair_weight_share_pct']:>7.1f}%{r['n_clusters_id30']:>6}"
          f"{r['largest_cluster_n']:>5}{r['auroc_layer26_esm']:>10.4f}"
          f"{r['auroc_composition']:>10.4f}{r['delta_auroc']:>+9.4f}")
print()
print(f"PRIMARY  ESM={pt['layer26_esm']:.6f}  composition={pt['composition']:.6f}")
print(f"         delta={delta:+.6f}  95% CI [{lo:+.6f}, {hi:+.6f}]  excludes_zero={int(lo>0 or hi<0)}")
print(f"         clusters={len(gs)}  undefined resamples={undef}  valid={len(diffs)}")
print(f"         degenerate-band draws: {dict(sorted(degenerate_bands.items()))}")
