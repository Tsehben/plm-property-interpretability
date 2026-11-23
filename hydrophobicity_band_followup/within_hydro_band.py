#!/usr/bin/env python3
"""Within-hydrophobicity-band comparison of the frozen layer-26 ESM-only probe and the
matched composition baseline on P4. Saved predictions only; no refitting, no extraction.

Identical analysis to the charge-band comparison: only extracellular-cytoplasmic pairs from
the SAME band contribute; band AUROCs are averaged with weights w_b = n_EX_b * n_CY_b; the
interval comes from whole 30%-identity cluster draws applied identically to both models and
all five bands, with pair weights recomputed inside each draw.

Usage: within_hydro_band.py <repo_root> <out_dir>
"""
import csv, json, os, sys, collections
import numpy as np

REPO, OUT = sys.argv[1], sys.argv[2]
B, SEED, WINDOW = 10000, 42, 30
KD = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5,
      'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
      'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}

BA = json.load(open(os.path.join(OUT, 'band_assignment.json')))
CUTS_STR = BA['cuts_exact_str']; CUTS = np.array(BA['cuts_float'])
ASSIGN = BA['assignment']
# Bands come from EXACT RATIONAL arithmetic, not float64. Kyte-Doolittle values are one
# decimal, so every first-30 mean is an exact multiple of 1/300 and every cut is an exact
# rational (-5/6, -41/100, 3/100, 173/300). In float64 the summed mean and the quantile
# carry different rounding error, so a protein sitting exactly ON a cut (e.g. true value
# 0.03 computed as 0.029999999999999995 against a cut of 0.03000000000000007) lands in the
# wrong band. Eight P4 proteins were affected. Fractions make the lower-inclusive /
# upper-exclusive rule exact at the boundary.

z = np.load(os.path.join(REPO, 'layer_sweep/v1/p4_predictions.npz'), allow_pickle=True)
acc = [str(a) for a in z['acc']]; y_all = z['y'].astype(int)
clu_all = np.array([str(c) for c in z['cluster']])
S_all = {'layer26_esm': z['esm_L26'].astype(np.float64),
         'composition': z['comp_L26'].astype(np.float64)}
assert len(acc) == 2366

missing = [a for a in acc if a not in ASSIGN]
keep = np.array([i for i, a in enumerate(acc) if a in ASSIGN])
y, clu = y_all[keep], clu_all[keep]
S = {k: v[keep] for k, v in S_all.items()}
band = np.array([ASSIGN[acc[i]] for i in keep], int)
cnt = collections.Counter(band.tolist())
assert len(keep) == 2365, f'expected 2,365 eligible, got {len(keep)}'
assert len(missing) == 1, f'expected 1 missing value, got {len(missing)}'
assert sum(cnt.values()) == 2365 and set(cnt) == {0, 1, 2, 3, 4}, 'every protein in exactly one band'

def auc(yy, ss):
    npos = int(yy.sum()); nneg = len(yy) - npos
    if npos == 0 or nneg == 0: return None
    order = np.argsort(ss, kind='mergesort'); s_sorted = ss[order]
    ranks = np.empty(len(ss), float); i = 0
    while i < len(ss):
        j = i
        while j + 1 < len(ss) and s_sorted[j + 1] == s_sorted[i]: j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return (ranks[yy == 1].sum() - npos * (npos + 1) / 2.0) / (npos * nneg)

def weighted(yy, ss, bb):
    num = den = 0.0
    for b in range(5):
        m = bb == b
        a = auc(yy[m], ss[m])
        if a is None: continue
        w = float(yy[m].sum()) * float((1 - yy[m]).sum())
        num += w * a; den += w
    return (num / den if den > 0 else None), den

rows = []
for b in range(5):
    m = band == b
    nex = int(y[m].sum()); ncy = int((1 - y[m]).sum())
    rows.append({'band': b + 1,
                 'lower': '' if b == 0 else f'{CUTS[b-1]:.10f}',
                 'upper': '' if b == 4 else f'{CUTS[b]:.10f}',
                 'n': int(m.sum()), 'n_extracellular': nex, 'n_cytoplasm': ncy,
                 'n_eligible_pairs': nex * ncy, 'pair_weight_share_pct': None,
                 'n_clusters_id30': len(set(clu[m])),
                 'largest_cluster_n': max(collections.Counter(clu[m]).values()),
                 'auroc_layer26_esm': round(auc(y[m], S['layer26_esm'][m]), 6),
                 'auroc_composition': round(auc(y[m], S['composition'][m]), 6),
                 'support': 'sparse' if nex < 20 else 'ok'})
tot = sum(r['n_eligible_pairs'] for r in rows)
assert tot == 90815, f'expected 90,815 within-band pairs, got {tot}'
for r in rows:
    r['pair_weight_share_pct'] = round(100 * r['n_eligible_pairs'] / tot, 3)
    r['delta_auroc'] = round(r['auroc_layer26_esm'] - r['auroc_composition'], 6)

pt = {}
for nm, s in S.items():
    v, w = weighted(y, s, band); pt[nm] = v
    assert abs(w - tot) < 1e-6
delta = pt['layer26_esm'] - pt['composition']

gs = sorted(set(clu.tolist())); gi = {g: np.where(clu == g)[0] for g in gs}
rng = np.random.default_rng(SEED)
diffs, undef, degen = [], 0, collections.Counter()
for _ in range(B):
    d = rng.choice(len(gs), size=len(gs), replace=True)
    sel = np.concatenate([gi[gs[x]] for x in d])
    yy, bb = y[sel], band[sel]
    a, _ = weighted(yy, S['layer26_esm'][sel], bb)
    c, _ = weighted(yy, S['composition'][sel], bb)
    if a is None or c is None: undef += 1; continue
    for b in range(5):
        m = bb == b
        if m.sum() == 0 or yy[m].sum() == 0 or (1 - yy[m]).sum() == 0: degen[b + 1] += 1
    diffs.append(a - c)
diffs = np.array(diffs); lo, hi = np.percentile(diffs, [2.5, 97.5])

res = {'label': 'FOLLOW-UP on previously examined P4 development data; not a fresh test set',
 'property': 'mean Kyte-Doolittle hydrophobicity over residues 1..min(30,L)',
 'saved_predictions_only': True, 'refit': False, 'extraction': False, 'P0_accessed': False,
 'band_definition': {'cuts_defined_on': 'discovery P1+P2 observed values only, quantiles .2/.4/.6/.8',
   'cuts_exact_rational': CUTS_STR,
   'cuts_float_view': [float(c) for c in CUTS],
   'arithmetic': 'EXACT rational (Fraction); float64 misassigns 8 P4 proteins that sit '
                 'exactly on a cut, because the summed mean and the quantile carry '
                 'different rounding error',
   'interval_rule': 'lower inclusive, upper exclusive; outer intervals unbounded',
   'n_P4_total': 2366, 'n_eligible': int(len(keep)),
   'n_missing_excluded': len(missing), 'missing_accessions': missing,
   'missing_reason': 'Kyte-Doolittle undefined for a residue in the first-30 window; the protein is kept in the cohort and excluded only from this comparison',
   'every_eligible_protein_in_exactly_one_band': True,
   'band_counts': {str(b + 1): rows[b]['n'] for b in range(5)}},
 'primary': {'definition': 'pair-weighted mean of within-band AUROCs; only extracellular-'
                           'cytoplasmic pairs from the SAME band contribute',
   'weight': 'w_b = n_extracellular_b * n_cytoplasmic_b',
   'total_eligible_pairs': int(tot),
   'auroc_layer26_esm': round(pt['layer26_esm'], 6),
   'auroc_composition': round(pt['composition'], 6),
   'delta_auroc_esm_minus_composition': round(delta, 6),
   'ci95_lo': round(float(lo), 6), 'ci95_hi': round(float(hi), 6),
   'ci95_excludes_zero': int(lo > 0 or hi < 0)},
 'bootstrap': {'unit': '30% identity cluster, drawn whole with replacement',
   'same_draw_for_both_models_and_all_bands': True, 'pair_weights_recomputed_each_draw': True,
   'B': B, 'seed': SEED, 'interval': 'percentile 95%', 'n_clusters_resampled': len(gs),
   'n_resamples_undefined_excluded': undef, 'n_valid_resamples': int(len(diffs)),
   'draws_in_which_a_band_was_degenerate': {str(k): v for k, v in sorted(degen.items())}},
 'bands': rows}
json.dump(res, open(os.path.join(OUT, 'hydro_band_results.json'), 'w'), indent=1)
cols = ['band','lower','upper','n','n_extracellular','n_cytoplasm','n_eligible_pairs',
        'pair_weight_share_pct','n_clusters_id30','largest_cluster_n',
        'auroc_layer26_esm','auroc_composition','delta_auroc','support']
with open(os.path.join(OUT, 'hydro_band_results.tsv'), 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in rows: f.write('\t'.join(str(r[c]) for c in cols) + '\n')
    f.write('\nPRIMARY (pair-weighted mean of within-band AUROCs)\nmodel\tauroc\n')
    f.write(f"layer26_esm\t{pt['layer26_esm']:.6f}\ncomposition\t{pt['composition']:.6f}\n")
    f.write(f"delta_esm_minus_composition\t{delta:.6f}\nci95_lo\t{lo:.6f}\nci95_hi\t{hi:.6f}\n")
    f.write(f"n_clusters\t{len(gs)}\nn_resamples\t{B}\nundefined_excluded\t{undef}\n")

print('cuts (exact):', CUTS_STR, '=', [float(c) for c in CUTS])
print(f'eligible={len(keep)}  missing={len(missing)} {missing}  within-band pairs={tot}')
print()
print(f"{'band':>5}{'n':>6}{'EX':>5}{'CY':>5}{'pairs':>9}{'wt%':>8}{'clus':>6}{'big':>5}"
      f"{'ESM':>10}{'comp':>10}{'delta':>9}  support")
for r in rows:
    print(f"{r['band']:>5}{r['n']:>6}{r['n_extracellular']:>5}{r['n_cytoplasm']:>5}"
          f"{r['n_eligible_pairs']:>9}{r['pair_weight_share_pct']:>7.1f}%{r['n_clusters_id30']:>6}"
          f"{r['largest_cluster_n']:>5}{r['auroc_layer26_esm']:>10.4f}"
          f"{r['auroc_composition']:>10.4f}{r['delta_auroc']:>+9.4f}  {r['support']}")
print()
print(f"PRIMARY  ESM={pt['layer26_esm']:.6f}  composition={pt['composition']:.6f}")
print(f"         delta={delta:+.6f}  95% CI [{lo:+.6f}, {hi:+.6f}]  excludes_zero={int(lo>0 or hi<0)}")
print(f"         clusters={len(gs)}  undefined={undef}  valid={len(diffs)}  degenerate-band draws={dict(sorted(degen.items()))}")
