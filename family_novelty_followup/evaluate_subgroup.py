#!/usr/bin/env python3
"""STEP 2: evaluate the FROZEN layer-26 ESM-only probe and the matched composition
baseline on the frozen subgroup. Saved predictions only - no refitting, no extraction.

Usage: evaluate_subgroup.py <repo_root> <out_dir>
"""
import csv, json, os, sys, collections
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score

REPO, OUT = sys.argv[1], sys.argv[2]
B, SEED = 10000, 42                       # identical to the sweep's frozen protocol

frozen = json.load(open(os.path.join(OUT, 'subgroup_frozen.json')))
mem = list(csv.DictReader(open(os.path.join(OUT, 'subgroup_membership.tsv')), delimiter='\t'))
assert len(mem) == frozen['n_members'] == 443

z = np.load(os.path.join(REPO, 'layer_sweep/v1/p4_predictions.npz'), allow_pickle=True)
acc = list(z['acc']); y_all = z['y']; clu_all = z['cluster']
pos = {a: i for i, a in enumerate(acc)}
assert len(pos) == len(acc) == 2366

idx = np.array([pos[r['ACC']] for r in mem])
y = y_all[idx]
clu = clu_all[idx]
# the frozen table's target and cluster must agree with the saved prediction file
assert all(int(r['target']) == int(y[i]) for i, r in enumerate(mem)), 'target disagreement'
assert all(r['cluster_id30'] == clu[i] for i, r in enumerate(mem)), 'cluster disagreement'
assert int(y.sum()) == 183 and int((1 - y).sum()) == 260

S = {'layer26_esm': z['esm_L26'][idx], 'composition': z['comp_L26'][idx]}
FULL = {'layer26_esm': z['esm_L26'], 'composition': z['comp_L26']}

def pboot(yy, a, b, g, n=B, seed=SEED):
    gs = sorted(set(g)); gi = {x: np.where(g == x)[0] for x in gs}
    rng = np.random.default_rng(seed)
    dau, dap, und = [], [], 0
    for _ in range(n):
        d = rng.choice(len(gs), size=len(gs), replace=True)
        s = np.concatenate([gi[gs[x]] for x in d])
        if len(set(yy[s])) < 2: und += 1; continue
        dau.append(roc_auc_score(yy[s], a[s]) - roc_auc_score(yy[s], b[s]))
        dap.append(average_precision_score(yy[s], a[s]) - average_precision_score(yy[s], b[s]))
    return np.array(dau), np.array(dap), und

cc = collections.Counter(clu); big, bign = cc.most_common(1)[0]
rows = []
for name, s in S.items():
    rows.append({'analysis': 'subgroup_443_family_unseen_vs_P1P2P3', 'model': name,
                 'n': len(y), 'n_extracellular': int(y.sum()), 'n_cytoplasm': int((1 - y).sum()),
                 'n_clusters_id30': len(cc), 'largest_cluster': big, 'largest_cluster_n': bign,
                 'largest_cluster_pct': round(100 * bign / len(y), 3),
                 'auroc_protein_weighted': round(float(roc_auc_score(y, s)), 6),
                 'ap_protein_weighted': round(float(average_precision_score(y, s)), 6)})
dau, dap, und = pboot(y, S['layer26_esm'], S['composition'], clu)
d_au = float(roc_auc_score(y, S['layer26_esm']) - roc_auc_score(y, S['composition']))
d_ap = float(average_precision_score(y, S['layer26_esm']) - average_precision_score(y, S['composition']))
diff = {'analysis': 'subgroup_443_family_unseen_vs_P1P2P3',
        'comparison': 'layer26_esm_minus_composition', 'n': len(y),
        'n_clusters_id30': len(cc), 'largest_cluster_pct': round(100 * bign / len(y), 3),
        'delta_auroc': round(d_au, 6),
        'auroc_ci95_lo': round(float(np.percentile(dau, 2.5)), 6),
        'auroc_ci95_hi': round(float(np.percentile(dau, 97.5)), 6),
        'auroc_ci95_excludes_zero': int(np.percentile(dau, 2.5) > 0 or np.percentile(dau, 97.5) < 0),
        'delta_ap': round(d_ap, 6),
        'ap_ci95_lo': round(float(np.percentile(dap, 2.5)), 6),
        'ap_ci95_hi': round(float(np.percentile(dap, 97.5)), 6),
        'ap_ci95_excludes_zero': int(np.percentile(dap, 2.5) > 0 or np.percentile(dap, 97.5) < 0),
        'n_resamples': B, 'n_resamples_undefined_excluded': und, 'seed': SEED}

# context: the existing full-P4 figures, recomputed from the same saved file
ccf = collections.Counter(clu_all); bigf, bignf = ccf.most_common(1)[0]
ctx = []
for name, s in FULL.items():
    ctx.append({'analysis': 'full_P4_context', 'model': name, 'n': len(y_all),
                'n_extracellular': int(y_all.sum()), 'n_cytoplasm': int((1 - y_all).sum()),
                'n_clusters_id30': len(ccf), 'largest_cluster': bigf, 'largest_cluster_n': bignf,
                'largest_cluster_pct': round(100 * bignf / len(y_all), 3),
                'auroc_protein_weighted': round(float(roc_auc_score(y_all, s)), 6),
                'ap_protein_weighted': round(float(average_precision_score(y_all, s)), 6)})
ctx.append({'analysis': 'full_P4_context', 'model': 'layer26_esm_minus_composition',
            'n': len(y_all), 'n_extracellular': int(y_all.sum()),
            'n_cytoplasm': int((1 - y_all).sum()), 'n_clusters_id30': len(ccf),
            'largest_cluster': bigf, 'largest_cluster_n': bignf,
            'largest_cluster_pct': round(100 * bignf / len(y_all), 3),
            'auroc_protein_weighted': round(float(
                roc_auc_score(y_all, FULL['layer26_esm']) - roc_auc_score(y_all, FULL['composition'])), 6),
            'ap_protein_weighted': round(float(
                average_precision_score(y_all, FULL['layer26_esm']) -
                average_precision_score(y_all, FULL['composition'])), 6)})

cols = list(rows[0])
with open(os.path.join(OUT, 'subgroup_results.tsv'), 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in rows + ctx: f.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')
    f.write('\n')
    dcols = list(diff)
    f.write('\t'.join(dcols) + '\n')
    f.write('\t'.join(str(diff[c]) for c in dcols) + '\n')
json.dump({'label': 'FOLLOW-UP analysis on previously examined P4 data; not a fresh test set',
           'refitting': False, 'extraction': False, 'P0_accessed': False,
           'subgroup': rows, 'difference': diff, 'full_P4_context': ctx,
           'frozen_membership_sha256': frozen['membership_sha256'],
           'bootstrap': {'unit': '30% identity cluster, drawn whole with replacement',
                         'pairing': 'same resample indices for both models',
                         'B': B, 'seed': SEED, 'interval': 'percentile 95%',
                         'mixed_class_clusters': 'kept intact, never split',
                         'stratified': False}},
          open(os.path.join(OUT, 'subgroup_results.json'), 'w'), indent=1)

print(f"{'model':<22}{'n':>5}{'EX':>5}{'CY':>5}{'clusters':>10}{'AUROC':>10}{'AP':>10}")
for r in rows:
    print(f"{r['model']:<22}{r['n']:>5}{r['n_extracellular']:>5}{r['n_cytoplasm']:>5}"
          f"{r['n_clusters_id30']:>10}{r['auroc_protein_weighted']:>10.4f}{r['ap_protein_weighted']:>10.4f}")
print()
print(f"paired difference (ESM - composition): AUROC {diff['delta_auroc']:+.6f} "
      f"[{diff['auroc_ci95_lo']:+.6f}, {diff['auroc_ci95_hi']:+.6f}] excl0={diff['auroc_ci95_excludes_zero']}")
print(f"                                       AP    {diff['delta_ap']:+.6f} "
      f"[{diff['ap_ci95_lo']:+.6f}, {diff['ap_ci95_hi']:+.6f}] excl0={diff['ap_ci95_excludes_zero']}")
print(f"clusters={diff['n_clusters_id30']}  largest share={diff['largest_cluster_pct']}%  "
      f"undefined resamples={diff['n_resamples_undefined_excluded']}")
print()
print('--- full P4 context ---')
for r in ctx:
    print(f"  {r['model']:<34} n={r['n']} AUROC={r['auroc_protein_weighted']:+.6f} AP={r['ap_protein_weighted']:+.6f}")
