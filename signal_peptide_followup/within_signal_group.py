#!/usr/bin/env python3
"""Within-signal-peptide-group comparison of the frozen layer-26 ESM-only probe and the
matched composition baseline on P4. Saved predictions only; no refitting, no extraction,
no P0 access. Rules frozen in frozen_rules.json before this script was run.

Primary: only extracellular-cytoplasmic pairs from the SAME annotation group
(annotated_present or not_annotated) contribute. Group AUROCs are averaged weighted by
w_g = n_extracellular_g * n_cytoplasm_g. 'unverifiable' is kept separate.

Interval: whole 30%-identity clusters resampled across all of P4, the SAME draw for both
models and both groups, pair weights recomputed in every draw (as in the charge and
hydrophobicity follow-ups).

Usage: within_signal_group.py <repo_root> <out_dir>
"""
import csv, json, os, sys, hashlib, collections
import numpy as np

REPO, OUT = sys.argv[1], sys.argv[2]
RULES = os.path.join(OUT, 'frozen_rules.json')
rules_sha = hashlib.sha256(open(RULES, 'rb').read()).hexdigest()
assert rules_sha == open(os.path.join(OUT, 'frozen_rules.sha256')).read().split()[0], 'rules changed'
B, SEED = 10000, 42
GROUPS = ['annotated_present', 'not_annotated']

z = np.load(os.path.join(REPO, 'layer_sweep/v1/p4_predictions.npz'), allow_pickle=True)
acc = [str(a) for a in z['acc']]; y = z['y'].astype(int); clu = np.array([str(c) for c in z['cluster']])
S = {'layer26_esm': z['esm_L26'].astype(np.float64), 'composition': z['comp_L26'].astype(np.float64)}
assert len(acc) == 2366 and len(set(acc)) == 2366

man = {r['ACC']: r for r in csv.DictReader(open(os.path.join(
    REPO, 'discovery_analysis/v1/inclusion_manifest.tsv')), delimiter='\t')}
for i, a in enumerate(acc):
    r = man[a]
    assert r['in_core'] == '1' and r['partition'] == '4', a
    assert int(r['target']) == y[i], a
status = np.array([man[a]['signal_status'] for a in acc])
status_counts = collections.Counter(status.tolist())
assert set(status_counts) <= set(GROUPS + ['unverifiable']), status_counts
gid = np.array([GROUPS.index(s) if s in GROUPS else -1 for s in status])

def auc(yy, ss):
    """Tie-aware AUROC via mid-ranks (identical to within_charge_band.py)."""
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

def weighted_within(yy, ss, gg):
    num = den = 0.0
    for g in range(len(GROUPS)):
        m = gg == g
        a = auc(yy[m], ss[m])
        if a is None: continue
        w = float(yy[m].sum()) * float((1 - yy[m]).sum())
        num += w * a; den += w
    return (num / den if den > 0 else None), den

def describe(m, name):
    nex = int(y[m].sum()); ncy = int((1 - y[m]).sum())
    a_e = auc(y[m], S['layer26_esm'][m]); a_c = auc(y[m], S['composition'][m])
    cex = len(set(clu[m & (y == 1)])); ccy = len(set(clu[m & (y == 0)]))
    return {'group': name, 'n': int(m.sum()), 'n_extracellular': nex, 'n_cytoplasm': ncy,
            'n_eligible_pairs': nex * ncy, 'pair_weight_share_pct': None,
            'n_clusters_id30': len(set(clu[m])), 'n_clusters_extracellular': cex,
            'n_clusters_cytoplasm': ccy,
            'largest_cluster_n': max(collections.Counter(clu[m]).values()) if m.any() else 0,
            'auroc_layer26_esm': None if a_e is None else round(a_e, 6),
            'auroc_composition': None if a_c is None else round(a_c, 6),
            'delta_auroc': None if a_e is None else round(a_e - a_c, 6),
            'sparse': bool(min(nex, ncy) < 30 or min(cex, ccy) < 20)}

rows = [describe(gid == g, GROUPS[g]) for g in range(len(GROUPS))]
tot = sum(r['n_eligible_pairs'] for r in rows)
for r in rows: r['pair_weight_share_pct'] = round(100 * r['n_eligible_pairs'] / tot, 3)
unv = describe(status == 'unverifiable', 'unverifiable (separate, not in primary)')
unv['pair_weight_share_pct'] = ''

pt = {nm: weighted_within(y, s, gid)[0] for nm, s in S.items()}
delta = pt['layer26_esm'] - pt['composition']

gs = sorted(set(clu.tolist())); gi = {g: np.where(clu == g)[0] for g in gs}
rng = np.random.default_rng(SEED)
diffs, undef, degenerate = [], 0, collections.Counter()
for _ in range(B):
    d = rng.choice(len(gs), size=len(gs), replace=True)
    sel = np.concatenate([gi[gs[x]] for x in d])
    yy, gg = y[sel], gid[sel]
    a, _ = weighted_within(yy, S['layer26_esm'][sel], gg)
    c, _ = weighted_within(yy, S['composition'][sel], gg)
    if a is None or c is None: undef += 1; continue
    for g in range(len(GROUPS)):
        m = gg == g
        if m.sum() == 0 or yy[m].sum() == 0 or (1 - yy[m]).sum() == 0: degenerate[GROUPS[g]] += 1
    diffs.append(a - c)
diffs = np.array(diffs); lo, hi = np.percentile(diffs, [2.5, 97.5])

res = {'label': 'FOLLOW-UP on previously analysed P4 development data; not a fresh test set',
       'frozen_rules_sha256': rules_sha,
       'saved_predictions_only': True, 'refit': False, 'extraction': False, 'P0_accessed': False,
       'status_counts_P4': dict(status_counts),
       'groups': rows, 'unverifiable': unv,
       'primary': {'auroc_layer26_esm': round(pt['layer26_esm'], 6),
                   'auroc_composition': round(pt['composition'], 6),
                   'delta_auroc_esm_minus_composition': round(delta, 6),
                   'ci95_lo': round(float(lo), 6), 'ci95_hi': round(float(hi), 6),
                   'ci95_excludes_zero': int(lo > 0 or hi < 0), 'total_eligible_pairs': int(tot)},
       'bootstrap': {'B': B, 'seed': SEED, 'n_clusters_resampled': len(gs),
                     'n_resamples_undefined_excluded': undef, 'n_valid_resamples': int(len(diffs)),
                     'draws_with_degenerate_group': dict(degenerate)}}
json.dump(res, open(os.path.join(OUT, 'signal_group_results.json'), 'w'), indent=1)
cols = ['group', 'n', 'n_extracellular', 'n_cytoplasm', 'n_eligible_pairs', 'pair_weight_share_pct',
        'n_clusters_id30', 'n_clusters_extracellular', 'n_clusters_cytoplasm', 'largest_cluster_n',
        'auroc_layer26_esm', 'auroc_composition', 'delta_auroc', 'sparse']
with open(os.path.join(OUT, 'signal_group_results.tsv'), 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in rows + [unv]: f.write('\t'.join(str(r[c]) for c in cols) + '\n')
    f.write(f"\nPRIMARY\tESM={pt['layer26_esm']:.6f}\tcomp={pt['composition']:.6f}\t"
            f"delta={delta:.6f}\tci95=[{lo:.6f},{hi:.6f}]\n")
print(json.dumps(res, indent=1))
