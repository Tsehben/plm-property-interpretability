#!/usr/bin/env python3
"""RETROSPECTIVE cluster-aware reanalysis of the four completed P4 paired contrasts.

Saved outcomes only. No model is rerun, no intervention is executed, no file in
evaluation_p4/ is modified. Original intervals are retained alongside the new ones.

Original design (evaluation_p4/v1/analysis_summary.json): BCa bootstrap, unit = protein,
stratified by localisation class, 10,000 resamples, seed 20261122, Bonferroni level 0.9875
for 4 primary contrasts.

This reanalysis: percentile bootstrap, unit = 30%-identity cluster, whole clusters drawn
with replacement and never split, UNSTRATIFIED (splitting a mixed-class cluster to obtain a
class-stratified resample is forbidden by the frozen protocol). Same B and the same
four-comparison Bonferroni level are retained.

Usage: reanalyse_interventions.py <repo_root> <out_dir>
"""
import os, sys, csv, json, collections
import numpy as np

REPO, OUT = sys.argv[1], sys.argv[2]
B, SEED, BONF = 10000, 42, 0.9875
os.makedirs(OUT, exist_ok=True)

clu = {r['ACC']: r['cluster_id30'] for r in csv.DictReader(
    open(os.path.join(REPO, 'homology_audit/v1/cluster_membership.tsv')), delimiter='\t')}
per = list(csv.DictReader(open(os.path.join(REPO,
    'evaluation_p4/v1/primary_per_protein.tsv')), delimiter='\t'))
orig = {(str(r['candidate']), str(float(r['multiplier']))): r
        for r in json.load(open(os.path.join(REPO,
            'evaluation_p4/v1/analysis_summary.json')))['primary']}

def boot(vals, groups, b=B, seed=SEED):
    """Paired cluster bootstrap of the mean. Whole clusters, with replacement."""
    rng = np.random.default_rng(seed)
    gs = sorted(set(groups))
    idx = {g: np.where(np.array(groups) == g)[0] for g in gs}
    out = np.empty(b)
    for i in range(b):
        drawn = rng.choice(len(gs), size=len(gs), replace=True)
        sel = np.concatenate([idx[gs[d]] for d in drawn])
        out[i] = vals[sel].mean()
    return out

rows = []
for (cand, mult), o in sorted(orig.items()):
    sub = [r for r in per if r['candidate'] == cand and float(r['multiplier']) == float(mult)]
    assert len(sub) == o['n_proteins'], (cand, mult, len(sub), o['n_proteins'])
    v = np.array([float(r['paired_kl_difference']) for r in sub])
    g = [clu[r['ACC']] for r in sub]
    cc = collections.Counter(g)
    big, bign = cc.most_common(1)[0]
    bigmean = v[np.array(g) == big].mean()
    d = boot(v, g)
    lo, hi = np.percentile(d, [2.5, 97.5])
    blo, bhi = np.percentile(d, [(1 - BONF) / 2 * 100, (1 - (1 - BONF) / 2) * 100])
    rows.append({
        'candidate': cand, 'multiplier': mult, 'n_proteins': len(sub),
        'n_Extracellular': o['n_Extracellular'], 'n_Cytoplasm': o['n_Cytoplasm'],
        'n_clusters_id30': len(cc),
        'largest_cluster': big, 'largest_cluster_n': bign,
        'largest_cluster_pct_of_contrast': round(100 * bign / len(sub), 2),
        'largest_cluster_mean_value': f'{bigmean:.6e}',
        'largest_cluster_share_of_estimate_pct':
            round(100 * (bign * bigmean) / (len(sub) * v.mean()), 2) if v.mean() != 0 else None,
        'n_mixed_class_clusters': sum(
            1 for c in cc if len({r['localisation_class'] for r in sub if clu[r['ACC']] == c}) > 1),
        'estimate': f"{v.mean():.6e}",
        'orig_ci95_lo': o['ci95_lo'], 'orig_ci95_hi': o['ci95_hi'],
        'orig_ci95_excludes_zero': o['ci95_excludes_zero'],
        'orig_ci9875_lo': o['ci9875_lo'], 'orig_ci9875_hi': o['ci9875_hi'],
        'orig_ci9875_excludes_zero': o['ci9875_excludes_zero'],
        'cluster_ci95_lo': f'{lo:.6e}', 'cluster_ci95_hi': f'{hi:.6e}',
        'cluster_ci95_excludes_zero': int(lo > 0 or hi < 0),
        'cluster_ci9875_lo': f'{blo:.6e}', 'cluster_ci9875_hi': f'{bhi:.6e}',
        'cluster_ci9875_excludes_zero': int(blo > 0 or bhi < 0),
        'interpretation_changed_at_95': int((lo > 0 or hi < 0) != bool(o['ci95_excludes_zero'])),
        'interpretation_changed_at_9875': int((blo > 0 or bhi < 0) != bool(o['ci9875_excludes_zero'])),
    })

cols = list(rows[0])
with open(os.path.join(OUT, 'intervention_cluster_sensitivity.tsv'), 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in rows: f.write('\t'.join(str(r[c]) for c in cols) + '\n')
json.dump({'rows': rows, 'B': B, 'seed': SEED, 'bonferroni_level': BONF,
           'n_primary_contrasts': 4,
           'original_design': {'unit': 'protein', 'method': 'BCa',
                               'stratified_by': 'localisation class',
                               'seed': 20261122, 'n_resamples': 10000},
           'reanalysis_design': {'unit': '30% identity cluster', 'method': 'percentile',
                                 'stratified': False,
                                 'stratification_note': 'unstratified: mixed-class clusters '
                                 'are never split to obtain class-stratified resamples',
                                 'whole_clusters': True},
           'label': 'RETROSPECTIVE sensitivity analysis. Original intervals retained. '
                    'No model rerun.'},
          open(os.path.join(OUT, 'intervention_cluster_sensitivity.json'), 'w'), indent=1)

print(f"{'cand':>5}{'mult':>6}{'n':>5}{'clus':>6}{'big':>5}{'big%':>7}"
      f"{'orig95_excl0':>14}{'clus95_excl0':>14}{'orig9875':>10}{'clus9875':>10}{'changed':>9}")
for r in rows:
    print(f"{r['candidate']:>5}{r['multiplier']:>6}{r['n_proteins']:>5}{r['n_clusters_id30']:>6}"
          f"{r['largest_cluster_n']:>5}{r['largest_cluster_pct_of_contrast']:>6.1f}%"
          f"{r['orig_ci95_excludes_zero']:>14}{r['cluster_ci95_excludes_zero']:>14}"
          f"{r['orig_ci9875_excludes_zero']:>10}{r['cluster_ci9875_excludes_zero']:>10}"
          f"{r['interpretation_changed_at_95']:>9}")
print('\nwrote intervention_cluster_sensitivity.{tsv,json}')
