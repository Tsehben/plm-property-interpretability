#!/usr/bin/env python3
"""STEP 1: freeze the subgroup membership. NO PREDICTION FILE IS OPENED HERE.

Subgroup: P4 (evaluation) proteins whose applicable Pfam FAMILY-type signatures were ALL
absent from P1+P2+P3 — i.e. classification 'all_unseen' relative to the probe's fitting
partitions (P1+P2) plus the layer-selection partition (P3).

"Unseen" refers strictly to recorded Pfam family IDs in those partitions under this
annotation snapshot. It is not a claim of evolutionary independence.

Usage: freeze_subgroup.py <repo_root> <out_dir>
"""
import csv, json, os, sys, hashlib, collections

REPO, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
csv.field_size_limit(10_000_000)
FAMILY_TYPE = 'family'

for forbidden in ('p4_predictions.npz', 'p3_predictions.npz', 'probe_parameters.npz'):
    assert forbidden not in ' '.join(sys.argv), 'no prediction file may be referenced here'

T = {r['ACC']: r for r in csv.DictReader(
    open(os.path.join(REPO, 'concept_table/concept_table.tsv')), delimiter='\t')}
PF = list(csv.DictReader(
    open(os.path.join(REPO, 'concept_table/concept_pfam_long.tsv')), delimiter='\t'))
man = {r['ACC']: r for r in csv.DictReader(
    open(os.path.join(REPO, 'discovery_analysis/v1/inclusion_manifest.tsv')), delimiter='\t')
    if r['in_core'] == '1'}
seqs = {}
with open(os.path.join(REPO, 'DeepLoc2/esm2_probe/deeploc_pairs.tsv')) as f:
    next(f)
    for line in f:
        k, s = line.rstrip('\n').split('\t'); seqs[k] = s

role = {a: r['role'] for a, r in T.items()}
applicable = {a: r['domain_annotation_applies_to'] == 'entry_matches_frozen_sequence'
              for a, r in T.items()}
fam_rows = [r for r in PF if r['pfam_entry_type'] == FAMILY_TYPE]
fam_by_acc = collections.defaultdict(list)
for r in fam_rows: fam_by_acc[r['ACC']].append(r['pfam'])
seen = {r['pfam'] for r in fam_rows if role[r['ACC']] in ('discovery', 'tuning')}

members = []
for a, r in T.items():
    if r['role'] != 'evaluation': continue
    if not applicable[a]: continue
    fams = fam_by_acc.get(a, [])
    if not fams: continue
    if any(f in seen for f in fams): continue
    members.append(a)
members.sort()

# ---- per-row verification against the frozen manifest and the actual sequence ----
bad = []
for a in members:
    r, m = T[a], man[a]
    if a not in man: bad.append((a, 'absent_from_core_manifest')); continue
    if r['frozen_sha256'] != m['frozen_sha256']: bad.append((a, 'sha_mismatch_vs_manifest'))
    if hashlib.sha256(seqs[a].encode()).hexdigest() != r['frozen_sha256']:
        bad.append((a, 'sha_mismatch_vs_sequence'))
    if r['target'] != m['target']: bad.append((a, 'target_mismatch'))
    if r['partition'] != m['partition'] or m['partition'] != '4': bad.append((a, 'partition_mismatch'))
    if r['role'] != 'evaluation': bad.append((a, 'role_mismatch'))
assert not bad, bad[:20]

ex = sum(1 for a in members if T[a]['target'] == '1')
cy = sum(1 for a in members if T[a]['target'] == '0')
clu = collections.Counter(T[a]['cluster_id30'] for a in members)
big, bign = clu.most_common(1)[0]

cols = ['ACC', 'frozen_sha256', 'partition', 'role', 'target', 'localisation_class',
        'cluster_id30', 'cluster_size_id30_full_cohort', 'n_pfam',
        'pfam_family_ids', 'n_pfam_family_ids', 'domain_annotation_applies_to']
path = os.path.join(OUT, 'subgroup_membership.tsv')
with open(path, 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for a in members:
        r = T[a]
        f.write('\t'.join([a, r['frozen_sha256'], r['partition'], r['role'], r['target'],
                           r['localisation_class'], r['cluster_id30'], r['cluster_size_id30'],
                           r['n_pfam'], ';'.join(fam_by_acc[a]), str(len(fam_by_acc[a])),
                           r['domain_annotation_applies_to']]) + '\n')

def fsha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''): h.update(b)
    return h.hexdigest()

rec = {
 'step': '1 - membership frozen BEFORE any prediction score was read',
 'prediction_files_opened': [],
 'definition': "P4 proteins with >=1 APPLICABLE Pfam family-type signature, none of which "
               "appears in P1+P2+P3",
 'unseen_means': "recorded Pfam family IDs absent from the probe's fitting partitions (P1+P2) "
                 "and the layer-selection partition (P3), under this annotation snapshot; "
                 "NOT a claim of evolutionary independence",
 'family_type_basis': "the Pfam signature's OWN entry type, not the integrated InterPro entry",
 'n_members': len(members),
 'n_extracellular': ex, 'n_cytoplasm': cy,
 'expected_183_260_443': bool(ex == 183 and cy == 260 and len(members) == 443),
 'verification': {'rows_checked': len(members), 'failures': 0,
                  'accession_in_core_manifest': True,
                  'frozen_sha256_matches_manifest': True,
                  'frozen_sha256_matches_actual_sequence': True,
                  'target_matches_manifest': True,
                  'partition_is_P4': True, 'role_is_evaluation': True},
 'n_clusters_id30': len(clu),
 'largest_cluster': big, 'largest_cluster_n': bign,
 'largest_cluster_pct': round(100 * bign / len(members), 3),
 'n_distinct_family_ids_in_subgroup': len({f for a in members for f in fam_by_acc[a]}),
 'P0_accessed': False,
 'membership_sha256': fsha(path),
}
json.dump(rec, open(os.path.join(OUT, 'subgroup_frozen.json'), 'w'), indent=1)
print(json.dumps(rec, indent=1))
