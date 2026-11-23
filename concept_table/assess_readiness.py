#!/usr/bin/env python3
"""Coverage and evaluation-feasibility assessment for the concept table.

Read-only over the frozen concept table. No probe training, no performance analysis.
Usage: assess_readiness.py <repo_root> <out_dir>
"""
import csv, json, os, sys, hashlib, collections
import numpy as np

REPO, OUT = sys.argv[1], sys.argv[2]
csv.field_size_limit(10_000_000)
QUANTILES = [.20, .40, .60, .80]          # same convention as the NESG within-property notebook
CONTINUOUS = ['hydrophobicity_kd_30', 'charge_balance_30']

T = list(csv.DictReader(open(os.path.join(OUT, 'concept_table.tsv')), delimiter='\t'))
PF = list(csv.DictReader(open(os.path.join(OUT, 'concept_pfam_long.tsv')), delimiter='\t'))
SG = list(csv.DictReader(open(os.path.join(OUT, 'concept_signal_long.tsv')), delimiter='\t'))
man = {r['ACC']: r for r in csv.DictReader(
    open(os.path.join(REPO, 'discovery_analysis/v1/inclusion_manifest.tsv')), delimiter='\t')
    if r['in_core'] == '1'}
seqs = {}
with open(os.path.join(REPO, 'DeepLoc2/esm2_probe/deeploc_pairs.tsv')) as f:
    next(f)
    for line in f:
        k, s = line.rstrip('\n').split('\t'); seqs[k] = s

R = {}
# ---------------------------------------------------------------- integrity
acc = [r['ACC'] for r in T]
ck = {
 'n_rows': len(T), 'n_distinct_accessions': len(set(acc)),
 'one_row_per_protein': len(T) == len(set(acc)) == 10104,
 'accession_set_equals_core_manifest': set(acc) == set(man),
 'pfam_long_accessions_subset_of_table': {r['ACC'] for r in PF} <= set(acc),
 'signal_long_accessions_subset_of_table': {r['ACC'] for r in SG} <= set(acc),
 'frozen_sha256_matches_manifest': all(r['frozen_sha256'] == man[r['ACC']]['frozen_sha256'] for r in T),
 'frozen_sha256_matches_actual_sequence': all(
     hashlib.sha256(seqs[r['ACC']].encode()).hexdigest() == r['frozen_sha256'] for r in T),
 'targets_unchanged': all(r['target'] == man[r['ACC']]['target'] for r in T),
 'partitions_unchanged': all(r['partition'] == man[r['ACC']]['partition'] for r in T),
 'roles_unchanged': all(r['role'] == man[r['ACC']]['role'] for r in T),
 'complete_accounting': len(T) == 10104 and set(acc) == set(man),
}
ck['pfam_long_join_is_many_to_one'] = all(
    collections.Counter((r['ACC'], r['pfam']) for r in PF)[k] == 1
    for k in set((r['ACC'], r['pfam']) for r in PF))
R['integrity_checks'] = ck

# ---------------------------------------------------------------- by role
by_role = {}
for role in ('discovery', 'tuning', 'evaluation'):
    g = [r for r in T if r['role'] == role]
    d = {'n': len(g),
         'partitions': sorted({r['partition'] for r in g}),
         'n_Extracellular': sum(1 for r in g if r['target'] == '1'),
         'n_Cytoplasm': sum(1 for r in g if r['target'] == '0')}
    for t in (30, 50, 80):
        d[f'n_clusters_id{t}'] = len({r[f'cluster_id{t}'] for r in g})
    for c in CONTINUOUS:
        n_missing = sum(1 for r in g if r[c] == '')
        d[f'{c}_missing'] = n_missing
        d[f'{c}_missing_pct'] = round(100 * n_missing / len(g), 3)
    d['pfam_no_ids_recorded'] = sum(1 for r in g if r['n_pfam'] == '0')
    d['pfam_recorded'] = sum(1 for r in g if r['pfam_annotation_status'] == 'recorded')
    d['pfam_recorded_applicability_uncertain'] = sum(
        1 for r in g if r['pfam_annotation_status'] == 'recorded_applicability_uncertain')
    d['signal_annotated_present'] = sum(1 for r in g if r['signal_status'] == 'annotated_present')
    d['signal_not_annotated'] = sum(1 for r in g if r['signal_status'] == 'not_annotated')
    d['signal_unverifiable'] = sum(1 for r in g if r['signal_status'] == 'unverifiable')
    d['signalp6_predictions'] = 0
    by_role[role] = d
R['by_role'] = by_role

# ---------------------------------------------------------------- Pfam seen/unseen
# CORRECTION: the family-level check uses the PFAM SIGNATURE'S OWN entry type
# (InterPro /api/entry/pfam/ -> metadata.type), NOT the type of the InterPro entry it is
# integrated into. A Pfam *domain* can be integrated into an InterPro entry labelled
# *family* (e.g. PF00001 is a Pfam domain inside InterPro IPR000276, a family), so using
# the integrated entry's label would silently promote domains to families.
FAMILY_TYPE = 'family'
role_of = {r['ACC']: r['role'] for r in T}
class_of = {r['ACC']: r['localisation_class'] for r in T}
applicable = {r['ACC']: r['domain_annotation_applies_to'] == 'entry_matches_frozen_sequence' for r in T}

fam_rows = [r for r in PF if r['pfam_entry_type'] == FAMILY_TYPE]
fam_by_acc = collections.defaultdict(list)
for r in fam_rows: fam_by_acc[r['ACC']].append(r['pfam'])
pf_by_acc_all = collections.defaultdict(list)
for r in PF: pf_by_acc_all[r['ACC']].append(r)

def ids_in(roles, rows):
    return {r['pfam'] for r in rows if role_of[r['ACC']] in roles}

# all Pfam signatures (any type), for the ID-level absence question
all_ids, fit_ids = {r['pfam'] for r in PF}, ids_in({'discovery'}, PF)
tun_ids, eval_ids = ids_in({'tuning'}, PF), ids_in({'evaluation'}, PF)
fit_tun_ids = fit_ids | tun_ids
# family-type signatures only
fam_all = {r['pfam'] for r in fam_rows}
fam_fit = ids_in({'discovery'}, fam_rows)
fam_tun = ids_in({'tuning'}, fam_rows)
fam_eval = ids_in({'evaluation'}, fam_rows)
fam_fit_tun = fam_fit | fam_tun

R['pfam_id_coverage'] = {
 'note': 'ID-level absence over ALL Pfam signature types',
 'n_distinct_pfam_all_roles': len(all_ids),
 'n_distinct_in_discovery_P1P2': len(fit_ids),
 'n_distinct_in_tuning_P3': len(tun_ids),
 'n_distinct_in_evaluation_P4': len(eval_ids),
 'n_absent_from_P1P2': len(all_ids - fit_ids),
 'n_absent_from_P1P2P3': len(all_ids - fit_tun_ids),
 'n_eval_ids_absent_from_P1P2': len(eval_ids - fit_ids),
 'n_eval_ids_absent_from_P1P2P3': len(eval_ids - fit_tun_ids),
}
R['pfam_family_type_id_coverage'] = {
 'note': "restricted to Pfam signatures whose OWN entry type is 'family'",
 'n_distinct_family_signatures_all_roles': len(fam_all),
 'n_distinct_in_discovery_P1P2': len(fam_fit),
 'n_distinct_in_tuning_P3': len(fam_tun),
 'n_distinct_in_evaluation_P4': len(fam_eval),
 'n_absent_from_P1P2': len(fam_all - fam_fit),
 'n_absent_from_P1P2P3': len(fam_all - fam_fit_tun),
 'n_eval_family_ids_absent_from_P1P2': len(fam_eval - fam_fit),
 'n_eval_family_ids_absent_from_P1P2P3': len(fam_eval - fam_fit_tun),
}
R['pfam_entry_type_distribution_in_cohort'] = dict(
    collections.Counter(r['pfam_entry_type'] for r in PF))

with open(os.path.join(OUT, 'pfam_ids_absent.tsv'), 'w') as f:
    f.write('pfam\tpfam_entry_type\tis_family_type\tabsent_from_P1P2\tabsent_from_P1P2P3\t'
            'n_proteins_total\tn_in_tuning\tn_in_evaluation\tentry_name\n')
    cnt = collections.Counter(r['pfam'] for r in PF)
    ct = collections.Counter((r['pfam'], role_of[r['ACC']]) for r in PF)
    meta = {r['pfam']: (r['pfam_entry_type'], r['pfam_name_interpro'] or r['pfam_entry_name_uniprot'])
            for r in PF}
    for pid in sorted(all_ids - fit_ids):
        t, nm = meta[pid]
        f.write(f"{pid}\t{t}\t{int(t == FAMILY_TYPE)}\t1\t{int(pid not in fit_tun_ids)}\t"
                f"{cnt[pid]}\t{ct[(pid,'tuning')]}\t{ct[(pid,'evaluation')]}\t{nm}\n")

# ---- protein-level classification, FAMILY-TYPE signatures only, explicit unknown bucket ----
def classify_family(a, seen):
    if not applicable[a]:
        return 'unknown:applicability_uncertain'
    fams = fam_by_acc.get(a, [])
    if fams:
        hit = [x in seen for x in fams]
        return 'all_seen' if all(hit) else ('all_unseen' if not any(hit) else 'mixed_seen_unseen')
    allp = pf_by_acc_all.get(a, [])
    if not allp:
        return 'unknown:no_pfam_recorded'
    if any(x['pfam_entry_type'] == 'unavailable' for x in allp):
        return 'unknown:entry_type_unavailable'
    return 'unknown:pfam_recorded_but_no_family_type_signature'

ev = [r for r in T if r['role'] == 'evaluation']
for lbl, seen in (('vs_P1P2', fam_fit), ('vs_P1P2P3', fam_fit_tun)):
    c = collections.Counter(classify_family(r['ACC'], seen) for r in ev)
    byclass = collections.defaultdict(collections.Counter)
    for r in ev:
        byclass[r['localisation_class']][classify_family(r['ACC'], seen)] += 1
    R[f'evaluation_protein_family_status_{lbl}'] = {
     'basis': "Pfam signatures whose OWN entry type is 'family'",
     'total_evaluation_proteins': len(ev),
     'overall': {k: c[k] for k in sorted(c)},
     'by_localisation_class': {k: dict(sorted(v.items())) for k, v in sorted(byclass.items())},
     'ONLY_UNSEEN_family_ids_by_class': {k: v.get('all_unseen', 0) for k, v in sorted(byclass.items())},
    }
with open(os.path.join(OUT, 'evaluation_family_status.tsv'), 'w') as f:
    f.write('reference_set\tlocalisation_class\tstatus\tn\n')
    for lbl, seen in (('P1P2', fam_fit), ('P1P2P3', fam_fit_tun)):
        byclass = collections.defaultdict(collections.Counter)
        for r in ev: byclass[r['localisation_class']][classify_family(r['ACC'], seen)] += 1
        for cls in sorted(byclass):
            for st in sorted(byclass[cls]):
                f.write(f'{lbl}\t{cls}\t{st}\t{byclass[cls][st]}\n')

# ---- Pfam-bearing counts: BEFORE vs AFTER the core exclusions ----
dom_all = collections.defaultdict(set)
for r in csv.DictReader(open(os.path.join(REPO, 'annotation_readiness/out/annotation_domains.tsv')),
                        delimiter='\t'):
    if r['database'] == 'Pfam': dom_all[r['ACC']].add(r['xref_id'])
full_manifest = {r['ACC'] for r in csv.DictReader(
    open(os.path.join(REPO, 'discovery_analysis/v1/inclusion_manifest.tsv')), delimiter='\t')}
core_set = set(man)
R['pfam_bearing_protein_counts'] = {
 'before_core_exclusions_n_cohort': len(full_manifest),
 'before_core_exclusions_n_with_pfam': len(set(dom_all) & full_manifest),
 'after_core_exclusions_n_cohort': len(core_set),
 'after_core_exclusions_n_with_pfam': len(set(dom_all) & core_set),
 'difference': len(set(dom_all) & full_manifest) - len(set(dom_all) & core_set),
 'clarification': 'the previously quoted 9,454 was computed over the FULL 10,139-protein '
                  'manifest, i.e. BEFORE the 35 core exclusions (8 E1 + 27 E2)',
 'n_with_at_least_one_family_type_signature_core': len(fam_by_acc),
 'n_with_pfam_but_no_family_type_signature_core': len(
     {r['ACC'] for r in PF} - set(fam_by_acc)),
}

# ---------------------------------------------------------------- prospective strata from P1+P2
strata = {}
for c in CONTINUOUS:
    tr = np.array([float(r[c]) for r in T if r['role'] == 'discovery' and r[c] != ''])
    cuts = np.unique(np.quantile(tr, QUANTILES, method='linear'))
    cuts = cuts[(cuts > tr.min()) & (cuts < tr.max())]
    out = {'defined_on': 'discovery P1+P2 observed values only', 'n_training_observed': int(tr.size),
           'quantiles': QUANTILES, 'cut_points': [round(float(x), 6) for x in cuts],
           'interval_rule': 'lower inclusive, upper exclusive; outer intervals unbounded; tied edges collapsed',
           'n_strata': int(len(cuts) + 1), 'bins': []}
    for i in range(len(cuts) + 1):
        b = {'stratum': i + 1,
             'lower': None if i == 0 else round(float(cuts[i - 1]), 6),
             'upper': None if i == len(cuts) else round(float(cuts[i]), 6)}
        for role in ('discovery', 'tuning', 'evaluation'):
            g = [r for r in T if r['role'] == role and r[c] != '']
            v = np.array([float(r[c]) for r in g])
            a = np.searchsorted(cuts, v, side='right')
            m = a == i
            ex = sum(1 for r, k in zip(g, m) if k and r['target'] == '1')
            cy = sum(1 for r, k in zip(g, m) if k and r['target'] == '0')
            b[role] = {'n': int(m.sum()), 'n_Extracellular': ex, 'n_Cytoplasm': cy,
                       'both_classes_present': bool(ex > 0 and cy > 0),
                       'n_clusters_id30': len({r['cluster_id30'] for r, k in zip(g, m) if k})}
        out['bins'].append(b)
    out['all_strata_have_both_classes_in_tuning'] = all(b['tuning']['both_classes_present'] for b in out['bins'])
    out['all_strata_have_both_classes_in_evaluation'] = all(b['evaluation']['both_classes_present'] for b in out['bins'])
    out['min_minority_count_evaluation'] = min(min(b['evaluation']['n_Extracellular'],
                                                   b['evaluation']['n_Cytoplasm']) for b in out['bins'])
    strata[c] = out
R['prospective_strata'] = strata
with open(os.path.join(OUT, 'prospective_strata.tsv'), 'w') as f:
    f.write('property\tstratum\tlower\tupper\trole\tn\tn_Extracellular\tn_Cytoplasm\t'
            'both_classes\tn_clusters_id30\n')
    for c, o in strata.items():
        for b in o['bins']:
            for role in ('discovery', 'tuning', 'evaluation'):
                d = b[role]
                f.write(f"{c}\t{b['stratum']}\t{'' if b['lower'] is None else b['lower']}\t"
                        f"{'' if b['upper'] is None else b['upper']}\t{role}\t{d['n']}\t"
                        f"{d['n_Extracellular']}\t{d['n_Cytoplasm']}\t{int(d['both_classes_present'])}\t"
                        f"{d['n_clusters_id30']}\n")

# family strata support (observed-in-discovery families)
fam = []
_mem_index = collections.defaultdict(set)
for r in PF: _mem_index[r['pfam']].add(r['ACC'])
_famtype = {r['pfam']: r['pfam_entry_type'] for r in PF}
for p in sorted(fit_ids):
    mem = _mem_index[p]
    row = {'pfam': p, 'pfam_entry_type': _famtype.get(p, 'unavailable'),
           'is_family_type': int(_famtype.get(p, '') == FAMILY_TYPE)}
    for role in ('discovery', 'tuning', 'evaluation'):
        g = [r for r in T if r['role'] == role and r['ACC'] in mem]
        row[f'{role}_n'] = len(g)
        row[f'{role}_ex'] = sum(1 for r in g if r['target'] == '1')
        row[f'{role}_cy'] = sum(1 for r in g if r['target'] == '0')
    fam.append(row)
_ft = [r for r in fam if r['is_family_type'] == 1]
R['family_strata'] = {
 'basis_all_signature_types': True,
 'n_signatures_observed_in_discovery_all_types': len(fam),
 'n_signatures_observed_in_discovery_family_type_only': len(_ft),
 'family_type_only': {
   'n_with_any_evaluation_member': sum(1 for r in _ft if r['evaluation_n'] > 0),
   'n_with_ge2_evaluation_members': sum(1 for r in _ft if r['evaluation_n'] >= 2),
   'n_with_both_classes_in_evaluation': sum(1 for r in _ft if r['evaluation_ex'] > 0 and r['evaluation_cy'] > 0),
   'n_with_both_classes_in_tuning': sum(1 for r in _ft if r['tuning_ex'] > 0 and r['tuning_cy'] > 0)},
 'n_families_observed_in_discovery': len(fam),
 'n_with_both_classes_in_evaluation': sum(1 for r in fam if r['evaluation_ex'] > 0 and r['evaluation_cy'] > 0),
 'n_with_any_evaluation_member': sum(1 for r in fam if r['evaluation_n'] > 0),
 'n_with_ge2_evaluation_members': sum(1 for r in fam if r['evaluation_n'] >= 2),
 'n_with_both_classes_in_tuning': sum(1 for r in fam if r['tuning_ex'] > 0 and r['tuning_cy'] > 0),
}
with open(os.path.join(OUT, 'family_strata_support.tsv'), 'w') as f:
    cols = list(fam[0]); f.write('\t'.join(cols) + '\n')
    for r in fam: f.write('\t'.join(str(r[c]) for c in cols) + '\n')

R['signalp6'] = {'installed': False, 'predictions_generated': 0,
                 'all_proteins_missing': len(T),
                 'substituted_with_uniprot': False}
json.dump(R, open(os.path.join(OUT, 'readiness_assessment.json'), 'w'), indent=1)

print('=== integrity ==='); [print(f'  {k}: {v}') for k, v in ck.items()]
print('\n=== by role ===')
for role, d in by_role.items():
    print(f"  {role:<11} n={d['n']:>5} EX={d['n_Extracellular']:>4} CY={d['n_Cytoplasm']:>4} "
          f"clus30={d['n_clusters_id30']:>5} pfam_none={d['pfam_no_ids_recorded']:>4} "
          f"sig_present={d['signal_annotated_present']:>4} "
          f"kd_missing={d['hydrophobicity_kd_30_missing']} chg_missing={d['charge_balance_30_missing']}")
print('\n=== pfam ID coverage (all signature types) ===')
for k, v in R['pfam_id_coverage'].items(): print(f'  {k}: {v}')
print('\n=== pfam FAMILY-TYPE signature coverage ===')
for k, v in R['pfam_family_type_id_coverage'].items(): print(f'  {k}: {v}')
print('\n=== pfam-bearing counts before/after core exclusions ===')
for k, v in R['pfam_bearing_protein_counts'].items(): print(f'  {k}: {v}')
print('\n=== evaluation protein family status (family-type signatures only) ===')
for lbl in ('vs_P1P2', 'vs_P1P2P3'):
    d = R[f'evaluation_protein_family_status_{lbl}']
    print(f"  --- {lbl} (n={d['total_evaluation_proteins']}) ---")
    for k, v in d['overall'].items(): print(f'      {k:<52} {v}')
    print(f"      ONLY-UNSEEN by class: {d['ONLY_UNSEEN_family_ids_by_class']}")
print('\n=== strata ===')
for c, o in strata.items():
    print(f"  {c}: cuts={o['cut_points']} strata={o['n_strata']} "
          f"both_classes_tuning={o['all_strata_have_both_classes_in_tuning']} "
          f"both_classes_eval={o['all_strata_have_both_classes_in_evaluation']} "
          f"min_minority_eval={o['min_minority_count_evaluation']}")
print('\n=== family strata ==='); [print(f'  {k}: {v}') for k, v in R['family_strata'].items()]
