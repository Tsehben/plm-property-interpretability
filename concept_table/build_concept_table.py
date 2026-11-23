#!/usr/bin/env python3
"""Build the accession-keyed biological-characteristic (concept) table.

Population: the existing 10,104-protein eligible core. Partitions, roles, targets and
frozen sequence hashes are carried through UNCHANGED from the frozen manifest; nothing
is downsampled, refitted or re-extracted.

DERIVATIONS (frozen here, before any property-conditioned performance analysis):

  hydrophobicity_kd_30
      Mean Kyte-Doolittle score over residues 1..min(30, L).
      KD is defined for the 20 standard amino acids only. If the window contains ANY
      residue outside that set (U, X, B, Z, O), the score is recorded as MISSING with
      the offending residues named. The protein is KEPT. Nothing is imputed.

  charge_balance_30
      (count(K) + count(R) - count(D) - count(E)) / window_length over residues
      1..min(30, L). This is a CHARGE BALANCE, not a pH-dependent net charge:
        * histidine is OMITTED (pKa ~6.0, partially protonated at physiological pH,
          so its contribution is not a fixed integer);
        * the free alpha-amino and alpha-carboxyl TERMINAL charges are OMITTED;
        * no pKa model, no ionic-strength or pH term is applied.
      Recorded as MISSING only if the window contains X, B or Z, because those codes are
      ambiguous with respect to membership of {K,R,D,E} (X = any residue, B = D/N,
      Z = E/Q) and the count would be undefined. U (selenocysteine) and O (pyrrolysine)
      are SPECIFIC residues that are not K, R, D or E, so each contributes exactly zero
      and NEITHER makes this score missing. Both still make the KD score missing.

  Pfam
      ALL recorded Pfam accessions are retained per protein, in a long table, with name
      and InterPro entry type. The first ID is NEVER taken as a single family label.
      Entry-type metadata is fetched by ACCESSION from the InterPro entry catalogue; no
      protein sequence is transmitted.

  Homology
      Existing MMseqs2 cluster IDs at 30/50/80% identity, reused unchanged from
      homology_audit/v1. No clustering is recomputed.

  Signal
      Existing UniProt status, intervals, evidence categories and the flag recording
      whether coordinates apply to the frozen sequence. No prediction is made here.

Usage: build_concept_table.py <repo_root> <out_dir>
"""
import csv, json, os, sys, hashlib, collections

REPO, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
csv.field_size_limit(10_000_000)

KD = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5,
      'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
      'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
AMBIG_FOR_CHARGE = set('XBZ')           # ambiguous w.r.t. membership of {K,R,D,E}
# U (selenocysteine) and O (pyrrolysine) are SPECIFIC, unambiguous residues that are not
# K, R, D or E, so each contributes exactly zero to the numerator of the literal
# (K+R-D-E)/window_length proxy. Neither makes the charge score missing. Both still make
# the Kyte-Doolittle score missing, because KD defines no value for them.
WINDOW = 30

def fsha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()

P = {
 'manifest': os.path.join(REPO, 'discovery_analysis/v1/inclusion_manifest.tsv'),
 'clusters': os.path.join(REPO, 'homology_audit/v1/cluster_membership.tsv'),
 'pairs':    os.path.join(REPO, 'DeepLoc2/esm2_probe/deeploc_pairs.tsv'),
 'primary':  os.path.join(REPO, 'annotation_readiness/out/annotation_primary.tsv'),
 'domains':  os.path.join(REPO, 'annotation_readiness/out/annotation_domains.tsv'),
 'features': os.path.join(REPO, 'annotation_readiness/out/annotation_features.tsv'),
 'pfam_types': os.path.join(OUT, 'cache', 'pfam_entry_types.json'),
}

man = [r for r in csv.DictReader(open(P['manifest']), delimiter='\t') if r['in_core'] == '1']
assert len(man) == 10104, len(man)
acc_order = [r['ACC'] for r in man]
assert len(set(acc_order)) == 10104, 'duplicate accession in core'

clu = {r['ACC']: r for r in csv.DictReader(open(P['clusters']), delimiter='\t')}
seqs = {}
with open(P['pairs']) as f:
    next(f)
    for line in f:
        k, s = line.rstrip('\n').split('\t'); seqs[k] = s
prim = {r['ACC']: r for r in csv.DictReader(open(P['primary']), delimiter='\t')}
pfam_types = json.load(open(P['pfam_types'])) if os.path.exists(P['pfam_types']) else {}

dom = collections.defaultdict(list)
for r in csv.DictReader(open(P['domains']), delimiter='\t'):
    if r['database'] == 'Pfam': dom[r['ACC']].append(r)
feat = collections.defaultdict(list)
for r in csv.DictReader(open(P['features']), delimiter='\t'):
    if r['feature_type'] == 'Signal': feat[r['ACC']].append(r)

def hydro(seq):
    w = seq[:WINDOW]
    bad = sorted({c for c in w if c not in KD})
    if bad:
        return None, f'KD_undefined_for_residues:{",".join(bad)}'
    return sum(KD[c] for c in w) / len(w), 'ok'

def charge(seq):
    w = seq[:WINDOW]
    bad = sorted({c for c in w if c in AMBIG_FOR_CHARGE})
    if bad:
        return None, f'ambiguous_residues_for_KRDE_membership:{",".join(bad)}'
    n = w.count('K') + w.count('R') - w.count('D') - w.count('E')
    return n / len(w), 'ok'

rows, pf_long, sg_long = [], [], []
miss = collections.Counter()
for r in man:
    a = r['ACC']; s = seqs[a]
    assert hashlib.sha256(s.encode()).hexdigest() == r['frozen_sha256'], a
    h, hs = hydro(s); c, cs = charge(s)
    if hs != 'ok': miss['hydrophobicity'] += 1
    if cs != 'ok': miss['charge_balance'] += 1
    pr = prim.get(a, {})
    pfs = dom.get(a, [])
    sgs = feat.get(a, [])
    for d in pfs:
        t = pfam_types.get(d['xref_id'], {})
        pf_long.append({'ACC': a, 'frozen_sha256': r['frozen_sha256'],
                        'partition': r['partition'], 'role': r['role'],
                        'target': r['target'], 'localisation_class': r['localisation_class'],
                        'pfam': d['xref_id'], 'pfam_entry_name_uniprot': d['entry_name'],
                        'pfam_name_interpro': t.get('name', ''),
                        'pfam_entry_type': t.get('type', 'unavailable'),
                        'pfam_integrated_interpro': t.get('integrated_interpro') or '',
                        'domain_annotation_applies_to': r['domain_annotation_applies_to'],
                        'coordinates_transferable': r['coordinates_transferable'],
                        'sequence_match': r['sequence_match']})
    for d in sgs:
        sg_long.append({'ACC': a, 'frozen_sha256': r['frozen_sha256'],
                        'partition': r['partition'], 'role': r['role'],
                        'target': r['target'], 'start': d['start'],
                        'start_modifier': d['start_modifier'], 'end': d['end'],
                        'end_modifier': d['end_modifier'],
                        'evidence_code': d['evidence_code'],
                        'evidence_category': d['evidence_category'],
                        'coordinates_apply_to_frozen_sequence': d['coordinates_apply_to_frozen_sequence']})
    types = [pfam_types.get(d['xref_id'], {}).get('type', 'unavailable') for d in pfs]
    rows.append({
        'ACC': a, 'frozen_sha256': r['frozen_sha256'], 'partition': r['partition'],
        'role': r['role'], 'target': r['target'], 'localisation_class': r['localisation_class'],
        'kingdom_deeploc': r['kingdom_deeploc'], 'frozen_len': r['frozen_len'],
        'is_isoform': r['is_isoform'],
        'cluster_id30': clu[a]['cluster_id30'], 'cluster_size_id30': clu[a]['cluster_size_id30'],
        'cluster_id50': clu[a]['cluster_id50'], 'cluster_size_id50': clu[a]['cluster_size_id50'],
        'cluster_id80': clu[a]['cluster_id80'], 'cluster_size_id80': clu[a]['cluster_size_id80'],
        'window_length': min(WINDOW, len(s)),
        'hydrophobicity_kd_30': '' if h is None else f'{h:.6f}',
        'hydrophobicity_kd_30_status': hs,
        'charge_balance_30': '' if c is None else f'{c:.6f}',
        'charge_balance_30_status': cs,
        'n_pfam': len(pfs),
        'pfam_ids': ';'.join(d['xref_id'] for d in pfs),
        'pfam_entry_names': ';'.join(d['entry_name'] for d in pfs),
        'pfam_entry_types': ';'.join(types),
        'n_pfam_type_family': sum(1 for t in types if t == 'family'),
        'n_pfam_type_domain': sum(1 for t in types if t == 'domain'),
        'n_pfam_type_other': sum(1 for t in types if t not in ('family', 'domain', 'unavailable')),
        'n_pfam_type_unavailable': sum(1 for t in types if t == 'unavailable'),
        'pfam_annotation_status': ('no_pfam_recorded' if not pfs else
                                   'recorded' if r['domain_annotation_applies_to'] == 'entry_matches_frozen_sequence'
                                   else 'recorded_applicability_uncertain'),
        'domain_annotation_applies_to': r['domain_annotation_applies_to'],
        'coordinates_transferable': r['coordinates_transferable'],
        'sequence_match': r['sequence_match'],
        'n_interpro': r['n_interpro'],
        'signal_status': r['signal_status'],
        'signal_status_reason': pr.get('signal_status_reason', ''),
        'n_signal_features': r['n_signal_features'],
        'signal_evidence_best': r['signal_evidence_best'],
        'signal_evidence_categories': pr.get('signal_evidence_categories', ''),
        'parent_entry_has_signal': pr.get('parent_entry_has_signal', ''),
        'signal_intervals': ';'.join(f"{d['start']}-{d['end']}" for d in sgs),
        'signal_coords_apply_to_frozen': ';'.join(d['coordinates_apply_to_frozen_sequence'] for d in sgs),
        'transit_status': r['transit_status'],
        'signalp6_predicted_class': '', 'signalp6_score': '',
        'signalp6_cleavage_site': '', 'signalp6_status': 'not_run_no_local_installation',
    })

def write(name, data):
    p = os.path.join(OUT, name)
    cols = list(data[0])
    with open(p, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for d in data: f.write('\t'.join(str(d[c]) for c in cols) + '\n')
    return p, len(data)

p1, n1 = write('concept_table.tsv', rows)
p2, n2 = write('concept_pfam_long.tsv', pf_long)
p3, n3 = write('concept_signal_long.tsv', sg_long)

prov = {
 'population': 'existing eligible core, in_core==1', 'n': len(rows),
 'partitions_and_roles': 'carried through UNCHANGED from the frozen manifest',
 'downsampled': False, 'probes_refitted': False, 'esm_rerun': False, 'P0_accessed': False,
 'frozen_before_property_conditioned_analysis': True,
 'window': WINDOW,
 'hydrophobicity': {'scale': 'Kyte-Doolittle', 'window': 'residues 1..min(30,L)',
   'defined_for': 'the 20 standard amino acids only',
   'missing_rule': 'any residue outside the 20 standard AAs makes the score missing; '
                   'the protein is kept and the offending residues are named; nothing imputed'},
 'charge_balance': {'formula': '(count(K)+count(R)-count(D)-count(E))/window_length',
   'window': 'residues 1..min(30,L)',
   'name': 'charge balance - NOT a pH-dependent net charge',
   'omissions': ['histidine (pKa ~6.0, partially protonated at physiological pH)',
                 'free alpha-amino and alpha-carboxyl terminal charges',
                 'any pKa model, pH term or ionic-strength term'],
   'missing_rule': 'missing only if the window contains X, B or Z (ambiguous w.r.t. '
                   'membership of {K,R,D,E}); U (selenocysteine) and O (pyrrolysine) are '
                   'specific residues outside {K,R,D,E}, each contributing exactly zero, '
                   'and NEITHER triggers missing'},
 'pfam': {'all_ids_retained': True, 'single_family_label_taken': False,
   'entry_type_source': 'InterPro entry catalogue, queried by Pfam accession',
   'sequences_transmitted': False,
   'n_entries_in_type_cache': len(pfam_types)},
 'homology': {'source': 'homology_audit/v1/cluster_membership.tsv',
   'thresholds': [30, 50, 80], 'recomputed': False},
 'signal': {'source': 'existing UniProt annotations',
   'fields': ['status', 'intervals', 'evidence categories', 'applicability flag'],
   'predictions_made_here': False},
 'signalp6': {'run': False, 'reason': 'no local installation found',
   'substituted_with_uniprot': False},
 'input_sha256': {k: fsha(v) for k, v in P.items() if os.path.exists(v)},
 'output_sha256': {},
 'missing_counts': dict(miss),
 'row_counts': {'concept_table.tsv': n1, 'concept_pfam_long.tsv': n2,
                'concept_signal_long.tsv': n3},
}
for p in (p1, p2, p3): prov['output_sha256'][os.path.basename(p)] = fsha(p)
json.dump(prov, open(os.path.join(OUT, 'concept_table_provenance.json'), 'w'), indent=1)

print(f'concept_table.tsv       {n1:,} rows')
print(f'concept_pfam_long.tsv   {n2:,} rows')
print(f'concept_signal_long.tsv {n3:,} rows')
print('missing:', dict(miss))
print('pfam entry types cached:', len(pfam_types))
