#!/usr/bin/env python3
"""Build the frozen FASTA input for the homology audit.

Audit population: ALL 10,139 P1-P4 development proteins, regardless of the
engineering (E1) / P0-base-collision (E2) / coordinate-transfer (E3) exclusion
flags. Those flags are CARRIED as columns, not applied: a protein excluded from
the scientific-analysis population is still a sequence that can be homologous to
an included one, so leaving it out would understate cross-partition relatedness.

Sequences are the ORIGINAL full-length DeepLoc sequences, not ESM-truncated.
Every sequence is verified byte-for-byte against the frozen sha256 in the
inclusion manifest before it is written. Nothing is substituted or dropped.

Non-standard residues (U/X/B) are written VERBATIM. MMseqs2 maps unknown
residue codes to X internally for scoring; we do not pre-edit the sequence.
Isoform accessions (-N suffix) are resolved by their own accession, exactly as
they appear in the DeepLoc CSV. The canonical sequence is never substituted.

Usage: build_input.py <repo_root> <out_dir>
"""
import csv, sys, os, json, hashlib, collections

REPO, OUT = sys.argv[1], sys.argv[2]
MANIFEST = os.path.join(REPO, 'discovery_analysis/v1/inclusion_manifest.tsv')
RAWCSV   = os.path.join(REPO, 'DeepLoc2/raw/Swissprot_Train_Validation_dataset.csv')
STD_AA   = set('ACDEFGHIKLMNPQRSTVWY')

def sha(s): return hashlib.sha256(s.encode()).hexdigest()
def filesha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()

man = list(csv.DictReader(open(MANIFEST), delimiter='\t'))
assert len(man) == 10139, len(man)
assert len({r['ACC'] for r in man}) == 10139, 'duplicate ACC in manifest'

raw = {r['ACC']: r for r in csv.DictReader(open(RAWCSV, newline='', encoding='utf-8'))}

recs, problems = [], []
for r in man:
    acc = r['ACC']
    if acc not in raw:
        problems.append((acc, 'absent_from_raw_csv')); continue
    seq = raw[acc]['Sequence'].strip()
    if sha(seq) != r['frozen_sha256']:
        problems.append((acc, 'sha256_mismatch_vs_frozen')); continue
    if len(seq) != int(r['frozen_len']):
        problems.append((acc, 'length_mismatch_vs_frozen')); continue
    recs.append((acc, seq, r))

assert not problems, problems[:20]
assert len(recs) == 10139

os.makedirs(OUT, exist_ok=True)
fa = os.path.join(OUT, 'development_10139.fasta')
with open(fa, 'w') as f:
    for acc, seq, _ in recs:          # manifest order; deterministic
        f.write(f'>{acc}\n{seq}\n')

nonstd = collections.Counter(c for _, s, _ in recs for c in s if c not in STD_AA)
n_with_nonstd = sum(1 for _, s, _ in recs if any(c not in STD_AA for c in s))
L = [len(s) for _, s, _ in recs]

report = {
    'audit_population': 'all 10,139 P1-P4 development proteins (exclusion flags carried, not applied)',
    'n_sequences': len(recs),
    'n_unique_accessions': len({a for a, _, _ in recs}),
    'n_unique_sequences': len({s for _, s, _ in recs}),
    'sequence_source': 'DeepLoc2/raw/Swissprot_Train_Validation_dataset.csv (original full length)',
    'truncated_to_1022': False,
    'all_verified_against_frozen_sha256': True,
    'n_sequences_dropped_or_substituted': 0,
    'length_min': min(L), 'length_median': sorted(L)[len(L) // 2], 'length_max': max(L),
    'n_proteins_with_nonstandard_residues': n_with_nonstd,
    'nonstandard_residue_counts': dict(sorted(nonstd.items())),
    'nonstandard_residue_handling': 'written verbatim; MMseqs2 maps unknown codes to X at scoring time; no pre-edit, no substitution, no drop',
    'n_isoform_accessions': sum(1 for _, _, r in recs if r['is_isoform'] == '1'),
    'isoform_handling': 'resolved by own accession; canonical never substituted',
    'n_engineering_protein_E1': sum(1 for _, _, r in recs if r['engineering_protein'] == '1'),
    'n_p0_base_collision_E2': sum(1 for _, _, r in recs if r['p0_base_collision'] != 'none'),
    'n_in_core': sum(1 for _, _, r in recs if r['in_core'] == '1'),
    'role_counts': dict(collections.Counter(r['role'] for _, _, r in recs)),
    'class_counts': dict(collections.Counter(r['localisation_class'] for _, _, r in recs)),
    'input_hashes': {
        'inclusion_manifest.tsv': filesha(MANIFEST),
        'Swissprot_Train_Validation_dataset.csv': filesha(RAWCSV),
        'development_10139.fasta': filesha(fa),
    },
}
json.dump(report, open(os.path.join(OUT, 'input_report.json'), 'w'), indent=1)
for k, v in report.items():
    if k != 'input_hashes': print(f'{k}: {v}')
print('\ninput hashes:')
for k, v in report['input_hashes'].items(): print(f'  {k}  {v}')
print(f'\nwrote {fa}')
