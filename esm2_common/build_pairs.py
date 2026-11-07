#!/usr/bin/env python3
"""Build <key>\t<sequence> inputs from the FROZEN manifests. Reserved data is never read."""
import os, sys, csv, json, hashlib, collections
ROOT = sys.argv[1]
U = os.path.join(ROOT, 'Violet969-PLM_Sol-3589a95')
D = os.path.join(ROOT, 'DeepLoc2')

def parse_fasta(p):
    recs, h, s = [], None, []
    for line in open(p):
        line = line.strip()
        if line.startswith('>'):
            if h is not None: recs.append((h, ''.join(s)))
            h, s = line[1:], []
        else: s.append(line)
    if h is not None: recs.append((h, ''.join(s)))
    return recs

# ---------- UESolDS: keyed on record_uid = "<split>:<fasta_record_index>" ----------
man = list(csv.DictReader(open(os.path.join(U, 'dev_cohort', 'dev_cohort_manifest.tsv')),
                          delimiter='\t'))
seqs = {}
for split, f in [('train', 'Train_dataset.fasta'), ('validation', 'validation_dataset.fasta')]:
    for i, (h, s) in enumerate(parse_fasta(os.path.join(U, 'embedding_dataset', f))):
        seqs[f'{split}:{i}'] = s
out = os.path.join(U, 'esm2_probe', 'uesolds_pairs.tsv')
miss = 0
with open(out, 'w') as fh:
    fh.write('key\tsequence\n')
    for r in man:
        s = seqs.get(r['record_uid'])
        if s is None: miss += 1; continue
        assert hashlib.sha256(s.encode()).hexdigest() == r['sha256'], \
            f"hash mismatch for {r['record_uid']}"
        assert len(s) == int(r['seq_len'])
        fh.write(f"{r['record_uid']}\t{s}\n")
print(f'UESolDS: manifest={len(man):,} written={len(man)-miss:,} missing={miss}  -> {out}')

# ---------- DeepLoc: development P1-4 only; P0 and HPA stay reserved ----------
dman = list(csv.DictReader(open(os.path.join(D, 'baseline', 'deeploc_binary_manifest.tsv')),
                           delimiter='\t'))
dev = [r for r in dman if r['split_role'] == 'development_P1_4']
rows = {r['ACC']: r['Sequence'].strip() for r in
        csv.DictReader(open(os.path.join(D, 'raw', 'Swissprot_Train_Validation_dataset.csv'),
                            newline='', encoding='utf-8'))}
out2 = os.path.join(D, 'esm2_probe', 'deeploc_pairs.tsv')
miss2 = 0
with open(out2, 'w') as fh:
    fh.write('key\tsequence\n')
    for r in dev:
        s = rows.get(r['ACC'])
        if s is None: miss2 += 1; continue
        assert hashlib.sha256(s.encode()).hexdigest() == r['sha256'], \
            f"hash mismatch for {r['ACC']}"
        fh.write(f"{r['ACC']}\t{s}\n")
print(f'DeepLoc dev(P1-4): manifest={len(dev):,} written={len(dev)-miss2:,} '
      f'missing={miss2}  -> {out2}')
print(f'  RESERVED and not written: P0={sum(1 for r in dman if r["split_role"]=="held_out_P0"):,}'
      f', hpa_testset.csv, UESolDS test_dataset.fasta')
