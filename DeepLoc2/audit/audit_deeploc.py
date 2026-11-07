#!/usr/bin/env python3
"""
DeepLoc 2.0 data audit (read-only).

Inputs (originals, untouched):
  DeepLoc2/raw/Swissprot_Train_Validation_dataset.csv   development set, published Partition 0-4
  DeepLoc2/raw/hpa_testset.csv                          HPA test set - HEADER ONLY here

Purpose: size the candidate binary contrast Extracellular vs Cytoplasm while
preserving the published Partition assignments. Test evaluation is reserved.

Usage: python3 audit_deeploc.py <DeepLoc2_dir> <out_dir>
"""
import os, sys, csv, json, hashlib, collections

DEV_FILE = 'Swissprot_Train_Validation_dataset.csv'
HPA_FILE = 'hpa_testset.csv'
A, Bc = 'Extracellular', 'Cytoplasm'          # the two target labels for the contrast
# The 10 localization columns in the DeepLoc 2.0 development set. 'Membrane' is a
# separate membrane-association flag, NOT a localization, so it is excluded here.
LOCS = ['Cytoplasm', 'Nucleus', 'Extracellular', 'Cell membrane', 'Mitochondrion',
        'Plastid', 'Endoplasmic reticulum', 'Lysosome/Vacuole', 'Golgi apparatus',
        'Peroxisome']

def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 20), b''): h.update(c)
    return h.hexdigest()

def as_binary(v):
    """Return 1, 0, or None if the value is missing/unrecognized."""
    s = (v or '').strip()
    if s == '': return None
    try: fv = float(s)
    except ValueError: return None
    if fv == 1.0: return 1
    if fv == 0.0: return 0
    return None

def main():
    DL, OUT = sys.argv[1], sys.argv[2]
    os.makedirs(OUT, exist_ok=True)
    dev = os.path.join(DL, 'raw', DEV_FILE)
    hpa = os.path.join(DL, 'raw', HPA_FILE)

    print('=' * 92)
    print('FILE PROVENANCE')
    print('=' * 92)
    for p, url in [(dev, 'https://services.healthtech.dtu.dk/services/DeepLoc-2.0/data/' + DEV_FILE),
                   (hpa, 'https://services.healthtech.dtu.dk/services/DeepLoc-2.0/data/' + HPA_FILE)]:
        print(f'  {os.path.basename(p)}\n    bytes  {os.path.getsize(p):,}\n'
              f'    sha256 {sha256(p)}\n    source {url}')

    # ---------------- development set ----------------
    with open(dev, newline='', encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    cols = list(rows[0].keys())
    print('\n' + '=' * 92); print('DEVELOPMENT SET - COLUMNS'); print('=' * 92)
    print(f'  rows={len(rows):,}  columns={len(cols)}')
    for c in cols: print(f'    - {c!r}')
    for t in (A, Bc, 'Partition'):
        assert t in cols, f'expected column {t!r} not found'

    print('\n  observed value sets (label columns):')
    for c in LOCS + ['Membrane']:
        if c in cols:
            vs = sorted({(r[c] or '').strip() for r in rows})
            print(f'    {c:<24}{vs if len(vs) <= 6 else str(vs[:6]) + " ..."}')
    print(f"    {'Partition':<24}{sorted({(r['Partition'] or '').strip() for r in rows})}")
    print(f"    {'Kingdom':<24}{sorted({(r['Kingdom'] or '').strip() for r in rows})}")

    # ---------------- partition x group table ----------------
    GROUPS = ['Extracellular_pos_Cytoplasm_neg', 'Cytoplasm_pos_Extracellular_neg',
              'both_positive', 'both_negative', 'missing_or_unrecognized']
    tab = collections.defaultdict(collections.Counter)
    other_loc = collections.defaultdict(collections.Counter)   # group -> partition -> count
    part_total = collections.Counter()
    bad_rows = []
    for r in rows:
        p = (r['Partition'] or '').strip()
        part_total[p] += 1
        a, b = as_binary(r[A]), as_binary(r[Bc])
        if a is None or b is None:
            g = 'missing_or_unrecognized'; bad_rows.append({'ACC': r.get('ACC'), 'Partition': p,
                                                            A: r[A], Bc: r[Bc]})
        elif a == 1 and b == 0: g = 'Extracellular_pos_Cytoplasm_neg'
        elif b == 1 and a == 0: g = 'Cytoplasm_pos_Extracellular_neg'
        elif a == 1 and b == 1: g = 'both_positive'
        else:                   g = 'both_negative'
        tab[g][p] += 1
        # does this protein carry ANY other positive localization label?
        if g in GROUPS[:2]:
            if any(as_binary(r[c]) == 1 for c in LOCS if c in cols and c not in (A, Bc)):
                other_loc[g][p] += 1

    parts = sorted(part_total)
    print('\n' + '=' * 92)
    print('PARTITION x TARGET-LABEL GROUP  (published Partition assignments preserved)')
    print('=' * 92)
    hdr = f"{'group':<36}" + ''.join(f'{"P"+p:>10}' for p in parts) + f'{"TOTAL":>10}'
    print(hdr); print('-' * len(hdr))
    for g in GROUPS:
        print(f'{g:<36}' + ''.join(f'{tab[g][p]:>10,}' for p in parts)
              + f'{sum(tab[g].values()):>10,}')
    print('-' * len(hdr))
    print(f"{'TOTAL rows':<36}" + ''.join(f'{part_total[p]:>10,}' for p in parts)
          + f'{len(rows):>10,}')

    print('\n  sum-check: groups must sum to the row count in every partition')
    ok = True
    for p in parts:
        s = sum(tab[g][p] for g in GROUPS)
        good = (s == part_total[p]); ok &= good
        print(f'    P{p}: groups={s:,}  rows={part_total[p]:,}  {"OK" if good else "MISMATCH"}')
    s = sum(sum(tab[g].values()) for g in GROUPS)
    ok &= (s == len(rows))
    print(f'    ALL: groups={s:,}  rows={len(rows):,}  {"OK" if s == len(rows) else "MISMATCH"}')
    assert ok, 'group counts do not reconcile with row counts'

    print('\n' + '=' * 92)
    print('MULTI-LOCALIZATION WITHIN THE TWO CANDIDATE GROUPS')
    print('(proteins ALSO positive for >=1 other localization; KEPT in their group)')
    print('=' * 92)
    hdr2 = f"{'group':<36}" + ''.join(f'{"P"+p:>10}' for p in parts) + f'{"TOTAL":>10}{"share":>9}'
    print(hdr2); print('-' * len(hdr2))
    for g in GROUPS[:2]:
        tot = sum(tab[g].values()); oth = sum(other_loc[g].values())
        print(f'{g:<36}' + ''.join(f'{other_loc[g][p]:>10,}' for p in parts)
              + f'{oth:>10,}{oth / tot * 100:>8.1f}%')

    cls = {g: sum(tab[g].values()) for g in GROUPS[:2]}
    print('\n  CANDIDATE CLASS SIZES for the binary contrast:')
    for g, n in cls.items(): print(f'    {g:<36}{n:>8,}')
    print(f'    {"combined":<36}{sum(cls.values()):>8,}'
          f'   ratio {cls[GROUPS[0]] / cls[GROUPS[1]]:.3f} : 1')

    # ---------------- duplicates & cross-partition overlap ----------------
    print('\n' + '=' * 92); print('DUPLICATE ACCESSIONS AND CROSS-PARTITION SEQUENCE OVERLAP'); print('=' * 92)
    accs = collections.Counter(r['ACC'] for r in rows)
    dup_acc = {k: v for k, v in accs.items() if v > 1}
    print(f'  rows={len(rows):,}  unique ACC={len(accs):,}  '
          f'ACCs appearing >1x={len(dup_acc):,}  extra rows={sum(v - 1 for v in dup_acc.values()):,}')
    if dup_acc:
        for k, v in sorted(dup_acc.items(), key=lambda kv: -kv[1])[:5]:
            ps = [r['Partition'] for r in rows if r['ACC'] == k]
            print(f'     {k} x{v} partitions={ps}')

    seqh = collections.defaultdict(list)
    for r in rows: seqh[hashlib.sha256(r['Sequence'].strip().encode()).hexdigest()].append(r)
    dup_seq = {k: v for k, v in seqh.items() if len(v) > 1}
    print(f'  unique sequences={len(seqh):,}  sequences appearing >1x={len(dup_seq):,}  '
          f'extra rows={sum(len(v) - 1 for v in dup_seq.values()):,}')
    cross = {k: v for k, v in dup_seq.items() if len({r['Partition'] for r in v}) > 1}
    print(f'  identical sequences spanning >1 PARTITION = {len(cross):,} '
          f'(rows involved: {sum(len(v) for v in cross.values()):,})')
    pair = collections.Counter()
    for v in cross.values():
        ps = sorted({r['Partition'] for r in v})
        for i in range(len(ps)):
            for j in range(i + 1, len(ps)): pair[(ps[i], ps[j])] += 1
    for (x, y), n in sorted(pair.items()): print(f'     P{x} <-> P{y}: {n}')
    conf = sum(1 for v in cross.values()
               if len({(as_binary(r[A]), as_binary(r[Bc])) for r in v}) > 1)
    print(f'  cross-partition duplicate sequences with CONFLICTING (Extracellular,Cytoplasm) = {conf}')

    # ---------------- HPA test set: header only ----------------
    print('\n' + '=' * 92); print('HPA TEST SET - HEADER INSPECTION ONLY (no rows read)'); print('=' * 92)
    with open(hpa, newline='', encoding='utf-8') as f:
        hcols = next(csv.reader(f))
    print(f'  columns ({len(hcols)}):')
    for c in hcols: print(f'    - {c!r}')
    avail = {t: (t in hcols) for t in (A, Bc)}
    print('\n  target-label availability:')
    for t, v in avail.items(): print(f'    {t:<16}{"PRESENT" if v else "*** ABSENT ***"}')

    # ---------------- outputs ----------------
    tsv = os.path.join(OUT, 'deeploc_partition_counts.tsv')
    with open(tsv, 'w') as f:
        f.write('group\t' + '\t'.join('P' + p for p in parts) + '\tTOTAL\n')
        for g in GROUPS:
            f.write(g + '\t' + '\t'.join(str(tab[g][p]) for p in parts)
                    + f'\t{sum(tab[g].values())}\n')
        f.write('TOTAL_rows\t' + '\t'.join(str(part_total[p]) for p in parts) + f'\t{len(rows)}\n')
        for g in GROUPS[:2]:
            f.write(g + '__also_other_localization\t'
                    + '\t'.join(str(other_loc[g][p]) for p in parts)
                    + f'\t{sum(other_loc[g].values())}\n')
    json.dump({'dev_file': DEV_FILE, 'dev_sha256': sha256(dev), 'dev_bytes': os.path.getsize(dev),
               'hpa_file': HPA_FILE, 'hpa_sha256': sha256(hpa), 'hpa_bytes': os.path.getsize(hpa),
               'n_rows': len(rows), 'columns': cols, 'partitions': parts,
               'partition_totals': dict(part_total),
               'groups': {g: {'by_partition': dict(tab[g]), 'total': sum(tab[g].values())}
                          for g in GROUPS},
               'also_other_localization': {g: {'by_partition': dict(other_loc[g]),
                                               'total': sum(other_loc[g].values())}
                                           for g in GROUPS[:2]},
               'candidate_class_sizes': cls,
               'duplicate_accessions': len(dup_acc),
               'duplicate_sequences': len(dup_seq),
               'cross_partition_duplicate_sequences': len(cross),
               'cross_partition_label_conflicts': conf,
               'hpa_columns': hcols, 'hpa_target_availability': avail},
              open(os.path.join(OUT, 'deeploc_audit_results.json'), 'w'), indent=1)
    print(f'\nwrote {tsv}')
    print(f"wrote {os.path.join(OUT, 'deeploc_audit_results.json')}")

if __name__ == '__main__':
    main()
