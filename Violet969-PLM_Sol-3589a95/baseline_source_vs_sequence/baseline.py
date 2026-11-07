#!/usr/bin/env python3
"""
UESolDS baseline: how much of the solubility label is explained by header-derived
"source proxies" alone, versus amino-acid composition, versus both.

Train = embedding_dataset/Train_dataset.fasta  (minus PDB_230301_28958)
Val   = embedding_dataset/validation_dataset.fasta
Test  = NOT TOUCHED (reserved for final evaluation).

Label convention (verified in the audit step):
  raw '1' = soluble, raw '0' = insoluble
  evidence: utils/general.py:10  SOLUBILITY = ['0','1']
            datasets/transforms.py:72-87  SolubilityToInt  -> index(); docstring
            "0 for membrane bound OR UNKNOWN or 1 for soluble"
"""
import os, re, sys, json, collections, warnings
import numpy as np

# numpy 2.0.2 built against Apple Accelerate emits spurious
# "divide by zero / overflow / invalid encountered in matmul" RuntimeWarnings even on
# clean, finite random input (reproduced independently on this machine). The feature
# matrices below are explicitly asserted finite, so these warnings are cosmetic.
warnings.filterwarnings('ignore', message='.*encountered in matmul.*',
                        category=RuntimeWarning)
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.metrics import roc_auc_score

BASE = sys.argv[1]
OUT  = sys.argv[2]
EXCLUDE_FROM_TRAIN = {'PDB_230301_28958'}   # the one train/test duplicate found in the audit
SEED = 42
C_FIXED = 1.0                                # fixed regularization for this first comparison
STD_AA = "ACDEFGHIKLMNPQRSTVWY"

# ---------------------------------------------------------------- source proxy rules
# These are the pilot's ID-class rules, applied unchanged. "Source proxy" = a group
# derived purely from the FASTA header identifier, standing in for the upstream
# database / consortium the record came from. It is NOT a biological annotation.
#
#   PDB_serial   ^PDB_230301_\d+$          authors' 2023-03-01 PDB pull, opaque serial
#   RefSeq       ^(NP|XP|YP|WP)_\d+        NCBI RefSeq accession
#   DNASU_clone  ^[A-Za-z]{2}CD\d+_pET     DNASU clone ID + expression vector suffix
#   TB_Rv        ^TB-Rv                    M. tuberculosis Rv locus tag
#   SG_target    ^(JCSG|MCSG|NESG|NYSGRC|SSGCID|SGPP|SECSG|EFI|BSGI|BSGC|MPP|CESG|SGC|APC)
#                                          structural-genomics consortium target ID
#   other        anything else
# Within SG_target the "center" is the substring before the first '-' (e.g. MCSG-APC84365
# -> MCSG). For every non-SG class the center is the class name itself, so that a single
# categorical column "source_center" is well defined for every record.
SG_CENTERS = r'JCSG|MCSG|NESG|NYSGRC|SSGCID|SGPP|SECSG|EFI|BSGI|BSGC|MPP|CESG|SGC|APC'

def id_class(i):
    if re.match(r'^PDB_230301_\d+$', i):        return 'PDB_serial'
    if re.match(r'^(NP|XP|YP|WP)_\d+', i):      return 'RefSeq'
    if re.match(r'^[A-Za-z]{2}CD\d+_pET', i):   return 'DNASU_clone'
    if i.startswith('TB-Rv'):                   return 'TB_Rv'
    if re.match(r'^(' + SG_CENTERS + r')', i):  return 'SG_target'
    return 'other'

def source_center(i, cls):
    if cls == 'SG_target':
        m = re.match(r'^(' + SG_CENTERS + r')', i)
        return m.group(1)
    return cls

def parse_fasta(path):
    """Multi-line-aware FASTA parser."""
    recs, hdr, seq = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            if hdr is not None: recs.append((hdr, ''.join(seq)))
            hdr, seq = line[1:], []
        else:
            seq.append(line)
    if hdr is not None: recs.append((hdr, ''.join(seq)))
    return recs

def load(path, drop=frozenset()):
    rows, dropped = [], 0
    for hdr, seq in parse_fasta(path):
        parts = hdr.split()
        # test-style headers carry an md5 first; ID is the token before the label field
        pid = parts[1] if len(parts) == 3 else parts[0]
        if pid in drop:
            dropped += 1
            continue
        cls = id_class(pid)
        rows.append({'id': pid, 'y': int(parts[-1].split('-')[-1]), 'seq': seq,
                     'len': len(seq), 'cls': cls, 'center': source_center(pid, cls)})
    return rows, dropped

# ---------------------------------------------------------------- features
def seq_features(rows):
    """20 amino-acid fractions + log(length).

    Non-standard residues (only 'B' occurs: 2 training records, 1 residue each;
    0 in validation) are counted in the denominator (full sequence length) but in
    none of the 20 numerators, so those two rows sum to slightly below 1.0. No
    record is dropped and no residue is imputed. With 2/70030 records affected at
    1 residue each, any other convention (e.g. splitting B as 0.5 D + 0.5 N) moves
    the feature matrix by <1e-5 in those rows only.
    """
    X = np.zeros((len(rows), len(STD_AA) + 1), dtype=np.float64)
    for r_i, r in enumerate(rows):
        n = len(r['seq'])
        c = collections.Counter(r['seq'])
        for a_i, a in enumerate(STD_AA):
            X[r_i, a_i] = c.get(a, 0) / n
        X[r_i, -1] = np.log(n)
    return X

SEQ_NAMES = [f'frac_{a}' for a in STD_AA] + ['log_len']

def main():
    tr, n_drop = load(os.path.join(BASE, 'embedding_dataset', 'Train_dataset.fasta'),
                      EXCLUDE_FROM_TRAIN)
    va, _      = load(os.path.join(BASE, 'embedding_dataset', 'validation_dataset.fasta'))
    print(f"train n={len(tr)} (excluded {n_drop}: {sorted(EXCLUDE_FROM_TRAIN)})   val n={len(va)}")

    # ---------- step 1: source proxy composition tables ----------
    def table(rows, title):
        print(f"\n### {title} — by ID class (source proxy)")
        g = collections.defaultdict(collections.Counter)
        for r in rows: g[r['cls']][r['y']] += 1
        print(f"{'id_class':<14}{'n':>7}{'insoluble(0)':>14}{'soluble(1)':>12}{'%soluble':>10}")
        for k in sorted(g, key=lambda k: -sum(g[k].values())):
            n = sum(g[k].values())
            print(f"{k:<14}{n:>7}{g[k][0]:>14}{g[k][1]:>12}{g[k][1]/n*100:>9.1f}%")
        n = len(rows); s = sum(r['y'] for r in rows)
        print(f"{'TOTAL':<14}{n:>7}{n-s:>14}{s:>12}{s/n*100:>9.1f}%")

        print(f"\n### {title} — SG_target by center prefix")
        g2 = collections.defaultdict(collections.Counter)
        for r in rows:
            if r['cls'] == 'SG_target': g2[r['center']][r['y']] += 1
        print(f"{'center':<14}{'n':>7}{'insoluble(0)':>14}{'soluble(1)':>12}{'%soluble':>10}")
        for k in sorted(g2, key=lambda k: -sum(g2[k].values())):
            n = sum(g2[k].values())
            print(f"{k:<14}{n:>7}{g2[k][0]:>14}{g2[k][1]:>12}{g2[k][1]/n*100:>9.1f}%")
        return g, g2

    table(tr, "TRAIN"); table(va, "VALIDATION")

    # ---------- step 2: three baselines ----------
    ytr = np.array([r['y'] for r in tr]); yva = np.array([r['y'] for r in va])

    # A. source only — one-hot of (id_class, source_center).
    # handle_unknown='ignore' => a validation category never seen in training becomes an
    # all-zero block, i.e. it falls back to the intercept rather than erroring.
    enc = OneHotEncoder(handle_unknown='ignore', sparse_output=False)
    Ctr = np.array([[r['cls'], r['center']] for r in tr], dtype=object)
    Cva = np.array([[r['cls'], r['center']] for r in va], dtype=object)
    enc.fit(Ctr)                                  # fit on TRAIN only
    Atr, Ava = enc.transform(Ctr), enc.transform(Cva)
    unseen = int((Ava.sum(axis=1) < 2).sum())
    print(f"\nsource one-hot dims={Atr.shape[1]}  validation rows with >=1 unseen category: {unseen}")

    # B. sequence only — standardized (fit on TRAIN only)
    Btr_raw, Bva_raw = seq_features(tr), seq_features(va)
    sc = StandardScaler().fit(Btr_raw)            # fit on TRAIN only
    Btr, Bva = sc.transform(Btr_raw), sc.transform(Bva_raw)

    # C. both
    Dtr, Dva = np.hstack([Atr, Btr]), np.hstack([Ava, Bva])

    for nm_, M in [('Atr',Atr),('Ava',Ava),('Btr',Btr),('Bva',Bva),('Dtr',Dtr),('Dva',Dva)]:
        assert np.isfinite(M).all(), f'non-finite values in {nm_}'
    print('feature matrices verified finite')

    models, scores = {}, {}
    for name, Xtr, Xva in [('A_source_only', Atr, Ava),
                           ('B_sequence_only', Btr, Bva),
                           ('C_source_plus_sequence', Dtr, Dva)]:
        m = LogisticRegression(penalty='l2', C=C_FIXED, solver='lbfgs',
                               max_iter=5000, random_state=SEED)
        m.fit(Xtr, ytr)
        models[name] = m
        scores[name] = m.predict_proba(Xva)[:, 1]
        print(f"fitted {name}: {Xtr.shape[1]} features, converged n_iter={m.n_iter_[0]}")

    # ---------- step 3: validation AUROC ----------
    def auroc(y, s):
        y = np.asarray(y)
        if len(np.unique(y)) < 2: return None
        return roc_auc_score(y, s)

    results = []
    def evaluate(mask, group, label):
        y = yva[mask]; n1 = int(y.sum()); n0 = int(len(y) - n1)
        row = {'group_type': group, 'group': label, 'n': int(len(y)),
               'n_insoluble_0': n0, 'n_soluble_1': n1}
        for name in models:
            a = auroc(y, scores[name][mask])
            row[name] = None if a is None else round(float(a), 4)
        row['auroc_defined'] = not (n0 == 0 or n1 == 0)
        results.append(row); return row

    evaluate(np.ones(len(va), bool), 'overall', 'ALL_validation')
    sg = np.array([r['cls'] == 'SG_target' for r in va])
    evaluate(sg, 'SG_target_only', 'SG_target')

    centers = collections.Counter()
    for r in va:
        if r['cls'] == 'SG_target': centers[r['center']] += 1
    per_center_included = []
    for c in sorted(centers):
        mask = np.array([r['cls'] == 'SG_target' and r['center'] == c for r in va])
        y = yva[mask]; n1 = int(y.sum()); n0 = int(len(y) - n1)
        if n1 >= 20 and n0 >= 20:
            evaluate(mask, 'SG_center', c); per_center_included.append(c)
        else:
            results.append({'group_type': 'SG_center_EXCLUDED', 'group': c,
                            'n': int(mask.sum()), 'n_insoluble_0': n0, 'n_soluble_1': n1,
                            'A_source_only': None, 'B_sequence_only': None,
                            'C_source_plus_sequence': None,
                            'auroc_defined': not (n0 == 0 or n1 == 0)})

    cols = ['group_type', 'group', 'n', 'n_insoluble_0', 'n_soluble_1',
            'A_source_only', 'B_sequence_only', 'C_source_plus_sequence', 'auroc_defined']
    os.makedirs(OUT, exist_ok=True)
    tsv = os.path.join(OUT, 'validation_auroc.tsv')
    with open(tsv, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in results:
            f.write('\t'.join('' if r.get(c) is None else
                              ('undefined_single_class' if c.startswith(('A_','B_','C_')) and r.get(c) is None
                               else str(r.get(c)))
                              for c in cols) + '\n')

    print("\n" + "=" * 104)
    print("VALIDATION AUROC  (centers shown when >=20 soluble AND >=20 insoluble)")
    print("=" * 104)
    print(f"{'group':<26}{'n':>6}{'insol':>7}{'sol':>6}{'A_source':>11}{'B_seq':>9}{'C_both':>9}")
    print("-" * 104)
    for r in results:
        if r['group_type'] == 'SG_center_EXCLUDED': continue
        def fmt(v):
            if v is None: return 'undefined' if not r['auroc_defined'] else 'n/a'
            return f'{v:.4f}'
        tag = '' if r['group_type'] != 'SG_center' else '  '
        print(f"{tag+r['group']:<26}{r['n']:>6}{r['n_insoluble_0']:>7}{r['n_soluble_1']:>6}"
              f"{fmt(r['A_source_only']):>11}{fmt(r['B_sequence_only']):>9}{fmt(r['C_source_plus_sequence']):>9}")
    print("-" * 104)
    print("\nSG centers NOT evaluated (fewer than 20 of either class in validation):")
    for r in results:
        if r['group_type'] == 'SG_center_EXCLUDED':
            note = '  <- single class, AUROC undefined' if not r['auroc_defined'] else ''
            print(f"   {r['group']:<12} n={r['n']:<5} insoluble={r['n_insoluble_0']:<5} soluble={r['n_soluble_1']:<5}{note}")

    # top source coefficients for interpretation
    nm = list(enc.get_feature_names_out(['id_class', 'center']))
    co = models['A_source_only'].coef_[0]
    order = np.argsort(co)
    print("\nModel A strongest source-proxy coefficients (log-odds of SOLUBLE):")
    for i in order[::-1][:6]: print(f"   {nm[i]:<28}{co[i]:+.3f}")
    print("   ...")
    for i in order[:5]: print(f"   {nm[i]:<28}{co[i]:+.3f}")

    json.dump({'results': results, 'C': C_FIXED, 'seed': SEED,
               'train_n': len(tr), 'val_n': len(va),
               'excluded_from_train': sorted(EXCLUDE_FROM_TRAIN),
               'centers_evaluated': per_center_included,
               'unseen_category_val_rows': unseen},
              open(os.path.join(OUT, 'results.json'), 'w'), indent=1)
    print(f"\nwrote {tsv}")
    print(f"wrote {os.path.join(OUT,'results.json')}")

if __name__ == '__main__':
    main()
