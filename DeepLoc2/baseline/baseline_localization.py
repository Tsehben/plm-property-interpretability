#!/usr/bin/env python3
"""
DeepLoc 2.0 - first localization baseline: Extracellular vs Cytoplasm.

COHORT (frozen here)
  positive (target=1): Extracellular==1 AND Cytoplasm==0
  negative (target=0): Cytoplasm==1 AND Extracellular==0
  Rows positive for both or neither are excluded from the contrast entirely.
  All other localization annotations are retained as METADATA only.

EVALUATION PROTOCOL (frozen here)
  Partition 0            -> RESERVED for final evaluation. Never read by any model.
  Partitions 1,2,3,4     -> development. Four-fold CV: hold out each published
                            partition once, train on the other three.
  hpa_testset.csv        -> RESERVED. Its header has no 'Extracellular' column
                            (verified in the audit step), so it cannot score this
                            contrast as published.

MODELS (L2 logistic regression, fixed C=1, no class re-weighting, so the observed
class frequencies are preserved)
  A  kingdom only
  B  whole-sequence amino-acid fractions + log(length)
  C  B + amino-acid fractions over the first 30 residues
  D  C + kingdom

Encoders and scalers are fit INSIDE each training fold. Accession, partition and
every localization annotation column are excluded from the predictors.

Usage: python3 baseline_localization.py <DeepLoc2_dir> <out_dir>
"""
import os, sys, csv, json, hashlib, collections, warnings
import numpy as np

# numpy 2.0.2 built against Apple Accelerate emits spurious "divide by zero /
# overflow / invalid encountered in matmul" RuntimeWarnings even on clean finite
# input (reproduced on this machine with random data). Every feature matrix below
# is explicitly asserted finite before fitting, so these warnings are cosmetic.
warnings.filterwarnings('ignore', message='.*encountered in matmul.*',
                        category=RuntimeWarning)
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score

DEV_FILE = 'Swissprot_Train_Validation_dataset.csv'
POS, NEG = 'Extracellular', 'Cytoplasm'
HELD_OUT_PARTITION = '0'
CV_PARTITIONS = ['1', '2', '3', '4']
NTERM = 30                      # prespecified N-terminal window
C_FIXED, SEED = 1.0, 42
STD_AA = 'ACDEFGHIKLMNPQRSTVWY'
# All ten localization columns. 'Membrane' is a membrane-association flag, NOT a
# localization, and is EXCLUDED when computing has_additional_localization.
LOCS = ['Cytoplasm', 'Nucleus', 'Extracellular', 'Cell membrane', 'Mitochondrion',
        'Plastid', 'Endoplasmic reticulum', 'Lysosome/Vacuole', 'Golgi apparatus',
        'Peroxisome']
OTHER_LOCS = [c for c in LOCS if c not in (POS, NEG)]

def as_bin(v):
    s = (v or '').strip()
    if s == '': return None
    try: f = float(s)
    except ValueError: return None
    return 1 if f == 1.0 else (0 if f == 0.0 else None)

def aa_fracs(seq):
    """Fractions of the 20 standard amino acids.

    Non-standard residues present in this cohort are U (selenocysteine, 59
    occurrences), X (unknown, 20) and B (Asp/Asn ambiguity, 1). They are counted
    in the DENOMINATOR (the full window length) but in none of the 20 numerators,
    so affected rows sum to slightly below 1.0. No residue is imputed and no
    protein is dropped.
    """
    n = len(seq)
    if n == 0: return [0.0] * len(STD_AA)
    c = collections.Counter(seq)
    return [c.get(a, 0) / n for a in STD_AA]

def main():
    DL, OUT = sys.argv[1], sys.argv[2]
    os.makedirs(OUT, exist_ok=True)
    src = os.path.join(DL, 'raw', DEV_FILE)
    with open(src, newline='', encoding='utf-8') as f:
        rows = list(csv.DictReader(f))

    # ---------------- 1. freeze the cohort ----------------
    coh = []
    for r in rows:
        p, n = as_bin(r[POS]), as_bin(r[NEG])
        if (p, n) == (1, 0):   t = 1
        elif (p, n) == (0, 1): t = 0
        else:                  continue
        seq = r['Sequence'].strip()
        coh.append({
            'ACC': r['ACC'], 'target': t, 'kingdom': r['Kingdom'].strip(),
            'partition': r['Partition'].strip(), 'seq': seq, 'len': len(seq),
            'sha256': hashlib.sha256(seq.encode()).hexdigest(),
            # flag EXCLUDES the separate 'Membrane' column by construction
            'has_additional_localization': int(any(as_bin(r[c]) == 1 for c in OTHER_LOCS)),
            'other_localizations': ';'.join(c for c in OTHER_LOCS if as_bin(r[c]) == 1),
        })
    print(f'source rows={len(rows):,}  frozen cohort={len(coh):,}'
          f'  (positive={sum(c["target"] for c in coh):,})')

    held = [c for c in coh if c['partition'] == HELD_OUT_PARTITION]
    dev  = [c for c in coh if c['partition'] in CV_PARTITIONS]
    print(f'\ndevelopment (P1-4) = {len(dev):,}   held out (P0, reserved) = {len(held):,}')
    assert len(dev) == 10139, f'expected 10,139 development records, got {len(dev)}'
    assert len(held) == 2712, f'expected 2,712 held-out records, got {len(held)}'
    assert len(dev) + len(held) == len(coh)
    print('  expected totals verified: 10,139 development / 2,712 held out')
    print(f"  RESERVED: partition {HELD_OUT_PARTITION}; and hpa_testset.csv, whose header "
          f"has no '{POS}' column\n           so it cannot score this contrast as published.")

    # manifest over the whole frozen cohort
    mcols = ['ACC', 'sha256', 'target', 'kingdom', 'partition', 'len',
             'has_additional_localization', 'other_localizations', 'split_role']
    man = os.path.join(OUT, 'deeploc_binary_manifest.tsv')
    with open(man, 'w') as f:
        f.write('\t'.join(mcols) + '\n')
        for c in coh:
            c['split_role'] = ('held_out_P0' if c['partition'] == HELD_OUT_PARTITION
                               else ('development_P1_4' if c['partition'] in CV_PARTITIONS
                                     else 'unused'))
            f.write('\t'.join(str(c[k]) for k in mcols) + '\n')
    hs = [c['sha256'] for c in coh]
    print(f'  manifest sequence hashes: n={len(hs):,} unique={len(set(hs)):,} '
          f'duplicates={len(hs) - len(set(hs))}')

    # ---------------- 2. counts by kingdom x target on P1-4 ----------------
    print('\n' + '=' * 88)
    print('DEVELOPMENT (P1-4): counts by kingdom and target')
    print('=' * 88)
    kt = collections.defaultdict(collections.Counter)
    for c in dev: kt[c['kingdom']][c['target']] += 1
    print(f"{'kingdom':<18}{'n':>8}{'Cytoplasm(0)':>15}{'Extracellular(1)':>19}{'%positive':>11}")
    for k in sorted(kt, key=lambda k: -sum(kt[k].values())):
        n = sum(kt[k].values())
        print(f'{k:<18}{n:>8,}{kt[k][0]:>15,}{kt[k][1]:>19,}{kt[k][1] / n * 100:>10.1f}%')
    n1 = sum(c['target'] for c in dev)
    print(f"{'TOTAL':<18}{len(dev):>8,}{len(dev) - n1:>15,}{n1:>19,}"
          f"{n1 / len(dev) * 100:>10.1f}%")

    print(f"\nper-fold (held-out partition) class counts:")
    print(f"{'partition':<12}{'n':>8}{'Cytoplasm(0)':>15}{'Extracellular(1)':>19}{'prevalence':>12}")
    for p in CV_PARTITIONS:
        v = [c for c in dev if c['partition'] == p]
        a = sum(c['target'] for c in v)
        print(f'P{p:<11}{len(v):>8,}{len(v) - a:>15,}{a:>19,}{a / len(v):>12.4f}')

    # ---------------- feature builders ----------------
    def feats_seq(recs):
        return np.array([aa_fracs(r['seq']) + [np.log(r['len'])] for r in recs], float)
    def feats_nterm(recs):
        # prespecified first-30-residue window; proteins shorter than 30 use what they have
        return np.array([aa_fracs(r['seq'][:NTERM]) for r in recs], float)
    short = sum(1 for r in dev if r['len'] < NTERM)
    print(f'\nN-terminal window = first {NTERM} residues; '
          f'development proteins shorter than {NTERM}: {short}')

    MODELS = ['A_kingdom', 'B_seq', 'C_seq_nterm', 'D_seq_nterm_kingdom']

    # ---------------- 3. four-fold CV, out-of-fold predictions ----------------
    oof = {m: np.full(len(dev), np.nan) for m in MODELS}
    idx = {p: np.array([i for i, c in enumerate(dev) if c['partition'] == p])
           for p in CV_PARTITIONS}
    y = np.array([c['target'] for c in dev])
    per_fold = []

    for p in CV_PARTITIONS:
        te = idx[p]
        tr = np.array(sorted(set(range(len(dev))) - set(te.tolist())))
        Rtr = [dev[i] for i in tr]; Rte = [dev[i] for i in te]

        # fit encoders/scalers INSIDE the training fold only
        enc = OneHotEncoder(handle_unknown='ignore', sparse_output=False).fit(
            np.array([[r['kingdom']] for r in Rtr], dtype=object))
        Ktr = enc.transform(np.array([[r['kingdom']] for r in Rtr], dtype=object))
        Kte = enc.transform(np.array([[r['kingdom']] for r in Rte], dtype=object))

        sS = StandardScaler().fit(feats_seq(Rtr))
        Str, Ste = sS.transform(feats_seq(Rtr)), sS.transform(feats_seq(Rte))
        sN = StandardScaler().fit(feats_nterm(Rtr))
        Ntr, Nte = sN.transform(feats_nterm(Rtr)), sN.transform(feats_nterm(Rte))

        blocks = {
            'A_kingdom':           (Ktr, Kte),
            'B_seq':               (Str, Ste),
            'C_seq_nterm':         (np.hstack([Str, Ntr]), np.hstack([Ste, Nte])),
            'D_seq_nterm_kingdom': (np.hstack([Str, Ntr, Ktr]), np.hstack([Ste, Nte, Kte])),
        }
        row = {'fold_heldout_partition': p, 'n_train': len(tr), 'n_test': len(te),
               'n_pos': int(y[te].sum()), 'n_neg': int(len(te) - y[te].sum()),
               'prevalence': round(float(y[te].mean()), 4)}
        for m in MODELS:
            Xtr, Xte = blocks[m]
            assert np.isfinite(Xtr).all() and np.isfinite(Xte).all(), f'non-finite in {m}'
            clf = LogisticRegression(penalty='l2', C=C_FIXED, solver='lbfgs',
                                     max_iter=5000, random_state=SEED).fit(Xtr, y[tr])
            s = clf.predict_proba(Xte)[:, 1]
            oof[m][te] = s
            row[m + '__auroc'] = round(float(roc_auc_score(y[te], s)), 4)
            row[m + '__ap']    = round(float(average_precision_score(y[te], s)), 4)
        per_fold.append(row)
        print(f'  fold P{p}: train={len(tr):,} test={len(te):,} '
              f'pos={row["n_pos"]:,} prev={row["prevalence"]:.4f}')

    for m in MODELS: assert not np.isnan(oof[m]).any(), f'missing OOF predictions for {m}'

    pred = os.path.join(OUT, 'deeploc_oof_predictions.tsv')
    pcols = ['ACC', 'partition', 'kingdom', 'target', 'has_additional_localization'] + MODELS
    with open(pred, 'w') as f:
        f.write('\t'.join(pcols) + '\n')
        for i, c in enumerate(dev):
            f.write('\t'.join([c['ACC'], c['partition'], c['kingdom'], str(c['target']),
                               str(c['has_additional_localization'])]
                              + [f'{oof[m][i]:.6f}' for m in MODELS]) + '\n')

    # ---------------- reporting ----------------
    def block(title, mask, results):
        yy = y[mask]; n1 = int(yy.sum()); n0 = int(len(yy) - n1)
        print('\n' + '=' * 110); print(title); print('=' * 110)
        print(f'  n={len(yy):,}   Cytoplasm(0)={n0:,}   Extracellular(1)={n1:,}   '
              f'positive prevalence={yy.mean():.4f}')
        print(f"\n{'model':<24}" + ''.join(f'{"P"+p:>9}' for p in CV_PARTITIONS)
              + f'{"MEAN":>10}{"|":>3}' + ''.join(f'{"P"+p:>9}' for p in CV_PARTITIONS)
              + f'{"MEAN":>10}')
        print(f"{'':<24}{'--- AUROC per held-out fold ---':^39}{'|':>3}"
              f"{'--- Average precision ---':^39}")
        print('-' * 110)
        for m in MODELS:
            au, ap = [], []
            for p in CV_PARTITIONS:
                fm = mask & np.isin(np.array([c['partition'] for c in dev]), [p])
                au.append(roc_auc_score(y[fm], oof[m][fm]))
                ap.append(average_precision_score(y[fm], oof[m][fm]))
            print(f'{m:<24}' + ''.join(f'{v:>9.4f}' for v in au) + f'{np.mean(au):>10.4f}{"|":>3}'
                  + ''.join(f'{v:>9.4f}' for v in ap) + f'{np.mean(ap):>10.4f}')
            results.append({'subset': title, 'model': m, 'n': int(len(yy)),
                            'n_neg_cytoplasm': n0, 'n_pos_extracellular': n1,
                            'prevalence': round(float(yy.mean()), 4),
                            **{f'auroc_P{p}': round(float(v), 4)
                               for p, v in zip(CV_PARTITIONS, au)},
                            'auroc_mean': round(float(np.mean(au)), 4),
                            **{f'ap_P{p}': round(float(v), 4)
                               for p, v in zip(CV_PARTITIONS, ap)},
                            'ap_mean': round(float(np.mean(ap)), 4)})
        # paired per-fold B -> C change
        dau, dap = [], []
        for p in CV_PARTITIONS:
            fm = mask & np.isin(np.array([c['partition'] for c in dev]), [p])
            dau.append(roc_auc_score(y[fm], oof['C_seq_nterm'][fm])
                       - roc_auc_score(y[fm], oof['B_seq'][fm]))
            dap.append(average_precision_score(y[fm], oof['C_seq_nterm'][fm])
                       - average_precision_score(y[fm], oof['B_seq'][fm]))
        print('-' * 110)
        print(f'{"delta C - B (paired)":<24}' + ''.join(f'{v:>+9.4f}' for v in dau)
              + f'{np.mean(dau):>+10.4f}{"|":>3}' + ''.join(f'{v:>+9.4f}' for v in dap)
              + f'{np.mean(dap):>+10.4f}')
        print(f'   folds improved by adding the N-terminal window: '
              f'AUROC {sum(1 for v in dau if v > 0)}/4, AP {sum(1 for v in dap if v > 0)}/4')
        results.append({'subset': title, 'model': 'delta_C_minus_B', 'n': int(len(yy)),
                        'n_neg_cytoplasm': n0, 'n_pos_extracellular': n1,
                        'prevalence': round(float(yy.mean()), 4),
                        **{f'auroc_P{p}': round(float(v), 4) for p, v in zip(CV_PARTITIONS, dau)},
                        'auroc_mean': round(float(np.mean(dau)), 4),
                        **{f'ap_P{p}': round(float(v), 4) for p, v in zip(CV_PARTITIONS, dap)},
                        'ap_mean': round(float(np.mean(dap)), 4)})
        return results

    results = []
    allm = np.ones(len(dev), bool)
    results = block('ALL development proteins (P1-4)', allm, results)
    single = np.array([c['has_additional_localization'] == 0 for c in dev])
    results = block('SENSITIVITY: proteins with NO additional localization annotation',
                    single, results)

    rcols = (['subset', 'model', 'n', 'n_neg_cytoplasm', 'n_pos_extracellular', 'prevalence']
             + [f'auroc_P{p}' for p in CV_PARTITIONS] + ['auroc_mean']
             + [f'ap_P{p}' for p in CV_PARTITIONS] + ['ap_mean'])
    res = os.path.join(OUT, 'deeploc_baseline_results.tsv')
    with open(res, 'w') as f:
        f.write('\t'.join(rcols) + '\n')
        for r in results: f.write('\t'.join(str(r[c]) for c in rcols) + '\n')

    json.dump({'results': results, 'per_fold_setup': per_fold, 'C': C_FIXED, 'seed': SEED,
               'n_cohort': len(coh), 'n_dev': len(dev), 'n_heldout_P0': len(held),
               'nterm_window': NTERM, 'models': MODELS,
               'kingdom_x_target': {k: dict(v) for k, v in kt.items()},
               'reserved': {'partition_0': len(held), 'hpa_testset':
                            'header lacks Extracellular column'}},
              open(os.path.join(OUT, 'deeploc_baseline_results.json'), 'w'), indent=1)
    print(f'\nwrote {man}\nwrote {pred}\nwrote {res}')
    print(f"wrote {os.path.join(OUT, 'deeploc_baseline_results.json')}")

if __name__ == '__main__':
    main()
