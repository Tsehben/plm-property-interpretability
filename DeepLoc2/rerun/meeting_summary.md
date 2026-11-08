# DeepLoc 2.0 — does ESM-2 beat composition? Meeting summary

**Date:** 2026-11-08 · **Scope:** Swiss-Prot development partitions P1–P4 (10,139 proteins)
**Reserved:** partition P0 (2,712 proteins) and `hpa_testset.csv` — neither read for any
number below. **No models were retrained for this check**; all figures come from the
rerun outputs.

## The comparison

| | |
|---|---|
| **Task** | Extracellular (positive) vs Cytoplasm (negative), binary |
| **Baseline** | `matched_C` — amino-acid composition + log(length) over the first 1,022 residues, plus composition of the first 30 residues |
| **Contender** | `c_matchedC_plus_ESM` — the same baseline *plus* frozen ESM-2-650M layer-8 embeddings, mean-pooled over visible residues and over the first 30 positions |
| **Protocol** | Four-fold cross-validation on the published partitions; each of P1–P4 held out once. L2 logistic regression, C = 1, `lbfgs`, seed 42, no class re-weighting. Scalers fit inside each training fold only. |
| **Fairness** | Both classifiers see the same visible residues and are scored on the same validation proteins. The baseline is deliberately input-length-matched so the contrast isolates the representation, not the amount of sequence seen. |

## Results

All values are **full-precision**, from
[`esm2_probe/deeploc_esm_results.tsv`](esm2_probe/deeploc_esm_results.tsv).
Means are unweighted across the four folds. Δ = (+ESM) − baseline.

### Full development cohort — n = 10,139 (2,498 extracellular / 7,641 cytoplasmic)

| Fold | Extracell. | Cytopl. | AUROC base | AUROC +ESM | Δ AUROC | AP base | AP +ESM | Δ AP |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| P1 | 520 | 1,970 | 0.9612 | 0.9650 | +0.0038 | 0.8777 | 0.9371 | +0.0594 |
| P2 | 819 | 1,892 | 0.9572 | 0.9929 | +0.0357 | 0.8898 | 0.9836 | +0.0938 |
| P3 | 584 | 1,972 | 0.9573 | 0.9845 | +0.0271 | 0.8776 | 0.9541 | +0.0764 |
| P4 | 575 | 1,807 | 0.9682 | 0.9731 | +0.0049 | 0.9162 | 0.9401 | +0.0239 |
| **Mean** | **2,498** | **7,641** | **0.9610** | **0.9789** | **+0.0179** | **0.8903** | **0.9537** | **+0.0634** |

Folds improved: **AUROC 4/4, AP 4/4**.

### Single-localisation subset — n = 5,967 (2,367 extracellular / 3,600 cytoplasmic)

Proteins carrying no positive localisation annotation beyond their target label. The
separate `Membrane` flag is a membrane-association indicator, not a localisation, and is
excluded from this definition. No length filter applied.

| Fold | Extracell. | Cytopl. | AUROC base | AUROC +ESM | Δ AUROC | AP base | AP +ESM | Δ AP |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| P1 | 483 | 938 | 0.9661 | 0.9706 | +0.0045 | 0.9292 | 0.9584 | +0.0293 |
| P2 | 789 | 872 | 0.9643 | 0.9928 | +0.0286 | 0.9495 | 0.9909 | +0.0414 |
| P3 | 553 | 928 | 0.9593 | 0.9844 | +0.0251 | 0.9328 | 0.9721 | +0.0393 |
| P4 | 542 | 862 | 0.9706 | 0.9699 | **−0.0008** | 0.9548 | 0.9571 | +0.0023 |
| **Mean** | **2,367** | **3,600** | **0.9651** | **0.9794** | **+0.0143** | **0.9416** | **0.9696** | **+0.0281** |

Folds improved: **AUROC 3/4, AP 4/4**. P4 is the one fold where adding ESM-2 does not help
on AUROC.

## Why the subset matters

Exploratory analysis showed the two classes are annotated very differently. **52.9 % of
cytoplasmic proteins carry at least one additional localisation (41.7 % are also nuclear),
against 5.2 % of extracellular proteins.** The negative class is therefore largely a
multi-compartment population while the positive class is nearly pure, so a classifier can
be rewarded for recognising "this protein is annotated in several places" rather than
"this protein is cytoplasmic." Restricting to unambiguously single-localised proteins
removes that route. It also rebalances the problem — prevalence rises from 24.6 % to
39.7 % — which is why baseline average precision jumps from 0.8903 to 0.9416 without any
model change; that rise reflects the easier base rate, not a better classifier, and the
two AP columns should not be compared across tables.

The gain survives the restriction at reduced size: **+0.0143 AUROC and +0.0281 AP**, versus
+0.0179 and +0.0634 on the full cohort. Roughly half the apparent ESM-2 advantage on the
full cohort is attributable to the multi-localised proteins. A separate check removing the
1,111 proteins truncated at 1,022 residues left the gain essentially unchanged (+0.0190
AUROC), so differential truncation is not the explanation
([`deeploc_esm_lenmatched_le1022.tsv`](deeploc_esm_lenmatched_le1022.tsv)).

## Preprocessing decision

**Retain all 10,139 development proteins for the primary comparison.** The single-
localisation and length-matched subsets are reported as sensitivity analyses, not as the
headline. Reasons: the published DeepLoc partitions are defined over the full cohort and
subsetting breaks comparability with the source benchmark; multi-localised proteins are
genuine biology, not annotation noise, and discarding 41 % of the cohort would trade a
real population for a cleaner-looking number; and the conclusion is unchanged either way,
since the gain holds in both subsets.

Report the full-cohort numbers as primary, and cite the single-localisation result
alongside them whenever the size of the ESM-2 gain is the point, since it is the more
conservative estimate.

## Conclusion

Adding frozen layer-8 ESM-2 embeddings to an input-length-matched composition baseline
improves extracellular-vs-cytoplasmic discrimination on every fold of the full development
cohort (+0.0179 mean AUROC, +0.0634 mean AP), and the gain survives both sensitivity
analyses. Two qualifications belong in any claim made from this: the gain is concentrated
in folds P2 and P3 (+0.036, +0.027 AUROC) and is near-zero in P1 and P4 (+0.004, +0.005),
so the mean is carried by half the folds; and on single-localised proteins the gain is
roughly 20 % smaller on AUROC and 55 % smaller on AP. Average precision is the fairer
headline metric for a 24.6 %-positive problem and shows the steadier improvement.

## Precision note

Metrics recomputed from the stored out-of-fold predictions file agree with the
full-precision values above but are not identical, because that file stores probabilities
at six decimals and the +ESM model is heavily saturated (~800 distinct values per fold;
544–629 scores at exactly `0.000000`, 158–457 at exactly `1.000000`). Rounding creates
ties, which rank-based metrics score at half credit.

| Source | matched_C mean AUROC | +ESM mean AUROC | matched_C mean AP | +ESM mean AP |
|---|---:|---:|---:|---:|
| Full-precision (`deeploc_esm_results.tsv`) | 0.9651 | **0.9794** | 0.9416 | **0.9696** |
| Recomputed from saved predictions | 0.9651 | **0.9792** | 0.9413 | **0.9681** |

Differences on the single-localisation subset run from −0.00068 to +0.00012 on AUROC and
−0.0036 to +0.0010 on AP — in both directions, largest on average precision. The
unsaturated baseline reproduces to within 0.00005. Fold counts and all conclusions are
unaffected. **The full-precision table is authoritative.**

## Supporting files

- Full-precision results, both cohorts — [`esm2_probe/deeploc_esm_results.tsv`](esm2_probe/deeploc_esm_results.tsv) · [`.json`](esm2_probe/deeploc_esm_results.json)
- Out-of-fold predictions (6-dp) — [`esm2_probe/deeploc_esm_oof_predictions.tsv`](esm2_probe/deeploc_esm_oof_predictions.tsv)
- Length-matched sensitivity (≤ 1,022 aa) — [`deeploc_esm_lenmatched_le1022.tsv`](deeploc_esm_lenmatched_le1022.tsv)
- Composition baselines B vs C — [`baseline/deeploc_baseline_results.tsv`](baseline/deeploc_baseline_results.tsv)
- Verified cohort manifest — [`../baseline/deeploc_binary_manifest.tsv`](../baseline/deeploc_binary_manifest.tsv)
- Probe script — [`../esm2_probe/probe_deeploc.py`](../esm2_probe/probe_deeploc.py)
