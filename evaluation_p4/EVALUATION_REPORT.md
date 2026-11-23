# Frozen P4 intervention evaluation — run once

Run once at the frozen sample size. No further P3 tuning, no new candidates, no altered
strengths, no sequence generation, no external predictor, nothing committed. P0 was never
read. All prior outputs are preserved.

| artefact | path |
|---|---|
| evaluation run | `evaluation_p4/run_eval.py` |
| analysis | `evaluation_p4/analyse_eval.py` |
| frozen configuration, code version, control IDs, sample | `v1/frozen_configuration.json`, `v1/eval_sample.tsv` |
| full-precision per-protein outputs (43,200 rows) | `v1/per_protein/*.json` |
| primary results and intervals | `v1/primary_results.tsv`, `v1/primary_by_class.tsv`, `v1/primary_per_protein.tsv`, `v1/primary_distributions.tsv` |
| coverage, all-sampled outcomes, secondary, top-1 | `v1/coverage.tsv`, `v1/all_sampled_candidate_outcomes.tsv`, `v1/secondary_descriptive.tsv`, `v1/top1_by_arm.tsv` |
| run and analysis records | `v1/run_report.json`, `v1/analysis_summary.json` |

Configuration hash `3ea85a0590681786…`, sample hash `9ecaa7e7c37c82cd…`, both written before
any P4 model output. Code version hashes for all five source files are in the frozen
configuration, as are the control IDs (`control_ids_sha256`), carried unchanged from the
pilot.

**Partition history, recorded as required.** P4 was used in the Check-in 1 analyses — probe
fitting and reporting — so it is **not a previously untouched final test set**. It was
excluded from the current candidate discovery (P1+P2) and from the intervention tuning pilot
(P3). P0 remains reserved from fitting, selection, tuning and scoring, and was not read here.

---

## 1. Sample and integrity

| check | result |
|---|---|
| proteins | **400 — 200 Extracellular / 200 Cytoplasm**, quotas asserted and met |
| selection | within class, ascending frozen-sequence SHA-256, first 200; no activation, annotation or outcome filter |
| eligibility | P4, `in_core == 1` (engineering and p0-collision already excluded), length 60–1022, standard 20 amino acids |
| sequence hashes | all 400 verified against the frozen manifest |
| positions sampled | exactly {10, 20, 30} |
| token-offset checks | **all passed** |
| **no-op preservation** | **5,895 no-op conditions, max ‖Δh‖ = 0.0, max KL (AA-20 and full) = 0.0** |
| latent coordinates changed per edit | **exactly 1** (max over 12,300 edits) |
| unintended top-K evictions / insertions | **0 / 0** |
| Δh confined to the edited token position | **all edits** |
| non-finite metric cells | **0** |
| **control norm matching** | **median relative error 2.91e-07, maximum 6.79e-07** against the frozen 1e-5 tolerance |
| conditions rejected for norm tolerance | **0** |
| conditions rejected for top-K violation | **0** |
| complete accounting | 43,200 rows, every one either edited or carrying an explicit unavailability status |
| runtime | **13,500 forward passes, 2,592.6 s (43.2 min)** |

No implementation failure occurred, and no setting was adjusted at any point during or after
the run.

---

## 2. Primary results — paired KL difference

**"Paired KL difference" is a neutral name.** It is the candidate's KL minus the mean of its
available matched controls' KLs, averaged over positions within protein, then averaged over
proteins. KL is KL(clean ‖ edited), AA-20 renormalised, in nats. **It is not a measure of
localisation specificity.**

Four primary contrasts: 2 candidates × 2 non-zero edit conditions. Protein-level resampling
stratified by class, 10,000 resamples, seed 20261122. Ordinary 95% BCa and
Bonferroni-adjusted 98.75% BCa. No interval required the pre-registered degenerate fallback;
all eight are BCa.

| candidate | edit | n proteins | n paired positions | **estimate (nats)** | 95% BCa | excl. 0 | 98.75% BCa | excl. 0 |
|---|---|---|---|---|---|---|---|---|
| **2564** | **ablation (×0)** | **61** | **61** | **+2.097e-05** | **[+1.13e-05, +3.50e-05]** | **yes** | **[+9.05e-06, +4.01e-05]** | **yes** |
| 2564 | amplification (×2) | 136 | 202 | +8.301e-05 | [−9.13e-06, +2.04e-04] | no | [−3.22e-05, +2.42e-04] | no |
| 3900 | ablation (×0) | 275 | 461 | +6.413e-06 | [−9.63e-06, +2.96e-05] | no | [−1.35e-05, +3.67e-05] | no |
| 3900 | amplification (×2) | 301 | 706 | +7.050e-05 | [−8.22e-06, +4.16e-04] | no | [−1.30e-05, +5.20e-04] | no |

**One of four contrasts excludes zero, and it survives Bonferroni adjustment: candidate 2564
under ablation.** The other three are inconclusive — their intervals contain zero at both
levels.

### Distribution of the protein-level differences

| candidate | edit | n | mean | median | p25 | p75 | positive / negative |
|---|---|---|---|---|---|---|---|
| 2564 | ×0 | 61 | +2.10e-05 | +2.34e-06 | −1.27e-06 | +3.79e-05 | 37 / 24 |
| 2564 | ×2 | 136 | +8.30e-05 | +1.67e-05 | −1.92e-05 | +9.38e-05 | 81 / 55 |
| 3900 | ×0 | 275 | +6.41e-06 | −2.03e-06 | −1.63e-05 | +6.64e-06 | 107 / 168 |
| 3900 | ×2 | 301 | +7.05e-05 | −2.68e-06 | −3.41e-05 | +1.17e-05 | 121 / 180 |

Histogram of the one contrast that excludes zero (2564, ablation, n = 61):

```
 [-6.0e-05,-3.8e-05)  ##                                   2
 [-3.8e-05,-1.5e-05)  ###                                  3
 [-1.5e-05,+7.2e-06)  #################################   33
 [+7.2e-06,+3.0e-05)  #####                                5
 [+3.0e-05,+5.2e-05)  #########                            9
 [+5.2e-05,+7.4e-05)  ##                                   2
 [+7.4e-05,+9.7e-05)  ##                                   2
 [+9.7e-05,+1.2e-04)  #                                    1
 [+1.2e-04,+1.4e-04)  ##                                   2
 [+1.4e-04,+1.6e-04)  #                                    1
 [+1.6e-04,+1.9e-04)                                       0
 [+1.9e-04,+2.1e-04)  #                                    1
```

**Read this honestly.** The interval excludes zero for the **mean**, and the mean is carried
by a right tail: the median protein difference is +2.3e-06, only 37 of 61 proteins are
positive, and the modal bin straddles zero. BCa is the right tool for a skewed statistic and
the result stands as stated — but the claim is about a mean, not about a typical protein.
3900's two contrasts have **more** proteins negative than positive while their means are
positive, which is the same skew without the interval support.

### Per class — reported alongside the pooled result

The **200/200 design is deliberately balanced and does not represent natural cohort
prevalence.** No post-result weighting is applied.

| candidate | edit | stratum | n | estimate | 95% BCa |
|---|---|---|---|---|---|
| 2564 | ×0 | Extracellular | **55** | +2.275e-05 | **[+1.17e-05, +3.79e-05]** |
| 2564 | ×0 | Cytoplasm | **6** | +4.644e-06 | [−2.99e-06, +2.49e-05] |
| 2564 | ×2 | Extracellular | 127 | +8.880e-05 | [−9.02e-06, +2.19e-04] |
| 2564 | ×2 | Cytoplasm | 9 | +1.281e-06 | [−6.61e-06, +1.59e-05] |
| 3900 | ×0 | Extracellular | 98 | −4.058e-06 | [−2.19e-05, +3.77e-06] |
| 3900 | ×0 | Cytoplasm | 177 | +1.221e-05 | [−1.19e-05, +4.93e-05] |
| 3900 | ×2 | Extracellular | 103 | −4.449e-06 | [−2.63e-05, +6.36e-06] |
| 3900 | ×2 | Cytoplasm | 198 | +1.095e-04 | [−1.08e-05, +6.18e-04] |

The 2564 ablation result is **effectively an Extracellular-only result**: 55 of its 61
proteins are Extracellular, because 2564 is active in 51.0% of Extracellular cells and 4.3%
of Cytoplasm cells. The six Cytoplasm proteins cannot support any statement.

---

## 3. Coverage, and what the paired contrasts describe

**Paired contrasts describe the subset of protein–position cells in which the candidate was
active *and* at least one matched control was feasible.** An inactive candidate is a genuine
no-edit outcome and is reported; an unavailable control comparison is **missing** and is
never imputed as zero.

| candidate | edit | cells | candidate active | % | paired | % of active | controls per paired cell |
|---|---|---|---|---|---|---|---|
| 2564 | ×0 | 1200 | 332 | 27.7 | **61** | **18.4** | 1.15 |
| 2564 | ×2 | 1200 | 332 | 27.7 | 202 | 60.8 | 1.27 |
| 3900 | ×0 | 1200 | 706 | 58.8 | 461 | 65.3 | 1.97 |
| 3900 | ×2 | 1200 | 706 | 58.8 | **706** | **100.0** | 3.33 |

The 2564 ablation cell has the **narrowest coverage in the design** — 61 of 332 active cells,
18.4%. The reason is structural: matching the candidate's ablation norm requires a large
negative latent change in the control, which frequently exceeds the control's own latent
value, and those conditions are marked unavailable rather than clipped. The loss is not
random with respect to activation magnitude: it removes exactly the conditions where the
candidate's latent was largest relative to the control's.

By position (all proteins):

| candidate | edit | pos 10 active / paired | pos 20 | pos 30 |
|---|---|---|---|---|
| 2564 | ×0 | 12 / 3 | 130 / **50** | 190 / 8 |
| 2564 | ×2 | 12 / 4 | 130 / 99 | 190 / 99 |
| 3900 | ×0 | 287 / 138 | 228 / 154 | 191 / 169 |
| 3900 | ×2 | 287 / 287 | 228 / 228 | 191 / 191 |

**50 of the 61 paired cells behind the significant contrast come from position 20.** That
contrast is therefore, in practice, a statement about ablating feature 2564 at residue 20 in
Extracellular proteins where a control was feasible. Positions are not independent sample
units and were averaged within protein before any summary; the protein is the unit throughout.

### All-sampled candidate outcomes, separate from the paired contrasts

Protein-level mean KL over active cells, no control subtraction, all 400 proteins in the
denominator:

| candidate | edit | active conditions / 1200 | proteins with ≥1 active position | mean KL | median KL | max KL |
|---|---|---|---|---|---|---|
| 1446 | ×0 | 373 (31.1%) | 202 (50.5%) | 2.98e-04 | 1.02e-04 | 9.90e-03 |
| 1446 | ×2 | 373 (31.1%) | 202 (50.5%) | 3.32e-04 | 1.19e-04 | 9.59e-03 |
| 2564 | ×0 | 332 (27.7%) | 203 (50.8%) | 3.21e-04 | 9.77e-05 | 5.32e-03 |
| 2564 | ×2 | 332 (27.7%) | 203 (50.8%) | 3.53e-04 | 1.07e-04 | 5.15e-03 |
| 3900 | ×0 | 706 (58.8%) | 301 (75.2%) | 1.60e-04 | 3.21e-05 | 1.14e-02 |
| 3900 | ×2 | 706 (58.8%) | 301 (75.2%) | 2.23e-04 | 3.59e-05 | 2.37e-02 |

**Candidate 1446 is descriptive only** — no eligible matched control exists under the
prespecified rule, so it appears here and in the secondary tables and in no contrast.

---

## 4. Secondary descriptive outputs

Not primary claims. Protein-level means over active candidate cells.

| candidate | edit | n | Δ expected KD | Δp(true AA) | Δ log p(true) | ‖Δh‖/‖h‖ | KL full vocab | Δ non-standard mass |
|---|---|---|---|---|---|---|---|---|
| 1446 | ×0 | 202 | **+8.02e-03** | +3.70e-05 | −4.11e-03 | 0.0049 | 2.98e-04 | +2.9e-07 |
| 1446 | ×2 | 202 | **−9.55e-03** | −1.70e-04 | +4.01e-03 | 0.0049 | 3.32e-04 | −1.5e-07 |
| 2564 | ×0 | 203 | −1.62e-03 | +1.75e-04 | +1.99e-03 | 0.0054 | 3.21e-04 | +3.8e-07 |
| 2564 | ×2 | 203 | −7.39e-04 | −1.10e-04 | −1.42e-03 | 0.0054 | 3.53e-04 | −3.1e-07 |
| 3900 | ×0 | 301 | −2.22e-03 | −1.38e-04 | −2.85e-04 | 0.0028 | 1.60e-04 | +1.6e-07 |
| 3900 | ×2 | 301 | **+2.53e-03** | −3.45e-05 | −1.33e-04 | 0.0028 | 2.23e-04 | −1.4e-07 |

Signed hydrophobicity takes approximately opposite values on the two sides of the no-op for
1446 (+8.0e-03 against −9.6e-03 KD units) and for 3900 (−2.2e-03 against +2.5e-03), and does
not for 2564 (−1.6e-03 against −7.4e-04, same sign). Magnitudes are unequal and no symmetry
was tested. Full-vocabulary KL tracks AA-20 KL to three digits throughout, and non-standard
token mass moves by at most 4e-07, so the renormalisation is not doing any work.

**Top-1 changes, denominators separated by arm, never pooled:**

| arm | group | edit conditions | top-1 AA-20 changed | % |
|---|---|---|---|---|
| candidate | 1446 | 746 | 13 | 1.74 |
| candidate | 2564 | 664 | 11 | 1.66 |
| candidate | 3900 | 1412 | 11 | 0.78 |
| control | matched to 2564 | 326 | 1 | 0.31 |
| control | matched to 3900 | 3257 | 28 | 0.86 |

---

## 5. Conclusion

**Evidence of differential model sensitivity — one contrast.** Ablating feature 2564 at a
masked position produces a larger change in the model's masked-token distribution than
ablating activation- and norm-matched control features, by +2.10e-05 nats on average, 95% BCa
[+1.13e-05, +3.50e-05], and the interval still excludes zero after Bonferroni adjustment for
the four primary contrasts, [+9.05e-06, +4.01e-05]. The edit machinery behind it is verified:
no-op exactly zero in 5,895 conditions, one latent coordinate changed, no support violations,
control norms matched to 3e-07 relative. **The licensed claim is exactly this: the model's
masked-token predictions are differentially sensitive to the identity of this feature
compared with matched alternatives.** It is bounded in four ways that belong in the same
sentence: it is a **mean** effect carried by a right tail (median +2.3e-06, 37 of 61 proteins
positive); it rests on **61 of 400** proteins, the narrowest coverage in the design; **55 of
those 61 are Extracellular**; and **50 of the 61 paired cells are at residue 20**.

**Inconclusive — three contrasts.** 2564 under amplification and 3900 under both ablation and
amplification all have intervals containing zero at 95% and at 98.75%. For 3900 the
protein-level medians are *negative* while the means are positive, on 275 and 301 proteins
with good coverage (65% and 100% of active cells paired). These are not evidence of absence;
they are results whose sign is not resolved at this sample size. Candidate 1446 has no
matched control and carries no contrast at all.

**Unanswered biological questions.** Nothing here speaks to protein localisation, secretion,
signal-peptide function or relocation. The readout is the model's own masked-token
distribution at one edited residue, measured in the same forward pass; it does not test
propagation to distant positions, to downstream layers' behaviour on unmasked inputs, to our
own classifier, or to any external predictor. The discovery-side association between these
features and annotated signal peptides remains an association, and this evaluation does not
convert it into a statement about what the features encode. Hydrophobicity shifts are
properties of the masked-token distribution, not localisation predictions.

**Carried-over limitation on uncertainty.** No homology clustering exists for this cohort, so
proteins are treated as independent when they are not known to be. The effective independent
sample size is unmeasured. I am not asserting a known effective n, and I am not asserting
that the intervals err in a known direction — only that an unmeasured dependence is present
and that the intervals do not account for it.

**Stopping here.** The sample was not expanded, nothing was adjusted after seeing results,
and no sequence-variant work is being initiated or cancelled on the basis of whether an
interval contains zero.
