# Follow-up: does ESM add predictive value among proteins of similar N-terminal hydrophobicity?

> **FOLLOW-UP analysis on previously examined P4 development data** — not a fresh test set.
> Saved predictions only; no refitting, no extraction, no layer re-selection.
> **P0 reserved and untouched.**

Results: [`hydro_band_results.tsv`](hydrophobicity_band_followup/hydro_band_results.tsv) ·
[`hydro_band_results.json`](hydrophobicity_band_followup/hydro_band_results.json) ·
code [`within_hydro_band.py`](hydrophobicity_band_followup/within_hydro_band.py).
Method identical to the charge-band comparison
([`CHARGE_BAND_FOLLOWUP.md`](CHARGE_BAND_FOLLOWUP.md)).

---

## 1. Bands, and a precision problem worth recording

Mean Kyte–Doolittle over residues 1..min(30, L); bands defined on **P1+P2 only** at quantiles
0.20/0.40/0.60/0.80, lower-inclusive, upper-exclusive, outer intervals unbounded.

**The cut points are exactly −5/6, −41/100, 3/100 and 173/300.**

⚠️ **float64 gets the boundary wrong here, and the exact arithmetic was needed to reproduce
your expected counts.** KD values carry one decimal, so every first-30 mean is an exact
multiple of 1/300 and every cut is an exact rational. In float64 the summed mean and the
`np.quantile` result carry *different* rounding error — a protein whose true value is exactly
0.03 computes as `0.029999999999999995` against a cut of `0.03000000000000007`, and so falls
below a cut it actually sits on. **Eight P4 proteins were misassigned this way**, giving 90,720
pairs instead of 90,815. Recomputing values and cuts as exact `Fraction`s fixes it. (The
charge-band analysis was checked for the same failure and had **0 disagreements**, so that
result stands.)

| Verification | Result |
|---|---|
| P4 total | 2,366 |
| **Eligible** | **2,365** |
| **Missing, excluded** | **1** — `P63301`, whose window `MALAVRVVYCGAUGYKPKYLQLKEKLEHEF` contains **U** (selenocysteine), for which KD defines no value |
| Every eligible protein in exactly one band | **Yes** — 443 + 501 + 512 + 383 + 526 = 2,365 |
| **Within-band EX–CY pairs** | **90,815** |

The missing protein is kept in the cohort and excluded only from this comparison; nothing is
imputed.

## 2. The five bands

| Band | n | EX | CY | Pairs | Weight | Clusters | Largest | **ESM** | **Composition** | Δ | Support |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 (< −5/6) | 443 | **17** | 426 | 7,242 | 8.0 % | 352 | 13 | 0.9119 | 0.9024 | +0.0095 | **sparse** |
| 2 [−5/6, −0.41) | 501 | **14** | 487 | 6,818 | 7.5 % | 396 | 5 | 0.9601 | 0.9877 | **−0.0276** | **sparse** |
| 3 [−0.41, 0.03) | 512 | **12** | 500 | 6,000 | 6.6 % | 417 | 5 | 0.9012 | 0.7800 | **+0.1212** | **sparse** |
| 4 [0.03, 173/300) | 383 | 101 | 282 | 28,482 | **31.4 %** | 291 | 7 | 0.9751 | 0.9086 | +0.0664 | ok |
| 5 (≥ 173/300) | 526 | 427 | 99 | 42,273 | **46.5 %** | 300 | 92 | 0.9405 | 0.8571 | +0.0834 | ok |

**Bands 1–3 are sparse** — 17, 14 and 12 extracellular proteins respectively. Their AUROCs rest
on fewer than 20 positives each and should be read as unstable: band 3's +0.1212 and band 2's
−0.0276 are both within the range such small samples produce. They are shown because hiding
them would misrepresent the spread, not because they are informative individually. Together
they carry **22.1 %** of the pair weight.

**Band 2 is the only band where composition beats ESM.** It is also the second-sparsest.

## 3. Primary comparison

Pair-weighted mean of the five within-band AUROCs, weights `w_b = n_EX_b × n_CY_b`, over
**90,815** within-band pairs:

| Model | Within-band AUROC |
|---|---|
| **Layer-26 ESM-only** | **0.947938** |
| **Composition baseline** | **0.881572** |
| **Difference (ESM − composition)** | **+0.066366** |
| **95 % CI** | **[+0.039119, +0.099529]** — **excludes zero** |

Interval from 10,000 whole-cluster draws over the 30 %-identity clusters, the same draw applied
to both models and all five bands, pair weights recomputed inside every draw; seed 42,
percentile.

| Diagnostics | |
|---|---|
| Clusters resampled | **1,546** (one fewer than the full-P4 1,547 — `P63301`'s cluster leaves with it) |
| **Undefined draws** | **0** of 10,000 |
| Valid draws | 10,000 |
| Draws in which a band lost a class | band 1 in **57**, band 2 in **1**; both re-weighted to zero within that draw rather than discarding the draw |

## 4. Interpretation

**Yes, and the advantage is resolved — more decisively than in any comparison so far.** Among
proteins sharing a hydrophobicity band, the layer-26 ESM probe reaches a pair-weighted AUROC of
**0.9479** against composition's **0.8816**, a difference of **+0.0664 with a 95 % interval of
[+0.0391, +0.0995]** that excludes zero by a wide margin.

**Bands 4 and 5 contribute most** — 31.4 % and 46.5 % of the pair weight, **77.9 % together** —
and both show large positive gaps (+0.0664, +0.0834). The primary number is substantially a
statement about those two bands. The three sparse bands supply the remaining 22.1 %, pull in
both directions, and do not drive the result.

**The contrast with the other two comparisons is the interesting part:**

| Comparison | Δ AUROC | 95 % CI |
|---|---|---|
| Unconditioned, all P4 pairs | +0.013298 | [+0.004658, +0.023729] |
| Within charge bands | +0.017653 | [+0.008132, +0.029287] |
| **Within hydrophobicity bands** | **+0.066366** | **[+0.039119, +0.099529]** |

Holding hydrophobicity fixed *widens* the gap roughly fivefold. The reason is visible in the
absolute numbers: restricted to within-band pairs the composition baseline falls to 0.8816
while ESM falls only to 0.9479. N-terminal hydrophobicity is evidently a large part of what the
41-feature composition model was using; take it away and that model loses much more than ESM
does.

⚠️ **Read those absolute numbers carefully.** Within-band AUROCs are computed over a *different
pair set* than the unconditioned AUROCs (0.9696 and 0.9829) — the easy cross-band pairs are
removed, so both models score lower by construction. That is a change of estimand, not a
measured drop in performance. **Only the paired differences in the table above are comparable
across rows.**

**What this does not establish.**

- **Band restriction is approximate.** The bands are wide: band 5 spans everything above
  0.5767 and holds 526 proteins. Holding a property inside a band narrows it; it does not
  eliminate it, and it does nothing to anything else correlated with it.
- **One descriptor at a time.** This conditions on hydrophobicity alone, as the charge analysis
  conditioned on charge alone. Neither is a joint control, and no matching or residualisation
  was performed.
- **Band 5 carries the toxin cluster** — 92 of its 526 proteins sit in one 30 % cluster, and
  band 5 supplies 46.5 % of the weight. The cluster bootstrap accounts for that grouping; it
  does not establish independence between clusters.
- **P4 is not a fresh test set**, and no multiplicity adjustment was applied across this
  comparison and the charge one.

> **Whether the model directly encodes localisation remains open.** Showing that the ESM
> advantage survives — indeed grows under — conditioning on N-terminal hydrophobicity rules out
> one specific confound and shows the composition baseline was leaning on it heavily. It does
> not show that ESM-2 represents subcellular localisation as such. A linear probe recovers what
> is linearly available in a pooled representation, and that remains equally consistent with
> the probe reading some other correlate — family identity among them — that hydrophobicity
> banding does not touch. §1.4 of
> [`LOCALISATION_EXPERIMENT_SPEC.md`](LOCALISATION_EXPERIMENT_SPEC.md) applies unchanged.
