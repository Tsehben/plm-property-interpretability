# P3 pilot correction — residual-space norm matching

v1 is preserved byte-for-byte. Corrected outputs are in `v2_normmatched/`. No P4 or P0
access, no new proteins, no new candidates, no window interventions, no sequence generation,
nothing committed.

| artefact | path |
|---|---|
| corrected control construction | `intervention_pilot/rerun_norm_matched_controls.py` |
| corrected paired analysis | `intervention_pilot/analyse_v2.py` |
| rebuilt control records | `v2_normmatched/per_protein/*.json` |
| availability, v1 against v2 | `v2_normmatched/control_availability_v2.tsv` |
| paired contrasts | `v2_normmatched/paired_contrasts_v2.tsv`, `v2_normmatched/paired_per_protein_v2.tsv` |
| top-1 denominators by arm | `v2_normmatched/top1_denominators_by_arm.tsv` |
| run and analysis records | `v2_normmatched/rerun_report.json`, `v2_normmatched/analysis_summary_v2.json` |

---

## 1. The corrected construction, and the criterion fixed before the rerun

v1 matched the **latent** delta. The requirement was to match the realised **residual-space**
perturbation. With the validated direct decoder-delta and one edited coordinate *f* at token
position *p*, the support is otherwise unchanged, so

    delta_h[p] = delta_a · w_dec[f] · std[p]        and        ‖delta_h[p]‖₂ = |delta_a| · ‖w_dec[f]‖₂ · std[p]

Given the candidate's **realised** target norm *T*, measured on the same masked input at the
same token, the control's required latent change is

    delta_a_control = sign(delta_a_candidate) · T / ( ‖w_dec[control]‖₂ · std[p] )

using the same per-token `std` the SAE itself computes. The realised norm is then recomputed
from the actual edited forward pass and checked.

Why v1 was off by so much for 3900: decoder row norms differ substantially. 3900's is
**0.599**, while its four controls are **0.424, 0.387, 0.554, 0.421** — so an equal latent
delta produced a residual perturbation up to 35% smaller. 2564's norm is 0.596 against
controls 0.563, 0.633, 0.587, 0.599, 0.623, which is why its v1 mismatch was only 1.5–6%.

### Numerical matching criterion — stated before the rerun

> **| ‖Δh_control‖ − T | / T ≤ 1e-5.**
>
> Justified from fp32 precision and the validated implementation, not from any observed
> effect. The intervention validation measured this operator against a float64 reference at
> median relative L2 error **3.367e-07** and maximum **9.829e-07**, with a minimum of **6.01**
> significant digits against fp32's ~7.2. A tolerance of 1e-5 sits about an order of
> magnitude above that worst case: it admits genuine fp32 rounding and rejects construction
> error.

**Achieved: median relative norm error 2.98e-07, maximum 5.96e-07** across every rebuilt
control condition — inside the criterion by a factor of ~17, and consistent with the
validated operator's own fp32 error. Mean candidate and control ‖Δh‖ agree to 7 digits
(e.g. 3900 at multiplier 2: 4.629 against 4.629, mean paired difference 1.0e-06).

Constraints held, with no silent repair: sign preserved; latents kept non-negative; **zero**
top-K support violations; **zero** conditions failing the norm tolerance. Nothing was clipped,
the matching rule was not widened, and the operator was not switched. A control that could
not reach *T* was marked unavailable.

---

## 2. What stayed fixed, and what became unavailable

Same 24 proteins, same masks, same candidates, same multipliers, same control IDs. Reuse was
verified by SHA-256 on the inclusion manifest, the feature ranking, the selected-feature
list, the DeepLoc manifest and the SAE checkpoint — all identical — and the sample was
rebuilt and asserted equal to the frozen accession list and sequence hashes. **The candidate
arm was therefore reused unchanged from v1**; only the control arm was recomputed.
884 forward passes, 182.5 s.

**1446 remains descriptive only.** No eligible matched control exists under the prespecified
rule, so no candidate-versus-control statement about it is possible. That is unchanged by
this correction.

Availability, v1 (latent-matched) against v2 (norm-matched):

| candidate | multiplier | available v1 | **available v2** | lost | reason lost |
|---|---|---|---|---|---|
| 2564 | 0.0 | 5 | **5** | 0 | — |
| 2564 | 0.5 | 7 | **7** | 0 | — |
| 2564 | 1.5 | 12 | **12** | 0 | — |
| 2564 | 2.0 | 12 | **12** | 0 | — |
| 3900 | 0.0 | 73 | **57** | **16** | required negative latent change now exceeds the control's latent value |
| 3900 | 0.5 | 124 | **109** | **15** | same |
| 3900 | 1.5 | 149 | **149** | 0 | — |
| 3900 | 2.0 | 149 | **149** | 0 | — |

The losses are a direct consequence of the correction: because 3900's controls have smaller
decoder norms, matching its perturbation requires a **larger** latent change, and on ablation
that larger negative change sometimes exceeds the control's own latent value. Those
conditions are marked unavailable rather than clipped. Ablation is where the loss falls,
which is exactly where the contrast matters most.

---

## 3. Corrected paired contrasts

Aggregation, in this order, reported at every level:

1. within each (protein, position, candidate, multiplier) cell, average the **available**
   control features' effects;
2. candidate minus that control average, at that position;
3. average those paired differences **within each protein**;
4. summarise **across proteins**.

Positions with no available control contribute nothing and are counted. Controls per position
average 1.0–1.1 for 2564 and 2.1–3.6 for 3900.

KL AA-20 (nats), protein-level paired differences:

| candidate | mult | n proteins | n paired positions | mean | SE | mean/SE | median | IQR | +/− proteins |
|---|---|---|---|---|---|---|---|---|---|
| 2564 | 0.0 | 5 | 5 | −3.09e-04 | 3.11e-04 | −0.99 | +5.48e-06 | [−3.5e-05, +1.1e-05] | 3 / 2 |
| 2564 | 0.5 | 7 | 7 | −5.58e-05 | 5.86e-05 | −0.95 | +1.94e-06 | [−3.6e-06, +4.6e-06] | 5 / 2 |
| 2564 | 1.5 | 8 | 11 | −4.20e-05 | 5.72e-05 | −0.74 | +2.66e-06 | [−3.8e-06, +8.1e-06] | 5 / 3 |
| 2564 | 2.0 | 8 | 11 | −1.73e-04 | 2.32e-04 | −0.75 | +1.01e-05 | [−2.0e-05, +3.4e-05] | 5 / 3 |
| 3900 | 0.0 | 15 | 27 | −5.44e-05 | 4.44e-05 | −1.22 | −2.89e-06 | [−1.6e-05, +1.2e-06] | 5 / 10 |
| 3900 | 0.5 | 17 | 42 | −4.93e-06 | 4.94e-06 | −1.00 | +4.03e-07 | [−3.5e-06, +1.7e-06] | 9 / 8 |
| 3900 | 1.5 | 17 | 42 | −3.91e-06 | 4.36e-06 | −0.90 | +3.14e-07 | [−2.8e-06, +1.1e-06] | 10 / 7 |
| 3900 | 2.0 | 17 | 42 | −1.43e-05 | 1.81e-05 | −0.79 | +1.59e-06 | [−1.1e-05, +4.6e-06] | 10 / 7 |
| **1446** | — | — | — | — | — | — | — | — | **no matched control exists** |

**Every cell has |mean/SE| < 1.3.** No contrast is distinguishable from zero at this sample
size. The full protein-level distribution for 3900 at multiplier 2 shows why:

    -2.8e-04, -5.5e-05, -5.4e-05, -1.6e-05, -1.1e-05, -2.8e-06, -4.8e-07,
    +1.5e-06, +1.6e-06, +3.0e-06, +3.5e-06, +4.6e-06, +4.6e-06, +1.5e-05,
    +3.4e-05, +4.2e-05, +6.7e-05

Ten of seventeen proteins are positive, and a single protein at −2.8e-04 — an order of
magnitude from the rest — pulls the mean negative. **Mean/median disagreement here is a
statement about a heavy-tailed distribution over 17 proteins, not evidence that the effect is
zero.** The between-protein spread (sd 7.5e-05) is several times the mean difference, which
is a statement about precision, not about the absence of an effect.

Hydrophobicity differences are reported separately and are signed:
3900 moves from mean −3.3e-04 KD units at ablation to +1.9e-04 at amplification, with
median paired differences of −8.6e-04 and +1.5e-03 and protein splits of 4/11 and 12/5.
2564 reverses: +2.7e-02 at ablation, −1.1e-02 at amplification, on 5 and 8 proteins.

**Top-1 changes, denominators separated by arm** (never pooled):

| arm | candidate | edit conditions | top-1 AA-20 changed | % |
|---|---|---|---|---|
| candidate | 2564 | 96 | 0 | 0.00 |
| candidate | 3900 | 168 | 3 | 1.79 |
| candidate | 1446 | 100 | 1 | 1.00 |
| control | (matched to 2564) | 36 | 0 | 0.00 |
| control | (matched to 3900) | 464 | 3 | 0.65 |

---

## 4. Interpretation corrections

**4.1 Small is not "at the resolution limit".** The v1 report said single-position edits sit
"near the floor of what this readout can resolve". That is not established and the evidence
points the other way: the no-op gives **exactly** 0.0 KL in all 252 conditions, the operator's
own fp32 relative error is ~3e-07, and the norm construction reproduces its target to 3e-07
relative. A median KL of 8e-06 nats is orders of magnitude above any numerical floor here.
Corrected statement: *the effects are small in magnitude and are resolved cleanly by the
arithmetic; what limits the contrasts is between-protein variance at n = 24, not numerical
precision.*

**4.2 Sparse activity at floor(L/2) does not establish an N-terminal feature.** The v1 report
inferred that activity "is strongly N-terminal". What was measured is that at four specific
positions — 10, 20, 30 and one sequence-dependent midpoint — the candidates' latents were
mostly zero at the midpoint. Four positions per protein, one of them at a length-dependent
location, cannot characterise a positional profile; the midpoints also vary in local sequence
context, and nothing was sampled between 30 and L/2 or after it. Corrected statement: *at the
four sampled positions the candidates were active mainly at 10, 20 and 30 and rarely at
floor(L/2). A positional profile was not measured.*

**4.3 KL cannot be antisymmetric; the signed quantities are described separately.** KL is
non-negative by construction, so it cannot change sign about the no-op, and any phrasing
suggesting otherwise was wrong. What the data show, stated as two separate facts: (a) KL
AA-20 is zero at the no-op and rises on both sides, approximately with the square of the
perturbation size; (b) the **signed** quantities — Δp(true amino acid), Δ log p(true) and
Δ expected Kyte–Doolittle hydropathy — take approximately opposite signs on the two sides of
the no-op, e.g. 1446 at +0.0091 KD on ablation against −0.0080 on amplification, and 2564 at
−0.0196 against +0.0101. "Approximately opposite", not antisymmetric: the magnitudes are not
equal and no symmetry was tested.

**4.4 Window interventions are a different experiment.** The v1 report suggested editing over
a window of positions as a remedy for small effects. That is not a remedy; it is a different
manipulation, with a different null, a different control construction (the norm match would
have to be defined over the whole window) and a different interpretation. It may well produce
larger effects, but larger effects from a different experiment do not answer the question
this one asks. It is recorded as a possible separate study, not a fix.

---

## 5. Are the corrected controls usable?

**The construction is correct and the controls are technically usable. The contrasts are not
informative at this sample size, and no amount of re-tuning this pilot will change that.**

- Technically: norm matched to 3e-07 relative, zero top-K violations, zero tolerance
  failures, availability fully accounted. The machinery is ready.
- Statistically: every paired cell has |mean/SE| < 1.3 over 5–17 proteins. Nothing is
  detectable, and the pilot was never sized to detect anything.
- 2564's cells (5–8 proteins) are too small to carry any contrast and should not be read as
  evidence either way.
- 1446 has no matched control at all and stays descriptive.

I am not proposing further tuning rounds. The next step is one frozen protocol, run once.

---

## 6. Recommended frozen P4 evaluation protocol — specified, not run

| element | specification |
|---|---|
| **Sample** | P4 core (2,366 after engineering and p0-collision exclusions). Draw a fixed **400 proteins, 200 per class**, by the same deterministic rule used here: within class, sort by frozen sequence SHA-256 ascending and take the first 200. **No activation, annotation or outcome filtering.** Length 60–1022, standard amino acids, applied before sorting. Accession list hashed before any P4 forward pass |
| **Positions** | one-based residues **10, 20, 30**. floor(L/2) is dropped: it yielded 0–2 usable cells here. No new positions are added |
| **Multipliers** | **0 (ablation), 1 (no-op), 2 (amplification)**. 0.5 and 1.5 are dropped; they add conditions without adding a distinct question |
| **Candidates** | **2564** and **3900**, each against its norm-matched controls (5 and 4 features). **1446 descriptive only**, reported with its activation and effect distributions and no control contrast |
| **Control construction** | the §1 formula, with the 1e-5 tolerance frozen as specified. Unavailable conditions marked, never clipped or substituted |
| **Primary outcome** | **ΔKL_spec = KL_candidate − mean(KL_available_controls)**, AA-20 renormalised, in nats, aggregated by the §3 four-step rule, with the **protein** as the unit. One outcome, fixed in advance |
| **Uncertainty** | BCa bootstrap over proteins, 10,000 resamples, 95% interval, seed recorded. Reported for every candidate × multiplier cell whether or not it excludes zero |
| **Denominators** | every table reports all-sampled conditions and active/paired conditions separately, with protein counts, position counts and controls-per-position. Top-1 denominators separated by arm |
| **Secondary, descriptive** | Δp(true), Δ log p(true), Δ expected Kyte–Doolittle, ‖Δh‖/‖h‖, realised norm error, availability. Reported two-sided with no directional hypothesis |
| **Stopping rule** | the protocol is run **once** at n = 400. **The sample is not increased if the result is unfavourable**, and no condition, candidate or outcome is added, dropped or reweighted after P4 data is seen. Any deviation is reported as a deviation |
| **Expected yield** | at the pilot's 25–44% active-cell rate, 400 proteins × 3 positions gives roughly 300–520 active cells per candidate, against 24–42 here — a 12× increase in paired units without changing the design |
| **Limitations carried into the report** | no homology clustering exists, so the effective independent sample size is unknown and below 400 and every interval is narrower than the truth warrants; the controls are activation- and norm-matched but differ in decoder direction, position profile, co-firing structure and encoding; the readout is the model's own masked-token distribution and speaks to local sensitivity, not to localisation, secretion or relocation |

**What a result would and would not license.** If ΔKL_spec excludes zero, the licensed claim
is: *the model's masked-token predictions at the edited position are differentially sensitive
to this feature compared with activation- and norm-matched alternatives.* It would not
license any statement containing the word "localisation". If it includes zero, that is a
reportable null for this design at this sample size, and the protocol ends there rather than
expanding.

**Stopping here.** No P4 contact, no evaluation run, nothing committed.
