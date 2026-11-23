# P3 intervention pilot — exploratory tuning

**This is a tuning pilot, not an evaluation.** P4 outcomes and P0 were never read. No
sequences generated, no external predictor, no feature screen, nothing committed. All prior
results and the hashed candidate list `070d0f43…` are unchanged.

| artefact | path |
|---|---|
| pilot run | `intervention_pilot/run_pilot.py` |
| aggregation | `intervention_pilot/analyse_pilot.py` |
| frozen configuration and sample | `v1/frozen_configuration.json`, `v1/pilot_sample.tsv` |
| control selection | `v1/control_selection.tsv` |
| raw per-condition records (5,760 rows) | `v1/per_protein/*.json` |
| aggregates | `v1/dose_response.tsv`, `v1/by_mask_position.tsv`, `v1/candidate_vs_control_paired.tsv`, `v1/realised_norm_match.tsv`, `v1/per_protein_summary.tsv`, `v1/condition_availability.tsv`, `v1/active_edit_denominators.tsv`, `v1/analysis_summary.json` |

Configuration hash `ad59f46e0356dba6a2067e5f758f0275c31513d5cd97f9eaa6ef040cd0e28cf1`.

---

## 1. Interpretations carried forward

- **Similar AUROCs do not establish interchangeable features.** The 19 Extracellular-directed
  candidates have discovery AUROCs within 0.04 of one another, and their decoder directions
  are near-orthogonal (median pairwise cosine 0.008). Equal detection performance is not
  equal function and not equal intervention behaviour.
- **Decoder directions differ, but they are not established biological mechanisms.** Distinct
  directions in a 4,096-dimensional dictionary are a geometric fact about the SAE. Nothing
  here shows any of them corresponds to a biological mechanism.
- **Signal-annotation counterexamples support an association, not an exclusive identity.**
  Activation tracks annotation status rather than localisation label, which is what a
  signal-associated feature looks like — and is equally consistent with a feature responding
  to correlated N-terminal biochemistry. It does not establish what the feature *is*.
- **Associations persist within broad length strata; length is not fully controlled.**
  Within-quartile AUROCs are essentially unchanged from unstratified ones, so length does not
  account for the associations. Conditioning on one coarse variable is not control: residual
  within-stratum variation, joint confounding with kingdom, and unmeasured confounders remain.
- **Effective independent sample size remains unknown.** No clustering exists. The 24 pilot
  proteins, like the 5,190 discovery proteins, are treated as independent and are not known
  to be.

---

## 2. Frozen configuration and sample integrity

Frozen **before** any outcome was observed.

| item | value |
|---|---|
| partition | P3 only |
| candidates | 1446, 2564, 3900 (fixed, from the hashed list) |
| sample | 24 proteins, **12 Extracellular / 12 Cytoplasm** |
| filters | original length 60–1022; standard 20 amino acids only |
| exclusions | engineering proteins and p0-collision records, inherited via `in_core == 1` |
| selection rule | deterministic: within class, sort by frozen sequence SHA-256 hex ascending, take the first 12. **No random seed.** No activation or annotation value enters the ordering |
| eligible pool | 2,287 proteins (559 Extracellular, 1,728 Cytoplasm) |
| pool exclusions | 261 outside the length window, 7 p0-collision, 1 engineering |
| **quota shortfalls** | **none** — both 12-protein quotas met |
| mask positions | one-based residues 10, 20, 30, floor(L/2) |
| duplicate/invalid positions | **none** in any of the 24 proteins |
| multipliers | 0, 0.5, **1.0 (no-op)**, 1.5, 2.0 |
| edit | validated direct decoder-delta, hook on `esm.encoder.layer[7]` (output = `hidden_states[8]`) |

Sequence SHA-256 was re-verified against the frozen manifest for all 24 proteins.
**Token-offset checks passed for every position of every protein**: `<cls>` at token 0, so
one-based residue *i* is token *i*, and the token at that index was asserted equal to the
sequence residue before masking.

---

## 3. Control selection and availability

Prespecified matching rule, **no AUROC restriction**: discovery (P1+P2) activation frequency
and mean active magnitude each within ±20% relative, excluding all 20 selected candidates,
ordered by matching distance then feature id.

| candidate | prevalence | mean active magnitude | controls requested | **found** | controls |
|---|---|---|---|---|---|
| **1446** | 0.420 | 2.948 | 5 | **0** | — |
| 2564 | 0.444 | 1.092 | 5 | **5** | 1634, 1574, 2, 2871, 1034 |
| 3900 | 0.951 | 3.261 | 5 | **4** | 431, 2255, 2326, 2914 |

**The shortfall is reported, not engineered away.** The rule was not widened. For 1446 the
binding constraint is magnitude: only **3** non-candidate eligible features have a mean
active magnitude within ±20% of 2.948, and none of those also matches prevalence. The
closest candidates miss magnitude by 56% (f1574) and 64% (f1634) while matching prevalence to
1–2%. For 3900, 4 features satisfy both criteria and a 5th does not exist. **Feature 1446
therefore ran with no matched control at all**, and no candidate-versus-control statement
about it is possible from this pilot.

### Realised perturbation

The control's **latent** delta was matched exactly in magnitude and sign to the candidate's.
That does **not** match the residual-space norm, because decoder rows differ in length:

| candidate | median relative ‖Δh‖ difference | p90 | max |
|---|---|---|---|
| 2564 | 0.015–0.054 | 0.054–0.059 | 0.062 |
| **3900** | **0.292** | 0.354 | 0.354 |

For 3900 the control edit injects a residual-space perturbation **~29% smaller** than the
candidate's at the same latent delta. Any candidate-minus-control difference for 3900 is
therefore partly a difference in edit size, not only in edit identity. This is a design
defect to fix before evaluation, and it has a closed form: with top-K support unchanged,
Δh = Δa · w_dec[f] · std, so rescaling the control's Δa by ‖w_dec[cand]‖ / ‖w_dec[ctrl]‖
matches the norms exactly. It was not applied here because the prespecified rule said to
match the latent delta, and changing the operator mid-run would have been a silent
substitution.

Controls were marked **unavailable** rather than substituted whenever the control latent was
zero at the position (insertion would change the top-K support) or a matched negative delta
would have driven the latent below zero. Availability counts are in
`v1/condition_availability.tsv`.

---

## 4. No-op and edit-integrity checks

| check | result |
|---|---|
| no-op conditions (multiplier 1.0) | 252 |
| **no-op Δh exactly zero** | **yes, all 252** — max ‖Δh‖ = 0.0 |
| no-op KL (AA-20 and full vocabulary) | **exactly 0.0** in every case |
| latent coordinates changed per edit | **exactly 1**, maximum over all 1,147 edits |
| unintended top-K support evictions | **0** |
| unintended top-K support insertions | **0** |
| Δh confined to the edited token position | **yes, all edits** (0 exceptions) |
| non-finite metric cells | **0** |
| token-offset checks | all passed |

**No numerical failures.** Runtime 337.7 s for 1,243 forward passes on MPS, fp32.

---

## 5. Active-edit denominators

The candidate is frequently **inactive** at the masked position, in which case no edit was
made and no position was substituted.

| candidate | protein–position cells | cells with the candidate active | % | proteins with ≥1 active cell |
|---|---|---|---|---|
| 1446 | 96 | 25 | 26.0% | **12 / 24** |
| 2564 | 96 | 24 | 25.0% | **13 / 24** |
| 3900 | 96 | 42 | 43.8% | **17 / 24** |

Activity is strongly N-terminal, which matters for the positional design: 1446 was active at
position 10 (12 cells) and 20 (12 cells) but only once at 30 and **never at floor(L/2)**;
3900 was active at 10, 20 and 30 and **never at floor(L/2)**; 2564 was active mainly at 30
(13 cells), at 20 (9) and at floor(L/2) only twice. **The floor(L/2) position is nearly
empty for all three candidates** and yields almost no usable conditions.

---

## 6. Dose response

Protein-level means (positions averaged within protein first; n = proteins with ≥1 active
cell). All-sampled and active-edit denominators are in §5.

**Candidate 1446** (n = 12 proteins)

| multiplier | ‖Δh‖/‖h‖ | KL AA-20 (nats) | Δp(true AA) | Δ expected KD |
|---|---|---|---|---|
| 0.0 | 0.00524 | 0.00038 | +0.00067 | +0.00909 |
| 0.5 | 0.00262 | 0.00010 | +0.00038 | +0.00456 |
| **1.0 (no-op)** | **0** | **0** | **0** | **0** |
| 1.5 | 0.00262 | 0.00012 | −0.00044 | −0.00423 |
| 2.0 | 0.00524 | 0.00050 | −0.00094 | −0.00795 |

**Candidate 2564** (n = 13)

| multiplier | ‖Δh‖/‖h‖ | KL AA-20 | Δp(true AA) | Δ expected KD |
|---|---|---|---|---|
| 0.0 | 0.00493 | 0.00052 | −0.00095 | −0.01958 |
| 0.5 | 0.00247 | 0.00012 | −0.00046 | −0.00869 |
| **1.0** | **0** | **0** | **0** | **0** |
| 1.5 | 0.00247 | 0.00010 | +0.00043 | +0.00620 |
| 2.0 | 0.00493 | 0.00036 | +0.00083 | +0.01008 |

**Candidate 3900** (n = 17)

| multiplier | ‖Δh‖/‖h‖ | KL AA-20 | Δp(true AA) | Δ expected KD |
|---|---|---|---|---|
| 0.0 | 0.00256 | 0.00004 | +0.00026 | −0.00166 |
| 0.5 | 0.00128 | 0.00001 | +0.00013 | −0.00083 |
| **1.0** | **0** | **0** | **0** | **0** |
| 1.5 | 0.00128 | 0.00001 | −0.00014 | +0.00082 |
| 2.0 | 0.00256 | 0.00005 | −0.00030 | +0.00160 |

Described, not required: the perturbation size is exactly linear in |multiplier − 1|, as the
construction guarantees; KL rises with it roughly quadratically and symmetrically; Δp(true)
and Δ expected KD are **signed and antisymmetric** about the no-op — ablation moves them one
way, amplification the other. 1446 and 2564 move expected hydrophobicity in **opposite
directions** under ablation (+0.009 against −0.020 KD units). No directional hypothesis was
imposed, and hydrophobicity here is a biochemical readout of the masked-token distribution,
not a localisation prediction.

**Absolute magnitudes are very small.** Over all 895 non-no-op active edits: median KL AA-20
**8.1e-06** nats, 99th percentile 1.7e-03, **maximum 8.2e-03**. Top-1 amino acid changed in
**5 of 895** conditions (same 5 for full vocabulary). ‖Δh‖/‖h‖ median 0.0016, max 0.0125.
Probability mass on non-standard and special tokens moved by at most 3.3e-06, so AA-20 and
full-vocabulary KL agree to five decimal places throughout.

---

## 7. Paired candidate-versus-control contrasts

Paired only on conditions where **both** arms were edited at the same protein, position and
multiplier. Controls are averaged within protein over positions **and** over control
features before pairing, so the unit is the protein — never the position, never the control
feature.

KL AA-20, mean paired difference (candidate − control):

| candidate | multiplier | n proteins | n conditions | mean difference | median difference | n proteins candidate greater |
|---|---|---|---|---|---|---|
| 2564 | 0.0 | 5 | 5 | −2.75e-04 | +6e-06 | 3 / 5 |
| 2564 | 0.5 | 7 | 7 | −4.9e-05 | +2e-06 | 5 / 7 |
| 2564 | 1.5 | 8 | 12 | −3.5e-05 | +3e-06 | 5 / 8 |
| 2564 | 2.0 | 8 | 12 | −1.45e-04 | +1.0e-05 | 5 / 8 |
| 3900 | 0.0 | 15 | 73 | −1.7e-05 | +2e-06 | 10 / 15 |
| 3900 | 0.5 | 17 | 124 | +1e-06 | +1e-06 | 13 / 17 |
| 3900 | 1.5 | 17 | 149 | +3e-06 | +2e-06 | 13 / 17 |
| 3900 | 2.0 | 17 | 149 | +1.6e-05 | +7e-06 | 13 / 17 |
| **1446** | — | **0** | **0** | — | — | **no matched control exists** |

**The mean and the median disagree in sign for every 2564 row and for 3900 at ablation.** The
typical protein shows a slightly larger candidate effect while the mean is dragged negative
by a few controls with large effects. With 5–17 proteins and no multiplicity control, this
is an inconsistent outcome, and it is reported as one rather than resolved in favour of
either statistic. For 3900 the contrast is additionally confounded by the 29% norm mismatch
in §3.

---

## 8. Recommendation: what is technically usable for a later frozen protocol

Three distinct judgements, kept apart.

### Numerical usability — what the machinery can deliver

**Usable, and demonstrated here:**

- The edit operator itself. No-op exactly zero in 252 conditions, exactly one latent
  coordinate changed, zero unintended support changes, Δh confined to the edited position,
  zero non-finite values. Nothing about the implementation needs further validation at this
  scale.
- Multipliers 0, 0.5, 1.5, 2.0 with 1.0 as the no-op, at positions **10, 20 and 30**.
- Candidates **2564** and **3900** against matched controls, with the norm correction below.

**Not usable as configured:**

- **Position floor(L/2).** Nearly all candidate latents are zero there — 0 active cells for
  1446 and 3900, 2 for 2564. It contributes almost nothing and should be dropped, or replaced
  by a position chosen to lie inside the active region, which is an N-terminal design.
- **Candidate 1446 against controls.** No matched control exists under the prespecified rule.
  Either the matching rule changes (prospectively, in writing, before evaluation) or 1446 is
  evaluated without a feature-matched control and that is stated.
- **Latent-delta matching for controls.** It leaves a 29% residual-norm gap for 3900.
  Replace it with the closed-form rescale ‖w_dec[cand]‖ / ‖w_dec[ctrl]‖, frozen before the
  evaluation run.
- **A 24-protein sample.** 5–17 paired proteins per cell is too few for any contrast. The
  active-cell rate (25–44%) means the usable n is a quarter to a half of the nominal n, and a
  frozen protocol must size on **active cells**, not proteins.

### Model sensitivity — what the model actually did

A single-feature, single-position edit at layer 8 moves the masked-token distribution by a
median of 8e-06 nats and at most 8e-03 nats, and changes the top-1 amino acid in 5 of 895
conditions. The dose response is clean, monotone in |multiplier − 1| and symmetric, so the
signal is real and well-behaved — it is simply very small. **This is a finding about the
design, not a failure**: editing one latent at one position is a minimal perturbation, and
the pilot shows it is near the floor of what this readout can resolve. A frozen protocol that
wants resolvable effects should edit over a **window of positions** rather than one, and
should expect to need a substantially larger active-cell count at single positions.

### Biological interpretation — what none of this shows

Nothing here speaks to localisation, secretion or relocation. The readout is the model's own
masked-token distribution at the edited position. Hydrophobicity shifts are a property of
that distribution. The opposite ablation directions of 1446 and 2564 are a description of two
features' local effects, not evidence of two mechanisms.

### Conditions I am *not* recommending on the basis of their outcome

The largest effects in this pilot sit at multiplier 2.0 on candidate 1446 at position 20.
That is **not** why it should be kept, and it should not be the basis of selection. The
recommendation above is made on availability, integrity and matching feasibility only, and it
retains the null and inconsistent outcomes: 3900's near-zero contrasts and the sign
disagreement between mean and median in §7 are both reasons to keep those conditions in a
frozen protocol, because they define what "no effect" looks like under this machinery.

**Stopping here.** No evaluation protocol frozen, no P4 contact, nothing committed.
