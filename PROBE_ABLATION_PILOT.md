# Bounded localisation-probe ablation pilot — feature 1446, P3

> **Exploratory P3 tuning pilot. Not an evaluation and not a held-out result.**
> **P4 and P0 were never read.** No masking was used anywhere — the proteins' original,
> unmasked sequences were run. Protocol frozen before any outcome was read:
> [`frozen_protocol.json`](probe_ablation_pilot/frozen_protocol.json) (SHA-256 `0516c6d1…`).

Tables: [`pilot_per_protein.tsv`](probe_ablation_pilot/pilot_per_protein.tsv) (24 rows) ·
[`pilot_conditions.tsv`](probe_ablation_pilot/pilot_conditions.tsv) (72 rows) ·
[`pilot_results.json`](probe_ablation_pilot/pilot_results.json) ·
code [`run_probe_ablation.py`](probe_ablation_pilot/run_probe_ablation.py).

---

## 1. What was run

24 P3 pilot proteins (12 extracellular / 12 cytoplasmic), three positions (10, 20, 30), one
position at a time, four conditions each — **clean**, **no-op**, **ablate 1446**, and a
**seeded random direction matched to the ablation's residual norm**. 72 condition rows, 240
forward passes.

The edit is the validated **direct decoder-delta**,
`delta_h = ((topK(a_edited) − topK(a)) @ w_dec) * std`, applied by a forward hook on
`encoder.layer[7]`. The readout is the **frozen layer-26 ESM-only localisation probe** with its
saved scaler and the sweep's pooling rules (residue tokens only; `<cls>`, `<eos>` and padding
excluded), reported as the change in **extracellular log-odds**.

## 2. Verification — all criteria passed

| Check | Criterion | Observed | |
|---|---|---|---|
| Hook site really is `hidden_states[8]` | exact equality | **max abs diff 0.000e+00**, 25 downstream blocks remain | ✅ |
| Clean reproduces the saved P3 predictions | max \|Δp\| ≤ 0.02 | **5.10e-05** (median 2.59e-08) | ✅ |
| " | median \|Δ log-odds\| ≤ 0.10 | **0.00513** | ✅ |
| " | max \|Δ log-odds\| ≤ 0.50 | **0.3387** | ✅ |
| No-op reproduces clean | \|Δ log-odds\| < 1e-6 | **exactly 0.0**, and `‖delta_h‖` **exactly 0.0** | ✅ |
| Random control norm-matched | — | max mismatch **2.0e-06** | ✅ |
| Non-finite values | none | **0** across all condition and per-protein rows | ✅ |

I verified the hook site for `EsmModel` rather than inheriting it from the `EsmForMaskedLM`
validation. The one-sided max log-odds gap of 0.339 is large on the log-odds scale but
corresponds to a probability difference of 5e-5 — it occurs on a saturated prediction, where
log-odds are numerically ill-conditioned. On the probability scale agreement is near-exact.

## 3. How often feature 1446 was active

| Position | Active / 24 |
|---|---|
| 10 | **12** |
| 20 | **12** |
| 30 | **2** |
| **All cases** | **26 / 72** |

**13 of 24 proteins** had 1446 active at ≥ 1 position; **2** had it active at all three. Active
latent values ranged 1.42 – 8.54 (median 3.75), giving ablation residual norms `‖delta_h‖` of
2.29 – 14.62 (median 7.08).

**The 11 proteins where 1446 was inactive at all three positions produced exactly zero edits
and exactly zero change**, for both the ablation and its norm-matched control — the inactive
policy behaving as specified. They are retained in the table and in the primary average.

> **A strong asymmetry worth recording: of the 13 proteins with 1446 active, 12 are
> extracellular and 1 is cytoplasmic.** In this pilot the feature fires almost exclusively in
> extracellular proteins at N-terminal positions.

## 4. Primary result

Change in extracellular log-odds (edited − clean), averaged over the three positions within
each protein, then paired candidate − control across proteins.

| | All 24 proteins | **13 with 1446 active** | 11 inactive |
|---|---|---|---|
| Mean Δ, ablation | +0.004534 | **+0.008370** | 0.000000 |
| Median Δ, ablation | +0.000124 | +0.004369 | 0.000000 |
| Mean Δ, random control | −0.001434 | −0.002648 | 0.000000 |
| Median Δ, random control | 0.000000 | +0.000291 | 0.000000 |
| **Paired candidate − control, mean** | **+0.005968** | **+0.011018** | 0.000000 |
| **Paired candidate − control, median** | **0.000000** | **+0.004768** | 0.000000 |
| SD / min / max | 0.019053 / −0.011202 / +0.087691 | — / −0.011202 / +0.087691 | 0 |
| Sign split (+ / − / 0) | 9 / 4 / 11 | **9 / 4 / 0** | 0 / 0 / 11 |

Within the 13 active proteins the paired difference splits 12 extracellular (mean +0.012114)
against the single cytoplasmic one (−0.002129).

**Descriptive AUROC over the 24 proteins: clean 1.0000, ablated 1.0000, change 0.0000.** This is
at ceiling and carries no information — the probe separates these 12-versus-12 proteins
perfectly before and after, and a 24-protein AUROC could not resolve an effect of this size in
any case.

## 5. Interpretation

**The effect is real in sign but very small, and the pilot does not resolve it.** Ablating 1446
moved the probe's extracellular log-odds more than a norm-matched random direction did, by a
mean of **+0.0110** among the 13 proteins where the feature was active (median +0.0048, 9 of 13
positive). Across all 24 the median is exactly **0**, because 11 proteins contribute structural
zeros. Against clean log-odds large enough to give a perfect AUROC, a shift of ~0.01 is
negligible: it changes no classification and moves no reported metric.

**The direction is the opposite of the naive story.** If 1446 simply carried
"extracellular-ness", removing it should *lower* the extracellular log-odds. It raised it, on
average, in both summaries. With effects this small, 13 informative proteins, and no interval
computed, that sign should not be interpreted — it is recorded so it is not quietly dropped if
a larger run points the other way.

**What this cannot support.**

- **This is Readout A — a manipulation check, not evidence about encoding.** Per §1.3 of the
  specification, a change in the score of a classifier *we fitted*, caused by editing a feature
  that was *selected for label association*, moves partly by construction. It confirms the edit
  reaches the representation the probe reads. It says nothing about what ESM-2 encodes.
- **The random control tests direction at matched magnitude only.** It is not an
  activation-matched SAE feature control, and **1446 still has no eligible one** — so "the
  candidate's identity matters more than a random direction of the same size" is the strongest
  available comparison, and it is weaker than the family-L control the specification prefers.
- **24 proteins, 3 positions, 13 informative cases.** No confidence interval was computed and
  none would be meaningful at this size. This is a pilot that establishes the machinery runs
  correctly and bounds the effect as small.
- **Local, N-terminal, single-layer.** Only positions 10, 20 and 30 were edited, only at layer
  8; the probe reads the representation 18 encoder blocks downstream of the layer-8 edit (25
  blocks remain to the end of the network). Nothing here
  tests other positions, other layers, or propagation.
- **No biological claim.** Per §1.4, nothing here can establish that a protein's location
  changed. The measurement is entirely inside a frozen language model and a linear probe.

**Bottom line for the next step.** The machinery is verified end to end — the hook site is
exactly `hidden_states[8]`, the no-op is bit-exact, clean scores reproduce the saved
predictions, and inactive cases behave as true zero edits. Feature 1446 is active in about a
third of the N-terminal cases tested and almost exclusively in extracellular proteins. Its
ablation produces a probe-score change that exceeds a norm-matched random direction in 9 of 13
informative proteins but is far too small, and far too weakly controlled, to support any claim
about encoding.
