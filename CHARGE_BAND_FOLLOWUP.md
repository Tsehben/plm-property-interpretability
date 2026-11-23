# Follow-up: does ESM add predictive value among proteins of similar N-terminal charge?

> **FOLLOW-UP analysis on previously examined P4 development data.** P4 has already carried
> the 400-protein intervention evaluation, the layer-sweep evaluation and the family-novelty
> follow-up. It is not a fresh test set. Saved predictions only — no refitting, no extraction,
> no layer re-selection. **P0 reserved and untouched.**

Results: [`charge_band_results.tsv`](charge_band_followup/charge_band_results.tsv) ·
[`charge_band_results.json`](charge_band_followup/charge_band_results.json) ·
code [`within_charge_band.py`](charge_band_followup/within_charge_band.py).

---

## 1. Bands and assignment

Charge balance = `(K+R−D−E)/window_length` over residues 1..min(30, L) — a charge balance, not
a pH-dependent net charge (histidine and terminal charges omitted). The five bands were defined
on **P1+P2 only**, at quantiles 0.20/0.40/0.60/0.80, lower-inclusive and upper-exclusive with
unbounded outer intervals.

The cut points are exactly **−2/30, −1/30, +1/30, +2/30**. Because every value is a multiple of
1/30, I recomputed both values and cuts at full precision rather than trusting the 6-decimal
stored column, where a protein sitting exactly on a boundary could land in the wrong band.

| Verification | Result |
|---|---|
| P4 proteins with a missing charge value | **0** |
| Every P4 protein in exactly one band | **Yes** — 475 + 248 + 638 + 342 + 663 = **2,366** |
| Full-precision vs 6-decimal assignment | **0 disagreements** across all 2,366 |
| Band counts vs the frozen `prospective_strata.tsv` | identical |

## 2. The five bands

Each band's AUROC uses only the extracellular–cytoplasmic pairs **inside that band**.

| Band | Range | n | EX | CY | Eligible pairs | Weight | Clusters | Largest | **ESM** | **Composition** | Δ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | < −2/30 | 475 | 26 | 449 | 11,674 | 5.3 % | 370 | 8 | 0.9980 | 0.9608 | **+0.0373** |
| 2 | [−2/30, −1/30) | 248 | 32 | 216 | 6,912 | 3.2 % | 211 | 4 | 0.9666 | 0.9514 | +0.0152 |
| 3 | [−1/30, +1/30) | 638 | 215 | 423 | 90,945 | **41.5 %** | 471 | 47 | 0.9846 | 0.9680 | +0.0165 |
| 4 | [+1/30, +2/30) | 342 | 135 | 207 | 27,945 | 12.8 % | 263 | 26 | 0.9776 | 0.9743 | **+0.0033** |
| 5 | ≥ +2/30 | 663 | 163 | 500 | 81,500 | **37.2 %** | 494 | 19 | 0.9790 | 0.9578 | +0.0212 |

**Bands 3 and 5 carry 78.7 % of the eligible pairs** and therefore dominate the primary number.
Bands 1 and 2 together contribute only 8.5 % despite showing the largest and second-largest
gaps — band 1's ESM AUROC of 0.9980 rests on just 26 extracellular proteins and should not be
read as a precise estimate.

## 3. Primary comparison

Pair-weighted mean of the five within-band AUROCs, weights `w_b = n_EX_b × n_CY_b`, total
**218,976** eligible within-band pairs:

| Model | Within-band AUROC |
|---|---|
| **Layer-26 ESM-only** | **0.981763** |
| **Composition baseline** | **0.964110** |
| **Difference (ESM − composition)** | **+0.017653** |
| **95 % CI** | **[+0.008132, +0.029287]** — **excludes zero** |

Interval by paired cluster bootstrap: whole 30 %-identity clusters resampled across all of P4,
the **same draw** applied to both models and all five bands, with **pair weights recomputed
inside every draw**; B = 10,000, seed 42, percentile.

| Diagnostics | |
|---|---|
| Clusters resampled | **1,547** |
| Undefined resamples excluded | **0** of 10,000 |
| Valid resamples | 10,000 |
| Draws in which any band lost a class | **0** |

For context, the **unconditioned** full-P4 difference between the same two frozen models is
**+0.013299**. Restricting to within-band pairs does not shrink the gap — it is slightly larger
at +0.0177, and the unconditioned value sits inside the within-band interval.

## 4. Interpretation

**Yes — ESM adds predictive value among proteins in similar charge ranges.** Comparing only
extracellular and cytoplasmic proteins that share a charge band, the layer-26 ESM probe reaches
a pair-weighted AUROC of 0.9818 against the composition baseline's 0.9641, a difference of
**+0.0177 with a 95 % cluster-bootstrap interval of [+0.0081, +0.0293]** that excludes zero.
The direction is consistent in **all five bands**, from +0.0033 to +0.0373, so the result is not
produced by one band.

That matters because N-terminal charge is one of the obvious ways a 41-feature composition
model could be separating the classes. Holding it approximately fixed does not remove the ESM
advantage, and does not even reduce it relative to the unconditioned comparison.

**What this does not establish.**

- **Range restriction is approximate, not a removal.** Each band still spans a range of charge
  values, and the bands are wide: band 3 holds 638 proteins across [−1/30, +1/30). Holding a
  property inside a band narrows it; it does not eliminate it, and it does nothing to anything
  else correlated with it — length, hydrophobicity, amino-acid composition beyond charge,
  taxonomy, or family.
- **One descriptor was controlled, not a set.** This conditions on charge balance alone. The
  companion hydrophobicity bands exist in the concept table but are a separate comparison and
  were not run here.
- **The weight is concentrated.** Bands 3 and 5 supply 78.7 % of the pairs, so the primary
  number is largely a statement about those two.
- **Clusters are an operational grouping.** Resampling them makes the interval account for that
  grouping; it does not establish independence between clusters.
- **P4 is not a fresh test set**, and no multiplicity adjustment was applied — this is one
  prespecified primary comparison alongside the five descriptive rows.

> **Whether the model directly encodes localisation remains open.** A linear probe recovers
> what is linearly available in a pooled representation. Showing that this signal survives
> conditioning on N-terminal charge rules out one specific confound; it does not show that
> ESM-2 represents subcellular localisation as such, and it is equally consistent with the
> probe reading some other correlate — family identity among them — that charge banding does
> not touch. §1.4 of [`LOCALISATION_EXPERIMENT_SPEC.md`](LOCALISATION_EXPERIMENT_SPEC.md)
> applies unchanged.
