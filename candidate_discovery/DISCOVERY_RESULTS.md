# Bounded exploratory candidate discovery — `candidate_discovery/v1`

Discovery data only (P1+P2 core, n = 5,190). No interventions, no sequence generation, no
external predictor, no layer sweep, nothing committed. P0, P3 and P4 were not read.

| artefact | path |
|---|---|
| discovery run (stages 1–3) | `candidate_discovery/run_discovery.py` |
| characterisation (stage 4) | `candidate_discovery/characterise_candidates.py` |
| complete eligible ranking (2,867 rows) | `v1/eligible_feature_ranking.tsv` |
| the 20 candidates | `v1/candidate_features.tsv` |
| hashed selected list | `v1/selected_feature_ids.json` |
| provenance, config, environment, hashes | `v1/discovery_provenance.json` |
| pooled per-protein summaries (97 MB, 21 shards) | `v1/pooled/` |
| annotation descriptions | `v1/activation_by_class.tsv`, `v1/activation_by_signal_status.tsv`, `v1/signal_interval_enrichment.tsv`, `v1/activation_by_kingdom.tsv`, `v1/activation_vs_length.tsv`, `v1/counterexamples.tsv` |

Configuration hash `33a3b89d20605f4a…`; selected-list SHA-256
`070d0f43a9aea7a4de9102e377d7eb70b3ca7e44578ff709fddd3beecddbb917`; eligible ranking
SHA-256 `ff0ab17f60ead16b867a97d53967637ce6bdba4206a49e5fab49aeec4d3c3443`.

---

## 1. Sample and extraction integrity

| check | result |
|---|---|
| proteins in sample | 5,190 = expected |
| accessions, labels, partitions vs the prepared manifest | asserted row by row, 0 mismatches |
| frozen sequence SHA-256 vs manifest | 5,190 / 5,190 match, 0 mismatches, 0 missing |
| partitions present | {1, 2} only — P0, P3, P4 absent by assertion |
| engineering and p0-collision exclusions | inherited from the prepared view (11 proteins already removed) |
| **reference SAE agreement** | **max abs difference 0.000e+00** on 4 proteins — the authors' path (SAE on the full (L+2) tensor, then slice) against ours (slice, then SAE) |
| key alignment, pooled shards ↔ sample | complete; duplicate-key and set-equality assertions passed |
| non-finite outputs | **0** proteins |
| SAE checkpoint hash | matches `1172503b…` |
| visible residues processed | 2,434,672 |
| runtime | 1,100 s (18.3 min), 4.72 seq/s, MPS, fp32 |
| truncated proteins (original length > 1,022) | 553 (50 Extracellular, 503 Cytoplasm) |

Model `facebook/esm2_t33_650M_UR50D` rev `08e4846e…`, `hidden_states[8]`, fp32, eval,
`inference_mode`, frozen. SAE `liambai/InterProt-ESM2-SAEs` rev `3ded2e51…`,
`esm2_plm1280_l8_sae4096.safetensors`, d_model 1280, d_hidden 4096, k 128, fp32 on CPU.
Tokenisation `<cls>` + L residue tokens + `<eos>`, pooled over `hidden[1:1+L]` only;
special tokens and padding excluded, asserted per batch. Python 3.9.6, torch 2.8.0,
transformers 4.57.6, numpy 2.0.2, scipy 1.13.1, macOS arm64.

**Truncation does not touch the primary summary.** The first-1,022-residue policy removes
only C-terminal residues, so `mean_nterm` (first 30 positions) is unaffected for every
protein. Truncation affects only `mean_visible`, which is descriptive and not used for
selection.

### Eligibility

Rule: a feature must be active (> 0) at least once within the first 30 positions in **≥5%**
of the 5,190 discovery proteins. The SAE's training-time dead list was **not** used.

| | n |
|---|---|
| features total | 4,096 |
| **eligible** | **2,867** |
| ineligible | 1,229 |
| never active in the window | 55 |
| never active anywhere in any visible residue | 23 |
| lowest prevalence among eligible | 0.05010 |

Per-partition prevalence for every eligible feature is in the ranking file
(`prevalence_nterm_P1`, `prevalence_nterm_P2`).

---

## 2. The 20 candidates

Primary summary `mean_nterm`; ranking statistic |AUROC − 0.5| against the frozen DeepLoc
binary label; ties broken by ascending feature id. Signed AUROC retained: > 0.5 means
higher activation associates with **Extracellular**.

| rank | feature | AUROC | \|AUROC−0.5\| | direction | AUROC P1 | AUROC P2 | rank P1 | rank P2 | prevalence | prev P1 / P2 | mean-active magnitude | mean Ext | mean Cyt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 3900 | 0.0230 | 0.4770 | **higher_in_Cytoplasm** | 0.0368 | 0.0145 | 3 | 1 | 0.951 | 0.968 / 0.936 | 3.261 | 0.576 | 3.979 |
| 2 | 642 | 0.9745 | 0.4745 | higher_in_Extracellular | 0.9635 | 0.9810 | 2 | 2 | 0.512 | 0.479 / 0.543 | 1.386 | 2.497 | 0.091 |
| 3 | 1446 | 0.9741 | 0.4741 | higher_in_Extracellular | 0.9651 | 0.9786 | 1 | 5 | 0.420 | 0.379 / 0.458 | 2.949 | 4.357 | 0.158 |
| 4 | 243 | 0.9717 | 0.4717 | higher_in_Extracellular | 0.9594 | 0.9792 | 4 | 4 | 0.509 | 0.464 / 0.550 | 2.272 | 4.121 | 0.128 |
| 5 | 963 | 0.9711 | 0.4711 | higher_in_Extracellular | 0.9593 | 0.9779 | 5 | 6 | 0.464 | 0.437 / 0.489 | 0.879 | 1.436 | 0.051 |
| 6 | 2040 | 0.9673 | 0.4673 | higher_in_Extracellular | 0.9480 | 0.9795 | 12 | 3 | 0.519 | 0.489 / 0.546 | 1.491 | 2.686 | 0.111 |
| 7 | 317 | 0.9670 | 0.4670 | higher_in_Extracellular | 0.9552 | 0.9743 | 8 | 7 | 0.499 | 0.459 / 0.536 | 0.894 | 1.479 | 0.088 |
| 8 | 250 | 0.9666 | 0.4666 | higher_in_Extracellular | 0.9559 | 0.9730 | 6 | 9 | 0.556 | 0.519 / 0.590 | 0.634 | 1.180 | 0.066 |
| 9 | 2737 | 0.9647 | 0.4647 | higher_in_Extracellular | 0.9559 | 0.9696 | 7 | 10 | 0.362 | 0.316 / 0.404 | 0.988 | 1.190 | 0.069 |
| 10 | 627 | 0.9630 | 0.4630 | higher_in_Extracellular | 0.9465 | 0.9736 | 14 | 8 | 0.587 | 0.564 / 0.609 | 0.235 | 0.412 | 0.043 |
| 11 | 1052 | 0.9611 | 0.4611 | higher_in_Extracellular | 0.9483 | 0.9689 | 11 | 11 | 0.455 | 0.424 / 0.484 | 0.655 | 1.051 | 0.038 |
| 12 | 50 | 0.9608 | 0.4608 | higher_in_Extracellular | 0.9522 | 0.9639 | 10 | 13 | 0.335 | 0.291 / 0.375 | 1.419 | 1.652 | 0.068 |
| 13 | 900 | 0.9607 | 0.4607 | higher_in_Extracellular | 0.9546 | 0.9624 | 9 | 15 | 0.486 | 0.453 / 0.516 | 1.453 | 2.326 | 0.144 |
| 14 | 2564 | 0.9571 | 0.4571 | higher_in_Extracellular | 0.9452 | 0.9643 | 15 | 12 | 0.444 | 0.403 / 0.481 | 1.092 | 1.675 | 0.072 |
| 15 | 3731 | 0.9547 | 0.4547 | higher_in_Extracellular | 0.9424 | 0.9625 | 17 | 14 | 0.282 | 0.239 / 0.322 | 0.463 | 0.480 | 0.009 |
| 16 | 1235 | 0.9538 | 0.4538 | higher_in_Extracellular | 0.9448 | 0.9580 | 16 | 19 | 0.427 | 0.372 / 0.477 | 0.624 | 0.886 | 0.051 |
| 17 | 999 | 0.9531 | 0.4531 | higher_in_Extracellular | 0.9384 | 0.9620 | 19 | 16 | 0.502 | 0.462 / 0.538 | 0.402 | 0.646 | 0.048 |
| 18 | 217 | 0.9524 | 0.4524 | higher_in_Extracellular | 0.9471 | 0.9545 | 13 | 21 | 0.419 | 0.382 / 0.453 | 0.317 | 0.396 | 0.042 |
| 19 | 1946 | 0.9523 | 0.4523 | higher_in_Extracellular | 0.9392 | 0.9607 | 18 | 17 | 0.619 | 0.603 / 0.635 | 0.400 | 0.699 | 0.091 |
| 20 | 657 | 0.9513 | 0.4513 | higher_in_Extracellular | 0.9372 | 0.9601 | 21 | 18 | 0.343 | 0.302 / 0.380 | 0.276 | 0.305 | 0.022 |

Full precision is preserved in `candidate_features.tsv` (17 significant digits).

### Stability, reported not enforced

Definitions, stated because the numbers are meaningless without them:

- **Top-20 overlap** = |top-20 by |AUROC−0.5| computed on **P1 alone** ∩ top-20 computed on
  **P2 alone**|, both restricted to the eligible set defined on P1+P2 together.
- **Rank correlation** = Spearman correlation across all 2,867 eligible features between the
  P1 and P2 values of the named statistic, over features, not proteins.

| quantity | value |
|---|---|
| top-20 P1 ∩ top-20 P2 | **19 / 20** |
| selected ∩ P1 top-20 | 19 |
| selected ∩ P2 top-20 | 19 |
| selected ∩ both | 18 |
| Spearman, signed AUROC, P1 vs P2 | **0.899** |
| Spearman, \|AUROC−0.5\|, P1 vs P2 | **0.737** |

No pass quota was applied and the ranking was not altered. These are exploratory discovery
statistics on the data the features were selected from, not held-out performance.

### Shape of the ranking

| \|AUROC−0.5\| ≥ | n eligible features |
|---|---|
| 0.45 | **21** |
| 0.40 | 28 |
| 0.35 | 33 |
| 0.30 | 55 |
| 0.25 | 84 |
| 0.20 | 171 |
| 0.10 | 801 |

**The cut at 20 is not a natural break.** Rank 20 is 0.4513 and rank 21 (feature 2344) is
0.4502. Twenty is the pre-registered budget, not a gap in the data. Twenty-five of the top
100 features are Cytoplasm-directed.

---

## 3. Annotation-based description

Produced **after** the candidate list was written and hashed; the hash was re-verified at the
end of the characterisation run and is unchanged. Nothing below altered the list.

### 3.1 Signal-peptide status (annotation-conditioned discovery, n = 5,009)

> `not_annotated` means a matching UniProt record was retrieved and carries no
> signal-peptide annotation. **It is not a confirmed biological negative.** It is never
> pooled with `unverifiable` and never read as absence.

Mean of `mean_nterm` (prevalence in brackets), selected features:

| feature | annotated_present (n=1,311) | not_annotated (n=3,698) | experimental evidence only (n=240) |
|---|---|---|---|
| 642 | 2.580 (0.99) | 0.062 (0.35) | 2.642 (1.00) |
| 1446 | 4.476 (1.00) | 0.100 (0.22) | 4.425 (1.00) |
| 243 | 4.282 (0.99) | 0.067 (0.34) | 4.516 (1.00) |
| 900 | 2.403 (0.99) | 0.108 (0.31) | 2.415 (1.00) |
| 3900 | 0.487 (0.81) | 4.023 (1.00) | 0.424 (0.76) |

The **experimentally evidenced subset (n = 240)** behaves like the full annotated set for
every candidate — the pattern is not an artefact of the 75% computationally evidenced
majority. Full table: `v1/activation_by_signal_status.tsv`.

### 3.2 Residue-level enrichment inside transferable signal intervals

Scope: 1,306 discovery proteins with an EXACT-boundary signal-peptide feature whose
coordinates are verified to apply to the frozen sequence (1,247 Extracellular, 59
Cytoplasm). 0 skipped for truncation. Comparison: mean activation over residues
[start, end] against all other visible residues, aggregated per protein.

| feature | mean inside | mean outside | mean per-protein difference | n inside > outside | n outside > inside |
|---|---|---|---|---|---|
| 1446 | 5.587 | 0.123 | +5.464 | 1,305 | 0 |
| 243 | 6.162 | 0.010 | +6.153 | 1,301 | 2 |
| 2040 | 3.855 | 0.046 | +3.809 | 1,296 | 4 |
| 900 | 3.342 | 0.022 | +3.320 | 1,302 | 3 |
| 642 | 2.145 | 0.292 | +1.853 | 1,289 | 16 |
| 1946 | 0.966 | 0.006 | +0.960 | 1,305 | 0 |
| 3900 | 0.495 | 0.046 | +0.449 | 903 | 303 |
| **2564** | **0.190** | **1.270** | **−1.079** | **166** | **1,137** |

Eighteen of twenty candidates are strongly enriched **inside** the annotated interval.
Full table: `v1/signal_interval_enrichment.tsv`. This is descriptive co-occurrence: it does
not show that any feature encodes a signal peptide, and it is not causal evidence.

### 3.3 Kingdom and length

Mean activation (prevalence):

| feature | Metazoa | Fungi | Viridiplantae | Other | Spearman vs original length |
|---|---|---|---|---|---|
| 1446 | 1.321 (0.45) | 0.982 (0.33) | 0.817 (0.35) | 2.856 (0.67) | −0.316 |
| 642 | 0.763 (0.52) | 0.560 (0.50) | 0.372 (0.41) | 1.878 (0.79) | −0.322 |
| 963 | 0.447 (0.45) | 0.277 (0.51) | 0.224 (0.35) | 1.168 (0.73) | −0.395 |
| 3900 | 2.990 (0.95) | 3.383 (0.96) | 3.552 (0.99) | 1.842 (0.86) | **+0.599** |

Every Extracellular-directed candidate activates markedly higher in the "Other" kingdom
bucket (342 cohort-wide, 161 in the annotation-conditioned discovery sample) than in
Metazoa, Fungi or Viridiplantae. Those are the Sar, Amoebozoa, Discoba and Metamonada
proteins, and within discovery **65.8%** of that bucket carries a signal annotation against
29.7% Metazoa, 17.5% Fungi and 12.5% Viridiplantae — so the kingdom pattern tracks
annotation prevalence rather than kingdom as such. That is a description, not a
decomposition: kingdom and annotation prevalence are not separated here. Full tables: `v1/activation_by_kingdom.tsv`, `v1/activation_vs_length.tsv`.

### 3.4 Counterexamples

- **59 Cytoplasm-labelled proteins carry a signal annotation, and the candidates fire on
  them at Extracellular levels.** Feature 1446: 4.448 on those 59 against 4.476 on the 1,311
  annotated proteins overall. The features track the annotated signal peptide, not the
  localisation label.
- **53 Extracellular proteins are `not_annotated`, and the candidates largely do not fire on
  them.** Feature 1446 averages 0.950 with prevalence 0.51 there, against 4.476 at
  prevalence 1.00 on annotated proteins. Either these are signal-peptide-independent
  secretion or the annotation is missing — `not_annotated` cannot distinguish the two.
- **Feature 2564 is anti-enriched**: higher outside the annotated interval in 1,137 of 1,306
  proteins. It reaches AUROC 0.957 without being a signal-peptide-interval feature, so
  whatever it tracks sits elsewhere in the N-terminal region.
- **Zero activation despite an annotation**: 3 proteins for feature 1446, 252 for 3900, 62
  for 3731.
- **Extracellular proteins with zero activation**: 19 (feature 1946) to 97 (feature 3731).
- **Cytoplasm proteins above the Extracellular median**: 9 (feature 963) to 110 (feature
  217), and 3,837 for 3900, whose direction is reversed.
- Highest-activating Cytoplasm proteins for feature 1446: `Q9LJU0` (8.94), `P09470` (7.74),
  `Q9VH19` (7.34), `Q9HDU9` (7.28). Lowest-activating Extracellular: `P31230`, `Q9VXN9`,
  `C1PGW0`, `A0A0P1AQV8`, all exactly 0.

---

## 4. Limitations and checks

**No check failed.** Every assertion passed: sample identity, sequence hashes, partition
reservation, token/residue alignment, key alignment, finiteness, SAE reference agreement
(exactly 0), SAE checkpoint hash, and the post-characterisation re-verification of the
selected-list hash.

Limitations, stated as limitations:

1. **These are screening statistics, not inference.** |AUROC−0.5| was computed for 2,867
   features on the data that selected them, with no multiplicity control and no held-out
   split. No significance is claimed and none can be read from the magnitudes.
2. **The candidates look like signal-peptide detectors.** 18 of 20 are enriched inside the
   annotated interval, they fire on signal-bearing Cytoplasm proteins, and they largely do
   not fire on Extracellular proteins without an annotation. The Check-in 1 concern "isn't
   this rediscovering signal peptides" is live and is not answered by this run.
3. **Annotation provenance.** 75% of the signal annotations in this cohort are
   computationally evidenced with an unidentified method. The experimental subset (n = 240)
   agrees, which weakens but does not remove the shared-bias concern.
4. **`not_annotated` is not absence.** The 53 low-activating Extracellular proteins may be
   missing annotations rather than annotation-free biology. Nothing here distinguishes them.
5. **Length is a confound, in one direction.** Raw original length alone reaches AUROC 0.239
   against this label (Cytoplasm median 470 residues, Extracellular 201). Feature 3900
   correlates with length at ρ = +0.599, so its class association is partly length-mediated.
   The 19 Extracellular-directed candidates correlate *negatively* (−0.26 to −0.40) and
   reach |AUROC−0.5| of 0.45–0.47 against length's 0.26, so length alone does not account
   for them — but it is not controlled.
6. **The top-20 boundary is a budget, not a gap** (rank 20 = 0.4513, rank 21 = 0.4502).
7. **No homology control.** No clustering exists for this cohort; discovery proteins may be
   homologous to one another, so the effective number of independent observations behind
   every statistic above is unknown and smaller than 5,190.
8. **Untyped InterPro entries are not treated as families** anywhere in this report, and no
   family analysis was performed here.

**Stopping point.** No intervention strengths chosen, no controls drawn, no evaluation on
P3 or P4.
