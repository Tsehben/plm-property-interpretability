# Layer-Sweep Protocol — frozen before any sweep outcome was computed

Machine-readable twin: [`v1/frozen_protocol.json`](v1/frozen_protocol.json)
(SHA-256 `0d393b69…2169177a`). Authorised by
[`LOCALISATION_EXPERIMENT_SPEC.md` Amendment A](../LOCALISATION_EXPERIMENT_SPEC.md#amendment-a--layer-wise-decodability-is-separated-from-intervention-effects).

## What this is not

- **Not a replication** of the Check-in 1 four-fold layer-8 result. That held each of P1–P4
  out once; this fits once on P1+P2. Different fitting and evaluation design — the numbers are
  not comparable as a reproduction, and are never presented as one.
- **P4 is not a fresh final test set.** Spec §2.4 already established that P1–P4 are
  development data, and P4 in particular has carried the 400-protein intervention evaluation.
  Every P4 figure here is a development evaluation, labelled as such.
- **Not an intervention sweep.** Amendment A authorises a *probe* sweep only. Editing at any
  layer other than 8 remains conditional exactly as §4.5 wrote it.

## Population and roles

**Eligible core, n = 10,104** — `inclusion_manifest.tsv` with `in_core == 1` (10,139 minus 8
`E1_engineering_protein` and 27 `E2_p0_base_collision`). `E3_coordinates_not_transferable` is
**not** excluded: E3 concerns transfer of UniProt residue coordinates, which this probe never
uses. **P0 is never read.**

| Role | Partitions | Use |
|---|---|---|
| Fit | P1 + P2 | Fit the scaler and the probe |
| Compare | P3 | Compare layers; the selection criterion is computed here |
| Evaluate | P4 | Touched **only after** the P3 comparison and selection rule are frozen |

## What is compared

Every encoder layer **1–33**. HF `hidden_states[i]` = output of transformer block *i* =
fair-esm `representations[i]`; `hidden_states[0]` is the pre-block embedding output and is not
swept. Three feature sets per layer:

| Set | Dims | Definition |
|---|---|---|
| `comp` | 41 | **Exactly** the existing `matched_C` — `fracs(seq[:1022])` + `log(min(len,1022))` + `fracs(seq[:30])`, verified against [`probe_deeploc.py:64-67`](../DeepLoc2/esm2_probe/probe_deeploc.py#L64-L67) |
| `esm` | 2,560 | `concat(` mean over all visible residue tokens, mean over first `min(30,L)` residue tokens `)` |
| `comp_plus_esm` | 2,601 | `concat(comp, esm)` — exactly the existing `c_matchedC_plus_ESM` construction |

Non-standard residues follow the project convention: counted in the denominator, in none of
the 20 numerators ([`baseline_localization.py:65-76`](../DeepLoc2/baseline/baseline_localization.py#L65-L76)).

**Held identical across all layers:** proteins, truncation at 1,022 residues, pooling,
classifier settings, training-only preprocessing, feature construction.

## Classifier

`LogisticRegression(penalty='l2', C=1.0, solver='lbfgs', max_iter=5000, random_state=42)`,
no class re-weighting — identical to [`probe_deeploc.py:87-88`](../DeepLoc2/esm2_probe/probe_deeploc.py#L87-L88).
`StandardScaler` fit on the **P1+P2 rows only** and applied unchanged to P3 and P4.

`n_iter_` and any `ConvergenceWarning` are recorded for every (layer, feature set).
**Settings are never altered selectively by layer.** If a fit fails to converge, the correction
is documented and applied identically to *all* layers and feature sets, and the pre-correction
run is retained.

## Layer selection — frozen before any P4 prediction

> **Criterion:** highest `comp_plus_esm` AUROC on **P3**. **Tie-break:** lower layer index.

The selected layer and **layer 8** are then evaluated on P4. **A replacement layer is never
chosen on the basis of P4.**

## Metrics, weighting and uncertainty — all fixed before outcomes

- **Reported:** AUROC, average precision, n, class counts, 30 % cluster count, largest-cluster
  share.
- **Primary weighting:** protein-weighted (every protein weight 1).
- **Sensitivity weighting:** *equal-cluster-weight* — every protein gets weight
  1 / (size of its 30 % cluster within the evaluated set), so each cluster contributes total
  weight 1. Explicitly labelled a sensitivity analysis, never the primary. This is a different
  question from cluster-held-out splitting (see `HOMOLOGY_AUDIT.md` §4.2).
- **Fitting-set-overlap sensitivity:** repeat the evaluation excluding P4 proteins whose 30 %
  cluster also contains a P1+P2 protein.
- **Toxin-cluster sensitivity:** the selected-layer vs layer-8 P4 comparison is also reported
  with the `Q2VBP4` cluster removed. **Sensitivity analysis, not a replacement primary.**
- **Mixed-class clusters are kept intact.** They are **never split** to obtain class-stratified
  resamples. Resampling is therefore unstratified at the cluster level, and each resample's
  realised class counts are recorded.
- **Undefined metrics:** a resample missing either class gives an undefined AUROC/AP. Such
  resamples are **excluded from the interval and counted in the report**. Never imputed, never
  silently dropped.

### Paired cluster bootstrap

30 % clusters present in the evaluated set are drawn **with replacement, whole** — every
protein of a drawn cluster enters the resample, and a cluster drawn twice contributes its
members twice. **The same resample indices are used for both models in a comparison**, so the
paired difference preserves per-protein pairing. B = 10,000, seed 42, percentile 95 % CI.

> Resampling clusters makes the interval account for the operational grouping of
> `HOMOLOGY_AUDIT.md` §2. It does **not** establish that clusters are independent of one
> another.

## Retrospective intervention reanalysis

The four completed P4 paired contrasts are re-interval-ed from **saved outcomes only** — no
model rerun, no file in `evaluation_p4/` modified. Original eligible cases, aggregation,
intervals and the four-comparison Bonferroni level (98.75 %) are all retained; cluster-aware
intervals and per-contrast cluster counts are added alongside. Labelled retrospective.

## Forbidden in this task

P0 access · new interventions · sequence generation · commits · pushes · altering previous
outputs · choosing a layer using P4 · claiming the best probe layer is the best intervention
layer · claiming any result establishes direct encoding.
