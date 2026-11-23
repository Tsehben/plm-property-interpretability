# Localisation Experiment Specification

**Status:** specification only. Nothing in this document has been run. No extraction was
launched, no model was fitted, no feature was selected, no intervention was executed, and
no sequence was submitted to any external service.

**Validated foundation.** The intervention mechanics this specification depends on are
already tested: see [`intervention_validation/FOLLOWUP.md`](intervention_validation/FOLLOWUP.md)
and [`intervention_validation/FOLLOWUP_ADDENDUM.md`](intervention_validation/FOLLOWUP_ADDENDUM.md).
The production edit is the **direct decoder-delta**

```
delta_h = ((topK(a_edited) - topK(a)) @ w_dec) * std
h' = h + delta_h
```

injected at `hidden_states[8]` = the output of `encoder.layer[7]` of a frozen
`facebook/esm2_t33_650M_UR50D`, with 25 transformer blocks remaining downstream. It passes
42/42 checks, agrees with a float64 reference to a median relative L2 error of 3.367e-07,
and is exactly zero under a no-op. That is an **engineering** result. It says nothing
about biology, and this specification is careful to keep the two apart.

---

## 1. Scientific question and contribution

### 1.1 The overarching question, unchanged

> Are solubility and subcellular localisation represented in ESM-2 through **distinct,
> property-associated features**, or through **biochemical composition, regional sequence
> patterns, and biological identity** (family, organism)?

Check-in 1 established *that* frozen layer-8 ESM-2 representations help predict
Extracellular vs Cytoplasm — mean validation AUROC 0.9610 → 0.9789 over a matched
composition baseline, 4/4 folds — and that the gain is **not** attributable to any
particular positions, because the N-terminal pooled vector attends to the whole sequence.
Probing cannot resolve this. A probe measures what is linearly *decodable* from a
representation; it is silent on whether the model itself *uses* that information.

### 1.2 The concrete experiment proposed here

Three stages, each with a defined output:

1. **Identify** a bounded set of candidate layer-8 SAE features by a prespecified,
   discovery-only association rule (§4).
2. **Characterise** what those features co-occur with — annotated signal peptides, family,
   kingdom, N-terminal hydrophobicity, length (§3). This is description, not evidence of
   mechanism.
3. **Test downstream effects** by editing each candidate inside the frozen network with the
   validated implementation and measuring what changes (§5).

Stage 3 is what takes this past linear probing. A probe answers "is it there?"; an
intervention answers "does the rest of the network's computation depend on it?"

### 1.3 Two readouts that must never be conflated

| | **Readout A — our classifier** | **Readout B — native ESM-2** |
|---|---|---|
| What is measured | Change in our localisation probe's score when it is fed the post-intervention pooled representation | Change in ESM-2's own masked-amino-acid distribution from the frozen LM head |
| What it can support | That the quantity our probe weights has moved | That ESM-2's own downstream computation depends on the edited direction |
| Circularity risk | **High.** If the feature was selected because it is label-associated, and the probe weights it, Readout A moves by construction | **Low w.r.t. the probe.** The LM head was never fitted by us and knows nothing about our labels |
| Weakness | Near-tautological on its own | Masked-residue prediction is a weak and indirect proxy for localisation |

Readout A is reported as a **manipulation check** — confirmation that the edit reached the
representation the probe reads — and never as evidence about encoding. Readout B is the
scientific readout. Both are reported; neither is promoted after the fact.

### 1.4 What this experiment cannot conclude

**Neither readout can establish that a protein's biological location changed.** A change in
our classifier is a change in a model we fitted. A change in ESM-2's masked-residue
distribution is a change in a language model's next-token statistics. Subcellular
localisation is determined by physical machinery — signal recognition particle binding,
translocon engagement, signal peptidase cleavage — none of which is in the causal path of
either measurement. The strongest permissible claim is:

> *A layer-8 SAE feature exists to which ESM-2's own masked-token predictions are
> differentially sensitive, in a position-specific way and with a biochemical direction
> matching a hypothesis frozen in advance, and which co-occurs with annotated targeting
> signals.*

That is consistent with a targeting-signal representation. It is not proof of one, and it
is certainly not proof of relocation.

---

## 2. Data roles

### 2.1 What is retained unchanged

- The frozen Swiss-Prot binary cohort: **10,139** proteins in partitions P1–P4
  (2,498 Extracellular / 7,641 Cytoplasm, prevalence 0.2464), verified by SHA-256 in
  [`DeepLoc2/baseline/deeploc_binary_manifest.tsv`](DeepLoc2/baseline/deeploc_binary_manifest.tsv).
- DeepLoc 2.0's **published partition assignments**, used exactly as distributed. No
  reshuffling, no re-partitioning, no stratified re-draw across partition boundaries.
- **P0 (2,712 proteins; 643 / 2,069) stays reserved** — from fitting, from feature
  selection, from hyperparameter tuning, and from scoring. Nothing in §3–§6 reads a P0 row
  for any of those purposes.
- **HPA remains excluded.** The HPA test file has no Extracellular column, so it cannot
  evaluate this contrast at all.

### 2.2 Where the eight engineering proteins fall

All eight proteins used in the intervention validation are in P1–P4 (`split_role =
development_P1_4`). They were selected by an outcome-independent rule — "sequence length
≤ 1022 and in P1–P4, sorted by accession ascending, take the first 8" — recorded in
`validation_report.json` as `selection_used_labels_or_predictions: false`.

| Accession | Partition | Target | Kingdom | Length | Additional localisation |
|---|---|---|---|---|---|
| O42873 | P1 | 0 (Cytoplasm) | Fungi | 570 | yes |
| K7ZGS2 | P2 | 1 (Extracellular) | Metazoa | 63 | no |
| O61866 | P2 | 1 (Extracellular) | Metazoa | 403 | yes |
| O43037 | P3 | 0 (Cytoplasm) | Fungi | 147 | yes |
| H2A0K6 | P4 | 1 (Extracellular) | Metazoa | 238 | no |
| O13747 | P4 | 0 (Cytoplasm) | Fungi | 547 | no |
| O43447 | P4 | 0 (Cytoplasm) | Metazoa | 177 | yes |
| O60462-6 | P4 | 1 (Extracellular) | Metazoa | 555 | no |

They are spread across all four partitions, with **four of the eight in P4** — the partition
proposed below for evaluation. That is the awkward case and it is handled explicitly in §2.4.

### 2.3 Recommended division within P1–P4

Three roles, each a whole partition, so no partition boundary is crossed or split:

| Role | Partitions | n | Extracellular | Cytoplasm | Prevalence | Single-localisation subset |
|---|---|---|---|---|---|---|
| **Discovery** — feature association, candidate and control selection | P1 + P2 | 5,201 | 1,339 | 3,862 | 0.2575 | 3,082 (1,272 pos) |
| **Tuning** — intervention strengths, thresholds, position schemes | P3 | 2,556 | 584 | 1,972 | 0.2285 | 1,481 (553 pos) |
| **Intervention evaluation** | P4 | 2,382 | 575 | 1,807 | 0.2414 | 1,404 (542 pos) |
| → after removing engineering proteins | P4 | **2,378** | **573** | **1,805** | 0.2410 | — |

Kingdom composition of each role:

| Role | Metazoa | Fungi | Viridiplantae | Other |
|---|---|---|---|---|
| Discovery (P1+P2) | 3,317 | 1,167 | 556 | 161 |
| Tuning (P3) | 1,556 | 614 | 268 | 118 |
| Evaluation (P4) | 1,561 | 554 | 204 | 63 |

Rationale for this split rather than any other:

- **Three roles, not two.** Separating strength tuning (P3) from evaluation (P4) is what
  stops the "turn the knob until something happens" failure mode. Without it, strength
  selection and outcome measurement share data and the evaluation is not an evaluation.
- **Whole partitions.** DeepLoc's partitions are the only homology control we have. Any
  split that cuts across them discards it.
- **P4 for evaluation** because it is the largest single partition after P3 and has the
  class balance closest to the cohort mean (0.2414 vs 0.2464).

### 2.4 P1–P4 are not a previously untouched final test

**This must be stated in any report of the result.** All four partitions have already been
used:

- The four-fold cross-validation behind the Check-in 1 probe result trained and scored on
  P1–P4, including P4 (`auroc_P4 = 0.9682` for `matched_C`, `0.9656` for `a_ESM_whole`).
- The composition baselines, the matched-length sensitivity analysis, and the
  single-localisation sensitivity analysis all ran on P1–P4.
- The eight engineering proteins were drawn from P1–P4 and have had activations, logits,
  and intervened forward passes computed on them.

So P4 is a **held-out-within-this-experiment** split, not a fresh test set. It controls
against selecting features on the proteins they are evaluated on; it does **not** control
against the accumulated analyst knowledge of P1–P4 built up across Check-in 1. The only
genuinely untouched Swiss-Prot data is **P0**, which stays reserved — and if a final
confirmatory number is ever wanted, P0 is the one place to get it, at the cost of spending
it once.

### 2.5 How the engineering proteins are excluded from confirmatory results

The four engineering proteins in P4 (H2A0K6, O13747, O43447, O60462-6) are removed from the
evaluation set by an explicit accession blocklist containing all eight accessions, applied
as a filter before any evaluation statistic is computed. The blocklist is written to a file
and hashed alongside the candidate-feature list (§4.5).

The other four (O42873, K7ZGS2, O61866, O43037) sit in the discovery and tuning splits. They
are left there — removing them would change nothing material (4 of 7,757) and would be an
unnecessary deviation — but they are flagged in the discovery output so that if any
candidate feature's association turns out to be driven by them, that is visible.

Reason for the exclusion: those eight proteins have had their intervened outputs inspected
in detail, including which feature is argmax at which position, so they are no longer
outcome-blind. Within the 2,716 token positions of those eight proteins, feature 619 was
active at 100% of positions and was the argmax feature at a mean of 99.21% (minimum 96.92%)
— a sample-specific observation that is now known to the analyst and would bias any
feature-selection or evaluation step those proteins entered.

### 2.6 Homology: what is known and what is not

**Available:** partition membership, exact-sequence duplicate status (zero cross-partition
exact duplicates, per the audit), kingdom, length.

**Not available:** cluster identifiers, pairwise sequence identities, family or domain
assignments, any measurement of residual homology between partitions.

> **Preserving the published partitions does not make the splits family-disjoint, and this
> specification does not claim it does.** DeepLoc 2.0's partitions were constructed by a
> homology-partitioning procedure, so some protection is inherited — but we have not
> measured its strength, we do not have their clustering, and "homology-partitioned by
> design" is a statement about their method, not a verified property of our splits.

§3 proposes the smallest work that would turn this from an assumption into a measurement,
and §5.6 specifies a cluster-level bootstrap that uses the result.

---

## 3. Annotation readiness

### 3.1 Inventory, by provenance

Four categories, deliberately kept apart because they carry different evidential weight.

#### A. Existing annotations — in hand, curator-derived

| Item | Source | Coverage over P1–P4 |
|---|---|---|
| Binary target (Extracellular = 1 / Cytoplasm = 0) | DeepLoc 2.0 Swiss-Prot CSV, 10 one-hot localisation columns | 10,139 / 10,139 |
| `has_additional_localization`, `other_localizations` | same, derived from the other 9 columns | 10,139 / 10,139 (5,967 single-localisation across P1–P4) |
| Kingdom (Metazoa / Fungi / Viridiplantae / Other) | DeepLoc 2.0 CSV `Kingdom` | 10,139 / 10,139 |
| Partition (0–4) | DeepLoc 2.0 CSV `Partition` | 12,851 / 12,851 including P0 |
| Sequence + SHA-256 | DeepLoc 2.0 CSV `Sequence`, hashed locally | 10,139 / 10,139 |

#### B. Computed quantities — deterministic functions of the sequence

These are **features we compute**, not annotations. They carry no independent biological
evidence; they cannot confirm anything the sequence does not already say.

| Item | Status | Note |
|---|---|---|
| N-terminal composition, first 30 residues | **already computed** — it is the `matched_C` baseline | The Check-in 1 positional finding: leucine +8.53 pp enriched in the first 30 residues of Extracellular proteins, against +0.02 pp across the whole sequence |
| Length, visible length (≤1022), length quartile | already in the manifest | 1,111 / 10,139 proteins truncated |
| Kyte–Doolittle hydrophobicity profile, max hydrophobic window, net charge in first 30 | **not yet computed**; trivial, no dependency | Needed for the directional readout in §5.5 |

#### C. Predictor-generated annotations — not ground truth

| Item | Status |
|---|---|
| SignalP (any version) predictions | **not installed, not run.** Would be a *prediction* of a signal peptide, not an annotation of one |
| WoLF PSORT calls | **not installed, not run.** See §6; it is not independent of our data |
| DeepLoc 2.0's own sorting-signal outputs | not obtained; would be the predictions of the model whose training labels we are using |

Any of these, if used, must be labelled as predictor output everywhere it appears. A
predicted signal peptide used to validate a feature that predicts signal peptides is
circular.

#### D. Missing — the real gaps

| Item | Status | Why it matters |
|---|---|---|
| **Signal-peptide annotations with evidence codes** | **zero coverage.** The DeepLoc Swiss-Prot CSV has no signal-peptide column; its columns are `ACC, Kingdom, Partition`, ten localisation one-hots, and `Sequence` | This is the single most important missing item. The Check-in 1 "next test we plan" — whether an internal feature overlaps annotated signal peptides — cannot run without it |
| **Transit-peptide annotations** | zero coverage | The natural negative-control targeting signal: mitochondrial/plastid transit peptides are also N-terminal but route elsewhere |
| **Family / domain annotations (Pfam, InterPro)** | **zero coverage for this cohort.** The 49/100 Pfam pilot reported at Check-in 1 was run on the *solubility* dataset, not on DeepLoc | Family is one of the three routes in the overarching question. Without it, the "biological identity" arm cannot be characterised at all |
| **Sequence-similarity groups / cluster IDs** | zero coverage; no clustering has been run | Needed both to characterise features and to treat correlated proteins honestly in §5.6 |

### 3.2 Missing signal-peptide annotation is not a verified negative

UniProt records a `SIGNAL` feature when curators have evidence for one. It does not record
a feature asserting the *absence* of a signal peptide. Therefore an accession with no
`SIGNAL` feature is **unannotated**, which pools together:

- proteins that genuinely have no signal peptide,
- proteins that have one that has not been curated,
- proteins whose entry predates the relevant curation effort,
- fragments and isoforms where the N-terminus is not the true N-terminus.

The variable must therefore be three-state, and must be reported as such:

| State | Meaning | Usable as |
|---|---|---|
| `annotated_present` | `SIGNAL` feature present, with its ECO evidence class recorded (experimental / by-similarity / predicted-by-UniProt-pipeline) | Positive class |
| `not_annotated` | No `SIGNAL` feature | **Not** a negative class. An unlabelled pool |
| `unverifiable` | The returned sequence's SHA-256 does not match our stored hash — the entry has changed since DeepLoc's snapshot | Excluded, counted, reported |

All signal-peptide analyses are scored as **`annotated_present` vs `not_annotated`**, with
that asymmetry named in the text, and a sensitivity analysis restricted to
experimental-evidence annotations only. No imputation. The denominator is reported, never
hidden.

### 3.3 The smallest additional annotation work

**One job, no GPU, no model.** UniProt REST over the 10,139 P1–P4 accessions. P0 rows are
never a request source. Fields requested:

```
accession, id, reviewed, protein_name,
organism_name, organism_id, lineage,
length, sequence,
ft_signal, ft_transit,
xref_pfam, xref_interpro,
cc_subcellular_location, cc_alternative_products, protein_families,
date_created, date_modified, date_sequence_modified, version
```

`lineage` and the three date/version fields are included for two downstream uses beyond
annotation itself: the organism-type resolution in §6.2, and the exposure bound in §6.3.

Mechanics the fetch must get right, because each has a silent-failure mode:

- **Isoform accessions** (268 of the 10,139 carry a `-N` suffix) are resolved by their own
  accession. The canonical sequence is **never** substituted for an isoform whose own
  sequence could not be retrieved.
- **Sequence verification** against the frozen `sha256` in the manifest, for every entry,
  before any annotation is attached to it.
- **No residue coordinate is transferred to a sequence that did not match.** UniProt feature
  intervals are defined on the canonical sequence, so an isoform whose sequence differs from
  canonical does not inherit them.
- **Nothing is dropped silently.** Sequence mismatches, merged and deleted accessions,
  missing records and unresolved isoforms are each a separate, counted state.
- Caching, resumability, and retry with exponential backoff; raw responses, retrieval
  timestamps and query URLs preserved.

> **Status: implemented and run.** Script
> [`annotation_readiness/fetch_annotations.py`](annotation_readiness/fetch_annotations.py);
> results in [`annotation_readiness/out/`](annotation_readiness/out/) and the coverage
> report [`annotation_readiness/ANNOTATION_COVERAGE_REPORT.md`](annotation_readiness/ANNOTATION_COVERAGE_REPORT.md).

**Second priority, not part of the minimum:** MMseqs2 `easy-cluster` at 30% identity /
80% coverage over the 10,139 sequences, producing a cluster ID per protein. Minutes on CPU,
no external dependency beyond the binary. This is what converts §2.6's assumption into a
measurement — it yields both the within-partition redundancy and the cross-partition
homology that is currently unmeasured — and it supplies the resampling unit for §5.6.

---

## 4. Candidate-feature selection

### 4.1 Pooling and residue-position handling, stated first

Per-residue SAE activations are **not stored anywhere** and will not be: at 128 non-zeros ×
4 bytes × 8.9 M residues that is ≈ 9 GB for the two cohorts, and nothing here needs them
kept. The discovery pass therefore computes per-residue activations on the fly and writes
only per-protein, per-feature summaries.

For protein $p$, feature $j$, over **visible residue positions only** (tokens 1 … min(L,
1022); `<cls>`, `<eos>` and padding excluded, matching
[`esm2_common/esm_extract.py:112`](esm2_common/esm_extract.py#L112)):

| Summary | Definition |
|---|---|
| `mean_whole[p,j]` | mean of $a_{p,i,j}$ over all visible positions $i$ |
| `mean_nterm[p,j]` | mean of $a_{p,i,j}$ over positions $1 \ldots \min(30, L)$ |
| `frac_active_nterm[p,j]` | fraction of positions $1 \ldots \min(30,L)$ at which $j$ is in that position's top-128 |
| `frac_active_whole[p,j]` | the same over all visible positions |
| `max_act_nterm[p,j]` | max of $a_{p,i,j}$ over positions $1 \ldots \min(30,L)$ |

Storage: 4,096 features × 5 summaries × 5,201 discovery proteins ≈ **426 MB** in fp32, or
213 MB in fp16. Acceptable; per-residue storage would not be.

**Why N-terminal pooling is the primary and whole-sequence is the control.** The Check-in 1
finding that motivates this experiment is positional. But the honest caveat from Check-in 1
applies to the SAE summaries exactly as it applied to the embeddings: *positions 1–30 at
layer 8 attend to the whole visible sequence*, so `mean_nterm` is not a local quantity. It
is "the representation at N-terminal positions", not "the representation of the N-terminal
residues". Both poolings are computed and both are reported; the contrast between them is
informative, and neither alone licenses a positional claim. The positional claim is
established by the **intervention** position schemes in §5.3, not by the pooling.

### 4.2 The discovery rule — one rule, prespecified, discovery data only

Computed on **P1 + P2 only** (5,201 proteins). P3 and P4 are not read.

1. **Eligibility.** Retain feature $j$ if `frac_active_nterm[p,j] > 0` for at least **5%**
   of discovery proteins. This is an activity filter measured on our own data.
   - It deliberately does **not** use the SAE's training-time dead-feature list. 338 of
     4,096 latents were dead at the end of training, but 4,010 were ever active on our
     pilot proteins, so at least 252 of those 338 fire here. Excluding features by the
     training log would have been an error; excluding them by measured inactivity on our
     cohort is not.
2. **Association statistic.** For each eligible $j$, compute the AUROC of `mean_nterm[·,j]`
   for Extracellular vs Cytoplasm across the discovery proteins. Record the same for
   `mean_whole[·,j]`.
3. **Ranking.** Rank eligible features by $|\text{AUROC}_{\text{nterm}} - 0.5|$, descending.
4. **Bounded candidate count: the top 20.** Fixed before any discovery number is looked at.
   Twenty is chosen so that the §5 evaluation stays inside the runtime budget in §7 while
   leaving room for the specificity structure; it is not tuned.
5. Record, for every candidate: AUROC (both poolings), activation frequency, mean and max
   activation when active, and the direction of association (which class it is higher in).

### 4.3 Stability check

Re-run steps 1–4 **independently on P1 alone and on P2 alone**. Report:

- The overlap of the two single-partition top-20 lists with each other, and of each with the
  combined list (counts and Jaccard index).
- The Spearman correlation of $|\text{AUROC} - 0.5|$ across all commonly-eligible features
  between P1 and P2.

**Prespecified bar: at least 10 of the 20 combined-list candidates must appear in both
single-partition top-20 lists.**

If fewer than 10 do, the candidate set is partition-unstable. In that case the instability
**is the finding** and is reported as such; the claim is reduced to "no stable candidate
feature at layer 8 under this rule". The rule is not re-tuned until it produces stability —
doing so would make the stability check vacuous. Any revision is documented as a deviation,
with the original rule's result reported alongside it.

### 4.4 Control features — chosen from discovery-only statistics

**Two control families, answering two different questions.** They are specified separately
and never merged, because a control that is conditioned on near-chance discovery AUROC is
not a sample from the dictionary at large.

#### 4.4a Control family L — **low-association controls** (the primary specificity control)

One matched control per candidate. A control for candidate $c$ is the feature $m$ that
minimises a matching distance subject to all three constraints:

| Constraint | Tolerance |
|---|---|
| Activation frequency `frac_active_nterm` matched to $c$ | within ±20% relative |
| Mean activation magnitude when active, matched to $c$ | within ±20% relative |
| **Low discovery association** | $\lvert\text{AUROC}_{\text{nterm}} - 0.5\rvert \le 0.02$ on P1+P2 |

Ties broken by lowest feature index, so the choice is deterministic and reproducible. If no
feature satisfies all three for some candidate, that candidate is reported as
**unmatchable** and dropped from the specificity analysis rather than matched loosely.

> **These are low-association controls. They are not an unbiased random-feature null.** The
> third constraint selects on a statistic measured on the discovery split, so family L is a
> conditioned subset of the dictionary, not a draw from it. Two consequences are stated
> wherever family L is used. First, the AUROC filter is itself estimated with sampling
> error, so a feature admitted at $\lvert\text{AUROC} - 0.5\rvert \le 0.02$ on P1+P2 is a
> feature whose *estimated* association is near chance, not one whose true association is
> zero; selecting on an estimate near a bound also selects, in part, on noise in that
> estimate. Second, the question family L answers is narrow and must be asked in those
> words: *does the candidate do more than a comparably-active feature that is not
> label-associated in discovery?* It does not answer what a typical feature of that
> activity profile would do, and a candidate-minus-L difference must never be reported as
> "versus a random feature".

#### 4.4b Control family M — **activation/norm-matched controls, no association filter**

A second control family, specified now and kept distinct: the same ±20% matching on
activation frequency and mean active magnitude, with **no constraint on AUROC at all**.
Selection among the features that satisfy the two matching constraints is by a fixed seed
recorded in the pre-registration, drawing $n_M$ controls per candidate (proposed $n_M = 5$),
so family M is a sample of the comparably-active dictionary rather than one chosen feature.

The two families answer different questions and are reported in separate rows, never pooled:

| Family | What is conditioned on | Question it answers |
|---|---|---|
| **L** (low-association) | activity, magnitude, **and** near-chance discovery AUROC | Does the candidate do more than a comparably-active feature that discovery says is not label-associated? |
| **M** (matched, unfiltered) | activity and magnitude only | Does the candidate do more than a comparably-active feature drawn without regard to label association? |

A candidate that beats L but not M, or M but not L, is an informative result and is reported
as such rather than resolved in favour of whichever is larger. Neither family is drawn and
no candidate feature is selected at this stage.

**Perturbation magnitudes for controls are the same multiplicative factors as for
candidates** (§5.1). Because activation magnitudes are matched to ±20%, the resulting
injected deltas are comparable in size — but *comparable* is not *equal*, so the realised
$\lVert \Delta h \rVert_2$ is recorded for every single condition and reported. The match is
then something a reader can check rather than something they must take on trust.

**The match is two scalars — plus, for family L only, one selection constraint.** A control
feature agrees with its candidate on activation frequency and mean active magnitude to ±20%;
a family-L control is additionally *selected* for near-chance discovery AUROC, which is a
conditioning step on a noisy estimate, not a matched quantity. It is a different feature in
every other respect: a different decoder
direction, a different set of proteins and positions it fires on, a different correlation
structure with the rest of the dictionary, and an unknown encoding. Those differences are
not controlled, and they are the standing alternative explanation for any
candidate-minus-control difference. Every claim built on this contrast names the matched
quantities (§5.3).

> **Intervention-evaluation outcomes are never used to choose features or strengths.** The
> candidate list, the control list, the eligibility set, and the engineering-protein
> blocklist are written to disk and SHA-256 hashed **before** any intervention runs on P4.
> The hashes go in the report. If any of them changes afterwards, that is a deviation and
> is reported as one.

### 4.5 Layer: 8 stays, and a sweep is optional

> ⚠️ **SUPERSEDED IN PART by [Amendment A](#amendment-a--layer-wise-decodability-is-separated-from-intervention-effects) (2026-11-22).**
> The paragraph below beginning *"A layer sweep is optional, not necessary"* made a layer
> sweep conditional on the layer-8 **intervention** first yielding a reproducible specific
> effect. Amendment A removes that condition for a **probe** sweep, on the grounds that the
> two measure different things. **The text below is preserved verbatim as the historical
> record and is not edited.** Everything in it about layer 8's validated status, about SAE
> checkpoints, and about interventions remains in force.

**Layer 8 is retained as the validated starting point**, for three independent reasons:

1. It is the only layer for which we have verified our activation tensor against the
   released InterProt SAE exactly — reference agreement 0.000e+00, active counts exactly
   k = 128 everywhere, alignment confirmed by both ID and hash on 128/128 pilot proteins.
2. It is the layer the Check-in 1 probe result was measured at, so the intervention and the
   probe evidence are about the same representation.
3. It is the layer the intervention implementation was validated at, with the hook site
   established empirically rather than assumed.

**A layer sweep is optional, not necessary**, and is not launched by this specification. It
answers a *different* question — where in the network the feature lives — and it costs a
separate SAE checkpoint per layer, each needing its own reference-agreement validation
before any result from it could be trusted. The condition under which it becomes warranted
is narrow: **only if the layer-8 experiment yields a candidate with a reproducible,
specific effect**, at which point "is the effect layer-local?" becomes a real question.

If one cross-layer check is ever made, the more valuable one is **layer 24**, not a sweep.
The layer-8 InterProt checkpoint ships without the layer-24 checkpoint's documentation and
validation, so a feature-level claim at layer 8 is worth checking against the
better-documented layer before it is relied on. This remains optional.

---

## 5. Intervention experiment

### 5.1 Strengths, and how tuning is separated from evaluation

**Candidate strength set, prespecified:** multiplicative scale $s$ applied to the candidate
feature's activation at the edit positions, $s \in \{0, 0.5, 2, 4\}$ — ablation, halving,
doubling, quadrupling — plus $s = 1$ as the no-op.

**Tuning happens on P3 and only on P3.** On the tuning split, measure for each candidate and
each $s$:

- the realised $\lVert \Delta h \rVert_2 / \lVert h \rVert_2$ at edit positions, and
- the number of **unintended** top-K evictions — features other than the edited one whose
  top-128 membership changes.

Two prespecified admissibility criteria, both engineering, neither referring to any outcome:

| Criterion | Threshold | Basis |
|---|---|---|
| Non-degenerate | $\lVert \Delta h \rVert_2 / \lVert h \rVert_2 > 10^{-3}$ | Three orders of magnitude above the implementation's own relative error (median 3.4e-07, max 9.8e-07), so the edit is unambiguously signal, not rounding |
| Support-clean | unintended evictions ≤ 1 per edited position, mean over the split | Keeps the edit interpretable as a single-feature change |

The admissible strength set is then **frozen in writing** before P4 is touched.

> **Strengths are never raised to obtain a larger effect.** If no admissible strength
> produces a measurable outcome on P4, that is the result. Escalating until something moves
> would convert a null into an artefact of the search.

**Top-K semantics, carried over from the validation and reported per condition.** The edit
operates inside a top-128 selection, and the three operations are not symmetric:

| Operation | Support effect |
|---|---|
| Amplify an already-active feature ($s>1$) | Support-preserving — still 128 non-zeros, no eviction |
| Ablate an active feature ($s=0$) | Leaves 127 non-zeros; no other feature is evicted |
| **Insert** an inactive feature | Evicts the weakest incumbent — a genuinely different operation |

Insertion conditions, if run at all, are reported as a separate arm with their eviction
counts, never pooled with amplification.

### 5.2 Edit positions and mask positions

**Edit position schemes** (which residue positions receive $\Delta h$):

| Scheme | Edit positions | Purpose |
|---|---|---|
| `E-nterm` | visible positions 1–30 | Hypothesis-relevant |
| `E-all` | all visible positions | Tests whether the feature matters only where it fires, or anywhere |

**Mask positions** (where the masked-residue distribution is read). For each protein, a
fixed, label-independent rule:

| Scheme | Mask positions |
|---|---|
| `M-nterm` | $\{5, 10, 15, 20, 25, 30\} \cap$ visible |
| `M-distal` | the mirrored positions from the C-terminus: $\{L', L'-5, \ldots, L'-25\}$ where $L' = \min(L, 1022)$ |

One mask per forward pass — a single position is masked, the rest of the sequence is intact.
Proteins shorter than 60 visible residues are excluded from the position analysis so that
`M-nterm` and `M-distal` do not overlap; the count of such proteins is reported.

The combination `E-nterm` × `M-distal` is the positional-specificity test: the edit is in
the N-terminus, the measurement is at the other end of the protein.

### 5.3 Controls

| # | Control | What it isolates |
|---|---|---|
| 1 | **No-op** ($s = 1$) | Implementation integrity. Must give *exactly* zero delta — verified exactly zero in both implementations, re-checked every run as a regression |
| 2L | **Low-association control feature** (§4.4a), same positions, same $s$ | **The primary specificity control.** Does the candidate's *identity* matter, or would a comparably-active feature that discovery says is not label-associated do as well? Conditioned on discovery AUROC — **not** a random-feature null |
| 2M | **Activation/norm-matched control features, no AUROC filter** (§4.4b), same positions, same $s$ | Does the candidate's identity matter against comparably-active features drawn *without* regard to label association? Reported separately from 2L, never pooled with it |
| 3 | **Random-direction control** — inject $\Delta h$ of the *same measured L2 norm* as the candidate's, in a direction drawn as a random unit combination of `w_dec` rows (fixed seed, recorded) | Does the *direction* matter, or would any perturbation of that size do? |
| 4 | **Reconstruction replacement** — replace $h$ with the SAE reconstruction, make no edit | A **diagnostic reference** on the SAE parameterisation, explained below. **Not** a specificity control |

**Control 4 is a diagnostic, not a specificity criterion.** Our edit is error-preserving:
$h' = h + \Delta h$, so it does **not** pass $h$ through the SAE and does **not** incur
reconstruction damage. Reconstruction replacement is a different manipulation from the edit:
it answers "how much downstream disturbance does the SAE's reconstruction error alone
cause?", which is a property of the dictionary, not of any feature. From the validation, on
8 proteins: mean KL(clean ‖ reconstructed) = **0.293 nats**, with the top-1 prediction
flipped in **2 of 8** proteins when all residues were reconstructed.

It is reported alongside the edit effects because it is informative about the SAE we are
using — a dictionary whose reconstruction error already perturbs downstream predictions
substantially is worth knowing about. It is **not** a threshold.

> **Specificity is assessed only against matched interventions — controls 2L, 2M and 3.** The
> comparison that licenses a specificity claim is candidate versus *magnitude-matched
> control feature* and candidate versus *norm-matched random direction*. What those
> controls match is stated exactly, and so is what they do not.
>
> - **Control 2L matches** `frac_active_nterm` and mean activation magnitude when active, to
>   within ±20% relative, and additionally *conditions* on $\lvert\text{AUROC}_{\text{nterm}}
>   - 0.5\rvert \le 0.02$ on the discovery split — a selection on a measured statistic, so
>   2L is a low-association control and not a random-feature null (§4.4a). **Control 2M
>   matches** the same two quantities with no association filter (§4.4b). **Neither matches**
>   the control feature's decoder direction, its position profile along the sequence, which
>   proteins it fires on, its correlation with the candidate or with other features, its
>   activation distribution beyond the matched mean, or whatever it encodes.
> - **Control 3 matches** the realised $\lVert \Delta h \rVert_2$. **It does not match**
>   direction, and a random unit combination of `w_dec` rows is not a sample from the
>   distribution of directions the model itself uses.
>
> So a candidate-minus-control difference is evidence that *this particular* feature or
> direction matters more than the matched alternative **along the matched quantities** — not
> evidence that the two interventions were identical in every respect but one. Matching
> activation frequency or perturbation norm does not make two features otherwise the same.
> The unmatched differences remain live alternative explanations, so the realised
> $\lVert \Delta h \rVert_2$, both matched statistics, and the per-condition values are
> recorded and reported for every condition, and any specificity claim is stated with the
> matched quantities named. Reconstruction replacement differs from the edit in the
> manipulation itself, so an edit effect being larger or smaller than it implies nothing
> about whether the edit was feature-specific. A small but reliably feature-specific effect
> is a real finding; a large but non-specific one is not. Size relative to control 4 is not
> evidence either way.

### 5.4 Primary outcome — one, defined operationally

> **ΔKL_spec.** For each evaluation protein $p$: let $K_{\text{cand}}(p)$ be the mean over
> that protein's `M-nterm` mask positions of
> $\mathrm{KL}(P_{\text{clean}} \,\|\, P_{\text{intervened}})$ over the **20 standard amino
> acids, renormalised** (token ids 4–23, re-softmaxed), in **nats**, under candidate
> ablation ($s = 0$) with edit scheme `E-nterm`. Let $K_{\text{ctrl}}(p)$ be the same
> quantity for that candidate's matched control feature. Then
> $$\Delta\mathrm{KL}_{\text{spec}}(p) = K_{\text{cand}}(p) - K_{\text{ctrl}}(p)$$
> and the primary statistic is the **mean of $\Delta\mathrm{KL}_{\text{spec}}(p)$ over
> evaluation proteins**, with a 95% confidence interval.

Fully specified: direction of the KL (clean first), vocabulary (AA-20 renormalised), units
(nats), aggregation (within-protein mean over mask positions, then across proteins),
condition (ablation, N-terminal edit), and the comparison it is a difference against.

The full 33-token vocabulary KL is computed and reported too, but AA-20 is primary: the
non-amino-acid tokens (`<cls>`, `<pad>`, `<mask>`, `<unk>`, `X`, `B`, `U`, `Z`, `O`, `.`,
`-`) carry probability mass that is not a biological prediction, and letting it into the
primary metric would let a shift in, say, `<unk>` mass register as a biological effect.

#### What ΔKL_spec actually measures

> **ΔKL_spec is a measure of differential sensitivity of masked-token predictions. It is not
> a measure of localisation specificity.**

Unpacked, because the distinction is easy to lose in the reporting:

- KL is **unsigned and content-free**. It records that the predicted residue distribution
  moved, not which way, not toward what, and not whether the movement has anything to do
  with localisation. A candidate feature whose ablation reshuffles probability among
  chemically unrelated residues scores the same as one that shifts it coherently.
- The candidate-minus-control difference establishes that **ESM-2's masked-token predictions
  are more sensitive to perturbing this feature than to perturbing a matched control**. That
  is a statement about the model's internal dependence structure at layer 8.
- Getting from there to localisation requires further, separate evidence: the positional
  structure of the effect (S2), its biochemical direction (S3), and the independent
  characterisation of what the feature co-occurs with (§3). **Each of those is reported as
  its own result, never folded into ΔKL_spec.**
- Consequently, S1 alone licenses only: *the model's masked-token predictions are
  differentially sensitive to this feature.* It does not license any sentence containing
  the word "localisation".

### 5.5 Secondary measurements — all pre-listed

None of these may be promoted to primary after the fact.

| Measurement | Definition | Note |
|---|---|---|
| Top-1 change rate | Fraction of mask positions where the argmax residue changes (AA-20 and full vocab) | **Low-powered.** In the validation, **0 of 64** edit conditions changed top-1 under either vocabulary. Reported for completeness, not relied on |
| Total variation distance | $\tfrac12 \sum_{a \in \text{AA20}} \lvert P_{\text{clean}}(a) - P_{\text{int}}(a)\rvert$ | Bounded in [0,1]; complements KL, which is unbounded |
| **Expected hydrophobicity shift** | $\sum_a P(a)\cdot \text{KD}(a)$ under clean vs intervened, at `M-nterm` positions | **A biochemical outcome**, signed, unlike KL. It measures a change in the hydrophobicity of the predicted residue distribution — nothing more. See the directional-hypothesis rule below |
| Probe-score shift (**Readout A**) | Change in the localisation probe's output on the post-intervention pooled representation | **Manipulation check only.** Circular if the candidate was selected by label association — which it was. Reported, never used as evidence of encoding |
| Realised $\lVert\Delta h\rVert_2$, unintended eviction count | Per condition | Engineering telemetry that makes the control matching auditable |

#### Hydrophobicity shifts are biochemical outcomes

The expected Kyte–Doolittle shift is a **biochemical** quantity: it says the model's
predicted residues became more or less hydrophobic at the measured positions. It is **not**
a localisation readout and **not** a signal-peptide readout. Hydrophobicity at N-terminal
positions is shared by signal peptides, transmembrane segments, transit peptides, and
ordinary hydrophobic N-termini with no targeting role at all. Reporting it as "the signal
peptide got stronger" would be an unlicensed substitution of a mechanism for a measurement.
It is reported as what it is, with its units and its sign.

#### How a directional hypothesis may be formed — and how it may not

> **A feature's intervention direction is never inferred from signal-peptide presence
> alone.** "This feature co-occurs with annotated signal peptides, therefore ablating it
> should reduce predicted hydrophobicity" is not a valid inference. Co-occurrence fixes
> neither the sign nor the existence of a downstream effect: the feature may encode the
> *boundary* of the region rather than its hydrophobic core, may fire on the cleavage-site
> motif, may be anti-correlated with the hydrophobic stretch within the annotated interval,
> or may have no signed effect on predictions at all.

A directional hypothesis is admissible only if it satisfies **all four** conditions:

1. **Justified** by one of exactly two permitted sources, named explicitly in the
   pre-registration:
   - **Discovery data** — the sign is derived from P1+P2 only. For example: within the
     discovery split, the per-position correlation between the feature's activation and
     local Kyte–Doolittle hydrophobicity, which directly estimates the sign of the
     association the intervention would perturb.
   - **Prior evidence** — a specific published finding, cited, about this or a comparable
     feature. Prior work reports signal-peptide-aligned latents at **layer 28** (L28/1541,
     L28/1470, L28/1555) from manual inspection of probe coefficients; those are a different
     layer, a different dictionary, and are not experimentally validated, so they can
     motivate a hypothesis but cannot by themselves fix a sign at layer 8.
2. **Frozen** in writing, with its sign and its justification, **before any evaluation-split
   data is touched** — in the same hashed pre-registration file as the candidate list, the
   control list, and the admissible strengths (§4.4).
3. **Labelled** in every table and every sentence of the report as either
   `direction: discovery-derived` or `direction: prior-evidence`, so a reader can see which
   kind of claim it is.
4. **Absent if neither source supports one.** If no sign can be justified, the hydrophobicity
   shift is reported **two-sided** as an exploratory, unsigned measurement, and **S3 is not
   run**. A directional test is only meaningful against a direction that was committed to in
   advance.

### 5.6 Aggregation, uncertainty, and correlated proteins

- **Unit of aggregation: the protein.** Mask positions within a protein are averaged first,
  so a long protein does not contribute more observations than a short one.
- **Uncertainty: BCa bootstrap over proteins**, 10,000 resamples, 95% CI, seed recorded.
- **Correlated proteins.** Proteins in the same homologous family are not independent
  observations, and a protein-level bootstrap would understate the interval. The honest
  treatment is a **cluster bootstrap resampling MMseqs2 clusters** (§3.3) rather than
  individual proteins. Both intervals are reported.
  - If the clustering is not run, the protein-level interval is reported **with the explicit
    statement that it assumes an independence we have not verified** — not quietly.
- **Multiplicity.** Up to 20 candidates are tested. The primary analysis is prespecified as
  the **best-performing candidate on the tuning split P3** (chosen on P3, not P4), giving one
  primary test. All 20 are reported with Benjamini–Hochberg-adjusted intervals as a
  secondary, clearly labelled, analysis.

### 5.7 Distinguishing general disruption from a specific effect

Three requirements, all operational, all measured on the evaluation split:

| # | Requirement | Operational test |
|---|---|---|
| **S1** | **Differential sensitivity to feature identity** | Mean $\Delta\mathrm{KL}_{\text{spec}} > 0$ with 95% CI excluding 0 — masked-token predictions are more sensitive to the candidate than to its magnitude-matched control. **Not** a localisation result on its own (§5.4) |
| **S2** | **Positional specificity** | (effect at `M-nterm`) − (effect at `M-distal`), same difference-against-control form, CI excluding 0 |
| **S3** | **Directional biochemical effect** | Expected-hydrophobicity shift matches the **frozen, labelled** directional hypothesis (§5.5) **and** exceeds the matched control's shift, CI excluding 0. **Not run** if no direction could be justified from discovery data or prior evidence |

A perturbation that moves every distribution by a similar amount regardless of which feature
was edited (**S1 fails**), where in the protein it was edited (**S2 fails**), and in no
consistent direction (**S3 fails**) is **general disruption**. The random-direction control
(§5.3 #3) is what makes S1 interpretable: without it, "ablating the candidate changes the
output" is satisfied by any perturbation of comparable magnitude.

### 5.8 What each possible result supports, and what stays unresolved

| Outcome | Supports | Still unresolved |
|---|---|---|
| **S1, S2, S3 all hold** | A layer-8 SAE feature to which ESM-2's masked-token predictions are differentially sensitive, position-specifically and with a biochemical direction matching a pre-frozen hypothesis. Combined with the §3 characterisation, this is **consistent with** a targeting-signal representation | Whether it is the *only* carrier; whether anything causal for the biological property was touched; whether a real protein would relocate; whether masked-token sensitivity relates to localisation at all beyond the co-occurrence evidence. **None of these is addressed** |
| **S1 holds; S2, S3 fail** | Masked-token predictions are differentially sensitive to the feature's identity, but not to where it fires or in any consistent biochemical direction — more consistent with a diffuse whole-protein correlate, i.e. the **indirect** arm of the overarching question | Whether a positional effect exists at another layer, or at a finer position resolution |
| **S1 fails** | No evidence that masked-token predictions are more sensitive to the candidate than to a magnitude-matched arbitrary direction. A **real, reportable negative**, and a direct answer to the indirect-encoding arm | Whether a different selection rule, pooling, or layer would find a used feature. The negative is about *this* rule, not about ESM-2 in general |
| **S1 holds but the absolute effect is small** | The same as S1 holding at any magnitude: the candidate's direction is used more than the matched controls. Effect size is reported, but smallness is not itself disqualifying | Whether an effect of this magnitude matters for anything downstream of masked-token prediction |
| **Candidate set fails the §4.3 stability bar** | That layer-8 N-terminal feature-label associations do not replicate across homology partitions — itself worth reporting | Whether a larger discovery split or different pooling would stabilise it |

> **Residual signal after controls does not prove direct encoding.** The control set is
> finite. One matched control feature and one random direction do not span the space of
> things correlated with the candidate — family, composition, length and kingdom all remain
> uncontrolled in the intervention itself. The strongest claim S1–S3 jointly support is:
> *the candidate's direction is used more than the specific controls we tested.* That is a
> bounded claim about a tested comparison, not a demonstration that encoding is direct.

---

## 6. Sequence variants and external evaluation

**Nothing in this section runs under this specification.** It is a proposed protocol whose
execution is gated on §7.

### 6.1 Bounded masking / sampling protocol

| Element | Specification |
|---|---|
| **Starting proteins** | **N = 50**, drawn from the P4 evaluation split with the engineering blocklist applied, stratified 25 Extracellular / 25 Cytoplasm, kingdom-stratified where counts allow, fixed seed, drawn **before** any variant is generated and by a rule that reads only the label and kingdom — never an outcome |
| **Matched across conditions** | **The same 50 proteins and the same edit positions are used in every condition.** A paired design: the only thing that varies across conditions is which intervention is active during sampling |
| **Equal edit budget** | Every condition changes **exactly B residues per protein**, $B \in \{1, 3, 5\}$, fixed in advance. Identical B in every condition — no condition gets more freedom to change the sequence than another |
| **Edit positions** | The protein's `M-nterm` positions (§5.2), taken in a fixed left-to-right order, first B of them |
| **Sampling** | Primary: greedy argmax over AA-20, deterministic. Secondary: temperature-1 sampling with 5 recorded seeds. At each position: mask it, forward pass **with the intervention active**, take the predicted residue, substitute, commit before moving to the next position |
| **Conditions** | 6 — no intervention; candidate ablate; candidate amplify; matched-control ablate; matched-control amplify; random-direction |
| **Total** | 6 conditions × 50 proteins × 3 budgets = **900 variants**, plus the 50 originals. Every variant stored with its SHA-256, condition, budget, positions and seed |

Equal budgets and shared starting proteins are what make the conditions comparable: without
them, a condition that happened to change more residues would look more effective for a
reason that has nothing to do with the feature.

### 6.2 WoLF PSORT — proposed role and its limits

**Proposed role:** a *second opinion* from a non-PLM localisation predictor on whether
variants differ from their own originals. A sanity readout. **Not** ground truth, **not** an
independent oracle, and **not** a validation of the intervention.

**What the primary documentation actually says** (fetched from `wolfpsort.hgc.jp`; the
`help` page returns 404 and could not be read):

| Item | What the site says |
|---|---|
| Organism types | Three options, labelled **"Animal"**, **"Plant"**, **"Fungi"**. One must be chosen per submission |
| Output format | Example result reads `extr: 27`, above a table headed **"27 Nearest Neighbors"**, all 27 of class `extr` |
| Output semantics | The site **never uses the words "vote", "score" or "probability"** for this number. The reading consistent with the page is a count of nearest neighbours carrying that site |
| Training / reference set | "over 12,000 animal sequences and more than 2,000 plant and fungi sequences respectively … gathered mainly from Uniprot", plus several hundred *Arabidopsis thaliana* sequences from Gene Ontology. **No version, no date** |
| Redundancy filtering | **Not mentioned anywhere on the pages we could read** |
| Accuracy claim | "Our experiments (paper in preparation) show that the overall prediction accuracy of WoLF PSORT is over 80%" |
| Download / batch / licence | **No download, source, licence, terms, or per-submission limit appears on any page we read.** Input is "From Text Area" or "From File", multi-FASTA |

**Two corrections to carry forward:**

1. The site's "paper in preparation" text is **stale**. A published method paper exists
   (Horton et al., *Nucleic Acids Research* 35:W585–W587, 2007). The per-organism k values
   — **k = 32 (animal), 27 (fungi), 14 (plant)** — come from that paper, not from the web
   pages, and should be cited to it.
2. **Do not quote the site's ">80% accuracy."** It is attached to an unpublished claim on a
   stale page, the published paper reports lower figures, and independent benchmarks report
   lower still.

**Score interpretation rule:** treat the number as a **neighbour count out of k**, not a
probability and not a calibrated score. Differences of one or two neighbours are not
meaningful. Report raw counts with k stated, never normalised into a pseudo-probability.

**Organism restriction, and how our cohort maps onto it.** Every submission must declare
Animal, Plant or Fungi. Metazoa → Animal, Viridiplantae → Plant, Fungi → Fungi are
straightforward.

The **342 "Other"** proteins in P1–P4 (**63** of them in the P4 evaluation split) need
**checking, not assuming**. "Other" is DeepLoc's residual kingdom bucket, not a statement
that the organism falls outside the predictor's three types: it collects whatever did not
fall into their three named kingdoms, and the NCBI taxonomy lineage retrieved in §3.3 is
what actually settles each case. Some of those organisms may sit inside a clade the
predictor's Animal, Plant or Fungi reference sets plausibly cover; others (bacteria,
archaea, viruses, and eukaryotes outside those three groups) will not.

The rule is therefore:

1. Resolve each "Other" protein's organism and lineage from the retrieved taxonomy.
2. Assign an organism type where the lineage **unambiguously** falls within one of the three.
3. Where it does not, or where the taxonomy could not be resolved, keep the protein as an
   **explicitly flagged record** — `organism_type = unsupported` or `unresolved` — carried
   through every table with its reason, and excluded from WoLF PSORT submission while
   remaining in every other analysis.

Dropping them silently, or forcing them into one of the three types, are both wrong. The
counts in each flag category are reported.

> **Checked.** The §3.3 fetch resolved all 342 lineages: 221 Sar, 98 Amoebozoa, 22 Discoba,
> 1 Metamonada, 1 unresolved (a deleted accession). All are eukaryotes outside Metazoa,
> Viridiplantae and Fungi, so in this instance the bucket and the support verdict coincide
> — but that is a measured result, not the assumption the rule above forbids. The 342 are
> retained in every table as flagged records.

### 6.3 Training overlap is unmeasured; independence is unestablished

Two distinct statements, and the difference between them matters:

| Statement | Status |
|---|---|
| WoLF PSORT's reference set **overlaps** our evaluation proteins | **Unmeasured.** Not demonstrated, not ruled out |
| WoLF PSORT is **independent** of our data | **Unestablished.** Not demonstrated either way |

**Shared provenance is not proof of overlap.** The about page says the reference set was
"gathered mainly from Uniprot"; our cohort is Swiss-Prot. That makes overlap *plausible* —
it is drawn from the same resource — but it does not establish that any particular protein
of ours is in their set. Their set is roughly 12,000 animal and ~2,000 each of plant and
fungal sequences against a Swiss-Prot of hundreds of thousands, it is undated, and we have
no listing of it. Asserting overlap would be as unfounded as asserting independence.

What we can state precisely:

- No version or date is published for their reference set, so we cannot bound overlap by date.
- No redundancy reduction or homology filtering *between query and reference* is documented
  on the pages we could read (the help page returns 404).
- We hold no copy of their reference set, so we cannot measure overlap directly.

**Consequence for interpretation.** Because independence is unestablished, agreement
between WoLF PSORT and our labels on *original* sequences cannot be read as independent
confirmation — it may reflect shared provenance, and we cannot tell. The predictor is used
as a second opinion whose independence is an open question, and every number from it is
reported with that caveat attached.

**Paired differences reduce one concern; they do not eliminate it.** The paired
within-protein comparison between a variant and its own original, under the same organism
type, removes a *constant* per-protein offset — including a constant lookup effect, if one
exists. It does **not** eliminate the overlap problem, and it does not remove predictor
bias:

- If a protein (or a close homolog) is in the reference set, its neighbour structure is
  anchored by that entry, and the *sensitivity* of the paired difference to small edits is
  itself shaped by that anchoring. A cancelled offset is not a cancelled dependency.
- The predictor's features and k-nearest-neighbour behaviour carry their own biases —
  toward N-terminal composition in particular, which is exactly what our variants perturb.
  A paired difference inherits those biases in full.
- Changing 1–5 residues may leave the nearest-neighbour set entirely unchanged, in which
  case the paired difference is *structurally* zero and carries no information about the
  intervention. This must be reported as a count ("in X of 50 proteins the neighbour set
  was unchanged"), never averaged into a mean that hides it.

**What could be measured, if this arm is ever run.** A partial bound is obtainable without
their dataset: report how many evaluation proteins are Swiss-Prot entries whose first-public
date precedes the plausible range of their snapshot, and report the kingdom composition
against their stated reference-set sizes. That bounds *exposure*, not overlap, and should be
labelled that way. The entry dates needed for it are collected in §3.3.

**The reportable statement** is the unembellished one: *we did not measure overlap between
our evaluation proteins and WoLF PSORT's reference set, we could not, and we therefore treat
its independence as unestablished rather than assuming it in either direction.*

### 6.4 Biological-concept annotations stay separate from external evaluation

A strict one-way wall:

- UniProt signal-peptide, transit-peptide and family annotations (§3) are used **only** to
  characterise features and to interpret §5 results. They are never used as inputs to WoLF
  PSORT, and never as a validation target for it.
- WoLF PSORT output is **never** used to select candidate features, choose strengths, pick
  proteins, or define any outcome in §5.

Without this wall, a concept annotation and an external predictor trained on the same
curated resource would mutually confirm each other with no new evidence entering.

---

## 7. Feasibility and stopping points

### 7.1 Hardware and measured runtimes

**Local hardware:** Apple M4 Pro, 24 GB unified memory, MPS backend, torch 2.8.0,
transformers 4.57.6, ESM-2 650M at revision `08e4846e537177426273712802403f7ba8261b6c`.

**Measured, not estimated:**

| Job | Scope | Measured cost |
|---|---|---|
| Layer-8 extraction, full P1–P4 | 10,139 proteins, 4,819,401 visible residues, fp16 | **2,469.6 s = 41.2 min**, 4.11 seq/s |
| Intervention validation | 8 proteins, 12 conditions, 3 repeat runs, fp32 | **32.3 s** |
| Edit-implementation comparison | 96 conditions with float64 references on CPU | **51.0 s** |
| SAE pilot | 64 proteins per dataset, reference-equivalence checks | within the pilot report |

### 7.2 Derived pilot costs, with the basis stated

| Stage | Basis | Derived cost |
|---|---|---|
| **Discovery pass** (P1+P2 per-residue activations → pooled summaries) | P1+P2 hold 2,439,159 visible residues = 50.6% of the cohort's 4,819,401; the full cohort took 41.2 min at fp16 | **≈ 21 min** for the ESM forward passes, **plus** the SAE encode (1280×4096 per residue, CPU), which is **not** reliably extrapolable from the validation's 32.3 s mixed-workload figure |
| → required first step | — | **A timed 200-protein pilot**, projecting the SAE-encode cost from measurement rather than from this table |
| Discovery storage | 4,096 features × 5 summaries × 5,201 proteins, fp32 | **≈ 426 MB** (213 MB fp16). Per-residue activations are **not** stored — that would be ≈ 9 GB |
| **Intervention, full P4** | 2,378 proteins × 6 conditions × 6 mask positions ≈ **85,600** masked forward passes, at ≥ 0.24 s each | **≈ 6–7 hours. This is the binding constraint.** |
| **Intervention, reduced P4** | A prespecified random **400** proteins from P4 (stratified by class, fixed seed, **drawn before any intervention runs**): 400 × 6 × 6 ≈ 14,400 passes | **≈ 1.2 h.** This is the minimum experiment; full P4 is an optional extension |
| **Variant generation** (§6) | 900 variants × up to 5 sequential masked passes ≈ 4,500 passes | **≈ 25 min** greedy; ×5 for the sampling seeds |
| **UniProt annotation fetch** (§3.3) | 10,139 accessions, rate-limited HTTP, no GPU | **< 1 h**, dominated by the API |
| **MMseqs2 clustering** (§3.3, second priority) | 10,139 sequences, CPU | **minutes** |

### 7.3 The minimum complete scientific experiment

Six steps. All six must run for the result to be reportable at all:

1. UniProt annotation fetch with SHA-256 verification, and the three-state signal-peptide
   coverage table (§3.3).
2. Timed 200-protein pilot, then the discovery pass on P1+P2, storing pooled summaries only
   (§4.1).
3. Candidate and control selection, with the P1-vs-P2 stability check; lists and blocklist
   written and hashed (§4.2–§4.5).
4. Strength admissibility on P3; admissible set frozen in writing (§5.1).
5. Intervention on the prespecified 400-protein P4 subset: candidate, matched control,
   random direction, no-op; `E-nterm` and `E-all` × `M-nterm` and `M-distal`; primary
   outcome ΔKL_spec plus the S1/S2/S3 specificity tests, with the reconstruction-replacement
   diagnostic reported alongside (§5).
6. The report — including, in full, the negative-result path.

Total compute on the minimum path: **roughly 2–3 hours** of local GPU-equivalent time, plus
the annotation fetch.

### 7.4 Optional extensions, in priority order

1. Full P4 instead of the 400-protein subset (≈ 6 extra hours).
2. MMseqs2 clustering and the cluster bootstrap (§5.6) — cheap, and it closes a real gap.
3. Sequence variants (§6.1).
4. WoLF PSORT evaluation (§6.2) — only after the §7.5 gate.
5. A layer-24 cross-check of any surviving candidate.
6. A layer sweep. Lowest priority, and only under the §4.5 condition.

### 7.5 Conditions for proceeding, revising, or stopping

**Proceed to intervention if, and only if, all three hold:**
- ≥ 1 candidate feature meets the stability bar (≥ 10 of 20 in both single-partition lists);
- that candidate has a matched control satisfying all three §4.4 tolerances;
- ≥ 1 strength on P3 is admissible under both §5.1 criteria.

**Proceed to sequence variants and WoLF PSORT if, and only if:**
- S1 holds on the evaluation split (CI for ΔKL_spec excludes 0).

Otherwise there is no established differential effect to generate variants from, and
generating them anyway would produce 900 sequences that differ for no established reason.
Note that this gate is on S1 — a *matched-control* comparison — not on absolute effect size
and not on any comparison with reconstruction replacement.

**Revise if:**
- the candidate set is partition-unstable. Report the instability as the finding; the
  permissible revisions are limited to pre-named alternatives (whole-sequence pooling
  instead of N-terminal, already specified in §4.1). Any revision is documented as a
  deviation **with the original rule's result reported alongside**.

**Stop if:**
- the primary outcome's CI includes 0. Report the negative. **Do not** escalate strengths,
  add candidates, switch pooling, or change layers in search of an effect;
- the annotation fetch shows signal-peptide coverage so low that the characterisation arm
  cannot run — in which case report the coverage, and reduce §1.2 stage 2 rather than
  substituting predictor output for annotation.

### 7.6 No requirement that interventions improve any score

Stated plainly, because the incentive runs the other way:

> **Success is a clear verdict, not a positive effect.** A null — "ablating the most
> label-associated layer-8 SAE features changes ESM-2's masked-residue predictions no more
> than a magnitude-matched arbitrary direction does" — is a direct answer to the overarching
> question and evidence for the **indirect**-encoding arm. It is as publishable within this
> project as a positive, and the design is symmetric by construction: every specificity test
> is a difference against a control that can come out at zero.

Nothing in this specification rewards a larger effect. Strengths are frozen before
evaluation, candidates are frozen before evaluation, and the stopping rules forbid searching
for an effect after a null.

---

## 8. Summary tables

### 8.1 The experiment, in one table

| | |
|---|---|
| **Inputs** | Frozen ESM-2 650M (`facebook/esm2_t33_650M_UR50D`, rev `08e4846e…`), eval mode, no gradients; InterProt layer-8 SAE (d=1280, 4,096 features, k=128); Swiss-Prot binary cohort, 10,139 proteins in P1–P4, sequence only |
| **Selection data** | **P1 + P2 only** (5,201; 1,339 Extracellular / 3,862 Cytoplasm). Top-20 features by \|AUROC−0.5\| of N-terminal pooled activation, among features active in ≥5% of discovery proteins; stability checked P1 vs P2. Strengths admitted on **P3** (2,556). **P4 and P0 unread at this stage** |
| **Intervention** | Direct decoder-delta at `hidden_states[8]`: $h' = h + ((\mathrm{topK}(a_{\text{ed}}) - \mathrm{topK}(a))\,W_{\text{dec}})\cdot\sigma$; candidate feature scaled by $s \in \{0, 0.5, 2, 4\}$ at edit positions (`E-nterm` = 1–30, or `E-all`) |
| **Controls** | (1) no-op $s{=}1$, must be exactly zero; (2) matched control feature — activity and magnitude matched ±20%, \|AUROC−0.5\|≤0.02; (3) random direction at matched L2 norm. **Specificity is assessed against (2) and (3) only.** (4) reconstruction replacement is a **diagnostic reference** on the SAE parameterisation (mean KL 0.293 nats, top-1 flipped 2/8 in validation), reported but never used as a threshold |
| **Evaluation data** | **P4 only**, minus the 8 engineering accessions (2,378 available), reduced to a prespecified seeded random 400 stratified by class. **P0 reserved. HPA excluded** |
| **Measurement** | **Primary:** ΔKL_spec = mean over proteins of [mean KL(clean ‖ intervened) over AA-20 renormalised, in nats, at `M-nterm` mask positions, candidate ablation] − [same for matched control], BCa bootstrap 95% CI over proteins and over homology clusters. **Specificity:** S1 differential sensitivity to feature identity, S2 positional, S3 directional biochemical effect against a frozen, labelled hypothesis (omitted if no direction could be justified). **ΔKL_spec measures differential sensitivity of masked-token predictions, not localisation specificity.** **Secondary:** top-1 change rate, total variation, expected-hydrophobicity shift (a biochemical outcome), probe-score shift (manipulation check only), realised ‖Δh‖₂, eviction counts |
| **Permissible conclusion** | At most: *a layer-8 SAE feature exists to which the frozen model's masked-token predictions are differentially sensitive — more so than to the specific matched controls tested — position-specifically, with a biochemical direction matching a pre-frozen hypothesis, and which co-occurs with annotated targeting signals.* **Never:** that a protein's biological localisation changed; that encoding is direct in general; that the effect is causal for the biological property; that a KL difference by itself shows localisation specificity; that a hydrophobicity shift shows a targeting signal strengthened; that WoLF PSORT independently confirmed anything |

### 8.2 Check-in 1 concerns → design choices

| # | Check-in 1 concern | Where raised | Design choice here |
|---|---|---|---|
| 1 | **Provenance leakage** (solubility dataset header text → 0.7792 pooled AUROC), and SAE features inherit it | Risks; SAE pilot §6 | This specification is **DeepLoc-only**. No pooled solubility claim is made, so the confound is out of scope rather than carried in silently |
| 2 | **Incomplete family annotations** (49/100 on the *other* dataset; **0/10,139 here**) | Risks; Gap 3 | §3.3 fetches Pfam/InterPro for this cohort in the minimum work; §3.1D states the coverage is currently zero; the denominator is reported, never imputed. Family is a **characterisation covariate**, not a control we claim to hold fixed |
| 3 | **Unverified homology overlap** in both datasets | Risks; Gap 2; Q5 | §2.6 **refuses** the family-disjointness claim outright; §3.3 adds MMseqs2 clustering to *measure* it; §5.6 specifies a cluster-level bootstrap and requires the independence assumption to be stated when it is not run |
| 4 | **Probes are correlational; only intervention tests causation** | Risks; Q4 | The whole of §5. §1.3 separates Readout A (our classifier, circular) from Readout B (native ESM-2); §1.4 forbids the relocation claim explicitly |
| 5 | **No external test set for the DeepLoc contrast; P0 is internal** | Risks; Gap 1; Q6 | §2.4 states plainly that P1–P4 are already used and P4 is not a fresh test; P0 stays reserved as the one untouched split; HPA stays excluded; §6.3 refuses to treat WoLF PSORT as the missing external check |
| 6 | **The gain cannot be attributed to the N-terminal embedding specifically** — positions 1–30 attend to the whole sequence | Check-in 1 §5; Q4 | §4.1 repeats the caveat for SAE summaries and declines to make a positional claim from pooling. The positional claim is tested by **intervention position schemes** (`E-nterm` vs `E-all`, `M-nterm` vs `M-distal`) in §5.2, which is the thing probing could not decide |
| 7 | **"Isn't this just rediscovering signal peptides?"** | Q4 | §5.7's S2 and S3 — positional and signed-directional specificity — are exactly the tests that distinguish "a feature co-occurs with the N-terminus" from "the model's computation depends on it there, in that direction" |
| 8 | **Layer 8 vs the layer-28 latents cited from prior work** | Q2 | §4.5 keeps layer 8 with three stated reasons, makes a sweep **optional** with a named trigger condition, and prefers a single layer-24 cross-check over a sweep if one is ever made |
| 9 | **338/4,096 SAE latents dead at train end; the layer-8 checkpoint is released without the layer-24 checkpoint's documentation** | SAE pilot §6 | §4.2 filters eligibility by **measured activity on our own cohort**, never by the training-time dead list — at least 252 of those 338 fire on our proteins. §4.5 names the layer-24 cross-check as the highest-value optional extension |
| 10 | **The SAE captures only ~40–48% of residue-to-residue variance at k=128** — the variation that would carry positional information is largely not captured | SAE pilot §6 | §5.3 control 4 turns this from a caveat into a **measured diagnostic** on the dictionary, reported beside every effect. It is explicitly *not* a specificity threshold: specificity is judged against the matched-feature and matched-norm interventions only |
| 11 | **Pooled mean-vector claims said nothing about per-residue encoding** | SAE pilot §6 | §4.1 defines five explicit pooled summaries and states what each is and is not; the per-residue question is answered by the intervention, not by the pooling |

---

## 9. The single next executable task

> **Fetch UniProt annotations for the 10,139 P1–P4 cohort accessions and emit a verified
> coverage table.**

Concretely: one script that queries the UniProt REST API for `accession, reviewed,
organism_id, sequence, length, ft_signal, ft_transit, ft_chain, ft_propep,
protein_families, xref_pfam, xref_interpro, cc_subcellular_location`; verifies
`sha256(returned_sequence)` against
[`DeepLoc2/baseline/deeploc_binary_manifest.tsv`](DeepLoc2/baseline/deeploc_binary_manifest.tsv)
for every entry; and writes one TSV plus a coverage table giving, broken down by partition,
class and kingdom:

- the **three-state** signal-peptide variable — `annotated_present` (with ECO evidence
  class), `not_annotated`, `unverifiable` (sequence hash mismatch) — never collapsed to a
  binary;
- transit-peptide coverage, on the same three-state basis;
- family/domain annotation coverage, with the denominator explicit.

**Reads P1–P4 accessions only. Does not query P0.** No model runs, no SAE, no intervention,
no sequence submitted anywhere except as accession lookups against UniProt.

**Why this one first.** Every other stage depends on it, and it is the only step whose
*outcome can change the design*: if signal-peptide coverage turns out to be very low, the
characterisation arm of §1.2 stage 2 shrinks and the experiment should be re-scoped before
any compute is spent. Spending 2–3 hours of GPU time before knowing that would be the wrong
order.

Note what this step does **not** decide. The directional hypothesis behind S3 is **not**
derived from signal-peptide annotations (§5.5): annotation coverage affects the
characterisation arm, not the sign of any intervention effect. A direction, if one is
claimed at all, comes from discovery-split statistics or cited prior evidence, and is frozen
before evaluation independently of what this fetch returns.

---

## Amendment A — layer-wise decodability is separated from intervention effects

**Date:** 2026-11-22. **Status:** in force. **Supersedes:** the launch condition in §4.5 only.
**Preserves:** every other clause of §4.5, and all of §1.3, §1.4, §5 and §6.

### A.1 What was amended, and the original wording

§4.5 stated:

> **A layer sweep is optional, not necessary**, and is not launched by this specification. […]
> The condition under which it becomes warranted is narrow: **only if the layer-8 experiment
> yields a candidate with a reproducible, specific effect**, at which point "is the effect
> layer-local?" becomes a real question.

That condition is **withdrawn for a probe sweep** and **retained for an intervention sweep**.
The original §4.5 text stays in the document, unedited, under a banner pointing here.

### A.2 Why — the condition conflated two different measurements

§1.3 already insists that two readouts must never be conflated. The original §4.5 condition
violated that principle in a third place: it gated a **decodability** question behind an
**intervention** result.

| | **Layer-wise decodability** (this amendment) | **Intervention effect** (§5, unchanged) |
|---|---|---|
| Question | At which layer is the Extracellular/Cytoplasm distinction most linearly readable from the pooled representation? | Does editing an SAE feature change ESM-2's own masked-token predictions? |
| Needs an SAE | **No.** A linear probe on pooled hidden states. | **Yes** — one validated checkpoint per layer |
| Needs a candidate feature | **No** | **Yes** |
| Cost | One forward pass, reused across all 33 layers | A separate SAE checkpoint and reference-agreement validation per layer |
| What a null result means | The information is not linearly decodable there under this pooling | The edit did not propagate |

The original condition made sense for an **intervention** sweep, whose per-layer cost is a
validated SAE checkpoint. It does not transfer to a **probe** sweep, whose marginal cost over
the existing layer-8 extraction is one forward pass with `output_hidden_states=True`. A probe
sweep is also informative precisely when the layer-8 intervention is null: "the effect was
null at layer 8" and "layer 8 is where the signal is most decodable" are separate facts, and
learning the second does not depend on the first.

**An intervention sweep remains conditional exactly as §4.5 wrote it.** Nothing here authorises
editing at any layer other than 8.

### A.3 The boundary this amendment must not be read across

A probe sweep locates where a linear readout works best. It does **not** locate where the model
"stores" localisation, and it does **not** nominate an intervention layer.

> **The best probe layer is not thereby the best intervention layer, and high decodability is
> not evidence of direct encoding.** A linear probe recovers what is linearly available in a
> pooled vector. That is consistent with a direct representation, with a diffuse correlate of
> amino-acid composition, and with the model having encoded something that merely co-varies
> with localisation, such as protein family. §1.4's limits apply unchanged, and the probe
> sweep adds no new licence.

### A.4 What the sweep may and may not be called

- It is **not a replication** of the Check-in 1 four-fold layer-8 result. That used four-fold
  cross-validation over P1–P4 with each partition held out once; this uses a single fit on
  P1+P2 with P3 for comparison and P4 for evaluation. Different fitting and evaluation design,
  different numbers, not comparable as a reproduction.
- P4 is **not a fresh final test set**. §2.4 already established that P1–P4 are not a
  previously untouched final test; P4 specifically has carried the 400-protein intervention
  evaluation. Any P4 figure here is a *development* evaluation and is labelled as such.
- The operative protocol is frozen in
  [`layer_sweep/LAYER_SWEEP_PROTOCOL.md`](layer_sweep/LAYER_SWEEP_PROTOCOL.md) and
  [`layer_sweep/v1/frozen_protocol.json`](layer_sweep/v1/frozen_protocol.json), both written
  before any sweep outcome was computed.

### A.5 Consequential update to §2.6 and §5.6

§2.6 recorded that cluster identifiers and cross-partition homology were unavailable. They are
now measured — [`HOMOLOGY_AUDIT.md`](HOMOLOGY_AUDIT.md). Accordingly:

- §5.6's fallback ("if the clustering is not run, the protein-level interval is reported with
  the explicit statement that it assumes an independence we have not verified") **no longer
  applies.** The cluster bootstrap is runnable and is used.
- The clusters are an **operational grouping**. Resampling them makes an interval account for
  that grouping; it does **not** establish that the clusters are independent of one another.
  Cluster counts are reported as cluster counts, never as effective sample sizes.
