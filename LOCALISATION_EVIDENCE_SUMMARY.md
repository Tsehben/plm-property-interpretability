# Localisation interpretability — evidence summary

Prepared from saved outputs only. No model was run, no variant generated, no
experimental setting changed, P0 never read, nothing committed or pushed. All
original results are preserved unchanged; §2 marks **corrections to
interpretation** and §§3–4 mark **new retrospective analyses**, which were not
pre-registered.

> **Correction log — 2026-11-22, second pass.** Five interpretation statements in
> this document were corrected. No number changed. (1) §1 and §6: feature 1446 is
> described as *strongly associated* with signal-peptide annotation, including in
> the label-discordant groups — the earlier "tracks the annotation, not the label"
> claimed exclusivity and implied an identification. (2) §1 and §6: the claim that
> the candidates are "redundant" and "not 20 independent mechanisms" is replaced
> by the measured overall and within-class correlations. (3) §3: the intervals are
> stated as conditional on an unverified independence assumption, with no asserted
> direction of interval error — the earlier text said dependence "generally makes
> such intervals too narrow." (4) §4: identical sequences across paired arms do not
> imply identical token distributions. (5) §6: the proposed divergence-conditioned
> follow-up sample is **withdrawn** as outcome-based conditioning.

---

## 1. The three completed stages

### Stage 1 — Discovery associations (P1 + P2)

| | |
|---|---|
| **Question** | Does any single SAE feature at ESM-2 layer 8 separate extracellular from cytoplasmic proteins from its activity over the first 30 residues? |
| **Comparison** | For each of 2,867 eligible SAE features (of 4,096), AUROC of the mean first-30 activation against the frozen DeepLoc binary label; ranked by \|AUROC − 0.5\|; top 20 distinct feature IDs retained. |
| **Observed** | Strong separations exist. Feature **1446**: AUROC **0.974** (P1 0.965, P2 0.979), higher in Extracellular; **2564** 0.957; **3900** 0.023, i.e. 0.977 toward Cytoplasm. 1446 holds across kingdoms (0.948–0.985) and length quartiles (0.962–0.982). Across the 190 candidate pairs: per-protein Spearman correlation **median 0.697 over the whole sample** (range −0.725 to 0.833; 91 pairs > 0.7, 0 pairs > 0.9), but **median 0.139 within Extracellular and 0.179 within Cytoplasm**; activation-set Jaccard median 0.536 (range 0.233–0.714); decoder-direction cosine **median 0.008** (range −0.098 to 0.305, 0 pairs > 0.5). |
| **Population represented** | 5,190 P1+P2 proteins — **1,336 Extracellular, 3,854 Cytoplasm** — after excluding 8 engineering proteins and 27 P0-collision records. Not a random sample of proteins; the DeepLoc cohort with its own inclusion history. |
| **Strongest supported conclusion** | Single layer-8 SAE features separate the two classes well on held-out-from-nothing screening data. Across the 190 candidate pairs, **per-protein behaviour is moderately correlated over the whole sample (median Spearman 0.697, range −0.725 to 0.833; 91 pairs above 0.7, none above 0.9) and weakly correlated within class (median 0.139 within Extracellular, 0.179 within Cytoplasm), while their decoder directions are close to orthogonal (median cosine 0.008, max 0.305)**. Nothing here is causal, and no significance is claimed from screening thousands of features. |

A further discovery result bears directly on what 1446 *is*. Conditioned on
signal-peptide annotation:

| 1446 activation | n | active | mean |
|---|---|---|---|
| Extracellular, signal annotated | 1,252 | 99.8% | 4.48 |
| **Cytoplasm, signal annotated** | **59** | **98.3%** | **4.45** |
| Extracellular, not annotated | 53 | 50.9% | 0.95 |
| Cytoplasm, not annotated | 3,645 | 21.5% | 0.09 |

**1446 is strongly associated with signal-peptide annotation, including in the
observed label-discordant groups.** The 59 cytoplasmic proteins carrying a signal
annotation activate it about as strongly as the 1,252 extracellular ones
(98.3% vs 99.8% active; mean 4.45 vs 4.48), and the 53 extracellular proteins
without an annotation activate it much less (50.9% active, mean 0.95).

Three limits on that sentence. The association is **not exclusive**: 21.5% of the
3,645 unannotated cytoplasmic proteins still show some activity, and 49.1% of
unannotated extracellular proteins show none, so annotation status does not
partition activation cleanly. It is **not a mechanistic identification**: nothing
here shows that 1446 detects a signal peptide, represents one, or is used by the
model as one. And `not_annotated` is **not confirmed biological absence** — it
means no signal-peptide feature was recorded for the frozen sequence.

### Stage 2 — P4 masked-token interventions

| | |
|---|---|
| **Question** | Does editing a candidate feature at a masked N-terminal residue change the model's own masked-token distribution more than editing norm-matched control features does? |
| **Comparison** | Paired KL difference: candidate KL(clean‖edited) minus the mean of available activation-matched, residual-norm-matched control features, at the same protein, position and multiplier. Four primary contrasts (2564, 3900 × multipliers 0, 2), protein-level stratified BCa, 10,000 resamples, Bonferroni 98.75%. |
| **Observed** | **One of four excludes zero and survives Bonferroni: 2564 ablation, +2.10e-05 nats, 95% [+1.13e-05, +3.50e-05], 98.75% [+9.05e-06, +4.01e-05].** Levels on the same paired cells: candidate 4.80e-05 vs control 2.70e-05. The other three include zero; for both of 3900's contrasts the *median* protein's controls move the distribution more than the candidate does. Integrity complete: 5,895 no-op conditions at exactly zero KL, 0 unintended top-K changes, control norm match max relative error 6.79e-07. |
| **Population represented** | 400 P4 proteins, 200/200 by design — **not natural prevalence**. The significant contrast describes a much narrower population: **61 of 400 proteins, 55 of them Extracellular, 50 of 61 paired cells at residue 20, 18.4% of candidate-active cells.** That loss is structural — matching an ablation norm needs a large negative control latent change, which removes exactly the cells where the candidate's latent was largest. P4 had been used in Check-in 1, so it is not an untouched test set. |
| **Strongest supported conclusion** | **Ablating 2564 at a masked N-terminal residue changes the masked-token distribution more than norm-matched control features do, conditional on a matched control being constructible.** The mean is carried by a right tail (median +2.3e-06; 37 of 61 proteins positive). Nothing here concerns localisation: the readout is the model's own token distribution at the edited residue in the same forward pass. **1446 was descriptive only — it had 0 of 5 eligible activation-matched controls.** |

### Stage 3 — P3 sequence-variant pilot

| | |
|---|---|
| **Question** | Does intervening on feature 1446 during masked-residue resampling change where a local predictor places the resulting sequence, relative to unmodified sampling and to magnitude-matched random directions? |
| **Comparison** | `margin(s) = score(extr) − score(cyto)` from local WoLF PSORT (animal model), `D = margin(variant) − margin(original)`; five arms (A unmodified, B 1446 ×2, C 1446 ×0, D random +u, E random −u) on 10 fixed positions, 4 replicates, shared uniforms; contrasts B−A, C−A, B−D, C−E formed within replicate, averaged within protein. |
| **Observed** | **All four contrasts have a median paired difference of exactly zero**, with 26–30 of 40 proteins exactly tied. Means +0.0375 to +0.1063 against protein-level sd 0.64–0.75. Arm A alone shifts the margin by +0.93. Integrity complete: 8,000/8,000 passes, arm A exactly a no-op, 0 unintended top-K changes, 0 latent coordinates touched in the random arms, 840/840 sequences scored with none lost. |
| **Population represented** | 40 Metazoa P3 proteins, 20/20 by design, length 60–1022, standard residues. **5 of the 40 (12.5%) are members of WoLF PSORT's own reference set** and are partly their own nearest neighbours. Metazoa only, so one organism model applies throughout. |
| **Strongest supported conclusion** | **At this budget, intervening on 1446 does not detectably move this predictor's extracellular-minus-cytoplasmic margin, relative to unmodified sampling or to magnitude-matched isotropic controls.** §4 shows why: in 83–87% of paired cells the two arms generated the *same sequence*. The pilot mainly establishes that the generation, control and scoring machinery works end to end and is clean. |

A method correction found during this stage and verified against the
implementation: **a class omitted from `runWolfPsortSummary`'s printed line does
not mean zero.** All 840 records omitted at least one non-zero class (mean 7.6,
max 15); `extr` was non-zero but omitted in 119 records and `cyto` in 89. Scores
are now recomputed for every class from the neighbour list, and that
recomputation reproduces **all 2,650** printed scores exactly.

---

## 2. Corrections to the sequence-pilot interpretation

These correct wording in
[sequence_variant_pilot/PILOT_RESULTS.md](sequence_variant_pilot/PILOT_RESULTS.md)
§8. **No result changes.**

**2.1 — What the random-direction comparisons actually test.**

> *Was:* "B − D and C − E answer 'does perturbation magnitude alone explain
> it?'"

That claims too much. B − D and C − E compare 1446's edit against **one specified
control family**: isotropic directions drawn uniformly on the unit sphere in
ℝ¹²⁸⁰, injected at the same token with the same realised norm. They do **not**
test every magnitude-only explanation. Untested magnitude-matched alternatives
include other SAE decoder directions, directions inside the model's high-variance
subspace, and any structured direction correlated with the residual stream — all
of which an isotropic draw is near-orthogonal to by construction (observed mean
cosine with the candidate's own Δh: 0.0016, max \|cos\| 0.070). The correct
statement is: *1446's edit was not distinguishable from magnitude-matched
isotropic random directions.* The activation-matched SAE control family, which
would test the sharper question, remains unavailable for 1446 (0 of 5 eligible).

**2.2 — Common random numbers retain experimental pairing.**

> *Was:* "Common random numbers are a coupling, not exact pairing… This is
> variance reduction only."

That understates the design. **The experimental pairing is retained throughout.**
Every contrast compares arms at the same protein, the same replicate index and
the same seed stream, and the uniform consumed at each (protein, replicate, step)
is identical across arms by construction. The paired-difference analysis is valid
as a paired analysis for the full ten-step trajectory. What divergence changes is
narrower: once the arms' sequences differ, the uniform is inverted against
*different conditional distributions*, so a step's residue difference is no
longer attributable to the distribution shift alone — it is partly inherited from
earlier divergence. Only step 1 is a clean single-step attribution. Pairing is
intact; single-step causal attribution is not.

**2.3 — The ordinary-resampling shift is context, not a threshold.**

> *Was:* "The dominant effect is masked resampling itself, not the
> intervention… the intervention effect is small against the sampling effect it
> sits in."

Arm A's +0.93 mean shift describes what masked resampling does to this margin on
these proteins. It is **not a benchmark the intervention has to exceed**, and the
contrasts are not diminished by being smaller than it. B − A, C − A, B − D and
C − E are *paired differences* in which the common resampling movement cancels;
that is exactly what pairing is for. Arm A's shift is useful context — it says the
readout is not static under neutral resampling, so a comparison against the
original alone would be badly confounded — but it sets no specificity threshold
and no effect size is judged against it.

---

## 3. NEW RETROSPECTIVE ANALYSIS — uncertainty for the four pilot contrasts

**Retrospective and exploratory. Not pre-registered.** The frozen P3 protocol
specified no uncertainty quantification; this was decided after the results were
seen. These intervals carry **no pre-registered error control**, no multiplicity
adjustment, and no confirmatory reading. Script:
[retrospective_analysis.py](sequence_variant_pilot/retrospective_analysis.py);
outputs in [v1_retrospective/](sequence_variant_pilot/v1_retrospective/).

| contrast | n | mean | median of paired differences | 95% interval | excludes 0 |
|---|---|---|---|---|---|
| B − A | 40 | +0.0375 | 0.0000 | [−0.277, +0.166] | no |
| C − A | 40 | +0.0375 | 0.0000 | [−0.263, +0.222] | no |
| B − D | 40 | +0.1000 | 0.0000 | [−0.191, +0.278] | no |
| C − E | 40 | +0.1063 | 0.0000 | [−0.031, +0.444] | no |

All four are genuine BCa; the pre-registered degenerate fallback was not needed.

**Resampling assumptions.** The unit is the protein (n = 40), replicates averaged
within protein first. Proteins are resampled with replacement, stratified by
localisation class (20/20), 10,000 resamples, seed 20261122; BCa with jackknife
acceleration. Replicates, steps and arms are **not** independent sample units.

**Unmeasured homology dependence.** The procedure assumes proteins are
exchangeable within class and independent. **They are not known to be**, and that
assumption is unverified: no homology clustering exists for this cohort, and 5 of
the 40 proteins are members of the predictor's own reference set. **These
intervals are therefore conditional on an unverified independence assumption.**
The dependence has not been measured, so no direction and no magnitude of
interval error is asserted here — not that they are too narrow, not that they are
too wide — and no effective sample size is claimed.

**An interval containing zero is not equivalence.** No equivalence bound was
pre-specified and none is claimed. These four intervals say only that this sample
at this width does not resolve the sign. C − E runs from −0.03 to +0.44 — it is
consistent with no effect and equally consistent with an effect several times its
own point estimate.

---

## 4. NEW RETROSPECTIVE ANALYSIS — what the ties are made of

**Retrospective and exploratory. Not pre-registered.** Computed from the saved
variant sequences and saved scores.

Three distinct things were being conflated. Per contrast, over 160 paired cells:

| contrast | tied cells | **tie: arms generated the SAME sequence** | **tie: different sequences, same score** | not tied |
|---|---|---|---|---|
| B − A | 142 (88.8%) | **133 (93.7% of ties)** | 9 (6.3%) | 18 |
| C − A | 144 (90.0%) | **139 (96.5%)** | 5 (3.5%) | 16 |
| B − D | 141 (88.1%) | **132 (93.6%)** | 9 (6.4%) | 19 |
| C − E | 144 (90.0%) | **137 (95.1%)** | 7 (4.9%) | 16 |

No cell had the same sequence with different scores, as it must not.

**Ties are overwhelmingly unchanged generation, not evaluator insensitivity.**
In 82.5–86.9% of paired cells the intervention did not alter a single sampled
residue: the shared uniform plus a perturbation too small to move the inverse-CDF
bin produced literally identical sequences. When the arms *did* diverge, they did
so by a median of **1 residue** (mean 1.64–2.05, max 5) out of ~7.2 substitutions.

**An identical sequence across paired arms does not mean the arms had identical
token distributions.** The edit was applied and verified on every active step: at
each such position the arms' AA-20 distributions genuinely differed. An identical
sequence means only that the *same uniform fell in the same inverse-CDF bin* in
both arms — the distributions shifted, but not far enough to move the sampled
residue across a bin boundary. So "unchanged generation" describes the realised
sequence, not the model's state: it is a statement about the sampling step's
resolution, not evidence that the intervention did nothing to the distribution.
The per-step distributions are saved and remain available for direct comparison.

**Evaluator insensitivity is real but secondary.** Among the cells where the
sequences genuinely differ, **24–33% still receive an identical predictor score**.
Across the whole run the predictor's resolution is visibly coarse: 800 variants
collapse to **210 distinct sequences**, a protein has on average 5.25 distinct
variant sequences but only **3.92 distinct margins**, and **15.4%** of
distinct-sequence pairs within a protein share a margin exactly.

**Both mechanisms are present; unchanged generation dominates.** The pilot's zero
medians therefore say more about the intervention's effect on *sampling* than
about the predictor's sensitivity — though the predictor is also coarse enough
that it would blunt a modest real effect. 35 of 800 variants were identical to
their own original.

---

## 5. Mapping to the Check-in 1 concerns

| Check-in 1 concern | What the completed work does | What remains unanswered |
|---|---|---|
| **"We have shown *that* ESM-2 helps; we have not shown *which biological features* explain the gain."** | Discovery identified specific layer-8 SAE features with strong class separation (1446 at AUROC 0.974), stable across kingdom and length strata. | **Unanswered.** A feature that separates classes is not a feature that explains the probe's gain. Nothing links these candidates to the Check-in 1 probe's +0.018 AUROC. |
| **"Comparison 2 adds the whole-input embedding *and* the first-30 embedding at the same time… leaves open which positions account for it."** | Discovery screened on a first-30 summary only, and interventions were applied at single N-terminal residues. | **Partly addressed, not resolved.** The screening window is N-terminal, but ESM-2 is bidirectional, so a first-30 activation still carries whole-sequence context. Position attribution is not established. |
| **"Our N-terminal embedding is pooled from positions that attend to the whole sequence."** | The ≤1022 stratification changed candidate AUROCs by ≤0.0051; length-quartile AUROCs for 1446 range 0.962–0.982. | **Still open.** Robustness to truncation and length is not the same as locality of information. |
| **"The next test: does an internal ESM-2 feature overlap with annotated signal peptides?"** | **Answered, and informatively.** 1446 is strongly associated with signal-peptide annotation, including in the 59 label-discordant cytoplasmic proteins that carry one (98.3% active, mean 4.45). | The association is **not exclusive** (21.5% of unannotated cytoplasmic proteins still show activity), and the overlap is with the **annotation**, much of which is itself computational (1,873 of 1,875 signal evidence rows carry no named source). Shared sequence-based bias between ESM-2 and the annotation tools is not excluded, and no mechanism is identified. |
| **"Overlap alone does not establish what ESM-2 encodes… testing causation requires intervention."** | Two interventions were run. P4: one of four contrasts shows differential model sensitivity, for 2564, narrowly scoped. P3: 1446's intervention did not detectably move an external predictor's margin. | **The causal question is not settled.** P4 shows an internal effect on the model's own token distribution, not on localisation. P3 shows no external effect at this budget, with 85% of cells not changing sequence at all. |
| **"The gain varied considerably across folds."** | Discovery stability was reported without a pass quota (1446: P1 0.965, P2 0.979). | **Generalisation unverified.** No homology clustering exists, so independent-group counts and family-controlled generalisation remain unmeasured — stated at every stage and still true. |

### Explicitly unanswered

**Biological mechanism.** Nothing shows that 1446, 2564 or 3900 *is* a signal
peptide detector in any mechanistic sense, that the model uses it to decide
localisation, or that changing it changes where a protein would go in a cell. The
P4 effect is on the model's own masked-token distribution at the edited residue
in the same forward pass — no propagation, no classifier, no cell.

**Evaluator dependence.** Every external claim rests on WoLF PSORT, a 2006 kNN
predictor over UniProt-derived annotations. §4 shows it resolves only 3.92
distinct margins per protein across 5.25 distinct sequences. 12.5% of the pilot
sample sits inside its reference set. No second evaluator was used, so predictor
dependence is total and unquantified.

**Generalisation.** Unmeasured homology dependence everywhere; one kingdom in the
P3 pilot; one organism model; P4 already used at Check-in 1; discovery,
tuning and evaluation all inside a single cohort with one inclusion history.

---

## 6. Bounded completion recommendation

### What can responsibly be reported now

1. **Discovery.** Single layer-8 SAE features separate extracellular from
   cytoplasmic proteins in this cohort, strongly and stably across kingdom and
   length strata. The best of them, 1446, is **strongly associated with
   signal-peptide annotation, including in the label-discordant groups** — the 59
   cytoplasmic proteins with a signal annotation activate it about as strongly as
   the 1,252 extracellular ones. Report as an association, with the
   annotation-provenance caveat, and without claiming it tracks annotation
   exclusively or that any mechanism has been identified.
2. **Candidate relationships — report the measurements, not a verdict.** Over the
   190 pairs: per-protein Spearman median **0.697** whole-sample but **0.139 /
   0.179** within Extracellular / Cytoplasm; activation-set Jaccard median
   **0.536**; decoder-direction cosine median **0.008**. Report those four
   numbers rather than calling the candidates interchangeable or redundant: the
   whole-sample correlation is substantially driven by the shared class signal,
   and nothing here establishes that any two of them are substitutable.
3. **P4 intervention.** One contrast of four shows differential model
   sensitivity: ablating 2564 changes the model's masked-token distribution more
   than norm-matched controls do, conditional on a matched control existing, and
   that conditioning restricts it to 61 of 400 proteins, mostly extracellular,
   mostly at residue 20. Report with those four bounds attached; it is an
   internal-sensitivity result, not a localisation result.
4. **P3 sequence pilot.** A negative, well-instrumented result: intervening on
   1446 did not detectably move WoLF PSORT's margin, and §4 shows the main reason
   is that the intervention rarely changed the sampled sequence at all. Report as
   a feasibility and sensitivity finding about the pipeline, not as evidence
   against the hypothesis.
5. **Method findings worth reporting in their own right.** The WoLF PSORT
   printed-line omission (omitted ≠ zero), the residual-space norm-matching
   correction, and the structural loss of ablation controls under norm matching.

### Would a separately frozen variant evaluation answer a useful remaining question?

> **WITHDRAWN.** An earlier version of this section recommended a follow-up
> evaluation conditioned on the intervention actually changing the generated
> sequence — a sample "drawn so that divergence is not a rare event." **That
> recommendation is withdrawn.** Intervention-induced sequence divergence is an
> outcome of the intervention, so selecting on it is outcome-based conditioning
> **even when it is applied before any external scoring**. The divergent cells
> are exactly those where the perturbation crossed an inverse-CDF boundary, which
> is not independent of how large the perturbation's effect was; conditioning on
> it would select cells for having responded and then measure whether they
> responded. Declaring it in advance does not fix that. No such sample should be
> drawn.

**The honest position is that no further intervention experiment is recommended
on the current design.** Conditioning on divergence is invalid for the reason
above; not conditioning on it reproduces the pilot, where 83–87% of paired cells
carry no sequence change. Changing the multiplier, the position schedule or the
candidate to force more divergence is a parameter search, which is also not
recommended.

**What would be legitimate is a change of instrument, not of sample.** Two
options, neither executed or scheduled here:

* **A second, independent evaluator.** Predictor dependence is currently total
  and unquantified (§5). Scoring the *existing* 840 saved sequences with a
  different localisation predictor requires no new generation, no new sample and
  no conditioning of any kind, and it would separate "the intervention had no
  effect" from "WoLF PSORT cannot see it."
* **A readout inside the model.** The saved per-step AA-20 distributions already
  record what the intervention did to the model's state at every edited position,
  on every protein and replicate, with no sampling-resolution loss. That is the
  quantity the pilot actually perturbed, and it is already on disk.

Both work from data that exists. Neither requires selecting a sample on an
outcome.

**Not recommended:** a larger parameter search over multipliers, positions,
candidates or layers. The pilot's limiting factor is not parameter choice, and no
such search is proposed here.

No further experiment is executed in this document.

---

## Sources

| stage | report | data |
|---|---|---|
| Discovery | [candidate_discovery/DISCOVERY_RESULTS.md](candidate_discovery/DISCOVERY_RESULTS.md), [addendum](candidate_discovery/DISCOVERY_RESULTS_ADDENDUM.md) | [v1/](candidate_discovery/v1/), [v2_redundancy/](candidate_discovery/v2_redundancy/) |
| P3 tuning pilot | [intervention_pilot/PILOT_REPORT.md](intervention_pilot/PILOT_REPORT.md), [norm correction](intervention_pilot/NORM_MATCHING_CORRECTION.md) | [v1/](intervention_pilot/v1/), [v2_normmatched/](intervention_pilot/v2_normmatched/) |
| P4 evaluation | [evaluation_p4/EVALUATION_REPORT.md](evaluation_p4/EVALUATION_REPORT.md), [addendum](evaluation_p4/EVALUATION_ADDENDUM.md) | [v1/](evaluation_p4/v1/), [v1_addendum/](evaluation_p4/v1_addendum/) |
| WoLF PSORT readiness | [wolfpsort_setup/WOLFPSORT_READINESS.md](wolfpsort_setup/WOLFPSORT_READINESS.md) | [verification/](wolfpsort_setup/verification/), [reference_audit/](wolfpsort_setup/reference_audit/) |
| P3 sequence pilot | [sequence_variant_pilot/PILOT_RESULTS.md](sequence_variant_pilot/PILOT_RESULTS.md), [spec](sequence_variant_pilot/PILOT_SPEC.md) | [v1/](sequence_variant_pilot/v1/), [v1_retrospective/](sequence_variant_pilot/v1_retrospective/) |
| Check-in 1 | [DEEPLOC_CHECKIN_1.md](DEEPLOC_CHECKIN_1.md), [CHECKIN_1.md](CHECKIN_1.md) | — |
