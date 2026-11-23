# Localisation interpretability — final results package

MIT 6.393 · Phase 2 · prepared 2026-11-22 from saved outputs only. No experiment
was run, no sequence generated, no sample selected, P0 never read, nothing
committed or pushed. Corrected interpretations are recorded in
[../LOCALISATION_EVIDENCE_SUMMARY.md](../LOCALISATION_EVIDENCE_SUMMARY.md); all
original result files are unchanged.

**One distinction runs through everything below.** We measured effects on a
*computational model* — ESM-2's internal features and its own token
distributions — and, separately, effects on a *sequence-based localisation
predictor*. Neither is a measurement of where a protein goes in a cell. Wherever
this document says "model", it means the model.

---

## Figures

| | |
|---|---|
| **Figure 1** | [Discovery: feature 1446 across localisation-label × signal-annotation groups](figures/fig1_discovery_1446_annotation_light.png) |
| **Figure 2** | [Internal interventions: the four P4 paired KL differences](figures/fig2_p4_internal_interventions_light.png) |
| **Figure 3** | [Sequence pilot: four paired score contrasts and what the ties are made of](figures/fig3_sequence_pilot_light.png) |

Light and dark variants, PNG (200 dpi) and SVG, in [figures/](figures/). Rendered
by [make_figures.py](make_figures.py), which reads the committed result tables and
recomputes nothing.

---

## Narrative

### Association (Figure 1)

Screening 2,867 eligible layer-8 SAE features over 5,190 P1+P2 proteins (1,336
Extracellular, 3,854 Cytoplasm) found features that separate the two classes
sharply. Feature **1446** reaches AUROC **0.974**, and holds across kingdoms
(0.948–0.985) and length quartiles (0.962–0.982).

Conditioning on signal-peptide annotation shows what the separation is riding on.
**1446 is strongly associated with the annotation, including in the
label-discordant groups**: the 59 *cytoplasmic* proteins that carry a signal
annotation activate it about as strongly (98.3% active, mean 4.45) as the 1,252
extracellular ones (99.8%, mean 4.48), while extracellular proteins without an
annotation activate it far less (50.9%, mean 0.95).

The association is not exclusive — 21.5% of the 3,645 unannotated cytoplasmic
proteins still show some activity, and 49.1% of unannotated extracellular
proteins show none — and `not_annotated` is not confirmed biological absence. The
annotation is also largely computational: 1,873 of 1,875 signal evidence rows
carry no named source, so ESM-2 and the annotation may share sequence-based
biases. **No mechanism is identified.**

On the 20 candidates as a set, we report the measurements rather than a verdict.
Over 190 pairs: per-protein Spearman median **0.697** across the whole sample
(range −0.725 to 0.833; 91 pairs above 0.7, none above 0.9) but only **0.139**
within Extracellular and **0.179** within Cytoplasm; activation-set Jaccard
median **0.536**; decoder-direction cosine median **0.008**. The high
whole-sample correlation is substantially the shared class signal. Nothing here
shows any two candidates are substitutable.

### Internal intervention (Figure 2)

Editing a candidate feature at a masked N-terminal residue and comparing against
activation- and norm-matched control features, over 400 P4 proteins (200/200 by
design), **one of four contrasts excludes zero and survives Bonferroni**:
ablating **2564** raises the paired KL difference by **+2.10×10⁻⁵ nats**, 95% BCa
[+1.13, +3.50]×10⁻⁵, 98.75% [+0.91, +4.01]×10⁻⁵.

Four bounds travel with that sentence. The mean is carried by a right tail
(median +0.23×10⁻⁵; 37 of 61 proteins positive). It rests on **61 of 400
proteins** — 18.4% of candidate-active cells, the narrowest coverage in the
design. **55 of those 61 are Extracellular**, and **50 of the 61 paired cells sit
at residue 20**. The coverage loss is structural, not random: matching an
ablation's realised norm requires a large negative control latent change, which
removes exactly the cells where the candidate's latent was largest.

The other three contrasts include zero. For both of 3900's, the *median* protein's
controls move the distribution more than the candidate does.

**The outcome here is the model's own masked-token distribution at the edited
residue in the same forward pass.** No propagation, no classifier, no cell. This
is evidence of differential model sensitivity, not of localisation.

### Sequence generation (Figure 3)

Generating variants by masked resampling under intervention — 40 Metazoa P3
proteins, feature 1446, 10 fixed positions, 4 replicates, five arms, shared
uniforms, 8,000 forward passes, 840 sequences scored locally — gives **a median
paired difference of exactly zero in all four contrasts** (means +0.037 to
+0.106; exploratory 95% intervals all containing zero).

The decomposition of those zeros is the substantive result. Of 160 paired cells
per contrast, **132–139 are cells where the two arms generated the identical
sequence.** Only 5–9 are cases of different sequences receiving the same
predictor score. **The pilot mostly did not change the sequence, so the zero
medians are not an evaluator verdict.**

Two readings that would be wrong. First, an identical sequence does **not** mean
identical token distributions: the edit was applied and verified at every active
step, and the arms' AA-20 distributions genuinely differed there — the shared
uniform simply fell in the same inverse-CDF bin. "Unchanged generation" describes
the realised sequence, not the model's state. Second, evaluator insensitivity is
real but secondary: among cells where the sequences *do* differ, 24–33% still
score identically, and across the run 800 variants collapse to 210 distinct
sequences but only 3.92 distinct margins per protein.

Instrumentation was clean throughout: arm A an exact no-op, zero unintended
top-K changes, zero latent coordinates touched in the random-direction arms,
control norm matching to 3.2×10⁻⁷ relative, 840/840 sequences scored with none
lost.

---

## Completed experimental work vs unresolved scientific questions

### Completed — and what each one settles

| Work | Settled |
|---|---|
| Feature screen over 2,867 features × 5,190 proteins | Which layer-8 SAE features separate the classes, and how stably across kingdom and length strata |
| Annotation-conditioned analysis | That 1446's activation is strongly associated with signal-peptide annotation, including in label-discordant groups, and that the association is not exclusive |
| Pairwise redundancy measurement (190 pairs) | The measured overall and within-class behavioural correlations and the decoder-direction geometry |
| P4 intervention, 400 proteins, 13,500 passes | That ablating 2564 changes the model's own masked-token distribution more than matched controls do, on the subset where a matched control is constructible |
| P3 sequence-variant pilot, 8,000 passes, 840 sequences | That the pipeline runs end to end and is clean; and that at this budget the intervention rarely changes the generated sequence |
| WoLF PSORT local install + 37-check verification | That a reproducible local evaluator exists, is deterministic, and that a class omitted from its printed line does **not** mean zero |
| Reference-overlap audit | That 5 of the 40 pilot proteins are inside the predictor's own reference set |

### Unresolved — and why

**Direct versus indirect localisation encoding remains unresolved.** This is the
central open question and nothing in this package decides it. Two accounts fit
every measurement we have: (a) ESM-2 represents subcellular localisation and
signal-peptide features are part of how it does so; (b) ESM-2 represents
signal-peptide-like sequence patterns, and localisation separation follows
*indirectly* because the two correlate in this cohort. Figure 1 is equally
consistent with both — a feature that fires on annotated cytoplasmic proteins is
what (b) predicts, and is also what (a) predicts if the model represents the
targeting signal as one input among several. **We cannot presently distinguish a
feature that encodes localisation from a feature that encodes a correlate of it.**

Also unresolved:

- **Biological mechanism.** No evidence that any candidate *is* a signal-peptide
  detector in a mechanistic sense, that the model uses it to decide localisation,
  or that changing it would change where a protein goes.
- **Evaluator dependence.** Every external claim rests on one 2006 kNN predictor,
  which resolves 3.92 distinct margins per protein across 5.25 distinct
  sequences, and 12.5% of the pilot sample sits inside its reference set. No
  second evaluator was used, so predictor dependence is total and unquantified.
- **Generalisation.** No homology clustering exists for this cohort, so
  independent-group counts and family-controlled generalisation are unmeasured.
  Uncertainty estimates are conditional on an unverified independence assumption;
  no direction or magnitude of interval error is asserted.
- **Attribution to the Check-in 1 gain.** Nothing links these candidates to the
  +0.018 AUROC that motivated the work.

---

## Mapping completed work to the Check-in 1 concerns

| Check-in 1 concern | Completed work | Status |
|---|---|---|
| "We have shown *that* ESM-2 helps; we have not shown *which biological features* explain the gain." | Feature screen; 1446 at AUROC 0.974, stable across strata | **Open.** A feature that separates classes is not a feature that explains the probe's gain. Nothing connects these candidates to the +0.018 AUROC. |
| "Comparison 2 adds the whole-input embedding *and* the first-30 embedding at the same time… leaves open which positions account for it." | Screening on a first-30 summary; interventions at single N-terminal residues | **Partly addressed.** The window is N-terminal, but ESM-2 is bidirectional, so a first-30 activation still carries whole-sequence context. Position attribution not established. |
| "Our N-terminal embedding is pooled from positions that attend to the whole sequence." | ≤1022 stratification (ΔAUROC ≤ 0.0051); length-quartile AUROCs 0.962–0.982 | **Open.** Robustness to truncation and length is not locality of information. |
| "The next test: does an internal ESM-2 feature overlap with annotated signal peptides?" | Annotation-conditioned analysis across four label × annotation groups (Figure 1) | **Answered.** 1446 is strongly associated with the annotation, including in the 59 label-discordant cytoplasmic proteins. Not exclusive; annotation is largely computational; no mechanism identified. |
| "Overlap alone does not establish what ESM-2 encodes… testing causation requires intervention." | P4 internal intervention (Figure 2); P3 sequence-variant pilot (Figure 3) | **Open.** P4 shows an effect on the model's own token distribution, not on localisation. P3 found no external effect at this budget, mostly because 83–87% of paired cells carried no sequence change. |
| "The gain varied considerably across folds." | Discovery stability reported without a pass quota (1446: P1 0.965, P2 0.979) | **Open.** No homology clustering, so independent-group counts and family-controlled generalisation remain unmeasured. |

---

## Presentation script (≈6 minutes)

> **Slide 0 — the frame.** "Everything I'm about to show is a measurement on a
> model. One part is ESM-2's own internal state; another is a sequence-based
> localisation predictor. Neither is a measurement of where a protein actually
> goes in a cell. I'll keep saying 'the model' because I mean the model."

> **Slide 1 — where we left off.** "At Check-in 1 we had shown *that* ESM-2's
> layer-8 features help predict extracellular versus cytoplasmic localisation —
> mean validation AUROC went from 0.961 to 0.979. We said the open question was
> *which* features. We also said the next test would be whether an internal
> feature lines up with annotated signal peptides, and that overlap alone would
> not settle causation. Phase 2 ran that test, and then two interventions."

> **Slide 2 — Figure 1. Association.** "We screened 2,867 sparse-autoencoder
> features. Feature 1446 separates the classes at AUROC 0.974. Now split the
> proteins four ways, by localisation label and by whether a signal peptide is
> annotated. Look at the top two rows: extracellular proteins with an annotation,
> and *cytoplasmic* proteins with an annotation. They are nearly identical — 99.8%
> versus 98.3% active, means 4.48 versus 4.45. The feature is strongly associated
> with the annotation, and that holds even where the annotation and the
> localisation label disagree.
>
> Three things it is not. It is not exclusive — about a fifth of unannotated
> cytoplasmic proteins still show some activity. 'Not annotated' is not confirmed
> absence; it means nothing was recorded. And this identifies no mechanism: it's
> an association measured on screening data."

> **Slide 3 — Figure 2. Internal intervention.** "Association isn't causation, so
> we intervened. We edited the feature at a masked N-terminal residue and asked
> whether the model's own prediction for that residue moves more than it does for
> matched control features. Four contrasts, 400 proteins. One excludes zero and
> survives correction: ablating feature 2564, about two parts in a hundred
> thousand nats.
>
> That is a real effect and a small one, and I want to be precise about what it
> covers. It rests on 61 of the 400 proteins — the narrowest coverage in the
> design — because you can only build a matched control for some cells, and the
> ones you lose are systematically the ones with the biggest edits. Fifty-five of
> those 61 proteins are extracellular, and 50 of the 61 cells are at a single
> residue position. The mean is carried by a tail.
>
> And the outcome is the model's own token distribution at the residue we edited,
> in the same forward pass. No propagation, no classifier, no cell. This is
> differential model sensitivity. It is not localisation."

> **Slide 4 — Figure 3. Sequence generation.** "So we tried to go all the way to
> sequences. Generate variants by masked resampling while intervening on 1446,
> then score them with a local predictor. Forty proteins, five arms, 8,000 forward
> passes, 840 sequences scored.
>
> Every contrast's median is exactly zero. The temptation is to read that as 'the
> intervention does nothing.' The right-hand panel is why that reading is wrong.
> Of 160 paired cells, about 135 are cells where the two arms produced the
> *identical sequence*. The intervention mostly didn't change what got sampled.
>
> Two careful points. Identical sequence does not mean identical distribution —
> the distributions did differ at every edited position; the shared random number
> just landed in the same bin. And the predictor is genuinely coarse on top of
> that: where sequences did differ, a quarter to a third still scored the same.
> So: mostly unchanged generation, with some evaluator insensitivity underneath."

> **Slide 5 — what's settled and what isn't.** "Settled: which features separate
> the classes; that 1446 is strongly associated with signal annotation including
> in discordant cases; that one internal intervention moves the model's own
> predictions more than matched controls, narrowly; and that the whole generation
> and scoring pipeline works and is clean.
>
> The central question is not settled, and I want to name it exactly. **We cannot
> distinguish whether ESM-2 encodes localisation directly, or encodes
> signal-peptide-like sequence patterns that correlate with localisation in this
> cohort.** Figure 1 is equally consistent with both. Everything downstream
> inherits that ambiguity.
>
> We also have total dependence on a single 2006 predictor, and no homology
> clustering anywhere, so our uncertainty estimates assume an independence we
> haven't verified."

> **Slide 6 — what we're not doing.** "We considered a follow-up that would sample
> only the cases where the intervention actually changed the sequence. We withdrew
> it: divergence is an outcome of the intervention, so selecting on it means
> picking the cells that responded and then measuring whether they responded —
> and declaring it in advance doesn't fix that. We're also not running a parameter
> search.
>
> What's legitimate is changing the instrument, not the sample. We have 840 saved
> sequences a second predictor could score, and we have the per-step token
> distributions already on disk, which is the quantity we actually perturbed. Both
> work from data that already exists."

---

## Provenance

| | |
|---|---|
| Discovery | config `33a3b89d…`, candidate list `070d0f43…`, 5,190 proteins |
| P4 evaluation | config `cc677954…`, sample `9ecaa7e7…`, 400 proteins, 13,500 passes, 2,592.6 s |
| P3 sequence pilot | protocol `2cffcda6…`, sample `c121ffa5…`, 40 proteins, 8,000 passes, 1,645 s; scoring 840 sequences, 266.7 s |
| WoLF PSORT | `wolfpsort:0.2`, commit `23786c0d`, linux/386 container, `--network none`, 37/37 checks |
| Figures | [make_figures.py](make_figures.py) — reads committed tables only |
