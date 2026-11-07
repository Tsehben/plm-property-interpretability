# Check-in 1 — Research Brief

**6.3930 Final Project · Ebenezer Tseh and Hoi · 2026-11-07**
**Format:** 7-minute presentation covering the whole project, then 4 minutes of questions.

**Submitted proposal:** *Mechanistic Interpretability of Property Representations in Protein
Language Models.* This brief tracks that plan directly — frozen ESM-2-650M, InterProt SAEs,
UESolDS, and the DeepLoc Extracellular-versus-Cytoplasm contrast. Phase 1 (property probes,
biological-concept probes, representation overlap) is partly complete; Phase 2 (feature
ablations, sequence variants) is planned.

> *Submission history per Ebenezer; the authoritative submitted text is not in the shared
> folder. A similarly titled local PDF describing steering of ProGen2/ProLLaMA is an **earlier
> draft**, superseded, and is not the submitted plan.*

**Status:** property probes complete for both properties; SAE feasibility verified; concept
probes and representation overlap not yet run.

**One change to report:** both labels come from DeepLoc's Swiss-Prot training/validation file,
which has them. The *separate* HPA test file has no Extracellular column, so it cannot evaluate
our contrast. We therefore reserve **DeepLoc partition 0** as an internal Swiss-Prot holdout.
External validation for this contrast **remains open** ([§3](#3-data-and-access)).

---

## 1. What we are trying to do, and why it matters

**Question:** are solubility and subcellular localization represented in ESM-2 through
*distinct, property-associated features*, or through *biochemical features, regional patterns,
and biological identity* — composition, local sequence regions, and family or organism?

A model can get a property right for three different reasons:

| Route | What the model is doing | Said plainly |
|---|---|---|
| **Biological identity** (family) | Recognising the evolutionary family and recalling its typical property | "This looks like a globin, and globins are soluble." |
| **Biochemical features** (composition) | Reading bulk amino-acid statistics and length | "This is 12 % lysine and short, so probably soluble." |
| **Regional pattern** (targeting motif) | Detecting a localised sequence element that *causes* the property | "The first 20 residues are a hydrophobic core flanked by positive charge — a signal peptide — so this is secreted." |

**Concrete example.** A secreted human protease. All three routes predict "extracellular,"
but only the third is the mechanism: the secretory machinery physically binds that N-terminal
signal peptide and cleaves it. Family and composition merely *travel with* secretion in the data.

**Why the distinction matters for unfamiliar proteins.** A family-based representation has
nothing to retrieve for a *de novo* designed protein with no relatives. A composition-based one
breaks on anything compositionally atypical for its class. A motif-based one transfers, because
the motif *is* the mechanism — attach a signal peptide to a novel sequence and it really is
secreted. **Two models with identical benchmark AUROC can behave completely differently the
moment you hand them a designed sequence.** Benchmark accuracy cannot tell them apart. That is
the gap this project targets, and it is the gap that matters if these models are used to screen
generated proteins.

## 2. Prior work and our contribution

- **Probe** — a simple classifier (L2 logistic regression) trained on *frozen* activations. It
  measures what is linearly *decodable*; it does **not** show the model uses that information.
- **Sparse autoencoder (SAE)** — a wide layer trained to reconstruct activations with few units
  active at a time, re-expressing each 1,280-dim vector as a handful of features from a large
  dictionary. It aims to undo *superposition* so individual features become interpretable.

| Work | One verified finding | Method | Question it leaves us |
|---|---|---|---|
| **Rives et al., PNAS 2021** [[1]](#r1) | Mean-pooling hidden states into one vector per protein gives a space whose **primary axes of variation are species and orthology**; orthologous genes cluster, and nearest-neighbour queries detect remote homologs at fold level comparably to HHblits. | Mean-pooled embeddings; t-SNE, PCA, nearest-neighbour retrieval, linear projections. | Mean pooling is *exactly* our pooling. If family is the dominant axis, how much of any property probe is family readout? |
| **InterPLM — Simon & Zou, Nat. Methods 2025** [[2]](#r2) | SAE features align with far more biology than neurons: up to **2,309** strongly concept-aligned features per layer, with detected concepts expanding from **15** (neurons) to **143** (SAE features). | SAEs on ESM-2-8M (6 layers) and 6 of ESM-2-650M's 33 layers (L1, L9, L18, L24, L30, L33); scored against Swiss-Prot by F1. | Features align with *annotations*. Does an interpretable feature carry the predictive signal for a property, or merely co-occur with it? |
| **InterProt — Adams et al., bioRxiv 2025** [[3]](#r3) | Linear probes on SAE latents are **competitive with probes on raw ESM embeddings** across layers. | TopK SAEs (k=128, 4,096 features) on the ESM-2-650M residual stream; per-layer weights released. | If SAE probes match embedding probes, the SAE buys interpretability rather than accuracy. What does that interpretability let us conclude? |

**InterProt already analysed localization signals — credit where due.** From probe coefficients
at **layer 28**, they report latents firing on an extracellular signal-peptide cleavage site
(**L28/1541**) and across the peptide (**L28/1470**, **L28/1555**), on bipartite nuclear
localization signals (**L28/2375**), and on mitochondrial transit peptides (**L28/3277**). We
are not first to the idea that SAE latents pick up targeting signals. Their interpretations come
from manual inspection and are **not experimentally validated**.

**What we add: controlled comparisons.** Prior work shows the information is *present*. It does
not separate the property from the concepts correlated with it. Our three controls, all built:

1. **Provenance control.** UESolDS header text alone predicts solubility at pooled AUROC
   **0.7792**, so we score *within* source groups where provenance is constant by construction.
2. **Matched-input control.** ESM features and composition are computed on the **same visible
   residues and same visible length**, so truncation cannot confound the comparison.
3. **Like-for-like increments.** Every gain is a paired per-fold difference against its own
   matched baseline, never against another experiment's number.

## 3. Data and access

Both datasets are downloaded, audited, hash-verified, and frozen. We have full access.

| | **UESolDS** (solubility) | **DeepLoc 2.0** (localization) |
|---|---|---|
| **Labels** | 1 = soluble, 0 = insoluble; verified in release code | Multi-label: 10 localizations as **separate 0/1 columns**. Positive = `Extracellular=1 AND Cytoplasm=0`; negative = `Cytoplasm=1 AND Extracellular=0`; rows with **both 1 or both 0 excluded**. Other localizations are allowed |
| **Release** | PLM_Sol V1.1 (Zenodo) | 28,303 records, SHA-256 `29a07b29…8384a2` |
| **Eligible** | 70,030 training records | **12,851** (3,141 Extracellular-pos + 9,710 Cytoplasm-pos; 160 both-positive and 15,292 both-negative excluded) |
| **Analysis cohort** | **10,000 train** (seed 42, stratified) + **4,000 validation** | **10,139** in partitions 1–4 (2,498 / 7,641; prevalence 0.2464) |
| **Evaluation** | Published validation, 2,000/2,000 balanced; **test reserved** | **4-fold CV** over published P1–P4; **P0 (2,712) reserved**; HPA reserved |
| **Main limitation** | **Provenance leakage** — PDB_serial (15,084/15,084) and RefSeq (567/567) are 100 % soluble in training, repeated in validation | **Imbalance (≈3:1), kingdom confounding, and no external test** for this contrast |

**Counts reconcile:** 28,303 released → 12,851 eligible → 2,712 reserved in P0 → **10,139**
in development. Other localization labels are **kept as metadata**, not discarded; a flag marks
additional positive localizations (the separate Membrane flag is excluded from it). In P1–4,
4,172 of 10,139 carry an additional localization, leaving a **5,967-protein single-localization
sensitivity subset** we re-score every model on.

**The evaluation change — two different files.** DeepLoc 2.0 ships a Swiss-Prot
training/validation set and a separate, independent HPA test set [[4]](#r4), and only one of
them carries our labels.

*The Swiss-Prot training/validation file **has both labels.*** Its header is
`, ACC, Kingdom, Partition, Membrane, Cytoplasm, Nucleus, Extracellular, Cell membrane,
Mitochondrion, Plastid, Endoplasmic reticulum, Lysosome/Vacuole, Golgi apparatus, Peroxisome,
Sequence` — ten localization columns including **Extracellular**. We built our binary cohort
from these annotations, and we had extracellular proteins throughout development. Swiss-Prot
contains **3,301** proteins labeled extracellular. Of these, **160** are also labeled
cytoplasmic. Excluding those leaves **3,141** positives for our binary comparison, with
**2,498** in development (P1–P4) and **643** in the reserved partition (P0).

*The separate HPA test file **lacks one of them.*** Its header is `sid, Cell membrane,
Cytoplasm, Endoplasmic reticulum, Golgi apparatus, Lysosome/Vacuole, Mitochondrion, Nucleus,
Peroxisome, Lengths, fasta` — eight localization columns, and **no Extracellular column**
(the paper's Table 2 scores six of them). Extracellular status is simply *unknown* for those
proteins; it is a missing annotation, not an absence of extracellular protein. A protein there
with `Cytoplasm=0` could be extracellular or could be anything else, so HPA cannot evaluate
this contrast either way.

We therefore reserve **partition 0 as an internal Swiss-Prot holdout** (2,712 proteins). It is
drawn from the same corpus as development, so it tests generalization across homology
partitions — DeepLoc partitions at ≤30 % sequence similarity — but **not** across data sources.
**External validation for this binary contrast remains an open problem.**

## 4. Preliminary evidence

**AUROC** = probability a random positive outranks a random negative (chance = 0.5; insensitive
to imbalance). **Average precision (AP)** = area under the precision–recall curve (chance =
prevalence, 0.2464 here; the more honest metric under imbalance).

### 4.1 UESolDS — the headline is a warning

| Model | **Pooled** validation AUROC | **Five-centre** mean AUROC |
|---|---|---|
| Source proxy only (header text, no sequence) | **0.7792** | 0.5000 (by construction) |
| Composition + length, matched inputs | 0.6563 | **0.4989** |
| Frozen ESM-2 layer 8, whole-input mean | 0.6759 | **0.5174** |

**0.7792 is pooled source-only AUROC** — a model given *only the header string* beats both
sequence models, because provenance is almost perfectly associated with the label. **0.4989
versus 0.5174 is the five-centre comparison** — composition against ESM-2 *within* source
groups, where provenance is held constant. These answer different questions: the first asks
"can you rank the whole validation set," the second "can you rank within one source."

Once provenance is controlled, composition sits at chance and ESM-2 barely moves above it.
**Source is strongly associated with the labels; whether ESM-2 *uses* source-related cues is
still open** — that is a Check-in 2 deliverable, not a claim we make now.

### 4.2 DeepLoc — a real, consistent signal (two distinct experiments)

| | Experiment | Metric | Result |
|---|---|---|---|
| **A** | Adding **first-30-residue composition** to the **full-length** baseline | **mean AUROC** | **0.8278 → 0.9605** (+0.1327, 4/4 folds) |
| **B** | Adding **ESM-2 features** to the **matched-input** baseline | **mean average precision** | **0.8903 → 0.9537** (+0.0634, 4/4 folds) |

> These are **separate experiments with different inputs and different metrics**. A's 0.9605 is
> full-length baseline C's AUROC; B's 0.8903 is *matched* baseline C's AP. They are not two
> points on one curve. Kingdom alone reaches only 0.6034 AUROC.

On the 5,967-protein single-localization subset, experiment B gives 0.9416 → 0.9696 AP — the
effect is not an artifact of multi-localized proteins.

**What this motivates.** A large gain from the first 30 residues *specifically* is what a
targeting signal would produce. But our N-terminal vector is pooled from positions that **attend
to the whole sequence**, so it can carry whole-protein context too. This is **predictive evidence
of an N-terminal contribution, not identification of a targeting signal.** Distinguishing them is
the Phase 1 overlap analysis.

### 4.3 SAE feasibility — verified

A 64 + 64 protein pilot confirms our layer-8 activations reproduce the authors' reference
InterProt implementation **exactly (max difference 0.000e+00)**, with alignment verified by ID
and sequence hash on 128/128 proteins. Full extraction is costed (≈1.3 h, ≈561 MB) and ready.

## 5. Plan, risks, and success

### 5.1 Plan, following the submitted phases

| Phase | Component | Status |
|---|---|---|
| **1** | **Property probes** — predict solubility / localization from frozen representations | ✅ **Complete**, both properties |
| **1** | **Biological-concept probes** — predict family, organism, kingdom, source centre from the *same* representations | ⬜ Check-in 2 |
| **1** | **Representation overlap** — do the features carrying the property coincide with those carrying the concepts? | ⬜ Check-in 2 |
| **2** | **Feature ablations** — remove identified features, measure the cost to each prediction | ⬜ Final |
| **2** | **Sequence variants** — evaluate constructed variants (e.g. signal-peptide transplants) | ⬜ Final |

**Baseline already in hand per property:**

| Property | Baseline | Evaluation group | Metric |
|---|---|---|---|
| Solubility | Matched composition + length — **0.4989** | 5 SG centres (EFI, JCSG, MCSG, SECSG, SGPP) | Unweighted mean AUROC |
| Localization | Matched composition + first-30 — **0.9610** AUROC / **0.8903** AP | 4 folds over published P1–P4 | AUROC and AP, Extracellular positive |

**The next property-versus-concept comparison (Check-in 2).** Run the full SAE extraction, then
fit the *same* fixed-C probes for (a) the property and (b) the biological concepts, on identical
splits, and measure overlap between the feature sets each relies on.

- **Supports distinct property features:** the property signal concentrates in features that the
  concept probes do not use, and for DeepLoc those features fire in residues 1–30 and overlap
  annotated signal peptides.
- **Supports indirect encoding:** the property and concept probes lean on the same features, or
  the property signal is diffuse with no positional preference.

### 5.2 Timeline

Submitted plan, with actual progress marked:

| Weeks | Planned | Actual |
|---|---|---|
| **1–2** | Data curation and representation extraction | ✅ Complete — both datasets audited, frozen, hash-verified; ESM-2 layer-8 embeddings for all 24,139 development records |
| **3–4** | Property probes, concept probes, representation overlap | 🟡 **In progress** — property probes complete and SAE feasibility verified *ahead of schedule*; concept probes and overlap are the **Check-in 2** deliverable |
| **5–6** | Feature ablations and interventions | ⬜ Planned |
| **Final** | Comparison across the two properties | ⬜ Planned |

### 5.3 Risks

| Risk | Status and mitigation |
|---|---|
| **Provenance** (UESolDS) | *Measured, severe, mitigated.* Header text alone gives pooled 0.7792. We never report pooled; we score within SG centres, where source-only is 0.5000 by construction. The honest cost: the leak-resistant metric shows near-chance performance. |
| **Incomplete family annotations** | *Measured, unresolved.* In the 100-protein pilot, **49 of 100** had a Pfam annotation and **51 of 100** did not, using the methods applied. Of those **51**, **23** had exact sequence matches but no Pfam in the queried records; the remaining **28** had no exact match at all — they are tagged, truncated expression constructs (His6, TEV sites, pET vector ends), and a tag-trimming retry recovered only 1 of 29. So family **cannot be fully controlled on UESolDS**, and the concept probe for family will run on partial coverage. |
| **Unverified homology overlap** | *Unmeasured in both datasets.* Exact-sequence overlap is zero everywhere, but homology-level overlap is not measured. DeepLoc's published partitions are homology-partitioned by design (partial protection); UESolDS has no equivalent. Leakage would inflate **every** model alike, so our **paired differences** are more robust than absolute numbers. |
| **Limits on causal interpretation** | *Structural.* Probes and SAE features are correlational. A latent firing on residues 1–30 is predictive evidence; only Phase 2 ablations and sequence variants test causation. We will not claim "ESM-2 detects signal peptides" on probe evidence. |
| **No external test for the DeepLoc contrast** | *Open.* The separate HPA test file lacks the Extracellular column, so extracellular status is unknown there. P0 is an internal Swiss-Prot holdout, not an independent source. Options if time allows: add Nucleus vs Cytoplasm as a secondary contrast, which HPA *can* evaluate. |

### 5.4 How a negative result still answers the question

The question is *to what extent*, not *whether*. A finding that layer-8 ESM-2 carries **no**
solubility information beyond provenance and composition is a direct answer — and we largely
have it (0.5174 vs 0.4989 once provenance is controlled). It is independently useful: it implies
published pooled AUROCs on this benchmark substantially measure **source database**, not
solubility. Equally, if the DeepLoc gain proves diffuse rather than motif-localized, that is
evidence for *indirect* encoding — the other arm of our own hypothesis. The design is symmetric:
both outcomes are reportable.

## 6. The eight course questions

| # | Question | One-sentence answer | Detail |
|---|---|---|---|
| **1** | **What are you trying to do?** | Determine whether ESM-2 represents solubility and subcellular localization through distinct property-associated features, or through biochemical composition, regional sequence patterns, and biological identity. | [§1](#1-what-we-are-trying-to-do-and-why-it-matters) |
| **2** | **How is it done today, and what are the limits?** | Today, properties are predicted by probing or fine-tuning PLM embeddings and reporting pooled benchmark accuracy — which cannot distinguish a model reading a causal motif from one recognising a family, and we show a UESolDS model given *only header text* scores 0.7792 pooled AUROC, beating both sequence models. | [§2](#2-prior-work-and-our-contribution), [§4.1](#41-uesolds--the-headline-is-a-warning) |
| **3** | **What is new in your approach?** | We hold the correlated concepts fixed rather than ignoring them — scoring within provenance groups, matching inputs between ESM and composition baselines, reporting paired per-fold increments, and then directly comparing which SAE features the property and concept probes rely on. | [§2](#2-prior-work-and-our-contribution), [§5.1](#51-plan-following-the-submitted-phases) |
| **4** | **Who cares, and what difference will it make?** | Anyone using a PLM to screen designed or unfamiliar proteins: family-based and motif-based representations give identical benchmark scores but fail differently on sequences with no evolutionary relatives, so knowing which one a model uses determines whether its predictions transfer. | [§1](#1-what-we-are-trying-to-do-and-why-it-matters) |
| **5** | **What are the risks?** | Provenance leakage in UESolDS (measured and mitigated by within-centre scoring), incomplete family annotation (49/100 Pfam coverage), unmeasured homology overlap in both datasets, the correlational limits of probes, and no external test set for the DeepLoc contrast. | [§5.3](#53-risks) |
| **6** | **What is your plan, and do you have the data?** | Phase 1 runs property probes, concept probes, and representation overlap; Phase 2 runs feature ablations and sequence variants — and **we have both datasets already**, downloaded, audited, hash-verified and frozen (UESolDS 10,000 + 4,000; DeepLoc 10,139 with P0 reserved). | [§3](#3-data-and-access), [§5.1](#51-plan-following-the-submitted-phases) |
| **7** | **Rough timeline through final presentations?** | Weeks 1–2 curation and extraction (**complete**), weeks 3–4 probes and overlap (**property probes and the SAE pilot complete; concept probes and overlap are the Check-in 2 deliverable**), weeks 5–6 interventions (planned), final presentation the comparison across properties (planned). | [§5.2](#52-timeline) |
| **8** | **Milestones and criteria for success?** | **Complete:** property probes for both properties and the SAE feasibility pilot. **At Check-in 2:** concept probes and the representation-overlap result — success is a clear verdict either way, since concentrated non-overlapping features support distinct property encoding and shared or diffuse features support indirect encoding. **At the end:** Phase 2 ablations and sequence variants, where success is that removing the identified features costs the property prediction specifically, plus the cross-property comparison. | [§5.1](#51-plan-following-the-submitted-phases), [§5.4](#54-how-a-negative-result-still-answers-the-question) |

---

# Appendix

## A. Supporting evidence and verification

**Headline numbers, with sources.** All read from saved tables, not from memory.

| Number | Meaning | Source |
|---|---|---|
| 0.7792 | Pooled source-only validation AUROC | [`dev_cohort_validation_auroc.tsv`](Violet969-PLM_Sol-3589a95/dev_cohort/dev_cohort_validation_auroc.tsv) |
| 0.4989 / 0.5174 | Five-centre mean AUROC, matched composition vs ESM-2 | [`uesolds_esm_results.tsv`](Violet969-PLM_Sol-3589a95/esm2_probe/uesolds_esm_results.tsv) |
| 0.8278 → 0.9605 | Mean AUROC, +first-30 composition on full-length baseline | [`deeploc_baseline_results.tsv`](DeepLoc2/baseline/deeploc_baseline_results.tsv) |
| 0.8903 → 0.9537 | Mean average precision, +ESM features on matched baseline | [`deeploc_esm_results.tsv`](DeepLoc2/esm2_probe/deeploc_esm_results.tsv) |
| 49 / 51 / 23 / 28 | Pfam and sequence-match coverage in the 100-protein pilot | [`uesolds_pilot_annotations.tsv`](Violet969-PLM_Sol-3589a95/pilot_annotation/uesolds_pilot_annotations.tsv) |
| 0.000e+00 | SAE reference-implementation agreement | [`sae_pilot_report.json`](esm2_common/sae_pilot_report.json) |

**Arithmetic checks.** 3,141 + 9,710 + 160 + 15,292 = 28,303 ✓ · 12,851 − 2,712 = 10,139 ✓ ·
10,139 = 2,498 + 7,641 ✓ · 10,139 − 4,172 = 5,967 ✓ · 49 + 51 = 100 ✓ · 23 + 28 = 51 ✓

**Five SG centres** (soluble / insoluble): EFI 21/71 · JCSG 621/269 · MCSG 160/367 ·
SECSG 92/208 · SGPP 94/361 — 2,264 validation proteins in total.

**Label verification.** Solubility encoding confirmed in release code:
[`utils/general.py:10`](Violet969-PLM_Sol-3589a95/utils/general.py),
[`datasets/transforms.py:72-87`](Violet969-PLM_Sol-3589a95/datasets/transforms.py).

**DeepLoc cohort re-verified against the raw release.** Rebuilding the eligibility rule
directly from [`Swissprot_Train_Validation_dataset.csv`](DeepLoc2/raw/Swissprot_Train_Validation_dataset.csv)
reproduces [`deeploc_binary_manifest.tsv`](DeepLoc2/baseline/deeploc_binary_manifest.tsv)
exactly: same 12,851 accessions, **0** label or partition mismatches, **0** sequence-SHA-256
mismatches, **0** length mismatches, and **0** mismatches on the additional-localization flag.
The 10,139 keys in [`deeploc_pairs.tsv`](DeepLoc2/esm2_probe/deeploc_pairs.tsv) are exactly the
P1–P4 development set, with sequences identical to the CSV. The file's 3,301 Extracellular-
positive rows match the 3,301 reported in the paper's Table 1 [[4]](#r4), confirming we hold the
published dataset. **No re-extraction is required.**

**Integrity.** DeepLoc: zero duplicate accessions, zero cross-partition duplicate sequences,
zero cross-partition label conflicts. UESolDS: `original_id` is **not** a primary key — 169
identifiers are reused across 374 training records as distinct construct variants — so the
manifest key is `record_uid = <split>:<fasta_record_index>`. All splits hash-verified and
byte-identical to the published originals.

**Layer choice.** Layer 8 was **prespecified for our first probe** and verified against a
released InterProt SAE (`esm2_plm1280_l8_sae4096.safetensors`, revision pinned, SHA-256
verified), with our activations reproducing the authors' reference implementation exactly. The
prior localization examples cited in [§2](#2-prior-work-and-our-contribution) are from
**layer 28**, which is where InterProt reports them. A comparison across layers is **optional
future work**, not a correction to the present design.

Full methods, extraction metadata, per-fold tables and all limitations: [`FINDINGS.md`](FINDINGS.md).

## B. Anticipated questions

**Q1. Why these two properties?**
They differ in exactly the way the hypothesis needs. Localization has a known, localised,
*causal* motif — the signal peptide — so "distinct property feature" is sharply decidable.
Solubility has no single causal motif and is a bulk biophysical property, so it tests the other
end. Running both lets the final presentation compare how a PLM represents a motif-driven
property against a diffuse one.

**Q2. Why layer 8, when the localization examples you cite are layer 28?**
Layer 8 was prespecified for the first probe because a layer-8 InterProt SAE is released and we
verified our tensor against it exactly (0.000e+00 difference). The layer-28 latents are cited as
*prior work*, which is where InterProt reports them — not as a benchmark we failed to match. A
later cross-layer comparison is optional future work.

**Q3. Your UESolDS result is essentially chance. Is the dataset usable?**
For pooled benchmarking, no — and that is itself a finding worth reporting. Within-centre it
still supports a valid, harder comparison: ESM 0.5174 vs composition 0.4989. We keep it as the
contrasting case to DeepLoc; a dataset where frozen ESM helps against one where it does not
survive a provenance control is the comparison our final presentation is built on.

**Q4. Isn't the N-terminal result just rediscovering signal peptides?**
That signal peptides sit at the N-terminus is textbook, and InterProt already reported
signal-peptide latents at layer 28. What is *not* established is whether the gain a frozen
mean-pooled representation provides is carried by those motif features or by family and
composition correlates. Our N-terminal vector attends to the whole sequence, so we cannot
currently tell — which is precisely the overlap experiment.

**Q5. How do you know the ESM gain is not homology leakage?**
We do not, fully, and we say so. Exact-sequence overlap is zero and DeepLoc's partitions are
homology-partitioned by design, but residual homology is unmeasured. Leakage would inflate
baseline and ESM models alike, so our **paired per-fold differences** are more robust than
absolute numbers — the +0.0634 AP gain held in 4/4 folds.

**Q6. What is the external validity story for DeepLoc?**
Weaker than we would like, and we flag it rather than paper over it. Our labels come from the
Swiss-Prot file; the separate HPA test file has no Extracellular column, so P0 is an internal
Swiss-Prot holdout — it tests generalization across homology
partitions, not across sources. If time allows, Nucleus vs Cytoplasm is a secondary contrast
that HPA *can* evaluate.

## C. Remaining gaps

1. **No external test set** for Extracellular vs Cytoplasm; P0 is internal to Swiss-Prot.
2. **Homology overlap unmeasured** in both datasets.
3. **Pfam coverage is partial** — 49 of 100 pilot proteins annotated; the family concept probe
   will run on incomplete coverage, and we will report the denominator rather than impute.
4. **The authoritative submitted proposal text is not in the shared folder**; the framing above
   is per Ebenezer. Worth dropping the real document in for the record.
5. **InterPLM count differs by version** — Nature Methods reports 2,309 concept-aligned features
   per layer, the arXiv preprint 2,548. We cite the published figure.

---

## References

<a name="r1"></a>**[1]** Rives A, Meier J, Sercu T, Goyal S, Lin Z, Liu J, Guo D, Ott M, Zitnick CL, Ma J, Fergus R.
*Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences.*
**PNAS** 118(15):e2016239118, 2021. [doi:10.1073/pnas.2016239118](https://doi.org/10.1073/pnas.2016239118) ·
[local PDF](../fall%202026/MIT/final_project/Biological%20structure%20and%20function%20emerge%20from%20scaling%20unsupervised%20learning.pdf)
*Quoted:* "training shapes the representation space so that orthologous genes are clustered";
"linear dimensionality reduction recovers species and orthology as primary axes of variation."

<a name="r2"></a>**[2]** Simon E, Zou J. *InterPLM: discovering interpretable features in protein language models via sparse autoencoders.*
**Nature Methods** 22:2107–2117, 2025. [doi:10.1038/s41592-025-02836-7](https://doi.org/10.1038/s41592-025-02836-7) ·
[local PDF](../fall%202026/MIT/final_project/InterPLM%20-%20discovering%20interpretable%20features%20in%20protein%20language%20models%20via%20sparse%20autoencoders.pdf)
*Quoted:* "up to 2,309 features with strong concept alignment per layer"; "expanding from 15
distinct concepts identified by neurons to 143 identified by SAE features."

<a name="r3"></a>**[3]** Adams E, Bai L, Lee M, Yu Y, AlQuraishi M. *From Mechanistic Interpretability to Mechanistic Biology:
Training, Evaluating, and Interpreting Sparse Autoencoders on Protein Language Models.*
**bioRxiv** 2025.02.06.636901, 2025. [doi:10.1101/2025.02.06.636901](https://doi.org/10.1101/2025.02.06.636901) ·
weights: [liambai/InterProt-ESM2-SAEs](https://huggingface.co/liambai/InterProt-ESM2-SAEs) ·
code: [etowahadams/interprot](https://github.com/etowahadams/interprot)
*Localization latents (layer 28):* L28/1541 signal-peptide cleavage site; L28/1470 and L28/1555
across the peptide; L28/2375 bipartite NLS; L28/3277 mitochondrial transit peptides.

<a name="r4"></a>**[4]** Thumuluri V, Almagro Armenteros JJ, Johansen AR, Nielsen H, Winther O.
*DeepLoc 2.0: multi-label subcellular localization prediction using protein language models.*
**Nucleic Acids Research** 50(W1):W228–W234, 2022.
[doi:10.1093/nar/gkac278](https://doi.org/10.1093/nar/gkac278) ·
[PMC9252801](https://pmc.ncbi.nlm.nih.gov/articles/PMC9252801/)
*Table 1* is the Swiss-Prot cross-validation dataset: ten localization categories, including
**Extracellular (3,301 proteins)**. *Table 2* is the independent Human Protein Atlas test set
and scores six categories — Cytoplasm, Nucleus, Cell membrane, Mitochondrion, Endoplasmic
reticulum, Golgi apparatus — with **no Extracellular**. HPA annotations come from confocal
microscopy and are "used for independent validation."

**Datasets.** PLM_Sol V1.1 / UESolDS (Zenodo) ·
DeepLoc 2.0 downloads, "Training and testing data sets":
[DTU Health Tech](https://services.healthtech.dtu.dk/services/DeepLoc-2.0/) ·
[data page](https://services.healthtech.dtu.dk/services/DeepLoc-2.0/4-Data.php).
The page lists the *Training/Validation set* and the *Test set HPA* as two separate downloads
and states that the Partition column marks five homology partitions (0–4) at "maximum 30%
sequence similarity." 
