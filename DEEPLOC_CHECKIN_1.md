# DeepLoc Check-in 1 — partner brief

**6.3930 Final Project · Ebenezer Tseh and Hoi · 2026-11-08**
Scope: the DeepLoc subcellular localisation contrast only. Solubility work is covered
separately in [`CHECKIN_1.md`](CHECKIN_1.md).

---

## 1. What are we asking?

Proteins are made in the cell but do not all stay there. Some are **secreted** — exported
outside the cell entirely — while others remain in the **cytoplasm**, the fluid interior.
Where a protein is found is its *subcellular localisation*. As **background biology**, the
cell's export machinery is known to recognise a short stretch at the very beginning of many
secreted proteins, bind it, and route the protein outward.

**What the labels mean.** `Extracellular` and `Cytoplasm` record **where a protein was
observed**. They do not specify how it got there. They are observations, not mechanisms, and
nothing in this brief establishes a route.

**Our question:** how much do frozen ESM-2 layer-8 representations help tell extracellular
proteins from cytoplasmic ones, *beyond* what plain amino-acid composition and sequence
length already provide?

**Two questions, kept separate.** We first ask whether a probe can recover localisation
information from ESM-2. Comparing it with the composition baseline tells us whether ESM-2
adds predictive value beyond those simple features. These are not the same question: a
representation can carry localisation information that the simple features already capture,
so failing to beat the baseline would not mean the information is absent.

Predicting the label is the **first step, not the goal**. Our project asks what ESM-2
internally *represents*. We cannot interpret a representation that carries no signal, so we
first establish that the signal is there and measure its size. Only then is it worth asking
*which* internal features carry it — the concept-probing and feature-overlap work planned
for Phase 1.

---

## 2. Which proteins did we use?

**Source:** `Swissprot_Train_Validation_dataset.csv`, the Swiss-Prot training/validation
file published with DeepLoc 2.0 (28,303 records; SHA-256 `29a07b29…84a2`).

**The two classes**, defined on that file's own columns:

| Class | Definition | n |
|---|---|---:|
| **Positive — extracellular** | `Extracellular` = 1 **and** `Cytoplasm` = 0 | 3,141 |
| **Negative — cytoplasmic** | `Cytoplasm` = 1 **and** `Extracellular` = 0 | 9,710 |

**Why two groups were excluded.** A record with *both* columns set to 1 (160 records) is
annotated in both compartments, so it has no single correct answer for a binary
extracellular-versus-cytoplasmic question. For a record with *both* set to 0 (15,292
records), **neither target location was annotated, so these records were excluded from the
binary comparison.** Including either group would ask the classifier to answer a question the
label does not pose.

**Why other localisation annotations were kept.** Many eligible proteins carry *additional*
labels beyond the target one — a cytoplasmic protein also marked `Nucleus`, for instance.
These are real biology, not annotation errors, and dropping them would discard a large,
non-random part of the cohort. We retained them and treat the restriction to
single-localisation proteins as a *sensitivity analysis* instead (§3, §4).

**Cohort and splits** — the published DeepLoc partitions were preserved exactly, never
re-split:

| | n | Role |
|---|---:|---|
| **Eligible proteins** | **12,851** | Total after the two exclusions above |
| **P1–P4** | **10,139** | Four-fold development; each partition held out once |
| **P0** | **2,712** | **Reserved from model fitting and scoring** — read only for split accounting; held for final evaluation |

Verification: 0 duplicate accessions, 0 identical sequences shared across partitions, 0
missing or unparseable labels.

**All performance in this brief is development performance** — four-fold cross-validation
*within* P1–P4, each partition held out once and scored by classifiers trained on the other
three. **P0 was reserved from model fitting and scoring** — it was read only for split
accounting, and no performance number below comes from it.

**The HPA file is excluded from this contrast.** `hpa_testset.csv` is a separate test file
shipped with DeepLoc 2.0, and **it has no `Extracellular` column** — its columns are `sid`,
`Cell membrane`, `Cytoplasm`, `Endoplasmic reticulum`, `Golgi apparatus`, `Lysosome/Vacuole`,
`Mitochondrion`, `Nucleus`, `Peroxisome`, `Lengths`, `fasta`. Without that column it cannot
express our positive class, so it **cannot score this contrast at all**. It was used for no
training, no validation, and no scoring here. **HPA was excluded because it lacks the
`Extracellular` label.** (Both raw DeepLoc CSVs were downloaded locally and are not tracked
in Git.) This is why P0 serves as our held-out test set, and why *external* validation for
this contrast remains an open gap.

---

## 3. What did the EDA tell us to change?

**Class balance.** Extracellular proteins are the minority at 24.6 % overall, and prevalence
varies by fold — P1 20.9 %, P2 30.2 %, P3 22.8 %, P4 24.1 %. Average precision must
therefore be read against each fold's own base rate, not compared across folds.

**Kingdom.** Prevalence ranges from 12.7 % (Viridiplantae) to 62.9 % (Other), with Metazoa
at 27.5 % and Fungi at 16.4 %. Kingdom alone reaches 0.6034 mean AUROC, so it is a real but
weak predictor on its own. Our reported comparisons use sequence-derived features and never
receive a kingdom label — but that does **not** make them kingdom-free. Sequence features can
carry information associated with kingdom regardless, so withholding the label is not a
control for it.

**Data integrity — clean.** Across all 10,139 development proteins: 0 missing sequences, 0
missing kingdom values, 0 missing embeddings, 0 embedding/sequence-hash mismatches, 0
duplicate accessions. Alignment was checked by both accession and sequence hash. Twenty-two
proteins contain non-standard residues (U, X, B); these are genuine UniProt content, counted
in the denominator but no numerator, with nothing imputed and no protein dropped.

**Truncation is class-correlated.** ESM-2 accepts 1,022 residues. Cytoplasmic proteins are
much longer (median 470 aa vs 221 aa), so **13.3 % of cytoplasmic versus 3.7 % of
extracellular** sequences get truncated — a difference that tracks the label.

**The finding that shaped our checks.** **An additional location was recorded for 52.9 % of
cytoplasmic proteins and 5.2 % of extracellular proteins** (41.7 % of the cytoplasmic class
also carries a `Nucleus` annotation). The two classes differ in how they are **annotated**,
not only in where the proteins are found — note this is a statement about recorded
annotations, not about every place these proteins occur. We retained all eligible proteins
and checked whether the gain persisted among proteins with one recorded location.

**What we did about it.** We did **not** change the cohort. We **retain all 10,139
development proteins for the primary comparison** — the published partitions are defined over
the full cohort, subsetting would break comparability with the source benchmark, and
multi-localised proteins are real biology rather than noise. Instead we **rescored the same
out-of-fold predictions on a subset**: the 5,967 proteins with exactly one recorded location
(2,367 extracellular / 3,600 cytoplasmic; the separate `Membrane` association flag is
excluded from this definition, and **no length filter is applied**). No model was retrained —
the same predictions are simply evaluated on fewer proteins — so the check is clean and costs
nothing. It asks one narrow question: does the gain persist when both classes carry a single
recorded location?

---

## 4. What did we find?

Four-fold cross-validation on P1–P4, each partition held out once. L2 logistic regression,
C = 1, scalers fit inside each training fold only. Means are unweighted across the four folds.
All figures are **development** performance; P0 is not scored.

| # | Comparison | Mean AUROC | Mean AP | Folds improved |
|---|---|---|---|---|
| **1** | Full-length composition + length **→ plus first-30-residue composition** | 0.8278 → **0.9605** | 0.6790 → **0.8887** | **4/4** AUROC, **4/4** AP |
| **2** | Input-length-matched composition baseline **→ plus ESM-2 layer-8** | 0.9610 → **0.9789** | 0.8903 → **0.9537** | **4/4** AUROC, **4/4** AP |
| **3** | Comparison 2 **rescored** on 5,967 single-localisation proteins | 0.9651 → **0.9794** | 0.9416 → **0.9696** | **3/4** AUROC, **4/4** AP |

Per-fold scores for all three: [`DeepLoc2/rerun/meeting_summary.md`](DeepLoc2/rerun/meeting_summary.md).

**Comparison 1 uses no ESM-2 features whatsoever.** Both classifiers count amino acids in
plain text. The second has the whole-protein features **plus**
first-30 counts, so the +0.1327 AUROC is the benefit of **adding those counts to whole-protein
features** — not the performance of the first 30 residues on their own. It is a statement
about *sequence*, not about any language model.

**Comparison 2 is input-matched.** Both classifiers used the same visible sequence, up to
**1,022 amino acids per protein**, and both were scored on the same validation proteins. This
lets us compare their predictions while holding the amount of sequence available constant.

**Comparison 3 evaluates the same two classifiers on a subset** — no refitting. The subset
is the 5,967 proteins with exactly one recorded localisation; the separate `Membrane`
association flag is excluded from that definition, and **no length filter is applied**. The
figures in the table are the rerun's full-precision values from
[`deeploc_esm_results.tsv`](DeepLoc2/rerun/esm2_probe/deeploc_esm_results.tsv).

**The mean gain persisted among proteins with one recorded location, although it was
smaller**: +0.0143 AUROC and +0.0281 AP, against +0.0179 and +0.0634 on the full cohort. The
subset changes two things at once — which proteins are evaluated, and the proportion labelled
extracellular, which rises from 24.6 % to 39.7 %. **The comparison therefore cannot assign a
fraction of the gain to multiple locations.** It is a sensitivity check, not a decomposition.

A further check — also a **rescore of the saved six-decimal out-of-fold predictions**, not a
refit — removed the 1,111 truncated proteins and evaluated the 9,028 that fit within the limit
in full, giving +0.0190 mean AUROC
([`deeploc_esm_lenmatched_le1022.tsv`](DeepLoc2/rerun/deeploc_esm_lenmatched_le1022.tsv)).
What this shows is that **the gain persisted among proteins that fit in full**. It leaves open
how truncation affected the 1,111 longer proteins, which the check removes rather than
examines.

**Where the numbers came from.** All three headline rows above, **including the
single-localisation values**, were computed from scores **before prediction scores were
rounded**, during the rerun itself. The predictions written to
[`deeploc_esm_oof_predictions.tsv`](DeepLoc2/rerun/esm2_probe/deeploc_esm_oof_predictions.tsv)
are **rounded to six decimals**. The later check rescored those saved six-decimal predictions
and **recovered the same qualitative result**, agreeing within 0.0002 mean AUROC and 0.0015
mean AP on the single-localisation subset; the small gap arises because rounding creates ties
that rank-based metrics then split. **The +0.0190 AUROC for the untruncated subset comes from
those saved predictions.** The full-precision table is authoritative for the headline rows.

---

## 5. What can we conclude at Check-in 1?

**Counting amino acids in the first 30 positions improved the composition baseline.** Mean
AUROC rose from 0.8278 to 0.9605 in 4/4 folds. One composition difference is strikingly
positional: **leucine** is
+8.53 percentage points enriched in the first 30 residues of extracellular proteins but only
+0.02 points across the whole sequence — a statement about leucine specifically, not about
every residue. **Leucine enrichment was concentrated in the first 30 positions. Whether ESM-2
recognises a targeting signal that helps explain localisation remains open.**

**Adding the two ESM-2 representations together produced a further gain** over the matched
baseline, improving 4/4 folds on the full cohort. Comparison 2 adds the whole-input embedding
**and** the first-30 embedding at the same time, so the gain cannot be attributed to the
first-30 embedding specifically: **this comparison leaves open which positions or biological
features account for it.** **The gain varied considerably across folds** — +0.036 and +0.027
AUROC in P2 and P3 against +0.004 and +0.005 in P1 and P4 — so the mean is carried by half the
folds. **It persisted among proteins with one recorded location, although it was smaller**
(+0.0143 AUROC, 3/4 folds; P4 is slightly negative at −0.0008).

**What remains open is the mechanism.** We have shown *that* ESM-2 helps; we have not shown
*which biological features* explain the gain. Our N-terminal embedding is pooled from
positions that attend to the whole sequence, so it can carry whole-protein context as well as
local signal. This is **predictive evidence of an N-terminal contribution, not identification
of a targeting signal.** Concept probing, representation-overlap analysis, and feature
ablations are future work.

**The next test we plan.** We will check whether an internal ESM-2 feature overlaps with
**annotated signal peptides** — the short N-terminal export tags the cell's machinery
physically recognises. If the same feature predicts both a signal peptide and extracellular
localisation, that would **support the hypothesis that targeting signals help explain the
localisation prediction**, and our observed N-terminal effect is the pattern that hypothesis
predicts.

What such overlap would *not* do is settle the question. **Overlap alone does not establish
what ESM-2 encodes, and it does not demonstrate a causal mechanism** — two quantities can
covary in the data without either driving the model's output, and **other explanations could
coexist** with a targeting-signal one. Testing causation requires intervention: the feature
ablations and sequence-variant work planned for Phase 2.

**Bottom line for the meeting.** Counting amino acids at the beginning of a protein already
tells us a lot about its localisation. When we gave our classifier that information and then
added ESM-2's layer-8 features, mean validation AUROC rose from 0.9610 to 0.9789. The gain
also persisted among proteins with only one recorded location. We have shown that ESM-2 helps
prediction; the next scientific question is which features explain why.

---

## Supporting files

- DeepLoc results summary, per-fold tables — [`DeepLoc2/rerun/meeting_summary.md`](DeepLoc2/rerun/meeting_summary.md)
- ESM-2 comparison, full precision — [`DeepLoc2/rerun/esm2_probe/deeploc_esm_results.tsv`](DeepLoc2/rerun/esm2_probe/deeploc_esm_results.tsv)
- Out-of-fold predictions, rounded to 6 dp — [`DeepLoc2/rerun/esm2_probe/deeploc_esm_oof_predictions.tsv`](DeepLoc2/rerun/esm2_probe/deeploc_esm_oof_predictions.tsv)
- Composition baselines (comparison 1) — [`DeepLoc2/rerun/baseline/deeploc_baseline_results.tsv`](DeepLoc2/rerun/baseline/deeploc_baseline_results.tsv)
- Untruncated-subset check — [`DeepLoc2/rerun/deeploc_esm_lenmatched_le1022.tsv`](DeepLoc2/rerun/deeploc_esm_lenmatched_le1022.tsv)
- Label and partition audit — [`DeepLoc2/audit/deeploc_partition_counts.tsv`](DeepLoc2/audit/deeploc_partition_counts.tsv) · [`deeploc_audit_results.json`](DeepLoc2/audit/deeploc_audit_results.json)
- Verified cohort manifest — [`DeepLoc2/baseline/deeploc_binary_manifest.tsv`](DeepLoc2/baseline/deeploc_binary_manifest.tsv)
- Full project context, both properties — [`CHECKIN_1.md`](CHECKIN_1.md)
- Detailed running log — [`FINDINGS.md`](FINDINGS.md)
