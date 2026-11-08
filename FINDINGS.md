# MIT 6.3930 — Protein Representation Project: Findings Log

Shared log for both datasets. Append new entries; do not rewrite existing ones.
All links are relative to this file.

- **Dataset 1 — UESolDS** (solubility), from the PLM_Sol V1.1 Zenodo release → [`Violet969-PLM_Sol-3589a95/`](Violet969-PLM_Sol-3589a95/)
- **Dataset 2 — DeepLoc 2.0** (subcellular localization), from DTU Health Tech → [`DeepLoc2/`](DeepLoc2/)

---

## Dataset 1 — UESolDS (PLM_Sol V1.1) — 2026-11-06

### Published splits and label encoding

| Split | File | Records |
|---|---|---|
| Train | [`Violet969-PLM_Sol-3589a95/embedding_dataset/Train_dataset.fasta`](Violet969-PLM_Sol-3589a95/embedding_dataset/Train_dataset.fasta) | 70,031 |
| Validation | [`Violet969-PLM_Sol-3589a95/embedding_dataset/validation_dataset.fasta`](Violet969-PLM_Sol-3589a95/embedding_dataset/validation_dataset.fasta) | 4,000 |
| Test | [`Violet969-PLM_Sol-3589a95/embedding_dataset/test_dataset.fasta`](Violet969-PLM_Sol-3589a95/embedding_dataset/test_dataset.fasta) | 4,000 — **reserved** |

Headers are `>ID <localization>-<solubility>`; localization is the constant placeholder `A`.
Test headers carry an extra leading MD5 field (verified: 4000/4000 equal the MD5 of their sequence).

**Evidence that `1` = soluble and `0` = insoluble** (code in this release):

- [`Violet969-PLM_Sol-3589a95/utils/general.py:10`](Violet969-PLM_Sol-3589a95/utils/general.py#L10) — `SOLUBILITY = ['0', '1']`
- [`Violet969-PLM_Sol-3589a95/datasets/transforms.py:72-87`](Violet969-PLM_Sol-3589a95/datasets/transforms.py#L72-L87) — class `SolubilityToInt`, docstring *"…solubility into 0 for membrane bound OR UNKNOWN or 1 for soluble"*, with `solubility = SOLUBILITY.index(solubility)`, so string `'1'` → integer 1 → soluble
- [`Violet969-PLM_Sol-3589a95/datasets/embeddings_dataset.py:40-47`](Violet969-PLM_Sol-3589a95/datasets/embeddings_dataset.py#L40-L47) — label extracted from the header via `.split('-')[-1]`

*Caveat:* this docstring is inherited from the upstream protein-localization/DeepLoc codebase, where the codes were `S`/`M`/`U`. The release ships no standalone data dictionary, so the direction is code-inferred, not author-stated.

### Train/test overlap and the eligible pool

The initial audit found exactly **one** exact sequence overlap between train and test: train `PDB_230301_28958` ≡ test `959117ec99c9f7d28a16d95836f58ef4 4O65_1` (166 aa, both label 1, no label conflict). We **excluded the training copy**, giving **70,030 eligible training records**. All baseline experiments evaluate on validation; test is untouched.

*Source: the initial audit step was reported, not written to a file — this log is the record of it. All later steps have saved artifacts, linked below.*

### Annotation pilot (100 training proteins, seed 42)

**69 of 100** proteins received a confident source-species assignment and **49 of 100** carried at least one Pfam annotation. 3 were organism-ambiguous (one sequence in ≥2 species), 28 had no exact sequence match. The 28 are tagged/truncated expression constructs, not unknown proteins — only 1 was recoverable by trimming tags, so the rest differ internally.

→ [`Violet969-PLM_Sol-3589a95/pilot_annotation/uesolds_pilot_annotations.tsv`](Violet969-PLM_Sol-3589a95/pilot_annotation/uesolds_pilot_annotations.tsv)

### Fixed development cohort

**10,000 sampled training proteins + all 4,000 published validation proteins = 14,000.**
Stratified by center × label (SG_target) and ID class × label (all others), largest-remainder proportional allocation, records sorted before sampling, seed 42. Reruns reproduce the manifest byte-for-byte.

- Manifest: [`Violet969-PLM_Sol-3589a95/dev_cohort/dev_cohort_manifest.tsv`](Violet969-PLM_Sol-3589a95/dev_cohort/dev_cohort_manifest.tsv)
- **Manifest SHA-256: `b830019cd27c271b8833f3614064f16a5a7394a8ad692e7c4959fbe3ece7b442`**
- Builder: [`Violet969-PLM_Sol-3589a95/dev_cohort/build_cohort.py`](Violet969-PLM_Sol-3589a95/dev_cohort/build_cohort.py)

**Why the key is `record_uid`, not `original_id`:** the published training FASTA **reuses 169 identifiers across 374 records** (205 extra rows). Each repeat carries a *different* sequence — alternative constructs of the same structural-genomics target (e.g. `SSGCID-EnhiA.18428.a` appears 6× at lengths 336/319/338/167/…). No repeated ID has conflicting labels. So one original identifier can refer to several distinct sequences and cannot serve as a primary key. The manifest is keyed on `record_uid` = `<split>:<fasta_record_index>`, with `original_id` retained as an attribute. Validation and test identifiers are unique (4,000/4,000 each).

| Split | n | insoluble (0) | soluble (1) |
|---|---|---|---|
| train (sampled) | 10,000 | 3,937 | 6,063 |
| validation (all) | 4,000 | 2,000 | 2,000 |
| **cohort** | **14,000** | **5,937** | **8,063** |

- **Zero exact train↔validation sequence overlaps** (SHA-256; 14,000/14,000 hashes unique cohort-wide).
- **Homology overlap remains unmeasured** — no clustering (CD-HIT / MMseqs2) has been run, so near-duplicate leakage between splits is still unknown.
- **218 cohort proteins exceed 1,022 residues** (1.56%; 157 train, 61 validation; 145 insoluble vs 73 soluble), relevant to ESM-2 truncation at the 1,024-token context.

### Baseline AUROC (validation) — measurements, not ESM-2 results

Models: L2 logistic regression, fixed C = 1.0, seed 42. **A** = one-hot source proxies (header-derived ID class + SG center); **B** = 20 amino-acid fractions + log length; **C** = both. Preprocessing fit on training only.

| Training set | Evaluation group | A source | B composition+length | C both |
|---|---|---|---|---|
| Full training (70,030) | Pooled validation (4,000) | 0.7792 | 0.6570 | 0.7996 |
| Full training (70,030) | SG_target only (2,959) | 0.6172 | 0.6260 | 0.6479 |
| **Fixed cohort (10,000)** | **Pooled validation (4,000)** | **0.7792** | **0.6565** | **0.7982** |
| **Fixed cohort (10,000)** | SG_target only (2,959) | 0.6172 | 0.6242 | 0.6460 |
| **Fixed cohort (10,000)** | **Unweighted mean, 5 eligible centers** | **0.5000** | **0.4995** | **0.5266** |

Eligible centers (≥20 of each class in validation): EFI, JCSG, MCSG, SECSG, SGPP. Model A is exactly 0.5000 within any single center by construction — its features are constant there.

→ [`Violet969-PLM_Sol-3589a95/baseline_source_vs_sequence/validation_auroc.tsv`](Violet969-PLM_Sol-3589a95/baseline_source_vs_sequence/validation_auroc.tsv) · [`Violet969-PLM_Sol-3589a95/dev_cohort/dev_cohort_validation_auroc.tsv`](Violet969-PLM_Sol-3589a95/dev_cohort/dev_cohort_validation_auroc.tsv)

**These are baseline measurements. ESM-2 has yet to be evaluated.**

### Key limitation carried forward

Header provenance nearly determines the label: in training, `PDB_serial` is 15,084/15,084 soluble and `RefSeq` 567/567 soluble (together ~22% of training, zero negatives); validation repeats this exactly. A model can score ~0.78 pooled by detecting *source database*, not solubility. Macro-averaging across centers removes almost all of it (0.4995 for composition). **Compare ESM-2 against the macro-averaged ~0.50, not the pooled 0.78.**

---

## Dataset 2 — DeepLoc 2.0 — 2026-11-06

### Acquisition

Downloaded from the official DTU data page, following the links labeled *Training/Validation set* and *Test set HPA*:
<https://services.healthtech.dtu.dk/services/DeepLoc-2.0/4-Data.php>
No existing official copies were present in the project beforehand. Original filenames preserved; files untouched.

| File | Source URL | Bytes | SHA-256 |
|---|---|---|---|
| [`DeepLoc2/raw/Swissprot_Train_Validation_dataset.csv`](DeepLoc2/raw/Swissprot_Train_Validation_dataset.csv) | `https://services.healthtech.dtu.dk/services/DeepLoc-2.0/data/Swissprot_Train_Validation_dataset.csv` | 17,724,293 | `29a07b293fed2994a966b70bdcd6bacc59915b8b01fa200cb2b07d8db18384a2` |
| [`DeepLoc2/raw/hpa_testset.csv`](DeepLoc2/raw/hpa_testset.csv) | `https://services.healthtech.dtu.dk/services/DeepLoc-2.0/data/hpa_testset.csv` | 570,089 | `557553bcd35e3dcb59bbfe90e1fd63341549f0a693c0e725d6c2c8387a463b37` |

### Development set structure

28,303 rows × 16 columns: unnamed index, `ACC`, `Kingdom`, `Partition`, `Membrane`, then ten localization columns (`Cytoplasm`, `Nucleus`, `Extracellular`, `Cell membrane`, `Mitochondrion`, `Plastid`, `Endoplasmic reticulum`, `Lysosome/Vacuole`, `Golgi apparatus`, `Peroxisome`), then `Sequence`.

Encoding is **multi-label one-hot**: every label column takes only `0.0` or `1.0` (verified exhaustively). `Partition` ∈ {0,1,2,3,4}; `Kingdom` ∈ {Fungi, Metazoa, Other, Viridiplantae}. `Membrane` is a membrane-association flag, **not** a localization, and is excluded from localization counts. Published `Partition` assignments are preserved untouched.

### Partition × target-label group (Extracellular vs Cytoplasm)

| Group | P0 | P1 | P2 | P3 | P4 | TOTAL |
|---|---|---|---|---|---|---|
| Extracellular positive, Cytoplasm negative | 643 | 520 | 819 | 584 | 575 | **3,141** |
| Cytoplasm positive, Extracellular negative | 2,069 | 1,970 | 1,892 | 1,972 | 1,807 | **9,710** |
| Both positive | 50 | 26 | 36 | 24 | 24 | 160 |
| Both negative | 3,201 | 2,935 | 2,984 | 3,116 | 3,056 | 15,292 |
| Missing or unrecognized | 0 | 0 | 0 | 0 | 0 | **0** |
| **TOTAL rows** | 5,963 | 5,451 | 5,731 | 5,696 | 5,462 | **28,303** |

Groups sum exactly to the row count in every partition and overall (verified in-script).

**Candidate class sizes for the binary contrast: 3,141 Extracellular vs 9,710 Cytoplasm** (12,851 combined, ratio 0.323 : 1 — roughly 3:1 imbalance).

**Multi-localization inside those two groups** (proteins also positive for ≥1 other localization; **kept** in their groups for this audit):

| Group | P0 | P1 | P2 | P3 | P4 | TOTAL | share |
|---|---|---|---|---|---|---|---|
| Extracellular positive, Cytoplasm negative | 66 | 37 | 30 | 31 | 33 | 197 | 6.3% |
| Cytoplasm positive, Extracellular negative | 1,087 | 1,032 | 1,020 | 1,044 | 945 | 5,128 | **52.8%** |

The asymmetry is large: over half the Cytoplasm class is multi-localized, against 6.3% of the Extracellular class.

### Duplicates and cross-partition overlap

- Accessions: 28,303 rows, **28,303 unique `ACC`** — no repeats.
- Sequences: **28,303 unique sequences** — no exact duplicates anywhere in the file.
- **Zero identical sequences spanning more than one partition**, so zero label conflicts across partitions.

The published partitions are clean at exact-match level. **Homology-level overlap is unmeasured here** (DeepLoc 2.0 states its partitions were built by homology partitioning, but we have not independently verified that).

### HPA test set — header only (no rows read, evaluation reserved)

11 columns: `sid`, `Cell membrane`, `Cytoplasm`, `Endoplasmic reticulum`, `Golgi apparatus`, `Lysosome/Vacuole`, `Mitochondrion`, `Nucleus`, `Peroxisome`, `Lengths`, `fasta`.

| Target label | Available in HPA test set? |
|---|---|
| `Cytoplasm` | **Present** |
| `Extracellular` | **Absent** |

**The HPA test set has no Extracellular column.** It covers 8 localizations; Extracellular and Plastid are not among them. The planned Extracellular-vs-Cytoplasm contrast therefore **cannot be evaluated on HPA as published**. This needs a decision before the first localization baseline (see below).

### Outputs

→ [`DeepLoc2/audit/audit_deeploc.py`](DeepLoc2/audit/audit_deeploc.py) · [`DeepLoc2/audit/deeploc_partition_counts.tsv`](DeepLoc2/audit/deeploc_partition_counts.tsv) · [`DeepLoc2/audit/deeploc_audit_results.json`](DeepLoc2/audit/deeploc_audit_results.json)

### Open decisions before the first localization baseline

1. **HPA cannot test Extracellular vs Cytoplasm.** Either (a) hold out one development partition as the test split and use HPA only as a secondary check on labels it does carry, (b) switch the contrast to one HPA supports (e.g. Nucleus vs Cytoplasm — both present), or (c) keep the contrast and accept cross-validation over partitions 0–4 with no independent test set.
2. **Multi-localized proteins.** Kept in their groups for this audit. 52.8% of the Cytoplasm class vs 6.3% of the Extracellular class are also positive elsewhere — a "pure single-localization" variant would shrink and rebalance the classes substantially.
3. **Class imbalance** of roughly 3:1 (9,710 vs 3,141) — decide between class weighting, subsampling, or reporting AUROC/AUPRC only.
4. **Kingdom as a confounder.** Not yet cross-tabulated against the two groups; it is the DeepLoc analogue of the provenance leak found in Dataset 1 and should be checked before trusting any pooled AUROC.

---

## Dataset 2 — DeepLoc 2.0 — first localization baseline — 2026-11-06

### Frozen cohort and evaluation protocol

Binary contrast, frozen here and not to be re-derived downstream:

- **Positive class (target = 1):** `Extracellular == 1` AND `Cytoplasm == 0`
- **Negative class (target = 0):** `Cytoplasm == 1` AND `Extracellular == 0`
- Rows positive for both (160) or neither (15,292) are excluded from the contrast. All other localization columns are retained as **metadata only** and never enter the predictors.

Frozen cohort = **12,851** rows; 12,851/12,851 sequence hashes unique, no duplicates.

| Role | Partitions | n | Cytoplasm (0) | Extracellular (1) |
|---|---|---|---|---|
| Development | 1, 2, 3, 4 | **10,139** | 7,641 | 2,498 |
| **Held out — reserved** | 0 | **2,712** | 2,069 | 643 |

Both expected totals verified in-script by assertion (10,139 / 2,712).

**Evaluation decision.** Partition 0 is reserved for final evaluation and is never read by any model. All baseline scoring runs on partitions 1–4 as **four-fold cross-validation**: each published partition is held out once and the model trains on the other three. Published partition assignments are preserved; no re-partitioning.

**`hpa_testset.csv` stays reserved and cannot score this contrast**: its header carries only `sid`, `Cell membrane`, `Cytoplasm`, `Endoplasmic reticulum`, `Golgi apparatus`, `Lysosome/Vacuole`, `Mitochondrion`, `Nucleus`, `Peroxisome`, `Lengths`, `fasta` — **no `Extracellular` column**. Partition 0 is therefore the only held-out evaluation available for Extracellular vs Cytoplasm.

Manifest (accession, sequence SHA-256, target, kingdom, published partition, length, additional-localization flag, the specific other localizations, split role): [`DeepLoc2/baseline/deeploc_binary_manifest.tsv`](DeepLoc2/baseline/deeploc_binary_manifest.tsv). The `has_additional_localization` flag is computed over the eight non-target localization columns and **excludes the separate `Membrane` flag**, which is a membrane-association indicator rather than a localization.

### Development composition (partitions 1–4)

| Kingdom | n | Cytoplasm (0) | Extracellular (1) | % positive |
|---|---|---|---|---|
| Metazoa | 6,434 | 4,666 | 1,768 | 27.5% |
| Fungi | 2,335 | 1,951 | 384 | 16.4% |
| Viridiplantae | 1,028 | 897 | 131 | 12.7% |
| Other | 342 | 127 | 215 | **62.9%** |
| **TOTAL** | **10,139** | **7,641** | **2,498** | 24.6% |

Per-fold prevalence varies: P1 0.2088, P2 0.3021, P3 0.2285, P4 0.2414.

### Models

L2 logistic regression, fixed **C = 1**, `lbfgs`, seed 42, **no class re-weighting** so observed class frequencies are preserved. One-hot encoders and scalers are fit **inside each training fold**. Accession, partition and every localization annotation column are excluded from the predictors.

- **A** kingdom only
- **B** whole-sequence amino-acid fractions + log(length)
- **C** B + amino-acid fractions over the **first 30 residues** (prespecified window; proteins shorter than 30 would use what they have — 0 such proteins in this cohort)
- **D** C + kingdom

**Non-standard residues.** The cohort contains U (selenocysteine, 59), X (unknown, 20) and B (Asp/Asn, 1). They count in the **denominator** (full window length) but in none of the 20 numerators, so affected rows sum slightly below 1.0. Nothing is imputed and no protein is dropped.

### Results — all development proteins (n = 10,139, prevalence 0.2464)

| Model | AUROC P1 | P2 | P3 | P4 | **mean** | AP P1 | P2 | P3 | P4 | **mean** |
|---|---|---|---|---|---|---|---|---|---|---|
| A kingdom | 0.5793 | 0.6222 | 0.6136 | 0.5985 | **0.6034** | 0.2754 | 0.3770 | 0.3056 | 0.2978 | **0.3139** |
| B sequence | 0.8196 | 0.8048 | 0.8286 | 0.8582 | **0.8278** | 0.5900 | 0.7277 | 0.6572 | 0.7410 | **0.6790** |
| C + N-term 30 | 0.9611 | 0.9554 | 0.9575 | 0.9679 | **0.9605** | 0.8765 | 0.8860 | 0.8764 | 0.9157 | **0.8887** |
| D + kingdom | 0.9661 | 0.9591 | 0.9647 | 0.9722 | **0.9655** | 0.8928 | 0.8952 | 0.8912 | 0.9276 | **0.9017** |
| **Δ C − B (paired)** | +0.1415 | +0.1506 | +0.1290 | +0.1097 | **+0.1327** | +0.2865 | +0.1583 | +0.2192 | +0.1746 | **+0.2097** |

The N-terminal window improves **4/4 folds** on both AUROC and AP.

### Sensitivity — no additional localization annotation (n = 5,967; Cytoplasm 3,600, Extracellular 2,367; prevalence 0.3967)

Same out-of-fold predictions, rescored on the subset.

| Model | AUROC P1 | P2 | P3 | P4 | **mean** | AP P1 | P2 | P3 | P4 | **mean** |
|---|---|---|---|---|---|---|---|---|---|---|
| A kingdom | 0.5788 | 0.6073 | 0.6080 | 0.5844 | **0.5946** | 0.4060 | 0.5411 | 0.4483 | 0.4388 | **0.4585** |
| B sequence | 0.8331 | 0.8106 | 0.8428 | 0.8709 | **0.8393** | 0.7251 | 0.8370 | 0.7836 | 0.8370 | **0.7957** |
| C + N-term 30 | 0.9657 | 0.9625 | 0.9596 | 0.9703 | **0.9645** | 0.9282 | 0.9472 | 0.9321 | 0.9548 | **0.9406** |
| D + kingdom | 0.9694 | 0.9633 | 0.9653 | 0.9733 | **0.9678** | 0.9364 | 0.9467 | 0.9396 | 0.9594 | **0.9455** |
| **Δ C − B (paired)** | +0.1325 | +0.1519 | +0.1167 | +0.0994 | **+0.1252** | +0.2031 | +0.1102 | +0.1485 | +0.1178 | **+0.1449** |

Conclusions are unchanged on single-localization proteins: C still gains over B in 4/4 folds. The subset is strongly asymmetric — 94.8% of positives (2,367/2,498) are single-localization against 47.1% of negatives (3,600/7,641) — which is why prevalence rises from 0.246 to 0.397.

### Interpretation

- **Kingdom alone is weak but not null** (AUROC 0.603). It is a real confounder — "Other" is 62.9% positive against 12.7% for Viridiplantae — so it must stay in the model comparison rather than be assumed harmless.
- **Composition alone is already strong** (0.828), and whole-sequence composition plus length carries most of a conventional signal.
- **The N-terminal window is the single largest gain: +0.133 AUROC and +0.210 AP, consistently in every fold.** The first 30 residues carry substantial predictive information about Extracellular vs Cytoplasm beyond whole-sequence composition. **This is predictive evidence only.** It is consistent with N-terminal secretion signals, but identifying a specific targeting signal requires a dedicated later test — the `SortingSignalsSwissprot.csv` set on the DTU page is one route.
- **Kingdom adds little once sequence is present** (D − C = +0.005 AUROC), so the kingdom effect is largely absorbed by composition.

**This is a baseline measurement. ESM-2 has not been evaluated on DeepLoc.** The bar to beat is **D = 0.9655 AUROC / 0.9017 AP** (all development), or C = 0.9605 / 0.8887 without kingdom.

### Limitations

1. **No independent test set for this contrast.** HPA lacks `Extracellular`, so partition 0 is the only held-out evaluation. Spending it is a one-time decision.
2. **Homology overlap unmeasured.** Exact duplicates are zero within and across partitions, but we have not independently verified DeepLoc 2.0's homology partitioning. Near-duplicate leakage across folds would inflate all four models.
3. **Class imbalance ≈3:1** preserved by design (no re-weighting). AP is reported alongside AUROC because of it.
4. **Multi-localized proteins kept** in the main analysis; the sensitivity table is the check, and it does not change the conclusions.
5. **Truncation untested.** Maximum sequence length in the cohort is 18,562 residues; models B–D use full-length composition, but any ESM-2 comparison will truncate at 1,022, so the comparison will not be strictly like-for-like.
6. **The N-terminal result is positional-composition evidence, not signal-peptide identification.**

### Outputs

→ [`DeepLoc2/baseline/baseline_localization.py`](DeepLoc2/baseline/baseline_localization.py) · [`DeepLoc2/baseline/deeploc_binary_manifest.tsv`](DeepLoc2/baseline/deeploc_binary_manifest.tsv) · [`DeepLoc2/baseline/deeploc_oof_predictions.tsv`](DeepLoc2/baseline/deeploc_oof_predictions.tsv) · [`DeepLoc2/baseline/deeploc_baseline_results.tsv`](DeepLoc2/baseline/deeploc_baseline_results.tsv) · [`DeepLoc2/baseline/deeploc_baseline_results.json`](DeepLoc2/baseline/deeploc_baseline_results.json)

Out-of-fold predictions: one row per development protein (10,139) with one column per model. Reruns reproduce them byte-for-byte (SHA-256 `f316c5dfac094fef0fc1ba5164f41afe6b9e0ea0924526c3c5dfe0e8a4a290aa`).

---

## Both datasets — first frozen ESM-2 representation experiment — 2026-11-07

### Model, implementation, layer indexing

| Item | Value |
|---|---|
| Checkpoint | `facebook/esm2_t33_650M_UR50D` (ESM-2, 33 layers, 1280-dim) |
| Implementation | HuggingFace `transformers.EsmModel` 4.57.6, torch 2.8.0 |
| Mode | `model.eval()` + `torch.inference_mode()`, all parameters frozen — no fine-tuning |
| **Layer** | **`hidden_states[8]`** |
| Layer indexing | HF `hidden_states[0]` is the embedding output *before* block 0; `hidden_states[i]` is the **output of transformer block i**. So `hidden_states[8]` is the output of block 8 — the same tensor fair-esm returns as `representations[8]` (`repr_layers=[8]`). Prespecified: the layer the supplied InterProt paper studied and the layer its released SAE was trained on. |
| Input | **Sequence alone.** No kingdom, header text, or annotation reaches the model. |
| Compute / storage dtype | float16 on Apple M4 Pro MPS / stored float32 |
| Truncation | first **1,022** residues |
| Pooling | mean over **residue tokens only**; `<cls>`, `<eos>` and padding excluded |

Reserved throughout and never read: **UESolDS test**, **DeepLoc partition 0** (2,712 records), **`hpa_testset.csv`**.

### Extraction pilot — verification before scaling

| Check | Result |
|---|---|
| One residue token per residue | Asserted every batch (`attention_mask.sum() == len+2`); 0 failures |
| Ambiguous residues | B, U, Z, O, X all have dedicated ESM-2 vocabulary ids — **zero `<unk>`**, asserted every batch |
| Padding / special tokens excluded from pooling | Re-ran single sequences **unpadded** and compared to batched vectors: cosine **1.00000000** at lengths 118, 459, 2843 |
| Key alignment | All 24,139 inputs SHA-256-verified against the frozen manifests *before* extraction; 0 missing |
| float16 fidelity | cosine **0.999999** vs float32 on identical records (whole and N-terminal), max rel. diff 6.3e-4 |

Measured throughput ~1,580–1,950 residues/s on MPS. GPU was available, so the full run proceeded rather than stopping at pilot. Pooled vectors are written as resumable shards (63 MB + 88 MB).

| Dataset | Sequences | Truncated > 1,022 | Wall time |
|---|---|---|---|
| UESolDS (10,000 train + 4,000 val) | 14,000 | **218** (157 train, 61 val) | 41 min |
| DeepLoc development (P1–4) | 10,139 | **1,111** | 41 min |

Both truncation counts match the frozen manifests exactly — an independent alignment check. Truncated keys are listed in the extraction metadata.

### UESolDS — validation AUROC (fixed cohort, L2 probe, C=1, training-only scaling)

| Model | pooled | EFI | JCSG | MCSG | SECSG | SGPP | **mean of 5** |
|---|---|---|---|---|---|---|---|
| **ESM-2 L8 whole-input mean** | **0.6759** | 0.4594 | 0.6171 | 0.5620 | 0.4631 | 0.4854 | **0.5174** |
| composition + length (matched, first 1022) | 0.6563 | 0.5057 | 0.5327 | 0.4976 | 0.3653 | 0.5931 | 0.4989 |
| composition + length (full length) | 0.6565 | 0.5044 | 0.5338 | 0.4996 | 0.3677 | 0.5919 | 0.4995 |

Per-center soluble/insoluble: EFI 21/71 · JCSG 621/269 · MCSG 160/367 · SECSG 92/208 · SGPP 94/361.

Truncation barely matters here (matched vs full-length differ by ≤0.001), as expected with only 218/14,000 affected.

**Frozen ESM-2 layer 8 beats composition but stays near chance once provenance is controlled**: +0.019 pooled (0.676 vs 0.657), and the unweighted five-centre mean moves only 0.4995 → **0.5174**. Both remain far below the 0.7792 that the **header-text-only** source-proxy baseline reaches. For reference, that source-only model scores exactly 0.5000 within any single center by construction.

### DeepLoc 2.0 — four-fold CV on P1–4 (fold-local scaling, same protocol as the baseline)

All development proteins (n = 10,139, prevalence 0.2464):

| Model | AUROC P1 | P2 | P3 | P4 | **mean** | AP P1 | P2 | P3 | P4 | **mean** |
|---|---|---|---|---|---|---|---|---|---|---|
| matched baseline C | 0.9612 | 0.9572 | 0.9573 | 0.9682 | **0.9610** | 0.8777 | 0.8898 | 0.8776 | 0.9162 | **0.8903** |
| full-length baseline C | 0.9611 | 0.9554 | 0.9575 | 0.9679 | 0.9605 | 0.8765 | 0.8860 | 0.8764 | 0.9157 | 0.8887 |
| **(a) ESM whole-input mean** | 0.9661 | 0.9878 | 0.9815 | 0.9656 | **0.9752** | 0.9392 | 0.9735 | 0.9527 | 0.9213 | **0.9467** |
| **(b) a + contextual N-term mean** | 0.9656 | 0.9928 | 0.9849 | 0.9746 | **0.9795** | 0.9377 | 0.9838 | 0.9546 | 0.9399 | **0.9540** |
| **(c) matched C + both ESM** | 0.9650 | 0.9929 | 0.9845 | 0.9731 | **0.9789** | 0.9371 | 0.9836 | 0.9541 | 0.9401 | **0.9537** |
| Δ a − matched C | +0.0049 | +0.0306 | +0.0242 | −0.0026 | **+0.0143** | +0.0615 | +0.0837 | +0.0750 | +0.0051 | **+0.0563** |
| Δ b − matched C | +0.0044 | +0.0356 | +0.0276 | +0.0064 | **+0.0185** | +0.0600 | +0.0941 | +0.0770 | +0.0237 | **+0.0637** |
| **Δ c − matched C** | +0.0038 | +0.0357 | +0.0271 | +0.0049 | **+0.0179** | +0.0594 | +0.0938 | +0.0764 | +0.0239 | **+0.0634** |

Model (c) beats matched baseline C in **4/4 folds on both AUROC and AP**.

Sensitivity — no additional localization annotation (n = 5,967; Cytoplasm 3,600, Extracellular 2,367; prevalence 0.3967), same out-of-fold predictions:

| Model | AUROC P1 | P2 | P3 | P4 | **mean** | AP P1 | P2 | P3 | P4 | **mean** |
|---|---|---|---|---|---|---|---|---|---|---|
| matched baseline C | 0.9661 | 0.9643 | 0.9593 | 0.9706 | **0.9651** | 0.9292 | 0.9495 | 0.9328 | 0.9548 | **0.9416** |
| (a) ESM whole | 0.9715 | 0.9877 | 0.9822 | 0.9626 | **0.9760** | 0.9637 | 0.9834 | 0.9715 | 0.9430 | **0.9654** |
| (b) a + contextual N-term | 0.9712 | 0.9931 | 0.9850 | 0.9715 | **0.9802** | 0.9593 | 0.9914 | 0.9726 | 0.9574 | **0.9702** |
| (c) matched C + both ESM | 0.9706 | 0.9928 | 0.9844 | 0.9699 | **0.9794** | 0.9584 | 0.9909 | 0.9721 | 0.9571 | **0.9696** |
| Δ c − matched C | +0.0045 | +0.0286 | +0.0251 | −0.0008 | **+0.0143** | +0.0293 | +0.0414 | +0.0393 | +0.0023 | **+0.0281** |

Conclusions hold on single-localization proteins (c beats matched C in 3/4 folds on AUROC, 4/4 on AP).

### Interpretation

- **The two datasets behave very differently, and that contrast is the finding.** On DeepLoc, frozen layer-8 ESM-2 adds a clear, consistent gain over hand-crafted composition (+0.0185 AUROC, +0.064 AP for model b) on top of an already-strong 0.961 baseline. On UESolDS, the same representation moves the provenance-controlled score from 0.4995 to 0.5174 — barely off chance.
- **Adding the contextual N-terminal mean helps on DeepLoc** (b > a by +0.0043 AUROC, +0.0073 AP), consistent with the earlier +0.133 finding for raw N-terminal composition. Positions 1–30 attend to the whole visible sequence, so this is contextual evidence, **not** signal-peptide identification.
- **Concatenating hand-crafted features adds nothing once ESM is present** (c ≈ b, −0.0006 AUROC): ESM-2 layer 8 already subsumes composition and length.
- **UESolDS remains dominated by provenance.** Header text alone scores 0.7792 pooled; frozen ESM-2 reaches 0.6759 pooled and ~0.52 macro-averaged. On this benchmark a strong pooled number still reflects source database more than solubility.

### Limitations and deviations

1. **Fold P4 is the dissenting fold on DeepLoc** — model (a) is slightly *worse* than matched C there (−0.0026 AUROC), and in the sensitivity subset (c) is marginally negative (−0.0008). Gains are real but not uniform.
2. **Single layer, single pooling.** Layer 8 was prespecified for SAE compatibility, not selected as optimal; mean pooling is the only scheme tried. No layer sweep was run.
3. **Truncation at 1,022** affects 1,111/10,139 DeepLoc records (11%) vs 218/14,000 UESolDS (1.6%), so the DeepLoc ESM models see less of each protein than the full-length baseline. The matched baselines control for this, and matched vs full-length differ by ≤0.002.
4. **float16 compute** on MPS; verified equivalent to float32 at cosine 0.999999.
5. **Homology leakage still unmeasured** in both datasets; it would inflate every model here equally.
6. **No test-set numbers.** UESolDS test, DeepLoc P0 and HPA remain untouched.
7. One transient filesystem timeout occurred when first reading the DeepLoc manifest (Desktop sync stall); the file was verified intact and the probe reran cleanly. No data was affected.

**Status: raw frozen-ESM probes complete for both datasets. SAE feature analysis comes next.**

### Outputs

Shared: [`esm2_common/esm_extract.py`](esm2_common/esm_extract.py) · [`esm2_common/build_pairs.py`](esm2_common/build_pairs.py)

UESolDS: [`Violet969-PLM_Sol-3589a95/esm2_probe/probe_uesolds.py`](Violet969-PLM_Sol-3589a95/esm2_probe/probe_uesolds.py) · [`uesolds_esm_results.tsv`](Violet969-PLM_Sol-3589a95/esm2_probe/uesolds_esm_results.tsv) · [`uesolds_esm_val_predictions.tsv`](Violet969-PLM_Sol-3589a95/esm2_probe/uesolds_esm_val_predictions.tsv) · [`extraction metadata`](Violet969-PLM_Sol-3589a95/esm2_probe/emb/uesolds_extraction_metadata.json)

DeepLoc: [`DeepLoc2/esm2_probe/probe_deeploc.py`](DeepLoc2/esm2_probe/probe_deeploc.py) · [`deeploc_esm_results.tsv`](DeepLoc2/esm2_probe/deeploc_esm_results.tsv) · [`deeploc_esm_oof_predictions.tsv`](DeepLoc2/esm2_probe/deeploc_esm_oof_predictions.tsv) · [`extraction metadata`](DeepLoc2/esm2_probe/emb/deeploc_extraction_metadata.json)

---

## Both datasets — layer-8 InterProt SAE compatibility pilot — 2026-11-07

**Status: pilot only. 64 + 64 proteins. No full SAE feature extraction has been run.**

### 1. The released SAE, and whether our tensor is the right one

A layer-8 SAE **does exist** in the authors' release, alongside layers 4, 12, 16, 20, 24,
28, 32 and 33. The repo README documents only layer 24 and states that *"all other layer
SAEs are trained with the same configurations"*.

| Item | Value |
|---|---|
| Repo | `liambai/InterProt-ESM2-SAEs` (Apache-2.0) |
| **Revision pinned** | `3ded2e51641df84d6b89f25430874e3a8eb42c24` (lastModified 2025-02-12) |
| File | `esm2_plm1280_l8_sae4096.safetensors`, 41,997,696 bytes |
| **SHA-256** | `1172503b90c950651e787a5625d91ff1c702a04ff3d851e8e95b700201e56ebd` (downloaded file verified byte-exact against the Hub blob record) |
| Reference code | [`interprot_ref/raw/sae_model.py`](interprot_ref/raw/sae_model.py), SHA-256 `1c5d4bc913d0f6896a2143e0fcc5df013d8aad390f8a5dcc6f6912117a7672fb`, from `etowahadams/interprot` |
| Class | `SparseAutoencoder(d_model=1280, d_hidden=4096)` — constructor defaults, exactly as the README instantiates it |
| **Sparsity** | **Top-K, k = 128** (constructor default; `k` is *not* stored in the checkpoint, so the default is what the authors' own inference path uses) |
| Encoder | `w_enc` (1280×4096), `b_enc` (4096) |
| Decoder | `w_dec` (4096×1280), `b_pre` (1280); decoder column norms measured 0.9914–1.0076, mean 1.00005 → trained unit-norm decoder, as expected |
| Dead features | `stats_last_nonzero` indicates **338 / 4096** features were dead at the end of training |

**Normalization (the authors', reproduced by calling their own class, not reimplemented):**

```
get_acts(x) = topK_128( (LN(x) - b_pre) @ w_enc + b_enc ),  ReLU on the kept values
LN(x): subtract mean over dim=-1, divide by (std over dim=-1 + 1e-5); torch default
       UNBIASED std; no learned scale or shift.
```

#### Activation hook — the correspondence is established by implementation, not by layer number

The authors' reference inference snippet is:

```python
esm_model = EsmModel.from_pretrained("facebook/esm2_t33_650M_UR50D")
outputs   = esm_model(**inputs, output_hidden_states=True)
esm_layer_acts = outputs.hidden_states[LAYER][0]      # (L+2, 1280)
sae_acts  = sae_model.get_acts(esm_layer_acts)
```

That is the **identical checkpoint, identical HuggingFace class, and identical indexing
convention** already recorded for our extraction (see the 2026-11-07 ESM-2 entry above).
The hook is the residual stream read straight off `hidden_states[i]` — no `_t33` fair-esm
re-indexing, no pre/post-LayerNorm variant, no attention- or MLP-output hook. So the
tensor we extract is the tensor this SAE was trained on. **Matching layer number and
dimension alone would not have settled this; matching the implementation path does.**

Two further points were verified numerically rather than assumed:

- **BoS/EoS exclusion is safe.** The reference feeds the full `(L+2, 1280)` tensor
  including `<cls>` and `<eos>`. Because the SAE's `LN` normalizes over `dim=-1`
  (per token, across the 1280 features) there is no cross-token mixing anywhere in
  `get_acts`. Applying the SAE to residue tokens only is therefore *token-identical* to
  applying it to the full sequence and slicing afterwards. Measured on a batch of 4 in
  each dataset: **max |reference − ours| = 0.000e+00, exactly, and the active feature set
  is identical on 4/4 proteins** in both datasets.
- **No unresolved issue blocks scaling.** Reported as required: none.

### 2. Pilot cohort

64 UESolDS **training** records and 64 DeepLoc **development (P1–4)** records, seed 42,
stratified as 2 labels × 4 length quartiles × 8. Sampler
[`esm2_common/sae_sample.py`](esm2_common/sae_sample.py); both samples reproduce
byte-identically on rerun.

| | n | label 0 / 1 | length min / median / max | truncated at 1022 |
|---|---|---|---|---|
| UESolDS | 64 | 32 / 32 | 81 / 232 / 997 | 0 |
| DeepLoc | 64 | 32 / 32 | 63 / 415 / 2616 | 8 |

UESolDS test, DeepLoc partition 0 and `hpa_testset.csv` were **not read**; the eligible
pools are restricted to the training and P1–4 development records and this is asserted
in the sampler.

### 3. Verification table — UESolDS

| Check | Result |
|---|---|
| Residue means reproduce the saved fp16 embeddings | cosine **0.9999995** mean, 0.9999989 min |
| — max abs difference / ‖saved‖ | 0.00108 (max abs diff 0.082) |
| Agreement with reference SAE implementation | **max abs diff 0.000e+00**; active sets identical 4/4 |
| Finite outputs | 64/64 proteins, 0 non-finite |
| Active features per residue | **128 / 128 / 128** (min/mean/max) — exactly k, so no top-128 pre-activation was clipped by ReLU anywhere |
| Reconstruction error, relative to ‖x‖² | 0.120 |
| **Reconstruction error, relative to a per-protein mean-vector baseline** | **0.601** |
| Latents active anywhere in the pilot | **3965 / 4096 = 96.8 %** |
| Pooled vector density (fraction non-zero) | 54.1 % |
| fp16 vs fp32 ESM: per-residue top-k overlap | 99.67 % mean, 99.58 % min |
| **fp16 vs fp32 ESM: pooled feature cosine** | **0.9999983** mean, 0.9999936 min |
| Alignment by record ID and sequence hash | 64/64, order preserved, no duplicates |
| Throughput (fp32, MPS) | 692 residues/s |

### 4. Verification table — DeepLoc

| Check | Result |
|---|---|
| Residue means reproduce the saved fp16 embeddings | cosine **0.9999994** mean, 0.9999989 min |
| — max abs difference / ‖saved‖ | 0.00098 (max abs diff 0.092) |
| Agreement with reference SAE implementation | **max abs diff 0.000e+00**; active sets identical 4/4 |
| Finite outputs | 64/64 proteins, 0 non-finite |
| Active features per residue | **128 / 128 / 128** (min/mean/max) |
| Reconstruction error, relative to ‖x‖² | 0.111 |
| **Reconstruction error, relative to a per-protein mean-vector baseline** | **0.518** |
| Latents active anywhere in the pilot | **4010 / 4096 = 97.9 %** |
| Pooled vector density (whole / N-terminal 1–30) | 60.0 % / 25.5 % |
| fp16 vs fp32 ESM: per-residue top-k overlap | 99.65 % mean, 99.52 % min |
| **fp16 vs fp32 ESM: pooled feature cosine** | **0.9999987** mean, 0.9999951 min |
| Alignment by record ID and sequence hash | 64/64, order preserved, no duplicates |
| Throughput (fp32, MPS) | 668 residues/s |

### 5. Full-run cost estimate

Measured on visible (post-truncation) residue counts read from the frozen manifests.

| | proteins | visible residues | est. runtime, fp32 | est. runtime, fp16 | pooled storage (fp32 dense) |
|---|---|---|---|---|---|
| UESolDS | 14,000 | 4,092,687 | 98.6 min | ~35 min | 229 MB (1 vector) |
| DeepLoc | 10,139 | 4,819,401 | 120.2 min | ~43 min | 332 MB (2 vectors) |
| **Total** | 24,139 | 8,912,088 | **3.6 h** | **~1.3 h** | **561 MB** |

Pooled SAE vectors are **54–60 % dense**, so sparse storage would not pay; dense fp32
is the right format. Per-residue activations are *not* stored — at 128 × 4 bytes ×
8.9 M residues that would be ≈ 9 GB even sparsely, and nothing downstream needs them.

### 6. What this does and does not change about the earlier interpretations

- **Provenance (UESolDS).** Unchanged, and it now propagates. SAE features are a
  deterministic function of the same layer-8 vectors, so the header-provenance confound
  documented above applies to them in full. A latent that appears to "detect solubility"
  on UESolDS may be detecting the source database. Any SAE feature claim on this dataset
  must be scored with the same per-centre macro-average, not pooled; the pooled number
  remains uninterpretable here.
- **Composition.** The earlier claim that layer 8 "subsumes composition and length"
  was a statement about the **mean-pooled** vector, and should be qualified as such.
  The SAE reconstructs only ~40 % (UESolDS) to ~48 % (DeepLoc) of the per-residue
  variance *around each protein's own mean vector*, while capturing ~89 % of raw norm
  energy. Most of a layer-8 residue vector's magnitude is a shared component that both
  mean-pooling and the SAE capture easily; the residue-to-residue variation that would
  carry positional or motif information is largely **not** captured at k=128. Mean
  pooling discards that variation entirely, so "subsumes composition" says nothing about
  what layer 8 encodes per residue.
- **Homology.** Still unmeasured in both datasets, and SAE features inherit it exactly —
  being a deterministic function of the embeddings, they cannot reduce homology leakage
  and offer no independent check on it. Any SAE result carries the same caveat.
- **N-terminal gain (DeepLoc).** Unchanged. The contextual positions 1–30 still attend to
  the whole visible sequence; the pilot shows the N-terminal pooled SAE vector is much
  sparser (25.5 % vs 60.0 % dense), consistent with a narrower feature set, but that is
  descriptive only and still is not signal-peptide identification.
- **A new caveat of its own.** 338/4096 latents were already dead at the end of training,
  and the layer-8 SAE is released without the layer-24 checkpoint's documentation or
  validation. Feature-level claims should be checked against the better-documented
  layer 24 before being relied on.

### 7. Recommendation

**Ready to scale, with one change.** Every correctness check passed, several exactly:
reference agreement is 0.000e+00, active counts are exactly k everywhere, outputs are
finite, and alignment holds by both ID and hash on 128/128 proteins.

The one decision the pilot settles is dtype. fp16 ESM flips ~0.35 % of per-residue top-k
memberships relative to fp32 — small, but not nothing, and top-k is discrete. At the
**pooled** level, which is what the probes consume, the two agree at cosine 0.999998
(min 0.999994). So: **run the scaled pooled-feature extraction in fp16** (~1.3 h total,
561 MB) and keep this measurement on record. If per-residue feature-level analysis is
added later — individual latent interpretation, max-activating-residue work — **re-extract
those in fp32**, where the 0.35 % jitter would otherwise fall exactly on the quantity
being examined.

### Outputs

- Pilot script — [`esm2_common/sae_pilot.py`](esm2_common/sae_pilot.py)
- Sampler — [`esm2_common/sae_sample.py`](esm2_common/sae_sample.py)
- Combined report — [`esm2_common/sae_pilot_report.json`](esm2_common/sae_pilot_report.json)
- Reference SAE code and checkpoint — [`interprot_ref/raw/`](interprot_ref/raw/)
- UESolDS sample / pooled features — [`sae_pilot_sample.tsv`](Violet969-PLM_Sol-3589a95/sae_pilot/sae_pilot_sample.tsv), [`sae_pilot_pooled.npz`](Violet969-PLM_Sol-3589a95/sae_pilot/sae_pilot_pooled.npz)
- DeepLoc sample / pooled features — [`sae_pilot_sample.tsv`](DeepLoc2/sae_pilot/sae_pilot_sample.tsv), [`sae_pilot_pooled.npz`](DeepLoc2/sae_pilot/sae_pilot_pooled.npz)

**Status: SAE compatibility pilot complete. Full pooled-feature extraction not started.**

---

## DeepLoc 2.0 — ESM-2 gain on untruncated proteins only — 2026-11-08

**Read-only check.** No refitting. Uses the saved out-of-fold predictions in
[`DeepLoc2/rerun/esm2_probe/`](DeepLoc2/rerun/esm2_probe/) joined to the verified manifest
by accession. P0 reserved, HPA excluded.

**Question.** The +ESM gain was measured on all 10,139 development proteins, 1,111 of
which exceed the 1,022-residue ESM-2 input limit and are truncated — and truncation is
class-correlated (13.3 % of cytoplasmic vs 3.7 % of extracellular). Does the gain survive
when the truncated proteins are removed, so that every protein's representation covers its
whole sequence?

**Yes.** On the 9,028 proteins whose *original* length is ≤ 1,022 (2,406 extracellular,
6,622 cytoplasmic), `c_matchedC_plus_ESM` beats `matched_C` in **4/4 folds on both
metrics**: mean AUROC 0.9587 → 0.9776 (**+0.0190**), mean AP 0.8935 → 0.9534
(**+0.0599**). These are within 0.001 of the full-cohort values (+0.0179 AUROC,
+0.0634 AP), so the gain is not an artifact of differential truncation.

Join integrity: 10,139/10,139 proteins with exactly one score per classifier, zero
duplicate accessions, zero label or fold mismatches, zero non-finite scores.

**One caveat on recomputed metrics.** Metrics recomputed from the saved predictions file
are slightly *below* the archived in-memory values for the +ESM model: mean AUROC 0.9788
vs the archived 0.9789, worst single fold P2 0.9924 vs 0.9929. Cause: the predictions file
stores probabilities at 6 decimals, and the +ESM model is heavily saturated — only ~800
distinct values per fold, with 544–629 scores at exactly `0.000000` and 158–457 at exactly
`1.000000`. AUROC is rank-based, so those rounding-induced ties take 0.5 credit and pull
the value down. `matched_C` is unsaturated (0 values at either extreme) and reproduces
exactly. The subset figures above therefore **understate** the +ESM column and the deltas
by up to ~0.0005 per fold; the direction and 4/4 counts are unaffected. Archived
in-memory values remain authoritative; recomputation from the stored file is bounded
below by tie-breaking.

### Outputs

- Per-fold table — [`DeepLoc2/rerun/deeploc_esm_lenmatched_le1022.tsv`](DeepLoc2/rerun/deeploc_esm_lenmatched_le1022.tsv)

**Status: check complete. Primary cohort and archived results unchanged.**

---

## DeepLoc 2.0 — single-localisation sensitivity + meeting summary — 2026-11-08

**Read-only.** Saved rerun predictions only, no retraining. P0 reserved, HPA excluded.

On the 5,967 development proteins with no additional localisation annotation (2,367
extracellular / 3,600 cytoplasmic; `Membrane` flag excluded from the definition, no length
filter), `c_matchedC_plus_ESM` beats `matched_C` by **+0.0143 mean AUROC** (0.9651 →
0.9794) and **+0.0281 mean AP** (0.9416 → 0.9696), improving **3/4 folds on AUROC** and
4/4 on AP — P4 is the exception at −0.0008. Against the full cohort's +0.0179 / +0.0634,
roughly half the apparent ESM-2 advantage traces to multi-localised proteins.

**Preprocessing decision: retain all 10,139 development proteins for the primary
comparison.** The single-localisation and ≤1,022-residue subsets stand as sensitivity
analyses. The published partitions are defined over the full cohort, multi-localised
proteins are real biology rather than annotation noise, and the conclusion holds in both
subsets.

Metrics recomputed from the 6-decimal predictions file match the full-precision table to
within −0.00068/+0.00012 AUROC and −0.0036/+0.0010 AP, differing in both directions and
most on AP, because the saturated +ESM scores round into ties. Full-precision values are
authoritative; see the summary's precision note.

### Outputs

- Meeting summary — [`DeepLoc2/rerun/meeting_summary.md`](DeepLoc2/rerun/meeting_summary.md)

**Status: DeepLoc single-localisation check complete. Full-cohort result unchanged.**
