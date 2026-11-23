# Biological-Characteristic Labels — Concept Table and Readiness Report

**Status:** complete through coverage and feasibility. **Stopped before any probe training or
property-conditioned performance analysis.** Partitions, roles, targets, the 10,104-protein
eligible core, the completed sweep, annotations, embeddings, probes and predictions are all
unchanged. No downsampling, no refitting, no ESM re-extraction, no interventions, no P0 access,
no commits, no pushes. `git status` reports **0 tracked files modified**.

---

## 1. Shared notebooks on `NESG-analysis`

The branch existed on the remote but not locally. `git fetch origin` brought the
remote-tracking ref only — **the working tree was never switched or modified** (still on
`deeploc`). Files were read with `git show origin/NESG-analysis:<path>`. **No notebook was
executed.**

All four are at `nesg_primary/esm_layer_probing/`:

| Notebook | Runtime | Reads | Writes |
|---|---|---|---|
| `ESM_Extraction_Colab.ipynb` | GPU (T4+) | `extraction_manifest.tsv`, `input_integrity.json`, `protocol.json` | pooled all-layer activations → private HF dataset `PLM-Interp/nesg-esm2-layerwise` |
| `ESM_Train_Validation_AUROC_Colab.ipynb` | CPU | the above + `local_labels.tsv`, `probe.py`, `run_layerwise_auroc.py` | 33 layer probes + composition baseline → HF |
| `ESM_Uncontrolled_Property_Probes_Colab.ipynb` | CPU | + `pfam_memberships.tsv`, `uncontrolled/properties.tsv`, `family_catalog.tsv` | per-layer property probes → HF |
| `ESM_Within_Property_Evaluation_Colab.ipynb` | CPU | + saved validation predictions | metrics recomputed **inside strata**; fits nothing |

**Scope check.** These notebooks target the **NESG solubility** cohort (2,361 train + 691
validation), not DeepLoc localisation. They mount Google Drive for code and use private
Hugging Face as the activation store.

**What can be reused for our localisation analysis, and what cannot:**

| Component | Reusable? |
|---|---|
| `ESM_Extraction_Colab` | **No.** We already hold pooled layers 1–33 for all 10,104 core proteins (`layer_sweep/v1/emb_store/`, 815/815 batches). Re-extracting is explicitly out of scope. |
| `ESM_Train_Validation_AUROC_Colab` | **No.** Our P3 layer curve and frozen layer-26 selection already exist. |
| `ESM_Uncontrolled_Property_Probes_Colab` | **Method only.** Its *design* (one probe per property per layer) maps onto our table, but running it would fit new probes — out of scope here. |
| **`ESM_Within_Property_Evaluation_Colab`** | **Yes, in principle — this is the relevant one.** It fits nothing: it recomputes metrics from *saved validation predictions* inside strata. Our saved P4 predictions ([`p4_predictions.npz`](layer_sweep/v1/p4_predictions.npz)) and probe parameters ([`probe_parameters.npz`](layer_sweep/v1/probe_parameters.npz)) are exactly the inputs that pattern needs. |

**Schema compatibility.** Our concept table deliberately mirrors theirs: a long
`pfam_memberships`-style table (`sha256`, `pfam`, name, `record_id`, `split`) and a wide
`properties`-style table. Their `charge_proxy` is the same `(K+R−D−E)/length` formula and is
likewise documented as *not* pH-dependent charge. One difference: their hydrophobicity is a
**whole-sequence** Kyte–Doolittle mean; ours is over **residues 1–30**, per this task's spec.

> Unrelated but worth knowing: **`origin/deeploc` has moved ahead of your local branch** —
> `8bcb17a "Revise title and scope description in DEEPLOC_CHECKIN_1.md"` sits on top of
> `c2fa725`. Nothing was pulled or merged.

---

## 2. The concept table

[`concept_table/concept_table.tsv`](concept_table/concept_table.tsv) — **10,104 rows, one per
protein**. Long companions: [`concept_pfam_long.tsv`](concept_table/concept_pfam_long.tsv)
(17,037 rows) and [`concept_signal_long.tsv`](concept_table/concept_signal_long.tsv) (2,674
rows). Provenance with input/output hashes:
[`concept_table_provenance.json`](concept_table/concept_table_provenance.json).

Carried through unchanged: accession, frozen sequence hash, original partition, analysis role,
localisation target and class, kingdom, length, isoform flag.

### Derivations, frozen before any property-conditioned analysis

**Hydrophobicity** — mean Kyte–Doolittle over residues 1..min(30, L). KD is defined for the 20
standard amino acids only; a window containing anything else makes the score **missing**, with
the offending residues named. **2 proteins affected** (1 discovery, 1 evaluation).

**Charge balance** — `(count(K)+count(R)−count(D)−count(E)) / window_length` over residues
1..min(30, L). Named **charge balance, not pH-dependent net charge**:

- **histidine is omitted** (pKa ≈ 6.0, only partially protonated at physiological pH, so its
  contribution is not a fixed integer);
- **terminal charges are omitted** — the free α-amino and α-carboxyl groups;
- no pKa model, no pH term, no ionic-strength term.

Missing **only** when the window contains **X, B or Z**, which are ambiguous with respect to
membership of {K,R,D,E} (X = any, B = D/N, Z = E/Q). **1 protein affected** (discovery).

> **Correction applied — O (pyrrolysine).** O was initially grouped with the ambiguous codes.
> That was wrong: pyrrolysine is a *specific* residue, and under this literal `count(K)`
> formula it is not K, R, D or E, so it contributes **exactly zero** to the numerator — the
> same status as U (selenocysteine). Neither now makes the charge score missing. Both still
> make the KD score missing, because KD defines no value for either.
>
> **Windows affected by the correction: 0.** No first-30 window in the 10,104-protein core
> contains O. Across the whole core, first-30 windows contain **U in 1 window, X in 1 window,
> and B, Z, O in none**. The correction changes the *rule*, not any value — reported because
> you asked for the count including if it is zero.

**Pfam** — all recorded accessions retained in long form; **the first ID is never taken as a
single family label**. Entry types come from the InterPro entry catalogue queried **by
accession**; no protein sequence was transmitted.

**Homology** — existing MMseqs2 cluster IDs at 30 / 50 / 80 % identity, reused unchanged from
`homology_audit/v1`. Nothing re-clustered.

**Signal** — existing UniProt status, intervals, evidence categories and the
coordinates-apply-to-frozen-sequence flag. 2,674 features; no prediction made here.

---

## 3. SignalP 6 readiness — **not available, predictions remain missing**

Checked: no `signalp` / `signalp6` / `signalp-6` executable on `PATH`; no `signalp` or
`signalp6` Python module; no conda environment; no matching files under `$HOME`.

**No predictions were generated.** All 10,104 rows carry
`signalp6_status = not_run_no_local_installation` with class, score and cleavage site empty.

**What has *not* been done**, per instruction: no model weights downloaded, no sequence
submitted to any external service, and **no UniProt annotation relabelled as a SignalP
prediction**. The UniProt signal columns remain explicitly named as UniProt annotations and are
a different kind of evidence — mostly predictor-derived themselves (75 % of annotated proteins
carry `sequence_analysis_predictor` as their strongest evidence).

**Precise requirement to unblock:**

1. **Authorised local installation.** SignalP 6.0 is distributed by DTU Health Tech under an
   academic-use licence; the package and its model weights are obtained by registering and
   downloading from their service, not from a public package index. An authorised copy must be
   installed locally, including weights — this is an access requirement, not a technical one.
2. **Configuration to record once available:** organism mode, model mode (the distilled "fast"
   model versus the full model), output format, and the N-terminal window the tool uses.
3. **Resources.** 10,104 sequences on CPU in fast mode is tractable; the full model is
   substantially heavier and benefits from a GPU. Current disk headroom is **3.85 GB**, which
   constrains how large a weight set can be stored — this should be re-measured at install
   time, since swap on this host has repeatedly consumed free space.
4. **Run shape when authorised:** one documented configuration, resumable sharded outputs, and
   per-protein status preserved alongside predicted class, score and cleavage site — the
   columns are already reserved in the table.

---

## 4. Coverage and evaluation feasibility

### 4.1 By analysis role

| | Discovery (P1+P2) | Tuning (P3) | Evaluation (P4) |
|---|---|---|---|
| n | 5,190 | 2,548 | 2,366 |
| Extracellular / Cytoplasm | 1,336 / 3,854 | 584 / 1,964 | 571 / 1,795 |
| Clusters (30 % id) | 3,304 | 1,637 | 1,547 |
| Hydrophobicity missing | 1 | 0 | 1 |
| Charge balance missing | 1 | 0 | 0 |
| No Pfam recorded | 337 | 185 | 159 |
| Signal: present / not annotated / unverifiable | 1,311 / 3,698 / 181 | 582 / 1,871 / 95 | 542 / 1,741 / 83 |
| SignalP 6 predictions | 0 | 0 | 0 |

### 4.2 Pfam-bearing counts — **before vs after the core exclusions**

> **The 9,454 figure quoted earlier was computed over the full 10,139-protein manifest, i.e.
> BEFORE the 35 core exclusions.** For the 10,104-protein core the count is **9,423** — a
> difference of 31 of the 35 excluded proteins carrying at least one Pfam assignment.

Within the core: **2,963** proteins carry ≥ 1 Pfam signature whose own entry type is
**family**; **6,460** carry Pfam signatures but **none of family type**.

### 4.3 Pfam IDs absent from the fitting partitions

[`pfam_ids_absent.tsv`](concept_table/pfam_ids_absent.tsv) lists every absent ID with its own
entry type.

| | All signature types | **Family-type signatures only** |
|---|---|---|
| Distinct IDs, all roles | 4,821 | 1,187 |
| In discovery P1+P2 | 2,933 | 665 |
| **Absent from P1+P2** | **1,888** | **522** |
| **Absent from P1+P2+P3** | **863** | **236** |
| Evaluation IDs absent from P1+P2 | 994 | 265 |
| Evaluation IDs absent from P1+P2+P3 | 863 | 236 |

### 4.4 Evaluation proteins classified by family status

> **Correction applied — entry type.** The classification uses the **Pfam signature's own
> entry type**, not the type of the InterPro entry it is integrated into. A Pfam *domain* can
> sit inside an InterPro entry labelled *family* — verified in our own data: **PF00001 is a
> Pfam `domain` integrated into InterPro IPR000276, which InterPro labels a `family`.** Using
> the integrated label would silently promote 36,530 catalogue-wide domains to families.
> Cohort-wide, only 1,187 of our 4,821 distinct signatures are family-type.

Proteins with no applicable family assignment are kept in **explicit unknown categories**, not
folded into "unseen":

| Status | vs P1+P2 | vs P1+P2+P3 |
|---|---|---|
| All recorded family IDs seen | 222 | 271 |
| **All recorded family IDs unseen** | **499** | **443** |
| Mixed seen / unseen | 12 | 19 |
| *unknown* — Pfam recorded but no family-type signature | 1,395 | 1,395 |
| *unknown* — no Pfam recorded | 154 | 154 |
| *unknown* — applicability uncertain | 83 | 83 |
| *unknown* — entry type unavailable | 1 | 1 |
| **Total** | **2,366** | **2,366** |

**Only-unseen family IDs, by localisation class** (the figure you asked for):

| Reference set | Extracellular | Cytoplasm | Total |
|---|---|---|---|
| vs **P1+P2** | **208** of 571 (36.4 %) | **291** of 1,795 (16.2 %) | 499 |
| vs **P1+P2+P3** | **183** of 571 (32.0 %) | **260** of 1,795 (14.5 %) | 443 |

Extracellular evaluation proteins are **more than twice as likely** to carry only unseen
family identifiers as cytoplasmic ones, under both reference sets.

> **These categories do not establish family independence.** "Unseen" means only that the
> recorded family-type Pfam accession does not appear in the reference partitions under this
> annotation snapshot. It does not mean the protein is evolutionarily unrelated to anything in
> those partitions: 1,633 of 2,366 evaluation proteins (69 %) fall in an *unknown* category
> where no family-type signature is applicable at all; Pfam coverage is incomplete by
> construction; two different family accessions can belong to the same clan; and
> domain-level or sub-threshold homology is invisible to this check. The MMseqs2 clusters in
> the table measure a different, also partial, notion of relatedness.

### 4.5 Prospective strata for the continuous characteristics

Cut points defined on **discovery (P1+P2) observed values only**, quantiles 0.20/0.40/0.60/0.80,
linear method, lower-inclusive / upper-exclusive, outer intervals unbounded, tied edges
collapsed — the same convention as the within-property notebook.
Full table: [`prospective_strata.tsv`](concept_table/prospective_strata.tsv).

**Hydrophobicity (KD, residues 1–30)** — cuts −0.8333, −0.4100, 0.0300, 0.5767 → 5 strata:

| Stratum | Tuning n (EX/CY) | Evaluation n (EX/CY) | Both classes |
|---|---|---|---|
| 1 (< −0.833) | 496 (4 / 492) | 443 (17 / 426) | yes |
| 2 | 544 (8 / 536) | 501 (14 / 487) | yes |
| 3 | 556 (29 / 527) | 512 (**12** / 500) | yes |
| 4 | 432 (111 / 321) | 383 (101 / 282) | yes |
| 5 (≥ 0.577) | 520 (432 / 88) | 526 (427 / 99) | yes |

**Charge balance (residues 1–30)** — cuts −0.0667, −0.0333, 0.0333, 0.0667 → 5 strata:

| Stratum | Tuning n (EX/CY) | Evaluation n (EX/CY) | Both classes |
|---|---|---|---|
| 1 (< −0.0667) | 495 (43 / 452) | 475 (**26** / 449) | yes |
| 2 | 234 (55 / 179) | 248 (32 / 216) | yes |
| 3 | 727 (214 / 513) | 638 (215 / 423) | yes |
| 4 | 361 (122 / 239) | 342 (135 / 207) | yes |
| 5 (≥ 0.0667) | 731 (150 / 581) | 663 (163 / 500) | yes |

**All 10 strata contain both localisation classes in both P3 and P4.** The binding constraint
is the minority class: the smallest evaluation minority cell is **12 extracellular proteins**
(hydrophobicity stratum 3), then **17** (stratum 1) and **14** (stratum 2). Three of five
hydrophobicity strata have fewer than 20 extracellular evaluation proteins, so subgroup
estimates there will be very imprecise, and a cluster-level interval in those cells would rest
on few clusters. Charge balance is better supported throughout (minimum 26).

### 4.6 Family strata

[`family_strata_support.tsv`](concept_table/family_strata_support.tsv).

| | All signature types | **Family-type only** |
|---|---|---|
| Signatures observed in discovery | 2,933 | 665 |
| …with any evaluation member | 688 | 112 |
| …with ≥ 2 evaluation members | 374 | 47 |
| …**with both classes in evaluation** | **26** | **2** |
| …with both classes in tuning | 23 | 0 |

**Within-family evaluation of the localisation contrast is close to infeasible.** Restricted to
family-type signatures, exactly **2** discovery-observed families contain both an extracellular
and a cytoplasmic protein in P4, and **none** does in P3. This is a structural property of the
contrast, not a gap in the table: localisation is strongly family-homogeneous, which is the
same asymmetry the homology audit found (only 0.66 % of 30 % clusters are class-mixed). Any
within-family comparison would be restricted to a couple of families and should not be
presented as a general control.

### 4.7 How the within-property notebook would use this

`ESM_Within_Property_Evaluation_Colab` **fits nothing**. It loads saved validation predictions
and recomputes metrics inside one restriction at a time. Transposed to our setting it would:

1. Load the **frozen** P4 predictions from [`p4_predictions.npz`](layer_sweep/v1/p4_predictions.npz)
   — the same layer-26 and layer-8 probes, unchanged, never refitted.
2. Apply one stratum mask at a time from §4.5 — the cuts being defined on P1+P2 only, so the
   strata are prospective with respect to P3 and P4.
3. Recompute AUROC and AP inside each stratum, with the cluster-aware paired bootstrap already
   specified in `frozen_protocol.json`, and report cluster counts and largest-cluster share per
   stratum alongside n and class counts.
4. Report cells with insufficient minority support as undefined rather than estimating them —
   relevant for the three hydrophobicity strata above.

> Range restriction is an **approximate** control. Holding a property inside a quantile band
> narrows it; it does not remove it, and it does not remove anything correlated with it.
> **No subgroup performance has been calculated. That is the next task, not this one.**

---

## 5. Integrity verification — all pass

| Check | Result |
|---|---|
| One row per protein | **Pass** — 10,104 rows, 10,104 distinct accessions |
| Accession set equals the core manifest | **Pass** |
| Complete accounting of all 10,104 | **Pass** |
| Frozen SHA-256 matches the manifest | **Pass**, all rows |
| Frozen SHA-256 matches the actual sequence | **Pass**, recomputed for all 10,104 |
| Targets unchanged | **Pass** |
| Partitions unchanged | **Pass** |
| Roles unchanged | **Pass** |
| Pfam long join is many-to-one, no duplicate (ACC, pfam) | **Pass** |
| Long-table accessions ⊆ concept table | **Pass**, both |

## 6. Artifacts

| File | Contents |
|---|---|
| [`concept_table.tsv`](concept_table/concept_table.tsv) | 10,104 rows, one per protein |
| [`concept_pfam_long.tsv`](concept_table/concept_pfam_long.tsv) | 17,037 (protein, Pfam) rows with the signature's own entry type |
| [`concept_signal_long.tsv`](concept_table/concept_signal_long.tsv) | 2,674 signal features with intervals, evidence and applicability |
| [`concept_table_provenance.json`](concept_table/concept_table_provenance.json) | derivations, omissions, missing-value rules, input/output hashes |
| [`readiness_assessment.json`](concept_table/readiness_assessment.json) | every number in §4 |
| [`pfam_ids_absent.tsv`](concept_table/pfam_ids_absent.tsv) · [`evaluation_family_status.tsv`](concept_table/evaluation_family_status.tsv) | absence and per-class status |
| [`prospective_strata.tsv`](concept_table/prospective_strata.tsv) · [`family_strata_support.tsv`](concept_table/family_strata_support.tsv) | stratum support |
| [`build_concept_table.py`](concept_table/build_concept_table.py) · [`assess_readiness.py`](concept_table/assess_readiness.py) · [`fetch_pfam_entry_types.py`](concept_table/fetch_pfam_entry_types.py) | code |
