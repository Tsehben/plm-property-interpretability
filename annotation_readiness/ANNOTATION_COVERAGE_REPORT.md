# Annotation coverage report — DeepLoc localisation development cohort

Generated from `out/annotation_primary.tsv` (10139 rows), `out/annotation_features.tsv` (2786 rows), `out/annotation_domains.tsv` (62382 rows).

Retrieval plan: 10017 base accessions in batches of 100, 268 isoform sequence requests, and 6 individual entry lookups for accessions the batch search either omitted or returned without a sequence. This run: 375 cache hits, 0 network requests, elapsed 1.5s — the fetch is cached and resumable, so a re-run reproduces the tables without re-querying UniProt. Provenance in `out/provenance.json`; raw responses retained under `cache/`.

---

## 1. Retrieval and sequence verification

| Retrieval status | n | share |
|---|---|---|
| ok | 10133 | 99.9% |
| resolved_via_merged_accession | 5 | 0.0% |
| inactive_deleted | 1 | 0.0% |

| Sequence vs frozen extraction input | n | share |
|---|---|---|
| exact | 10019 | 98.8% |
| mismatch | 117 | 1.2% |
| unverifiable | 3 | 0.0% |

**Isoforms.** 268 of 10139 accessions carry a `-N` suffix. 262 matched the frozen sequence exactly; 2 could not have their own isoform sequence retrieved. The canonical sequence was never substituted for an unresolved isoform.

**Merged accessions (5).** These cohort accessions are no longer primary; UniProt resolves them to a successor entry. The successor's record is used, and the substitution is flagged rather than applied silently:

| cohort accession | resolved to | successor id | sequence vs frozen | signal status |
|---|---|---|---|---|
| A0A140G945 | P02489 | CRYAA_HUMAN | exact | not_annotated |
| Q2PE51 | Q90Y12 | BNP_CRODU | mismatch | unverifiable |
| Q9D5U9 | Q8BML1 | MICA2_MOUSE | mismatch | unverifiable |
| G0Y276 | G2WWV6 | PEVD1_VERDV | exact | annotated_present |
| Q0ZHI9 | G2X0L1 | EIX2_VERDV | mismatch | unverifiable |

**Deleted accessions (1).** `Q9Y6N3` (Deleted from Swiss-Prot). No sequence and no annotation are available; these rows are `unverifiable` throughout and are retained in the cohort.

**Coordinate transfer.** 381 of 10139 records are marked `coordinates_transferable = 0`: no UniProt residue interval is applied to them, because the sequence they would be applied to is not the one we verified.

Of those, **261 are isoforms whose own sequence matched the frozen input exactly but differs from the canonical sequence** that UniProt's feature coordinates are defined on. Their parent entries are not blank: **50** carry a signal-peptide annotation, recorded in the `parent_entry_has_signal` column but deliberately **not** transferred. Recovering them would require aligning each isoform to its canonical sequence and mapping the interval across — a bounded, well-defined follow-up that is not performed here, because transferring an interval without that alignment is exactly the error this column exists to prevent.

**Failures, itemised** (122 rows; full detail in `out/annotation_failures.tsv`):

| retrieval_status | sequence_match | n | roles affected |
|---|---|---|---|
| ok | mismatch | 114 | discovery:64, evaluation:18, tuning:32 |
| resolved_via_merged_accession | mismatch | 3 | discovery:3 |
| resolved_via_merged_accession | exact | 2 | discovery:2 |
| ok | unverifiable | 2 | discovery:2 |
| inactive_deleted | unverifiable | 1 | discovery:1 |

---

## 2. Signal-peptide annotation coverage

> `not_annotated` means a matching record was retrieved and carries no signal-peptide annotation. **It is not a confirmed biological negative.** UniProt records evidence for presence; it does not assert absence.

| role | class | n | annotated_present | not_annotated | unverifiable |
|---|---|---|---|---|---|
| discovery | Extracellular | 1339 | 1254 (93.7%) | 53 (4.0%) | 32 (2.4%) |
| discovery | Cytoplasm | 3862 | 59 (1.5%) | 3648 (94.5%) | 155 (4.0%) |
| tuning | Extracellular | 584 | 551 (94.3%) | 22 (3.8%) | 11 (1.9%) |
| tuning | Cytoplasm | 1972 | 31 (1.6%) | 1851 (93.9%) | 90 (4.6%) |
| evaluation | Extracellular | 575 | 515 (89.6%) | 41 (7.1%) | 19 (3.3%) |
| evaluation | Cytoplasm | 1807 | 29 (1.6%) | 1704 (94.3%) | 74 (4.1%) |

**By kingdom (DeepLoc bucket), all roles pooled:**

| kingdom | n | annotated_present | not_annotated | unverifiable |
|---|---|---|---|---|
| Fungi | 2335 | 375 (16.1%) | 1925 (82.4%) | 35 (1.5%) |
| Metazoa | 6434 | 1718 (26.7%) | 4380 (68.1%) | 336 (5.2%) |
| Other | 342 | 214 (62.6%) | 125 (36.5%) | 3 (0.9%) |
| Viridiplantae | 1028 | 132 (12.8%) | 889 (86.5%) | 7 (0.7%) |

Full three-way breakdown (role × class × kingdom) with every counter: `out/coverage_by_role_class_kingdom.tsv`.

---

## 3. Evidence strength, recorded separately from presence

> An annotation is not automatically experimental evidence. These are the ECO evidence categories attached to the signal-peptide features themselves.

**Signal-peptide feature-evidence rows** (2680 rows over 2507 proteins; a feature may carry several evidence attributions):

| evidence category | n evidence rows | n distinct proteins |
|---|---|---|
| experimental | 642 | 491 |
| curator_inference | 24 | 22 |
| author_statement | 4 | 4 |
| by_similarity | 89 | 89 |
| sequence_analysis_predictor | 1875 | 1875 |
| imported | 9 | 9 |
| unspecified | 37 | 37 |

**Strongest evidence category per protein**, over the 2439 proteins with `annotated_present`:

| strongest category | n proteins | share of annotated | Extracellular | Cytoplasm |
|---|---|---|---|---|
| experimental | 479 | 19.6% | 468 | 11 |
| curator_inference | 21 | 0.9% | 21 | 0 |
| author_statement | 4 | 0.2% | 4 | 0 |
| by_similarity | 85 | 3.5% | 83 | 2 |
| sequence_analysis_predictor | 1820 | 74.6% | 1714 | 106 |
| unspecified | 30 | 1.2% | 30 | 0 |

**Experimental-evidence subset: 479 proteins** (19.6% of annotated, 4.7% of the cohort). This is the denominator for any sensitivity analysis restricted to experimentally-evidenced signal peptides.

**Transit peptides:** 73 proteins carry an annotated transit peptide (106 evidence rows). Evidence breakdown: sequence_analysis_predictor 73, experimental 17, by_similarity 7, combinatorial_experimental 4, curator_inference 3, unspecified 2.

---

## 4. Family and domain coverage

> **A Pfam domain is not a whole-protein family.** Pfam and InterPro entries are region-level annotations: a multidomain protein carries several, and two proteins sharing one domain are not thereby in the same family. The counts below are *annotation coverage*, not a family assignment and not a clustering.

> **No homology clusters have been computed.** No MMseqs2 or CD-HIT run has been performed, so no sequence-similarity groups exist for this cohort. The residual-homology question of §2.6 of the specification remains unmeasured.

| quantity | value |
|---|---|
| proteins with ≥1 Pfam cross-reference | 9454 (93.2%) |
| proteins with ≥1 InterPro cross-reference | 9681 (95.5%) |
| proteins with a UniProt SIMILARITY (family) comment | 7680 (75.7%) |
| distinct Pfam accessions observed | 4828 |
| distinct InterPro accessions observed | 10239 |
| total Pfam cross-reference rows | 17118 |
| total InterPro cross-reference rows | 45264 |
| max Pfam domains on one protein | 15 |
| proteins with >1 Pfam domain | 4224 |

| role | n | ≥1 Pfam | ≥1 InterPro | SIMILARITY comment |
|---|---|---|---|---|
| discovery | 5201 | 93.5% | 95.5% | 76.2% |
| tuning | 2556 | 92.8% | 94.9% | 77.3% |
| evaluation | 2382 | 93.2% | 96.1% | 73.0% |

The SIMILARITY comment is the closest thing UniProt gives to a whole-protein family statement, and it is free text, not a controlled identifier. It is reported as its own column rather than being conflated with the domain counts.

---

## 5. Organism resolution (for the optional §6 external-evaluation arm)

| organism type | n | share |
|---|---|---|
| animal | 6433 | 63.4% |
| fungi | 2335 | 23.0% |
| plant | 1028 | 10.1% |
| unsupported | 342 | 3.4% |
| unresolved | 1 | 0.0% |

**DeepLoc's "Other" kingdom bucket resolved by taxonomy** (342 proteins). This bucket was *not* assumed unsupported; each lineage was checked:

| organism type | n |
|---|---|
| unsupported | 342 |

Every one of these resolved to a **eukaryote outside Metazoa, Viridiplantae and Fungi** — so the DeepLoc bucket and the predictor-support verdict happen to coincide here, but only because the lineages were checked. What they actually are, by major eukaryotic group (second lineage rank):

| eukaryotic group | n proteins | example organism |
|---|---|---|
| Sar | 221 | Phytophthora sojae (strain P6497) |
| Amoebozoa | 98 | Entamoeba histolytica (strain ATCC 30459 / HM-1:IMSS / ABRM) |
| Discoba | 22 | Trypanosoma brucei brucei (strain 927/4 GUTat10.1) |
| Metamonada | 1 | Giardia intestinalis (strain ATCC 50803 / WB clone C6) |
| (unresolved) | 1 |  |

Unsupported and unresolved records are retained in every output table with their reason in `wolfpsort_organism_reason`, excluded only from a WoLF PSORT submission, and never dropped from the cohort or from any other analysis.

---

## 6. Integrity checks

| check | result |
|---|---|
| primary rows | 10139 (expected 10139) |
| unique accessions | 10139 |
| exactly one row per development protein | PASS |
| labels and partitions unchanged vs manifest | PASS |
| frozen sequences re-verified against manifest sha256 | 10139 match, 0 mismatch, 0 missing |
| P0 rows requested as annotation targets | 0 |
| HPA files read | 0 |
| every row has an explicit status (no silent drops) | PASS |
| engineering proteins flagged | 8 |
| role counts | discovery 5201, evaluation 2382, tuning 2556 |

### P0 boundary — one case that needs stating plainly

No P0 row was requested as an annotation target, and no P0 row has an annotation row in any output table. However, **11 development isoform accessions have a parent accession that is itself a P0 row**, and **8 more share a base accession with a P0 isoform**. Retrieving the parent entry is unavoidable if those development isoforms are to be annotated at all.

Parent accessions that are themselves P0 rows: `O35426`, `O75410`, `O95278`, `O95831`, `P60484`, `Q12955`, `Q6NYC8`, `Q9BXL6`, `Q9HDD0`, `Q9QZU4`, `Q9UMS6`

Every affected development row carries `p0_base_collision != none`, so these rows can be excluded from any later analysis without refetching, if the reviewer prefers that. **This is flagged for a decision, not decided here.** The P0 reservation as written covers fitting, feature selection, tuning and scoring, none of which this step performs.

---

## 7. What these annotations support, and what they do not

*Assessment of analysis readiness. The experiment is not redesigned here.*

### Supported

- **Signal-peptide co-occurrence characterisation (§1.2 stage 2).** 2439 proteins carry an annotation (24.1% of the cohort): 2320/2498 Extracellular (92.9%) against 119/7641 Cytoplasm (1.6%). Both arms are populated, so the characterisation can run.

- **Per-residue interval analysis** on the 9758 records whose coordinates are transferable — the annotated intervals can be aligned to the exact extraction-input sequence without re-indexing.

- **Domain-level covariates.** Pfam/InterPro coverage is high enough to use as a descriptive covariate in feature characterisation.

- **Organism typing** for the optional external-evaluation arm, resolved from taxonomy rather than from the DeepLoc kingdom bucket.

- **Entry dates** for the exposure bound described in §6.3 of the specification.


### Not supported

- **Any analysis treating `not_annotated` as a biological negative.** The 7319 `not_annotated` proteins are an unlabelled pool, not a negative class.

- **An experimentally-evidenced-only analysis at full scale.** Only 479 proteins have experimental or combinatorial-experimental evidence (4.7% of the cohort), so that sensitivity analysis is small and must report its denominator.

- **Family-level control or family-disjointness claims.** Pfam/InterPro are region annotations, not family assignments, and no clustering has been run. The residual-homology question stays open.

- **Cluster-level uncertainty (§5.6 cluster bootstrap).** Blocked until a clustering run exists. Until then a protein-level bootstrap assumes an independence that remains unverified.

- **Treating any annotation as independent of a predictor.** A large share of signal-peptide annotations carry `sequence_analysis_predictor` (ECO:0000255) evidence, i.e. they are themselves predictor-derived. Validating a feature that predicts signal peptides against predictor-derived annotations is partly circular, and any such comparison must be reported split by evidence category.

