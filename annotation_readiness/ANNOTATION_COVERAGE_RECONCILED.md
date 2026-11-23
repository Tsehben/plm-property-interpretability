# Annotation-readiness coverage — reconciled

**This file does not replace `ANNOTATION_COVERAGE_REPORT.md`.** That report is preserved
unchanged. This one records the corrected figures and the corrected wording, produced by
`reconcile_coverage.py` from the saved tables only: no network requests, no model runs, no
annotation expansion, no clustering, and no change to `out/` or to the primary cohort.

Inputs read: `out/annotation_primary.tsv`, `out/annotation_features.tsv`,
`out/annotation_domains.tsv`, `DeepLoc2/baseline/deeploc_binary_manifest.tsv`.
Corrected outputs: `out_reconciled/`.

---

## 1. Reconciled coverage table

Counts are **proteins**, one row per development protein. All three signal-peptide states
are shown, so each row sums to its own denominator and nothing is hidden inside a
percentage.

| role | class | n (denominator) | annotated_present | not_annotated | unverifiable | of which experimental evidence |
|---|---|---|---|---|---|---|
| discovery | Extracellular | 1339 | 1254 (93.7%) | 53 (4.0%) | 32 (2.4%) | 236 (17.6% of 1339) |
| discovery | Cytoplasm | 3862 | 59 (1.5%) | 3648 (94.5%) | 155 (4.0%) | 4 (0.1%) |
| tuning | Extracellular | 584 | 551 (94.3%) | 22 (3.8%) | 11 (1.9%) | 94 (16.1%) |
| tuning | Cytoplasm | 1972 | 31 (1.6%) | 1851 (93.9%) | 90 (4.6%) | 3 (0.2%) |
| evaluation | Extracellular | 575 | 515 (89.6%) | 41 (7.1%) | 19 (3.3%) | 138 (24.0%) |
| evaluation | Cytoplasm | 1807 | 29 (1.6%) | 1704 (94.3%) | 74 (4.1%) | 4 (0.2%) |
| **all** | **Extracellular** | **2498** | **2320 (92.9%)** | 116 (4.6%) | 62 (2.5%) | 468 (18.7%) |
| **all** | **Cytoplasm** | **7641** | **119 (1.6%)** | 7203 (94.3%) | 319 (4.2%) | 11 (0.1%) |
| **all** | **all** | **10139** | **2439 (24.1%)** | 7319 (72.2%) | 381 (3.8%) | 479 (4.7%) |

Machine-readable: `out_reconciled/reconciled_coverage.tsv`.

---

## 2. Corrected passages

### 2.1 Organism categories — §5 of the original report

The original report's eukaryotic-group table is headed as a breakdown of the 342
DeepLoc-"Other" proteins but appends a row `(unresolved) | 1`, giving 221 + 98 + 22 + 1 + 1
= 343 against a denominator of 342. The unresolved record is **not** one of the 342. It is
`Q9Y6N3`, DeepLoc kingdom *Metazoa*, whose UniProt entry is deleted and for which no lineage
was retrieved at all. "Resolved to an organism outside WoLF PSORT's supported groups" and
"taxonomy not resolved" are different states and must not appear in one series.

**Replaces the two organism tables in §5:**

> **Organism categories, mutually exclusive, over all 10,139 development proteins.** A
> protein is *resolved* if a UniProt lineage was retrieved for it; resolved proteins are
> then either inside or outside the three organism types WoLF PSORT supports.
>
> | category | n | example |
> |---|---|---|
> | resolved, supported — animal | 6433 | *Bos taurus* |
> | resolved, supported — fungi | 2335 | *Trichophyton rubrum* |
> | resolved, supported — plant | 1028 | *Oryza sativa* subsp. *japonica* |
> | resolved, **outside** the supported groups — Sar | 221 | *Phytophthora sojae* (strain P6497) |
> | resolved, **outside** the supported groups — Amoebozoa | 98 | *Entamoeba histolytica* (ATCC 30459 / HM-1:IMSS / ABRM) |
> | resolved, **outside** the supported groups — Discoba | 22 | *Trypanosoma brucei brucei* (927/4 GUTat10.1) |
> | resolved, **outside** the supported groups — Metamonada | 1 | *Giardia intestinalis* (ATCC 50803 / WB clone C6) |
> | **taxonomy unresolved** (no UniProt record retrieved) | 1 | `Q9Y6N3` |
> | **total** | **10139** | |
>
> The 342 outside-the-supported-groups proteins are exactly DeepLoc's "Other" kingdom
> bucket: 221 + 98 + 22 + 1 = 342, every one resolved. The bucket and the support verdict
> coincide, but only because each lineage was checked. The single unresolved protein is
> Metazoa by DeepLoc kingdom and sits outside that bucket; for it, organism support is not
> *unsupported* — it is *unknown*, and is recorded as unknown.
>
> Outside-group and unresolved records are retained in every output table with their reason
> in `wolfpsort_organism_reason`. They are excluded only from a WoLF PSORT submission,
> never from the cohort or from any other analysis.

Table: `out_reconciled/organism_categories.tsv`.

### 2.2 Status accounting — §1 of the original report

"**Failures, itemised** (122 rows)" mixes dimensions. 122 is the **union** of *accession not
current* (6) and *sequence not exact* (120), which overlap in 4 rows; it is not a sum, and
it is not a count of failed annotations. Two of the 122 are merged accessions whose
retrieved sequence is SHA-256 identical to the frozen input — those are provenance caveats,
and their annotations are usable.

**Replaces the "Failures, itemised" block:**

> **Four status dimensions, each one exhaustive over all 10,139 proteins.** They are
> reported separately and never added together: a protein that is non-exact in dimension B
> is non-transferable in dimension C for the *same* reason, so summing across dimensions
> double-counts.
>
> | dimension | category | n | % of 10,139 |
> |---|---|---|---|
> | **A. Retrieval / accession status** | resolved as requested (current accession) | 10131 | 99.9% |
> | | resolved via merged/secondary accession | 5 | 0.0% |
> | | active entry, isoform sequence unresolved (HTTP 404) | 2 | 0.0% |
> | | inactive, deleted — no record | 1 | 0.0% |
> | **B. Agreement with frozen input sequence** | exact (SHA-256 identical) | 10019 | 98.8% |
> | | mismatch (entry re-sequenced since the DeepLoc snapshot) | 117 | 1.2% |
> | | unverifiable (no sequence retrieved) | 3 | 0.0% |
> | **C. Residue-coordinate transferability** | transferable | 9758 | 96.2% |
> | | not transferable: isoform matches frozen input, differs from canonical | 261 | 2.6% |
> | | not transferable: retrieved sequence differs from frozen input | 117 | 1.2% |
> | | not transferable: no sequence retrieved | 3 | 0.0% |
> | **D. Signal-peptide status** | annotated_present | 2439 | 24.1% |
> | | not_annotated (**not** a confirmed biological negative) | 7319 | 72.2% |
> | | unverifiable | 381 | 3.8% |
>
> **How the dimensions overlap.** 6 proteins have a non-current accession and 120 have a
> non-exact sequence; 4 proteins are in both, so the union is 122, not 126. Reclassifying
> that union by *annotation outcome* rather than by status label:
>
> | outcome | n | accessions |
> |---|---|---|
> | usable, with a provenance caveat — merged accession, sequence identical to the frozen input | 2 | `A0A140G945` (→ P02489), `G0Y276` (→ G2WWV6) |
> | annotation not transferable — retrieved sequence differs from the frozen input | 117 | `out_reconciled/non_exact_records_reclassified.tsv` |
> | no annotation available — no sequence retrieved | 3 | `P22455-3`, `Q11179-3` (isoform 404), `Q9Y6N3` (deleted) |
>
> So **3** proteins have no annotation at all, **117** have an entry but no transferable
> residue coordinates, and **2** are fully annotated under a different primary accession.

Tables: `out_reconciled/status_dimensions.tsv`,
`out_reconciled/non_exact_records_reclassified.tsv`.

### 2.3 What is counted, and the 261 isoforms — §2–3 of the original report

**Added before the coverage table:**

> **These are protein counts.** 2,439 is the number of development proteins whose signal
> status is `annotated_present`; each carries exactly one signal feature. Three other
> numbers appear in these outputs and are not interchangeable with it: **2,507** signal
> *feature records* across all proteins (2,439 on annotated proteins plus 68 on proteins
> whose status is `unverifiable`), **2,680** rows in `annotation_features.tsv`, which holds
> one row per *(feature, evidence)* pair, and **2,498** Extracellular proteins, a cohort
> denominator. The original report's evidence table uses 2,507 as its denominator and its
> strongest-evidence table uses 2,439; both are correct for what they count, and they are
> labelled here so the two are not read as the same population.

**Replaces the paragraph on the 261 isoforms:**

> **The 261 exact-matching isoforms are classified `unverifiable`, not `not_annotated`.**
> Their own sequence is SHA-256 identical to the frozen extraction input, but UniProt's
> feature coordinates are defined on the canonical sequence, which differs. We therefore
> know neither that a signal peptide is present at a stated interval in *this* sequence nor
> that one is absent. `unverifiable` is the only state consistent with that, and the
> distinction matters because **50** of the 261 have a parent entry that does carry a
> signal-peptide annotation (`parent_entry_has_signal = 1`), held back rather than
> transferred. Classifying any of them as `not_annotated` would convert an unresolved
> coordinate question into a claimed biological negative. By role and class:
>
> | role | class | n | parent entry has a signal annotation |
> |---|---|---|---|
> | discovery | Extracellular | 19 | 17 |
> | discovery | Cytoplasm | 98 | 5 |
> | tuning | Extracellular | 7 | 4 |
> | tuning | Cytoplasm | 62 | 6 |
> | evaluation | Extracellular | 17 | 16 |
> | evaluation | Cytoplasm | 58 | 2 |
>
> The other 120 unverifiable proteins are the 117 sequence mismatches and the 3 with no
> sequence; the full breakdown is `out_reconciled/unverifiable_breakdown.tsv`. Recovering
> the 261 requires aligning each isoform to its canonical sequence and mapping the interval
> across — a bounded follow-up, not performed here.

### 2.4 Evidence categories and tool attribution — §3 of the original report

**Added to §3:**

> **Evidence categories overlap; the strongest-evidence collapse does not.** A signal
> feature may carry several evidence attributions, and **20** annotated proteins carry more
> than one category, so the overlapping set sizes sum to more than 2,439. Both views are
> reported and must not be mixed:
>
> | evidence category | overlapping set membership (proteins) | mutually exclusive, strongest per protein |
> |---|---|---|
> | sequence_analysis_predictor | 1830 (75.0%) | 1820 (74.6%) |
> | experimental | 479 (19.6%) | 479 (19.6%) |
> | by_similarity | 85 (3.5%) | 85 (3.5%) |
> | unspecified | 30 (1.2%) | 30 (1.2%) |
> | curator_inference | 22 (0.9%) | 21 (0.9%) |
> | imported | 9 (0.4%) | — (never strongest) |
> | author_statement | 4 (0.2%) | 4 (0.2%) |
>
> Percentages are of the 2,439 annotated proteins. Only the right-hand column partitions
> them.
>
> **An evidence code is not a tool.** ECO:0000255 is "match to sequence model evidence used
> in manual assertion". It names a *kind* of evidence, not a program. In this cohort the
> evidence source field is empty for **1,873 of the 1,875** ECO:0000255 signal-peptide
> evidence rows; the remaining 2 name HAMAP-Rule (MF_03060, MF_03108). No signal-peptide
> annotation in this cohort identifies SignalP or any other named predictor. The generating
> method is therefore **unidentified**, and must not be reported as any particular tool.

**Replaces the "partly circular" bullet in §7:**

> - **Treating a computational annotation as an independent reference.** 1,830 of the 2,439
>   signal-peptide annotations (75.0%) carry ECO:0000255 evidence, and the method behind
>   them is not identified in the record. Computational annotations and a sequence model's
>   features may share sequence-based biases — N-terminal hydrophobicity and charge, in
>   particular — so agreement between them is weaker evidence than agreement with an
>   experimentally evidenced annotation. Whether this amounts to circularity depends on two
>   things that are under our control and must be stated for each analysis: the
>   **provenance** of the annotations used, and the **role** they play. Using them only to
>   *describe* a cohort carries little risk; using them to *select* candidate features and
>   then to *evaluate* those same features is where the dependency becomes a confound. Any
>   comparison between a feature and these annotations is reported split by evidence
>   category, with the 479 experimentally evidenced proteins (19.6% of annotated, 4.7% of
>   the cohort) as the low-dependency subset.

### 2.5 P0-adjacent retrievals — §6 of the original report

The original says "**11** development isoform accessions have a parent accession that is
itself a P0 row, and **8** more share a base accession with a P0 isoform." Those are counts
of *base accessions*, reported as if they were counts of development proteins, and the two
groups are described as if comparable units. The affected development proteins are **27**.
`integrity_checks.json` also reports `p0_rows_requested_as_annotation_targets: 0`, which is
true, but the report should not be read as saying no P0 record was fetched — 11 were.

**Replaces the "P0 boundary" section:**

> ### P0 boundary — stated in full
>
> Three separate facts, none of which implies the others:
>
> 1. **No P0 cohort row was an annotation target.** No P0 accession has a primary row, a
>    feature row or a domain row in any output table. `p0_rows_requested_as_annotation_targets
>    = 0`.
> 2. **11 P0 accessions were nonetheless fetched from UniProt**, as the parent entries of
>    development isoforms: `O35426`, `O75410`, `O95278`, `O95831`, `P60484`, `Q12955`,
>    `Q6NYC8`, `Q9BXL6`, `Q9HDD0`, `Q9QZU4`, `Q9UMS6`. Retrieving the parent entry is
>    unavoidable if the corresponding development isoform is to be annotated at all. It is
>    not correct to describe this step as having made zero P0-associated retrievals.
> 3. **27 unique development proteins are affected**, in two groups that do not overlap:
>
>    | group | development proteins | unique base accessions |
>    |---|---|---|
>    | base accession **is** a P0 cohort row | **15** | 11 |
>    | base accession **shares an entry** with a P0 isoform | **12** | 8 |
>    | **total, unique** | **27** | 19 |
>
>    The 15: `O35426-2`, `O75410-10`, `O95278-4`, `O95278-5`, `O95278-7`, `O95831-4`,
>    `O95831-5`, `P60484-2`, `Q12955-7`, `Q6NYC8-4`, `Q9BXL6-2`, `Q9BXL6-3`, `Q9HDD0-2`,
>    `Q9QZU4-2`, `Q9UMS6-5`. The 12: `F1MAB7`, `F1MAB7-2`, `O00165`, `O00165-3`, `O00165-5`,
>    `P17693-6`, `P30279`, `P55957`, `Q03001`, `Q03001-11`, `Q86W56-4`, `Q8WU39`. In the
>    second group the fetched record is a development protein's own entry; the P0 row is a
>    *different isoform* of it, and no P0 accession was fetched for those 8.
>
> **What P0 information reached development rows.** The DeepLoc P0 **label** was never read:
> the manifest's `target` column was used only to assert that no P0 row is an annotation
> target. However, the entry-level fields on the 15 development rows in group 1 come from a
> record that *is* a P0 cohort row: `organism_*`, `entry_version`, the three entry dates,
> `n_pfam`, `n_interpro`, `pfam_ids`, `interpro_ids`, `family_similarity_text`,
> `parent_entry_has_signal`, and — the one that matters for the holdout —
> `subcellular_location_text`, populated for all 15. That is UniProt subcellular-location
> annotation about a held-out protein appearing on a development row. It is reported here
> rather than claimed absent. No additional P0 record was inspected to establish this.
>
> Every affected row keeps `p0_base_collision != none`. Nothing was deleted, the holdout was
> not altered, and the decision of whether these 27 rows may carry parent-derived entry
> fields is left to the reviewer. The P0 reservation as written covers fitting, feature
> selection, tuning and scoring — none of which this step performs.

Table: `out_reconciled/reconciliation_summary.json`, key `p0_accounting`.

### 2.6 Family and domain coverage

Unchanged and reconfirmed from the saved tables: **9,454** proteins (93.2%) with ≥1 Pfam,
**9,681** (95.5%) with ≥1 InterPro, **7,680** (75.7%) with a SIMILARITY comment; 4,828
distinct Pfam and 10,239 distinct InterPro accessions over **17,118** Pfam and **45,264**
InterPro cross-reference rows; 4,224 proteins with >1 Pfam domain. (Figures quoted in
earlier correspondence as 9,449 / 9,676 / 7,676 / 17,108 / 45,242 / 10,238 / 4,222 were
transcription errors in that correspondence; the saved report and tables have always read
the values above.) By role: `out_reconciled/family_domain_coverage.tsv`. A Pfam or InterPro entry
remains a region-level annotation, not a whole-protein family assignment; no homology
clustering has been computed and no family-disjointness is claimed.

---

## 3. Remaining data-integrity blocker

**One blocker.**

**Entry-level P0 annotation on 15 development rows (§2.5, group 1).** The holdout as
written reserves P0 from fitting, selection, tuning and scoring. Annotation retrieval is
none of those, so nothing in the holdout has been violated yet. But 15 development rows now
carry UniProt entry annotation — including `subcellular_location_text` — sourced from a
record that is itself a P0 cohort row, and any later analysis that conditions on those
fields would be using held-out-protein annotation. This needs a decision before annotations
are used in feature selection, and it is not a decision to make unilaterally. The options
are to accept it and say so, to blank the parent-derived fields on those 15 rows, or to drop
the 15 rows from annotation-conditioned analyses. All three are one filter on
`p0_base_collision` and need no refetching.

**Two items that are not blockers.**

- **50 withheld parent-entry signal annotations** among the 261 exact-matching isoforms
  (§2.3). Bounded, well-defined follow-up: isoform-to-canonical alignment and interval
  mapping. Until then they stay `unverifiable`, which is the correct state.
- **1,830 annotations (75.0%) with unidentified computational evidence** (§2.4). A
  limitation on what the annotations can validate, handled by splitting every
  feature-versus-annotation comparison by evidence category. Not an integrity defect in the
  data.

---

## 4. Scope

No network requests. No model runs. No annotation expansion, no clustering, no change to
the primary cohort, labels, partitions, roles, engineering-protein flags or collision flags.
`out/` is byte-identical to the original run; all corrected outputs are in `out_reconciled/`.
Nothing committed or pushed.
