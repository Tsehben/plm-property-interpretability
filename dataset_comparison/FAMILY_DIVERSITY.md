# Pfam diversity and within-family feasibility

Pfam 37.0 / PyHMMER 0.11.2 annotation is complete for every unique input sequence, with existing UESolDS hits reused by exact SHA256. Full-length source sequences are preserved; NetSolP's already-trimmed constructs are annotated as delivered. All reported family statistics are descriptive, not selected by predictive performance.

## Diversity and independent support

| policy | total | annotated | annotation_percent | families | clans | groups | mixed_families | families_proteins_ge3 | families_proteins_ge5 | families_groups_ge3 | families_groups_ge5 | families_exclusive_groups_ge3 | families_exclusive_groups_ge5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PSI_original | 12216 | 11273 | 92.28 | 2093 | 316 | 3530 | 996 | 415 | 220 | 187 | 77 | 71 | 30 |
| PSI_NetSolP | 11226 | 10385 | 92.51 | 2075 | 314 | 3423 | 962 | 394 | 208 | 181 | 71 | 66 | 32 |
| NESG_NetSolP | 1323 | 1169 | 88.36 | 852 | 191 | 985 | 157 | 25 | 6 | 13 | 2 | 7 | 2 |
| UESolDS_DNASU_proxy | 460 | 400 | 86.96 | 405 | 109 | 424 | 30 | 0 | 0 | 0 | 0 | 0 | 0 |
| NESG_raw_scores | 9703 | 8971 | 92.46 | 2101 | 328 | 3734 | — | — | — | — | — | — | — |
| NESG_unique_scores | 9272 | 8560 | 92.32 | 2086 | 328 | 3697 | — | — | — | — | — | — | — |
| NESG_unique_E45_S45_vs0 | 3762 | 3510 | 93.30 | 1383 | 243 | 1820 | 390 | 120 | 64 | 55 | 28 | 35 | 13 |

Family and clan counts are distinct concepts. “Groups” are connected components of direct MMseqs2 hits at ≥30% identity, ≥80% coverage of each sequence and E ≤ 0.001 across the entire audit union. They are a reproducible operational homology-control criterion, not proof of evolutionary independence. Added sequences can bridge components; therefore UESolDS audit-union group counts differ from its original frozen cluster IDs, which are preserved separately and unchanged.

A protein contributes once to each distinct Pfam accession. Multiple domains may put it in multiple families; family prevalences therefore do not sum to 100%. Counts of proteins in qualifying families use the union of record IDs. ≥3/≥5 **group** thresholds require that many groups containing each class; a mixed-label group can contribute to both class counts. The stricter `exclusive_groups` columns require class-exclusive groups on each side, so the same group cannot supply both classes. Pfam sharing does not establish identical full-protein architecture.

| policy | multifamily_proteins | proteins_in_groups_ge3 | proteins_in_groups_ge5 |
| --- | --- | --- | --- |
| PSI_original | 3634 | 4205 | 2649 |
| PSI_NetSolP | 3495 | 3793 | 2367 |
| NESG_NetSolP | 300 | 95 | 15 |
| UESolDS_DNASU_proxy | 105 | 0 | 0 |
| NESG_raw_scores | 3228 | — | — |
| NESG_unique_scores | 3120 | — | — |
| NESG_unique_E45_S45_vs0 | 1337 | 883 | 698 |

## Annotation ambiguity sensitivity

| policy | filter | families | mixed_families | families_groups_ge3 | families_groups_ge5 |
| --- | --- | --- | --- | --- | --- |
| PSI_NetSolP | all_GA | 2075 | 962 | 181 | 71 |
| PSI_NetSolP | no_competing_overlap | 1664 | 692 | 75 | 21 |
| PSI_NetSolP | HMM_coverage80 | 1733 | 754 | 117 | 40 |
| NESG_NetSolP | all_GA | 852 | 157 | 13 | 2 |
| NESG_NetSolP | no_competing_overlap | 645 | 86 | 2 | 0 |
| NESG_NetSolP | HMM_coverage80 | 685 | 104 | 6 | 2 |
| NESG_unique_E45_S45_vs0 | all_GA | 1383 | 390 | 55 | 28 |
| NESG_unique_E45_S45_vs0 | no_competing_overlap | 1079 | 245 | 16 | 3 |
| NESG_unique_E45_S45_vs0 | HMM_coverage80 | 1168 | 282 | 27 | 11 |

`no_competing_overlap` removes both hits when different Pfam families from the same clan overlap in aligned coordinates. It is intentionally conservative, not official PfamScan adjudication. `HMM_coverage80` requires ≥80% of the profile aligned. Neither removes a protein from heterogeneous analysis A. The substantial reduction warns against interpreting multiple related overlapping Pfams as multiple independent biological contrasts. All qualifying-family identities and counts are in `family_support.tsv.gz`.

## Homology-aware evaluation feasibility

The prospective diagnostic split assigns each union component to one of five folds by SHA256(`dataset-audit-v1:` + lexicographically smallest sequence hash) modulo 5. Fold 0 is the prespecified evaluation fold; others are development. It is label-blind, without family or class balancing. This is a newly audited design, **not a pristine historical external holdout**. Original UESolDS roles and published NetSolP partitions remain intact in the canonical records.

| policy | train_groups_ge3_eval_both | train_groups_ge5_eval_both | train_groups_ge3_eval_ge3 | train_groups_ge5_eval_ge5 |
| --- | --- | --- | --- | --- |
| PSI_original | 61 | 39 | 18 | 4 |
| PSI_NetSolP | 58 | 40 | 17 | 3 |
| NESG_NetSolP | 2 | 1 | 0 | 0 |
| UESolDS_DNASU_proxy | 0 | 0 | 0 | 0 |
| NESG_raw_scores | — | — | — | — |
| NESG_unique_scores | — | — | — | — |
| NESG_unique_E45_S45_vs0 | 21 | 15 | 5 | 0 |

The first two columns require ≥3/≥5 development groups per class and at least one evaluation group per class in the same family. The last two require ≥3/≥5 in evaluation as well. A single group per evaluation class establishes minimal evaluability, not a well-powered family-specific estimate. `five_fold_feasibility.tsv` reports every fold rather than selecting the most favorable. Published-PSI fold support is separately recorded in `published_PSI_family_feasibility.tsv`.

For strict NESG, A is biologically heterogeneous and B is feasible for a limited panel; after conservative Pfam overlap removal, only 16 families meet ≥3 groups/class and 3 meet ≥5. NetSolP PSI supports a broader B panel (75 and 21 under the same conservative rule). NetSolP NESG alone is weak for a primary controlled study (2 and 0 after that rule).

Within NetSolP PSI, pET15 has 346 mixed-label families and 48/11 with ≥3/≥5 groups per class; pET21 has 677 mixed families and 116/48 meeting those thresholds. These are same-vector strata, not guaranteed matched experimental batches.

PSI's two vector strata are recorded in `PSI_family_vector_support.tsv`; this permits same-vector contrasts without flattening natural family proportions. NESG's common study pipeline gives better study-level comparability, but complete record-level vector, condition and replicate metadata are absent. Family agreement cannot repair that gap.

## Actual named biological composition

The figures `figures/top20_family.png` and `figures/top20_clan.png` show actual biological names. `top20_named_families_clans.tsv` includes top-20 combined ranks and top-20 ranks within each class (ties included), with exact counts, denominators and percentages. Columns `_0` are insoluble and `_1` soluble. Families absent in one class have zero prevalence there.

### PSI_original: family

| accession | name | description | N_0 | N_1 | percent_0 | percent_1 |
| --- | --- | --- | --- | --- | --- | --- |
| PF13450 | NAD_binding_8 | NAD(P)-binding Rossmann-like domain | 73 | 91 | 1.84 | 1.10 |
| PF08241 | Methyltransf_11 | Methyltransferase domain | 58 | 102 | 1.46 | 1.24 |
| PF03816 | LytR_cpsA_psr | LytR_cpsA_psr family | 25 | 134 | 0.63 | 1.63 |
| PF00535 | Glycos_transf_2 | Glycosyl transferase family 2 | 47 | 111 | 1.18 | 1.35 |
| PF13649 | Methyltransf_25 | Methyltransferase domain | 50 | 101 | 1.26 | 1.23 |
| PF00069 | Pkinase | Protein kinase domain | 117 | 30 | 2.94 | 0.36 |
| PF01266 | DAO | FAD dependent oxidoreductase | 75 | 68 | 1.89 | 0.83 |
| PF04055 | Radical_SAM | Radical SAM superfamily | 70 | 72 | 1.76 | 0.87 |
| PF07714 | PK_Tyr_Ser-Thr | Protein tyrosine and serine/threonine kinase | 113 | 29 | 2.84 | 0.35 |
| PF13489 | Methyltransf_23 | Methyltransferase domain | 51 | 83 | 1.28 | 1.01 |

### PSI_original: clan

| accession | name | description | N_0 | N_1 | percent_0 | percent_1 |
| --- | --- | --- | --- | --- | --- | --- |
| CL0063 | NADP_Rossmann | FAD/NAD(P)-binding Rossmann fold Superfamily | 306 | 540 | 7.69 | 6.55 |
| CL0023 | P-loop_NTPase | P-loop containing nucleoside triphosphate hydrolase superfamily | 245 | 385 | 6.16 | 4.67 |
| CL0123 | HTH | Helix-turn-helix clan | 106 | 355 | 2.66 | 4.31 |
| CL0021 | OB | OB fold | 44 | 258 | 1.11 | 3.13 |
| CL0036 | TIM_barrel | Common phosphate binding-site TIM barrel superfamily | 150 | 147 | 3.77 | 1.78 |
| CL0110 | GT-A | Glycosyl transferase clan GT-A | 86 | 177 | 2.16 | 2.15 |
| CL0072 | Ubiquitin | Ubiquitin superfamily | 60 | 203 | 1.51 | 2.46 |
| CL0016 | PKinase | Protein kinase superfamily | 167 | 62 | 4.20 | 0.75 |
| CL0186 | Beta_propeller | Beta propeller clan | 99 | 101 | 2.49 | 1.23 |
| CL0020 | TPR | Tetratrico peptide repeat superfamily | 58 | 130 | 1.46 | 1.58 |

### PSI_NetSolP: family

| accession | name | description | N_0 | N_1 | percent_0 | percent_1 |
| --- | --- | --- | --- | --- | --- | --- |
| PF13450 | NAD_binding_8 | NAD(P)-binding Rossmann-like domain | 73 | 91 | 1.96 | 1.21 |
| PF03816 | LytR_cpsA_psr | LytR_cpsA_psr family | 25 | 134 | 0.67 | 1.79 |
| PF00535 | Glycos_transf_2 | Glycosyl transferase family 2 | 47 | 111 | 1.26 | 1.48 |
| PF08241 | Methyltransf_11 | Methyltransferase domain | 50 | 95 | 1.34 | 1.27 |
| PF00069 | Pkinase | Protein kinase domain | 116 | 29 | 3.11 | 0.39 |
| PF01266 | DAO | FAD dependent oxidoreductase | 74 | 68 | 1.99 | 0.91 |
| PF07714 | PK_Tyr_Ser-Thr | Protein tyrosine and serine/threonine kinase | 112 | 28 | 3.01 | 0.37 |
| PF04055 | Radical_SAM | Radical SAM superfamily | 70 | 70 | 1.88 | 0.93 |
| PF13649 | Methyltransf_25 | Methyltransferase domain | 45 | 94 | 1.21 | 1.25 |
| PF13489 | Methyltransf_23 | Methyltransferase domain | 44 | 78 | 1.18 | 1.04 |

### PSI_NetSolP: clan

| accession | name | description | N_0 | N_1 | percent_0 | percent_1 |
| --- | --- | --- | --- | --- | --- | --- |
| CL0063 | NADP_Rossmann | FAD/NAD(P)-binding Rossmann fold Superfamily | 293 | 523 | 7.86 | 6.97 |
| CL0023 | P-loop_NTPase | P-loop containing nucleoside triphosphate hydrolase superfamily | 235 | 367 | 6.31 | 4.89 |
| CL0123 | HTH | Helix-turn-helix clan | 97 | 320 | 2.60 | 4.27 |
| CL0036 | TIM_barrel | Common phosphate binding-site TIM barrel superfamily | 146 | 142 | 3.92 | 1.89 |
| CL0021 | OB | OB fold | 44 | 240 | 1.18 | 3.20 |
| CL0110 | GT-A | Glycosyl transferase clan GT-A | 86 | 174 | 2.31 | 2.32 |
| CL0072 | Ubiquitin | Ubiquitin superfamily | 57 | 195 | 1.53 | 2.60 |
| CL0016 | PKinase | Protein kinase superfamily | 166 | 61 | 4.46 | 0.81 |
| CL0186 | Beta_propeller | Beta propeller clan | 97 | 99 | 2.60 | 1.32 |
| CL0010 | SH3 | Src homology-3 domain | 20 | 154 | 0.54 | 2.05 |

### NESG_unique_E45_S45_vs0: family

| accession | name | description | N_0 | N_1 | percent_0 | percent_1 |
| --- | --- | --- | --- | --- | --- | --- |
| PF13450 | NAD_binding_8 | NAD(P)-binding Rossmann-like domain | 76 | 43 | 6.62 | 1.64 |
| PF01266 | DAO | FAD dependent oxidoreductase | 70 | 25 | 6.10 | 0.96 |
| PF01494 | FAD_binding_3 | FAD binding domain | 57 | 35 | 4.97 | 1.34 |
| PF00535 | Glycos_transf_2 | Glycosyl transferase family 2 | 18 | 59 | 1.57 | 2.26 |
| PF00890 | FAD_binding_2 | FAD binding domain | 51 | 25 | 4.44 | 0.96 |
| PF00575 | S1 | S1 RNA binding domain | 3 | 61 | 0.26 | 2.33 |
| PF08220 | HTH_DeoR | DeoR-like helix-turn-helix domain | 25 | 32 | 2.18 | 1.22 |
| PF00455 | DeoRC | DeoR C terminal sensor domain | 26 | 31 | 2.26 | 1.19 |
| PF01656 | CbiA | CobQ/CobB/MinD/ParA nucleotide binding domain | 19 | 38 | 1.66 | 1.45 |
| PF04055 | Radical_SAM | Radical SAM superfamily | 22 | 33 | 1.92 | 1.26 |

### NESG_unique_E45_S45_vs0: clan

| accession | name | description | N_0 | N_1 | percent_0 | percent_1 |
| --- | --- | --- | --- | --- | --- | --- |
| CL0063 | NADP_Rossmann | FAD/NAD(P)-binding Rossmann fold Superfamily | 228 | 222 | 19.86 | 8.49 |
| CL0023 | P-loop_NTPase | P-loop containing nucleoside triphosphate hydrolase superfamily | 121 | 213 | 10.54 | 8.15 |
| CL0021 | OB | OB fold | 22 | 150 | 1.92 | 5.74 |
| CL0123 | HTH | Helix-turn-helix clan | 32 | 110 | 2.79 | 4.21 |
| CL0110 | GT-A | Glycosyl transferase clan GT-A | 32 | 87 | 2.79 | 3.33 |
| CL0036 | TIM_barrel | Common phosphate binding-site TIM barrel superfamily | 39 | 67 | 3.40 | 2.56 |
| CL0105 | Hybrid | Barrel sandwich hybrid superfamily | 10 | 70 | 0.87 | 2.68 |
| CL0167 | Zn_Beta_Ribbon | Zinc beta-ribbon | 16 | 62 | 1.39 | 2.37 |
| CL0050 | HotDog | HotDog superfamily | 22 | 44 | 1.92 | 1.68 |
| CL0219 | RNase_H | Ribonuclease H-like superfamily | 2 | 51 | 0.17 | 1.95 |

### NESG_NetSolP: family

| accession | name | description | N_0 | N_1 | percent_0 | percent_1 |
| --- | --- | --- | --- | --- | --- | --- |
| PF00535 | Glycos_transf_2 | Glycosyl transferase family 2 | 5 | 16 | 1.04 | 1.90 |
| PF00696 | AA_kinase | Amino acid kinase family | 5 | 11 | 1.04 | 1.31 |
| PF20143 | NAD_kinase_C | ATP-NAD kinase C-terminal domain | 5 | 11 | 1.04 | 1.31 |
| PF08241 | Methyltransf_11 | Methyltransferase domain | 8 | 7 | 1.66 | 0.83 |
| PF13649 | Methyltransf_25 | Methyltransferase domain | 6 | 7 | 1.25 | 0.83 |
| PF04909 | Amidohydro_2 | Amidohydrolase | 10 | 3 | 2.08 | 0.36 |
| PF13450 | NAD_binding_8 | NAD(P)-binding Rossmann-like domain | 9 | 4 | 1.87 | 0.48 |
| PF01513 | NAD_kinase | ATP-NAD kinase N-terminal domain | 5 | 8 | 1.04 | 0.95 |
| PF13641 | Glyco_tranf_2_3 | Glycosyltransferase like family 2 | 4 | 9 | 0.83 | 1.07 |
| PF00265 | TK | Thymidine kinase | 0 | 12 | 0 | 1.43 |

### NESG_NetSolP: clan

| accession | name | description | N_0 | N_1 | percent_0 | percent_1 |
| --- | --- | --- | --- | --- | --- | --- |
| CL0063 | NADP_Rossmann | FAD/NAD(P)-binding Rossmann fold Superfamily | 42 | 48 | 8.73 | 5.70 |
| CL0036 | TIM_barrel | Common phosphate binding-site TIM barrel superfamily | 24 | 28 | 4.99 | 3.33 |
| CL0023 | P-loop_NTPase | P-loop containing nucleoside triphosphate hydrolase superfamily | 16 | 36 | 3.33 | 4.28 |
| CL0110 | GT-A | Glycosyl transferase clan GT-A | 7 | 23 | 1.46 | 2.73 |
| CL0123 | HTH | Helix-turn-helix clan | 7 | 17 | 1.46 | 2.02 |
| CL0029 | Cupin | Cupin fold | 6 | 17 | 1.25 | 2.02 |
| CL0105 | Hybrid | Barrel sandwich hybrid superfamily | 3 | 17 | 0.62 | 2.02 |
| CL0701 | NAD_kinase_C | NAD-kinase beta sandwich domain superfamily | 8 | 11 | 1.66 | 1.31 |
| CL0240 | PFK | PFK-like superfamily | 8 | 8 | 1.66 | 0.95 |
| CL0050 | HotDog | HotDog superfamily | 6 | 10 | 1.25 | 1.19 |

The composition differs descriptively: in strict NESG, NAD(P)-binding Rossmann-like PF13450 occurs in 76/1,148 insoluble proteins (6.62%) versus 43/2,614 soluble proteins (1.64%); clan CL0063 occurs in 19.86% versus 8.49%. In NetSolP PSI, protein kinase PF00069 occurs in 116/3,726 insoluble proteins (3.11%) versus 29/7,500 soluble proteins (0.39%). These natural biological associations are retained, not normalized away. Related overlapping hits and broader architectures must be considered before attributing an effect to an individual family.

## Long tail and multidomain accounting

| policy | families | singleton_families | families_at_most5 | largest_family |
| --- | --- | --- | --- | --- |
| PSI_original | 2093 | 472 | 1288 | 164 |
| PSI_NetSolP | 2075 | 489 | 1312 | 164 |
| NESG_raw_scores | 2101 | 894 | 1537 | 301 |
| NESG_unique_scores | 2086 | 900 | 1541 | 297 |
| NESG_unique_E45_S45_vs0 | 1383 | 640 | 1132 | 119 |
| NESG_NetSolP | 852 | 556 | 793 | 21 |
| UESolDS_DNASU_proxy | 405 | 311 | 398 | 12 |

`representative_distributions.tsv` provides an optional mutually exclusive 100% distribution: choose the longest aligned domain, then greatest domain bit score, accession and coordinate as deterministic tie-breaks. This is label-independent but loses multidomain information; the main prevalence analysis uses all memberships. No-hit and annotated-without-clan categories are separate. `family_long_tail.png` displays the full ranked frequency distribution.

Length, all 20 amino-acid fractions, charge proxy and mean Kyte–Doolittle hydropathy are summarized by class in `biochemical_summary.tsv`. Unknown-residue sequences retain missing physicochemical summaries rather than silently imputing amino acids. These descriptions are not trained predictors or mechanistic conclusions.
