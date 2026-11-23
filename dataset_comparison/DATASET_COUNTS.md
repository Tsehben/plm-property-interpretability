# Dataset sizes and outcome policies

All numbers come from the pinned downloaded files. `Soluble`/`insoluble` are deliberately blank for raw ordinal-score cohorts and the usability task. “Unclassified” means no binary outcome was assigned, not that the experimental score is missing. No class balancing was applied.

## Exact counts

| policy | total | soluble | insoluble | percent_soluble | unclassified | parent_N | excluded | outcome |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PSI_original | 12216 | 8238 | 3978 | 67.44 | 0 | 12216 | 0 | solubility |
| PSI_NetSolP | 11226 | 7500 | 3726 | 66.81 | 0 | 11226 | 0 | solubility |
| NESG_NetSolP | 1323 | 842 | 481 | 63.64 | 0 | 1323 | 0 | solubility |
| UESolDS_development | 13880 | 7996 | 5884 | 57.61 | 0 | 13880 | 0 | solubility |
| Price_usability | 7259 | — | — | — | 0 | 7259 | 0 | usability_NOT_solubility |
| UESolDS_DNASU_proxy | 460 | 409 | 51 | 88.91 | 0 | 460 | 0 | solubility |
| NESG_raw_scores | 9703 | — | — | — | 9703 | 9703 | 0 | E_and_S_ordinal_scores_NO_binary_label |
| NESG_unique_scores | 9272 | — | — | — | 9272 | 9703 | 431 | E_and_S_ordinal_scores_NO_binary_label |
| NESG_consistent_scores | 9328 | — | — | — | 9328 | 9703 | 375 | E_and_S_ordinal_scores_NO_binary_label |
| NESG_raw_E45_S45_vs0 | 3943 | 2753 | 1190 | 69.82 | 0 | 9703 | 5760 | solubility |
| NESG_raw_E345_S45_vs0 | 4641 | 3149 | 1492 | 67.85 | 0 | 9703 | 5062 | solubility |
| NESG_raw_E45_S345_vs0 | 4672 | 3482 | 1190 | 74.53 | 0 | 9703 | 5031 | solubility |
| NESG_raw_E345_S345_vs0 | 5710 | 4218 | 1492 | 73.87 | 0 | 9703 | 3993 | solubility |
| NESG_raw_E45_S45_vs01 | 4227 | 2753 | 1474 | 65.13 | 0 | 9703 | 5476 | solubility |
| NESG_raw_expressed_ordinal | 7543 | — | — | — | 7543 | 9703 | 2160 | S_ordinal_conditioned_on_E_positive |
| NESG_unique_E45_S45_vs0 | 3762 | 2614 | 1148 | 69.48 | 0 | 9272 | 5510 | solubility |
| NESG_unique_E345_S45_vs0 | 4436 | 2992 | 1444 | 67.45 | 0 | 9272 | 4836 | solubility |
| NESG_unique_E45_S345_vs0 | 4458 | 3310 | 1148 | 74.25 | 0 | 9272 | 4814 | solubility |
| NESG_unique_E345_S345_vs0 | 5457 | 4013 | 1444 | 73.54 | 0 | 9272 | 3815 | solubility |
| NESG_unique_E45_S45_vs01 | 4036 | 2614 | 1422 | 64.77 | 0 | 9272 | 5236 | solubility |
| NESG_unique_expressed_ordinal | 7200 | — | — | — | 7200 | 9272 | 2072 | S_ordinal_conditioned_on_E_positive |
| NESG_consistent_E45_S45_vs0 | 3792 | 2638 | 1154 | 69.57 | 0 | 9328 | 5536 | solubility |
| NESG_consistent_E345_S45_vs0 | 4467 | 3017 | 1450 | 67.54 | 0 | 9328 | 4861 | solubility |
| NESG_consistent_E45_S345_vs0 | 4490 | 3336 | 1154 | 74.30 | 0 | 9328 | 4838 | solubility |
| NESG_consistent_E345_S345_vs0 | 5491 | 4041 | 1450 | 73.59 | 0 | 9328 | 3837 | solubility |
| NESG_consistent_E45_S45_vs01 | 4066 | 2638 | 1428 | 64.88 | 0 | 9328 | 5262 | solubility |
| NESG_consistent_expressed_ordinal | 7234 | — | — | — | 7234 | 9328 | 2094 | S_ordinal_conditioned_on_E_positive |

`parent_N` is the input to that row's filtering policy; published releases are initially audited as delivered. NetSolP PSI excludes 990 original SoDoPe records (738 positive, 252 negative). Original Price publication N = 9,644 is not silently substituted for the recovered archive N = 9,703; its exact row membership and binary counts remain unavailable.

## Additional eligibility views

These are separate manifests, not changes to the originals. Consistent-deduplicated PSI excludes all exact-sequence binary conflicts, then counts one lexicographically selected record per agreeing sequence. This cannot resolve differences in experimental constructs removed by NetSolP tag trimming.

| policy | total | soluble | insoluble | families | mixed_families | groups | families_groups_ge3 | families_groups_ge5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NESG_NetSolP_integer_rule | 1320 | 841 | 479 | 849 | 156 | 983 | 13 | 2 |
| NESG_NetSolP_no_PSI_component | 1107 | 698 | 409 | 717 | 113 | 839 | 6 | 2 |
| NESG_unique_E45_S45_vs0_no_PSI_component | 1176 | 741 | 435 | 770 | 120 | 899 | 7 | 2 |
| PSI_NetSolP_consistent_deduplicated | 11110 | 7402 | 3708 | 2075 | 959 | 3422 | 181 | 71 |
| PSI_original_consistent_deduplicated | 12147 | 8185 | 3962 | 2093 | 995 | 3529 | 187 | 77 |

## Complete E × S table, as distributed

| E | 0 | 1 | 2 | 3 | 4 | 5 |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 2092 | 10 | 6 | 18 | 3 | 22 |
| 1 | 74 | 71 | 20 | 25 | 7 | 20 |
| 1.50 | 1 | 0 | 0 | 0 | 0 | 0 |
| 2 | 143 | 80 | 170 | 70 | 34 | 96 |
| 2.50 | 0 | 1 | 0 | 0 | 0 | 0 |
| 3 | 302 | 78 | 160 | 340 | 116 | 280 |
| 3.50 | 2 | 0 | 0 | 1 | 0 | 2 |
| 4 | 301 | 79 | 161 | 237 | 390 | 334 |
| 4.50 | 2 | 0 | 0 | 0 | 0 | 0 |
| 5 | 889 | 205 | 340 | 492 | 376 | 1653 |

E = 0 rows are never treated as experimentally demonstrated insolubility. Nine fractional E records remain unrounded. The 59 E = 0 / S > 0 records need source clarification. Their count matching the publication/mirror size difference does not establish their historical identity.

## Prespecified sensitivity rules

- `E45_S45_vs0`: integer E in {4,5}; positive S in {4,5}; negative S = 0. Primary high-expression extreme-recovery policy.
- `E345_S45_vs0`: additionally admits E = 3; tests support gained at lower expression and potentially noisier recovery assessment.
- `E45_S345_vs0`: additionally calls S = 3 positive; changes the estimand to moderate-or-high recovery versus none.
- `E345_S345_vs0`: combines those relaxations; an exploratory sensitivity, not evidence that moderate scores are unambiguous extremes.
- `E45_S45_vs01`: includes S = 1 in the low-recovery class. These are **low-recovery**, not literally no-soluble-recovery negatives.
- `expressed_ordinal`: retains integer E > 0 with all S values; no binary relabeling. E and S are ordinal, not precise continuous concentrations.

`NESG_exclusion_reasons.tsv` partitions each raw/unique policy's exclusions into noninteger score, E = 0, expression below threshold, and intermediate S. `NESG_repeated_records.tsv` preserves every repeat and conflict. A missing/ambiguous measurement is never used as a negative. The 1,323-record NetSolP release is not identical to literal score-rule filtering: three fractional-E exceptions remain in the released file. The integer-rule view contains 1,320 records (841 positive / 479 negative).

## Direct pET15 outcome cross-check

| label | expression_status | solubility_status | purification_status | N |
| --- | --- | --- | --- | --- |
| 0 | Protein_Confirmed | Tested_Not_Soluble | Not_Tested | 1540 |
| 1 | Protein_Confirmed | Protein_Soluble | Not_Tested | 1124 |
| 1 | Protein_Confirmed | Protein_Soluble | Protein_Purified | 686 |
| 1 | Protein_Confirmed | Protein_Soluble | Tested_Not_Purified | 86 |

Class mappings are audited against the released raw status fields. pET21 has author classes and sequence files, but not the equivalent complete raw outcome-status export.
