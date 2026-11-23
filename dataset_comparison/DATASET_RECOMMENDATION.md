# Recommendation for direct-versus-indirect encoding experiments

**Use the strict NESG score cohort as the primary experimental solubility dataset: 3,762 proteins, 2,614 soluble and 1,148 insoluble.** Preserve its full 9,272 unique-sequence score archive alongside it; use the broader expressed-protein ordinal view for a separately declared sensitivity analysis. The strongest labels are the high-expression extremes, not all S = 0 entries in the raw archive.

This is a qualified recommendation for the pinned public NESG/SoluProt archive, not a claim that the missing original 9,644-record release was recovered or that every row's assay protocol is independently verified. It is based on interpretable measured outcomes and sufficient heterogeneous family coverage, not on model performance. No ESM/SAE probes were trained.

## Candidate comparison

| policy | total | families | mixed_families | groups | families_groups_ge3 | families_groups_ge5 |
| --- | --- | --- | --- | --- | --- | --- |
| PSI_original | 12216 | 2093 | 996 | 3530 | 187 | 77 |
| PSI_NetSolP | 11226 | 2075 | 962 | 3423 | 181 | 71 |
| NESG_NetSolP | 1323 | 852 | 157 | 985 | 13 | 2 |
| UESolDS_DNASU_proxy | 460 | 405 | 30 | 424 | 0 | 0 |
| NESG_raw_scores | 9703 | 2101 | — | 3734 | — | — |
| NESG_unique_scores | 9272 | 2086 | — | 3697 | — | — |
| NESG_unique_E45_S45_vs0 | 3762 | 1383 | 390 | 1820 | 55 | 28 |

- **Primary A:** Retain all 3,762 strict-score proteins, including single-class families, no-hit proteins and the natural 69.48% soluble prevalence. This retains 1,383 Pfam families, 243 clans and 1,820 audit homology components. Do not balance family distributions or restrict A to mixed families.
- **Controlled B, from the same A:** 390 families have both labels; 55 have ≥3 groups per class and 28 have ≥5. These cover 883 and 698 proteins respectively (union counts). Only 35/13 families meet the stricter class-exclusive-group thresholds. Under the prespecified split, 21/15 families have ≥3/≥5 development groups/class plus both labels in evaluation. These are candidate contrasts, not 55/28 independent biological mechanisms. Competing Pfam assignments reduce conservative support to 16/3 families; select a small biologically interpretable, architecture-checked panel using development data before any probes.
- **Broader biological comparison:** NetSolP PSI is the stronger option for family diversity and broad B coverage: 2,075 families, 962 mixed families, 181/71 with ≥3/≥5 groups per class. Retain its original DNASU records for auditing. If used for primary binary modeling, flag the 12 conflicting exact-sequence groups and avoid counting repeated sequences as independent. The consistent-deduplicated view has 11,110 proteins (7,402 soluble / 3,708 insoluble), 959 mixed families, and preserves 181/71 group-supported families. It is a credible secondary dataset, but its binary statuses give less quantitative control over expression/recovery than NESG scores.
- **NetSolP NESG alone:** Excellent for reproducing a published external benchmark, insufficient as the main controlled-family dataset: 13/2 group-supported families, falling to 2/0 under conservative overlap handling. The three fractional-E exceptions also need explicit treatment. Its strict integer-rule view is 1,320, not 1,323.
- **DNASU proxy in UESolDS:** 460 records / 409 positive / 51 negative; no family meets ≥3 proteins of each class. It cannot substitute for full PSI or provide the desired B analysis.

## Ordinal NESG analysis

The 7,200 unique-sequence records with integer E > 0 provide an ordinal S analysis spanning intermediate outcomes. Keep E explicit, report expression strata, and do not call S precise continuous solubility. The alternative E ≥3/extreme-S policy yields 4,436 proteins (2,992 / 1,444), 1,522 families and 66/31 group-supported families. Its larger support comes with a changed expression-conditioned population, so it remains a prespecified sensitivity rather than an automatic replacement for the primary rule.

## Training/development versus external evaluation

Do not simply train on PSI and call the full NESG archive independent: the strict cohorts share 1,138 exact sequences and extensive homology. Conversely, fitting on full NESG prevents treating all PSI as an independent external evaluation. A conservative PSI-development / NESG-external design can retain 1,176 strict NESG proteins (741 / 435) after removing components shared with original PSI, but has only seven families with ≥3 groups/class. It supports an aggregate external check better than broad within-family validation. See the explicit eligibility manifest and overlap report before choosing direction.

Do not pool PSI and NESG binary labels into one undifferentiated training target. Their related history, duplicate experiments, vector differences and binary-versus-score definitions would recreate the provenance problem we are trying to address. They can be compared as separate outcome strata in a shared annotation inventory, with shared sequence components assigned consistently.

## Remaining information and smallest useful next steps

1. Obtain the original 9,644-ID E/S export (or a dated mapping to the public 9,703 records), including the meaning of fractional E and E = 0 / S > 0 entries. The present strict integer/high-E policy avoids those anomalies but cannot establish historical membership.
2. Recover replicate/construct/condition details for the four exact PSI–strict-NESG binary disagreements and the 12 conflicting PSI sequence groups. This small targeted audit is more informative than collecting more embeddings.
3. For a development-selected B panel, inspect overlapping Pfam calls and complete architectures, and verify assay/construct comparability. Freeze a homology-aware evaluation protocol before training. Our saved five-fold feasibility split is descriptive and must not be portrayed as the original untouched external test.

The evidence supports a useful experimental-solubility representation study with a conditional target: **soluble recovery of a recombinant construct in the documented E. coli production setting**. It does not turn sequence associations into causal computations. SAE interventions on the model would test effects inside the model; experimental changes in protein solubility would still require wet-lab validation.
