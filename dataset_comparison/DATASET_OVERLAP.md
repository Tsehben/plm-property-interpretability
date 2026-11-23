# Dataset overlap, conflicting outcomes and independence

These are releases of related experiments, not six independent cohorts. Exact protein hashes, original IDs and direct homology hits are recorded separately. No sequence was modified to manufacture an exact match.

## Direct overlap

| A | B | exact_sequences | exact_conflicting_sequences | A_direct_30 | B_direct_30 | A_direct_90 | B_direct_90 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| PSI_NetSolP | NESG_raw_scores | 2289 | 0 | 4585 | 6115 | 3942 | 4105 |
| PSI_NetSolP | NESG_unique_E45_S45_vs0 | 1138 | 4 | 3527 | 2514 | 2189 | 1921 |
| PSI_NetSolP | NESG_NetSolP | 0 | 0 | 504 | 191 | 0 | 0 |
| PSI_NetSolP | UESolDS_development | 99 | 3 | 2357 | 852 | 960 | 548 |
| NESG_unique_E45_S45_vs0 | UESolDS_development | 24 | 8 | 736 | 411 | 189 | 178 |

`A_direct_*` / `B_direct_*` count records on that side with a direct hit to the other at the stated identity and ≥80% bilateral coverage. Full pairwise results, including ≥95%, are in `dataset_overlap.tsv`. Exact conflict counts concern assigned binary outcomes; raw NESG E/S rows have no binary labels and thus their zero binary-conflict count is **not** evidence of agreement. Score conflicts are retained separately.

Original tagged PSI and trimmed NetSolP PSI have zero exact full-sequence overlaps, despite 11,226 accession-matched records with identical labels and fully explained terminal offsets. All 460 UESolDS DNASU-proxy IDs map to original PSI records with agreeing labels, and every cohort sequence is contained in the original construct; none are full-sequence identical. This demonstrates why exact hashes alone understate shared biological/experimental provenance.

PSI NetSolP and strict NESG share 1,138 exact sequences; four have different binary labels. These may reflect different experimental records, constructs, thresholds or updates. Do not overwrite either outcome. At ≥90% identity, additional near-duplicate pairs disagree; their counts are in the overlap table. Common-source history means the remaining proteins are not necessarily independent assays even when sequence leakage is removed.

## Duplicates

| dataset | records | unique_sequences | repeated_sequences | repeated_records | conflicting_binary_sequences | conflicting_binary_records |
| --- | --- | --- | --- | --- | --- | --- |
| NESG_NetSolP | 1323 | 1323 | 0 | 0 | 0 | 0 |
| NESG_raw | 9703 | 9478 | 206 | 431 | 0 | 0 |
| PSI_NetSolP | 11226 | 11122 | 98 | 202 | 12 | 26 |
| PSI_original | 12216 | 12157 | 53 | 112 | 10 | 22 |
| Price_usability | 7259 | 7259 | 0 | 0 | 0 | 0 |
| UESolDS_development | 13880 | 13880 | 0 | 0 | 0 | 0 |

Binary-conflict fields are not defined as ordinal E/S disagreement. For raw NESG, 150 exact-sequence groups / 316 records have differing E or S scores. Every repeated record is retained in `NESG_repeated_records.tsv`; PSI and cross-dataset binary conflicts are retained in `duplicate_records.tsv` and `exact_conflicts.tsv`. Group support prevents repeat records from masquerading as independent examples.

## Published versus new partitions

The released PSI partitions keep all exact duplicates within one fold. Two cross-fold direct MMseqs pairs pass our local ≥30%/80%-coverage rule (maximum identity 30.7%). This is a different alignment/normalization criterion from the published ≤25% global ggsearch rule, so these hits do not by themselves refute that claim. We verified the delivered folds with our specified criterion; we did not rerun proprietary USEARCH or reproduce the full global-identity graph partitioner.

The 1,323-record NetSolP NESG release has zero exact or ≥90% direct matches to PSI, but 191 records have direct ≥30%/80% hits to curated PSI. Therefore its published external designation is not sufficient under our operational criterion. A conservative removal of components intersecting **original PSI** leaves 1,107 released NESG records (698 soluble / 409 insoluble). This view still preserves the published fractional-score labels; use the integer-score restriction as well for a strict experimental readout.

| policy | total | soluble | insoluble | families | mixed_families | groups | families_groups_ge3 | families_groups_ge5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NESG_NetSolP_no_PSI_component | 1107 | 698 | 409 | 717 | 113 | 839 | 6 | 2 |
| NESG_unique_E45_S45_vs0_no_PSI_component | 1176 | 741 | 435 | 770 | 120 | 899 | 7 | 2 |

Component exclusion is more conservative than direct hit exclusion: intermediates from any audit dataset can bridge a component. All eligibility decisions are separate flags/manifests; originals remain intact. Because we inspected these releases for feasibility, they are not described as untouched prospective holdouts. Before future modeling, freeze the selected component-level split and isolate its evaluation labels.

## UESolDS boundary

Only the frozen 13,880 operational proteins were included in the sequence/Pfam/homology union. Exact-hash overlap with the larger original train/validation manifest is separately available in `full_training_validation_hash_overlap.tsv`; no additional UESolDS sequences were extracted. **The locked UESolDS independent test was never accessed.** Existing records, annotation files, embeddings and original cluster assignments were not edited.
