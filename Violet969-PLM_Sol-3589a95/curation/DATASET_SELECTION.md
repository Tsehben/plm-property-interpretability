# Primary dataset decision — 2026-11-12

**Recommendation: policy C, one canonical heterogeneous cohort with separate evidence and analysis-eligibility flags. Retain all 13,880 proteins for biological-concept work, but do not yet designate any complete source as a validated primary binary solubility cohort.** Explicit DNASU outcomes are the first recovery priority; currently no reliable, comparable binary property core has been established.

This differs from treating all-source A as a scientifically validated solubility benchmark. A remains useful for exploratory published-label prediction and broad biological concepts. A good pooled score, source matching or restriction to a single centre cannot repair unsupported negative labels.

## Candidate policies: measured coverage

| policy | n | n0 | n1 | sources | families | clans | groups | mixed_families | families_both_splits_ge3 | family_source_pairs_both_splits_ge3 | families_lost | groups_lost |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_all | 13880 | 5884 | 7996 | 17 | 6043 | 557 | 13716 | 1660 | 69 | 4 | 0 | 0 |
| B_DNASU_candidate | 460 | 51 | 409 | 1 | 405 | 109 | 456 | 30 | 0 | 0 | 5638 | 13260 |
| C_canonical | 13880 | 5884 | 7996 | 17 | 6043 | 557 | 13716 | 1660 | 69 | 4 | 0 | 0 |
| C_verified_property_core | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 6043 | 13716 |
| D_JCSG | 2529 | 1101 | 1428 | 1 | 1547 | 285 | 2520 | 297 | 3 | 3 | 4496 | 11196 |
| D_MCSG | 2245 | 928 | 1317 | 1 | 1405 | 264 | 2233 | 273 | 1 | 1 | 4638 | 11483 |
| D_JCSG_MCSG | 4774 | 2029 | 2745 | 2 | 2503 | 370 | 4743 | 582 | 10 | 4 | 3540 | 8973 |
| D_training_supported_centres | 10058 | 5601 | 4457 | 8 | 4251 | 483 | 9960 | 1223 | 30 | 4 | 1792 | 3756 |

“Both splits ≥3” requires at least three proteins per label independently in training and validation. Pfam families, clans and homology groups are different units. A “mixed” family in the pooled cohort is weaker than a supported held-out within-family comparison. C_verified_property_core = 0 means that no comparable binary assay cohort was established; it does not deny the two recovered MCSG soluble statuses or the explicit DNASU repository statuses.

Original split/class counts for each policy:

| policy | split | n | n0 | n1 |
| --- | --- | --- | --- | --- |
| A_all | train | 9880 | 3884 | 5996 |
| A_all | validation | 4000 | 2000 | 2000 |
| B_DNASU_candidate | train | 334 | 35 | 299 |
| B_DNASU_candidate | validation | 126 | 16 | 110 |
| C_canonical | train | 9880 | 3884 | 5996 |
| C_canonical | validation | 4000 | 2000 | 2000 |
| C_verified_property_core | train | 0 | 0 | 0 |
| C_verified_property_core | validation | 0 | 0 | 0 |
| D_JCSG | train | 1639 | 832 | 807 |
| D_JCSG | validation | 890 | 269 | 621 |
| D_MCSG | train | 1718 | 561 | 1157 |
| D_MCSG | validation | 527 | 367 | 160 |
| D_JCSG_MCSG | train | 3357 | 1393 | 1964 |
| D_JCSG_MCSG | validation | 1417 | 636 | 781 |
| D_training_supported_centres | train | 7191 | 3788 | 3403 |
| D_training_supported_centres | validation | 2867 | 1813 | 1054 |

Exact source composition for every policy is in [policy_source_counts.tsv](source_audit_results/policy_source_counts.tsv); the [canonical manifest](source_audit_results/policy_manifest.tsv) retains original IDs, sequence hashes, source, split, label and homology group and adds eligibility/evidence columns only. Existing datasets and prior results are not overwritten.

## A — all-source heterogeneous property cohort

Retains all biological diversity and observed family–label correlations. Suitable for descriptive concept associations, source diagnostics and exploratory prediction of **published UESolDS labels**. Unsuitable as a verified solubility benchmark: PDB positives are structure-selected, TargetTrack negatives can reflect selection/cloning, and eSOL linkage is unresolved. Interventions changing this label readout would not establish a direct biochemical solubility mechanism.

## B — selected-source heterogeneous property cohort

DNASU is the only currently identified repository stratum with directly recoverable explicit positive and negative outcome fields. Its 460 proteins (409 positive, 51 negative), 405 families and 456 groups offer a candidate, not certified, property subset. Twelve audited statuses agree with released labels, but construct reconciliation is incomplete for all twelve and assay comparability is unresolved. It loses 5,638 families relative to A and has no ≥3/class family contrast in both splits. Thus it is not yet a substitute for the heterogeneous direct-versus-indirect study, even if its published labels are more plausible.

No alternative selected-source list has demonstrated better negative-outcome evidence. Excluding PDB/RefSeq solely because they are positive-only sacrifices diversity without repairing retained negatives. Selection should follow recoverable measurement definitions, never an AUROC ranking.

## C — canonical cohort plus evidence-qualified views (preferred)

Keep every original protein and family association. `concept_eligible=True` for all records; annotation-specific availability determines which concepts can be analyzed without imputing biology. `published_label_exploratory_only=True` preserves the released outcome as a descriptive target. `verified_property_eligible=False` until comparable construct-linked binary evidence is adjudicated. Separate flags record exact TargetTrack trial matches, explicit soluble stages and matching DNASU statuses. These flags distinguish recovery from certification.

For primary property labels, prioritize reconciled DNASU Protein_Soluble/Tested_Not_Soluble outcomes under a documented assay definition. Add JCSG, MCSG or other centres only when individual assayed constructs have explicit outcomes—including actual negatives—and compatible conditions. Do not use selected/cloned states, absent scores, centre names or successful structure deposition alone as primary binary truth. PDB and other extension proteins remain useful for family, structural and biochemical concept probing.

This is one dataset with transparent analysis eligibility, not a balanced replacement dataset. The primary label subset may ultimately have a narrower scope than the concept cohort; report that scope rather than claiming direct generalization to all proteins. If family × source support stays sparse after evidence recovery, controlled analysis B must remain exploratory or require new measured examples.

## D — single or limited centres

JCSG and MCSG provide complementary families and ten shared detected groups. Their union retains 2,503/6,043 families, losing 3,540, with four supported family/source cells. Neither has recovered explicit negative solubility measurements in the audited sample. The eight training-supported centres retain 10,058 proteins but still lose 1,792 families and do not increase the four controlled cells. Source comparability improves only at the proxy level; assay and negative-label comparability remain unproven. These are diagnostic views, not preferred primary cohorts.

## Suitability for the planned research

| Use | Current defensible scope | What remains necessary |
|---|---|---|
| Biological-concept probes | Heterogeneous canonical cohort with sequence/Pfam and documented predictor annotations | Respect annotation missingness, domain overlap and model-derived evidence |
| Correlational published-label probes | Exploratory all-source results with matched centre/family strata | Explicitly name the outcome as published UESolDS label, not validated solubility |
| Primary solubility readout | Not yet established | Construct-linked positive **and negative** outcomes with a consistent assay definition |
| SAE intervention feasibility | Can test technical reconstruction/ablation machinery separately | Valid property readout, token-level activations, matched controls, model-output effects; probe change alone is not causal biochemical evidence |
| Controlled within-family conclusions | Limited exploratory views | More independent families/groups with reliable two-class outcomes and fixed evaluation support |

## Final decision and smallest next step

**We cannot yet certify a reliable, biologically heterogeneous binary solubility cohort from the available UESolDS evidence.** We can preserve a useful heterogeneous concept cohort and now have a concrete route to stronger labels. Start with the fixed twelve DNASU records already audited: reconcile the original fusion sequences and obtain their expression/solubility score definitions and protocol versions. Then extend that exact rule to all 51 DNASU negatives and corresponding positives. Ask for the pre-cleaning UESolDS sequence-to-assay crosswalk or the original NESG/DNASU assay export; the current terminal differences are directly measurable, not hypothetical.

For structural-genomics negatives, request centre LIMS results for an actual soluble-expression test on the matching construct; a TargetTrack selection/cloning record cannot resolve the question. If those results cannot be recovered, a small matched experimental panel with measured soluble fractions is the next reliable option. No filtering, class balancing, label rewriting or new embeddings can substitute for that missing evidence.

Supporting evidence: [source validity](SOURCE_AUDIT.md), [diversity](SOURCE_DIVERSITY.md), [sequence-source and property probes](SOURCE_SHORTCUTS.md). The independent test remains unopened, all original records and previous artifacts are preserved, and cohort policy was not selected to maximize AUROC.
