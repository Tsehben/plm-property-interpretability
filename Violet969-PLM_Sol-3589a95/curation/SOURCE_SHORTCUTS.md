# Sequence-accessible source information and solubility sensitivity — 2026-11-12

**Source-associated information is accessible from sequence and frozen ESM activations.** It does not follow that the solubility readout uses a particular source cue or that a centre effect is technical rather than biological. The original-record audit independently finds weak evidence for many negative labels.

## Fixed evaluation protocol

The models use exactly 9,880 training and 4,000 validation records, with original roles, labels and detected homology groups preserved. No detected group, exact sequence or reused identifier spans these splits. This protects against detected leakage at the existing MMseqs threshold, not all remote homology or domain overlap. Validation has already been inspected repeatedly; all results are exploratory development evidence. The independent test remains locked.

Source prediction uses all 17 frozen header-derived centre/repository classes. Fifteen occur in validation (BSGC and CESG do not). Composition is 20 AA fractions plus log length on the first 1,022 residues; frozen ESM-2 layer-8 activations are cached 1,280-dimensional means on the same visible residues. Pfam uses full-sequence multi-hot memberships, a training-only family vocabulary and an explicit no-hit indicator; it has a modest full-length information advantage over ESM on long proteins. Validation-only families are ignored as unseen features, not counted as absent biological annotation.

L2 logistic regression, C=1, lbfgs, tolerance 1e-4, maximum 5,000 iterations, no class weights or hyperparameter search. Dense sequence inputs are standardized using only the relevant fitting training subset. Encoders see training only. All fits converged (source: 74–425 iterations; solubility transfer: 7–209). No ESM weights were loaded or updated. The existing all-source solubility predictions are reused, including their saved precision; new fits are limited to training-source sensitivities.

Approximate 95% intervals use 200 resamples of validation homology groups, keeping group members together. Paired differences use identical resamples. These are conditional on fitted models and do not include training-fit or label uncertainty; there is no multiplicity adjustment. AUROC is undefined for single-class strata. AP is interpreted with observed positive prevalence, and is trivially 1 for all-positive strata. For tiny minorities—especially SSGCID with one negative—bootstrap intervals cannot support reliable population inference; undefined bootstrap draws are counted and excluded.

## Source prediction

| model | accuracy | balanced_accuracy | macro_f1 | macro_ovr_auroc | macro_ap | logloss |
| --- | --- | --- | --- | --- | --- | --- |
| training_prior | 0.187 [0.176, 0.199] | 0.067 [0.067, 0.067] | 0.021 [0.020, 0.022] | 0.500 [0.500, 0.500] | 0.067 [0.067, 0.067] | 2.343 [2.309, 2.376] |
| composition | 0.292 [0.279, 0.307] | 0.142 [0.136, 0.150] | 0.130 [0.123, 0.138] | 0.712 [0.698, 0.728] | 0.149 [0.145, 0.159] | 2.117 [2.079, 2.153] |
| Pfam | 0.296 [0.284, 0.309] | 0.157 [0.149, 0.168] | 0.149 [0.137, 0.161] | 0.658 [0.646, 0.674] | 0.167 [0.159, 0.182] | 2.146 [2.108, 2.191] |
| ESM_layer8 | 0.370 [0.357, 0.384] | 0.227 [0.215, 0.239] | 0.222 [0.212, 0.233] | 0.745 [0.732, 0.760] | 0.229 [0.220, 0.246] | 4.269 [4.109, 4.461] |
| composition_plus_ESM | 0.369 [0.355, 0.383] | 0.231 [0.219, 0.243] | 0.227 [0.215, 0.237] | 0.747 [0.733, 0.761] | 0.229 [0.220, 0.247] | 4.300 [4.135, 4.494] |

The training-prior predictor always selects the most frequent training source (PDB_serial): 18.7% validation accuracy. Balanced accuracy is mean recall over the 15 represented validation classes; macro F1 uses those classes too. Macro one-vs-rest AUROC/AP omit absent validation classes, and log loss retains all 17 probability columns. The explicit prior model accounts for unequal class sizes; chance accuracy is not assumed to be 1/17.

ESM improves classification accuracy over composition by **7.80 percentage points** (paired 95% CI 6.25–9.33), and balanced accuracy by about 8.43 points. However, ESM has much worse log loss (4.269 versus 2.117), indicating poor confidence calibration under this fixed high-dimensional probe and source shift. The model has not been calibrated or tuned after seeing validation. AUROC/accuracy cannot conceal that failure.

Adding ESM to composition improves accuracy over composition by **7.72 points** (6.02–9.46). This demonstrates incremental linear predictive information under the specified models; it is not a causal adjustment or proof of independence from all nonlinear composition effects. Source-target identity itself mixes laboratory, repository and naming conventions.

[Per-centre precision, recall, support, one-vs-rest AUROC/AP and intervals](source_audit_results/source_prediction_per_centre.tsv), [confusion matrices](source_audit_results/source_confusion.tsv), [paired increments](source_audit_results/source_paired_differences.tsv) and all probability tables are retained.

## Source distinguishability within families

| model | families | accuracy | balanced_accuracy | macro_ovr_auroc | logloss |
| --- | --- | --- | --- | --- | --- |
| training_prior | 27 | 0.3517 | 0.4130 | 0.5000 | 1.1665 |
| composition | 27 | 0.4719 | 0.4836 | 0.6814 | 0.9959 |
| Pfam | 27 | 0.4498 | 0.4105 | 0.5333 | 1.0527 |
| ESM_layer8 | 27 | 0.5814 | 0.5744 | 0.6620 | 1.9275 |
| composition_plus_ESM | 27 | 0.5780 | 0.5771 | 0.6847 | 1.9202 |
| family_training_prior | 27 | 0.4664 | 0.4130 | 0.5000 | 0.9642 |

There are **27 qualifying families** with at least two sources having ≥3 records in both training and validation. Restrict candidate sources to those supported in that family and renormalize every model equally. A family-specific training-frequency prior provides a stronger baseline than a global prior. These are equal-family descriptive means; memberships overlap, and they are not independent replicates or a validation-selected optimal subset. Exact per-family counts and metrics are in source_within_family.tsv and source_within_family_support.tsv.

Per-family intervals are saved in [source_within_family_uncertainty.tsv](source_audit_results/source_within_family_uncertainty.tsv). They resample validation homology groups within each fixed stratum; draws missing a candidate source are excluded and counted. There is no independence claim across overlapping families and no macro-mean confidence interval.

Source discrimination survives family conditioning: composition and ESM both show conditional ranking signal. ESM does not dominate every metric: mean conditional one-vs-rest AUROC is about 0.662 for ESM and 0.681 for composition. Sharing one Pfam does not hold all domains, architecture, organism or chemistry fixed. These results do not identify an experimental artefact.

## Solubility: pooled and within-centre, same fitted predictions

| evaluation | n | n0 | n1 | prevalence | composition_AUROC_CI | composition_AP | ESM_layer8_AUROC_CI | ESM_layer8_AP |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ALL | 4000 | 2000 | 2000 | 0.5000 | 0.656 [0.639, 0.671] | 0.64408698930867869 | 0.677 [0.659, 0.691] | 0.67926182928558776 |
| BSGI | 35 | 35 | 0 | 0.0000 | undefined |  | undefined |  |
| DNASU_clone | 126 | 16 | 110 | 0.8730 | 0.627 [0.477, 0.760] | 0.91960653141488224 | 0.619 [0.506, 0.741] | 0.92848464105997364 |
| EFI | 92 | 71 | 21 | 0.2283 | 0.499 [0.348, 0.693] | 0.29131306814578345 | 0.482 [0.337, 0.636] | 0.2619593494115926 |
| JCSG | 890 | 269 | 621 | 0.6978 | 0.534 [0.491, 0.571] | 0.75332838516673006 | 0.619 [0.581, 0.656] | 0.80602453428900511 |
| MCSG | 527 | 367 | 160 | 0.3036 | 0.497 [0.452, 0.553] | 0.29304702609618283 | 0.561 [0.506, 0.612] | 0.35062902480609603 |
| MPP | 57 | 57 | 0 | 0.0000 | undefined |  | undefined |  |
| NESG | 152 | 152 | 0 | 0.0000 | undefined |  | undefined |  |
| NYSGRC | 384 | 384 | 0 | 0.0000 | undefined |  | undefined |  |
| PDB_serial | 748 | 0 | 748 | 1.0000 | undefined | 1 | undefined | 1 |
| RefSeq | 32 | 0 | 32 | 1.0000 | undefined | 1 | undefined | 1 |
| SECSG | 300 | 208 | 92 | 0.3067 | 0.365 [0.305, 0.433] | 0.24034245982337557 | 0.461 [0.393, 0.528] | 0.27096228378161358 |
| SGPP | 455 | 361 | 94 | 0.2066 | 0.595 [0.530, 0.656] | 0.26019760483202781 | 0.491 [0.428, 0.552] | 0.21260237936083135 |
| SSGCID | 67 | 1 | 66 | 0.9851 | 0.121 [0.061, 0.212] | 0.97029720979684886 | 0.576 [0.462, 0.688] | 0.99188167623575973 |
| TB_Rv | 100 | 53 | 47 | 0.4700 | 0.523 [0.407, 0.635] | 0.51278867149885954 | 0.487 [0.365, 0.614] | 0.46165086690104401 |
| other | 35 | 26 | 9 | 0.2571 | 0.829 [0.637, 0.971] | 0.62401394901394913 | 0.573 [0.315, 0.828] | 0.47993567797489367 |

Pooled ESM AUROC is 0.6769 versus 0.6565 for matched composition; the paired difference is 0.0204 [0.0030, 0.0360]. The metadata-only source model is 0.7791 [0.7664, 0.7915]. These compare the same 4,000 records. They do **not** share an estimand with a five-centre mean, and subtracting a pooled score from a within-centre average is not a causal source-reliance estimate.

Per-centre paired ESM minus composition AUROC:

| evaluation | n | difference | low | high | valid_bootstraps |
| --- | --- | --- | --- | --- | --- |
| ALL | 4000 | 0.0204 | 0.0030 | 0.0360 | 200 |
| DNASU_clone | 126 | -0.0085 | -0.1716 | 0.1595 | 200 |
| EFI | 92 | -0.0168 | -0.1865 | 0.1786 | 200 |
| JCSG | 890 | 0.0856 | 0.0417 | 0.1264 | 200 |
| MCSG | 527 | 0.0639 | 0.0116 | 0.1104 | 200 |
| SECSG | 300 | 0.0966 | 0.0172 | 0.1679 | 200 |
| SGPP | 455 | -0.1037 | -0.1766 | -0.0370 | 200 |
| SSGCID | 67 | 0.4545 | 0.3177 | 0.5765 | 119 |
| TB_Rv | 100 | -0.0357 | -0.1813 | 0.0987 | 200 |
| other | 35 | -0.2564 | -0.5043 | 0.0146 | 200 |

ESM has higher ranking performance than composition within JCSG and MCSG, but is weaker in SGPP. SECSG improves relative to an inverted composition signal while remaining near/below chance. Centre-wise variation, tiny strata and possible label-generation differences prevent interpreting these as universal biochemical laws. A constant metadata source model scores 0.5 within a centre by construction; that does not exclude finer source-related sequence signals.

## Cross-centre training

| protocol | model | evaluation | n0 | n1 | AUROC_95CI | ap | prevalence |
| --- | --- | --- | --- | --- | --- | --- | --- |
| only_JCSG | ESM_layer8 | JCSG | 269 | 621 | 0.645 [0.608, 0.679] | 0.81293936889531926 | 0.6978 |
| only_JCSG | ESM_layer8 | MCSG | 367 | 160 | 0.569 [0.514, 0.625] | 0.36896849748318539 | 0.3036 |
| only_JCSG | composition | JCSG | 269 | 621 | 0.640 [0.605, 0.673] | 0.81671394171149014 | 0.6978 |
| only_JCSG | composition | MCSG | 367 | 160 | 0.531 [0.483, 0.587] | 0.32254715232599895 | 0.3036 |
| only_MCSG | ESM_layer8 | JCSG | 269 | 621 | 0.552 [0.512, 0.587] | 0.76006531164532198 | 0.6978 |
| only_MCSG | ESM_layer8 | MCSG | 367 | 160 | 0.481 [0.433, 0.530] | 0.2892090893286845 | 0.3036 |
| only_MCSG | composition | JCSG | 269 | 621 | 0.511 [0.474, 0.548] | 0.72711757343337746 | 0.6978 |
| only_MCSG | composition | MCSG | 367 | 160 | 0.434 [0.385, 0.480] | 0.26345364207303423 | 0.3036 |

Each row trains only on the named centre’s original training proteins and evaluates the named, unchanged validation stratum. JCSG→MCSG and MCSG→JCSG are separate transfer tasks, not paired population comparisons. Their results do not upgrade either centre’s label validity. The complete matrix includes all validation sources, with intervals in cross_centre_uncertainty.tsv.

## Leave-one-source-out training: fixed 4,000-protein validation population

| protocol | model | n | AUROC_95CI | ap |
| --- | --- | --- | --- | --- |
| without_BSGC | composition | 4000 | 0.656 [0.639, 0.671] | 0.6440198212499757 |
| without_BSGC | ESM_layer8 | 4000 | 0.677 [0.660, 0.692] | 0.68004288749672526 |
| without_BSGI | composition | 4000 | 0.656 [0.639, 0.671] | 0.64392358419963713 |
| without_BSGI | ESM_layer8 | 4000 | 0.676 [0.659, 0.690] | 0.67841093863129576 |
| without_CESG | composition | 4000 | 0.656 [0.639, 0.671] | 0.64404976042486317 |
| without_CESG | ESM_layer8 | 4000 | 0.676 [0.659, 0.691] | 0.67964091492866241 |
| without_DNASU_clone | composition | 4000 | 0.657 [0.640, 0.672] | 0.64555309497162772 |
| without_DNASU_clone | ESM_layer8 | 4000 | 0.676 [0.658, 0.691] | 0.68040535150980919 |
| without_EFI | composition | 4000 | 0.657 [0.640, 0.672] | 0.64322826714862313 |
| without_EFI | ESM_layer8 | 4000 | 0.676 [0.659, 0.690] | 0.6770072382740363 |
| without_JCSG | composition | 4000 | 0.647 [0.630, 0.662] | 0.62857247932130045 |
| without_JCSG | ESM_layer8 | 4000 | 0.660 [0.642, 0.674] | 0.64912015227519659 |
| without_MCSG | composition | 4000 | 0.655 [0.637, 0.670] | 0.64528145536025516 |
| without_MCSG | ESM_layer8 | 4000 | 0.673 [0.657, 0.687] | 0.67475835097813852 |
| without_MPP | composition | 4000 | 0.656 [0.639, 0.671] | 0.64401372365087606 |
| without_MPP | ESM_layer8 | 4000 | 0.678 [0.660, 0.692] | 0.6798122650778492 |
| without_NESG | composition | 4000 | 0.657 [0.639, 0.672] | 0.64447059140399687 |
| without_NESG | ESM_layer8 | 4000 | 0.676 [0.660, 0.690] | 0.67736758851671997 |
| without_NYSGRC | composition | 4000 | 0.653 [0.635, 0.668] | 0.64609290834448707 |
| without_NYSGRC | ESM_layer8 | 4000 | 0.664 [0.646, 0.679] | 0.66757211156959606 |
| without_PDB_serial | composition | 4000 | 0.648 [0.631, 0.664] | 0.63031498402114794 |
| without_PDB_serial | ESM_layer8 | 4000 | 0.634 [0.614, 0.648] | 0.63228242999515016 |
| without_RefSeq | composition | 4000 | 0.657 [0.640, 0.672] | 0.64479875378712093 |
| without_RefSeq | ESM_layer8 | 4000 | 0.678 [0.661, 0.692] | 0.68135450042878087 |
| without_SECSG | composition | 4000 | 0.659 [0.642, 0.675] | 0.64774599140620603 |
| without_SECSG | ESM_layer8 | 4000 | 0.686 [0.669, 0.701] | 0.69308503083846118 |
| without_SGPP | composition | 4000 | 0.657 [0.640, 0.674] | 0.64327970001809598 |
| without_SGPP | ESM_layer8 | 4000 | 0.667 [0.651, 0.680] | 0.67278225030357164 |
| without_SSGCID | composition | 4000 | 0.656 [0.639, 0.671] | 0.64438011288968411 |
| without_SSGCID | ESM_layer8 | 4000 | 0.678 [0.661, 0.692] | 0.68036482609338478 |
| without_TB_Rv | composition | 4000 | 0.657 [0.639, 0.671] | 0.64414298699225947 |
| without_TB_Rv | ESM_layer8 | 4000 | 0.676 [0.659, 0.690] | 0.67921756375754105 |
| without_other | composition | 4000 | 0.656 [0.639, 0.671] | 0.64387584939884634 |
| without_other | ESM_layer8 | 4000 | 0.677 [0.660, 0.692] | 0.67923381003437489 |

All training-source omissions above evaluate the **same** validation proteins; no removed source is also silently removed from evaluation. Separate omitted-centre-only results and paired AUROC/AP changes are saved. Removing PDB training lowers pooled ESM AUROC to about 0.634; this changes training size, class composition and biology simultaneously. It is a sensitivity result, not a causal estimate of a PDB shortcut or a reason to discard PDB proteins.

## Family-overlapping versus family × source solubility evaluation

| kind | model | strata | mean_auroc | mean_ap |
| --- | --- | --- | --- | --- |
| family | ESM_layer8 | 69 | 0.6583 | 0.6430 |
| family | composition | 69 | 0.6818 | 0.6220 |
| family | source | 69 | 0.7761 | 0.7950 |
| family_source | ESM_layer8 | 4 | 0.6317 | 0.7660 |
| family_source | composition | 4 | 0.5119 | 0.6522 |
| family_source | source | 4 | 0.5000 | 0.5179 |

These use global all-source fitted predictions, not separate family-tuned fits. Training support is ≥3/class, and the reported strata also have ≥3/class in validation: 69 Pfam families versus only four family/source cells. The source model still ranks labels within families (mean AUROC 0.776), showing that family matching alone leaves source differences. Within family/source cells the source model is constant, while sequence probes can vary. Only three related methyltransferase accessions supply the four cells; this tiny, overlapping set cannot establish a general within-family mechanism. [Individual metrics and group-bootstrap intervals](source_audit_results/solubility_family_strata.tsv) accompany the means. No family associations were erased or rebalanced.

## Interpretation and feasibility

High pooled AUROC does not validate labels. Near-chance within-centre performance does not establish that ESM lacks solubility information. Predictable source information does not show that the solubility probe uses it. Predictive readouts or overlapping probe weights do not establish causal computation; ablation must distinguish a change in the probe from a biologically meaningful output change.

The cached cohort suffices for this audit. Expanding embeddings is not the next priority: previous analysis identified 59,368 extra clean training proteins, but their uniform Pfam annotations are unknown and unchanged validation limits controlled contrasts. First recover verified construct-specific outcomes. Mean-pooled caches suffice for these logistic probes, but future SAE interventions require appropriate token-level forward passes and reconstruction controls; no SAE analysis was performed here.

Reproduce with `source_validity_audit.py --stage probes`, `source_audit_uncertainty.py`, then `render_source_audit.py`. Protocol, fit diagnostics, predictions and input hashes are saved in source_audit_results/. The procedure does not rerun annotations, clustering or ESM extraction.
