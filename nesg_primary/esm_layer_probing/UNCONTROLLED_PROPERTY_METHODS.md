# NESG uncontrolled property probes

## Scope and execution

Open `ESM_Uncontrolled_Property_Probes_Colab.ipynb` in Colab and run all cells. A CPU runtime is sufficient. This is a new notebook; the executed extraction and solubility notebooks are preserved. It uses all 33 mean-pooled ESM-2 layers, not a selected subset. Controlled comparisons, SAE discovery, interventions, model fine-tuning and reserved-test evaluation are outside this run.

The notebook reads code and frozen development inputs from the existing Drive project. `HF_TOKEN_PLM` is read securely from Colab Secrets. Persistent activations/results use the existing private `PLM-Interp/nesg-esm2-layerwise` dataset. No repository is created or made public. A completed analysis is checked before tensor downloads; incomplete runs restore only pooled tensors, never residue maps or ESM weights.

## Populations, targets and provenance

The frozen development cohort is 2,361 training and 691 validation proteins. Protein IDs, sequence SHA-256, split assignments, solubility labels and homology groups remain unchanged. The reserved test was historically inspected for family-label composition, but this workflow does not access it or use it for tuning.

| Target | Definition | Observed train / validation |
|---|---|---:|
| Solubility | E ∈ {4,5}; S ∈ {4,5} versus S = 0 | 2,361 / 691 |
| Hydrophobicity | Existing mean Kyte–Doolittle hydropathy | 2,359 / 690 |
| Charge proxy | (K+R−D−E)/original length; not estimated net charge at a pH | 2,359 / 690 |
| Log length | Natural logarithm of original residue count | 2,361 / 691 |
| Aggregation propensity | Authors' AGGRESCAN server Na4vSS | 2,359 / 690 |
| APR burden | Authors' predicted APR residue fraction | 2,359 / 690 |
| Pfam membership | Presence of ≥1 GA-passing Pfam 37.0 assignment, per family | 2,361 / 691 per target |

Solubility is an experimental soluble-recovery phenotype under the NESG assay. Hydropathy, charge proxy and length are deterministic sequence descriptors. AGGRESCAN outputs are **predicted aggregation properties**, not experimentally measured aggregation. APR fraction is in [0,1]; genuine zero-hotspot results remain zero. The three X-containing sequences remain in the cohort; unsupported targets are missing and are not imputed. No additional disorder or secondary-structure labels are invented.

`property_inputs.py` verifies existing `results/development_records.tsv.gz`, `data/pfam_memberships.tsv`, and the development AGGRESCAN export. It reads only the `train_*` / `validation_*` job files listed by that export, checking sequence alignment and recorded output hashes. It never opens the combined aggregation table or reserved-test jobs. The frozen output is `data/uncontrolled/`; its `integrity.json` records exact source and generated-file hashes. Existing source annotation methods remain authoritative: no Pfam, clustering or server annotation is repeated.

All GA-passing family assignments are retained, including competing overlapping assignments. Multiple assignments do not duplicate a protein within a binary family target. Pfam negatives mean no detected assignment to that family and can include unannotated proteins; they do not establish experimentally verified absence of the family/domain. Domain membership does not establish identical full-protein architecture.

## Family task eligibility

Family targets require ≥5 distinct positive and ≥5 distinct negative **training** homology groups, and both membership classes in each of the five held-out training folds and its fitting complement. There are 48 families meeting the group threshold and **21** meeting the additional frozen-fold support requirement. The full catalogue records every exclusion reason. Validation labels/memberships do not determine eligibility.

This chooses family-prediction tasks, not a restricted dataset. Every development protein is included in every eligible membership task. Sixteen eligible targets have validation positives; five have none and retain undefined validation discrimination. Several of the sixteen have only one positive. These thresholds do not describe soluble-versus-insoluble support within a family. Membership is constant inside a single-family subset, so that later controlled subset cannot support a discrimination test for its defining family.

## Models and tuning

Each target is fitted independently at every layer. StandardScaler is learned separately inside every fitting fold. Binary tasks use L2 logistic regression (`lbfgs`, unweighted classes, max_iter=3000), with C ∈ {0.0001, 0.001, 0.01, 0.1, 1, 10}. Mean held-out-fold AUROC chooses C; exact ties prefer stronger regularization. A convergence warning stops the fit instead of silently accepting an unconverged model.

Continuous tasks use ridge regression, alpha ∈ {0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000}, chosen by lowest mean fold MSE; exact ties prefer larger alpha. Targets remain in their stated units. A float64 SVD of standardized fitting examples is reused across ridge candidates/targets sharing the same observed-example mask. The implementation is tested numerically against scikit-learn Ridge. Held-out examples never contribute to scaling or the SVD. Final coefficients and intercepts are refitted using all observed training examples.

For each task, the baseline uses the 20 standard residue fractions divided by original length, plus natural log length, with the same tuning and evaluation protocol. X is not substituted; the denominator stays the original length. Training-mean/prevalence constant references are also recorded. Hydropathy and charge are linear functions of composition; log length is itself a baseline input. Their near-perfect baseline performance is an expected implementation check.

The same frozen training folds select hyperparameters and provide reported training-CV scores: these are tuning estimates, not unbiased nested-CV estimates. Validation fits/tuning are prohibited. No family or class balancing occurs. No model or layer is selected automatically for subsequent experiments.

Compatible existing solubility bundles are identified using frozen input, model-activation, protocol and probe-code hashes. Selection among compatible bundles is deterministic, never based on performance. Existing predictions/parameters are reused and their ID/hash alignment is rechecked. If no compatible completed bundle exists, the runner announces a new fit under the original protocol. Other target fits are new.

## Metrics, figures and uncertainty

Binary targets report AUROC and average precision. Continuous targets report R², MAE and RMSE. Constant-target R² and single-class binary discrimination are undefined and retained as missing metrics. Negative R² is allowed. Validation AP must be interpreted relative to target prevalence.

The primary plots show training mean held-out-fold scores and validation scores across layers 1–33. Training bands show ±1 fold SD; validation bands use 2,000 homology-group bootstrap resamples (seed 20261118) conditional on fitted probes. The same resamples are used across layers and composition baselines, with target-specific missingness masks. `paired_vs_composition.tsv` reports ESM minus baseline differences and paired intervals. Negative MAE/RMSE differences favor ESM; positive AUROC/AP/R² differences favor ESM.

Family macro summaries use a fixed set of training-eligible families whose validation membership has both classes. All eligible tasks, including undefined validation metrics, are still saved and listed. Macro bootstrap draws require all members of that fixed set to have defined discrimination. Valid draw counts are reported; an interval is suppressed if fewer than 90% of draws remain valid, avoiding a heavily conditioned macro interval. No varying-family-set bootstrap average is substituted. Individual-family intervals remain exploratory with sparse positives; exact support tables and a named AUROC heatmap accompany them.

Intervals do not include fitting, hyperparameter, layer-selection or target-selection uncertainty. All-layer curves are descriptive, not 33 independent confirmatory tests. Strong linear decoding does not distinguish direct from indirect encoding; poor linear decoding does not demonstrate absence of information. The later controlled analysis needs an explicitly defined estimand and support checks: a restriction-induced AUROC change alone is not a causal confounding estimate.

## Files, checkpoints and reproduction

- `property_inputs.py`: build hash-verified development target inputs.
- `property_probes.py`: fit/reuse probes and save coefficients, scalers, predictions, CV and bootstrap metrics.
- `property_report.py`: complete-run tables and PNG/SVG plots.
- `run_property_probes.py`: CPU runner, cached-result checks and private HF checkpointing.
- `build_property_notebook.py`: regenerate this new notebook only (clears its outputs; do not run over an executed copy you want to preserve).
- `test_property_probes.py`: numerical, masking, training-only selection, missing-validation-class, privacy, rendering and cache-integrity tests.

The notebook reuses the pinned CPU package installer: huggingface-hub 0.36.0, SciPy 1.15.3, scikit-learn 1.7.2, joblib 1.5.2, threadpoolctl 3.6.0 and matplotlib 3.10.7. Working NumPy/Pandas versions are preserved and recorded, not silently downgraded. The existing frozen ESM model/revision/indexing remains in `protocol.json`. No ESM installation is needed for this CPU analysis.

HF checkpoints use `esm_stage1/uncontrolled_property_probing_v1/<input-code-software-fingerprint>/layer_XX/`. Layer 00 denotes composition. Each checkpoint contains `results.zip` plus a SHA-256 completion index, committed atomically; HTTP 429 respects server cooldowns. Reopening downloads completed layers instead of refitting them. Within an unchanged runtime, completed target fits inside the unfinished layer are reused by hash. The final `summary/results.zip` contains aggregate tables, the methods-aware report and figures. A completed summary skips all fitting and activation downloads.

Prediction tables contain stable ID, sequence hash, split, predicted value/decision score and (for binary tasks) probability. Observed per-protein targets, sequences, E/S scores, family memberships and homology/CV assignments remain in separate local manifests. They can be joined later without recomputing embeddings. Fitted standardized-coordinate coefficients can be applied as `((X - mean) / scale) @ coef + intercept`; use the sigmoid for logistic probabilities. Imported solubility parameters use the original saved array shapes.

Local checks, from this directory:

```sh
python -m unittest test_property_probes -v
python run_property_probes.py --baseline-only --output /tmp/nesg-uncontrolled-validation --threads 2
```

These require the already prepared input snapshot. Only source-annotation rebuilding needs `python property_inputs.py`. Full all-layer execution requires authenticated access to the private pooled cache and is performed by the Colab notebook. No full new ESM-property run is claimed from local baseline checks.
