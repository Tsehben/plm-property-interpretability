# Within-property evaluation of unchanged NESG probes

## Purpose and order

Run `ESM_Within_Property_Evaluation_Colab.ipynb` after the uncontrolled-property notebook finishes and its plots have been inspected. Keep all 33 layers for the primary descriptive analysis. The expensive global probes already exist: this notebook evaluates their saved predictions on subsets. Selecting layers 4/16/32 may simplify a later illustration but is unnecessary for computation and must not replace the complete curves or be chosen to maximize an observed drop.

There is no model fitting, hyperparameter selection, activation extraction, SAE analysis, intervention, cohort modification or reserved-test access. The 2,361 training proteins define restrictions/support; metrics use the unchanged 691 validation proteins, with target-specific missingness. The three X-containing proteins are not removed from the canonical cohort. No training-fit AUROC is presented as validation performance.

## Questions evaluated

1. Experimental soluble-recovery prediction within ranges of hydropathy, charge proxy, log length, predicted AGGRESCAN Na4vSS and predicted APR fraction.
2. Soluble-recovery prediction within each training-observed Pfam family, using all GA-passing memberships and a separate noncompeting-assignment sensitivity.
3. The five continuous targets within the same numeric/family subsets; the existing 21 family-membership targets within numeric subsets.

Every subset uses the same globally trained probe at that layer. No family-specific model is trained. Family membership is constant inside its defining family; family-membership probes are therefore not evaluated inside family subsets. Global family probes are evaluated within numeric ranges where membership may vary.

## Restrictions and support

For numeric properties, cut points are the 20th, 40th, 60th and 80th percentiles of observed training values, using NumPy's linear quantile method. Duplicate cut points and endpoint cuts are collapsed. Intervals are lower-inclusive, upper-exclusive, with unbounded outer intervals. Missing descriptors belong to no numeric range. Training quantiles are not selected using labels, validation values or model performance. These broad ranges are approximate restrictions, not exact matching or a claim that all confounding is removed. Only one property is restricted at a time.

Family identities come from training-observed Pfam 37.0 GA assignments. There is no ≥3/≥5 training-solubility-support threshold because no new probe is fitted. All 949 training-observed families are in the support output. A family/assignment variant needs at least two validation proteins for point calculations; binary discrimination additionally needs both labels. Unsupported families are retained in support tables, rather than assigned chance scores. Validation-only families are outside this training-defined family catalogue.

Membership is protein-level: repeated hits do not duplicate a protein. A multidomain protein may participate in separate evaluations for multiple families. Those evaluations overlap and must not be pooled into a duplicated protein population or interpreted as independent replicates. In the sensitivity variant, memberships marked `competing_overlap` are excluded for that family only. Primary all-GA results remain intact. Shared domain membership does not imply identical architecture, organism, experimental construct or assay.

The local support audit found 46 all-GA and 20 noncompeting family strata with both validation solubility labels. These are **estimable point comparisons, not 46/20 adequately powered controlled experiments**. Many have only one protein or independent group for a label. Exact counts are displayed beside each family row and in `all_training_family_support.tsv` / `stratum_target_support.tsv`.

## Metrics and interpretation

Binary targets use AUROC and average precision; continuous targets use R², MAE and RMSE. Single-class binary metrics are undefined. Constant-target R² is undefined; negative R² is retained. Continuous target SD, minimum and maximum are reported because reducing target variation can reduce R² without increasing absolute error. Hydropathy/charge/length are sequence descriptors; AGGRESCAN outputs are predictions rather than aggregation measurements.

Primary solubility intervals use 2,000 homology-group bootstrap draws from the fixed validation cohort (seed 20261118), preserving groups and using the same draws across layers/baseline. The stratum mask is applied within each draw. An interval is reported only with at least three distinct groups per solubility class and at least 90% valid draws. These are minimum reporting safeguards, not validated confirmatory-power thresholds. Sparse point estimates remain explicitly exploratory. Intervals are conditional on fitted models and omit model-fitting, hyperparameter-selection and multiple-inspection uncertainty. Secondary targets have point estimates only.

ESM-versus-composition differences and intervals compare exactly the same proteins and target values within a stratum. Overall-versus-restricted curves describe different populations: an AUROC change is not a paired treatment effect or a causal estimate of dependence on a property. Range restriction changes difficulty, biological variation, class composition and sampling uncertainty. Inspecting many layers/strata and reporting the largest drop is exploratory. This analysis alone cannot establish direct/indirect encoding or causal computation.

## Reuse, storage and reproduction

The runner requires a complete source under private HF `esm_stage1/uncontrolled_property_probing_v1/`, including all 34 checkpoints (composition + 33 layers) and the final summary. Input/protocol/probe-code hashes must match. An optional explicit `SOURCE_RUN` selects a prefix; automatic selection uses deterministic prefix order among compatible completed runs, never performance. An unfinished source produces a clear stop without training or downloading activations.

Existing local prediction files are reused after source-identity and file-hash checks. Otherwise original probe archives are downloaded and verified by SHA-256; only validation prediction TSVs are extracted. Archives also contain fitted-probe artifacts, but these are not loaded or applied. No ESM model weights, pooled activations or full residue maps are fetched. Every prediction table must have precisely the expected unique IDs, sequence hashes, split and finite scores for that target. Observed outcomes remain in separate frozen local manifests.

`run_within_property.py` orchestrates reuse/evaluation. `within_property_evaluation.py` defines strata, metrics and plots. `build_within_property_notebook.py` creates only the new notebook and preserves the running uncontrolled notebook. CPU dependencies and pinned versions reuse the prior installer; installed pins are checked before any installation. Current software, code hashes, source checkpoint hashes and restriction protocol are recorded in outputs.

Local metric calculations checkpoint by layer. Once complete, tables/plots/provenance are uploaded to the existing private repository under `esm_stage1/within_property_evaluation_v1/<fingerprint>/`. A completed result is reused on reopening. A deleted runtime before final upload can require recomputing cheap metrics, not refitting models. The current uncontrolled run uses a separate work/output directory and remains unchanged. Drive supplies code/frozen inputs; HF persists research results.

The notebook includes solubility range curves, named family heatmaps with protein/group counts, competing-hit sensitivity, continuous-property heatmaps, complete mixed-family support tables and an exact-result viewer for any target/subset. All numeric results are retained, including values beyond heatmap color limits and unsupported tasks. No synthetic performance curves are included in the notebook.

## Validation and current execution status

Eight local tests cover training-only cut points, partition ties/missingness, unique family counting, prediction alignment/corruption rejection, same-population baseline comparison with fitting disabled, incomplete-source handling, selective archive extraction, notebook syntax and all-layer plotting. Temporary synthetic fixtures are engineering tests and are not retained as research results.

`validation_checks/within_property/` records the actual development-only support audit and source hashes. It contains no new property-performance results. Full conditional performance remains pending the authenticated Colab run after uncontrolled probing completes. The reserved NESG test and locked UESolDS independent test were not accessed.

```sh
python -m unittest test_within_property -v
python build_within_property_notebook.py
```

The second command rebuilds the unexecuted new notebook; do not run it over a later executed copy whose outputs you want to keep.
