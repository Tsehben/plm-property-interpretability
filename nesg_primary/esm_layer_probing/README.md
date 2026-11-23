# NESG ESM-2 analysis workflows

Use the notebooks below in order. Saved outputs confirm completion of extraction, all 33 solubility probes and the uncontrolled property probes. Within-property evaluation is the next workflow; its notebook began receiving execution outputs through Drive sync during cleanup. This status is based on the saved Colab outputs; private HF artifacts were not downloaded again during repository cleanup.

| Notebook | Purpose | Status |
|---|---|---|
| [ESM_Extraction_Colab.ipynb](ESM_Extraction_Colab.ipynb) | Extract frozen ESM-2 representations; GPU runtime | Executed; pooled activations and float16 full maps completed |
| [ESM_Train_Validation_AUROC_Colab.ipynb](ESM_Train_Validation_AUROC_Colab.ipynb) | Solubility probes across all 33 layers; CPU runtime | Executed, with saved figures and metrics |
| [ESM_Uncontrolled_Property_Probes_Colab.ipynb](ESM_Uncontrolled_Property_Probes_Colab.ipynb) | Solubility, five continuous properties and 21 Pfam membership targets; CPU runtime | Executed; compatible solubility probes reused |
| [ESM_Within_Property_Evaluation_Colab.ipynb](ESM_Within_Property_Evaluation_Colab.ipynb) | Evaluate saved global predictions within property ranges and families; CPU runtime | Execution started; completion not yet verified |

The solubility notebook reports its highest observed validation AUROC at layer 32: **0.8146**, with homology-group bootstrap 95% CI **0.7697–0.8567**. This is a descriptive maximum over 33 layers, not an independently confirmed best layer. Composition-baseline AUROC is **0.7216**. Consult the executed notebooks for all layers and uncertainty.

## Open or resume in Colab

1. Open the relevant notebook from the existing Drive project and choose **Run all**. Use a GPU only for extraction or full-map numerical checks; the probing/evaluation notebooks use CPU.
2. Enable notebook access to **HF_TOKEN_PLM** in Colab Secrets. Do not put the credential in code, files or outputs.
3. Keep the configured private dataset **PLM-Interp/nesg-esm2-layerwise**. Existing-repository checks verify privacy. Repository creation is only attempted through explicit configuration; use `CREATE_PRIVATE_REPO=False` for an existing dataset.
4. Let setup copy current code and required development inputs from Drive into temporary `/content` storage. There is no manual ZIP upload.
5. Matching completed results and extraction shards are checked and reused. A new runtime may need downloads and dependency setup again; it does not require repeating compatible completed experiments.

HF is the persistent activation/result store. Drive holds project code, manifests and saved notebooks. Runtime caches under `/content` are disposable. The standalone solubility runner uploads its completed run as one bundle; the uncontrolled runner checkpoints each layer. Within-property evaluation checks the completed uncontrolled run before using its predictions and downloads no ESM activations.

## Methods and frozen inputs

Read [UNCONTROLLED_PROPERTY_METHODS.md](UNCONTROLLED_PROPERTY_METHODS.md) and [WITHIN_PROPERTY_METHODS.md](WITHIN_PROPERTY_METHODS.md) for target definitions, eligibility, missingness, subset construction, uncertainty and interpretation.

- Development cohort: **2,361 training + 691 validation** proteins. Labels, sequences, hashes, homology groups and folds are frozen. All three X-containing proteins remain in extraction; undefined biochemical annotations remain missing.
- `data/extraction_manifest.tsv` contains development IDs, hashes, splits and sequences. `data/local_labels.tsv` contains experimental labels and homology/CV links. Pfam memberships preserve multiple assignments and competing-overlap information.
- Reserved-test sequences and labels are excluded from these pipelines. The test's family-label composition was inspected historically; its existing lock and assignments are retained.
- Model: `facebook/esm2_t33_650M_UR50D`, pinned revision `08e4846e537177426273712802403f7ba8261b6c`. No fine-tuning.
- `hidden_states[1]` through `[33]` are the corresponding transformer outputs; `[0]` is the embedding output. Layer 33 includes the final encoder normalization. Layers are not concatenated.
- Pooled representations have shape `[3052, 33, 1280]`, stored in float32. Residue pooling excludes BOS, EOS and padding. Original lengths are 50–810; no truncation.
- Solubility probes use fold-fitted standardization and L2 logistic regression, selecting C from `[0.0001, 0.001, 0.01, 0.1, 1, 10]` using the five frozen homology-aware training folds. Final fits use training proteins only. Validation does not tune C.
- Reported training CV reuses hyperparameter-selection folds; it is not nested CV. Validation intervals resample homology groups and condition on fitted models. They do not include training or layer-selection uncertainty.
- Within-property evaluation uses the same global probes on restricted validation populations. A change in AUROC or R² alone does not establish causal confounding or direct encoding.

## Code map

| Files | Responsibility |
|---|---|
| `common.py`, `prepare_data.py`, `protocol.json` | Frozen inputs, identities, shared utilities and original solubility protocol |
| `extract.py`, `install_colab_dependencies.py`, `requirements-colab.txt` | Extraction, numerical pilot, isolated Colab GPU environment |
| `hf_storage.py` | Private repository checks, validated shard uploads and restoration |
| `save_full_maps.py`, `pilot_full_map_precision.py` | Full residue maps, storage-precision checks and shared HF retry helper |
| `probe.py`, `plot_train_validation.py`, `run_layerwise_auroc.py` | Solubility fitting, plots, CPU environment and reusable result archives |
| `property_inputs.py`, `property_probes.py`, `property_report.py`, `run_property_probes.py` | Uncontrolled targets, fitting, reporting and per-layer checkpoints |
| `within_property_evaluation.py`, `run_within_property.py` | Subset support, saved-prediction evaluation and private result storage |
| `build_notebooks.py`, `build_*_notebook.py` | Notebook sources and shared construction/setup helpers |
| `test_*.py` | Numerical, alignment, access, resumption and workflow regression checks |

Notebook builders overwrite outputs when executed. Do not run them merely to reopen an analysis. Builders import shared helpers from one another; they are retained intentionally. Computational source hashes participate in HF cache identity, so cosmetic edits to those scripts can invalidate result reuse.

## Full residue maps and precision

`save_full_maps.py` writes `[33, sequence_length, 1280]` tensors per protein, excluding special tokens. Float16 storage totals approximately 63.5 GB; float32 would double that. Maps are separate from the approximately 516 MB pooled float32 arrays. The completed run used float16 storage.

The recorded numerical pilot tested 11 proteins at nine pretrained SAE layers and passed its configured tolerances: worst protein/layer mean active-set Jaccard **0.998480**, maximum latent relative L2 **0.010227**. This supports the working storage format on that sample; it does not establish biological equivalence for every feature or protein. Load maps as float32 for the float32 SAE. No SAE discovery or intervention analysis has been performed here.

Released InterProt SAE layers are **4, 8, 12, 16, 20, 24, 28, 32, 33**. The precision pilot records checkpoint/source checksums and uses the reference normalization and top-k settings. The all-layer probes and within-property evaluation retain all 33 ESM layers independently of SAE checkpoint availability.

## Storage, recovery and local reproduction

Private HF paths:

- `esm_stage1/`: pooled shards, checksums and extraction metadata.
- `esm_stage1/full_maps_float16_v1/`: per-protein residue maps and checkpoint manifests.
- `esm_stage1/precision_pilot/`: numerical pilot reports.
- `esm_stage1/layerwise_probing_v1/<fingerprint>/`: complete solubility result bundles.
- `esm_stage1/uncontrolled_property_probing_v1/<fingerprint>/`: composition/layer checkpoints and final report.
- `esm_stage1/within_property_evaluation_v1/<fingerprint>/`: within-property results once executed.

Published prediction bundles omit per-protein experimental labels and sequences; rejoin the local manifest by stable ID, sequence hash and split. Uploads require an accessible private destination and dataset write permission. Model-download access alone is insufficient.

For access diagnostics, run `hf_storage.py diagnose --repo PLM-Interp/nesg-esm2-layerwise` with the secret supplied through the environment. Errors redact raw responses and credentials. HTTP 429 triggers bounded retries honoring server reset hints; retain the same destination and arguments when resuming. Full-map batches are skipped only after their remote checkpoints validate. A transfer progress bar reaching 100% does not itself prove the repository commit completed.

Local Python 3.12 dependencies are recorded in `requirements.txt`. Colab uses its dedicated installers and preserves compatible runtime packages. For local pooled reuse:

```sh
python nesg_primary/esm_layer_probing/hf_storage.py download --repo PLM-Interp/nesg-esm2-layerwise --folder /tmp/nesg_activations
python nesg_primary/esm_layer_probing/extract.py validate --output /tmp/nesg_activations
```

Supply `HF_TOKEN_PLM` securely before authenticated commands. Use each runner's `--help` for explicit local paths and options; the notebooks provide the complete configured workflows.

## Historical local baseline

[ESM_Layerwise_Solubility.ipynb](ESM_Layerwise_Solubility.ipynb), [ESM_LAYERWISE_SUMMARY.md](ESM_LAYERWISE_SUMMARY.md), `results/` and `figures/composition_baseline.*` preserve the initial local baseline execution. Its pending-ESM statements describe that earlier run, not the current Colab status. `report.py`, `execute_analysis.py` and `audit_outputs.py` retain its local reproduction path. These are actual scientific outputs, not discarded plotting drafts.

Repository cleanup is recorded in [CLEANUP_2026-11-22.json](../../CLEANUP_2026-11-22.json). Executed notebooks, current pipeline source, frozen data and HF artifacts were preserved.
