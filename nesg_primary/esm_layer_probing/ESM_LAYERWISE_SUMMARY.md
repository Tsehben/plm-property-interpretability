> **Historical local baseline snapshot (2026-11-18).** Pending-extraction statements below describe that run. Extraction and the 33-layer Colab probes have since completed; see the [current workflow index](README.md) and [executed AUROC notebook](ESM_Train_Validation_AUROC_Colab.ipynb). The original results below are retained for provenance.

# NESG Stage 1 — ESM layer-wise probing

**Status: Partial: composition complete; Colab ESM extraction pending.** No Stage 2 analysis or interventions were run.

## Verified development data

| Split | Proteins | Soluble | Insoluble | Homology groups |
|---|---:|---:|---:|---:|
| Training | 2,361 | 1,645 | 716 | 1,090 |
| Validation | 691 | 491 | 200 | 364 |

All original IDs, sequence hashes, labels, E/S scores, split assignments and CV folds match frozen safe inputs. All three X-containing proteins are retained; lengths are 50–810, with no truncation. Test count 710 and cross-split group disjointness were checked using the frozen label-free assignment table only. Reserved-test sequences and labels were not accessed. The test was historically inspected in prior family-label reports and remains prospectively locked.

## Actual results

Composition + log-length, selected C = 0.1:
- Validation AUROC **0.7216** (homology-group bootstrap 95% CI **0.6669–0.7739**).
- Average precision **0.8406** (95% CI 0.7814–0.8919); soluble prevalence 491/691 = 0.7106.
- Training-CV AUROCs: 0.6182, 0.7353, 0.7112, 0.6980, 0.7176; mean 0.6961, SD 0.0456.

No ESM activations, ESM pilot results, or ESM probe metrics have been produced locally. The user chose Colab Secrets for HF_TOKEN_PLM. No best layer or layer-wise pattern can yet be reported. Layers 8, 16 and 32 are prespecified comparison anchors, not evidence-based SAE recommendations.

CV scores reuse the folds used to select C and therefore are tuning estimates, not unbiased nested-CV estimates. Bootstrap intervals resample validation homology groups, retain all group members, and condition on fitted models; they do not include training variability. The same draws are used for paired model comparisons. Family correlations are preserved. No decodability score alone establishes direct or indirect encoding or causal computation.

## Execution and remaining work

The analysis notebook is executed for available results and explicitly marks unavailable layer plots. The Colab notebook is a runnable extraction workflow, not an executed model experiment. Engineering tests cover hash alignment, pooling, X tokenization, group bootstrap, shard corruption and private-repository safeguards; synthetic tests are not ESM correctness-pilot results.

Run the Colab notebook with HF_TOKEN_PLM in Colab Secrets: pinned-model FP32/mixed-precision pilot, all-layer extraction, and complete-shard validation. Each shard uploads immediately to an explicitly configured existing private Hugging Face dataset; Google Drive supplies existing code and input sequences only, while activation storage stays on HF. Restore interrupted extraction from HF, or download verified shards to a disposable working cache and run `probe.py --activations PATH`, then rerun the analysis notebook. Its four ESM figures and this report populate only from actual completed probes. No HF upload has occurred in the local run, and the configured destination `PLM-Interp/nesg-esm2-layerwise` has not been privacy-verified. The HF-only Colab workflow requires an existing private destination, refuses public or incompatible destinations, and creates repositories only through an explicit `init --create-private` request (disabled by default).

## Files and provenance

- [Executed analysis notebook](ESM_Layerwise_Solubility.ipynb)
- [Colab extraction notebook](ESM_Extraction_Colab.ipynb)
- [Reproduction instructions](README.md)
- [Metrics](results/metrics.tsv), [all CV candidates](results/all_cv.tsv), [baseline predictions](results/composition_validation_predictions.tsv)
- [Protein labels/group links](data/local_labels.tsv), [multi-family/domain annotations](data/pfam_domains.tsv.gz)
- Activations: canonical copy in the configured private HF dataset; Colab uses `/content/nesg_esm_stage1_scratch/activations/` as temporary processing storage only. No extraction/upload has run yet.

Model: `facebook/esm2_t33_650M_UR50D`, revision `08e4846e537177426273712802403f7ba8261b6c`. The [pinned ESM checkpoint](https://huggingface.co/facebook/esm2_t33_650M_UR50D/commit/08e4846e537177426273712802403f7ba8261b6c) and [InterProt model card](https://huggingface.co/liambai/InterProt-ESM2-SAEs/blob/3ded2e51641df84d6b89f25430874e3a8eb42c24/README.md) establish model identity and matching `hidden_states[LAYER]` indexing. Released SAE layers: 4, 8, 12, 16, 20, 24, 28, 32, 33. No SAE weights were downloaded.
