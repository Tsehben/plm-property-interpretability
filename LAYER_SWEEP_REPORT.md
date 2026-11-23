# Layer Sweep and Cluster-Aware Sensitivity — Report

**Status:** complete. No P0 access, no new interventions, no sequence generation, no
reconstruction pass, not committed, not pushed. Previous outputs preserved.

**Authorised by** [`LOCALISATION_EXPERIMENT_SPEC.md` Amendment A](LOCALISATION_EXPERIMENT_SPEC.md#amendment-a--layer-wise-decodability-is-separated-from-intervention-effects),
which separates layer-wise **decodability** from **intervention effects** and withdraws the
old requirement that layer 8 first show a reproducible specific effect. The original §4.5 text
is preserved verbatim under a banner. Protocol frozen before any outcome in
[`layer_sweep/LAYER_SWEEP_PROTOCOL.md`](layer_sweep/LAYER_SWEEP_PROTOCOL.md) and
[`frozen_protocol.json`](layer_sweep/v1/frozen_protocol.json) (SHA-256 `0d393b69…`).

> **This is not a replication** of the Check-in 1 four-fold layer-8 result. That held each of
> P1–P4 out once; this fits once on P1+P2. Different design, not comparable as a reproduction.
> **P4 is not a fresh final test set** — it is development data that has already carried the
> 400-protein intervention evaluation. Every P4 figure here is a development evaluation.

---

## 1. What was run

| | |
|---|---|
| Cohort | eligible core, `in_core == 1`, **10,104** (10,139 − 8 E1 − 27 E2) |
| Roles | fit P1+P2 = 5,190 · compare P3 = 2,548 · evaluate P4 = 2,366 |
| Classes | 2,491 Extracellular / 7,613 Cytoplasm |
| Layers | **1–33**, `hidden_states[i]` = output of block *i*; `hidden_states[0]` not swept |
| Feature sets | `comp` (41) · `esm` (2,560) · `comp_plus_esm` (2,601) |
| Classifier | `l2`, `C=1.0`, `lbfgs`, `max_iter=5000`, `random_state=42`, no class weighting |
| Preprocessing | `StandardScaler` fitted on P1+P2 only |
| Truncation | 1,022 residues; **1,108** of 10,104 affected |
| Convergence | **0 non-converged fits**; max `n_iter` = 159 of 5,000 |

99 probe fits for the P3 curve, 6 more for P4. **Settings were never altered by layer.**

**Numerical note.** The log carries `overflow / divide by zero / invalid encountered in matmul`
RuntimeWarnings from `sklearn.utils.extmath`. These are a documented pre-existing condition in
this project — [`baseline_localization.py:35`](DeepLoc2/baseline/baseline_localization.py#L35)
records that they occur "even on clean finite" inputs and suppresses them; v2 did not suppress
them, so they surfaced. All outputs were checked: **0 non-finite values across 397 stored
arrays**, predictions in [0, 1], 0 failed convergence.

**Provenance and artifacts.** Per-layer SHA-256 of both pooled arrays for all 33 layers in
[`pooled_embedding_provenance.json`](layer_sweep/v1/pooled_embedding_provenance.json).
Probe coefficients, intercepts and training-fitted scaler parameters for all 73 probes in
[`probe_parameters.npz`](layer_sweep/v1/probe_parameters.npz), tagged `original_run_object` —
captured during fitting, **not reconstructed**.

**Execution.** Extraction ran in two segments (2,212.5 s + 148.6 s) because a disk guard
stopped it at 770/815 batches and it resumed from its checkpoint. The checkpoint/resume path
was verified beforehand on a 48-sequence subset: an interrupted-and-resumed extraction was
**byte-identical** to an uninterrupted one on all 33 layers and both poolings, and all 9
compatibility mutations were refused
([`resume_test_report.json`](layer_sweep/v1/resume_test_report.json)).

---

## 2. P3 layer curve — [`p3_layer_curves.tsv`](layer_sweep/v1/p3_layer_curves.tsv)

AUROC on P3 (n = 2,548; 584 Extracellular / 1,964 Cytoplasm). `comp` is layer-independent at
**0.9590** and is the flat reference line.

| Layer | `esm` | `comp_plus_esm` | | Layer | `esm` | `comp_plus_esm` |
|---|---|---|---|---|---|---|
| 1 | 0.9441 | **0.9439** | | 18 | 0.9795 | 0.9798 |
| 2 | 0.9710 | 0.9708 | | 19 | 0.9832 | 0.9839 |
| 3 | 0.9786 | 0.9786 | | 20 | 0.9837 | 0.9843 |
| 4 | 0.9829 | 0.9828 | | 21 | 0.9841 | 0.9829 |
| 5 | 0.9826 | 0.9821 | | 22 | 0.9832 | 0.9840 |
| 6 | 0.9826 | 0.9825 | | 23 | 0.9859 | 0.9860 |
| 7 | 0.9813 | 0.9810 | | 24 | 0.9876 | 0.9866 |
| **8** | 0.9837 | **0.9838** | | 25 | 0.9865 | 0.9853 |
| 9 | 0.9847 | 0.9840 | | **26** | 0.9893 | **0.9888** ← max |
| 10 | 0.9801 | 0.9796 | | 27 | 0.9851 | 0.9843 |
| 11 | 0.9843 | 0.9841 | | 28 | 0.9840 | 0.9838 |
| 12 | 0.9795 | 0.9792 | | 29 | 0.9886 | 0.9877 |
| 13 | 0.9814 | 0.9822 | | 30 | 0.9881 | 0.9876 |
| 14 | 0.9771 | 0.9782 | | 31 | 0.9885 | 0.9879 |
| 15 | 0.9793 | 0.9800 | | 32 | 0.9885 | 0.9883 |
| 16 | 0.9812 | 0.9823 | | 33 | 0.9882 | 0.9880 |
| 17 | 0.9767 | 0.9765 | | | | |

Shape: layer 1 is **below** the composition baseline (0.9439 vs 0.9590); the curve rises
steeply through layers 2–4, sits on a broad plateau around 0.978–0.984 through the middle of
the network, and reaches its maximum at **layer 26**. Layers 23–33 are all ≥ 0.9838. The whole
range above layer 3 spans 0.0123 AUROC — a shallow landscape, not a sharp peak.

### Selection — frozen before any P4 prediction

[`layer_selection.json`](layer_sweep/v1/layer_selection.json), written and `fsync`'d before
any P4 feature matrix was built:

> **Selected layer 26** (P3 `comp_plus_esm` AUROC **0.988768**). Layer 8 = **0.983810**.
> Criterion: highest P3 `comp_plus_esm` AUROC, lower index breaking exact ties.

The P3 margin of layer 26 over layer 8 is **+0.00496**.

---

## 3. P4 evaluation — [`p4_results.tsv`](layer_sweep/v1/p4_results.tsv)

**Development evaluation.** n = 2,366 (571 Extracellular / 1,795 Cytoplasm), **1,547 clusters**
at 30 % identity, largest `Q2VBP4` = 93 proteins (3.93 %), 22 proteins sharing a cluster with
the fitting set.

### Protein-weighted (primary)

| Layer | Features | AUROC | AP |
|---|---|---|---|
| — | `comp` | 0.9696 | 0.9182 |
| 8 | `esm` | 0.9612 | 0.9398 |
| 8 | `comp_plus_esm` | 0.9628 | 0.9406 |
| 26 | `esm` | **0.9829** | **0.9637** |
| 26 | `comp_plus_esm` | 0.9829 | 0.9624 |

### Equal-cluster-weight (explicitly labelled sensitivity, never primary)

Each protein weighted 1 / (its cluster's size within the evaluated set), so every cluster
contributes total weight 1.

| Layer | Features | AUROC | AP |
|---|---|---|---|
| — | `comp` | 0.9638 | 0.8833 |
| 8 | `esm` | 0.9737 | 0.9394 |
| 8 | `comp_plus_esm` | 0.9746 | 0.9406 |
| 26 | `esm` | 0.9784 | 0.9422 |
| 26 | `comp_plus_esm` | 0.9794 | 0.9420 |

### Paired cluster bootstrap — [`p4_cluster_bootstrap.tsv`](layer_sweep/v1/p4_cluster_bootstrap.tsv)

Whole 30 % clusters drawn with replacement, same resample indices for both models so pairing is
preserved; B = 10,000, seed 42, percentile 95 % CI. **Mixed-class clusters kept intact and never
split**, so resampling is unstratified. **0 resamples were undefined** in any comparison.

| Subset | Comparison | ΔAUROC | 95 % CI | ≠0 | ΔAP | 95 % CI | ≠0 |
|---|---|---|---|---|---|---|---|
| All | L26 vs L8 (`comp_plus_esm`) | +0.0201 | [−0.0035, +0.0630] | no | +0.0219 | [−0.0066, +0.0657] | no |
| All | L8 `comp_plus_esm` vs `comp` | −0.0068 | [−0.0446, +0.0178] | no | +0.0224 | [−0.0131, +0.0723] | no |
| All | **L26 `comp_plus_esm` vs `comp`** | **+0.0133** | **[+0.0046, +0.0237]** | **yes** | **+0.0442** | **[+0.0201, +0.0849]** | **yes** |
| All | L8 `esm` vs `comp` | −0.0084 | [−0.0486, +0.0172] | no | +0.0216 | [−0.0143, +0.0708] | no |
| All | **L26 `esm` vs `comp`** | **+0.0133** | **[+0.0047, +0.0237]** | **yes** | **+0.0455** | **[+0.0209, +0.0865]** | **yes** |
| Excl. fitting-set overlap (n=2,344) | L26 vs L8 | +0.0213 | [−0.0024, +0.0644] | no | +0.0222 | [−0.0071, +0.0667] | no |
| Excl. fitting-set overlap | **L26 `comp_plus_esm` vs `comp`** | **+0.0147** | **[+0.0064, +0.0249]** | **yes** | **+0.0451** | **[+0.0208, +0.0855]** | **yes** |
| Excl. fitting-set overlap | L8 `comp_plus_esm` vs `comp` | −0.0065 | [−0.0451, +0.0181] | no | +0.0229 | [−0.0136, +0.0737] | no |
| Excl. toxin cluster (n=2,273) | L26 vs L8 | +0.0237 | [−0.0045, +0.0696] | no | +0.0246 | [−0.0102, +0.0755] | no |
| Excl. toxin cluster | **L26 `comp_plus_esm` vs `comp`** | **+0.0157** | **[+0.0054, +0.0258]** | **yes** | **+0.0654** | **[+0.0402, +0.0926]** | **yes** |
| Excl. toxin cluster | L8 `comp_plus_esm` vs `comp` | −0.0080 | [−0.0513, +0.0199] | no | +0.0408 | [−0.0078, +0.0817] | no |

The toxin-cluster row is a **sensitivity analysis, not a replacement primary result.**

---

## 4. Retrospective intervention reanalysis — [`intervention_cluster_sensitivity.tsv`](layer_sweep/v1/intervention_cluster_sensitivity.tsv)

Saved outcomes only; no model rerun; `evaluation_p4/` untouched. Original eligible cases,
aggregation, BCa intervals and the four-comparison Bonferroni level (98.75 %) all retained.
Resampling unit changed from protein to 30 % cluster, whole clusters, unstratified.

| Contrast | n | Clusters | Largest cluster | Original 98.75 % ≠0 | **Cluster 98.75 % ≠0** |
|---|---|---|---|---|---|
| 2564 × 0.0 | 61 | 38 | `Q2VBP4` 22 (36.1 %) | **yes** | **no** |
| 2564 × 2.0 | 136 | 86 | `Q2VBP4` 31 (22.8 %) | no | no |
| 3900 × 0.0 | 275 | 243 | `O00755` 10 (3.6 %) | no | no |
| 3900 × 2.0 | 301 | 266 | `O00755` 10 (3.3 %) | no | no |

> **None of the four adjusted cluster intervals excludes zero.** The only contrast that was
> significant under the original adjusted BCa interval, 2564 × 0.0, is not under cluster
> resampling: its adjusted interval becomes **[−1.9 × 10⁻⁶, 3.5 × 10⁻⁵]**. The 22 `Q2VBP4`
> proteins in it are 36.1 % of its cases and contribute **73.4 % of its point estimate**.
> At the unadjusted 95 % level that one contrast still excludes zero
> ([4.7 × 10⁻⁷, 3.3 × 10⁻⁵]); the adjusted level is the one the original analysis
> prespecified. **The interpretation changes: the retrospective analysis does not support a
> differential-sensitivity effect for any of the four contrasts.** This is a sensitivity
> analysis; the original intervals stand alongside it.

---

## 5. Conclusions, kept separate

### 5.1 Where localisation is most decodable under this probe protocol

Most decodable in the **upper-middle to late layers**. On P3 the maximum is **layer 26**
(0.9888), and layers 23–33 all sit at ≥ 0.9838. The landscape is **broad, not peaked**: above
layer 3 the entire 31-layer range spans 0.0123 AUROC. Decodability rises steeply over layers
1–4 — layer 1 is *below* the composition baseline — then plateaus.

**Layer 26 is the argmax of the frozen P3 criterion. It is not established as better than
layer 8.** On P4 the L26-vs-L8 difference is +0.0201 AUROC with a 95 % cluster-bootstrap CI of
[−0.0035, +0.0630], which includes zero, and it still includes zero in both sensitivity
subsets. The honest statement is that the sweep *locates* a maximum under a prespecified rule,
not that it *resolves* one layer over another.

### 5.2 Whether ESM adds information beyond composition

**Layer-dependent, and that is the main finding.**

- **At layer 26: yes.** `esm` beats `comp` by **+0.0133 AUROC [+0.0047, +0.0237]** and
  **+0.0455 AP [+0.0209, +0.0865]**, both excluding zero, and both hold in the
  fitting-set-overlap and toxin-excluded subsets (AP strengthens to +0.0670 without the toxin
  cluster). Adding composition to the layer-26 embedding changes almost nothing
  (0.9829 either way), i.e. the embedding already carries what composition supplies.
- **At layer 8: not established.** The point estimate is *negative* for AUROC
  (−0.0084 [−0.0486, +0.0172]) and positive for AP (+0.0216 [−0.0143, +0.0708]); both CIs
  include zero. Under this single-fit design, a layer-8 ESM probe does not beat a
  41-feature amino-acid-composition baseline on P4 AUROC.

This does **not** contradict the Check-in 1 four-fold result; it is a different fitting and
evaluation design (one fit on P1+P2 versus four folds) and the two are not comparable as a
reproduction.

### 5.3 How measured clustering affects performance and intervention uncertainty

**The weighting choice changes a conclusion.** Moving from protein-weighted to
equal-cluster-weight *reverses the sign* of the layer-8 ESM-versus-composition comparison: the
composition baseline falls (AUROC 0.9696 → 0.9638; AP 0.9182 → 0.8833) while layer-8 `esm`
rises (0.9612 → 0.9737). The composition baseline is the model that benefits most from
redundant clusters being counted protein by protein. Both weightings are reported; neither is
the "corrected" one, because they answer different questions — per protein, or per family-like
group.

**One cluster carries a visible share of the baseline.** Removing `Q2VABP4`'s 93 proteins
(3.93 % of P4, all extracellular) drops the composition baseline's AP from 0.9182 to 0.8808,
and widens layer 26's AP advantage from +0.0442 to +0.0654.

**Cross-partition overlap is negligible here.** Only 22 of 2,366 P4 proteins share a 30 %
cluster with the fitting set, and excluding them moves nothing materially.

**For the interventions, cluster-aware resampling removes the one positive result** (§4).

> **Three limits on all of the above.** (i) These clusters are an *operational grouping*;
> resampling them accounts for that grouping without establishing independence between
> clusters, and a cluster interval is not guaranteed to be wider than a protein-level one.
> (ii) **The best probe layer is not thereby the best intervention layer** — decodability and
> intervention effects are different measurements, and nothing here nominates an intervention
> layer. (iii) **No result here establishes direct encoding.** A linear probe recovers what is
> linearly available in a pooled vector, which is equally consistent with a direct
> representation, a diffuse compositional correlate, or something that merely co-varies with
> localisation such as protein family. §1.4 of the specification applies unchanged.

---

## 6. Artifacts

| File | Contents |
|---|---|
| [`frozen_protocol.json`](layer_sweep/v1/frozen_protocol.json) · [`LAYER_SWEEP_PROTOCOL.md`](layer_sweep/LAYER_SWEEP_PROTOCOL.md) | protocol, frozen pre-run |
| [`pooled_embedding_provenance.json`](layer_sweep/v1/pooled_embedding_provenance.json) | model, versions, dtypes, pooling, truncation, per-layer SHA-256 ×33, 16-field compat key |
| [`p3_layer_curves.tsv`](layer_sweep/v1/p3_layer_curves.tsv) | 99 rows: AUROC, AP, `n_iter`, convergence |
| [`p3_predictions.npz`](layer_sweep/v1/p3_predictions.npz) | all-layer P3 predictions (99 arrays) |
| [`layer_selection.json`](layer_sweep/v1/layer_selection.json) | frozen selection |
| [`p4_results.tsv`](layer_sweep/v1/p4_results.tsv) · [`p4_results.json`](layer_sweep/v1/p4_results.json) | both weightings × 3 subsets |
| [`p4_cluster_bootstrap.tsv`](layer_sweep/v1/p4_cluster_bootstrap.tsv) | 15 paired comparisons |
| [`p4_predictions.npz`](layer_sweep/v1/p4_predictions.npz) | P4 predictions, cluster ids, overlap flags |
| [`probe_parameters.npz`](layer_sweep/v1/probe_parameters.npz) | 73 probes: coef, intercept, scaler mean/scale — **original-run objects** |
| [`intervention_cluster_sensitivity.tsv`](layer_sweep/v1/intervention_cluster_sensitivity.tsv) | retrospective reanalysis |
| [`resume_test_report.json`](layer_sweep/v1/resume_test_report.json) | checkpoint/resume verification |
| [`resource_check.json`](layer_sweep/v1/resource_check.json) · [`resource_plan.json`](layer_sweep/v1/resource_plan.json) · [`resource_recheck.json`](layer_sweep/v1/resource_recheck.json) · [`resource_stop_report.json`](layer_sweep/v1/resource_stop_report.json) · [`stop_record.json`](layer_sweep/v1/stop_record.json) · [`stop_reconciliation.json`](layer_sweep/v1/stop_reconciliation.json) | resource history |
| `emb_store/` | 33 per-layer memmaps, 1.707 GB, 815/815 batches committed |
