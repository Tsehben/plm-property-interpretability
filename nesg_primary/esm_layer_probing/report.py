"""Figures and a concise report; never draw placeholder layer-performance curves."""

import matplotlib.pyplot as plt
import pandas as pd
from common import HERE, PROTOCOL, SAE_LAYERS

NAVY = "#182c49"
PURPLE = "#8162ba"
plt.rcParams.update(
    {
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 11,
        "axes.labelcolor": NAVY,
        "text.color": NAVY,
        "svg.fonttype": "none",
        "figure.dpi": 110,
    }
)


def save(fig, name):
    path = HERE / "figures"
    path.mkdir(exist_ok=True)
    fig.tight_layout()
    for ext in ["png", "svg"]:
        fig.savefig(path / f"{name}.{ext}", dpi=220, bbox_inches="tight")
    return fig


def baseline_figure(metrics, cv):
    row = metrics[metrics.model == "composition"].iloc[0]
    folds = cv[(cv.model == "composition") & (cv.C == row.C)].sort_values("fold")
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.scatter(folds.fold, folds.auroc, color=NAVY, label="Training CV (selected C)")
    ax.errorbar(
        [6],
        [row.auroc],
        yerr=[[row.auroc - row.auroc_low], [row.auroc_high - row.auroc]],
        fmt="o",
        color=PURPLE,
        capsize=5,
        label="Fixed validation, group-bootstrap 95% CI",
    )
    ax.set(
        xticks=[0, 1, 2, 3, 4, 6],
        xticklabels=["CV 1", "CV 2", "CV 3", "CV 4", "CV 5", "Validation"],
        ylabel="AUROC",
        ylim=(0.45, 0.9),
    )
    ax.axhline(0.5, color="gray", ls=":", lw=1)
    ax.legend(fontsize=8, loc="upper left")
    return save(fig, "composition_baseline")


def layer_figure(which, metrics, cv, prevalence):
    layers = metrics[metrics.layer > 0].sort_values("layer")
    base = metrics[metrics.model == "composition"].iloc[0]
    assert layers.layer.tolist() == list(range(1, 34)), (
        "All 33 completed layers are required"
    )
    fig, ax = plt.subplots(figsize=(9, 4))
    if which == 1:
        ax.plot(
            layers.layer,
            layers.auroc,
            color=PURPLE,
            marker=".",
            label="ESM linear probe",
        )
        ax.fill_between(
            layers.layer,
            layers.auroc_low,
            layers.auroc_high,
            color=PURPLE,
            alpha=0.18,
            label="Pointwise group-bootstrap 95% CI",
        )
        ax.axhline(base.auroc, color=NAVY, ls="--", label="Composition baseline")
        ax.axhline(0.5, color="gray", ls=":", label="Chance")
        ax.set(
            xlabel="ESM transformer layer",
            ylabel="Validation AUROC",
            xticks=[1, 4, 8, 12, 16, 20, 24, 28, 33],
        )
        name = "layerwise_auroc"
    elif which == 2:
        chosen = cv[(cv.model != "composition") & (cv.C == cv.selected_C)].copy()
        chosen["layer"] = chosen.model.str.split("_").str[-1].astype(int)
        for fold, g in chosen.groupby("fold"):
            g = g.sort_values("layer")
            ax.plot(g.layer, g.auroc, marker=".", lw=1, label=f"CV fold {fold + 1}")
        ax.set(
            xlabel="ESM transformer layer",
            ylabel="Training-CV AUROC",
            xticks=[1, 4, 8, 12, 16, 20, 24, 28, 33],
        )
        name = "training_cv_consistency"
    elif which == 3:
        ax.plot(
            layers.layer,
            layers.average_precision,
            color=PURPLE,
            marker=".",
            label="ESM linear probe",
        )
        ax.fill_between(
            layers.layer,
            layers.ap_low,
            layers.ap_high,
            color=PURPLE,
            alpha=0.18,
            label="Pointwise group-bootstrap 95% CI",
        )
        ax.axhline(
            prevalence,
            color="gray",
            ls=":",
            label=f"Validation soluble prevalence ({prevalence:.3f})",
        )
        ax.axhline(
            base.average_precision, color=NAVY, ls="--", label="Composition baseline"
        )
        ax.set(
            xlabel="ESM transformer layer",
            ylabel="Validation average precision",
            xticks=[1, 4, 8, 12, 16, 20, 24, 28, 33],
        )
        name = "layerwise_average_precision"
    else:
        plt.close(fig)
        fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 4))
        chosen = layers[layers.layer.isin(PROTOCOL["representative_layers"])]
        ax.errorbar(
            chosen.layer.astype(str),
            chosen.auroc,
            yerr=[chosen.auroc - chosen.auroc_low, chosen.auroc_high - chosen.auroc],
            fmt="o",
            capsize=5,
            color=PURPLE,
            label="Validation AUROC",
        )
        ax.axhline(base.auroc, color=NAVY, ls="--")
        ax.set(xlabel="Prespecified representative layer", ylabel="AUROC")
        paired = pd.read_csv(
            HERE / "results/representative_paired_differences.tsv", sep="\t"
        )
        for i, r in paired.iterrows():
            ax2.errorbar(
                r.delta,
                i,
                xerr=[[r.delta - r.low], [r.high - r.delta]],
                fmt="o",
                color=PURPLE,
                capsize=4,
            )
        ax2.axvline(0, color="gray", ls=":")
        ax2.set(
            yticks=range(len(paired)),
            yticklabels=[
                f"L{int(r.layer_a)} − L{int(r.layer_b)}" for r in paired.itertuples()
            ],
            xlabel="Paired validation AUROC difference",
        )
        name = "representative_layer_comparison"
    ax.legend(fontsize=8)
    return save(fig, name)


def recommendation(metrics):
    eligible = metrics[metrics.layer.isin(SAE_LAYERS)]
    picked = []
    for low, high in [(1, 11), (12, 22), (23, 33)]:
        choices = eligible[eligible.layer.between(low, high)]
        if len(choices):
            # Within one CV SD of the highest training-CV mean, prefer stable CV.
            best = choices.loc[choices.cv_mean.idxmax()]
            near = choices[choices.cv_mean >= best.cv_mean - best.cv_sd]
            picked.append(int(near.sort_values(["cv_sd", "layer"]).iloc[0].layer))
    return picked


def write_summary(metrics, cv):
    baseline = metrics[metrics.model == "composition"].iloc[0]
    layers = metrics[metrics.layer > 0]
    status = (
        "Complete"
        if len(layers) == 33
        else "Partial: composition complete; Colab ESM extraction pending"
    )
    layer_text = "No ESM activations, ESM pilot results, or ESM probe metrics have been produced locally. The user chose Colab Secrets for HF_TOKEN_PLM. No best layer or layer-wise pattern can yet be reported. Layers 8, 16 and 32 are prespecified comparison anchors, not evidence-based SAE recommendations."
    if len(layers) == 33:
        best = layers.loc[layers.auroc.idxmax()]
        layer_text = f"Highest observed validation AUROC: layer {int(best.layer)}, {best.auroc:.4f} (pointwise group-bootstrap 95% CI {best.auroc_low:.4f}–{best.auroc_high:.4f}). This maximum is selected across 33 comparisons; its interval is not selection-adjusted. Candidate SAE layers: {recommendation(metrics)}, chosen within early/middle/late bins using training-CV mean and stability, restricted to released InterProt checkpoints. Review paired intervals before interpreting layer differences."
    folds = cv[(cv.model == "composition") & (cv.C == baseline.C)].sort_values("fold")
    summary = f"""# NESG Stage 1 — ESM layer-wise probing

**Status: {status}.** No Stage 2 analysis or interventions were run.

## Verified development data

| Split | Proteins | Soluble | Insoluble | Homology groups |
|---|---:|---:|---:|---:|
| Training | 2,361 | 1,645 | 716 | 1,090 |
| Validation | 691 | 491 | 200 | 364 |

All original IDs, sequence hashes, labels, E/S scores, split assignments and CV folds match frozen safe inputs. All three X-containing proteins are retained; lengths are 50–810, with no truncation. Test count 710 and cross-split group disjointness were checked using the frozen label-free assignment table only. Reserved-test sequences and labels were not accessed. The test was historically inspected in prior family-label reports and remains prospectively locked.

## Actual results

Composition + log-length, selected C = {baseline.C:g}:
- Validation AUROC **{baseline.auroc:.4f}** (homology-group bootstrap 95% CI **{baseline.auroc_low:.4f}–{baseline.auroc_high:.4f}**).
- Average precision **{baseline.average_precision:.4f}** (95% CI {baseline.ap_low:.4f}–{baseline.ap_high:.4f}); soluble prevalence 491/691 = {491 / 691:.4f}.
- Training-CV AUROCs: {", ".join(f"{v:.4f}" for v in folds.auroc)}; mean {baseline.cv_mean:.4f}, SD {baseline.cv_sd:.4f}.

{layer_text}

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

Model: `{PROTOCOL["model"]}`, revision `{PROTOCOL["revision"]}`. The [pinned ESM checkpoint](https://huggingface.co/facebook/esm2_t33_650M_UR50D/commit/{PROTOCOL["revision"]}) and [InterProt model card](https://huggingface.co/liambai/InterProt-ESM2-SAEs/blob/{PROTOCOL["interprot_revision"]}/README.md) establish model identity and matching `hidden_states[LAYER]` indexing. Released SAE layers: {", ".join(map(str, SAE_LAYERS))}. No SAE weights were downloaded.
"""
    (HERE / "ESM_LAYERWISE_SUMMARY.md").write_text(summary)
    return summary


def main():
    metrics = pd.read_csv(HERE / "results/metrics.tsv", sep="\t")
    cv = pd.read_csv(HERE / "results/all_cv.tsv", sep="\t")
    baseline_figure(metrics, cv)
    if (metrics.layer > 0).sum() == 33:
        for i in range(1, 5):
            layer_figure(i, metrics, cv, 491 / 691)
    write_summary(metrics, cv)


if __name__ == "__main__":
    main()
