"""Plot only completed NESG probes, with explicit training evaluation definitions."""

import argparse
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def plot(results):
    results = Path(results)
    table = pd.read_csv(results / "metrics.tsv", sep="\t")
    layers = table[table.layer > 0].sort_values("layer")
    assert layers.layer.tolist() == list(range(1, 34)), (
        "All 33 probes must finish before plotting."
    )
    assert layers.validation_n.eq(691).all() and layers.training_n.eq(2361).all()
    base = table[table.model == "composition"].iloc[0]
    figures = results / "figures"
    figures.mkdir(exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.fonttype": "none",
        }
    )
    for training, filename, title in [
        (
            "cv_mean",
            "train_validation_auroc",
            "Training cross-validation and fixed validation",
        ),
        (
            "training_fit_auroc",
            "train_fit_validation_auroc",
            "Training fit and fixed validation (overfitting diagnostic)",
        ),
    ]:
        fig, ax = plt.subplots(figsize=(10, 5.5))
        x = layers.layer.to_numpy()
        label = (
            "Train: mean held-out-fold AUROC"
            if training == "cv_mean"
            else "Train: fitted-model AUROC (in-sample)"
        )
        ax.plot(x, layers[training], "o-", ms=3.5, color="#236f95", label=label)
        if training == "cv_mean":
            ax.fill_between(
                x,
                (layers.cv_mean - layers.cv_sd).clip(0, 1),
                (layers.cv_mean + layers.cv_sd).clip(0, 1),
                color="#236f95",
                alpha=0.12,
                label="Train: ±1 SD across 5 homology-aware folds",
            )
        ax.plot(
            x,
            layers.auroc,
            "o-",
            ms=3.5,
            color="#8555ac",
            label="Validation: 691 proteins",
        )
        ax.fill_between(
            x,
            layers.auroc_low,
            layers.auroc_high,
            color="#8555ac",
            alpha=0.17,
            label="Validation: 95% homology-group bootstrap CI",
        )
        ax.axhline(
            base.auroc,
            color="#8555ac",
            ls="--",
            lw=1,
            label=f"Composition validation: {base.auroc:.3f}",
        )
        ax.axhline(
            base[training],
            color="#236f95",
            ls="--",
            lw=1,
            label=f"Composition train: {base[training]:.3f}",
        )
        ax.axhline(0.5, color="#888888", ls=":", lw=1, label="Chance")
        ax.set(
            xlabel="ESM-2 transformer layer",
            ylabel="AUROC",
            xlim=(0.5, 33.5),
            ylim=(0, 1.02),
            xticks=[1, 4, 8, 12, 16, 20, 24, 28, 33],
            title=title,
        )
        ax.legend(fontsize=8, loc="lower right")
        note = (
            "Train: 2,361 proteins; C selected using these same folds (not nested CV). "
            "Validation never used to fit or tune probes.\n"
            "Shading has different meanings: training fold SD versus conditional validation CI. Reserved test not evaluated."
            if training == "cv_mean"
            else "Training-fit AUROC is optimistic: these proteins fitted the probe. Use the CV figure to assess generalization.\n"
            "Validation CI resamples homology groups; it does not include training or layer-selection uncertainty."
        )
        fig.text(0.08, 0.01, note, fontsize=8, color="#444444")
        fig.tight_layout(rect=(0, 0.075, 1, 1))
        for extension in ("png", "svg"):
            fig.savefig(
                figures / f"{filename}.{extension}", dpi=250, bbox_inches="tight"
            )
        plt.close(fig)
    columns = [
        "layer",
        "C",
        "training_n",
        "cv_mean",
        "cv_sd",
        "training_oof_auroc",
        "training_fit_auroc",
        "validation_n",
        "auroc",
        "auroc_low",
        "auroc_high",
        "average_precision",
    ]
    layers[columns].to_csv(
        results / "train_validation_auroc.tsv", sep="\t", index=False
    )
    best = layers.loc[layers.auroc.idxmax()]
    summary = (
        "# NESG layer-wise solubility AUROC\n\n"
        "All 33 frozen ESM-2 layer probes completed. Training: 2,361; validation: 691; reserved test not evaluated.\n\n"
        f"Highest observed validation AUROC: layer {int(best.layer)}, {best.auroc:.4f} "
        f"(group-bootstrap 95% CI {best.auroc_low:.4f}–{best.auroc_high:.4f}). "
        "This is a descriptive maximum across 33 layers, not an independently confirmed best layer.\n\n"
        f"Composition baseline: train CV {base.cv_mean:.4f} ± {base.cv_sd:.4f} fold SD; validation {base.auroc:.4f}.\n\n"
        f"Layer-wise train CV AUROC range: {layers.cv_mean.min():.4f}–{layers.cv_mean.max():.4f}; "
        f"validation range: {layers.auroc.min():.4f}–{layers.auroc.max():.4f}.\n\n"
        "The main train curve averages the five held-out-fold AUROCs at training-selected C. "
        "Because those folds also select C, these are tuning estimates, not unbiased nested-CV estimates. "
        "Pooled OOF AUROC is included in the table, but can be affected by fold-to-fold score calibration. "
        "The second figure shows explicitly in-sample training fit as an overfitting diagnostic. "
        "Natural class/family distributions and frozen splits are preserved. AUROC establishes linear decodability, not direct encoding or causality.\n"
    )
    (results / "TRAIN_VALIDATION_AUROC.md").write_text(summary)
    print(summary)
    return layers[columns]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    plot(parser.parse_args().results)
