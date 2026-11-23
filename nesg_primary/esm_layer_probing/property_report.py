"""Summarize actual complete uncontrolled runs; no synthetic/fallback figures."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import write_json
from property_inputs import CONTINUOUS, SPEC

NAMES = dict(hydrophobicity="Mean Kyte–Doolittle hydropathy", charge_proxy="Charge proxy (K+R−D−E)/L",
             log_length="Natural log sequence length", aggregation_Na4vSS="AGGRESCAN Na4vSS (predicted)",
             apr_fraction="Predicted APR residue fraction", solubility="Experimental soluble recovery")


def render(root, catalog):
    root = Path(root)
    targets = ["solubility", *CONTINUOUS, *catalog.loc[catalog.eligible, "pfam"]]
    rows, deltas, constants = [], [], []
    family_boots = {}
    for layer in range(34):
        directory = root/f"layer_{layer:02d}"
        for target in targets:
            folder = directory/target
            table = pd.read_csv(folder/"metrics.tsv", sep="\t")
            assert len(table) == (3 if target in CONTINUOUS else 2)
            table["target"], table["layer"] = target, layer
            rows.append(table)
            con = pd.read_csv(folder/"constant.tsv", sep="\t")
            con["target"], con["layer"] = target, layer
            constants.append(con)
            with np.load(folder/"bootstrap.npz") as boot:
                if target.startswith("PF"):
                    family_boots[layer, target] = {k: boot[k].copy() for k in boot.files}
                if layer:
                    with np.load(root/"layer_00"/target/"bootstrap.npz") as baseline:
                        for metric in boot.files:
                            difference = boot[metric]-baseline[metric]
                            good = difference[np.isfinite(difference)]
                            lo, hi = np.quantile(good, [.025, .975]) if len(good) else (np.nan, np.nan)
                            point = table.set_index("metric").loc[metric, "validation"]
                            base = pd.read_csv(root/"layer_00"/target/"metrics.tsv", sep="\t").set_index("metric").loc[metric, "validation"]
                            deltas.append(dict(layer=layer, target=target, metric=metric, delta=point-base,
                                               low=lo, high=hi, valid_bootstrap=len(good)))
    metrics = pd.concat(rows, ignore_index=True)
    metrics.to_csv(root/"metrics.tsv", sep="\t", index=False)
    pd.concat(constants).to_csv(root/"constant_baselines.tsv", sep="\t", index=False)
    pd.DataFrame(deltas).to_csv(root/"paired_vs_composition.tsv", sep="\t", index=False)
    # Validate per-target population consistency before displaying any curves.
    assert metrics.groupby("target").train_n.nunique().eq(1).all()
    assert metrics.groupby("target").validation_n.nunique().eq(1).all()
    family_rows = metrics[metrics.target.str.startswith("PF")]
    eligible = catalog[catalog.eligible].copy()
    support = family_rows[family_rows.layer.eq(0) & family_rows.metric.eq("auroc")]
    eligible = eligible.merge(support[["target", "validation_positive", "validation_n", "validation"]],
                              left_on="pfam", right_on="target", validate="one_to_one")
    eligible["validation_negative"] = eligible.validation_n-eligible.validation_positive
    eligible["validation_estimable"] = np.isfinite(eligible.validation)
    eligible.to_csv(root/"family_support.tsv", sep="\t", index=False)
    evaluable = eligible.loc[eligible.validation_estimable, "pfam"].tolist()
    macro = []
    for layer in range(34):
        for metric in ("auroc", "average_precision"):
            table = family_rows[family_rows.layer.eq(layer) & family_rows.metric.eq(metric)].set_index("target").loc[evaluable]
            if len(evaluable):
                draws = np.stack([family_boots[layer, t][metric] for t in evaluable])
                valid = np.isfinite(draws).all(axis=0)
                samples = draws[:, valid].mean(axis=0)
                adequate = valid.mean() >= SPEC["family_macro_ci_min_valid_fraction"]
                lo, hi = np.quantile(samples, [.025, .975]) if adequate else (np.nan, np.nan)
            else:
                samples = np.array([])
                lo = hi = np.nan
            macro.append(dict(layer=layer, metric=metric, validation=table.validation.mean(), low=lo, high=hi,
                              cv_mean=table.cv_mean.mean(), families=len(evaluable), bootstrap_valid=len(samples),
                              ci_reported=bool(np.isfinite(lo))))
    macro = pd.DataFrame(macro)
    macro.to_csv(root/"family_macro.tsv", sep="\t", index=False)
    figdir = root/"figures"
    figdir.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

    def save(fig, name):
        fig.tight_layout()
        for ext in ("png", "svg"):
            fig.savefig(figdir/f"{name}.{ext}", dpi=180, bbox_inches="tight")
        plt.close(fig)

    def curves(ax, target, metric):
        sub = metrics[(metrics.target == target) & (metrics.metric == metric)].sort_values("layer")
        base, d = sub[sub.layer.eq(0)].iloc[0], sub[sub.layer.gt(0)]
        ax.plot(d.layer, d.cv_mean, color="#24718d", label="Train: mean held-out-fold score")
        ax.fill_between(d.layer, d.cv_mean-d.cv_sd, d.cv_mean+d.cv_sd, color="#24718d", alpha=.10)
        ax.plot(d.layer, d.validation, color="#8054a6", label="Validation")
        ax.fill_between(d.layer, d.low, d.high, color="#8054a6", alpha=.17)
        ax.axhline(base.validation, color="#8054a6", ls="--", label="Composition: validation")
        ax.axhline(base.cv_mean, color="#24718d", ls="--", label="Composition: train CV")
        ax.set(xlabel="ESM-2 transformer layer", ylabel=metric.upper(), xticks=[1,8,16,24,33])
        ax.set_title(f"{NAMES[target]}\ntrain n={int(base.train_n):,}; validation n={int(base.validation_n):,}")

    for metric in ("r2", "mae"):
        fig, axes = plt.subplots(3,2,figsize=(13,12))
        for ax, target in zip(axes.flat, CONTINUOUS):
            curves(ax,target,metric)
            if metric == "r2":
                ax.axhline(0, color="gray", ls=":")
        axes.flat[-1].axis("off")
        handles, labels = axes.flat[0].get_legend_handles_labels()
        axes.flat[-1].legend(handles,labels,loc="center",frameon=False)
        fig.suptitle("Uncontrolled continuous-property decoding" + (" — lower MAE is better" if metric == "mae" else " — R² can be negative"))
        save(fig, f"continuous_{metric}")
    fig, ax = plt.subplots(figsize=(11,5))
    curves(ax,"solubility","auroc")
    ax.axhline(.5,color="gray",ls=":")
    ax.set_ylim(0,1)
    ax.legend(fontsize=8)
    save(fig,"solubility_auroc")
    fig, axes = plt.subplots(1,2,figsize=(12,4.5))
    for ax, metric in zip(axes,("auroc","average_precision")):
        sub = macro[macro.metric.eq(metric)].sort_values("layer")
        d = sub[sub.layer.gt(0)]
        ax.plot(d.layer,d.validation,color="#8054a6",label="Validation macro mean")
        ax.fill_between(d.layer,d.low,d.high,color="#8054a6",alpha=.17)
        ax.plot(d.layer,d.cv_mean,color="#24718d",label="Train CV macro mean")
        ax.axhline(sub.iloc[0].validation,color="#8054a6",ls="--",label="Composition validation")
        ax.set(xlabel="ESM-2 layer",ylabel=metric,ylim=(0,1),title=f"Pfam membership: {len(evaluable)} evaluable families")
        if metric == "auroc":
            ax.axhline(.5,color="gray",ls=":")
        else:
            support = eligible[eligible.validation_estimable]
            ax.axhline((support.validation_positive/support.validation_n).mean(),color="gray",ls=":",label="Mean validation prevalence")
        ax.legend(fontsize=8)
    save(fig,"family_macro")
    heat = family_rows[family_rows.layer.gt(0) & family_rows.metric.eq("auroc")].pivot(index="target",columns="layer",values="validation")
    heat = heat.reindex(eligible.sort_values(["train_positive","pfam"],ascending=[False,True]).pfam)
    names = catalog.set_index("pfam")["name"]
    fig, ax = plt.subplots(figsize=(13, max(5,len(heat)*.29)))
    image = ax.imshow(np.ma.masked_invalid(heat.to_numpy()),vmin=0,vmax=1,aspect="auto",cmap="viridis")
    ax.set(yticks=range(len(heat)),yticklabels=[f"{p}: {names[p]}" for p in heat.index],
           xticks=[0,7,15,23,32],xticklabels=[1,8,16,24,33],xlabel="ESM-2 layer",title="Validation family-membership AUROC (blank = undefined)")
    fig.colorbar(image,ax=ax,label="AUROC")
    save(fig,"family_heatmap")
    text = f"""# Uncontrolled property-probe results

All 33 layers evaluated on the unchanged development cohort: 2,361 training and 691 validation proteins. No reserved-test records were accessed. Continuous targets with undefined X-residue descriptors use 2,359/690 proteins; length and binary targets use all 2,361/691. Missing targets are never zero-filled.

There are {len(eligible)} training-eligible Pfam membership targets; {len(evaluable)} have both membership classes in validation. Family qualification uses training only. All proteins remain eligible examples for each membership probe. Pfam presence is annotation-defined; a negative can mean no detected hit, not verified biological absence. Overlapping assignments are retained.

R² and MAE describe continuous-target decoding; AUROC/AP describe binary targets. These scales are not comparable. Hydrophobicity and charge are linear functions of composition; log length is a baseline input. Their composition baseline should be near-perfect and is a diagnostic reference, not an independent biological discovery. AGGRESCAN targets are predictions, not aggregation measurements.

Training bands show ±1 fold SD. Validation bands are pointwise 95% homology-group bootstrap intervals conditional on fitted models. The same training folds select regularization; CV scores are tuning estimates, not nested-CV estimates. Intervals omit model-fitting and layer/target selection uncertainty. Family macro intervals retain the same families in every valid resample; valid draw counts are reported. A macro interval is suppressed if fewer than 90% of resamples support every included family, rather than presenting a heavily conditioned interval. Sparse individual-family intervals also need caution.

`metrics.tsv` contains every target/layer result; `paired_vs_composition.tsv` uses the same group resamples for paired differences (ESM minus composition; negative MAE/RMSE differences favor ESM). `family_support.tsv` reports exact denominators. `constant_baselines.tsv` gives training-mean/prevalence references. Each layer directory stores coefficients, scalers, tuning scores and ID/hash-aligned predictions. Observed per-protein targets remain in local inputs.

Review these unrestricted results before specifying controlled comparisons. Strong decoding does not establish direct encoding or a causal computation; weak decoding does not establish absence of information. No restriction, matching, family-controlled evaluation, SAE feature analysis or intervention was performed.
"""
    (root/"UNCONTROLLED_RESULTS.md").write_text(text)
    write_json(root/"analysis_integrity.json",dict(layers=33,targets=len(targets),families=len(eligible),
               families_validation_estimable=len(evaluable),same_target_population_across_layers=True,reserved_test_accessed=False))
