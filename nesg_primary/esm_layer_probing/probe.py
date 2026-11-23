"""Identical frozen-fold linear-probe protocol for composition and ESM layers."""

import argparse
import importlib.metadata
import json
import warnings
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
from scipy.special import expit
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits
from common import HERE, PROTOCOL, composition, read_local, sha, write_json


def bootstrap_weights(groups, repetitions=2000, seed=20261118):
    unique, inverse = np.unique(groups, return_inverse=True)
    rng = np.random.default_rng(seed)
    multiplicities = rng.multinomial(
        len(unique), np.full(len(unique), 1 / len(unique)), size=repetitions
    )
    return multiplicities[:, inverse]


def weighted_metrics(y, scores, weights):
    """Vectorized group bootstrap AUROC/AP, including exact score ties."""
    order = np.argsort(scores, kind="stable")
    s = np.asarray(scores)[order]
    starts = np.r_[0, np.flatnonzero(np.diff(s)) + 1]
    positive = np.add.reduceat(weights[:, order] * np.asarray(y)[order], starts, axis=1)
    negative = np.add.reduceat(
        weights[:, order] * (1 - np.asarray(y)[order]), starts, axis=1
    )
    P = positive.sum(axis=1)
    N = negative.sum(axis=1)
    valid = (P > 0) & (N > 0)
    auc = np.full(len(weights), np.nan)
    ap = auc.copy()
    auc[valid] = (positive * (np.cumsum(negative, axis=1) - 0.5 * negative)).sum(
        axis=1
    )[valid] / (P[valid] * N[valid])
    p = positive[:, ::-1]
    n = negative[:, ::-1]
    cp = p.cumsum(axis=1)
    den = cp + n.cumsum(axis=1)
    precision = np.divide(cp, den, out=np.zeros_like(cp, dtype=float), where=den > 0)
    ap[valid] = (p * precision).sum(axis=1)[valid] / P[valid]
    return auc, ap


def fit_probe(X, data, name, destination, weights):
    tr = (data.split == "train").to_numpy()
    va = ~tr
    y = data.label.to_numpy()
    folds = data.inner_cv_fold.to_numpy()
    assert np.isfinite(X).all() and X.shape[0] == len(data)
    cv = []
    candidates = {}
    # Scaling is learned anew on the fitting fold for every candidate.
    for C in PROTOCOL["C_grid"]:
        oof = np.full(tr.sum(), np.nan)
        for fold in range(5):
            fit = tr & (folds != fold)
            held = tr & (folds == fold)
            assert not set(data.loc[fit, "audit_cluster"]) & set(
                data.loc[held, "audit_cluster"]
            )
            pipeline = make_pipeline(
                StandardScaler(),
                LogisticRegression(
                    C=C,
                    penalty="l2",
                    solver="lbfgs",
                    max_iter=PROTOCOL["max_iter"],
                    random_state=PROTOCOL["seed"],
                ),
            )
            with warnings.catch_warnings():
                warnings.simplefilter("error", ConvergenceWarning)
                pipeline.fit(X[fit], y[fit])
            score = pipeline.decision_function(X[held])
            oof[folds[tr] == fold] = score
            cv.append(
                dict(
                    model=name,
                    C=C,
                    fold=fold,
                    n=int(held.sum()),
                    auroc=roc_auc_score(y[held], score),
                    average_precision=average_precision_score(y[held], score),
                )
            )
        candidates[C] = oof
    table = pd.DataFrame(cv)
    ranked = (
        table.groupby("C")
        .auroc.mean()
        .sort_index()
        .sort_values(ascending=False, kind="stable")
    )
    C = float(ranked.index[0])
    oof = candidates[C]
    pipeline = make_pipeline(
        StandardScaler(),
        LogisticRegression(
            C=C,
            penalty="l2",
            solver="lbfgs",
            max_iter=PROTOCOL["max_iter"],
            random_state=PROTOCOL["seed"],
        ),
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        pipeline.fit(X[tr], y[tr])
    scores = pipeline.decision_function(X[va])
    training_fit_scores = pipeline.decision_function(X[tr])
    probabilities = pipeline.predict_proba(X[va])[:, 1]
    au, ap = weighted_metrics(y[va], scores, weights)
    model_dir = destination / "models"
    model_dir.mkdir(exist_ok=True)
    joblib.dump(pipeline, model_dir / f"{name}.joblib")
    scaler = pipeline.named_steps["standardscaler"]
    lr = pipeline.named_steps["logisticregression"]
    np.savez(
        model_dir / f"{name}_parameters.npz",
        mean=scaler.mean_,
        scale=scaler.scale_,
        coef=lr.coef_,
        intercept=lr.intercept_,
        classes=lr.classes_,
    )
    table["selected_C"] = C
    table.to_csv(destination / f"{name}_cv.tsv", sep="\t", index=False)
    cols = ["record_id", "sha256", "split", "label", "audit_cluster", "inner_cv_fold"]
    pred = data.loc[va, cols].copy()
    pred["model"] = name
    pred["layer"] = int(name.split("_")[1]) if name.startswith("layer_") else 0
    pred["decision_score"] = scores
    pred["probability"] = probabilities
    pred.to_csv(
        destination / f"{name}_validation_predictions.tsv", sep="\t", index=False
    )
    trainpred = data.loc[tr, cols].copy()
    trainpred["decision_score"] = oof
    trainpred["probability"] = expit(oof)
    trainpred.to_csv(destination / f"{name}_training_oof.tsv", sep="\t", index=False)
    fitted_train = data.loc[tr, cols].copy()
    fitted_train["decision_score"] = training_fit_scores
    fitted_train["probability"] = expit(training_fit_scores)
    fitted_train.to_csv(destination / f"{name}_training_fit.tsv", sep="\t", index=False)
    np.savez(destination / f"{name}_bootstrap.npz", auroc=au, average_precision=ap)
    selected = table[table.C == C]
    metric = dict(
        model=name,
        layer=int(name.split("_")[1]) if name.startswith("layer_") else 0,
        C=C,
        training_n=int(tr.sum()),
        training_fit_auroc=roc_auc_score(y[tr], training_fit_scores),
        training_oof_auroc=roc_auc_score(y[tr], oof),
        validation_n=int(va.sum()),
        auroc=roc_auc_score(y[va], scores),
        average_precision=average_precision_score(y[va], scores),
        auroc_low=float(np.nanquantile(au, 0.025)),
        auroc_high=float(np.nanquantile(au, 0.975)),
        ap_low=float(np.nanquantile(ap, 0.025)),
        ap_high=float(np.nanquantile(ap, 0.975)),
        cv_mean=float(selected.auroc.mean()),
        cv_sd=float(selected.auroc.std()),
        cv_min=float(selected.auroc.min()),
        cv_max=float(selected.auroc.max()),
        bootstrap_valid=int(np.isfinite(au).sum()),
    )
    write_json(destination / f"{name}_metrics.json", metric)
    return metric


def fit_or_reuse(X, data, name, destination, weights, resume=False):
    """A completed model is reusable only when every saved artifact still hashes."""
    marker = destination / f"{name}_complete.json"
    if resume and marker.is_file():
        saved = json.loads(marker.read_text())
        assert saved["protocol_sha256"] == sha(destination / "probe_protocol.json")
        for filename, digest in saved["files"].items():
            assert ".." not in Path(filename).parts and not Path(filename).is_absolute()
            assert sha(destination / filename) == digest, (
                f"Corrupt probe output: {filename}"
            )
        print(f"Reusing verified {name} probe; no fitting.", flush=True)
        return json.loads((destination / f"{name}_metrics.json").read_text())
    result = fit_probe(X, data, name, destination, weights)
    if resume:
        paths = (
            list(destination.glob(f"{name}_*"))
            + list((destination / "models").glob(f"{name}.*"))
            + list((destination / "models").glob(f"{name}_*"))
        )
        write_json(
            marker,
            dict(
                protocol_sha256=sha(destination / "probe_protocol.json"),
                files={
                    str(p.relative_to(destination)): sha(p)
                    for p in paths
                    if p != marker
                },
            ),
        )
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--activations", type=Path)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--output", type=Path, default=HERE / "results")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    data = read_local()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    if args.resume:
        identity = dict(
            protocol_sha256=sha(HERE / "protocol.json"),
            script_sha256=sha(__file__),
            input_sha256=sha(HERE / "data/extraction_manifest.tsv"),
            labels_sha256=sha(HERE / "data/local_labels.tsv"),
            activation_completion_sha256=sha(args.activations / "completed.json")
            if args.activations
            else None,
            software={
                p: importlib.metadata.version(p)
                for p in ["numpy", "scipy", "scikit-learn", "joblib"]
            },
        )
        path = out / "probe_protocol.json"
        if path.exists():
            assert json.loads(path.read_text()) == identity, (
                "Probe settings/environment changed; use another output folder."
            )
        else:
            assert not list(out.glob("*_metrics.json")), (
                "Use a new folder for resumable probes."
            )
            write_json(path, identity)
    weights = bootstrap_weights(
        data.loc[data.split == "validation", "audit_cluster"].to_numpy(),
        PROTOCOL["bootstrap_replicates"],
        PROTOCOL["seed"],
    )
    metrics = []
    with threadpool_limits(limits=args.threads):
        print("Fitting composition baseline with frozen training CV", flush=True)
        metrics.append(
            fit_or_reuse(
                composition(data.sequence),
                data,
                "composition",
                out,
                weights,
                args.resume,
            )
        )
        if args.activations:
            from extract import load_all

            X = load_all(
                args.activations,
                data[["record_id", "sha256", "split", "length", "sequence"]],
            )
            for layer in range(1, 34):
                print(f"Fitting layer {layer}/33", flush=True)
                metrics.append(
                    fit_or_reuse(
                        X[:, layer - 1, :],
                        data,
                        f"layer_{layer:02}",
                        out,
                        weights,
                        args.resume,
                    )
                )
    table = pd.DataFrame(metrics)
    base = np.load(out / "composition_bootstrap.npz")["auroc"]
    for row in metrics:
        draws = np.load(out / f"{row['model']}_bootstrap.npz")["auroc"] - base
        row["delta_vs_composition"] = row["auroc"] - metrics[0]["auroc"]
        row["delta_low"] = float(np.nanquantile(draws, 0.025))
        row["delta_high"] = float(np.nanquantile(draws, 0.975))
    table = pd.DataFrame(metrics)
    table.to_csv(out / "metrics.tsv", sep="\t", index=False)
    pd.concat(
        [pd.read_csv(out / f"{m['model']}_cv.tsv", sep="\t") for m in metrics]
    ).to_csv(out / "all_cv.tsv", sep="\t", index=False)
    pair = []
    for a, b in [(16, 8), (32, 8), (32, 16)]:
        if len(metrics) == 34:
            diff = (
                np.load(out / f"layer_{a:02}_bootstrap.npz")["auroc"]
                - np.load(out / f"layer_{b:02}_bootstrap.npz")["auroc"]
            )
            pair.append(
                dict(
                    layer_a=a,
                    layer_b=b,
                    delta=float(
                        table.set_index("layer").loc[a, "auroc"]
                        - table.set_index("layer").loc[b, "auroc"]
                    ),
                    low=float(np.nanquantile(diff, 0.025)),
                    high=float(np.nanquantile(diff, 0.975)),
                )
            )
    pd.DataFrame(pair, columns=["layer_a", "layer_b", "delta", "low", "high"]).to_csv(
        out / "representative_paired_differences.tsv", sep="\t", index=False
    )
    write_json(
        out / "run_status.json",
        dict(
            ESM_layers_completed=len(metrics) - 1,
            baseline_completed=True,
            test_evaluated=False,
            extraction_directory=str(args.activations) if args.activations else None,
            protocol_sha256=sha(HERE / "protocol.json"),
            input_manifest_sha256=sha(HERE / "data/extraction_manifest.tsv"),
            local_labels_sha256=sha(HERE / "data/local_labels.tsv"),
            CV_interpretation="Selected-C fold estimates reuse hyperparameter selection folds; not unbiased nested-CV estimates",
            uncertainty="2000 paired homology-group bootstrap draws conditional on fitted models; no training-resampling or layer-selection correction",
        ),
    )
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
