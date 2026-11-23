"""Uncontrolled linear property probes; target masks never change the cohort."""
import json
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.linalg import svd
from scipy.special import expit
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.exceptions import ConvergenceWarning
from common import sha, write_json
from property_inputs import SPEC, CONTINUOUS
from probe import bootstrap_weights, weighted_metrics


def regression_metrics(y, p, weights):
    """Vectorized weighted metrics; constant-target R² is undefined, not forced to 0/1."""
    w = np.asarray(weights, dtype=float)
    n = w.sum(axis=1)
    safe = np.where(n > 0, n, np.nan)
    mean = (w @ y) / safe
    sst = (w @ (y*y)) - safe*mean*mean
    sse = w @ ((y-p)**2)
    r2 = np.full(len(w), np.nan)
    valid = (n > 0) & (sst > np.finfo(float).eps * np.maximum(1., w @ (y*y))*16)
    r2[valid] = 1 - sse[valid]/sst[valid]
    return dict(r2=r2, mae=(w @ np.abs(y-p))/safe, rmse=np.sqrt(sse/safe))


def binary_metrics(y, score, weights):
    a, p = weighted_metrics(y, score, weights)
    return dict(auroc=a, average_precision=p)


def scalar_metrics(y, p, binary):
    values = (binary_metrics if binary else regression_metrics)(y, p, np.ones((1, len(y))))
    return {k: float(v[0]) for k, v in values.items()}


class RidgeBasis:
    """One training-only standardized SVD reused across alpha candidates/targets."""
    def __init__(self, X):
        self.scaler = StandardScaler().fit(X)
        z = self.scaler.transform(X)
        self.u, self.s, self.vt = svd(z, full_matrices=False, check_finite=False)

    def parameters(self, y, alpha):
        intercept = float(np.mean(y))
        coef = self.vt.T @ ((self.s/(self.s*self.s + alpha)) * (self.u.T @ (y-intercept)))
        return coef, intercept


def target_list(data, properties, memberships, catalog):
    result = [("solubility", True, data.label.to_numpy(dtype=float))]
    result += [(name, False, properties[name].to_numpy(dtype=float)) for name in CONTINUOUS]
    for family in catalog.loc[catalog.eligible, "pfam"]:
        ids = memberships.loc[memberships.pfam.eq(family), "record_id"]
        result.append((family, True, data.record_id.isin(ids).to_numpy(dtype=float)))
    return result


def valid_checkpoint(folder):
    marker = folder / "target_complete.json"
    if not marker.is_file():
        return False
    saved = json.loads(marker.read_text())
    for name, digest in saved["files"].items():
        assert Path(name).name == name
        assert sha(folder/name) == digest, f"Damaged property checkpoint: {name}"
    return True


def finish(folder):
    write_json(folder/"target_complete.json", {"files": {
        p.name: sha(p) for p in sorted(folder.iterdir()) if p.is_file() and p.name != "target_complete.json"
    }})


def save_results(folder, data, y, tr, va, oof, fitted, pred, cv, chosen, binary, weights, source):
    """Published per-protein tables contain predictions/IDs, never observed targets."""
    cv = pd.DataFrame(cv)
    cv["selected"] = cv.parameter.eq(chosen)
    cv.to_csv(folder/"cv.tsv", sep="\t", index=False)
    selected = cv[cv.selected]
    boots = (binary_metrics if binary else regression_metrics)(y[va], pred, weights)
    point = scalar_metrics(y[va], pred, binary)
    trainpoint = scalar_metrics(y[tr], fitted, binary)
    rows = []
    for name, samples in boots.items():
        good = samples[np.isfinite(samples)]
        lo, hi = np.quantile(good, [.025, .975]) if len(good) else (np.nan, np.nan)
        rows.append(dict(metric=name, validation=point[name], low=lo, high=hi,
                         cv_mean=selected[name].mean(), cv_sd=selected[name].std(),
                         training_fit=trainpoint[name], parameter=chosen,
                         train_n=int(tr.sum()), validation_n=int(va.sum()),
                         train_positive=int(y[tr].sum()) if binary else np.nan,
                         validation_positive=int(y[va].sum()) if binary else np.nan,
                         bootstrap_valid=len(good), source=source))
    pd.DataFrame(rows).to_csv(folder/"metrics.tsv", sep="\t", index=False)
    np.savez_compressed(folder/"bootstrap.npz", **boots)
    for label, mask, values in [("validation", va, pred), ("training_oof", tr, oof), ("training_fit", tr, fitted)]:
        table = data.loc[mask, ["record_id", "sha256", "split"]].copy()
        table["prediction"] = values
        if binary:
            table["probability"] = expit(values)
        table.to_csv(folder/f"{label}_predictions.tsv", sep="\t", index=False)
    # Constant reference is fitted using observed training outcomes only.
    constant = np.full(va.sum(), y[tr].mean())
    constants = scalar_metrics(y[va], constant, binary)
    pd.DataFrame([dict(metric=k, validation=v, training_constant=y[tr].mean())
                  for k, v in constants.items()]).to_csv(folder/"constant.tsv", sep="\t", index=False)
    finish(folder)


def fit_target(X, data, target, binary, y, folder, bases, weights):
    folder.mkdir(parents=True, exist_ok=True)
    if valid_checkpoint(folder):
        print(f"  Reusing {target}", flush=True)
        return
    observed = np.isfinite(y)
    tr = data.split.eq("train").to_numpy() & observed
    va = data.split.eq("validation").to_numpy() & observed
    folds = data.inner_cv_fold.to_numpy()
    w = weights[:, observed[data.split.eq("validation").to_numpy()]]
    assert np.isfinite(X).all() and np.isfinite(y[observed]).all()
    grid = SPEC["classification_C_grid" if binary else "regression_alpha_grid"]
    cv, candidates = [], {}
    for parameter in grid:
        oof = np.full(tr.sum(), np.nan)
        for f in range(5):
            fit, held = tr & (folds != f), tr & (folds == f)
            assert set(data.loc[fit, "audit_cluster"]).isdisjoint(data.loc[held, "audit_cluster"])
            coef, intercept, scaler = fit_parameters(X, y, fit, parameter, binary, bases)
            p = scaler.transform(X[held]) @ coef + intercept
            oof[folds[tr] == f] = p
            metrics = scalar_metrics(y[held], p, binary)
            cv.append(dict(parameter=parameter, fold=f, n=int(held.sum()),
                           mse=float(np.mean((p-y[held])**2)) if not binary else np.nan, **metrics))
        assert np.isfinite(oof).all()
        candidates[parameter] = oof
    scores = pd.DataFrame(cv).groupby("parameter")["auroc" if binary else "mse"].mean()
    chosen = float(min(grid, key=lambda p: (-scores[p], p) if binary else (scores[p], -p)))
    coef, intercept, scaler = fit_parameters(X, y, tr, chosen, binary, bases)
    np.savez(folder/"parameters.npz", mean=scaler.mean_, scale=scaler.scale_, coef=coef,
             intercept=intercept, parameter=chosen, binary=binary)
    save_results(folder, data, y, tr, va, candidates[chosen],
                 scaler.transform(X[tr]) @ coef+intercept, scaler.transform(X[va]) @ coef+intercept,
                 cv, chosen, binary, w, "new_training_only_fit")
    print(f"  Completed {target}", flush=True)


def fit_parameters(X, y, fit, parameter, binary, bases):
    if binary:
        assert np.unique(y[fit]).size == 2
        scaler = StandardScaler().fit(X[fit])
        model = LogisticRegression(C=parameter, max_iter=SPEC["max_iter"],
                                   solver="lbfgs", random_state=SPEC["seed"])
        with warnings.catch_warnings():
            warnings.simplefilter("error", ConvergenceWarning)
            model.fit(scaler.transform(X[fit]), y[fit])
        return model.coef_[0], float(model.intercept_[0]), scaler
    key = np.packbits(fit).tobytes()
    if key not in bases:
        bases[key] = RidgeBasis(X[fit])
    basis = bases[key]
    coef, intercept = basis.parameters(y[fit], parameter)
    return coef, intercept, basis.scaler


def reuse_solubility(old, model, data, folder, weights):
    """Normalize a previously verified solubility fit; never refit it."""
    import shutil
    folder.mkdir(parents=True, exist_ok=True)
    if valid_checkpoint(folder):
        return
    tr, va = data.split.eq("train").to_numpy(), data.split.eq("validation").to_numpy()
    arrays = []
    for suffix, mask in [("training_oof", tr), ("training_fit", tr), ("validation_predictions", va)]:
        table = pd.read_csv(old/f"{model}_{suffix}.tsv", sep="\t", float_precision="round_trip",
                            dtype={"record_id":str,"sha256":str,"split":str})
        assert table.record_id.is_unique and set(table.record_id) == set(data.loc[mask, "record_id"])
        table = table.set_index("record_id").loc[data.loc[mask, "record_id"]]
        assert table.sha256.tolist() == data.loc[mask, "sha256"].tolist()
        arrays.append(table.decision_score.to_numpy())
    cv = pd.read_csv(old/f"{model}_cv.tsv", sep="\t").rename(columns={"C": "parameter"})
    chosen = float(json.loads((old/f"{model}_metrics.json").read_text())["C"])
    shutil.copyfile(old/"models"/f"{model}_parameters.npz", folder/"parameters.npz")
    save_results(folder, data, data.label.to_numpy(), tr, va, *arrays,
                 cv.to_dict("records"), chosen, True, weights, "reused_existing_solubility_probe")


def run_layer(X, data, targets, folder, old=None, model=None):
    X = np.asarray(X, dtype=np.float64)
    bases = {}
    weights = bootstrap_weights(data.loc[data.split.eq("validation"), "audit_cluster"],
                                SPEC["bootstrap_repetitions"], SPEC["seed"])
    for target, binary, y in targets:
        if target == "solubility" and old is not None:
            reuse_solubility(old, model, data, folder/target, weights)
        else:
            fit_target(X, data, target, binary, y, folder/target, bases, weights)
