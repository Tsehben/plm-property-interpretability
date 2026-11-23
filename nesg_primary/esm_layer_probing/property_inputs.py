"""Build development-only property inputs from completed, hash-verified annotations.

Deliberately never opens combined AGGRESCAN tables: they contain reserved-test
records. Read only train_*/validation_* jobs listed in the development export.
"""

import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from common import HERE, AA, PROTOCOL, read_local, sha, write_json

CONTINUOUS = ["hydrophobicity", "charge_proxy", "log_length", "aggregation_Na4vSS", "apr_fraction"]
SPEC = dict(
    version=1, scope="uncontrolled development property probes only",
    layers=list(range(1, 34)), regression_alpha_grid=[.0001, .001, .01, .1, 1., 10., 100., 1000., 10000.],
    regression_selection="minimum unweighted mean held-out-fold MSE; ties prefer larger alpha",
    classification_C_grid=PROTOCOL["C_grid"], classification_selection="maximum mean held-out-fold AUROC; ties prefer smaller C",
    family_min_train_groups_per_class=5, family_each_cv_partition_requires_both_classes=True,
    family_target="presence of >=1 GA-passing Pfam 37.0 assignment; multiple memberships retained",
    family_macro="fixed set of training-eligible families with both classes in validation; undefined targets still reported",
    family_macro_ci_min_valid_fraction=0.9,
    cv="five frozen training-only homology-group folds", class_weight=None,
    bootstrap_repetitions=2000, seed=20261118, max_iter=3000,
    continuous_targets=CONTINUOUS, test_accessed=False,
    target_missingness="target-specific observed records only; no target imputation or primary-cohort filtering",
)


def family_catalog(data, memberships):
    tr = data.split.eq("train")
    rows = []
    for accession, group in memberships.groupby("pfam", sort=True):
        y = data.record_id.isin(group.record_id)
        pos = data.loc[tr & y, "audit_cluster"].nunique()
        neg = data.loc[tr & ~y, "audit_cluster"].nunique()
        estimable = all(y[tr & data.inner_cv_fold.eq(f)].nunique() == 2 and
                        y[tr & data.inner_cv_fold.ne(f)].nunique() == 2 for f in range(5))
        eligible = pos >= 5 and neg >= 5 and estimable
        first = group.iloc[0]
        rows.append(dict(pfam=accession, name=first["name"], description=first.description,
                         train_positive=int((tr & y).sum()), train_negative=int((tr & ~y).sum()),
                         train_positive_groups=int(pos), train_negative_groups=int(neg),
                         all_cv_partitions_both_classes=estimable, eligible=eligible,
                         reason="eligible" if eligible else "insufficient_training_group_support" if min(pos, neg) < 5 else "single_class_CV_partition"))
    return pd.DataFrame(rows)


def prepare():
    data = read_local()
    original_pin = json.loads((HERE / "data/input_integrity.json").read_text())
    dev = HERE.parent / "results/development_records.tsv.gz"
    assert sha(dev) == original_pin["original_safe_inputs"][dev.name]
    records = pd.read_csv(dev, sep="\t").set_index("record_id").loc[data.record_id]
    assert records.sha256.tolist() == data.sha256.tolist()
    memberships_path = HERE / "data/pfam_memberships.tsv"
    assert sha(memberships_path) == original_pin["pfam_memberships_sha256"]
    memberships = pd.read_csv(memberships_path, sep="\t", keep_default_na=False)
    assert set(memberships.record_id) <= set(data.record_id)
    mapped = memberships.merge(data[["record_id", "sha256", "split"]], on="record_id", suffixes=("", "_frozen"), validate="many_to_one")
    assert mapped.sha256.eq(mapped.sha256_frozen).all() and mapped.split.eq(mapped.split_frozen).all()
    exports = HERE.parent / "server_uploads/aggrescan"
    manifest_path = exports / "upload_manifest.tsv"
    provenance = json.loads((exports / "export_provenance.json").read_text())
    assert sha(manifest_path) == provenance["outputs_sha256"]["upload_manifest.tsv"]
    manifest = pd.read_csv(manifest_path, sep="\t", keep_default_na=False)
    assert set(manifest.split) == {"train", "validation"} and len(manifest) == 3052
    manifest = manifest.set_index("record_id").loc[data.record_id]
    assert manifest.sha256.tolist() == data.sha256.tolist()
    annotations = []
    inputs = {str(dev.relative_to(HERE.parent)): sha(dev),
              str(manifest_path.relative_to(HERE.parent)): sha(manifest_path),
              "esm_layer_probing/data/pfam_memberships.tsv": sha(memberships_path)}
    for batch in sorted(set(manifest.batch_file) - {""}):
        name = Path(batch).stem
        assert name.startswith(("train_", "validation_")) and Path(batch).parent == Path("batches")
        folder = HERE.parent / "aggregation_server/jobs" / name
        state_path, result_path = folder / "job.json", folder / "protein_results.json"
        state = json.loads(state_path.read_text())
        assert state["status"] == "complete" and state["split"] in {"train", "validation"}
        assert state["sequence_hashes_verified"] and state["residue_coordinates_verified"]
        assert sha(result_path) == state["outputs_sha256"]["protein_results.json"]
        rows = json.loads(result_path.read_text())
        assert {r["record_id"] for r in rows} == set(manifest.index[manifest.batch_file.eq(batch)])
        annotations.extend(rows)
        for path in (state_path, result_path):
            inputs[str(path.relative_to(HERE.parent))] = sha(path)
    agg = pd.DataFrame(annotations).set_index("record_id")
    assert agg.index.is_unique and len(agg) == 3049
    assert agg.sha256.equals(data.set_index("record_id").loc[agg.index, "sha256"])
    out = data[["record_id", "sha256", "split"]].copy()
    out["hydrophobicity"] = records.hydropathy.to_numpy()
    out["charge_proxy"] = records.charge_proxy.to_numpy()
    out["log_length"] = np.log(data.length)
    for target, source in [("aggregation_Na4vSS", "Na4vSS"), ("apr_fraction", "apr_fraction")]:
        out[target] = data.record_id.map(agg[source])
    out["aggregation_status"] = data.record_id.map(agg.status).fillna("unsupported_input")
    standard = data.sequence.map(lambda s: set(s) <= set(AA))
    assert standard.sum() == 3049
    for target in CONTINUOUS:
        assert out[target].notna().equals(pd.Series(True, index=data.index) if target == "log_length" else standard)
    assert out.apr_fraction.dropna().between(0, 1).all()
    # Recheck the documented charge proxy; no pH-dependent charge is implied.
    expected = data.loc[standard, "sequence"].map(lambda s: (s.count("K") + s.count("R") - s.count("D") - s.count("E")) / len(s))
    assert np.allclose(out.loc[standard, "charge_proxy"], expected)
    catalog = family_catalog(data, memberships.drop_duplicates(["record_id", "pfam"]))
    dest = HERE / "data/uncontrolled"
    dest.mkdir(parents=True, exist_ok=True)
    out.to_csv(dest / "properties.tsv", sep="\t", index=False)
    catalog.to_csv(dest / "family_catalog.tsv", sep="\t", index=False)
    write_json(dest / "protocol.json", SPEC)
    coverage = []
    for split in ("train", "validation"):
        for target in CONTINUOUS:
            mask = data.split.eq(split)
            observed = out[target].notna() & mask
            coverage.append(dict(split=split, target=target, total=int(mask.sum()), observed=int(observed.sum()),
                                 missing=int((mask & ~observed).sum())))
    pd.DataFrame(coverage).to_csv(dest / "coverage.tsv", sep="\t", index=False)
    write_json(dest / "integrity.json", dict(
        source_sha256=inputs, extraction_manifest_sha256=sha(HERE / "data/extraction_manifest.tsv"),
        local_labels_sha256=sha(HERE / "data/local_labels.tsv"),
        files={p.name: sha(p) for p in sorted(dest.iterdir()) if p.name != "integrity.json"},
        records=3052, eligible_families=int(catalog.eligible.sum()),
        original_sequences_splits_labels_groups_unchanged=True, reserved_test_opened=False,
    ))
    print(pd.DataFrame(coverage).to_string(index=False))
    print("Training-eligible Pfam targets:", int(catalog.eligible.sum()))


def load_inputs():
    data = read_local()
    folder = HERE / "data/uncontrolled"
    pin = json.loads((folder / "integrity.json").read_text())
    assert pin["extraction_manifest_sha256"] == sha(HERE / "data/extraction_manifest.tsv")
    assert pin["local_labels_sha256"] == sha(HERE / "data/local_labels.tsv")
    for name, digest in pin["files"].items():
        assert Path(name).name == name and sha(folder / name) == digest
    properties = pd.read_csv(folder / "properties.tsv", sep="\t")
    assert properties[["record_id", "sha256", "split"]].equals(data[["record_id", "sha256", "split"]])
    original_pin = json.loads((HERE / "data/input_integrity.json").read_text())
    assert sha(HERE / "data/pfam_memberships.tsv") == original_pin["pfam_memberships_sha256"]
    memberships = pd.read_csv(HERE / "data/pfam_memberships.tsv", sep="\t").drop_duplicates(["record_id", "pfam"])
    catalog = pd.read_csv(folder / "family_catalog.tsv", sep="\t")
    assert json.loads((folder / "protocol.json").read_text()) == SPEC
    return data, properties, memberships, catalog


if __name__ == "__main__":
    argparse.ArgumentParser(description=__doc__).parse_args()
    prepare()
