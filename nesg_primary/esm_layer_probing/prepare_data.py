"""Export verified development-only inputs; preserve domain multi-memberships."""

import hashlib
import json
import pandas as pd
from common import HERE, PROTOCOL, sha, write_json, read_local


def main():
    src = HERE.parent / "results"
    dest = HERE / "data"
    dest.mkdir(exist_ok=True)
    lock = json.loads((src / "test_lock.json").read_text())
    names = [
        "development_records.tsv.gz",
        "development_domains.tsv.gz",
        "training_manifest.tsv",
        "validation_manifest.tsv",
        "training_cv_manifest.tsv",
        "split_assignments.tsv",
        "initial_split_counts.tsv",
    ]
    for name in names:
        assert sha(src / name) == lock["rerun_safe_checksums"][name], name
    d = pd.read_csv(src / "development_records.tsv.gz", sep="\t")
    assert d.split.value_counts().to_dict() == PROTOCOL["counts"]
    assert d.record_id.is_unique and d.sha256.is_unique
    assert d.sequence.map(lambda s: hashlib.sha256(s.encode()).hexdigest()).equals(
        d.sha256
    )
    assert d.exp.isin([4, 5]).all() and d.sol.isin([0, 4, 5]).all()
    assert (d.label == d.sol.isin([4, 5]).astype(int)).all()
    # This file has identifiers/groups/splits only, not reserved-test sequences or labels.
    assignments = pd.read_csv(src / "split_assignments.tsv", sep="\t")
    assert assignments.split.value_counts().to_dict() == {
        "train": 2361,
        "validation": 691,
        "reserved_test": 710,
    }
    assert assignments.groupby("audit_cluster").split.nunique().max() == 1
    assert assignments.groupby("sha256").split.nunique().max() == 1
    fields = ["record_id", "sha256", "split", "audit_cluster", "inner_cv_fold"]
    pd.testing.assert_frame_equal(
        d[fields].sort_values("record_id").reset_index(drop=True),
        assignments[assignments.split.isin(["train", "validation"])][fields]
        .sort_values("record_id")
        .reset_index(drop=True),
        check_dtype=False,
    )
    for split, name in [
        ("train", "training_manifest.tsv"),
        ("validation", "validation_manifest.tsv"),
    ]:
        original = pd.read_csv(src / name, sep="\t").set_index("record_id").sort_index()
        current = d[d.split == split].set_index("record_id").sort_index()
        pd.testing.assert_frame_equal(
            current[original.columns], original, check_dtype=False
        )
    cv = (
        pd.read_csv(src / "training_cv_manifest.tsv", sep="\t")
        .sort_values("record_id")
        .reset_index(drop=True)
    )
    pd.testing.assert_frame_equal(
        d[d.split == "train"][cv.columns]
        .sort_values("record_id")
        .reset_index(drop=True),
        cv,
        check_dtype=False,
    )
    d = d.sort_values(["split", "record_id"], kind="stable").reset_index(drop=True)
    d[["record_id", "sha256", "split", "length", "sequence"]].to_csv(
        dest / "extraction_manifest.tsv", sep="\t", index=False
    )
    d[
        [
            "record_id",
            "sha256",
            "split",
            "sid",
            "label",
            "exp",
            "sol",
            "audit_cluster",
            "inner_cv_fold",
        ]
    ].to_csv(dest / "local_labels.tsv", sep="\t", index=False)
    domains = pd.read_csv(src / "development_domains.tsv.gz", sep="\t")
    assert set(domains.sha256) <= set(d.sha256)
    assert (domains.sequence_score >= domains.sequence_ga).all() and (
        domains.domain_score >= domains.domain_ga
    ).all()
    domains.merge(
        d[["record_id", "sha256", "split"]], on="sha256", validate="many_to_one"
    ).to_csv(
        dest / "pfam_domains.tsv.gz",
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    domains[
        [
            "sha256",
            "pfam",
            "name",
            "description",
            "clan",
            "clan_name",
            "competing_overlap",
        ]
    ].drop_duplicates().merge(
        d[["record_id", "sha256", "split"]], on="sha256", validate="many_to_one"
    ).to_csv(dest / "pfam_memberships.tsv", sep="\t", index=False)
    pin = {
        "original_safe_inputs": {name: sha(src / name) for name in names},
        "extraction_manifest_sha256": sha(dest / "extraction_manifest.tsv"),
        "local_labels_sha256": sha(dest / "local_labels.tsv"),
        "pfam_domains_sha256": sha(dest / "pfam_domains.tsv.gz"),
        "pfam_memberships_sha256": sha(dest / "pfam_memberships.tsv"),
    }
    write_json(dest / "input_integrity.json", pin)
    p = HERE / "protocol.json"
    if p.exists():
        assert json.loads(p.read_text()) == PROTOCOL, (
            "Protocol changed; review explicitly"
        )
    else:
        write_json(p, PROTOCOL)
    read_local()
    checks = dict(
        proteins=len(d),
        counts=d.groupby(["split", "label"]).size().unstack().to_dict("index"),
        groups=d.groupby("split").audit_cluster.nunique().to_dict(),
        min_length=int(d.length.min()),
        max_length=int(d.length.max()),
        X_proteins=int(d.sequence.str.contains("X").sum()),
        reserved_count_verified_from_label_free_assignments=710,
        reserved_sequences_or_labels_opened=False,
        all_input_hashes_verified=True,
        sequence_hashes_verified=True,
        split_and_fold_integrity=True,
        GA_hits_verified=True,
    )
    write_json(dest / "checks.json", checks)
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
