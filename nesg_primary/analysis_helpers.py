"""Cohort validation, family counts and guarded data access for the NESG notebook."""

from pathlib import Path
import hashlib
import json
import sys
import zipfile
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display, Markdown

AA = "ACDEFGHIKLMNPQRSTVWY"
COLORS = {0: "#9aabba", 1: "#8168aa"}
LABELS = {0: "Insoluble (S=0)", 1: "Soluble (S=4–5)"}
ORDER = ["train", "validation", "reserved_test"]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def seqhash(s):
    return hashlib.sha256(s.encode()).hexdigest()


def read(path):
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def save(df, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(
        path,
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None,
    )


def md(s):
    display(Markdown(s))


class ReadGuard:
    """Kernel audit hook prevents reopening test/upstream project data after freeze."""

    def __init__(self, repo, root, locked):
        self.repo = Path(repo).resolve()
        self.root = Path(root).resolve()
        self.locked = locked
        self.accesses = []
        sys.addaudithook(self.check)

    def check(self, event, args):
        if event != "open" or not isinstance(args[0], (str, bytes)):
            return
        p = Path(args[0].decode() if isinstance(args[0], bytes) else args[0]).absolute()
        if not p.is_relative_to(self.repo):
            return
        rel = str(p.relative_to(self.repo))
        self.accesses.append({"path": rel, "mode": str(args[1])})
        restricted = p.is_relative_to(
            self.root / "results/reserved_test"
        ) or p.is_relative_to(self.root / "results/restricted_prelock")
        if self.locked and (restricted or not p.is_relative_to(self.root)):
            raise PermissionError(
                "Prospective test lock forbids reopening project input: " + rel
            )


def show(fig, name, folder):
    fig.savefig(Path(folder) / (name + ".png"), dpi=170, bbox_inches="tight")
    display(fig)
    plt.close(fig)


def prevalence(q, column):
    return (
        q.groupby(column, dropna=False)
        .agg(N=("label", "size"), soluble=("label", "sum"))
        .assign(
            insoluble=lambda d: d.N - d.soluble,
            percent_soluble=lambda d: 100 * d.soluble / d.N,
        )
        .reset_index()
    )


def support(records, members, split_col="split"):
    x = records[["record_id", "sha256", "label", "audit_cluster", split_col]].merge(
        members.drop_duplicates(["sha256", "pfam"]),
        on="sha256",
        validate="many_to_many",
    )
    f = x.groupby("pfam").agg(
        name=("family_name", "first"),
        description=("description", "first"),
        clan=("clan", "first"),
        N=("record_id", "nunique"),
    )

    def counts(df, field, prefix):
        t = (
            df.groupby(["pfam", "label"])[field]
            .nunique()
            .unstack(fill_value=0)
            .reindex(index=f.index, columns=[0, 1], fill_value=0)
            .fillna(0)
        )
        for lab in [0, 1]:
            f[prefix + str(lab)] = t[lab].astype(int)

    counts(x, "record_id", "n")
    counts(x, "audit_cluster", "g")
    both = (
        x.groupby(["pfam", "audit_cluster"])
        .label.nunique()
        .rename("labels_in_group")
        .reset_index()
    )
    a = x.merge(both, on=["pfam", "audit_cluster"])
    counts(a[a.labels_in_group == 1], "audit_cluster", "exclusive_g")
    for split in records[split_col].unique():
        z = x[x[split_col] == split]
        counts(z, "record_id", split + "_n")
        counts(z, "audit_cluster", split + "_g")
    return f.reset_index(), x


def component_ids(ids, edges):
    parent = {s: s for s in ids}

    def root(s):
        while parent[s] != s:
            parent[s] = parent[parent[s]]
            s = parent[s]
        return s

    for a, b in edges:
        a, b = root(a), root(b)
        if a != b:
            parent[max(a, b)] = min(a, b)
    return {s: root(s) for s in parent}


def verify_frozen_snapshot(results):
    """Verify public aggregate snapshots and development files, never open reserved-test rows."""
    lock = json.loads((results / "test_lock.json").read_text())
    for rel, h in lock["rerun_safe_checksums"].items():
        assert digest(results / rel) == h, rel
    return lock


def initialize_inputs(repo, results):
    """Only called before lock creation. Narrow explicit artifact allowlist."""
    assert not (results / "test_lock.json").exists(), (
        "Initial full-cohort audit is closed."
    )
    base = repo / "dataset_comparison"
    p = base / "results"
    files = [
        "records.tsv.gz",
        "policy_groups.tsv.gz",
        "domains_with_clans.tsv.gz",
        "sequence_features.tsv.gz",
        "sequence_groups.tsv.gz",
        "homology_edges.tsv.gz",
        "pfam_run.json",
        "homology_run.json",
    ]
    pins = json.loads((p / "output_checksums.json").read_text())
    checks = {}
    for name in files:
        path = p / name
        h = digest(path)
        assert h == pins["results/" + name], f"Prior artifact changed: {name}"
        checks[str(path.relative_to(repo))] = h
    allr = read(p / "records.tsv.gz")
    raw = allr[allr.dataset == "NESG_raw"].copy()
    curated = allr[allr.dataset == "NESG_NetSolP"].copy()
    del allr
    for x in [raw, curated]:
        for col in ["exp", "sol"]:
            x[col] = pd.to_numeric(x[col])
    assert len(raw) == 9703 and raw.record_id.is_unique and raw.sid.is_unique
    assert (raw.sequence.map(seqhash) == raw.sha256).all()
    archive = base / "raw/nesg.zip"
    pin = next(
        x
        for x in json.loads((base / "source_manifest.json").read_text())
        if x["local_path"] == "raw/nesg.zip"
    )
    assert digest(archive) == pin["sha256"]
    checks[str(archive.relative_to(repo))] = pin["sha256"]
    with zipfile.ZipFile(archive) as z:
        scores = pd.read_csv(z.open("nesg.csv")).set_index("id").sort_index()
        seqs = {
            block.splitlines()[0]: "".join(block.splitlines()[1:])
            for block in z.read("nesg.fasta").decode().split(">")
            if block.strip()
        }
    original = raw.set_index("sid").sort_index()
    assert original.index.equals(scores.index)
    assert (
        original[["exp", "sol"]].to_numpy() == scores[["exp", "sol"]].to_numpy()
    ).all()
    assert original.sequence.to_dict() == seqs
    repeats = raw.sha256.duplicated(False)
    scorecounts = raw.groupby("sha256")[["exp", "sol"]].nunique()
    conflicts = set(scorecounts.index[scorecounts.max(axis=1) > 1])
    raw["repeated_sequence"] = repeats
    raw["conflicting_scores"] = raw.sha256.isin(conflicts)
    raw["exclusion_reason"] = "eligible"
    raw.loc[~raw.sol.isin([0, 4, 5]), "exclusion_reason"] = "intermediate_S"
    raw.loc[~raw.exp.isin([4, 5]), "exclusion_reason"] = "expression_not_4_or_5"
    raw.loc[raw.exp == 0, "exclusion_reason"] = "E_zero_not_interpretable_as_solubility"
    raw.loc[~raw.exp.isin(range(6)) | ~raw.sol.isin(range(6)), "exclusion_reason"] = (
        "noninteger_score"
    )
    raw.loc[repeats, "exclusion_reason"] = "repeated_sequence_consistent_scores"
    raw.loc[raw.conflicting_scores, "exclusion_reason"] = (
        "repeated_sequence_conflicting_scores"
    )
    cohort = raw[raw.exclusion_reason == "eligible"].copy()
    cohort["label"] = (cohort.sol >= 4).astype(int)
    policies = read(p / "policy_groups.tsv.gz")
    old = policies[policies.policy == "NESG_unique_E45_S45_vs0"]
    assert set(cohort.record_id) == set(old.record_id)
    cohort = cohort.drop(
        columns=["audit_cluster", "audit_fold"], errors="ignore"
    ).merge(
        old[["record_id", "audit_cluster", "audit_fold", "label"]].rename(
            columns={"label": "old_label"}
        ),
        on="record_id",
        validate="one_to_one",
    )
    assert (cohort.label == pd.to_numeric(cohort.old_label)).all()
    cohort = cohort.drop(columns="old_label")
    features = read(p / "sequence_features.tsv.gz")
    featurecols = [
        "sha256",
        "charge_proxy",
        "hydropathy",
        "terminal_his",
        "pfam_annotated",
    ] + ["frac_" + a for a in AA]
    cohort = cohort.merge(features[featurecols], on="sha256", validate="one_to_one")
    assert len(cohort) == 3762 and cohort.label.value_counts().to_dict() == {
        1: 2614,
        0: 1148,
    }
    for col in ["charge_proxy", "hydropathy"] + ["frac_" + a for a in AA]:
        cohort[col] = pd.to_numeric(cohort[col], errors="coerce")
    domains = read(p / "domains_with_clans.tsv.gz")
    domains = domains[domains.sha256.isin(cohort.sha256)].copy()
    assert (
        (domains.sequence_score >= domains.sequence_ga)
        & (domains.domain_score >= domains.domain_ga)
    ).all()
    assert ((domains.ali_from >= 1) & (domains.ali_to <= domains.sequence_length)).all()
    edges = pd.read_csv(
        p / "homology_edges.tsv.gz",
        sep="\t",
        names=["query", "target", "identity", "qcov", "tcov", "evalue", "bits"],
    )
    assert (
        (edges.identity >= 0.30)
        & (edges.qcov >= 0.8)
        & (edges.tcov >= 0.8)
        & (edges.evalue <= 0.001)
    ).all()
    save(
        raw.drop(columns="sequence"),
        results / "restricted_prelock/raw_record_decisions.tsv.gz",
    )
    (results / "input_checksums.json").write_text(json.dumps(checks, indent=2) + "\n")
    return raw, curated, cohort, domains, edges


def freeze(cohort, domains, edges, results):
    """One-time split creation/audit; materialize development-only inputs, then write the lock."""
    assert not (results / "test_lock.json").exists(), (
        "Never overwrite finalized splits."
    )
    c = cohort.copy()
    c["split"] = c.audit_fold.map(
        {0: "reserved_test", 1: "validation", 2: "train", 3: "train", 4: "train"}
    )
    assert c.split.notna().all() and c.sha256.is_unique
    assert c.groupby("audit_cluster").split.nunique().max() == 1
    # Training-only inner CV, predetermined hash namespace; no search or balancing.
    c["inner_cv_fold"] = -1
    train = c.split == "train"
    c.loc[train, "inner_cv_fold"] = c.loc[train, "audit_cluster"].map(
        lambda s: int(seqhash("nesg-primary-inner-v1:" + s)[:16], 16) % 5
    )
    assert c.loc[train].groupby("audit_cluster").inner_cv_fold.nunique().max() == 1
    lookup = c.set_index("sha256").split
    e = edges[edges["query"].isin(lookup.index) & edges.target.isin(lookup.index)]
    cross = e["query"].map(lookup) != e.target.map(lookup)
    assert not cross.any()
    assert c.groupby("split").label.nunique().eq(2).all()
    # Confirm original component mapping against all recorded edges, including paths through outside proteins.
    old = (
        read(
            results.parent.parent / "dataset_comparison/results/sequence_groups.tsv.gz"
        )
        .set_index("sha256")
        .audit_cluster
    )
    assert (edges["query"].map(old) == edges.target.map(old)).all()
    cols = [
        "record_id",
        "sid",
        "sha256",
        "length",
        "exp",
        "sol",
        "label",
        "audit_cluster",
        "audit_fold",
        "split",
        "inner_cv_fold",
    ]
    restricted = results / "reserved_test"
    restricted.mkdir(exist_ok=True)
    save(
        c[cols].sort_values("record_id"),
        restricted / "final_cohort_and_split_manifest.tsv",
    )
    save(
        c[c.split == "reserved_test"][cols].sort_values("record_id"),
        restricted / "test_manifest.tsv",
    )
    save(
        c[c.split == "reserved_test"][["record_id", "sha256", "sequence"]],
        restricted / "test_sequences.tsv.gz",
    )
    # Public assignment file intentionally omits test phenotype and sequence.
    save(
        c[
            ["record_id", "sha256", "audit_cluster", "split", "inner_cv_fold"]
        ].sort_values("record_id"),
        results / "split_assignments.tsv",
    )
    dev = c[c.split != "reserved_test"].copy()
    save(dev, results / "development_records.tsv.gz")
    save(dev[cols], results / "development_manifest.tsv")
    save(dev[dev.split == "train"][cols], results / "training_manifest.tsv")
    save(dev[dev.split == "validation"][cols], results / "validation_manifest.tsv")
    save(
        dev[dev.split == "train"][["record_id", "audit_cluster", "inner_cv_fold"]],
        results / "training_cv_manifest.tsv",
    )
    dd = domains[domains.sha256.isin(dev.sha256)]
    save(dd, results / "development_domains.tsv.gz")
    memberships = domains[
        ["sha256", "pfam", "clan", "family_name", "description"]
    ].drop_duplicates(["sha256", "pfam"])
    rows = []
    for split, q in c.groupby("split"):
        mm = q[["sha256"]].merge(memberships, on="sha256")
        rows.append(
            dict(
                split=split,
                N=len(q),
                soluble=int(q.label.sum()),
                insoluble=int((q.label == 0).sum()),
                percent_soluble=100 * q.label.mean(),
                groups=q.audit_cluster.nunique(),
                families=mm.pfam.nunique(),
                clans=mm.loc[mm.clan != "", "clan"].nunique(),
                annotated=q.sha256.isin(mm.sha256).sum(),
            )
        )
    save(
        pd.DataFrame(rows).set_index("split").loc[ORDER].reset_index(),
        results / "initial_split_counts.tsv",
    )
    sets = {s: set(c[c.split == s].sha256) for s in ORDER}
    famsets = {
        s: set(memberships.loc[memberships.sha256.isin(sets[s]), "pfam"]) for s in ORDER
    }
    familyrep = []
    for s in ORDER:
        familyrep.append(
            dict(
                split=s,
                families=len(famsets[s]),
                shared_with_train=len(famsets[s] & famsets["train"]),
                not_in_train=len(famsets[s] - famsets["train"]),
            )
        )
    save(pd.DataFrame(familyrep), results / "initial_split_family_representation.tsv")
    cv = (
        dev[dev.split == "train"]
        .groupby("inner_cv_fold")
        .agg(
            N=("label", "size"),
            soluble=("label", "sum"),
            groups=("audit_cluster", "nunique"),
        )
        .reset_index()
    )
    cv["insoluble"] = cv.N - cv.soluble
    assert (cv[["soluble", "insoluble"]] > 0).all().all()
    save(cv, results / "inner_cv_counts.tsv")
    integrity = dict(
        cohort_records=len(c),
        unique_sequences=c.sha256.nunique(),
        soluble=int(c.label.sum()),
        insoluble=int((c.label == 0).sum()),
        exact_cross_split_sequences=0,
        cross_split_groups=0,
        cohort_direct_edges=len(e),
        cross_split_detected_edges=int(cross.sum()),
        all_union_edges_verified=len(edges),
        training_cv_group_leakage=0,
        original_strict_membership_and_scores_unchanged=True,
    )
    (results / "integrity_checks.json").write_text(
        json.dumps(integrity, indent=2) + "\n"
    )
    frozen = [p for p in restricted.iterdir() if p.is_file()] + list(
        (results / "restricted_prelock").glob("*")
    )
    safe = [p for p in results.iterdir() if p.is_file() and p.name != "test_lock.json"]
    lock = dict(
        status="LOCKED",
        created_UTC=pd.Timestamp.now(tz="UTC").isoformat(),
        policy="Existing audit_fold 0=test, 1=validation, 2/3/4=train; no candidate search",
        test_fraction=len(c[c.split == "reserved_test"]) / len(c),
        previous_exposure="Full-cohort descriptive audit previously performed. This is a prospective lock, not a historically untouched holdout.",
        restricted_checksums={str(p.relative_to(results)): digest(p) for p in frozen},
        rerun_safe_checksums={str(p.relative_to(results)): digest(p) for p in safe},
    )
    (results / "test_lock.json").write_text(json.dumps(lock, indent=2) + "\n")
    return dev, dd, lock
