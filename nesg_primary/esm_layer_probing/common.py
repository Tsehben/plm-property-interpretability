"""Shared Stage 1 protocol. No test sequences or test labels are read."""

from pathlib import Path
import hashlib
import json
import os
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
MODEL = "facebook/esm2_t33_650M_UR50D"
REVISION = "08e4846e537177426273712802403f7ba8261b6c"
WEIGHTS_SHA256 = "a08adabb949fa67ad3c14b509d04fd60368b35007b0095e3358f81200c4f4db0"
INTERPROT_REVISION = "3ded2e51641df84d6b89f25430874e3a8eb42c24"
SAE_LAYERS = [4, 8, 12, 16, 20, 24, 28, 32, 33]
PROTOCOL = dict(
    model=MODEL,
    revision=REVISION,
    layers=list(range(1, 34)),
    width=1280,
    indexing="hidden_states[L] for L=1..33; layer 0 is embedding output",
    pooling="float32 mean of residue positions 1..length; no BOS, EOS or padding",
    C_grid=[0.0001, 0.001, 0.01, 0.1, 1.0, 10.0],
    seed=20261118,
    bootstrap_replicates=2000,
    representative_layers=[8, 16, 32],
    class_weight=None,
    solver="lbfgs",
    max_iter=3000,
    counts={"train": 2361, "validation": 691},
    max_residues=1022,
    baseline="20 standard amino-acid counts / original length plus natural log length; X not imputed",
    SAE_available_layers=SAE_LAYERS,
    interprot_revision=INTERPROT_REVISION,
)
AA = "ACDEFGHIKLMNPQRSTVWY"


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")
    temp.replace(path)


def token():
    value = os.environ.get("HF_TOKEN_PLM")
    if not value:
        raise RuntimeError(
            "HF_TOKEN_PLM is unavailable. Set it securely from Colab Secrets."
        )
    return value


def read_extraction(path=None):
    path = Path(path or HERE / "data/extraction_manifest.tsv")
    pin = json.loads((path.parent / "input_integrity.json").read_text())
    assert sha(path) == pin["extraction_manifest_sha256"], "Input manifest changed"
    d = pd.read_csv(path, sep="\t", keep_default_na=False)
    assert d.split.value_counts().to_dict() == PROTOCOL["counts"]
    assert len(d) == 3052 and d.record_id.is_unique and d.sha256.is_unique
    assert d.sequence.map(lambda s: hashlib.sha256(s.encode()).hexdigest()).equals(
        d.sha256
    )
    assert (d.sequence.str.len() == d.length).all() and d.length.between(1, 1022).all()
    assert d.sequence.map(lambda s: set(s) <= set(AA + "X")).all()
    assert d.sequence.str.contains("X").sum() == 3
    return d


def read_local():
    d = read_extraction()
    p = HERE / "data/local_labels.tsv"
    pin = json.loads((HERE / "data/input_integrity.json").read_text())
    assert sha(p) == pin["local_labels_sha256"]
    labels = pd.read_csv(p, sep="\t")
    out = d.merge(
        labels, on=["record_id", "sha256", "split"], validate="one_to_one", sort=False
    )
    assert len(out) == 3052 and out.groupby("audit_cluster").split.nunique().max() == 1
    tr = out[out.split == "train"]
    assert tr.groupby("audit_cluster").inner_cv_fold.nunique().max() == 1
    return out


def composition(sequences):
    return np.array(
        [[s.count(a) / len(s) for a in AA] + [np.log(len(s))] for s in sequences],
        dtype=np.float64,
    )
