"""Check deliverables and record checksums, without reopening upstream/test records."""

import ast
import json
import re
import subprocess
import sys
import nbformat
import pandas as pd
import numpy as np
from common import HERE, read_local, sha, write_json


def main():
    data = read_local()
    integrity = json.loads((HERE / "data/input_integrity.json").read_text())
    for name, key in [
        ("pfam_domains.tsv.gz", "pfam_domains_sha256"),
        ("pfam_memberships.tsv", "pfam_memberships_sha256"),
    ]:
        assert sha(HERE / "data" / name) == integrity[key]
    domains = pd.read_csv(HERE / "data/pfam_domains.tsv.gz", sep="\t")
    assert set(domains.record_id) <= set(data.record_id)
    assert set(domains.sha256) <= set(data.sha256)
    checked = domains.merge(
        data[["record_id", "sha256", "split"]],
        on=["record_id", "sha256", "split"],
        validate="many_to_one",
    )
    assert len(checked) == len(domains)
    for name in ["ESM_Extraction_Colab.ipynb", "ESM_Layerwise_Solubility.ipynb"]:
        nb = nbformat.read(HERE / name, as_version=4)
        nbformat.validate(nb)
        for cell in nb.cells:
            if cell.cell_type == "code":
                ast.parse(cell.source)
            assert not any(o.output_type == "error" for o in cell.get("outputs", []))
    colab = nbformat.read(HERE / "ESM_Extraction_Colab.ipynb", as_version=4)
    cells = "\n".join(c.source for c in colab.cells if c.cell_type == "code")
    assert "drive.mount('/content/drive')" in cells
    assert "files.upload" not in cells and "colab_bundle.zip" not in cells
    assert "OUTPUT = Path('/content/nesg_esm_stage1_scratch/activations')" in cells
    assert "'--hf-repo', DESTINATION_REPO" in cells
    metrics = pd.read_csv(HERE / "results/metrics.tsv", sep="\t")
    expected = data[data.split == "validation"][
        ["record_id", "sha256", "split", "label", "audit_cluster", "inner_cv_fold"]
    ].reset_index(drop=True)
    for model in metrics.model:
        pred = pd.read_csv(
            HERE / f"results/{model}_validation_predictions.tsv", sep="\t"
        )
        pd.testing.assert_frame_equal(
            pred[expected.columns], expected, check_dtype=False
        )
        assert np.isfinite(pred[["decision_score", "probability"]]).all().all()
    tests = subprocess.run(
        [sys.executable, str(HERE / "test_stage1.py")], capture_output=True, text=True
    )
    (HERE / "results/test_results.txt").write_text(tests.stdout + tests.stderr)
    assert tests.returncode == 0, "See results/test_results.txt"
    freeze = subprocess.run(
        [sys.executable, "-m", "pip", "freeze"],
        capture_output=True,
        text=True,
        check=True,
    )
    (HERE / "results/environment_freeze.txt").write_text(freeze.stdout)
    output = {}
    for path in sorted(HERE.rglob("*")):
        if (
            not path.is_file()
            or "__pycache__" in path.parts
            or ".ipynb_checkpoints" in path.parts
            or path.name == "artifact_checksums.json"
        ):
            continue
        if path.suffix in {".py", ".md", ".tsv", ".json", ".ipynb", ".txt"}:
            assert not re.search(r"hf_[A-Za-z0-9]{20,}", path.read_text()), (
                f"Possible credential found in {path.name}"
            )
        output[str(path.relative_to(HERE))] = sha(path)
    write_json(HERE / "artifact_checksums.json", output)
    print(
        f"Passed: {len(data)} aligned development proteins, {len(domains)} domain rows, notebooks, Drive-input/HF-output configuration, predictions, tests and secret-pattern scan. Hashed {len(output)} files."
    )


if __name__ == "__main__":
    main()
