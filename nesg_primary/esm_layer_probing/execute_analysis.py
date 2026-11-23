"""Execute the analysis notebook with this interpreter, retaining real outputs."""

import importlib.metadata
import json
import os
from pathlib import Path
import sys
import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpec
from common import HERE, write_json

for key, suffix in [
    ("MPLCONFIGDIR", "matplotlib"),
    ("IPYTHONDIR", "ipython"),
    ("JUPYTER_RUNTIME_DIR", "jupyter"),
]:
    os.environ[key] = f"/tmp/nesg-esm-{suffix}"
    Path(os.environ[key]).mkdir(exist_ok=True)
os.environ["OPENBLAS_NUM_THREADS"] = "2"
km = KernelManager(kernel_name="python3")
km._kernel_spec = KernelSpec(
    argv=[sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
    display_name="NESG Stage 1",
    language="python",
)
p = HERE / "ESM_Layerwise_Solubility.ipynb"
nb = nbformat.read(p, as_version=4)
try:
    NotebookClient(
        nb, km=km, timeout=7200, resources={"metadata": {"path": str(HERE)}}
    ).execute()
finally:
    nbformat.write(nb, p)
errors = [
    (i, o.get("ename"))
    for i, c in enumerate(nb.cells)
    for o in c.get("outputs", [])
    if o.output_type == "error"
]
assert not errors, errors
report = dict(
    code_cells=sum(c.cell_type == "code" for c in nb.cells),
    executed_cells=sum(
        c.cell_type == "code" and c.execution_count is not None for c in nb.cells
    ),
    inline_png_outputs=sum(
        "image/png" in o.get("data", {}) for c in nb.cells for o in c.get("outputs", [])
    ),
    errors=errors,
    python=sys.version.split()[0],
    versions={
        p: importlib.metadata.version(p)
        for p in [
            "torch",
            "transformers",
            "huggingface-hub",
            "numpy",
            "pandas",
            "scikit-learn",
            "scipy",
            "matplotlib",
            "nbformat",
            "nbclient",
        ]
    },
)
write_json(HERE / "results/notebook_execution.json", report)
print(json.dumps(report, indent=2))
