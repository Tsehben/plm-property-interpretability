"""Execute the real notebook with this Python environment; preserve errors for inspection."""

from pathlib import Path
import os
import sys
import json
import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpec

ROOT = Path(__file__).resolve().parent
os.environ["MPLCONFIGDIR"] = "/tmp/nesg-validation-matplotlib"
os.environ["IPYTHONDIR"] = "/tmp/nesg-validation-ipython"
os.environ["JUPYTER_RUNTIME_DIR"] = "/tmp/nesg-validation-jupyter"
os.environ["OPENBLAS_NUM_THREADS"] = "2"
for p in ["MPLCONFIGDIR", "IPYTHONDIR", "JUPYTER_RUNTIME_DIR"]:
    Path(os.environ[p]).mkdir(exist_ok=True)
km = KernelManager(kernel_name="python3")
km._kernel_spec = KernelSpec(
    argv=[sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
    display_name="NESG validation",
    language="python",
)
p = ROOT / "NESG_Dataset_Validation.ipynb"
nb = nbformat.read(p, as_version=4)
client = NotebookClient(
    nb, km=km, timeout=300, resources={"metadata": {"path": str(ROOT)}}
)
try:
    client.execute()
finally:
    nbformat.write(nb, p)
errors = [
    (i, o.get("ename"), o.get("evalue"))
    for i, c in enumerate(nb.cells)
    for o in c.get("outputs", [])
    if o.output_type == "error"
]
images = sum(
    "image/png" in o.get("data", {}) for c in nb.cells for o in c.get("outputs", [])
)
assert not errors, errors
assert images >= 10, images
report = dict(
    code_cells=sum(c.cell_type == "code" for c in nb.cells),
    executed_cells=sum(
        c.cell_type == "code" and c.execution_count is not None for c in nb.cells
    ),
    inline_png_outputs=images,
    errors=errors,
    locked=(ROOT / "results/test_lock.json").exists(),
)
(ROOT / "results/notebook_execution.json").write_text(
    json.dumps(report, indent=2) + "\n"
)
print(json.dumps(report, indent=2))
