"""Preflight the GPU, then install extraction packages into a separate venv.

The venv reads Colab's existing GPU torch/numerical stack through system site
packages. Pip writes only into that venv; it never installs into the notebook
interpreter, and --only-binary prevents source compilation.
"""

import argparse
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import venv

HERE = Path(__file__).resolve().parent
DEFAULT_ENV = Path("/content/nesg_esm_env")
PROTECTED = {
    "torch",
    "torchvision",
    "torchaudio",
    "triton",
    "numpy",
    "pandas",
    "scipy",
    "scikit-learn",
    "psutil",
}
PREFLIGHT = r"""
import json, sys
from importlib import metadata
import torch
if torch.version.cuda is None or not torch.cuda.is_available():
    raise RuntimeError(f'GPU required BEFORE installing anything: torch={torch.__version__}, compiled CUDA={torch.version.cuda}, CUDA available={torch.cuda.is_available()}. Select a fresh GPU runtime; do not repair this CPU installation in place.')
# Actually allocate on the device, rather than trusting a package suffix.
torch.zeros(1, device='cuda'); torch.cuda.synchronize()
import numpy, pandas, psutil
if sys.version_info >= (3, 13) and int(numpy.__version__.split('.')[0]) < 2:
    raise RuntimeError('Python 3.13 with NumPy 1.x is not a supported base for this workflow. Use a fresh GPU runtime.')
for name, module in [('torch',torch), ('numpy',numpy), ('pandas',pandas), ('psutil',psutil)]:
    recorded = metadata.version(name)
    if recorded.split('+')[0] != module.__version__.split('+')[0]:
        raise RuntimeError(f'Inconsistent package metadata for {name}: recorded={recorded}, imported={module.__version__}. Use a fresh GPU runtime.')
print(json.dumps({'python':sys.version.split()[0], 'torch':torch.__version__, 'cuda':torch.version.cuda, 'gpu':torch.cuda.get_device_name(), 'numpy':numpy.__version__, 'pandas':pandas.__version__, 'psutil':psutil.__version__}))
"""


def child_environment():
    env = os.environ.copy()
    # ESM uses torch only; prevent incidental TensorFlow/Flax imports.
    env.update(USE_TORCH="1", USE_TF="0", USE_FLAX="0", PYTHONNOUSERSITE="1")
    return env


def gpu_preflight():
    print(
        "[1/3] Checking the existing GPU runtime BEFORE package changes...", flush=True
    )
    result = subprocess.run(
        [sys.executable, "-u", "-c", PREFLIGHT],
        capture_output=True,
        text=True,
        timeout=120,
        env=child_environment(),
    )
    if result.returncode:
        print(result.stderr, file=sys.stderr, flush=True)
        raise RuntimeError(
            "GPU/runtime preflight failed. No packages were installed. Disconnect/delete this runtime and start a fresh GPU runtime."
        )
    info = json.loads(result.stdout.strip().splitlines()[-1])
    print("GPU preflight passed:", info, flush=True)
    return info


def active_versions(imported_versions=None):
    """Preserve priority-resolved versions, checking core packages against imports.

    Colab can expose several dist-info entries for one package. Enumeration order
    is not import precedence: an inactive older entry must not become a pin or
    force a runtime reset. Do not delete or rewrite shared package metadata.
    """
    imported_versions = imported_versions or {}
    versions = {}
    seen = {}
    for distribution in metadata.distributions():
        name = distribution.metadata.get("Name", "").lower().replace("_", "-")
        if name in PROTECTED or name.startswith("nvidia-"):
            seen.setdefault(name, set()).add(distribution.version)
    for name, candidates in seen.items():
        selected = metadata.version(name)
        imported = imported_versions.get(name)
        if imported is not None and selected.split("+")[0] != imported.split("+")[0]:
            raise RuntimeError(
                f"Active metadata/import mismatch for {name}: metadata={selected}, imported={imported}. "
                "Stopping before installation; duplicate entries alone are not an error."
            )
        if len(candidates) > 1:
            print(
                f"Multiple metadata entries for {name}: {sorted(candidates)}; "
                f"keeping active version {selected}"
                + (
                    " (matches imported version)."
                    if imported is not None
                    else " (resolved by interpreter metadata lookup)."
                ),
                flush=True,
            )
        versions[name] = selected
    return versions


def pip_command(python, constraints):
    return [
        sys.executable,
        "-m",
        "pip",
        "--python",
        str(python),
        "install",
        "--disable-pip-version-check",
        "--only-binary=:all:",
        "--timeout",
        "30",
        "--retries",
        "2",
        "--upgrade-strategy",
        "only-if-needed",
        "--constraint",
        str(constraints),
        "-r",
        str(HERE / "requirements-colab.txt"),
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--env-dir", type=Path, default=DEFAULT_ENV)
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    started = time.monotonic()
    # Nothing that mutates a Python environment may precede these checks.
    base_info = gpu_preflight()
    base_versions = active_versions(base_info)
    if args.preflight_only:
        print("Preflight complete; no environment created or packages installed.")
        return
    env_dir = args.env_dir.resolve()
    assert env_dir != Path(sys.prefix).resolve(), (
        "Refusing the active notebook environment"
    )
    marker = env_dir / "nesg_environment.json"
    identity = dict(
        base_python=str(Path(sys.executable).resolve()),
        base_info=base_info,
        protected_versions=base_versions,
    )
    if env_dir.exists():
        if not marker.exists() or json.loads(marker.read_text()) != identity:
            raise RuntimeError(
                "Existing extraction environment belongs to another runtime. Start a fresh runtime or choose another --env-dir."
            )
    else:
        # No ensurepip or apt needed: the host pip manages this interpreter with --python.
        venv.EnvBuilder(system_site_packages=True, with_pip=False).create(env_dir)
        marker.write_text(json.dumps(identity, indent=2) + "\n")
    python = env_dir / "bin/python"
    assert python.exists()
    constraints = env_dir / "keep-base-versions.txt"
    constraints.write_text(
        "\n".join(
            f"{name}=={version}" for name, version in sorted(base_versions.items())
        )
        + "\n"
    )
    # Query the actual dedicated interpreter, not a success flag or the host's
    # metadata. If pins already match, no pip resolver/network invocation is needed.
    check_pins = """
from importlib.metadata import version
from pathlib import Path
import sys
for line in Path(sys.argv[1]).read_text().splitlines():
    if line.strip() and not line.lstrip().startswith('#'):
        name, expected = line.strip().split('==')
        assert version(name) == expected
"""
    pinned = (
        subprocess.run(
            [str(python), "-c", check_pins, str(HERE / "requirements-colab.txt")],
            capture_output=True,
            timeout=30,
            env=child_environment(),
        ).returncode
        == 0
    )
    print(
        "[2/3] "
        + (
            "SKIP pip: pinned packages already installed in "
            if pinned
            else "Installing binary extraction packages into "
        )
        + str(env_dir),
        flush=True,
    )
    print(
        "The notebook Python environment is not an installation target; no source builds are allowed.",
        flush=True,
    )
    if not pinned:
        subprocess.run(
            pip_command(python, constraints), check=True, env=child_environment()
        )
    assert active_versions(base_info) == base_versions, (
        "Base package metadata changed unexpectedly"
    )
    print(
        f"[3/3] Testing ESM imports in the extraction environment (elapsed {time.monotonic() - started:.1f}s)...",
        flush=True,
    )
    check = r"""
import json, sys
from pathlib import Path
from importlib import metadata
import torch, numpy, pandas, psutil
from transformers import EsmModel, AutoTokenizer
from huggingface_hub import HfApi
assert torch.cuda.is_available() and torch.version.cuda is not None
assert sys.prefix != sys.base_prefix, 'Must run inside the extraction venv'
packages = ['torch','numpy','pandas','psutil','transformers','huggingface-hub','safetensors','tokenizers']
versions = {name:metadata.version(name) for name in packages}
for line in Path(sys.argv[1]).read_text().splitlines():
    if line.strip() and not line.lstrip().startswith('#'):
        name,expected = line.strip().split('=='); assert versions[name] == expected
Path(sys.argv[2]).write_text(json.dumps(versions, indent=2)+'\n')
print('ESM imports and CUDA passed:', versions, flush=True)
"""
    subprocess.run(
        [
            str(python),
            "-u",
            "-c",
            check,
            str(HERE / "requirements-colab.txt"),
            str(HERE / "colab_environment.json"),
        ],
        check=True,
        timeout=180,
        env=child_environment(),
    )
    print("Ready. All extraction/HF commands must use:", python, flush=True)


if __name__ == "__main__":
    main()
