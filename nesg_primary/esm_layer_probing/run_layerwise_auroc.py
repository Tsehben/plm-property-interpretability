"""Colab entry point: CPU probes from existing pooled ESM, plots, private HF results.

Never loads ESM weights or full residue maps. Inputs/labels remain separate local
files; published prediction tables omit labels, scores E/S, homology and CV groups.
"""

import argparse
import hashlib
import importlib.metadata as metadata
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import zipfile

PINS = {
    "huggingface-hub": "0.36.0",
    "scipy": "1.15.3",
    "scikit-learn": "1.7.2",
    "joblib": "1.5.2",
    "threadpoolctl": "3.6.0",
    "matplotlib": "3.10.7",
}


def ensure_packages():
    """Preserve working Torch/NumPy; install only missing probing wheels in venv."""
    assert sys.prefix != sys.base_prefix, "Use EXTRACTION_PYTHON via run_stage1."
    missing = []
    for name, version in PINS.items():
        try:
            installed = metadata.version(name)
        except metadata.PackageNotFoundError:
            installed = None
        if installed != version:
            missing.append(f"{name}=={version}")
    if missing:
        import numpy
        import pandas

        protected = dict(numpy=numpy.__version__, pandas=pandas.__version__)
        try:
            protected["torch"] = metadata.version("torch")
        except metadata.PackageNotFoundError:
            pass  # CPU probing and pooled-cache validation do not require Torch.
        constraints = Path(sys.prefix) / "nesg_probe_constraints.txt"
        constraints.write_text(
            "\n".join(f"{k}=={v}" for k, v in protected.items()) + "\n"
        )
        print("Installing only probing/plotting packages:", missing, flush=True)
        subprocess.run(
            [
                sys._base_executable,
                "-m",
                "pip",
                "--python",
                sys.executable,
                "install",
                "--disable-pip-version-check",
                "--only-binary=:all:",
                "--constraint",
                str(constraints),
                "--timeout",
                "30",
                "--retries",
                "2",
                *missing,
            ],
            check=True,
        )
    else:
        print("Probing packages already installed; skipping pip.", flush=True)


def make_bundle(results):
    """Publish enough to reuse predictions/probes, without uploading experimental labels."""
    import pandas as pd

    results = Path(results)
    buffer = io.BytesIO()
    checksums = {}
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(results.rglob("*")):
            if not path.is_file() or path.suffix not in {
                ".json",
                ".tsv",
                ".npz",
                ".joblib",
                ".png",
                ".svg",
                ".md",
            }:
                continue
            # Local completion markers hash local labelled prediction tables;
            # these are not the public bundle's identity/validation mechanism.
            if path.name.endswith("_complete.json"):
                continue
            name = str(path.relative_to(results))
            raw = path.read_bytes()
            if path.suffix == ".tsv":
                table = pd.read_csv(path, sep="\t", float_precision="round_trip")
                table = table.drop(
                    columns=[
                        "label",
                        "audit_cluster",
                        "inner_cv_fold",
                        "sequence",
                        "E",
                        "S",
                    ],
                    errors="ignore",
                )
                raw = table.to_csv(sep="\t", index=False).encode()
            checksums[name] = hashlib.sha256(raw).hexdigest()
            archive.writestr(
                zipfile.ZipInfo(name), raw, compress_type=zipfile.ZIP_DEFLATED
            )
        archive.writestr(
            zipfile.ZipInfo("bundle_checksums.json"),
            json.dumps(checksums, sort_keys=True),
        )
    return buffer.getvalue()


def unpack_bundle(raw, destination):
    destination = Path(destination)
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        names = archive.namelist()
        assert len(names) == len(set(names)), "Duplicate ZIP members"
        hashes = json.loads(archive.read("bundle_checksums.json"))
        assert set(names) == set(hashes) | {"bundle_checksums.json"}
        for name, digest in hashes.items():
            assert not Path(name).is_absolute() and ".." not in Path(name).parts
            content = archive.read(name)
            assert hashlib.sha256(content).hexdigest() == digest
            target = destination / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(content)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--activations", type=Path)
    parser.add_argument("--repo")
    parser.add_argument("--setup-only", action="store_true")
    parser.add_argument(
        "--output", type=Path, default=Path("/content/nesg_layerwise_probes")
    )
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    ensure_packages()
    if args.setup_only:
        print(
            "CPU probing environment ready; no model weights or activations loaded.",
            flush=True,
        )
        return
    if not args.activations or not args.repo:
        parser.error("--activations and --repo are required for probing")
    # Re-exec imports after any installs, so the numerical stack cannot be stale.
    command = [
        sys.executable,
        "-u",
        str(Path(__file__).resolve()),
        "--worker",
        "--activations",
        str(args.activations),
        "--repo",
        args.repo,
        "--output",
        str(args.output),
        "--threads",
        str(args.threads),
    ]
    subprocess.run(command, check=True)


def worker():
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--activations", required=True, type=Path)
    parser.add_argument("--repo", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--threads", type=int, required=True)
    args = parser.parse_args()
    from common import HERE, read_local, sha, token, write_json
    from hf_storage import private_api, verified_files
    from huggingface_hub import CommitOperationAdd, hf_hub_download
    from save_full_maps import remote

    data = read_local()
    assert data.split.value_counts().to_dict() == {"train": 2361, "validation": 691}
    assert set(data.loc[data.split == "train", "audit_cluster"]).isdisjoint(
        data.loc[data.split == "validation", "audit_cluster"]
    )
    verified_files(args.activations)
    identity = dict(
        inputs={
            p: sha(HERE / p)
            for p in [
                "common.py",
                "probe.py",
                "plot_train_validation.py",
                "run_layerwise_auroc.py",
                "protocol.json",
                "data/extraction_manifest.tsv",
                "data/local_labels.tsv",
                "data/input_integrity.json",
            ]
        },
        activation_completion_sha256=sha(args.activations / "completed.json"),
        software={p: metadata.version(p) for p in ["numpy", "pandas", *PINS]},
        train=2361,
        validation=691,
        test_evaluated=False,
    )
    fingerprint = hashlib.sha256(
        json.dumps(identity, sort_keys=True).encode()
    ).hexdigest()
    prefix = f"esm_stage1/layerwise_probing_v1/{fingerprint}"
    destination = args.output / fingerprint
    destination.mkdir(parents=True, exist_ok=True)
    api, info = private_api(args.repo)
    names = set(
        remote(api.list_repo_files, args.repo, repo_type="dataset", revision=info.sha)
    )
    if prefix + "/completed.json" in names:
        index = Path(
            remote(
                hf_hub_download,
                args.repo,
                prefix + "/completed.json",
                repo_type="dataset",
                revision=info.sha,
                token=token(),
            )
        )
        complete = json.loads(index.read_text())
        assert complete["identity"] == identity and complete["models"] == 34
        archive = Path(
            remote(
                hf_hub_download,
                args.repo,
                prefix + "/results.zip",
                repo_type="dataset",
                revision=info.sha,
                token=token(),
            )
        )
        assert sha(archive) == complete["sha256"]
        unpack_bundle(archive.read_bytes(), destination / "published")
        print(
            "SKIP fitting: complete matching probes and figures restored from private HF.",
            flush=True,
        )
    else:
        results = destination / "local_results"
        # CPU linear probes only. Existing per-layer checkpoints survive a
        # notebook reopen in the same runtime and are verified before reuse.
        subprocess.run(
            [
                sys.executable,
                "-u",
                str(HERE / "probe.py"),
                "--activations",
                str(args.activations),
                "--output",
                str(results),
                "--resume",
                "--threads",
                str(args.threads),
            ],
            check=True,
        )
        from plot_train_validation import plot

        plot(results)
        write_json(results / "reproducibility.json", identity)
        raw = make_bundle(results)
        archive = destination / "results.zip"
        archive.write_bytes(raw)
        complete = dict(
            identity=identity,
            models=34,
            sha256=sha(archive),
            label_policy="Experimental labels/sequences and homology/CV group IDs omitted from published protein predictions; join by stable ID/hash to local manifest.",
        )
        write_json(destination / "completed.json", complete)
        current = remote(api.repo_info, args.repo, repo_type="dataset")
        assert current.private is True
        remote(
            api.create_commit,
            repo_id=args.repo,
            repo_type="dataset",
            parent_commit=current.sha,
            operations=[
                CommitOperationAdd(
                    path_in_repo=prefix + "/results.zip", path_or_fileobj=str(archive)
                ),
                CommitOperationAdd(
                    path_in_repo=prefix + "/completed.json",
                    path_or_fileobj=str(destination / "completed.json"),
                ),
            ],
            commit_message="Save NESG 33-layer train/validation linear probes and AUROC figures",
        )
        unpack_bundle(raw, destination / "published")
        print(
            "Complete: 33 ESM probes + composition baseline; results committed to private HF.",
            flush=True,
        )
    published = destination / "published"
    write_json(
        args.output / "latest.json",
        dict(results=str(published), repository=args.repo, hf_prefix=prefix),
    )
    print((published / "TRAIN_VALIDATION_AUROC.md").read_text())
    print("Figures:", published / "figures", flush=True)
    print("Private HF results:", f"{args.repo}/{prefix}/results.zip", flush=True)


if __name__ == "__main__":
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/nesg-probe-matplotlib")
    if "--worker" in sys.argv:
        worker()
    else:
        main()
