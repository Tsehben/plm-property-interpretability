"""Save all 33 residue-level ESM maps for the frozen 3052-protein development cohort.

Run after the pooled extraction. Uploads one safetensors file per protein to a
separate private HF prefix, with atomic batches of eight and verified resume.
No sequences, solubility labels, or reserved-test records are uploaded.
Float32 requires about 127 GB; optional float16 requires about 63.5 GB.
"""

import argparse
from pathlib import Path
import math
import re
import time
from email.utils import parsedate_to_datetime


def retry_delay(error, attempt):
    """Respect server reset hints; otherwise exponential backoff, 60s to 15min."""
    headers = getattr(getattr(error, "response", None), "headers", {}) or {}
    headers = {k.lower(): v for k, v in headers.items()}
    delays = []
    value = headers.get("retry-after", "")
    try:
        delays.append(float(value))
    except (ValueError, TypeError):
        try:
            delays.append(parsedate_to_datetime(value).timestamp() - time.time())
        except (ValueError, TypeError, OverflowError):
            pass
    # HF's standardized RateLimit header includes t=<seconds until reset>.
    for value in re.findall(
        r"(?:^|[;,])\s*t\s*=\s*(\d+)", headers.get("ratelimit", "")
    ):
        delays.append(float(value))
    delays = [d for d in delays if math.isfinite(d) and d >= 0]
    return max(1, math.ceil(max(delays)) + 2) if delays else min(60 * 2**attempt, 900)


def remote(function, *args, **kwargs):
    """Retry rate limits without printing exception bodies, URLs or credentials."""
    max_retries = 12
    for attempt in range(max_retries + 1):
        try:
            return function(*args, **kwargs)
        except Exception as error:
            status = getattr(getattr(error, "response", None), "status_code", None)
            if status != 429 or attempt == max_retries:
                raise RuntimeError(
                    f"HF {function.__name__} failed (HTTP {status}); check access, quota or connectivity. "
                    "Completed batches remain saved; rerun to resume."
                ) from None
            remaining = retry_delay(error, attempt)
            print(
                f"HF rate limit (429) during {function.__name__}; waiting {remaining}s "
                f"before retry {attempt + 1}/{max_retries}. Completed batches are saved.",
                flush=True,
            )
            while remaining > 0:
                step = min(60, remaining)
                time.sleep(step)
                remaining -= step
                if remaining > 0:
                    print(f"Rate-limit cooldown: {remaining}s remaining.", flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, help="Existing private HF dataset")
    parser.add_argument(
        "--pooled-dir",
        type=Path,
        required=True,
        help="Completed, locally cached pooled extraction directory",
    )
    parser.add_argument(
        "--cache",
        type=Path,
        default=Path.home() / ".cache/nesg-esm",
        help="Existing ESM model cache; no model weights saved in the repository",
    )
    parser.add_argument(
        "--dtype",
        choices=["float32", "float16"],
        default="float32",
        help="Saved map precision; inference matches the pooled extraction",
    )
    parser.add_argument(
        "--scratch-dir",
        type=Path,
        default=Path("/content/nesg_full_maps_scratch"),
        help="Temporary batch files; removed only after successful HF commits",
    )
    args = parser.parse_args(argv)
    import json
    from contextlib import nullcontext
    import numpy as np
    import torch
    from safetensors.torch import save_file
    from huggingface_hub import CommitOperationAdd, hf_hub_download
    from common import HERE, MODEL, REVISION, read_extraction, sha, token, write_json
    from extract import load_model, tokenize, versions, load_all, verify_config
    from hf_storage import private_api, git_blob_sha

    repo = args.repo
    pooled_dir = args.pooled_dir
    cache = str(args.cache)
    storage_dtype = args.dtype
    assert storage_dtype in ("float16", "float32")
    assert (pooled_dir / "completed.json").exists(), (
        "Finish the pooled extraction first."
    )
    assert torch.cuda.is_available()
    data = (
        read_extraction()
    )  # Verifies all 3052 IDs, sequences, hashes and train/val counts.
    pooled_config = verify_config(pooled_dir, data)
    software = versions()
    assert all(software.get(k) == v for k, v in pooled_config["software"].items()), (
        "Use the same extraction environment."
    )
    reference = load_all(pooled_dir, data)  # Also verifies every pooled shard checksum.
    precision = pooled_config["extraction_dtype"]
    assert precision in ("float32", "float16", "bfloat16")
    save_dtype = getattr(torch, storage_dtype)
    prefix = f"esm_stage1/full_maps_{storage_dtype}_v1"
    scratch = args.scratch_dir / storage_dtype
    scratch.mkdir(parents=True, exist_ok=True)
    api, initial = private_api(repo)
    credential = token()

    def fetch_json(path):
        local = remote(
            hf_hub_download,
            repo,
            path,
            repo_type="dataset",
            revision=initial.sha,
            token=credential,
            cache_dir=str(scratch / "metadata_cache"),
        )
        return json.loads(Path(local).read_text())

    def commit(files):
        current = remote(api.repo_info, repo, repo_type="dataset")
        assert current.private is True, "Destination must remain private."
        found = {
            r.path: r
            for r in remote(
                api.get_paths_info,
                repo,
                list(files),
                repo_type="dataset",
                revision=current.sha,
            )
        }
        operations = []
        for name, path in files.items():
            if name in found:
                r = found[name]
                same = (
                    r.lfs.sha256 == sha(path)
                    if r.lfs
                    else r.blob_id == git_blob_sha(path)
                )
                assert same, f"Refusing to overwrite a different artifact: {name}"
            else:
                operations.append(
                    CommitOperationAdd(path_in_repo=name, path_or_fileobj=str(path))
                )
        if operations:
            remote(
                api.create_commit,
                repo_id=repo,
                repo_type="dataset",
                parent_commit=current.sha,
                operations=operations,
                commit_message="Checkpoint NESG full residue activations",
            )

    config = dict(
        format_version=1,
        model=MODEL,
        revision=REVISION,
        inference_precision=precision,
        storage_dtype=storage_dtype,
        axes=["layer", "residue", "feature"],
        layers=list(range(1, 34)),
        width=1280,
        residue_indexing="index 0 = original residue 1; no BOS/EOS/padding",
        batch_size=1,
        proteins=3052,
        splits={"train": 2361, "validation": 691},
        manifest_sha256=sha(HERE / "data/extraction_manifest.tsv"),
        pooled_completion_sha256=sha(pooled_dir / "completed.json"),
        software=software,
    )
    write_json(scratch / "run.json", config)
    config_hash = sha(scratch / "run.json")
    existing = set(
        remote(api.list_repo_files, repo, repo_type="dataset", revision=initial.sha)
    )
    if prefix + "/run.json" in existing:
        assert fetch_json(prefix + "/run.json") == config, (
            "Different existing run; do not mix outputs."
        )
    else:
        assert not any(p.startswith(prefix + "/") for p in existing), (
            "Existing prefix lacks its run metadata."
        )
        commit({prefix + "/run.json": scratch / "run.json"})

    model = tokenizer = None
    all_entries = []
    for start in range(0, len(data), 8):
        indices = list(range(start, min(start + 8, len(data))))
        batch_name = f"{prefix}/batches/{start // 8:04}.json"
        if batch_name in existing:
            saved = fetch_json(batch_name)
            assert saved["config_sha256"] == config_hash
            entries = saved["proteins"]
            assert len(entries) == len(indices)
            records = {
                r.path: r
                for r in remote(
                    api.get_paths_info,
                    repo,
                    [e["path"] for e in entries],
                    repo_type="dataset",
                    revision=initial.sha,
                )
            }
            for i, e in zip(indices, entries):
                row = data.iloc[i]
                assert (e["record_id"], e["sha256"], e["split"], e["shape"]) == (
                    row.record_id,
                    row.sha256,
                    row.split,
                    [33, int(row.length), 1280],
                )
                r = records[e["path"]]
                assert r.size == e["bytes"]
                assert (
                    (r.lfs.sha256 == e["file_sha256"])
                    if r.lfs
                    else (r.blob_id == e["git_blob_sha1"])
                )
            all_entries.extend(entries)
            print(
                f"Reused and verified {len(all_entries)}/{len(data)} proteins",
                flush=True,
            )
            continue

        if model is None:
            model, tokenizer = load_model("cuda", cache)
        entries, files = [], {}
        for i in indices:
            row = data.iloc[i]
            length = int(row.length)
            inputs = {
                k: v.to("cuda") for k, v in tokenize(tokenizer, [row.sequence]).items()
            }
            context = (
                nullcontext()
                if precision == "float32"
                else torch.autocast("cuda", dtype=getattr(torch, precision))
            )
            with torch.inference_mode(), context:
                output = model(**inputs, output_hidden_states=True, return_dict=True)
                assert len(output.hidden_states) == 34
                values = (
                    torch.stack(
                        [
                            h[0, 1 : length + 1, :].to(save_dtype)
                            for h in output.hidden_states[1:]
                        ]
                    )
                    .cpu()
                    .contiguous()
                )
            del output, inputs
            assert (
                tuple(values.shape) == (33, length, 1280)
                and torch.isfinite(values).all().item()
            )
            means = values.float().mean(dim=1).numpy()
            ref = reference[i]
            relative = np.linalg.norm(means - ref, axis=1) / np.maximum(
                np.linalg.norm(ref, axis=1), 1e-12
            )
            cosine = (means * ref).sum(axis=1) / np.maximum(
                np.linalg.norm(means, axis=1) * np.linalg.norm(ref, axis=1), 1e-12
            )
            assert relative.max() <= 0.05 and cosine.min() >= 0.999, (
                f"Mean alignment check failed for {row.record_id}"
            )
            path = scratch / f"{row.sha256}.safetensors"
            save_file(
                {"activations": values},
                str(path),
                metadata={
                    "record_id": row.record_id,
                    "sequence_sha256": row.sha256,
                    "split": row.split,
                    "axes": "layer,residue,feature",
                    "layers": "1-33",
                    "residue_indexing": "0 = original residue 1",
                    "config_sha256": config_hash,
                },
            )
            destination = f"{prefix}/proteins/{path.name}"
            entries.append(
                dict(
                    record_id=row.record_id,
                    sha256=row.sha256,
                    split=row.split,
                    shape=[33, length, 1280],
                    path=destination,
                    bytes=path.stat().st_size,
                    file_sha256=sha(path),
                    git_blob_sha1=git_blob_sha(path),
                    max_relative_mean_difference=float(relative.max()),
                    min_mean_cosine=float(cosine.min()),
                )
            )
            files[destination] = path
            del values
        manifest = scratch / "batch.json"
        write_json(manifest, dict(config_sha256=config_hash, proteins=entries))
        files[batch_name] = manifest
        commit(files)  # Atomic: eight protein files plus their checksum/ID manifest.
        all_entries.extend(entries)
        for path in files.values():
            path.unlink()  # Remove only this batch's temporary files, after successful HF commit.
        print(f"Uploaded {len(all_entries)}/{len(data)} full protein maps", flush=True)

    assert len(all_entries) == 3052
    assert len({e["record_id"] for e in all_entries}) == 3052
    assert {e["record_id"] for e in all_entries} == set(data.record_id)
    write_json(
        scratch / "completed.json",
        dict(config_sha256=config_hash, proteins=all_entries),
    )
    commit({prefix + "/completed.json": scratch / "completed.json"})
    print(
        f"COMPLETE: {repo}/{prefix} — 3052 train/validation proteins, all 33 layers; no test proteins.",
        flush=True,
    )


if __name__ == "__main__":
    main()
