"""Frozen, resumable all-layer ESM extraction. Model access uses HF_TOKEN_PLM only."""

import argparse
from contextlib import nullcontext
import importlib.metadata
import json
import os
from pathlib import Path
import time
import numpy as np
from common import (
    HERE,
    MODEL,
    REVISION,
    PROTOCOL,
    WEIGHTS_SHA256,
    read_extraction,
    sha,
    token,
    write_json,
)

KEYS = {"activations", "record_id", "sha256", "split"}


def versions():
    return {
        p: importlib.metadata.version(p)
        for p in [
            "torch",
            "transformers",
            "huggingface-hub",
            "numpy",
            "pandas",
            "safetensors",
            "tokenizers",
            "psutil",
        ]
    }


def device_name(requested):
    import torch

    if requested != "auto":
        return requested
    return (
        "cuda"
        if torch.cuda.is_available()
        else ("mps" if torch.backends.mps.is_available() else "cpu")
    )


def memory(device):
    import torch
    import psutil

    out = {
        "device": device,
        "RAM_available_GB": round(psutil.virtual_memory().available / 1e9, 2),
    }
    if device == "cuda":
        out.update(
            GPU=torch.cuda.get_device_name(),
            allocated_GB=torch.cuda.memory_allocated() / 1e9,
            peak_GB=torch.cuda.max_memory_allocated() / 1e9,
            total_GB=torch.cuda.get_device_properties(0).total_memory / 1e9,
        )
    elif device == "mps":
        out["MPS_allocated_GB"] = torch.mps.current_allocated_memory() / 1e9
    return out


def load_model(device, cache):
    from transformers import AutoTokenizer, EsmModel
    from transformers.utils import logging as log
    from huggingface_hub import HfApi, hf_hub_download
    from huggingface_hub.utils import logging as hflog

    credential = token()
    log.set_verbosity_error()
    hflog.set_verbosity_error()
    # Never propagate network exception messages that might include request credentials.
    try:
        api = HfApi(token=credential)
        api.whoami()
        assert api.model_info(MODEL, revision=REVISION).sha == REVISION
        weights = hf_hub_download(
            MODEL,
            "model.safetensors",
            revision=REVISION,
            token=credential,
            cache_dir=cache,
        )
        assert sha(weights) == WEIGHTS_SHA256
        tokenizer = AutoTokenizer.from_pretrained(
            MODEL, revision=REVISION, token=credential, cache_dir=cache
        )
        model = EsmModel.from_pretrained(
            MODEL,
            revision=REVISION,
            token=credential,
            cache_dir=cache,
            use_safetensors=True,
            add_pooling_layer=False,
            attn_implementation="eager",
        )
    except Exception:
        raise RuntimeError(
            "Authenticated pinned ESM download failed; check HF_TOKEN_PLM permissions, network and cache. Details suppressed to protect credentials."
        ) from None
    assert model.config.num_hidden_layers == 33 and model.config.hidden_size == 1280
    model.requires_grad_(False)
    model.eval()
    model.to(device)
    return model, tokenizer


def tokenize(tokenizer, sequences):
    assert all(0 < len(s) <= 1022 for s in sequences)
    encoded = tokenizer(
        list(sequences), padding=True, truncation=False, return_tensors="pt"
    )
    for i, s in enumerate(sequences):
        count = int(encoded["attention_mask"][i].sum())
        assert count == len(s) + 2
        ids = encoded["input_ids"][i].tolist()
        assert (
            ids[0] == tokenizer.cls_token_id
            and ids[len(s) + 1] == tokenizer.eos_token_id
        )
        expected = tokenizer.convert_tokens_to_ids(list(s))
        assert (
            ids[1 : len(s) + 1] == expected and tokenizer.unk_token_id not in expected
        )
        assert all(t == tokenizer.pad_token_id for t in ids[len(s) + 2 :])
    return encoded


def pool_hidden(hidden, lengths):
    """Works for test tensors too; excludes BOS/EOS/padding by explicit positions."""
    import torch

    return torch.stack(
        [
            torch.stack(
                [h[i, 1 : int(n) + 1].float().mean(0) for i, n in enumerate(lengths)]
            )
            for h in hidden
        ],
        dim=1,
    )


def forward(model, tokenizer, sequences, device, precision):
    import torch

    encoded = {k: v.to(device) for k, v in tokenize(tokenizer, sequences).items()}
    context = (
        nullcontext()
        if precision == "float32"
        else torch.autocast(
            device_type=device,
            dtype={"float16": torch.float16, "bfloat16": torch.bfloat16}[precision],
        )
    )
    with torch.inference_mode(), context:
        result = model(**encoded, output_hidden_states=True, return_dict=True)
        assert len(result.hidden_states) == 34
        pooled = pool_hidden(result.hidden_states[1:], list(map(len, sequences)))
        assert tuple(pooled.shape) == (len(sequences), 33, 1280)
        output = pooled.cpu().numpy().astype(np.float32, copy=False)
    del result, pooled, encoded
    assert np.isfinite(output).all()
    return output


def plans(data, shard_size):
    order = np.lexsort((data.record_id.to_numpy(), data.length.to_numpy()))
    return [order[i : i + shard_size] for i in range(0, len(order), shard_size)]


def batch_indices(lengths, batch_size, max_tokens):
    out = []
    for i, n in enumerate(lengths):
        proposed = out + [i]
        if out and (
            len(proposed) > batch_size
            or (max(lengths[j] for j in proposed) + 2) * len(proposed) > max_tokens
        ):
            yield out
            out = []
        out.append(i)
    if out:
        yield out


def collate_records(items):
    return items


def run_frame(
    model, tokenizer, frame, device, precision, batch_size, max_tokens, workers
):
    from torch.utils.data import DataLoader

    records = [(i, s) for i, s in enumerate(frame.sequence)]
    batches = list(batch_indices(frame.length.tolist(), batch_size, max_tokens))
    loader = DataLoader(
        records, batch_sampler=batches, num_workers=workers, collate_fn=collate_records
    )
    result = np.empty((len(frame), 33, 1280), dtype=np.float32)
    for items in loader:
        indices, sequences = zip(*items)
        result[list(indices)] = forward(model, tokenizer, sequences, device, precision)
    return result


def verify_shard(folder, name, expected):
    file = Path(folder) / f"{name}.npz"
    meta = Path(folder) / f"{name}.json"
    info = json.loads(meta.read_text())
    assert sha(file) == info["sha256"]
    with np.load(file, allow_pickle=False) as a:
        assert set(a.files) == KEYS
        for col in ["record_id", "sha256", "split"]:
            assert np.array_equal(a[col], expected[col].to_numpy(dtype=str)), (
                f"{name}: {col} mismatch"
            )
        values = a["activations"]
        assert values.shape == (len(expected), 33, 1280) and values.dtype == np.float32
        assert np.isfinite(values).all()
    assert info["shape"] == [len(expected), 33, 1280]
    return values


def save_shard(folder, name, frame, values):
    folder = Path(folder)
    temp = folder / f"{name}.npz.tmp"
    file = folder / f"{name}.npz"
    assert (
        values.shape == (len(frame), 33, 1280)
        and values.dtype == np.float32
        and np.isfinite(values).all()
    )
    with open(temp, "wb") as f:
        np.savez(
            f,
            activations=values,
            **{
                col: frame[col].to_numpy(dtype=str)
                for col in ["record_id", "sha256", "split"]
            },
        )
    temp.replace(file)
    write_json(
        folder / f"{name}.json",
        dict(sha256=sha(file), shape=list(values.shape), dtype="float32"),
    )
    verify_shard(folder, name, frame)


def verify_config(folder, data):
    config = json.loads((Path(folder) / "extraction.json").read_text())
    assert config["model"] == MODEL and config["revision"] == REVISION
    assert config["layers"] == list(range(1, 34)) and config["width"] == 1280
    assert config["saved_dtype"] == "float32"
    assert config["extraction_dtype"] in {"float32", "float16", "bfloat16"}
    assert len(data) == 3052 and data.record_id.is_unique and data.sha256.is_unique
    assert data.split.value_counts().to_dict() == PROTOCOL["counts"]
    assert isinstance(config["shard_size"], int) and config["shard_size"] > 0
    assert (
        config["pooling"] == PROTOCOL["pooling"]
        and config["indexing"] == PROTOCOL["indexing"]
    )
    identity = data[["record_id", "sha256", "split"]].to_dict("records")
    assert config["proteins"] == identity
    return config


def load_all(folder, data=None):
    data = read_extraction() if data is None else data
    config = verify_config(folder, data)
    array = np.empty((3052, 33, 1280), dtype=np.float32)
    for i, idx in enumerate(plans(data, config["shard_size"])):
        array[idx] = verify_shard(folder, f"shard_{i:04}", data.iloc[idx])
    assert np.isfinite(array).all()
    return array


def pilot(model, tokenizer, data, args, device, precision):
    ordered = data.sort_values(["length", "record_id"])
    indices = set(
        ordered.iloc[np.linspace(0, len(data) - 1, 5).astype(int)].index.tolist()
    )
    indices.update(data[data.sequence.str.contains("X")].index)
    frame = data.loc[sorted(indices)].copy()
    # FP32 single-protein reference also checks mixed-batch padding invariance.
    ref = run_frame(model, tokenizer, frame, device, "float32", 1, args.max_tokens, 0)
    fp = run_frame(
        model,
        tokenizer,
        frame,
        device,
        "float32",
        args.batch_size,
        args.max_tokens,
        args.workers,
    )
    padding_error = float(np.max(np.abs(ref - fp)))
    assert np.allclose(ref, fp, atol=2e-4, rtol=2e-3), "Batch/padding invariance failed"
    report = dict(
        model=MODEL,
        revision=REVISION,
        proteins=frame.record_id.tolist(),
        lengths=frame.length.tolist(),
        X_proteins=int(frame.sequence.str.contains("X").sum()),
        padding_max_abs_error=padding_error,
        precision=precision,
        software=versions(),
        device=device,
        passed=True,
    )
    if precision != "float32":
        mixed = run_frame(
            model,
            tokenizer,
            frame,
            device,
            precision,
            args.batch_size,
            args.max_tokens,
            args.workers,
        )
        cosine = np.sum(ref * mixed, axis=2) / (
            np.linalg.norm(ref, axis=2) * np.linalg.norm(mixed, axis=2)
        )
        relative = np.linalg.norm(ref - mixed, axis=2) / np.maximum(
            np.linalg.norm(ref, axis=2), 1e-12
        )
        report.update(
            mixed_min_cosine=float(cosine.min()),
            mixed_max_relative_L2=float(relative.max()),
            mixed_max_abs_error=float(np.max(np.abs(ref - mixed))),
        )
        assert cosine.min() >= 0.999 and relative.max() <= 0.05, (
            "Mixed precision pilot failed; use float32, do not relax thresholds"
        )
    else:
        report["mixed_precision_comparison"] = (
            "not performed; use CUDA float16 or bfloat16 pilot to compare"
        )
    args.output.mkdir(parents=True, exist_ok=True)
    write_json(args.output / "pilot.json", report)
    print(json.dumps(report, indent=2))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("action", choices=["pilot", "extract", "validate"])
    p.add_argument(
        "--manifest", type=Path, default=HERE / "data/extraction_manifest.tsv"
    )
    p.add_argument("--output", type=Path, default=HERE / "activations")
    p.add_argument("--device", choices=["auto", "cuda", "mps", "cpu"], default="auto")
    p.add_argument(
        "--precision",
        choices=["auto", "float32", "float16", "bfloat16"],
        default="auto",
    )
    p.add_argument("--batch-size", type=int, default=2)
    p.add_argument("--max-tokens", type=int, default=2048)
    p.add_argument("--shard-size", type=int, default=64)
    p.add_argument("--workers", type=int, default=0)
    p.add_argument(
        "--hf-repo",
        help="Existing private dataset; upload every shard before continuing",
    )
    p.add_argument(
        "--cache",
        default=os.environ.get("NESG_ESM_CACHE", str(Path.home() / ".cache/nesg-esm")),
    )
    args = p.parse_args()
    data = read_extraction(args.manifest)
    assert (
        args.batch_size > 0
        and args.shard_size > 0
        and args.workers >= 0
        and args.max_tokens >= int(data.length.max()) + 2
    )
    if args.action == "validate":
        array = load_all(args.output, data)
        print("Verified finite aligned activations", array.shape)
        return
    if args.hf_repo:
        from hf_storage import private_api

        private_api(args.hf_repo)
    device = device_name(args.device)
    precision = (
        ("float16" if device == "cuda" else "float32")
        if args.precision == "auto"
        else args.precision
    )
    assert precision == "float32" or device == "cuda", (
        "Mixed precision is supported here only on CUDA"
    )
    if (args.output / "completed.json").is_file():
        from hf_storage import verified_files

        config = verify_config(args.output, data)
        assert config["extraction_dtype"] == precision, (
            "Completed cache uses another precision."
        )
        verified_files(args.output)
        print(
            f"SKIP {args.action}: all 3052 pooled activations already complete and verified; no ESM model loaded.",
            flush=True,
        )
        return
    pilot_file = args.output / "pilot.json"
    if args.action == "pilot" and pilot_file.is_file():
        check = json.loads(pilot_file.read_text())
        if (
            check.get("passed")
            and check.get("model") == MODEL
            and check.get("revision") == REVISION
            and check.get("precision") == precision
            and check.get("software") == versions()
            and check.get("device") == device
        ):
            print(
                "SKIP pilot: matching successful correctness pilot already saved.",
                flush=True,
            )
            return
    print(json.dumps(memory(device)), flush=True)
    model, tokenizer = load_model(device, args.cache)
    if args.action == "pilot":
        pilot(model, tokenizer, data, args, device, precision)
        return
    pilot_file = args.output / "pilot.json"
    assert pilot_file.exists(), "Run the real ESM correctness pilot first"
    check = json.loads(pilot_file.read_text())
    assert (
        check["passed"]
        and check["model"] == MODEL
        and check["revision"] == REVISION
        and check["precision"] == precision
        and check["software"] == versions()
        and check["device"] == device
    )
    args.output.mkdir(parents=True, exist_ok=True)
    config = dict(
        model=MODEL,
        revision=REVISION,
        layers=list(range(1, 34)),
        width=1280,
        pooling=PROTOCOL["pooling"],
        indexing=PROTOCOL["indexing"],
        extraction_dtype=precision,
        saved_dtype="float32",
        software=versions(),
        shard_size=args.shard_size,
        manifest_sha256=sha(args.manifest),
        proteins=data[["record_id", "sha256", "split"]].to_dict("records"),
    )
    path = args.output / "extraction.json"
    if path.exists():
        assert json.loads(path.read_text()) == config, (
            "Incompatible existing extraction; use a separate output directory"
        )
    else:
        write_json(path, config)
    for i, idx in enumerate(plans(data, args.shard_size)):
        name = f"shard_{i:04}"
        frame = data.iloc[idx]
        if (args.output / f"{name}.json").exists():
            verify_shard(args.output, name, frame)
            if args.hf_repo:
                from hf_storage import upload_checkpoint

                upload_checkpoint(args.hf_repo, args.output, name)
            print(name, "verified; reused", flush=True)
            continue
        start = time.monotonic()
        values = run_frame(
            model,
            tokenizer,
            frame,
            device,
            precision,
            args.batch_size,
            args.max_tokens,
            args.workers,
        )
        save_shard(args.output, name, frame, values)
        if args.hf_repo:
            from hf_storage import upload_checkpoint

            upload_checkpoint(args.hf_repo, args.output, name)
        print(
            name,
            "completed",
            len(frame),
            "proteins",
            round(time.monotonic() - start, 1),
            "seconds",
            memory(device),
            flush=True,
        )
    del model
    full = load_all(args.output, data)
    files = {
        f.name: sha(f)
        for f in sorted(args.output.glob("shard_*.*"))
        if f.suffix in [".npz", ".json"]
    }
    files["extraction.json"] = sha(path)
    write_json(
        args.output / "completed.json",
        dict(shape=list(full.shape), files=files, model=MODEL, revision=REVISION),
    )
    if args.hf_repo:
        from hf_storage import upload

        upload(args.hf_repo, args.output)
    print("Complete:", full.shape, flush=True)


if __name__ == "__main__":
    main()
