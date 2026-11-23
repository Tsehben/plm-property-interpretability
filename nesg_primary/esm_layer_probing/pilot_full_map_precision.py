"""Label-blind NESG float16 storage/InterProt compatibility pilot; no probe fitting.

Run in the extraction environment after pooled extraction. Tests actual safetensors
round trips against unrounded values from the SAME ESM forward pass, isolating
storage rounding from inference precision. Only a compact report is uploaded.
"""

import argparse
import importlib.util
import json
from pathlib import Path
from contextlib import nullcontext
from tempfile import TemporaryDirectory

import numpy as np
import torch

from common import (
    HERE,
    MODEL,
    REVISION,
    INTERPROT_REVISION,
    SAE_LAYERS,
    read_extraction,
    sha,
    token,
    write_json,
)

SAE_REPO = "liambai/InterProt-ESM2-SAEs"
SOURCE_SHA = "1c5d4bc913d0f6896a2143e0fcc5df013d8aad390f8a5dcc6f6912117a7672fb"
# Engineering tolerances specified before examining NESG SAE results. These are
# not established biological-equivalence thresholds or guarantees for every feature.
LIMITS = dict(
    raw_relative_l2_max=0.001,
    mean_active_jaccard_min=0.99,
    latent_relative_l2_max=0.02,
    reconstruction_relative_l2_max=0.01,
)

# Previous reviewed implementation used for the user's successful Colab pilot.
PRIOR_PILOT_SHA = "11257c02be0ecd1e2327a5002ae5f3c5a8e9ec2d1a37097f9df200b843fe3008"


def report_matches(report, sample, config, completion_sha):
    expected = {
        (r.record_id, r.sha256, r.split, int(r.length), layer)
        for r in sample.itertuples()
        for layer in SAE_LAYERS
    }
    try:
        observed = {
            (r["record_id"], r["sequence_sha256"], r["split"], r["length"], r["layer"])
            for r in report["results"]
        }
        return (
            report["pilot_script_sha256"] in {PRIOR_PILOT_SHA, sha(__file__)}
            and report["model"] == MODEL
            and report["model_revision"] == REVISION
            and report["sae_repo"] == SAE_REPO
            and report["sae_revision"] == INTERPROT_REVISION
            and report["sae_source_sha256"] == SOURCE_SHA
            and report["sae_topk"] == 128
            and report["inference_precision"] == config["extraction_dtype"]
            and report["software"] == config["software"]
            and report["manifest_sha256"] == sha(HERE / "data/extraction_manifest.tsv")
            and report["pooled_completion_sha256"] == completion_sha
            and report["engineering_tolerances"] == LIMITS
            and report["proteins"] == len(sample)
            and observed == expected
            and len(report["results"]) == len(expected)
        )
    except (KeyError, TypeError):
        return False


def reuse_remote_report(api, info, args, sample, config):
    """Read tiny immutable reports only; never download full-map tensors here."""
    from huggingface_hub import hf_hub_download

    try:
        names = api.list_repo_files(args.repo, repo_type="dataset", revision=info.sha)
        for name in sorted(names):
            if not name.startswith("esm_stage1/precision_pilot/") or not name.endswith(
                ".json"
            ):
                continue
            path = Path(
                hf_hub_download(
                    args.repo,
                    name,
                    repo_type="dataset",
                    revision=info.sha,
                    token=token(),
                )
            )
            assert sha(path) == Path(name).stem, "Precision report checksum differs."
            report = json.loads(path.read_text())
            if not report_matches(
                report, sample, config, sha(args.pooled_dir / "completed.json")
            ):
                continue
            write_json(args.output, report)
            if report["status"] != "PASS" or not all(
                r["within_engineering_tolerances"] for r in report["results"]
            ):
                raise SystemExit(
                    "Existing matching precision pilot requires review; not rerunning it automatically."
                )
            print(
                f"SKIP precision pilot: verified previous PASS for {len(sample)} proteins / {len(SAE_LAYERS)} SAE layers in private HF. No models downloaded or run.",
                flush=True,
            )
            return True
    except Exception:
        raise RuntimeError(
            "Could not verify saved precision pilot; check HF access/network. No pilot was rerun."
        ) from None
    return False


def select_sample(data):
    """Four length ranks per split, plus every X-containing sequence (<=11)."""
    ids = set(data.loc[data.sequence.str.contains("X"), "record_id"])
    for split in ("train", "validation"):
        group = data[data.split == split].sort_values(["length", "record_id"])
        positions = np.linspace(0, len(group) - 1, 4).round().astype(int)
        ids.update(group.iloc[positions].record_id)
    return data[data.record_id.isin(ids)].sort_values(["length", "record_id"])


def relative_l2(a, b):
    return float(
        torch.linalg.vector_norm(a - b) / torch.linalg.vector_norm(a).clamp_min(1e-12)
    )


def compare(sae, original, rounded):
    """Both inputs and SAE arithmetic are float32; positive activity defines support."""
    assert original.dtype == rounded.dtype == torch.float32
    assert original.shape == rounded.shape
    with torch.inference_mode():
        z0, z1 = sae.get_acts(original), sae.get_acts(rounded)
        r0, r1 = sae.forward_val(original), sae.forward_val(rounded)
    assert all(
        torch.isfinite(t).all().item() for t in (original, rounded, z0, z1, r0, r1)
    )
    s0, s1 = z0 > 0, z1 > 0
    union = (s0 | s1).sum(-1)
    intersection = (s0 & s1).sum(-1)
    jaccard = torch.where(union > 0, intersection.float() / union.clamp_min(1), 1.0)
    metrics = dict(
        raw_relative_l2=relative_l2(original, rounded),
        max_absolute_rounding=float((original - rounded).abs().max()),
        mean_active_jaccard=float(jaccard.mean()),
        minimum_residue_active_jaccard=float(jaccard.min()),
        changed_support_residue_fraction=float((s0 != s1).any(-1).float().mean()),
        latent_relative_l2=relative_l2(z0, z1),
        reconstruction_relative_l2=relative_l2(r0, r1),
        original_reconstruction_relative_error=relative_l2(original, r0),
        rounded_reconstruction_relative_error=relative_l2(original, r1),
    )
    metrics["within_engineering_tolerances"] = (
        metrics["raw_relative_l2"] <= LIMITS["raw_relative_l2_max"]
        and metrics["mean_active_jaccard"] >= LIMITS["mean_active_jaccard_min"]
        and metrics["latent_relative_l2"] <= LIMITS["latent_relative_l2_max"]
        and metrics["reconstruction_relative_l2"]
        <= LIMITS["reconstruction_relative_l2_max"]
    )
    return metrics


def roundtrip(values, directory, metadata):
    """Use the same [33,L,1280] tensor/file layout as save_full_maps.py."""
    from safetensors import safe_open
    from safetensors.torch import save_file, load_file

    assert (
        values.dtype == torch.float32
        and values.shape[0] == 33
        and values.shape[2] == 1280
    )
    arrays = []
    for dtype in (torch.float32, torch.float16):
        path = directory / f"{dtype}.safetensors"
        save_file(
            {"activations": values.to(dtype).contiguous()}, str(path), metadata=metadata
        )
        with safe_open(str(path), framework="pt") as handle:
            assert handle.metadata() == metadata
        saved = load_file(str(path))["activations"]
        assert saved.dtype == dtype and saved.shape == values.shape
        assert torch.isfinite(saved).all().item()
        arrays.append(saved.float())
        path.unlink()
    assert torch.equal(arrays[0], values)
    return arrays


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pooled-dir", required=True, type=Path)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--sae-source", required=True, type=Path)
    parser.add_argument(
        "--repo",
        required=True,
        help="Existing private destination for compact report only",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("/content/nesg_full_map_precision_pilot.json"),
    )
    args = parser.parse_args()
    from extract import load_model, tokenize, verify_config, versions, load_all
    from hf_storage import private_api
    from huggingface_hub import hf_hub_download
    from safetensors.torch import load_file

    assert torch.cuda.is_available(), "Use the existing GPU extraction runtime."
    assert (args.pooled_dir / "completed.json").is_file(), (
        "Complete pooled extraction first."
    )
    data = read_extraction()
    config = verify_config(args.pooled_dir, data)
    assert versions() == config["software"], "Use the original extraction environment."
    reference = load_all(args.pooled_dir, data)
    sample = select_sample(data)
    assert set(sample.split) == {"train", "validation"}
    assert sample.sequence.str.contains("X").sum() == 3
    assert sha(args.sae_source) == SOURCE_SHA, (
        "SAE implementation differs from reviewed reference."
    )
    spec = importlib.util.spec_from_file_location(
        "pinned_interprot_sae", args.sae_source
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    api, info = private_api(args.repo)
    if reuse_remote_report(api, info, args, sample, config):
        return
    # Disable TF32 for this precision comparison; the ESM inference mode below
    # otherwise matches the full-map exporter. No independent full-fp32 ESM run.
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    checkpoints, saes = {}, {}
    for layer in SAE_LAYERS:
        filename = f"esm2_plm1280_l{layer}_sae4096.safetensors"
        try:
            path = hf_hub_download(
                SAE_REPO,
                filename,
                revision=INTERPROT_REVISION,
                token=token(),
                cache_dir=str(args.cache),
            )
        except Exception:
            raise RuntimeError(
                f"Pinned SAE download failed for layer {layer}; check access/network."
            ) from None
        sae = module.SparseAutoencoder(1280, 4096, k=128)
        sae.load_state_dict(load_file(path), strict=True)
        saes[layer] = sae.eval().requires_grad_(False).to("cuda", dtype=torch.float32)
        checkpoints[str(layer)] = dict(filename=filename, sha256=sha(path))
    model, tokenizer = load_model("cuda", str(args.cache))
    precision = config["extraction_dtype"]
    assert precision in ("float32", "float16", "bfloat16")
    results = []
    print(
        f"Testing {len(sample)} proteins at {len(SAE_LAYERS)} SAE layers; no bulk-map upload.",
        flush=True,
    )
    for number, row in enumerate(sample.itertuples(), 1):
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
            values = torch.stack(
                [
                    h[0, 1 : row.length + 1].float().cpu()
                    for h in output.hidden_states[1:]
                ]
            ).contiguous()
        del output, inputs
        metadata = dict(
            record_id=row.record_id,
            sequence_sha256=row.sha256,
            split=row.split,
            axes="layer,residue,feature",
            layers="1-33",
            residue_indexing="0 = original residue 1",
        )
        with TemporaryDirectory(prefix="nesg_precision_") as tmp:
            original, rounded = roundtrip(values, Path(tmp), metadata)
        means = rounded.mean(1).numpy()
        ref = reference[data.index.get_loc(row.Index)]
        cosine = (means * ref).sum(1) / np.maximum(
            np.linalg.norm(means, axis=1) * np.linalg.norm(ref, axis=1), 1e-12
        )
        relative = np.linalg.norm(means - ref, axis=1) / np.maximum(
            np.linalg.norm(ref, axis=1), 1e-12
        )
        assert cosine.min() >= 0.999 and relative.max() <= 0.05, (
            "Pooled alignment check failed."
        )
        for layer, sae in saes.items():
            # Stored layer index 7 is ESM layer 8, NOT layer index 8.
            metrics = compare(
                sae, original[layer - 1].cuda(), rounded[layer - 1].cuda()
            )
            results.append(
                dict(
                    record_id=row.record_id,
                    sequence_sha256=row.sha256,
                    split=row.split,
                    length=int(row.length),
                    layer=layer,
                    min_pooled_cosine=float(cosine.min()),
                    **metrics,
                )
            )
        print(
            f"Checked {number}/{len(sample)}: {row.record_id}, length={row.length}",
            flush=True,
        )
    passed = all(r["within_engineering_tolerances"] for r in results)
    report = dict(
        status="PASS" if passed else "REVIEW_REQUIRED",
        proteins=len(sample),
        selection="Four evenly spaced length ranks within each development split plus all X sequences; no labels read",
        model=MODEL,
        model_revision=REVISION,
        inference_precision=precision,
        comparison="Same forward pass: float32 storage versus float16 storage reloaded as float32; SAE arithmetic float32",
        sae_repo=SAE_REPO,
        sae_revision=INTERPROT_REVISION,
        sae_source_sha256=SOURCE_SHA,
        sae_checkpoints=checkpoints,
        sae_topk=128,
        software=versions(),
        manifest_sha256=sha(HERE / "data/extraction_manifest.tsv"),
        pilot_script_sha256=sha(__file__),
        pooled_completion_sha256=sha(args.pooled_dir / "completed.json"),
        engineering_tolerances=LIMITS,
        results=results,
        limitations="Small numerical pilot, not biological validation. Does not establish float32-inference equivalence, individual-feature invariance, or intervention validity. No labels, probes, feature associations or test data used.",
    )
    write_json(args.output, report)
    print(
        f"{report['status']}: worst protein/layer mean active-set Jaccard={min(r['mean_active_jaccard'] for r in results):.6f}; "
        f"max latent relative L2={max(r['latent_relative_l2'] for r in results):.6f}",
        flush=True,
    )
    destination = f"esm_stage1/precision_pilot/{sha(args.output)}.json"
    try:
        current = api.repo_info(args.repo, repo_type="dataset")
        assert current.private is True
        api.upload_file(
            path_or_fileobj=str(args.output),
            path_in_repo=destination,
            repo_id=args.repo,
            repo_type="dataset",
            parent_commit=current.sha,
            commit_message="Record label-blind full-map storage precision pilot",
        )
    except Exception:
        raise RuntimeError(
            f"Pilot report saved at {args.output}, but private HF upload failed."
        ) from None
    print(f"Report saved in private HF: {destination}", flush=True)
    if not passed:
        raise SystemExit(
            "Review precision differences before launching the bulk float16 export."
        )


if __name__ == "__main__":
    main()
