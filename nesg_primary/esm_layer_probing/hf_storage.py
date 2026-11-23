"""Private storage; repository creation requires explicit init --create-private."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
from common import HERE, read_extraction, sha, token, write_json
from extract import load_all, plans, verify_config, verify_shard


def status_code(error):
    # Never serialize exception messages, request headers or response bodies.
    status = getattr(getattr(error, "response", None), "status_code", None)
    return status if isinstance(status, int) else None


def authenticated_api(repo):
    from huggingface_hub import HfApi
    from huggingface_hub.utils import logging, validate_repo_id

    if not repo or "/" not in repo:
        raise ValueError("Explicit owner/private-dataset destination required")
    try:
        validate_repo_id(repo)
    except Exception:
        raise ValueError("Invalid dataset ID; use organization/dataset-name.") from None
    logging.set_verbosity_error()
    api = HfApi(token=token())
    try:
        api.whoami()
    except Exception as error:
        status = status_code(error)
        if status in (401, 403):
            raise RuntimeError(
                "HF_TOKEN_PLM authentication failed (HTTP 401/403). Check that the Colab Secret contains a valid, non-revoked Hugging Face token."
            ) from None
        raise RuntimeError(
            "Could not verify HF_TOKEN_PLM with Hugging Face. Check connectivity and retry; no repository was modified."
        ) from None
    return api


def repository_error(error):
    status = status_code(error)
    if status in (401, 404):
        return RuntimeError(
            "Token authentication succeeded, but the configured dataset is missing or hidden from this token (HTTP 401/404). Check the exact repository name and organization/token access. If it has not been created, explicitly run hf_storage.py init --repo ORGANIZATION/DATASET --create-private."
        )
    if status == 403:
        return RuntimeError(
            "Token authentication succeeded, but dataset access is forbidden (HTTP 403). Grant this token access to the private dataset and check PLM-Interp organization membership/resource-group permissions."
        )
    return RuntimeError(
        "Token authentication succeeded, but dataset lookup failed. Check connectivity/service availability and retry; no repository was modified."
    )


def require_private(info):
    if info.private is not True:
        raise ValueError(
            "Refusing non-private destination. Choose a private dataset; this command will not change an existing repository visibility."
        )


def private_api(repo):
    api = authenticated_api(repo)
    try:
        info = api.repo_info(repo, repo_type="dataset")
    except Exception as error:
        raise repository_error(error) from None
    require_private(info)
    return api, info


def initialize_private(repo, create_private=False):
    """Only this explicitly requested command can create a new private dataset."""
    api = authenticated_api(repo)
    try:
        info = api.repo_info(repo, repo_type="dataset")
    except Exception as error:
        if not create_private or status_code(error) not in (401, 404):
            raise repository_error(error) from None
        try:
            api.create_repo(
                repo_id=repo, repo_type="dataset", private=True, exist_ok=False
            )
        except Exception as creation_error:
            status = status_code(creation_error)
            if status == 409:
                raise RuntimeError(
                    "Private dataset creation failed (HTTP 409): a repository already exists under this name. No overwrite was performed; check access and retry the check command."
                ) from None
            if status in (401, 403):
                advice = "Token authentication succeeded, but creation was denied. Check token write scope for PLM-Interp and organization contributor/write/admin access."
            elif status == 429:
                advice = "Hugging Face rate-limited the request. Wait and retry; changing runtime packages will not help."
            elif status is not None and status >= 500:
                advice = "Hugging Face returned a service error. Retry later."
            elif status in (400, 404, 422):
                advice = "Check the exact organization/name and try the Hugging Face New Dataset form for any naming, policy or account requirements."
            else:
                advice = "Check connectivity and retry. The server outcome is unknown; check whether the repository exists before retrying creation."
            raise RuntimeError(
                f"Private dataset creation failed (HTTP {status if status is not None else 'status unavailable'}). {advice} No public repository or upload was requested."
            ) from None
        try:
            info = api.repo_info(repo, repo_type="dataset")
        except Exception as verification_error:
            status = status_code(verification_error)
            raise RuntimeError(
                f"Private dataset creation succeeded, but subsequent access verification failed (HTTP {status if status is not None else 'status unavailable'}). "
                "The repository may already exist: check it in the HF website and grant this token access to that dataset, then run check. No upload was attempted."
            ) from None
    require_private(info)
    print("Token authentication and private dataset access verified:", repo, flush=True)
    return api, info


def diagnose_access(repo):
    """Read-only diagnosis: print only allowlisted account/permission fields."""
    api = authenticated_api(repo)
    try:
        identity = api.whoami()
    except Exception:
        raise RuntimeError(
            "Could not read account permissions; retry. Server details suppressed."
        ) from None
    organization = repo.split("/")[0]
    account = identity.get("name", "not reported")
    access_token = identity.get("auth", {}).get("accessToken", {})
    role = access_token.get("role")
    token_type = (
        role
        if role in {"read", "write", "fineGrained", "fine-grained"}
        else "not reported"
    )
    membership = next(
        (org for org in identity.get("orgs", []) if org.get("name") == organization), {}
    )
    role = membership.get("roleInOrg")
    org_role = (
        role
        if role in {"admin", "write", "contributor", "read", "no_access"}
        else "not reported (does not establish non-membership)"
    )
    print(
        json.dumps(
            {
                "authenticated_account": account,
                "token_type": token_type,
                "organization": organization,
                "organization_role": org_role,
            },
            indent=2,
        ),
        flush=True,
    )
    try:
        info = api.repo_info(repo, repo_type="dataset")
    except Exception as error:
        print(
            f"Dataset lookup HTTP status: {status_code(error)}. {repository_error(error)}",
            flush=True,
        )
        return
    print("Dataset accessible; private:", info.private is True, flush=True)
    print(
        "This read-only diagnosis does not verify write permission and creates/uploads nothing.",
        flush=True,
    )


def verified_files(folder):
    folder = Path(folder)
    # Revalidate every array and identifier before allowing any upload.
    load_all(folder, read_extraction())
    complete = json.loads((folder / "completed.json").read_text())
    assert complete["shape"] == [3052, 33, 1280]
    files = {}
    for name, digest in complete["files"].items():
        assert Path(name).name == name
        assert name == "extraction.json" or (
            name.startswith("shard_") and Path(name).suffix in [".npz", ".json"]
        )
        path = folder / name
        assert sha(path) == digest
        files[name] = path
    safe_config(folder)
    files["completed.json"] = folder / "completed.json"
    return files


def safe_config(folder):
    """Only this schema may leave the runtime; labels/sequences are excluded."""
    folder = Path(folder)
    config = json.loads((folder / "extraction.json").read_text())
    allowed = {
        "model",
        "revision",
        "layers",
        "width",
        "pooling",
        "indexing",
        "extraction_dtype",
        "saved_dtype",
        "software",
        "shard_size",
        "manifest_sha256",
        "proteins",
    }
    assert set(config) == allowed
    assert all(
        set(row) == {"record_id", "sha256", "split"} for row in config["proteins"]
    )
    return config


def git_blob_sha(path):
    raw = path.read_bytes()
    return hashlib.sha1(f"blob {len(raw)}\0".encode() + raw).hexdigest()


def commit_files(repo, files):
    """Upload new files atomically; never replace incompatible research artifacts."""
    from huggingface_hub import CommitOperationAdd

    api, info = private_api(repo)
    try:
        existing = {
            r.path: r
            for r in (
                api.get_paths_info(
                    repo,
                    [f"esm_stage1/{name}" for name in files],
                    repo_type="dataset",
                    revision=info.sha,
                )
                if info.sha
                else []
            )
        }
        changed = []
        for name, path in files.items():
            remote = existing.get("esm_stage1/" + name)
            if remote:
                lfs = getattr(remote, "lfs", None)
                same = (
                    lfs.sha256 == sha(path)
                    if lfs
                    else remote.blob_id == git_blob_sha(path)
                )
                if not same:
                    raise ValueError(
                        "Destination contains a different artifact; use a different private repository."
                    )
            else:
                changed.append(
                    CommitOperationAdd(
                        path_in_repo="esm_stage1/" + name, path_or_fileobj=str(path)
                    )
                )
        assert api.repo_info(repo, repo_type="dataset").private is True
        commit = (
            api.create_commit(
                repo_id=repo,
                repo_type="dataset",
                operations=changed,
                parent_commit=info.sha,
                commit_message="Checkpoint verified NESG development ESM activations",
            )
            if changed
            else None
        )
    except ValueError:
        # No remote exception text is propagated, even for ValueError.
        raise RuntimeError(
            "Private upload refused: an existing artifact differs. Use a fresh private repository; no overwrite was performed."
        ) from None
    except Exception:
        raise RuntimeError(
            "Private upload failed; check write permission/network and retry. Extraction stops until the checkpoint is durable."
        ) from None
    return dict(
        repository=repo,
        private_verified=True,
        uploaded_files=len(changed),
        commit=commit.oid if commit else info.sha,
        files={name: sha(path) for name, path in files.items()},
    )


def upload(repo, folder):
    receipt = commit_files(repo, verified_files(folder))
    write_json(Path(folder) / "upload_receipt.json", receipt)
    print(
        "Private repository verified; uploaded changed files:",
        receipt["uploaded_files"],
    )


def upload_checkpoint(repo, folder, name):
    """Called after every completed shard, including reused local shards."""
    folder = Path(folder)
    data = read_extraction()
    config = verify_config(folder, data)
    safe_config(folder)
    expected = {
        f"shard_{i:04}": data.iloc[idx]
        for i, idx in enumerate(plans(data, config["shard_size"]))
    }
    assert name in expected
    verify_shard(folder, name, expected[name])
    files = {
        "extraction.json": folder / "extraction.json",
        name + ".npz": folder / (name + ".npz"),
        name + ".json": folder / (name + ".json"),
    }
    receipt = commit_files(repo, files)
    write_json(folder / "upload_receipt.json", receipt)
    print(name, "checkpoint durable in private HF repository", flush=True)


def restore(repo, folder):
    """Recover partial checkpoints from one pinned remote commit, without model inference."""
    from huggingface_hub import hf_hub_download

    api, info = private_api(repo)
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    credential = token()
    if not info.sha:
        print(
            "Private destination verified; empty repository, no checkpoints to restore."
        )
        return
    # A reopened notebook commonly still has the complete local cache. Verify it
    # and its remote completion hash before doing any listing or downloading.
    if (folder / "completed.json").is_file():
        try:
            verified_files(folder)
        except (AssertionError, OSError, ValueError, KeyError):
            print("Local cache is incomplete or damaged; checking HF for repair.")
        else:
            try:
                entries = api.get_paths_info(
                    repo,
                    ["esm_stage1/completed.json"],
                    repo_type="dataset",
                    revision=info.sha,
                )
                matching = any(
                    e.path == "esm_stage1/completed.json"
                    and (
                        e.lfs.sha256 == sha(folder / "completed.json")
                        if e.lfs
                        else e.blob_id == git_blob_sha(folder / "completed.json")
                    )
                    for e in entries
                )
            except Exception:
                raise RuntimeError(
                    "Cannot verify the private HF completion record; retry when connected."
                ) from None
            if matching:
                print(
                    "SKIP restore: all 3052 local pooled activations verified against private HF; nothing downloaded."
                )
                return
    try:
        names = set(api.list_repo_files(repo, repo_type="dataset", revision=info.sha))
        if "esm_stage1/extraction.json" not in names:
            assert not any(n.startswith("esm_stage1/") for n in names), (
                "Incomplete remote metadata"
            )
            print("Private destination verified; no previous extraction checkpoint.")
            return

        def fetch(name):
            return Path(
                hf_hub_download(
                    repo,
                    "esm_stage1/" + name,
                    repo_type="dataset",
                    revision=info.sha,
                    token=credential,
                )
            )

        if "esm_stage1/completed.json" in names:
            # One small index replaces a separate HTTP request for every sidecar.
            index = fetch("completed.json")
            complete = json.loads(index.read_text())
            assert complete["shape"] == [3052, 33, 1280]
            downloaded_count = 0
            for name, digest in complete["files"].items():
                assert Path(name).name == name and (
                    name == "extraction.json"
                    or (
                        name.startswith("shard_")
                        and Path(name).suffix in [".npz", ".json"]
                    )
                )
                target = folder / name
                if target.is_file() and sha(target) == digest:
                    continue
                cached = fetch(name)
                assert sha(cached) == digest
                shutil.copyfile(cached, target)
                downloaded_count += 1
            shutil.copyfile(index, folder / "completed.json")
            verified_files(folder)
            print(
                f"Pooled cache ready: downloaded {downloaded_count} missing/damaged files; reused the rest. No ESM computation."
            )
            return

        remote_config = fetch("extraction.json")
        target = folder / "extraction.json"
        if target.exists():
            assert json.loads(target.read_text()) == json.loads(
                remote_config.read_text()
            )
        else:
            shutil.copyfile(remote_config, target)
        data = read_extraction()
        config = verify_config(folder, data)
        safe_config(folder)
        restored = 0
        for i, idx in enumerate(plans(data, config["shard_size"])):
            name = f"shard_{i:04}"
            present = [
                f"esm_stage1/{name}{suffix}" in names for suffix in [".npz", ".json"]
            ]
            assert present[0] == present[1], (
                "Remote shard must be an atomic array/metadata pair"
            )
            if not all(present):
                continue
            meta = fetch(name + ".json")
            digest = json.loads(meta.read_text())["sha256"]
            array = folder / (name + ".npz")
            if not array.exists() or sha(array) != digest:
                downloaded = fetch(name + ".npz")
                assert sha(downloaded) == digest
                shutil.copyfile(downloaded, array)
            shutil.copyfile(meta, folder / (name + ".json"))
            verify_shard(folder, name, data.iloc[idx])
            restored += 1
        if "esm_stage1/completed.json" in names:
            shutil.copyfile(fetch("completed.json"), folder / "completed.json")
            verified_files(folder)
    except Exception:
        raise RuntimeError(
            "Checkpoint restore failed or incompatible data found; remote data were not modified. Check destination and credentials."
        ) from None
    print(
        "Restored and verified", restored, "shards from private HF revision", info.sha
    )


def download(repo, folder):
    from huggingface_hub import hf_hub_download

    api, info = private_api(repo)
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    credential = token()
    try:
        index = hf_hub_download(
            repo,
            "esm_stage1/completed.json",
            repo_type="dataset",
            revision=info.sha,
            token=credential,
        )
        completed = json.loads(Path(index).read_text())
        for name, digest in completed["files"].items():
            assert Path(name).name == name and (
                name == "extraction.json"
                or (
                    name.startswith("shard_") and Path(name).suffix in [".json", ".npz"]
                )
            )
            target = folder / name
            if target.exists() and sha(target) == digest:
                continue
            cached = hf_hub_download(
                repo,
                "esm_stage1/" + name,
                repo_type="dataset",
                revision=info.sha,
                token=credential,
            )
            assert sha(cached) == digest
            shutil.copyfile(cached, target)
        shutil.copyfile(index, folder / "completed.json")
    except Exception:
        raise RuntimeError(
            "Private download or checksum verification failed; detailed network errors suppressed."
        ) from None
    verified_files(folder)
    print("Downloaded and validated all 3052 development proteins; revision", info.sha)


def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "action", choices=["upload", "download", "restore", "check", "init", "diagnose"]
    )
    p.add_argument("--repo", required=True)
    p.add_argument("--folder", type=Path, default=HERE / "activations")
    p.add_argument(
        "--create-private",
        action="store_true",
        help="With init only: explicitly create a missing private dataset",
    )
    args = p.parse_args()
    if args.create_private and args.action != "init":
        p.error("--create-private is only valid with init")
    if args.action == "diagnose":
        diagnose_access(args.repo)
        return
    if args.action in ("check", "init"):
        initialize_private(args.repo, create_private=args.create_private)
        return
    {"upload": upload, "download": download, "restore": restore}[args.action](
        args.repo, args.folder
    )


if __name__ == "__main__":
    main()
