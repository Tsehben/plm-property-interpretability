"""Offline reopening tests; remote services and real model computation are mocked."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest import TestCase, main
from unittest.mock import MagicMock, patch

import extract
import hf_storage as storage
import pilot_full_map_precision as precision
from common import MODEL, REVISION, read_extraction, sha, write_json


class ReopenTests(TestCase):
    def test_complete_local_restore_downloads_nothing(self):
        with TemporaryDirectory() as tmp:
            folder = Path(tmp)
            write_json(folder / "completed.json", {})
            api = MagicMock()
            api.get_paths_info.return_value = [
                SimpleNamespace(
                    path="esm_stage1/completed.json",
                    lfs=None,
                    blob_id=storage.git_blob_sha(folder / "completed.json"),
                )
            ]
            with (
                patch.object(
                    storage,
                    "private_api",
                    return_value=(api, SimpleNamespace(sha="pin")),
                ),
                patch.object(storage, "token", return_value="synthetic"),
                patch.object(storage, "verified_files") as verify,
                patch("huggingface_hub.hf_hub_download") as download,
            ):
                storage.restore("owner/private", folder)
                verify.assert_called_once_with(folder)
                download.assert_not_called()
                api.list_repo_files.assert_not_called()

    def test_new_runtime_and_damaged_cache_fetch_only_needed_files(self):
        with TemporaryDirectory() as src, TemporaryDirectory() as dest:
            src, dest = Path(src), Path(dest)
            for name in ("extraction.json", "shard_0000.json", "shard_0000.npz"):
                (src / name).write_bytes(b"synthetic:" + name.encode())
            write_json(
                src / "completed.json",
                dict(
                    shape=[3052, 33, 1280],
                    files={p.name: sha(p) for p in src.iterdir()},
                ),
            )
            api = MagicMock()
            api.list_repo_files.return_value = [
                "esm_stage1/" + p.name for p in src.iterdir()
            ]
            fetched = []

            def fetch(repo, name, **kwargs):
                self.assertEqual(kwargs["revision"], "pin")
                fetched.append(Path(name).name)
                return str(src / Path(name).name)

            with (
                patch.object(
                    storage,
                    "private_api",
                    return_value=(api, SimpleNamespace(sha="pin")),
                ),
                patch.object(storage, "token", return_value="synthetic"),
                patch.object(storage, "verified_files") as verify,
                patch("huggingface_hub.hf_hub_download", side_effect=fetch),
            ):
                storage.restore("owner/private", dest)
                self.assertEqual(set(fetched), {p.name for p in src.iterdir()})
                fetched.clear()
                (dest / "shard_0000.npz").write_bytes(b"damaged")
                verify.side_effect = [AssertionError("damaged"), {}]
                storage.restore("owner/private", dest)
                self.assertEqual(fetched, ["completed.json", "shard_0000.npz"])

    def test_complete_extraction_and_pilot_never_load_model(self):
        with TemporaryDirectory() as tmp:
            write_json(Path(tmp) / "completed.json", {})
            for action in ("pilot", "extract"):
                with (
                    patch(
                        "sys.argv",
                        [
                            "extract",
                            action,
                            "--output",
                            tmp,
                            "--device",
                            "cpu",
                            "--precision",
                            "float32",
                        ],
                    ),
                    patch.object(
                        extract,
                        "verify_config",
                        return_value={"extraction_dtype": "float32"},
                    ),
                    patch.object(storage, "verified_files") as verify,
                    patch.object(extract, "load_model") as model,
                ):
                    extract.main()
                    verify.assert_called_once()
                    model.assert_not_called()
                    verify.side_effect = AssertionError("bad cache")
                    with self.assertRaises(AssertionError):
                        extract.main()
                    model.assert_not_called()

    def test_saved_extraction_pilot_reused_only_for_matching_environment(self):
        with TemporaryDirectory() as tmp:
            write_json(
                Path(tmp) / "pilot.json",
                dict(
                    passed=True,
                    model=MODEL,
                    revision=REVISION,
                    precision="float32",
                    device="cpu",
                    software={"version": "pinned"},
                ),
            )
            with (
                patch(
                    "sys.argv", ["extract", "pilot", "--output", tmp, "--device", "cpu"]
                ),
                patch.object(extract, "versions", return_value={"version": "pinned"}),
                patch.object(extract, "load_model") as model,
            ):
                extract.main()
                model.assert_not_called()

    def test_previous_sae_pilot_uses_exact_cohort_precision_and_code(self):
        sample = precision.select_sample(read_extraction())
        config = dict(extraction_dtype="float16", software={"test": "version"})
        report = dict(
            pilot_script_sha256=precision.PRIOR_PILOT_SHA,
            model=MODEL,
            model_revision=REVISION,
            sae_repo=precision.SAE_REPO,
            sae_revision=precision.INTERPROT_REVISION,
            sae_source_sha256=precision.SOURCE_SHA,
            sae_topk=128,
            inference_precision="float16",
            software=config["software"],
            manifest_sha256=sha(precision.HERE / "data/extraction_manifest.tsv"),
            pooled_completion_sha256="completion",
            engineering_tolerances=precision.LIMITS,
            proteins=len(sample),
            results=[
                dict(
                    record_id=r.record_id,
                    sequence_sha256=r.sha256,
                    split=r.split,
                    length=int(r.length),
                    layer=layer,
                )
                for r in sample.itertuples()
                for layer in precision.SAE_LAYERS
            ],
        )
        self.assertTrue(precision.report_matches(report, sample, config, "completion"))
        self.assertFalse(precision.report_matches(report, sample, config, "other"))
        for key, value in (
            ("inference_precision", "float32"),
            ("pilot_script_sha256", "unknown-code"),
            ("sae_source_sha256", "different"),
            ("results", report["results"][:-1]),
        ):
            changed = json.loads(json.dumps(report))
            changed[key] = value
            self.assertFalse(
                precision.report_matches(changed, sample, config, "completion")
            )


if __name__ == "__main__":
    main()
