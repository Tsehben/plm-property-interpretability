"""Regression tests for the CPU/duplicate-metadata/source-build installation failure."""

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import tempfile
import unittest
import importlib.util
import subprocess
import sys
import venv
import zipfile
import install_colab_dependencies as setup


class ColabSetupTests(unittest.TestCase):
    def test_real_pip_writes_only_to_disposable_target_environment(self):
        self.assertIsNone(importlib.util.find_spec("nesg_setup_check"))
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            environment = root / "venv"
            venv.EnvBuilder(system_site_packages=True, with_pip=False).create(
                environment
            )
            python = environment / "bin/python"
            wheel = root / "nesg_setup_check-0.0.0-py3-none-any.whl"
            prefix = "nesg_setup_check-0.0.0.dist-info/"
            files = {
                "nesg_setup_check.py": "VALUE = 7\n",
                prefix
                + "METADATA": "Metadata-Version: 2.1\nName: nesg-setup-check\nVersion: 0.0.0\n",
                prefix
                + "WHEEL": "Wheel-Version: 1.0\nGenerator: test\nRoot-Is-Purelib: true\nTag: py3-none-any\n",
            }
            files[prefix + "RECORD"] = "".join(
                f"{name},,\n" for name in [*files, prefix + "RECORD"]
            )
            with zipfile.ZipFile(wheel, "w") as archive:
                for name, value in files.items():
                    archive.writestr(name, value)
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "pip",
                    "--python",
                    str(python),
                    "install",
                    "--no-index",
                    "--no-deps",
                    "--only-binary=:all:",
                    "--disable-pip-version-check",
                    str(wheel),
                ],
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )
            subprocess.run(
                [
                    str(python),
                    "-c",
                    "import nesg_setup_check; assert nesg_setup_check.VALUE == 7",
                ],
                check=True,
                capture_output=True,
                text=True,
                timeout=30,
            )
            self.assertIsNone(importlib.util.find_spec("nesg_setup_check"))

    def test_cpu_preflight_stops_before_environment_mutation(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "environment"
            with (
                patch("sys.argv", ["setup", "--env-dir", str(target)]),
                patch.object(
                    setup, "gpu_preflight", side_effect=RuntimeError("CPU only")
                ),
                patch.object(setup.venv, "EnvBuilder") as builder,
                patch.object(setup.subprocess, "run") as run,
            ):
                with self.assertRaisesRegex(RuntimeError, "CPU only"):
                    setup.main()
                builder.assert_not_called()
                run.assert_not_called()
                self.assertFalse(target.exists())

    def test_inactive_numpy_metadata_does_not_override_verified_version(self):
        dist = [
            SimpleNamespace(metadata={"Name": "numpy"}, version=v)
            for v in ["1.26.4", "2.4.6"]
        ]
        for entries in [dist, dist[::-1]]:
            with (
                patch.object(setup.metadata, "distributions", return_value=entries),
                patch.object(setup.metadata, "version", return_value="2.4.6"),
            ):
                self.assertEqual(
                    setup.active_versions({"numpy": "2.4.6"}), {"numpy": "2.4.6"}
                )

    def test_actual_import_metadata_disagreement_still_fails(self):
        dist = [SimpleNamespace(metadata={"Name": "numpy"}, version="1.26.4")]
        with (
            patch.object(setup.metadata, "distributions", return_value=dist),
            patch.object(setup.metadata, "version", return_value="1.26.4"),
        ):
            with self.assertRaisesRegex(
                RuntimeError, "Active metadata/import mismatch for numpy"
            ):
                setup.active_versions({"numpy": "2.4.6"})

    def test_pip_targets_separate_interpreter_and_forbids_source_builds(self):
        command = setup.pip_command(
            Path("/content/nesg_esm_env/bin/python"), Path("/tmp/constraints")
        )
        self.assertEqual(
            command[command.index("--python") + 1], "/content/nesg_esm_env/bin/python"
        )
        self.assertIn("--only-binary=:all:", command)
        self.assertNotIn("--ignore-installed", command)
        self.assertEqual(setup.child_environment()["USE_TF"], "0")
        self.assertEqual(setup.child_environment()["USE_FLAX"], "0")

    def test_restart_reuses_managed_env_and_keeps_correct_command_order(self):
        info = {"torch": "synthetic-gpu", "cuda": "synthetic", "gpu": "synthetic"}
        with tempfile.TemporaryDirectory() as tmp:
            target = (Path(tmp) / "environment").resolve()

            def fake_create(directory):
                (directory / "bin").mkdir(parents=True)
                (directory / "bin/python").touch()

            with (
                patch("sys.argv", ["setup", "--env-dir", str(target)]),
                patch.object(setup, "gpu_preflight", return_value=info),
                patch.object(setup, "active_versions", return_value={"torch": "2.8.0"}),
                patch.object(setup.venv, "EnvBuilder") as builder,
                patch.object(setup.subprocess, "run") as run,
            ):
                builder.return_value.create.side_effect = fake_create
                # First run: missing pins -> install -> imports. Second: pins
                # match -> imports only (no pip or environment recreation).
                run.side_effect = [
                    SimpleNamespace(returncode=n) for n in [1, 0, 0, 0, 0]
                ]
                setup.main()
                self.assertEqual(
                    run.call_args_list[1].args[0][4], str(target / "bin/python")
                )
                self.assertEqual(
                    run.call_args_list[2].args[0][0], str(target / "bin/python")
                )
                setup.main()
                builder.return_value.create.assert_called_once()
                self.assertEqual(run.call_count, 5)
                self.assertEqual(sum("pip" in c.args[0] for c in run.call_args_list), 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
