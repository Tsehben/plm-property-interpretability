"""No-network tests for authenticated private repository setup and diagnostics."""

from types import SimpleNamespace
from unittest.mock import MagicMock, patch
import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
import hf_storage as hf


def http_error(status):
    error = RuntimeError("synthetic-sensitive-response-not-to-be-printed")
    error.response = SimpleNamespace(status_code=status)
    return error


class HFAccessTests(unittest.TestCase):
    def test_creation_denial_reports_status_without_response_body(self):
        api = MagicMock()
        api.repo_info.side_effect = http_error(404)
        api.create_repo.side_effect = http_error(403)
        with patch.object(hf, "authenticated_api", return_value=api):
            with self.assertRaisesRegex(
                RuntimeError, "creation failed.*HTTP 403"
            ) as error:
                hf.initialize_private("PLM-Interp/test", create_private=True)
        self.assertNotIn("synthetic-sensitive", str(error.exception))
        self.assertIn("creation was denied", str(error.exception))

    def test_successful_creation_followed_by_failed_verification_is_distinct(self):
        api = MagicMock()
        api.repo_info.side_effect = [http_error(404), http_error(403)]
        with patch.object(hf, "authenticated_api", return_value=api):
            with self.assertRaisesRegex(
                RuntimeError, "creation succeeded.*verification failed.*HTTP 403"
            ):
                hf.initialize_private("PLM-Interp/test", create_private=True)
        api.create_repo.assert_called_once()

    def test_diagnosis_is_read_only_and_never_prints_token_details(self):
        api = MagicMock()
        api.whoami.return_value = {
            "name": "researcher",
            "auth": {
                "accessToken": {
                    "role": "fineGrained",
                    "token": "synthetic-secret",
                    "displayName": "private-token-name",
                }
            },
            "orgs": [{"name": "PLM-Interp", "roleInOrg": "contributor"}],
        }
        api.repo_info.side_effect = http_error(404)
        output = StringIO()
        with (
            patch.object(hf, "authenticated_api", return_value=api),
            redirect_stdout(output),
        ):
            hf.diagnose_access("PLM-Interp/test")
        text = output.getvalue()
        self.assertIn("fineGrained", text)
        self.assertIn("contributor", text)
        self.assertIn("404", text)
        for sensitive in [
            "synthetic-secret",
            "private-token-name",
            "synthetic-sensitive",
        ]:
            self.assertNotIn(sensitive, text)
        api.create_repo.assert_not_called()
        api.create_commit.assert_not_called()

    def test_bad_token_is_distinguished_without_disclosing_response(self):
        with (
            patch.object(hf, "token", return_value="synthetic"),
            patch("huggingface_hub.HfApi") as factory,
        ):
            api = factory.return_value
            api.whoami.side_effect = http_error(401)
            with self.assertRaisesRegex(RuntimeError, "authentication failed") as error:
                hf.private_api("PLM-Interp/test")
            self.assertNotIn("synthetic-sensitive", str(error.exception))
            api.repo_info.assert_not_called()
            api.create_repo.assert_not_called()

    def test_missing_repo_does_not_trigger_automatic_creation(self):
        api = MagicMock()
        api.repo_info.side_effect = http_error(404)
        with patch.object(hf, "authenticated_api", return_value=api):
            with self.assertRaisesRegex(RuntimeError, "missing or hidden"):
                hf.initialize_private("PLM-Interp/test")
        api.create_repo.assert_not_called()

    def test_explicit_creation_is_private_dataset_only(self):
        api = MagicMock()
        api.repo_info.side_effect = [
            http_error(404),
            SimpleNamespace(private=True, sha=None),
        ]
        with patch.object(hf, "authenticated_api", return_value=api):
            hf.initialize_private("PLM-Interp/test", create_private=True)
        api.create_repo.assert_called_once_with(
            repo_id="PLM-Interp/test", repo_type="dataset", private=True, exist_ok=False
        )
        api.create_commit.assert_not_called()

    def test_public_repo_is_never_reconfigured_or_used(self):
        api = MagicMock()
        api.repo_info.return_value = SimpleNamespace(private=False)
        with patch.object(hf, "authenticated_api", return_value=api):
            with self.assertRaisesRegex(ValueError, "non-private"):
                hf.initialize_private("PLM-Interp/test", create_private=True)
        api.create_repo.assert_not_called()
        api.update_repo_settings.assert_not_called()

    def test_forbidden_repo_does_not_trigger_creation(self):
        api = MagicMock()
        api.repo_info.side_effect = http_error(403)
        with patch.object(hf, "authenticated_api", return_value=api):
            with self.assertRaisesRegex(RuntimeError, "forbidden"):
                hf.initialize_private("PLM-Interp/test", create_private=True)
        api.create_repo.assert_not_called()

    def test_empty_private_repo_can_start_without_a_revision(self):
        api = MagicMock()
        with (
            tempfile.TemporaryDirectory() as folder,
            patch.object(
                hf,
                "private_api",
                return_value=(api, SimpleNamespace(private=True, sha=None)),
            ),
            patch.object(hf, "token", return_value="synthetic"),
        ):
            hf.restore("PLM-Interp/test", folder)
        api.list_repo_files.assert_not_called()
        api.create_repo.assert_not_called()


if __name__ == "__main__":
    unittest.main(verbosity=2)
