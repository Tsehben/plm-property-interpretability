"""Offline rate-limit recovery checks; no network requests or real sleeps."""

from contextlib import redirect_stdout
from io import StringIO
from types import SimpleNamespace
from unittest import TestCase, main
from unittest.mock import Mock, patch

from save_full_maps import remote, retry_delay


def failure(status=429, headers=None):
    error = Exception("SENSITIVE_REQUEST_DETAILS")
    error.response = SimpleNamespace(status_code=status, headers=headers or {})
    return error


class RetryTests(TestCase):
    def test_honors_server_headers(self):
        self.assertEqual(retry_delay(failure(headers={"Retry-After": "120"}), 0), 122)
        self.assertEqual(
            retry_delay(failure(headers={"RateLimit": '"api";r=0;t=302'}), 0), 304
        )
        self.assertEqual(
            retry_delay(
                failure(headers={"Retry-After": "120", "RateLimit": '"api";r=0;t=302'}),
                0,
            ),
            304,
        )
        with patch("save_full_maps.time.time", return_value=0):
            self.assertEqual(
                retry_delay(
                    failure(headers={"Retry-After": "Thu, 01 Jan 1970 00:02:00 GMT"}), 0
                ),
                122,
            )

    def test_fallback_and_invalid_headers(self):
        self.assertEqual(retry_delay(failure(), 0), 60)
        self.assertEqual(retry_delay(failure(), 1), 120)
        self.assertEqual(retry_delay(failure(), 12), 900)
        self.assertEqual(retry_delay(failure(headers={"Retry-After": "nan"}), 0), 60)

    def test_retry_keeps_commit_arguments_and_sanitizes_output(self):
        call = Mock(side_effect=[failure(headers={"Retry-After": "61"}), "committed"])
        call.__name__ = "create_commit"
        output = StringIO()
        with patch("save_full_maps.time.sleep") as sleep, redirect_stdout(output):
            result = remote(call, repo_id="private/repo", parent_commit="same-parent")
        self.assertEqual(result, "committed")
        self.assertEqual(call.call_count, 2)
        self.assertEqual(call.call_args_list[0], call.call_args_list[1])
        self.assertEqual([c.args[0] for c in sleep.call_args_list], [60, 3])
        self.assertNotIn("SENSITIVE", output.getvalue())

    def test_no_retry_of_auth_conflict_or_unknown_commit_failure(self):
        for status in (401, 403, 409, 500, None):
            call = Mock(side_effect=failure(status))
            call.__name__ = "create_commit"
            with (
                patch("save_full_maps.time.sleep") as sleep,
                self.assertRaises(RuntimeError) as caught,
            ):
                remote(call)
            self.assertEqual(call.call_count, 1)
            sleep.assert_not_called()
            self.assertNotIn("SENSITIVE", str(caught.exception))

    def test_retry_limit_is_finite(self):
        call = Mock(side_effect=failure())
        call.__name__ = "create_commit"
        with (
            patch("save_full_maps.time.sleep"),
            redirect_stdout(StringIO()),
            self.assertRaises(RuntimeError),
        ):
            remote(call)
        self.assertEqual(call.call_count, 13)


if __name__ == "__main__":
    main()
