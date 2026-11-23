"""Parser checks against the saved authors' server pilot; no network or test labels."""

import unittest

from run_aggrescan_server import WORK, parse_result, read_fasta, verify_spreadsheet


class AggrescanServerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pilot = WORK / "pilot"
        cls.expected = read_fasta(pilot / "input.fasta")
        cls.html = (pilot / "crlf_submission.html").read_text()
        cls.spreadsheet = (pilot / "results.txt.xls").read_bytes()

    def test_pilot_identity_metrics_and_coordinates(self):
        proteins, residues, hotspots = parse_result(self.html, self.expected)
        verify_spreadsheet(self.spreadsheet, proteins)
        by_id = {r["record_id"]: r for r in proteins}
        self.assertEqual(by_id["NESG_raw:100"]["Na4vSS"], -3.7)
        self.assertEqual(by_id["NESG_raw:100"]["nHS"], 7)
        self.assertEqual(by_id["NESG_raw:1005"]["nHS"], 9)
        self.assertEqual(len(residues), 475)
        self.assertEqual(len(hotspots), 16)

    def test_rejects_correct_name_with_wrong_sequence(self):
        expected = self.expected.copy()
        uid = next(iter(expected))
        expected[uid] = "A" + expected[uid][1:]
        with self.assertRaisesRegex(AssertionError, "Sequence mismatch"):
            parse_result(self.html, expected)

    def test_rejects_missing_protein(self):
        expected = self.expected | {"NESG_raw:999999": "AAAAA"}
        with self.assertRaisesRegex(AssertionError, "Missing result IDs"):
            parse_result(self.html, expected)

    def test_rejects_spreadsheet_disagreement(self):
        proteins, _, _ = parse_result(self.html, self.expected)
        wrong = self.spreadsheet.replace(b"-3.7\t", b"-9.9\t")
        with self.assertRaises(AssertionError):
            verify_spreadsheet(wrong, proteins)

    def test_rejects_server_error_page(self):
        with self.assertRaisesRegex(AssertionError, "Missing result IDs"):
            parse_result("<html>UNIDENTIFIED ERROR</html>", self.expected)


if __name__ == "__main__":
    unittest.main()
