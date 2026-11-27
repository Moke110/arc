"""Check that the formal publication refuses a tampered scored table."""

from hashlib import sha256
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dose_formal import checked_report, REPORT_TABLES


class FormalReportTests(unittest.TestCase):
    def test_report_hash_is_checked_before_publication(self):
        with TemporaryDirectory() as work:
            report = Path(work)
            for file in REPORT_TABLES:
                (report / file).write_text("original\n", encoding="utf-8")
            record = {
                "schema": "carps-vae-dose-report-v1", "bounded": False,
                "primary_metrics": ["Δ PCC", "Δ RMSE"], "scores": 900,
                "outputs": {name: sha256(b"original\n").hexdigest() for name in REPORT_TABLES},
            }
            (report / "manifest.json").write_text(json.dumps(record), encoding="utf-8")
            self.assertEqual(checked_report(report), record)
            (report / "foldwise-summary.csv").write_text("changed\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "formal report"):
                checked_report(report)


if __name__ == "__main__":
    unittest.main()
