"""Check formal publication conflicts and artifact checks."""

from hashlib import sha256
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dose_formal import REPORT_TABLES, build, checked_report


class FormalReportTests(unittest.TestCase):
    def test_publish_cannot_overlap_the_local_report(self):
        with TemporaryDirectory() as work:
            bundle = Path(work) / "bundle"
            bundle.mkdir()
            report = bundle / "report"
            alias = bundle / "report-alias"
            alias.symlink_to(report, target_is_directory=True)
            for publication in (report, report / "subdir", alias):
                with self.subTest(publication=publication.name):
                    with self.assertRaisesRegex(ValueError, "publication.*overlaps.*report"):
                        build(bundle, publication, bundle / "unused-scratch", bundle)
                    self.assertFalse((bundle / ".report.publication.lock").exists())

    def test_report_hash_is_checked_before_publication(self):
        with TemporaryDirectory() as work:
            report = Path(work)
            for file in REPORT_TABLES:
                (report / file).write_text("original\n", encoding="utf-8")
            record = {
                "schema": "cpst-dose-report-v1", "bounded": False,
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
