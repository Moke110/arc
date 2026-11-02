"""Record + report tests: the tracked record derives from frozen predictions only."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import _support  # noqa: E402
import report  # noqa: E402
import write_record  # noqa: E402


class RecordTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.root = _support.clean_dir(_support.tmp_root() / "record")
        cls.pop = _support.write_synthetic_population(
            cls.root, drugs_per_fold=4, advantage=0.006
        )
        cls.baselines = _support.write_synthetic_baselines(
            cls.root, drugs_by_fold=cls.pop["drugs_by_fold"]
        )
        cls.target = cls.root / "record_out"
        cls.result = write_record.build_record(
            phase="test",
            units_path=cls.pop["units_path"],
            out=cls.pop["out"],
            data_dir=cls.pop["data"]["dir"],
            baselines_dir=cls.baselines,
            allow_subset=False,
            target=cls.target,
            make_figures=False,
            tahoe_seen_keys=set(next(iter(cls.pop["drugs_by_fold"].values()))),
        )

    def test_record_files_and_summary(self):
        for name in (
            "unit_results.csv", "fold_seed_table.csv", "per_drug_dose.csv",
            "paired_by_drug.csv", "drug_gains.csv", "strata.csv", "summary.json",
            "manifest.json", "record_check.json",
        ):
            self.assertTrue((self.target / name).exists(), f"missing {name}")
        summary = self.result["summary"]
        self.assertFalse(summary["dry_run"])
        self.assertEqual(summary["units_manifest"]["n_scored"], 30)
        self.assertEqual(summary["test_population"]["n_conditions"], 60)
        self.assertAlmostEqual(
            summary["models"]["tahoe"]["macro_auprc_200_over_folds"]
            - summary["models"]["random"]["macro_auprc_200_over_folds"],
            0.006,
            places=10,
        )

    def test_primary_bar_and_questions(self):
        summary = self.result["summary"]
        self.assertAlmostEqual(summary["primary_bar"]["auprc_200"], 0.15)
        q_a = summary["questions"]["q_a_tahoe_vs_random"]
        self.assertEqual(q_a["verdict"], "positive")
        self.assertGreater(q_a["overall"]["ci_low"], 0)
        self.assertIn("q_b_tahoe_vs_primary_baseline", summary["questions"])
        self.assertEqual(
            summary["questions"]["q_b_tahoe_vs_primary_baseline"]["verdict"], "negative"
        )
        self.assertEqual(summary["strata"]["response_label"], "unknown")
        self.assertEqual(summary["strata"]["toxicity_label"], "unknown")

    def test_gates_pass_and_record_hashes_present(self):
        check = self.result["check"]
        failed = [key for key, value in check["gates"].items() if isinstance(value, bool) and not value]
        self.assertEqual(failed, [])
        self.assertTrue(check["passed"])
        manifest = self.result["manifest"]
        self.assertIn("summary.json", manifest["record_sha256"])
        self.assertEqual(len(manifest["record_sha256"]), len(
            [name for name in manifest["record_sha256"]]
        ))

    def test_record_rederives_byte_identically(self):
        second = write_record.build_record(
            phase="test",
            units_path=self.pop["units_path"],
            out=self.pop["out"],
            data_dir=self.pop["data"]["dir"],
            baselines_dir=self.baselines,
            allow_subset=False,
            target=self.root / "record_again",
            make_figures=False,
            tahoe_seen_keys=set(next(iter(self.pop["drugs_by_fold"].values()))),
        )
        import common

        for name in ("unit_results.csv", "paired_by_drug.csv", "drug_gains.csv", "strata.csv"):
            self.assertEqual(
                common.sha256_file(self.target / name),
                common.sha256_file(Path(second["target"]) / name),
                f"{name} changed on re-derivation",
            )
        first = {
            key: value
            for key, value in self.result["summary"].items()
            if key not in ("created_utc", "figures", "tables")
        }
        again = {
            key: value
            for key, value in second["summary"].items()
            if key not in ("created_utc", "figures", "tables")
        }
        self.assertEqual(first, again)

    def test_dry_run_cannot_be_written_into_the_tracked_record_dir(self):
        with self.assertRaises(write_record.RecordError):
            write_record.build_record(
                phase="valid",
                units_path=self.pop["units_path"],
                out=self.pop["out"],
                data_dir=self.pop["data"]["dir"],
                baselines_dir=self.baselines,
                allow_subset=True,
                target=self.pop["out"] / "record",
                make_figures=False,
            )

    def test_report_section_renders_and_splices(self):
        section = report.render(self.target)
        for needle in (
            "Headline test results",
            "Acceptance question (a)",
            "delta_mean",
            "unknown",
            "Negative-result clause",
            "Direction AUPRC@200" if False else "direction AUPRC@200",
        ):
            self.assertIn(needle, section)
        readme = self.root / "README.md"
        readme.write_text(
            f"# test\n\n{report.BEGIN_MARKER}\nold\n{report.END_MARKER}\ntail\n",
            encoding="utf-8",
        )
        changed = report.splice_readme(readme, section)
        self.assertTrue(changed)
        text = readme.read_text(encoding="utf-8")
        self.assertIn("Headline test results", text)
        self.assertNotIn("\nold\n", text)
        self.assertFalse(report.splice_readme(readme, section))

    def test_missing_markers_raise(self):
        readme = self.root / "README_nomarkers.md"
        readme.write_text("# no markers\n", encoding="utf-8")
        with self.assertRaises(report.EvaluationError):
            report.splice_readme(readme, "section")


class DryRunRecordTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.root = _support.clean_dir(_support.tmp_root() / "record_dry")
        cls.pop = _support.write_synthetic_population(
            cls.root, drugs_per_fold=2, inits=("random",), folds=(0,), phase="valid"
        )
        cls.target = cls.root / "dry"
        cls.result = write_record.build_record(
            phase="valid",
            units_path=cls.pop["units_path"],
            out=cls.pop["out"],
            data_dir=cls.pop["data"]["dir"],
            baselines_dir=None,
            allow_subset=True,
            target=cls.target,
            make_figures=False,
        )

    def test_dry_run_is_labelled_and_has_no_baseline_claim(self):
        summary = self.result["summary"]
        self.assertTrue(summary["dry_run"])
        self.assertIsNone(summary["primary_bar"])
        self.assertIn("HARNESS DRY RUN", report.render(self.target))


if __name__ == "__main__":
    unittest.main()
