"""Record, fold coverage and determinism of the five-fold baseline artifact.

Runs the full ticket-#8 record on the synthetic fixture (fold 0, all four
baselines) and checks the properties #9 pairs on:

* every test observation appears exactly once and every test condition is scored
  exactly once; the test population and ``c*`` are identical across baselines;
* the record files exist with the documented schema and the manifest carries the
  input hashes + git SHA + config;
* re-running a fold produces byte-identical predictions, metrics and aggregates
  (the baseline outputs are closed-form, so this must be exact);
* the resume path reuses a completed fold without changing its files.
"""

from __future__ import annotations

import json
import unittest
from pathlib import Path

import numpy as np

import _support
import baselines


class FullRecordTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out_dir = _support.tmp_root() / "record" / "06.8_baselines" / "output"
        cls.result = baselines.run_all(
            data_dir=_support.fixture_dir(),
            out_dir=cls.out_dir,
            folds=[0],
            baselines=list(baselines.BASELINES),
            command=["test_record"],
        )
        cls.summary = cls.result["summary"]

    def test_tracked_record_files_exist(self):
        for name in (
            "summary_table.csv",
            "metrics.csv",
            "per_drug_dose.csv",
            "fold_coverage.csv",
            "manifest.json",
        ):
            self.assertTrue((self.out_dir / name).is_file(), name)

    def test_summary_table_has_per_fold_and_macro_rows(self):
        baselines_in_table = set(self.summary["baseline"].unique())
        self.assertEqual(baselines_in_table, set(baselines.BASELINES))
        per_fold = self.summary[self.summary["scope"] == "per_fold"]
        macro = self.summary[self.summary["scope"] == "macro_over_folds"]
        self.assertEqual(len(per_fold), len(baselines.BASELINES))
        self.assertEqual(len(macro), len(baselines.BASELINES))
        # a fold is a single seed: the macro row equals the one per-fold value
        for _, row in macro.iterrows():
            value = float(per_fold[per_fold["baseline"] == row["baseline"]]["auprc"].iloc[0])
            self.assertAlmostEqual(float(row["auprc"]), value, places=9)

    def test_every_test_condition_scored_exactly_once(self):
        coverage = self.result["manifest"]["folds"]["0"]
        self.assertTrue(coverage["scored_conditions_unique"])
        self.assertEqual(coverage["n_scored_conditions"], coverage["n_test_conditions"])
        self.assertEqual(coverage["n_scored_replicates"], coverage["n_test_rows"])
        self.assertEqual(coverage["n_test_conditions_seen_in_train"], 0)
        import pandas as pd

        frame = pd.read_csv(self.out_dir / "per_drug_dose.csv")
        for baseline in baselines.BASELINES:
            sub = frame[frame["baseline"] == baseline]
            self.assertEqual(sub["condition"].nunique(), coverage["n_test_conditions"])
            self.assertEqual(int(sub["n_replicates"].sum()), coverage["n_test_rows"])
        # identical test population + c* across baselines (the pairing guarantee)
        invariants = self.result["manifest"]["invariants"]
        self.assertTrue(invariants["same_test_population_across_baselines"])
        self.assertTrue(invariants["c_star_identical_across_baselines"])
        self.assertTrue(invariants["every_main_drug_tested_exactly_once"])

    def test_no_change_metrics_match_the_constant_boundary(self):
        import pandas as pd

        frame = pd.read_csv(self.out_dir / "per_drug_dose.csv")
        sub = frame[frame["baseline"] == "no_change"]
        self.assertTrue(bool((sub["pred_constant"]).all()))
        coverage = self.result["manifest"]["folds"]["0"]
        n_genes = int(self.result["manifest"]["data"]["gene_axis"]["n_genes"])
        expected_chance = float(
            np.mean((sub["n_up_200"].to_numpy() + sub["n_down_200"].to_numpy()) / (2.0 * n_genes))
        )
        summary = pd.read_csv(self.out_dir / "summary_table.csv")
        row = summary[(summary["baseline"] == "no_change") & (summary["scope"] == "per_fold")].iloc[0]
        self.assertAlmostEqual(float(row["auprc"]), expected_chance, places=6)
        self.assertGreater(int(coverage["n_test_rows"]), 0)

    def test_manifest_carries_input_git_and_config(self):
        manifest = self.result["manifest"]
        self.assertEqual(manifest["schema"], baselines.SCHEMA)
        self.assertEqual(manifest["ticket"], 8)
        self.assertEqual(len(manifest["input"]["artifact_sha256"]), 9)
        self.assertFalse(any(k.endswith("manifest.json") for k in manifest["input"]["artifact_sha256"]))
        self.assertIn("sha256", manifest["input"]["locked_config"])
        self.assertTrue(manifest["input"]["gene_axis"]["n_genes"] > 0)
        self.assertIn("git", manifest)
        self.assertEqual(manifest["config"]["ridge_alpha"], baselines.RIDGE_ALPHA)
        self.assertEqual(manifest["config"]["metric"]["headline_k"], 200)
        for key in ("summary_table", "metrics", "per_drug_dose", "fold_coverage"):
            self.assertIn(key, manifest["record_sha256"])


class DeterminismTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.frozen = _support.open_frozen()
        cls.folds = _support.folds_table(cls.frozen)
        cls.tmp = _support.tmp_root() / "determinism"
        cls.tmp.mkdir(parents=True, exist_ok=True)
        cls.dir_a = cls.tmp / "a"
        cls.dir_b = cls.tmp / "b"

    def _file_hashes(self, root: Path) -> dict[str, str]:
        from ptlib import manifest as manifest_mod

        out = {}
        for path in sorted(root.rglob("*")):
            if not path.is_file():
                continue
            relative = path.relative_to(root).as_posix()
            if relative.endswith("manifest.json") or relative.endswith("fold_record.json"):
                continue  # timestamps differ by design
            if relative.endswith(".npz"):
                # npz containers embed zip timestamps; hash the arrays instead
                with np.load(path) as payload:
                    out[relative] = baselines.arrays_sha256(
                        {key: payload[key] for key in payload.files}
                    )
                continue
            out[relative] = manifest_mod.sha256_file(path)
        return out

    def test_rerunning_a_fold_is_byte_identical(self):
        for fold in (0,):
            baselines.run_fold(
                self.frozen, self.folds, fold, self.dir_a,
                list(baselines.BASELINES), command=["test_determinism"],
            )
            baselines.run_fold(
                self.frozen, self.folds, fold, self.dir_b,
                list(baselines.BASELINES), command=["test_determinism"],
            )
            hashes_a = self._file_hashes(self.dir_a)
            hashes_b = self._file_hashes(self.dir_b)
            self.assertEqual(
                set(hashes_a), set(hashes_b),
                "a rerun produced a different set of record files",
            )
            for relative, digest in hashes_a.items():
                if relative.endswith("run.log"):
                    continue
                self.assertEqual(
                    digest, hashes_b[relative], f"{relative} is not deterministic"
                )
            # and the metrics themselves agree exactly
            for baseline in baselines.BASELINES:
                aggregate_a = json.loads(
                    (self.dir_a / "folds" / f"fold{fold}" / baseline / "eval" / "aggregate.json")
                    .read_text(encoding="utf-8")
                )
                aggregate_b = json.loads(
                    (self.dir_b / "folds" / f"fold{fold}" / baseline / "eval" / "aggregate.json")
                    .read_text(encoding="utf-8")
                )
                self.assertEqual(aggregate_a["by_k"]["200"], aggregate_b["by_k"]["200"])

    def test_resume_reuses_without_changing_files(self):
        from ptlib import manifest as manifest_mod

        dir_c = self.tmp / "resume"
        first = baselines.run_fold(
            self.frozen, self.folds, 0, dir_c, list(baselines.BASELINES),
            command=["test_resume"],
        )
        hashes_first = self._file_hashes(dir_c)
        second = baselines.run_fold(
            self.frozen, self.folds, 0, dir_c, list(baselines.BASELINES),
            command=["test_resume"],
        )
        self.assertTrue(all(record["prediction"].get("reused") for record in second["baselines"]))
        hashes_second = self._file_hashes(dir_c)
        for relative, digest in hashes_first.items():
            self.assertEqual(digest, hashes_second[relative], relative)
        for record in second["baselines"]:
            self.assertTrue(record["eval_manifest_sha256"])
        self.assertEqual(first["coverage"], second["coverage"])


class CheckInputTests(unittest.TestCase):
    def test_check_reports_every_fold_without_writing(self):
        report = baselines.check_input(data_dir=_support.fixture_dir())
        self.assertEqual(report["status"], "ok")
        self.assertEqual(report["n_genes"], 18164)
        self.assertEqual(report["n_folds"], 5)
        self.assertEqual(report["headline_k"], 200)
        self.assertEqual(sorted(report["folds"]), ["0", "1", "2", "3", "4"])
        for fold, record in report["folds"].items():
            self.assertEqual(record["n_test_conditions_seen_in_train"], 0)
            self.assertGreater(record["n_train_rows"], 0)
            self.assertGreater(record["n_test_rows"], 0)
            self.assertGreater(record["n_test_drugs"], 0)


if __name__ == "__main__":
    unittest.main()
