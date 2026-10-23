"""Offline harness tests for the Stage 2 benchmark (ticket 12).

No GPU.  The heavy case bundle is built once per test class under
`/tmp/opencode`; the full no-change dry run is the smoke path from ticket 05.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
CODE = HERE.parents[1]
BENCH = CODE.parent
DATA = BENCH / "data"
sys.path.insert(0, str(CODE))

import baselines  # noqa: E402
import evaluate  # noqa: E402
import runner  # noqa: E402
from cases import BenchData, FORBIDDEN_TEST_KEYS, META_KEYS  # noqa: E402


class DryRunTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tmp = tempfile.TemporaryDirectory(prefix="bench-dry-", dir="/tmp/opencode")
        cls.root = Path(cls.tmp.name)
        cls.output = cls.root / "output"
        cls.results = cls.root / "results"
        cls.data = BenchData(DATA)
        cls.n_genes = int(cls.data.manifest["n_genes"])
        cls.cases, cls.case_info = runner.ensure_cases(
            cls.data, 0, cls.output, force=True
        )

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmp.cleanup()

    def _run(self) -> dict:
        args = argparse.Namespace(
            dry_run=True,
            model=["no_change"],
            fold=[0],
            output_dir=self.output,
            results_dir=self.results,
            timeout=120.0,
            seed=20261021,
            force=False,
        )
        return runner.run(args)

    def test_case_bundles_are_leak_free(self) -> None:
        with np.load(self.cases["test"], allow_pickle=False) as test:
            keys = set(test.files)
            self.assertFalse(keys & FORBIDDEN_TEST_KEYS)
            self.assertTrue(set(META_KEYS) <= keys)
            self.assertEqual(len(test["well_id"]), 837)
            test_wells = set(map(str, test["well_id"]))
        with np.load(self.cases["train"], allow_pickle=False) as train:
            self.assertEqual(train["X"].shape, (2523, self.n_genes))
            self.assertEqual(train["delta"].shape, (2523, self.n_genes))
            self.assertFalse(test_wells & set(map(str, train["well_id"])))
        with np.load(self.cases["valid"], allow_pickle=False) as valid:
            self.assertEqual(valid["X"].shape, (846, self.n_genes))
            self.assertFalse(test_wells & set(map(str, valid["well_id"])))
        self.assertEqual(self.case_info["n_test_groups"], 276)

    def test_full_no_change_dry_run(self) -> None:
        result = self._run()
        statuses = result["statuses"]
        self.assertEqual(len(statuses), 1)
        self.assertEqual(statuses[0]["status"], "completed")
        leaderboard = pd.read_csv(self.results / "leaderboard.csv")
        self.assertEqual(len(leaderboard), 1)
        row = leaderboard.iloc[0]
        self.assertEqual(row["model"], "no_change")
        self.assertAlmostEqual(float(row["headline_auprc"]), 200 / self.n_genes, places=12)
        self.assertTrue(math.isnan(float(row["pcc"])))
        self.assertTrue(np.isfinite(float(row["r2"])))
        self.assertFalse(bool(row["winner_eligible"]))
        self.assertEqual(int(row["n_folds"]), 1)
        pred = np.load(self.output / "no_change" / "fold_0" / "pred_test.npy")
        self.assertEqual(pred.shape, (837, self.n_genes))
        self.assertTrue(np.array_equal(pred, np.zeros_like(pred)))

    def test_result_artifacts_and_schema(self) -> None:
        if not (self.results / "leaderboard.csv").exists():
            self._run()
        metrics = pd.read_csv(self.results / "metrics.csv")
        self.assertTrue(len(metrics) > 0)
        self.assertEqual(
            list(metrics.columns),
            [
                "model", "fold", "level", "subset", "group_id", "drug_key",
                "dose_um", "k", "direction", "metric", "value", "nan_reason",
            ],
        )
        self.assertIn("headline", set(metrics["level"]))
        self.assertIn("atc-moa", set(metrics["subset"]))
        per_drug = pd.read_csv(self.results / "per_drug.csv")
        self.assertEqual(len(per_drug), 92)
        with np.load(self.results / "bootstrap.npz", allow_pickle=False) as boot:
            self.assertEqual(boot["drug_idx"].shape, (2000, 5, 92))
            self.assertEqual(boot["fold_auprc"].shape, (2000, 5, 1))
            self.assertEqual(boot["headline"].shape, (2000, 1))
            self.assertEqual(int(boot["seed"]), 20261022)
        self.assertTrue((self.output / "run_status.csv").exists())


class PerfectEvaluationTest(unittest.TestCase):
    def test_perfect_prediction_leaderboard(self) -> None:
        data = BenchData(DATA)
        with tempfile.TemporaryDirectory(prefix="bench-perfect-", dir="/tmp/opencode") as tmp:
            root = Path(tmp)
            output = root / "output"
            results = root / "results"
            run_dir = output / "perfect" / "fold_0"
            run_dir.mkdir(parents=True)
            # Use the shared case bundle when present, otherwise build it.
            cases, _ = runner.ensure_cases(data, 0, output, force=True)
            with np.load(cases["test"], allow_pickle=False) as test:
                group_ids = np.asarray(test["group_id"]).astype(str)
            true = np.stack([data.group_delta[data.group_row[g]] for g in group_ids]).astype(np.float32)
            np.save(run_dir / "pred_test.npy", true)
            (run_dir / "group_ids.json").write_text(
                json.dumps({"valid": [], "test": group_ids.tolist()}) + "\n"
            )
            (run_dir / "status.json").write_text(json.dumps({"status": "completed"}) + "\n")
            result = evaluate.evaluate_run(
                output,
                results,
                models=["perfect"],
                folds=[0],
                registry={"perfect": {"ideal_fit": True}},
                data_dir=DATA,
            )
            row = result["leaderboard"].iloc[0]
            self.assertAlmostEqual(float(row["headline_auprc"]), 1.0, places=12)
            self.assertAlmostEqual(float(row["pcc"]), 1.0, places=12)
            self.assertAlmostEqual(float(row["spearman"]), 1.0, places=12)
            self.assertAlmostEqual(float(row["r2"]), 1.0, places=12)
            self.assertAlmostEqual(float(row["rmse"]), 0.0, places=12)
            self.assertTrue(bool(row["ideal_fit"]))
            self.assertEqual(int(row["nan_counts"]), 0)


class IdealFitLabelTest(unittest.TestCase):
    def test_string_labels_are_preserved(self) -> None:
        self.assertEqual(evaluate._ideal_fit("partial"), "partial")
        self.assertEqual(evaluate._ideal_fit("full"), "full")
        self.assertEqual(evaluate._ideal_fit(True), "full")
        self.assertEqual(evaluate._ideal_fit(False), "none")


class BaselineCoreTest(unittest.TestCase):
    def test_all_baselines_fit_predict(self) -> None:
        rng = np.random.default_rng(7)
        n_genes = 5000
        groups = ["D1|0.2", "D1|1", "D1|10", "D2|0.2", "D2|1", "D2|10"]
        drugs = ["D1", "D1", "D1", "D2", "D2", "D2"]
        doses = [0.2, 1.0, 10.0, 0.2, 1.0, 10.0]
        smiles = ["CCO", "CCO", "CCO", "c1ccccc1", "c1ccccc1", "c1ccccc1"]
        well_group = np.repeat(groups, 2)
        well_drug = np.repeat(drugs, 2)
        well_dose = np.repeat(doses, 2)
        well_smiles = np.repeat(smiles, 2)
        well_delta = rng.normal(size=(len(well_group), n_genes)).astype(np.float32)
        with tempfile.TemporaryDirectory(prefix="bench-baselines-", dir="/tmp/opencode") as tmp:
            root = Path(tmp)
            train_path = root / "train.npz"
            test_path = root / "test.npz"
            np.savez(
                train_path,
                group_id=well_group.astype(str),
                drug_key=well_drug.astype(str),
                dose_um=np.asarray(well_dose, dtype=np.float64),
                smiles=well_smiles.astype(str),
                delta=well_delta,
            )
            np.savez(
                test_path,
                well_id=np.asarray([f"w{i}" for i in range(len(groups))], dtype=str),
                group_id=np.asarray(groups, dtype=str),
                drug_key=np.asarray(drugs, dtype=str),
                dose_um=np.asarray(doses, dtype=np.float64),
                smiles=np.asarray(smiles, dtype=str),
            )
            group_means = np.stack(
                [well_delta[well_group == g].mean(axis=0) for g in groups]
            )
            expected_mean = group_means.mean(axis=0)
            expected_dose = {
                d: group_means[np.asarray(doses) == d].mean(axis=0)
                for d in (0.2, 1.0, 10.0)
            }
            config = {
                "baseline_params": {"fp_bits": 2048, "fp_radius": 2, "ridge_alpha": 1.0}
            }
            for name in ("delta_mean", "delta_mean_dose", "ridge_fp", "nearest_drug"):
                run_dir = root / name
                checkpoint = baselines.fit_baseline(name, train_path, train_path, config, run_dir)
                out = root / f"{name}.npy"
                pred = baselines.predict_baseline(name, test_path, checkpoint, out)
                self.assertEqual(pred.shape, (len(groups), n_genes), name)
                self.assertTrue(np.isfinite(pred).all(), name)
                if name == "delta_mean":
                    self.assertTrue(np.allclose(pred, expected_mean, atol=1e-6), name)
                if name == "delta_mean_dose":
                    for i, dose in enumerate(doses):
                        self.assertTrue(np.allclose(pred[i], expected_dose[dose], atol=1e-6), name)


if __name__ == "__main__":
    unittest.main()
