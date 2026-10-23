"""Offline metric sanity cases for ticket 07.

Runs the five fixed cases against the first group of the frozen group table
(Oxazepam 0.2 µM) at k = 100, 200, 500.  No GPU, no build.

Run:
    python -m unittest discover -s stage2/03_bench_models/temp_env_all_5k/scripts/tests -v
"""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
CODE = HERE.parents[1]
BENCH = CODE.parent
DATA = BENCH / "data"
sys.path.insert(0, str(CODE))

import build_data  # noqa: E402
import metrics  # noqa: E402


class MetricSanityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.group_delta = np.load(DATA / "group_delta.npy", mmap_mode="r")
        cls.groups = pd.read_parquet(DATA / "groups.parquet")
        cls.ids = cls.groups["group_id"].astype(str).to_numpy()
        cls.i = 0
        cls.true = np.asarray(cls.group_delta[cls.i], dtype=np.float64)
        cls.n = cls.true.size
        cls.true_sets = {
            k: np.load(DATA / f"true_sets_k{k}.npz", allow_pickle=False)
            for k in (100, 200, 500)
        }

    def test_first_group_matches_contract(self) -> None:
        self.assertEqual(self.ids[self.i], "ADIMAYPTOBDMTL-UHFFFAOYSA-N|0.2")
        self.assertTrue(np.isfinite(self.true).all())
        self.assertGreater(float(self.true.std()), 0.0)

    def test_perfect_prediction(self) -> None:
        for k in (100, 200, 500):
            m = metrics.group_metrics(self.true, self.true, k)
            self.assertAlmostEqual(m["auprc"], 1.0, places=12)
            self.assertAlmostEqual(m["pcc"], 1.0, places=12)
            self.assertAlmostEqual(m["spearman"], 1.0, places=12)
            self.assertAlmostEqual(m["r2"], 1.0, places=12)
            self.assertAlmostEqual(m["rmse"], 0.0, places=12)

    def test_constant_prediction(self) -> None:
        for k in (100, 200, 500):
            m = metrics.group_metrics(self.true, np.zeros_like(self.true), k)
            self.assertAlmostEqual(m["auprc"], k / self.n, places=12)
            self.assertTrue(math.isnan(m["pcc"]))
            self.assertTrue(math.isnan(m["spearman"]))
            expected_r2 = 1.0 - float(np.mean(self.true**2)) / float(np.var(self.true))
            expected_rmse = float(np.sqrt(np.mean(self.true**2)))
            self.assertAlmostEqual(m["r2"], expected_r2, places=12)
            self.assertAlmostEqual(m["rmse"], expected_rmse, places=12)

    def test_sign_flipped_prediction(self) -> None:
        for k in (100, 200, 500):
            m = metrics.group_metrics(self.true, -self.true, k)
            self.assertAlmostEqual(m["auprc_up"], m["auprc_down"], places=12)
            golden = build_data.SANITY_GOLDEN[k]["reverse_auprc"]
            self.assertAlmostEqual(m["auprc"], golden, places=9)
            self.assertAlmostEqual(m["pcc"], -1.0, places=12)
            self.assertAlmostEqual(m["spearman"], -1.0, places=12)

    def test_random_prediction(self) -> None:
        rng = np.random.default_rng(20261021)
        pred = rng.normal(size=self.n)
        for k in (100, 200, 500):
            m = metrics.group_metrics(self.true, pred, k)
            golden = build_data.SANITY_GOLDEN[k]["random_auprc"]
            self.assertAlmostEqual(m["auprc"], golden, places=9)

    def test_top_bottom_only_prediction(self) -> None:
        for k in (100, 200, 500):
            labels = self.true_sets[k]["labels"][self.i]
            pred = np.where(labels != 0, self.true, 0.0)
            m = metrics.group_metrics(self.true, pred, k)
            self.assertAlmostEqual(m["auprc"], 1.0, places=12)
            for metric in ("pcc", "spearman", "r2"):
                self.assertAlmostEqual(
                    m[metric], build_data.SANITY_GOLDEN[k][f"masked_{metric}"], places=9
                )

    def test_persisted_labels_match_metric_code(self) -> None:
        for k in (100, 200, 500):
            labels = self.true_sets[k]["labels"][self.i]
            up = np.flatnonzero(labels == 1)
            down = np.flatnonzero(labels == -1)
            ref_up, ref_down = metrics.true_set_indices(self.true, k)
            self.assertTrue(np.array_equal(np.sort(up), np.sort(ref_up)))
            self.assertTrue(np.array_equal(np.sort(down), np.sort(ref_down)))

    def test_float32_input_tolerance(self) -> None:
        true32 = self.true.astype(np.float32)
        pred32 = (self.true + 1e-3).astype(np.float32)
        m32 = metrics.group_metrics(true32, pred32, 200)
        m64 = metrics.group_metrics(
            true32.astype(np.float64), pred32.astype(np.float64), 200
        )
        self.assertLessEqual(abs(m32["pcc"] - m64["pcc"]), 1e-12)
        self.assertLessEqual(abs(m32["auprc"] - m64["auprc"]), 1e-12)
        # A float32 round-trip is within the contract's 1e-6 PCC tolerance.
        m_ref = metrics.group_metrics(self.true, self.true + 1e-3, 200)
        self.assertLessEqual(abs(m32["pcc"] - m_ref["pcc"]), 1e-6)


if __name__ == "__main__":
    unittest.main()
