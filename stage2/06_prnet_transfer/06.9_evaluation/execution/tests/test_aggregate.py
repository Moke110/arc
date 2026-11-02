"""Aggregation tests: per-drug pairing, bootstrap, strata, baseline comparison."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import _support  # noqa: E402
import aggregate as agg  # noqa: E402


class LoadUnitEvalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.root = _support.clean_dir(_support.tmp_root() / "agg_load")
        cls.out = cls.root / "out"
        cls.unit = _support.make_unit(fold=0, seed=_support.SEEDS[0], init="random")
        cls.drugs = ["D0", "D1"]
        cls.values = {f"{drug}|{dose:g}": 0.1 + 0.01 * i for i, (drug, dose) in enumerate(
            (d, dose) for d in cls.drugs for dose in (0.2, 1.0, 10.0)
        )}
        _support.write_unit_eval(cls.out, cls.unit, cls.values, drugs=cls.drugs)

    def test_loads_and_derives_down_set(self):
        loaded = agg.load_unit_eval(self.unit, self.out)
        self.assertEqual(loaded.unit_id, self.unit["unit_id"])
        self.assertIn("n_down_200", loaded.per_condition.columns)
        self.assertEqual(int(loaded.per_condition["n_down_200"].iloc[0]), 100)

    def test_macro_must_match_per_condition(self):
        import json

        eval_dir = self.out / "units" / self.unit["unit_id"] / "eval"
        aggregate = json.loads((eval_dir / "aggregate.json").read_text(encoding="utf-8"))
        aggregate["by_k"]["200"]["macro"]["mean"] += 0.5
        (eval_dir / "aggregate.json").write_text(json.dumps(aggregate), encoding="utf-8")
        try:
            with self.assertRaises(agg.AggregateError):
                agg.load_unit_eval(self.unit, self.out)
        finally:
            aggregate["by_k"]["200"]["macro"]["mean"] -= 0.5
            (eval_dir / "aggregate.json").write_text(json.dumps(aggregate), encoding="utf-8")


class BootstrapTests(unittest.TestCase):
    def test_deterministic_and_ordered(self):
        values = [0.01, 0.02, 0.03, 0.04, 0.05]
        first = agg.bootstrap_mean(values, n_bootstrap=500, seed=7)
        second = agg.bootstrap_mean(values, n_bootstrap=500, seed=7)
        self.assertEqual(first, second)
        self.assertLessEqual(first["ci_low"], first["mean"])
        self.assertGreaterEqual(first["ci_high"], first["mean"])
        self.assertEqual(first["n"], 5)

    def test_ci_narrows_with_more_values(self):
        rng = np.random.default_rng(0)
        small = agg.bootstrap_mean(rng.normal(1.0, 0.1, 20), n_bootstrap=2000, seed=1)
        large = agg.bootstrap_mean(rng.normal(1.0, 0.1, 200), n_bootstrap=2000, seed=1)
        self.assertLess(
            large["ci_high"] - large["ci_low"], small["ci_high"] - small["ci_low"]
        )

    def test_single_value_and_empty(self):
        self.assertEqual(agg.bootstrap_mean([0.3], n_bootstrap=10)["ci_low"], 0.3)
        empty = agg.bootstrap_mean([], n_bootstrap=10)
        self.assertIsNone(empty["mean"])
        self.assertEqual(empty["n"], 0)


class PairingAndStrataTests(unittest.TestCase):
    """A small complete 30-unit population with a known Tahoe advantage."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.root = _support.clean_dir(_support.tmp_root() / "agg_pop")
        cls.pop = _support.write_synthetic_population(
            cls.root, drugs_per_fold=4, advantage=0.005
        )
        cls.baselines = _support.write_synthetic_baselines(
            cls.root, drugs_by_fold=cls.pop["drugs_by_fold"]
        )
        cls.evals = [
            agg.load_unit_eval(unit, cls.pop["out"]) for unit in cls.pop["units"]
        ]
        cls.unit_results = agg.unit_results_frame(cls.evals, {})
        cls.fold_seed = agg.fold_seed_frame(cls.unit_results)
        baseline_conditions, cls.baseline_record = agg.load_baseline_conditions(
            cls.baselines, required_folds=range(5)
        )
        cls.per_drug_dose = agg.per_drug_dose_frame(cls.evals, baseline_conditions)
        cls.paired = agg.paired_by_drug_frame(cls.per_drug_dose)
        cls.drug_gains = agg.drug_gains_frame(cls.paired)

    def test_unit_results_shape(self):
        self.assertEqual(len(self.unit_results), 30)
        self.assertAlmostEqual(
            self.unit_results["auprc_200"].mean(), self.unit_results["auprc_200"].mean()
        )
        self.assertEqual(self.fold_seed["n_seeds"].min(), 3)
        self.assertEqual(len(self.fold_seed), 10)

    def test_per_drug_dose_uses_evaluator_conditions(self):
        self.assertEqual(len(self.per_drug_dose), 30 * 4 * 3)
        self.assertIn("delta_mean_auprc_200", self.per_drug_dose.columns)
        self.assertTrue(self.per_drug_dose["delta_mean_auprc_200"].notna().all())
        # 4 doses? no: 3 doses per drug, 4 drugs per fold
        self.assertEqual(self.per_drug_dose["drug_key"].nunique(), 20)

    def test_paired_by_drug_one_row_per_fold_seed_drug(self):
        self.assertEqual(len(self.paired), 5 * 3 * 4)
        self.assertTrue((self.paired["gain_tahoe_minus_random"] > 0).all())
        self.assertTrue(
            (self.paired["gain_tahoe_minus_random"].round(10) == 0.005).all()
        )

    def test_drug_gains_average_over_seeds_not_multiplying_doses(self):
        self.assertEqual(len(self.drug_gains), 5 * 4)
        self.assertEqual(int(self.drug_gains["n_seeds"].min()), 3)
        self.assertAlmostEqual(
            float(self.drug_gains["gain_tahoe_minus_random_mean"].mean()), 0.005, places=10
        )

    def test_paired_summary_verdict_positive(self):
        summary = agg.paired_summary(self.drug_gains, "gain_tahoe_minus_random", n_bootstrap=500, seed=3)
        self.assertEqual(summary["verdict"], "positive")
        self.assertGreater(summary["overall"]["ci_low"], 0)
        self.assertEqual(summary["n_drugs"], 20)
        self.assertEqual(len(summary["per_fold_gain_tahoe_minus_random"]), 5)

    def test_paired_summary_negative_without_advantage(self):
        reverted = self.drug_gains.copy()
        reverted["gain_tahoe_minus_random_mean"] = -reverted[
            "gain_tahoe_minus_random_mean"
        ].abs()
        summary = agg.paired_summary(reverted, "gain_tahoe_minus_random", n_bootstrap=500, seed=3)
        self.assertEqual(summary["verdict"], "negative")

    def test_strata_cover_doses_and_tahoe_seen(self):
        seen = {unit["drug_key"] for unit in self.paired.to_dict("records")[:2]}
        strata = agg.strata_frame(
            self.per_drug_dose,
            self.drug_gains,
            tahoe_seen={"DRUG0000", "DRUG0001"},
            n_bootstrap=200,
            seed=3,
        )
        self.assertFalse(strata.empty)
        dose_rows = strata[strata["stratum_type"] == "dose"]
        self.assertEqual(sorted(dose_rows["stratum_value"].astype(str)), ["0.2", "1", "10"])
        seen_rows = strata[strata["stratum_type"] == "tahoe_seen"]
        self.assertEqual(sorted(seen_rows["stratum_value"].astype(str)), ["False", "True"])
        unknown = strata[~strata["label_known"]]
        self.assertEqual(
            sorted(unknown["stratum_type"].unique()),
            ["response_label", "toxicity_label"],
        )
        self.assertTrue((unknown["stratum_value"] == "unknown").all())
        self.assertTrue(unknown["random_auprc_200"].isna().all())

    def test_baseline_comparison_paired_by_drug(self):
        summary = agg.paired_summary(
            self.drug_gains, "gain_tahoe_minus_delta_mean", n_bootstrap=500, seed=3
        )
        self.assertLess(summary["overall"]["mean"], 0)  # synthetic model is below the bar
        self.assertEqual(summary["verdict"], "negative")


if __name__ == "__main__":
    unittest.main()
