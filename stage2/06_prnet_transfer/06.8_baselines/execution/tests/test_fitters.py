"""Pure fitter/predictor logic (no frozen input): pooling rules, fallbacks, encoding.

These tests pin the *definitions* that #9 pairs against:

* ``delta_mean`` pools equal weight per training condition (replicate aggregation
  happened before pooling), uses an exact condition lookup when available, and
  otherwise the documented pooled fallback;
* ``delta_mean_dose`` is a dose pool with the pooled-mean fallback;
* ``ridge_fp`` uses the study's own dose encoding (rFCFP4 × log10(dose+1)) and a
  fixed a-priori alpha, and its prediction is ``X·coef + intercept``.
"""

from __future__ import annotations

import types
import unittest

import numpy as np

import _support  # noqa: F401  (sets sys.path)
import baselines


def make_table(condition_delta, drugs, doses, fold: int = 0) -> baselines.FitTable:
    condition_delta = np.asarray(condition_delta, dtype=np.float64)
    drugs = np.asarray(drugs, dtype=object)
    doses = np.asarray(doses, dtype=np.float64)
    keys = np.asarray([f"{d}|{z:g}" for d, z in zip(drugs, doses)], dtype=object)
    n_rows = condition_delta.shape[0]
    return baselines.FitTable(
        fold=fold,
        n_genes=int(condition_delta.shape[1]),
        train_rows=np.arange(n_rows, dtype=np.int64),
        row_conditions=np.arange(n_rows, dtype=np.int64),
        unique_conditions=keys,
        condition_delta=condition_delta,
        condition_drug=drugs,
        condition_dose=doses,
        pooled_mean=condition_delta.mean(axis=0),
        counts={"n_train_rows": n_rows, "n_train_conditions": n_rows},
        train_rows_sha256="train-rows-sha",
        train_delta_sha256="train-delta-sha",
    )


def make_context(drugs, doses, smiles) -> dict:
    drugs = np.asarray(drugs, dtype=object)
    doses = np.asarray(doses, dtype=np.float64)
    smiles = np.asarray(smiles, dtype=object)
    keys = np.asarray([f"{d}|{z:g}" for d, z in zip(drugs, doses)], dtype=object)
    return {
        "rows": np.arange(len(drugs), dtype=np.int64),
        "condition_keys": keys,
        "row_conditions": np.arange(len(drugs), dtype=np.int64),
        "condition_drug": drugs,
        "condition_dose": doses,
        "condition_smiles": smiles,
    }


def fake_frozen(n_genes: int, drug_dim: int = 16):
    return types.SimpleNamespace(n_genes=n_genes, structure={"drug_dim": drug_dim})


def fake_structures(mapping):
    return types.SimpleNamespace(smiles_for=lambda key: mapping.get(str(key), ""))


class NoChangeTests(unittest.TestCase):
    def test_no_change_is_zero_delta(self):
        table = make_table([[1.0, -1.0]], ["A"], [1.0])
        params, record = baselines.fit_baseline("no_change", table, None, fake_frozen(2))
        self.assertEqual(record["role"], "primary")
        np.testing.assert_array_equal(params["delta"], np.zeros(2, dtype=np.float32))
        result = baselines.predict_baseline(
            "no_change", params, table, make_context(["A"], [1.0], ["C"]), fake_frozen(2)
        )
        np.testing.assert_array_equal(result.delta, np.zeros((1, 2), dtype=np.float32))


class DeltaMeanTests(unittest.TestCase):
    def test_pooling_is_equal_weight_per_condition(self):
        # Condition means (not raw rows): the 2-condition drug and the 1-condition
        # drug must contribute equally to the pooled vector.
        table = make_table(
            condition_delta=[[1.0, 2.0], [3.0, 4.0], [10.0, 20.0]],
            drugs=["A", "A", "B"],
            doses=[1.0, 1.0, 0.2],
        )
        params, record = baselines.fit_baseline("delta_mean", table, None, fake_frozen(2))
        np.testing.assert_allclose(
            params["pooled"], np.array([14.0 / 3.0, 26.0 / 3.0]), rtol=1e-6
        )
        self.assertIn("equal weight per training condition", record["pooling"])
        self.assertEqual(record["role"], "primary")

    def test_lookup_first_then_documented_fallback(self):
        table = make_table(
            condition_delta=[[1.0, 0.0], [0.0, 2.0]],
            drugs=["A", "B"],
            doses=[1.0, 1.0],
        )
        params, _ = baselines.fit_baseline("delta_mean", table, None, fake_frozen(2))
        context = make_context(["A", "C"], [1.0, 1.0], ["C", "CC"])
        result = baselines.predict_baseline("delta_mean", params, table, context, fake_frozen(2))
        np.testing.assert_allclose(result.delta[0], table.condition_delta[0].astype(np.float32))
        np.testing.assert_allclose(result.delta[1], params["pooled"])
        self.assertEqual(result.info["n_conditions_seen_in_train"], 1)


class DeltaMeanDoseTests(unittest.TestCase):
    def test_dose_pool_and_unseen_dose_fallback(self):
        table = make_table(
            condition_delta=[[1.0], [3.0], [10.0]],
            drugs=["A", "B", "C"],
            doses=[0.2, 0.2, 1.0],
        )
        params, record = baselines.fit_baseline("delta_mean_dose", table, None, fake_frozen(1))
        self.assertEqual(record["role"], "auxiliary")
        self.assertEqual(sorted(record["doses"]), ["0.2", "1"])
        context = make_context(["X", "Y"], [0.2, 10.0], ["C", "CC"])
        result = baselines.predict_baseline(
            "delta_mean_dose", params, table, context, fake_frozen(1)
        )
        np.testing.assert_allclose(result.delta[0], np.array([2.0], dtype=np.float32))
        np.testing.assert_allclose(result.delta[1], params["pooled"])
        self.assertEqual(result.info["n_conditions_with_unseen_dose"], 1)


class RidgeFeatureTests(unittest.TestCase):
    def test_features_are_fingerprint_times_log_dose(self):
        features = baselines.ridge_features(["C", "C", "C"], [0.0, 9.0, 99.0], 16)
        self.assertEqual(features.shape, (3, 16))
        np.testing.assert_allclose(features[0], np.zeros(16))
        # log10(9+1) / log10(99+1) = 0.5 for the same molecule.
        np.testing.assert_allclose(features[1], 0.5 * features[2], rtol=1e-12)

    def test_alpha_is_fixed_a_priori(self):
        self.assertEqual(baselines.RIDGE_ALPHA, 1.0)


class RidgeFitTests(unittest.TestCase):
    def test_zero_dose_fit_predicts_the_intercept(self):
        n_genes, n_bits = 3, 16
        # All conditions at dose 0 -> the feature matrix is all zeros, so the
        # ridge coefficients are zero and the prediction is the training mean.
        table = make_table(
            condition_delta=[[1.0, 2.0, 3.0], [3.0, 4.0, 5.0]],
            drugs=["C", "CC"],
            doses=[0.0, 0.0],
        )
        frozen = fake_frozen(n_genes, n_bits)
        structures = fake_structures({"C": "C", "CC": "CC"})
        params, record = baselines.fit_baseline("ridge_fp", table, structures, frozen)
        self.assertEqual(record["role"], "auxiliary")
        self.assertEqual(record["alpha"], baselines.RIDGE_ALPHA)
        self.assertEqual(record["n_bits"], n_bits)
        np.testing.assert_allclose(params["coef"], np.zeros_like(params["coef"]), atol=1e-8)
        np.testing.assert_allclose(params["intercept"], table.pooled_mean.astype(np.float32))
        result = baselines.predict_baseline(
            "ridge_fp", params, table, make_context(["C"], [0.0], ["C"]), frozen
        )
        np.testing.assert_allclose(result.delta[0], table.pooled_mean.astype(np.float32))


class HashAndErrorTests(unittest.TestCase):
    def test_params_hash_is_content_based_and_stable(self):
        a = {"x": np.arange(4, dtype=np.float32)}
        b = {"x": np.arange(4, dtype=np.float32)}
        self.assertEqual(baselines.arrays_sha256(a), baselines.arrays_sha256(b))
        b["x"][0] = 99
        self.assertNotEqual(baselines.arrays_sha256(a), baselines.arrays_sha256(b))

    def test_unknown_baseline_rejected(self):
        table = make_table([[0.0, 1.0]], ["C"], [0.0])
        with self.assertRaises(baselines.BaselineError):
            baselines.fit_baseline("nope", table, None, fake_frozen(2))


if __name__ == "__main__":
    unittest.main()
