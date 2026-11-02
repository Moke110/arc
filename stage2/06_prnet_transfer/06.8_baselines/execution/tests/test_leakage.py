"""Leakage: only the fold's TRAINING treated responses may influence a fit.

Attacks the boundary empirically rather than trusting the row-selection code:
poison every held-out (valid / test) treated well in a copy of the fixture, refit
the baselines on the same fold, and require the fitted parameters and the
predictions to be byte-identical.  A fit that read a held-out response — even to
compute a control mean — would change.

Also checks the row sets themselves: train ∩ test = ∅, and the fit table's rows
are exactly the fold's ``role=train`` treated rows.
"""

from __future__ import annotations

import unittest

import numpy as np

import _support
import baselines
from ptlib import data as data_mod


def _fit_all(frozen, fold_table, fold):
    structures = frozen.structures()
    table = baselines.build_fit_table(frozen, fold_table, fold, structures=structures)
    fits = {}
    for name in baselines.BASELINES:
        params, record = baselines.fit_baseline(name, table, structures, frozen)
        fits[name] = {
            "params_sha256": baselines.arrays_sha256(params),
            "fit": record,
            "params": params,
        }
    return table, fits


class LeakageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fold = 0
        cls.base_dir = _support.copy_fixture_to(
            _support.fixture_dir(), _support.tmp_root() / "leakage_base"
        )
        cls.frozen = _support.open_frozen(cls.base_dir)
        cls.folds = _support.folds_table(cls.frozen)
        cls.structures = cls.frozen.structures()
        cls.source = _support.source_for(cls.frozen, cls.fold)
        cls.train_rows, cls.valid_rows, cls.test_rows = _support.test_rows_for_fold(
            cls.base_dir, cls.fold
        )
        cls.base_table, cls.base_fits = _fit_all(cls.frozen, cls.folds, cls.fold)

    def test_train_rows_are_exactly_the_fold_train_role(self):
        np.testing.assert_array_equal(self.base_table.train_rows, self.train_rows)
        self.assertEqual(int(np.intersect1d(self.train_rows, self.test_rows).size), 0)
        self.assertEqual(int(np.intersect1d(self.train_rows, self.valid_rows).size), 0)

    def test_poisoning_test_targets_does_not_change_the_fit(self):
        poisoned_dir = _support.poison_cm_targets(
            self.base_dir, _support.tmp_root() / "leakage_poison_test", self.test_rows, delta=5.0
        )
        frozen = _support.open_frozen(poisoned_dir)
        table, fits = _fit_all(frozen, _support.folds_table(frozen), self.fold)

        self.assertEqual(table.train_delta_sha256, self.base_table.train_delta_sha256)
        for name in baselines.BASELINES:
            self.assertEqual(
                fits[name]["params_sha256"],
                self.base_fits[name]["params_sha256"],
                f"{name} changed after poisoning the test targets",
            )
        # and the predictions (target-free by construction) are unchanged too
        source = _support.source_for(frozen, self.fold)
        rows = np.asarray(source.rows_for("test"), dtype=np.int64)
        context = baselines.condition_context(frozen, source, self.structures, rows)
        base_context = baselines.condition_context(
            self.frozen, self.source, self.structures, rows
        )
        for name in baselines.BASELINES:
            predicted = baselines.predict_baseline(
                name, fits[name]["params"], table, context, frozen
            )
            base_predicted = baselines.predict_baseline(
                name,
                self.base_fits[name]["params"],
                self.base_table,
                base_context,
                self.frozen,
            )
            self.assertEqual(
                predicted.delta.tobytes(),
                base_predicted.delta.tobytes(),
                f"{name} prediction changed after poisoning the test targets",
            )

    def test_poisoning_valid_targets_does_not_change_the_fit(self):
        poisoned_dir = _support.poison_cm_targets(
            self.base_dir, _support.tmp_root() / "leakage_poison_valid", self.valid_rows, delta=3.0
        )
        frozen = _support.open_frozen(poisoned_dir)
        table, fits = _fit_all(frozen, _support.folds_table(frozen), self.fold)
        self.assertEqual(table.train_delta_sha256, self.base_table.train_delta_sha256)
        for name in baselines.BASELINES:
            self.assertEqual(
                fits[name]["params_sha256"],
                self.base_fits[name]["params_sha256"],
                f"{name} changed after poisoning the validation targets",
            )

    def test_fit_does_not_depend_on_test_control_means(self):
        # The truth Δ of a *training* well uses that well's plate DMSO mean; the
        # matched-control pool must be the QC-pass DMSO wells, not any treated row.
        obs = self.frozen.cm_obs.reset_index(drop=True)
        control_rows = np.flatnonzero(
            obs["is_control"].to_numpy(dtype=bool) & obs["qc_pass"].to_numpy(dtype=bool)
        ).astype(np.int64)
        self.assertGreater(control_rows.size, 0)
        self.assertEqual(int(np.intersect1d(control_rows, self.test_rows).size), 0)
        self.assertEqual(int(np.intersect1d(control_rows, self.valid_rows).size), 0)


if __name__ == "__main__":
    unittest.main()
