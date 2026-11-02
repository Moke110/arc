"""The evaluator's 口径 on known ranking situations — using the real evaluator.

The spec requires the pairing reference to be checkable before #9 trusts it.  On
the synthetic frozen fixture (real 18,164-gene axis, real frozen direction
AUPRC@200) this module writes prediction directories in the pipeline format and
scores them with ``ptlib.evaluate.run_evaluation`` — the exact code path #9 uses:

* a **perfect** Δ (the condition-mean truth Δ) ranks the top-200 up and the
  top-200 down genes first → ``auprc_up = auprc_down = auprc = 1``;
* a **constant** Δ (``no_change``) cannot rank → AP is the chance level
  ``k / n_genes`` for each direction (all scores tied);
* a **sign-flipped** Δ (``-truth``) ranks the true sets last → both directions
  fall below chance.

The evaluator's own guards are exercised too: a prediction directory built
against a different baseline than ``c*`` (or with a different row order) is
refused, which is what makes the #9 pairing safe.
"""

from __future__ import annotations

import json
import unittest
from pathlib import Path

import numpy as np

import _support
import baselines
from ptlib import data as data_mod
from ptlib import evaluate as evaluate_mod


class EvaluatorSanityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.frozen = _support.open_frozen()
        cls.folds = _support.folds_table(cls.frozen)
        cls.structures = cls.frozen.structures()
        cls.source = _support.source_for(cls.frozen, 0)
        cls.fold = 0
        cls.rows = np.asarray(cls.source.rows_for("test"), dtype=np.int64)
        cls.meta = cls.source.meta(cls.rows)
        cls.c_star = np.asarray(cls.source.c_star, dtype=np.float32)
        cls.k = int(cls.frozen.config["metric"]["k"])
        cls.n_genes = int(cls.frozen.n_genes)
        (
            cls.per_row_truth,
            cls.condition_truth,
            cls.condition_keys,
            cls.condition_inverse,
            cls.truth_info,
        ) = _support.condition_level_truth(cls.frozen, cls.source, cls.rows)
        cls.tmp = _support.tmp_root() / "evaluator-sanity"
        cls.tmp.mkdir(parents=True, exist_ok=True)

    def _score(self, name: str, delta: np.ndarray) -> dict:
        pred_dir = self.tmp / name / "pred"
        eval_dir = self.tmp / name / "eval"
        baselines.write_prediction_dir(
            pred_dir,
            baseline=name,
            role="primary",
            fold=self.fold,
            meta=self.meta,
            c_star=self.c_star,
            delta=delta,
            fit_record={"test_case": name},
            frozen=self.frozen,
            repeats=2,
        )
        evaluate_mod.run_evaluation(
            evaluate_mod.EvaluateOptions(
                stage="finetune",
                pred_dir=pred_dir,
                out_dir=eval_dir,
                data_dir=self.frozen.data_dir,
                role="test",
                fold=self.fold,
                command=["test_evaluator_sanity"],
            )
        )
        return json.loads((eval_dir / "aggregate.json").read_text(encoding="utf-8"))

    def _per_condition(self, name: str) -> np.ndarray:
        import pandas as pd

        frame = pd.read_csv(self.tmp / name / "eval" / "per_condition.csv")
        return frame[f"n_up_{self.k}"].to_numpy(dtype=np.int64)

    def _chance_levels(self, name: str) -> tuple[float, float]:
        """Chance-level AP per direction for tied scores: ``n_true_set / n_genes``.

        True-set sizes come from the evaluator's own ``deltas.npz`` so the
        fallback (fewer than k non-zero genes in a direction) is respected.
        """
        payload = np.load(self.tmp / name / "eval" / "deltas.npz")
        truth = np.asarray(payload["true_delta"], dtype=np.float64)
        n_pos = (truth > 0).sum(axis=1)
        n_neg = (truth < 0).sum(axis=1)
        up = np.where(n_pos >= self.k, self.k, n_pos)
        down = np.where(n_neg >= self.k, self.k, n_neg)
        return (
            float(np.mean(up / self.n_genes)),
            float(np.mean(down / self.n_genes)),
        )

    def test_perfect_prediction_ranks_both_directions(self):
        aggregate = self._score("perfect", self.per_row_truth.astype(np.float32))
        headline = aggregate["by_k"][str(self.k)]
        self.assertGreaterEqual(headline["up"]["mean"], 0.999)
        self.assertGreaterEqual(headline["down"]["mean"], 0.999)
        self.assertGreaterEqual(headline["macro"]["mean"], 0.999)
        # the fixture has at least k non-zero genes in each direction
        self.assertGreater(int(np.min(self._per_condition("perfect"))), 0)

    def test_constant_prediction_is_chance_level(self):
        aggregate = self._score(
            "constant", np.zeros((self.rows.size, self.n_genes), dtype=np.float32)
        )
        headline = aggregate["by_k"][str(self.k)]
        # chance for tied scores is n_true_set / n_genes per direction
        expected_up, expected_down = self._chance_levels("constant")
        self.assertAlmostEqual(headline["up"]["mean"], expected_up, places=9)
        self.assertAlmostEqual(headline["down"]["mean"], expected_down, places=9)
        self.assertAlmostEqual(
            headline["macro"]["mean"], float(np.mean([expected_up, expected_down])), places=9
        )
        # and the chance level is far below a useful ranking
        self.assertLess(headline["macro"]["mean"], 0.05)

    def test_sign_flipped_prediction_is_below_chance(self):
        aggregate = self._score(
            "sign_flipped", (-self.per_row_truth).astype(np.float32)
        )
        headline = aggregate["by_k"][str(self.k)]
        chance_up, chance_down = self._chance_levels("constant")
        self.assertLess(headline["up"]["mean"], chance_up)
        self.assertLess(headline["down"]["mean"], chance_down)
        self.assertLess(headline["macro"]["mean"], float(np.mean([chance_up, chance_down])))

    def test_delta_is_prediction_minus_c_star(self):
        # write_prediction_dir stores Δ = pred_expr − c*; the evaluator refuses a
        # mismatched baseline, so a wrong-baseline directory must be rejected.
        pred_dir = self.tmp / "wrong_baseline" / "pred"
        baselines.write_prediction_dir(
            pred_dir,
            baseline="wrong_baseline",
            role="primary",
            fold=self.fold,
            meta=self.meta,
            c_star=self.c_star,
            delta=np.zeros((self.rows.size, self.n_genes), dtype=np.float32),
            fit_record={},
            frozen=self.frozen,
        )
        np.save(
            pred_dir / "pred_expr.npy",
            (self.c_star.astype(np.float64)[None, :] + 0.5).astype(np.float32),
        )
        np.save(
            pred_dir / "delta.npy",
            np.full((self.rows.size, self.n_genes), 0.5, dtype=np.float32),
        )
        with self.assertRaises(data_mod.DataError):
            evaluate_mod.run_evaluation(
                evaluate_mod.EvaluateOptions(
                    stage="finetune",
                    pred_dir=pred_dir,
                    out_dir=self.tmp / "wrong_baseline" / "eval",
                    data_dir=self.frozen.data_dir,
                    role="test",
                    fold=self.fold,
                    command=["test_evaluator_sanity"],
                )
            )

    def test_row_order_mismatch_is_rejected(self):
        pred_dir = self.tmp / "bad_order" / "pred"
        baselines.write_prediction_dir(
            pred_dir,
            baseline="bad_order",
            role="primary",
            fold=self.fold,
            meta=self.meta,
            c_star=self.c_star,
            delta=np.zeros((self.rows.size, self.n_genes), dtype=np.float32),
            fit_record={},
            frozen=self.frozen,
        )
        meta = self.meta.copy()
        meta[["row"]] = meta[["row"]].iloc[::-1].to_numpy()
        meta.to_csv(pred_dir / "meta.csv", index=False)
        with self.assertRaises(data_mod.DataError):
            evaluate_mod.run_evaluation(
                evaluate_mod.EvaluateOptions(
                    stage="finetune",
                    pred_dir=pred_dir,
                    out_dir=self.tmp / "bad_order" / "eval",
                    data_dir=self.frozen.data_dir,
                    role="test",
                    fold=self.fold,
                    command=["test_evaluator_sanity"],
                )
            )


if __name__ == "__main__":
    unittest.main()
