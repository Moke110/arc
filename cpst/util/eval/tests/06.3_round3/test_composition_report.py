"""Test the Step 06.3 paired report from independently specified scores."""

from pathlib import Path
import sys
import unittest

import pandas as pd

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
sys.path.insert(0, str(REPO / "cpst/06_tsg_shape/06.3_round3/execution"))
from composition_scores import paired_report


SHAPES = ("0_original", "1_output-delta", "2_output-delta-no-x2",
          "3_z-interact_FiLM", "4_1+3", "5_2+3")
METRICS = ("Δ PCC", "Δ RMSE")


class PairedReportTest(unittest.TestCase):
    def setUp(self):
        self.expected = [(shape, fold, seed) for shape in SHAPES for fold in (0, 1) for seed in (17, 42)]
        self.scores = pd.DataFrame([
            {"shape": shape, "fold": fold, "seed": seed, "metric": metric,
             "value": 4 * fold + (seed == 42) + (1, 3, 2, 4, 8, 7)[SHAPES.index(shape)],
             "n_conditions": 1, "n_eligible": 1, "n_undefined": 0}
            for shape, fold, seed in reversed(self.expected)
            for metric in METRICS
        ])

    def test_pairing_uses_fold_and_seed_and_reports_all_contrasts(self):
        self.scores.loc[self.scores["shape"].eq("4_1+3") & self.scores.seed.eq(17), "value"] += 2
        pairs, summary, contrasts, contrast_summary = paired_report(self.scores, self.expected, METRICS)
        combo = pairs.query("shape == '4_1+3' and metric == 'Δ PCC'")
        self.assertEqual(combo.control_score.tolist(), [1, 2, 5, 6])
        self.assertEqual(combo.score_difference.tolist(), [9, 7, 9, 7])
        self.assertEqual(summary.query("shape == '4_1+3' and metric == 'Δ PCC'").score_difference.item(), 8)
        one_coordinate = contrasts.query("fold == 1 and seed == 42 and metric == 'Δ RMSE'")
        self.assertEqual(one_coordinate.set_index("contrast").value.to_dict(),
                         {"4−1": 5, "4−3": 4, "5−2": 5, "5−3": 3, "4−1−3+0": 2, "5−2−3+0": 2})
        interaction = contrast_summary.query("contrast == '4−1−3+0' and metric == 'Δ PCC'")
        self.assertEqual(interaction.value.item(), 3)
        self.assertAlmostEqual(interaction.sample_sd.item(), 1.1547005383792515)
        self.assertEqual(interaction.n_defined.item(), 4)

    def test_missing_duplicate_and_mismatched_coordinates_fail(self):
        invalid = [
            self.scores.iloc[1:],
            pd.concat([self.scores, self.scores.iloc[[0]]], ignore_index=True),
            self.scores.assign(seed=self.scores.seed.where(self.scores.seed != 42, 43)),
            self.scores.assign(fold=self.scores.fold.where(self.scores.fold != 1, 2)),
        ]
        for scores in invalid:
            with self.subTest(rows=len(scores), seeds=set(scores.seed), folds=set(scores.fold)):
                with self.assertRaises(ValueError):
                    paired_report(scores, self.expected, METRICS)

    def test_undefined_scores_do_not_create_secondary_estimates(self):
        missing = (self.scores["shape"].eq("3_z-interact_FiLM") & self.scores.seed.eq(17)
                   & self.scores.fold.eq(0) & self.scores.metric.eq("Δ PCC"))
        self.scores.loc[missing, ["value", "n_eligible", "n_undefined"]] = [float("nan"), 0, 1]
        _, _, _, summary = paired_report(self.scores, self.expected, METRICS)
        self.assertEqual(summary.query("contrast == '4−3' and metric == 'Δ PCC'").n_defined.item(), 3)
        self.assertEqual(summary.query("contrast == '4−1' and metric == 'Δ PCC'").n_defined.item(), 4)
