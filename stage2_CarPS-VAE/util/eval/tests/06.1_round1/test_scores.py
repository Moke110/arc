import math
from pathlib import Path
import sys
import unittest

import pandas as pd

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
sys.path.insert(0, str(REPO / "stage2_CarPS-VAE/06_tsg_shape/06.1_round1/execution"))
from evaluate_shape import pair_scores


class PairScoresTest(unittest.TestCase):
    def test_differences_average_only_folds_with_two_defined_scores(self):
        fold_scores = pd.DataFrame([
            {"shape": shape, "fold": fold, "metric": "Δ PCC", "value": value,
             "n_conditions": 1, "n_eligible": int(math.isfinite(value)),
             "n_undefined": int(not math.isfinite(value))}
            for shape, values in (("original", (-0.25, 0.75, -0.1)),
                                  ("widen_z-130", (0.25, math.nan, 0.2)))
            for fold, value in enumerate(values)
        ])
        pairs, summary = pair_scores(fold_scores)
        rows = pairs.loc[pairs["shape"] == "widen_z-130"].set_index("fold")
        self.assertEqual(rows.loc[0, "score_difference"], 0.5)
        self.assertTrue(math.isnan(rows.loc[1, "score_difference"]))
        self.assertAlmostEqual(rows.loc[2, "score_difference"], 0.3)
        candidate = summary.set_index("shape").loc["widen_z-130"]
        self.assertAlmostEqual(candidate.value, 0.225)
        self.assertAlmostEqual(candidate.score_difference, 0.4)
        self.assertEqual(candidate.n_paired_folds, 2)
