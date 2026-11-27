"""Verify the public full-axis scoring and pairing report seam."""

from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dose_report import summarize_scores


class PairingTests(unittest.TestCase):
    def test_reference_is_absolute_and_carps_contrasts_pair_fold_and_seed(self):
        roster = tuple((arm, 0, seed) for arm in ("TranSiGen", "no-dose", "Xmol-FiLM") for seed in (17, 42))
        scores = pd.DataFrame([
            dict(arm=arm, fold=fold, seed=seed, metric="Δ PCC", value=value, width=5000,
                 n_conditions=2, n_eligible=2, n_undefined=0)
            for (arm, fold, seed), value in zip(roster, [0.1, 0.2, 0.3, 0.8, 0.6, 0.7], strict=True)
        ])
        paired, _, summary = summarize_scores(scores, roster, metrics=("Δ PCC",))
        self.assertTrue(paired.loc[paired.arm.eq("TranSiGen"), "paired_difference"].isna().all())
        self.assertEqual(summary.set_index("arm").loc["TranSiGen", "n_paired"], 0)
        np.testing.assert_allclose(paired.loc[paired.arm.eq("Xmol-FiLM"), "paired_difference"], [0.3, -0.1])
        self.assertAlmostEqual(summary.set_index("arm").loc["Xmol-FiLM", "paired_mean"], 0.1)


if __name__ == "__main__":
    unittest.main()
