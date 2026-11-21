"""Behavior of the bounded development analysis."""
from pathlib import Path
from types import SimpleNamespace
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'output/api'))
from tsg_transfer.development import delta_scores, split_validation


class DevelopmentTests(unittest.TestCase):
    def test_validation_keeps_all_profiles_of_a_parent_together(self):
        data = SimpleNamespace(obs=pd.DataFrame(
            {'drug_parent_inchikey': list('AABBCCDD')},
        ))
        roles = split_validation(data, [0, 1], [2, 3, 4, 5, 6, 7])
        selector = set(data.obs.iloc[roles['selection']].drug_parent_inchikey)
        assessment = set(data.obs.iloc[roles['assessment']].drug_parent_inchikey)
        self.assertTrue(selector and assessment)
        self.assertFalse(selector & assessment)
        self.assertEqual(selector | assessment, {'B', 'C', 'D'})

    def test_delta_scores_use_matched_control_and_preserve_undefined_pcc(self):
        truth = [[3, 4, 6], [3, 3, 3]]
        prediction = [[5, 4, 4], [1, 2, 3]]
        control = [[1, 2, 3], [2, 2, 2]]
        pcc, rmse = delta_scores(truth, prediction, control)
        self.assertAlmostEqual(pcc[0], -4 / (2 * np.sqrt(7)))
        self.assertTrue(np.isnan(pcc[1]))
        np.testing.assert_allclose(rmse, np.sqrt([8 / 3, 5 / 3]))
