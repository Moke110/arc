import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'axis_local'))
from baselines import mean_cell


class BaselineTest(unittest.TestCase):
    def test_mean_uses_only_matching_fit_rows_and_unique_vehicles(self):
        data = SimpleNamespace(
            obs=pd.DataFrame({'cell_line': ['a', 'a', 'b', 'a']}),
            expression=np.array([[2.], [6.], [900.], [800.], [10.], [14.], [700.]], np.float32),
            controls={0: [4, 5], 1: [4], 2: [6], 3: [6]},
        )
        predictions = mean_cell(data, np.array([0, 1, 2]), np.array([3]))
        np.testing.assert_array_equal(predictions, [[8.]])
        with self.assertRaisesRegex(ValueError, 'training cell coverage'):
            mean_cell(data, np.array([2]), np.array([3]))
