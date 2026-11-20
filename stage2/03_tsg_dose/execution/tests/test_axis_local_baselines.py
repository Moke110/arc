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


class BaselineReportTest(unittest.TestCase):
    def test_axis_break_keeps_fold_spreads_and_dose_trajectories_intact(self):
        import runpy

        figure = Path(__file__).resolve().parents[2] / 'figure'
        successor = runpy.run_path(str(figure / 'baselines/script/render.py'))
        base = successor['load_reference'](figure / 'script/render_axis_local.py')
        choose = successor['separated_metric_axis']
        axis = choose(base, [(-13, -9), (-0.25, .15), (.15, .20)])
        self.assertEqual(len(axis.segments), 2)
        self.assertLess(axis.segments[0]['upper'], -8)
        self.assertGreater(axis.segments[1]['lower'], -1)
        self.assertTrue(all(np.isfinite(axis.map([-13, -9, -.25, 0, .2]))))
        self.assertIsNone(choose(base, [(-13, -9), (-10, 0), (-.25, .2)]))
