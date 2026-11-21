import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "axis_local"))
from control_reliability import paired_control_scores


class ControlReliabilityTest(unittest.TestCase):
    def test_same_wells_and_equal_size_shared_or_disjoint_vehicle_pools(self):
        first = np.array([2., -2., 0., 0.])
        second = np.array([0., 0., 2., -2.])
        controls = np.array([[0., 0., 0., 0.], [1., -1., 1., -1.]])
        shared, independent = paired_control_scores(first, second, controls)
        self.assertAlmostEqual(shared, -.5)
        self.assertAlmostEqual(independent, -1 / np.sqrt(2))
