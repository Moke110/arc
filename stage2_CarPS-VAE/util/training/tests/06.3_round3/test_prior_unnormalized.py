"""Require the complete historical non-FiLM coordinate set before reuse."""

from itertools import product
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str((Path(__file__).resolve().parents[4] / "06_tsg_shape/06.3_round3/execution/tests" / Path(__file__).name).parents[1]))
from prior_unnormalized import select_reusable


class PriorUnnormalizedTest(unittest.TestCase):
    def setUp(self):
        self.retained = list(product(
            ("0_original", "1_output-delta", "2_output-delta-no-x2"), range(3), (17, 2026, 43),
        ))
        self.other = ("3_z-interact_FiLM", 0, 17)
        self.completed = [{"coordinate": list(key)} for key in [*self.retained, self.other]]

    def test_only_unchanged_arms_are_reused(self):
        self.assertEqual(list(select_reusable(self.completed)), self.retained)

    def test_missing_retained_or_duplicate_historical_coordinates_fail(self):
        for invalid in (self.completed[1:], [*self.completed, self.completed[-1]]):
            with self.subTest(rows=len(invalid)), self.assertRaises(ValueError):
                select_reusable(invalid)
