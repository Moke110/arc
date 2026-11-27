"""Verify the Step 10 matrix uses the frozen crossed tenv coordinates."""

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dose_protocol import ARMS, coordinates


class CoordinateTests(unittest.TestCase):
    def test_complete_formal_roster(self):
        roster = coordinates()
        self.assertEqual(len(roster), 270)
        self.assertEqual(len(set(roster)), 270)
        self.assertEqual([entry[1:] for entry in roster if entry[0] == ARMS[0]],
                         [(fold, seed) for fold in range(9) for seed in (17, 42, 2026)])

    def test_bounded_roster_exercises_a_new_crossed_fold(self):
        self.assertEqual(coordinates(bounded=True), tuple((arm, 1, 17) for arm in ARMS))
