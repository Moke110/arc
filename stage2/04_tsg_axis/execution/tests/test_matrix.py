"""The formal matrix preflight and registration boundary."""
import sys
import unittest
from pathlib import Path

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP / "execution"))

from matrix import match_rosters


class FormalMatrixTest(unittest.TestCase):
    def test_rejects_reordered_held_out_rows_before_registration(self):
        fold = {"train": ["a"], "valid": ["b"], "test": ["c", "d"]}
        self.assertEqual(match_rosters({"axis-1": [fold], "axis-2": [fold]}), [fold])
        changed = {**fold, "test": ["d", "c"]}
        with self.assertRaisesRegex(ValueError, "axis-2.*fold 0.*test"):
            match_rosters({"axis-1": [fold], "axis-2": [changed]})
