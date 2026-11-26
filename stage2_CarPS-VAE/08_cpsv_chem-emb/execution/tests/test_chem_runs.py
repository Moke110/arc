"""Test frozen coordinate and chemical input boundaries."""

from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from chem_runs import aligned_features, coordinates


class VectorTable:
    """Representations in a source order unlike the condition order."""

    def vectors(self, arm, parents):
        width = {"KPGT": 2304, "E3FP": 1024}[arm]
        return np.broadcast_to(np.array([{"a": 11, "b": 22}[parent] for parent in parents],
                                        dtype=np.float32)[:, None], (len(parents), width))


class ChemicalRunsTest(unittest.TestCase):
    def test_all_arms_use_the_frozen_formal_coordinates(self):
        roster = coordinates()
        self.assertEqual(len(set(roster)), 63)
        arms = {"ECFP4C", "ECFP6C", "FCFP4C", "FCFP6C", "E3FP", "KPGT", "UniMol"}
        self.assertEqual(set(roster), {(arm, fold, seed) for arm in arms for fold in range(3)
                                     for seed in (17, 42, 2026)})
        self.assertEqual(set(coordinates(smoke=True)), {(arm, 0, 17) for arm in arms})

    def test_native_width_vectors_follow_condition_order(self):
        features = aligned_features(VectorTable(), "E3FP", np.array(["b", "a", "b"]))
        self.assertEqual(features.shape, (3, 1024))
        np.testing.assert_array_equal(features[:, 0], [22, 11, 22])
        self.assertEqual(features.dtype, np.float32)


if __name__ == "__main__":
    unittest.main()
