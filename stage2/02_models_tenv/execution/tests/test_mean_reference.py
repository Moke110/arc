"""The deterministic ``Mean`` reference (#165).

The retired ``delta_mean`` baseline is the authority for the exact arithmetic:
one global prediction vector, the equal-weight mean over training-partition
drug-dose groups of each group's mean Δ, copied to every test well.  A cold-drug
test group is absent from the training partition by construction and receives
the same vector — there is no per-group lookup and no fallback.
"""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

STEP_ROOT = Path(__file__).resolve().parents[1].parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "script"))
import mean_reference  # noqa: E402
from evaluation_driver import validate_run  # noqa: E402


class GroupMeanTests(unittest.TestCase):
    def test_equal_weight_group_mean_matches_the_hand_computed_value(self):
        rows = [{"sample_id": "a1", "group_id": "A"},
                {"sample_id": "a2", "group_id": "A"},
                {"sample_id": "b1", "group_id": "B"}]
        responses = np.array([[1.0, 2.0], [3.0, 4.0], [10.0, 20.0]])
        prediction = mean_reference.group_mean_prediction(rows, responses)
        # group A mean [2, 3], group B mean [10, 20]; equal weight => [6, 11.5]
        np.testing.assert_allclose(prediction, np.array([6.0, 11.5], dtype=np.float32))

    def test_cold_drug_groups_receive_the_same_global_reference(self):
        training = [{"sample_id": f"w{index}", "group_id": "seen"} for index in range(3)]
        responses = np.tile(np.array([1.0, 3.0, 5.0]), (3, 1)).T
        responses = np.repeat(responses, 18164 // 3 + 1, axis=1)[:, :18164]
        delta = mean_reference.group_mean_prediction(training, responses)
        self.assertEqual(delta.shape, (18164,))
        # An unseen (cold-drug) test group is never looked up; every test well is
        # copied the global vector.
        cases = ["test_seen", "test_cold_drug", "test_cold_drug_2"]
        with tempfile.TemporaryDirectory() as tmp:
            mean_reference.publish(Path(tmp), fold=0, seed=17, case_ids=cases,
                                   prediction=np.broadcast_to(delta, (3, 18164)),
                                   source={"fixture": True})
            prediction = np.load(Path(tmp) / "prediction.npy")
            self.assertEqual(prediction.shape, (3, 18164))
            np.testing.assert_array_equal(prediction[0], prediction[1])
            np.testing.assert_array_equal(prediction[1], prediction[2])
            np.testing.assert_allclose(prediction[0], delta)

    def test_empty_training_groups_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "training"):
            mean_reference.group_mean_prediction([], np.zeros((0, 18164)))


class MeanRunArtifactTests(unittest.TestCase):
    def run_dir(self, directory, cases):
        delta = np.linspace(-1.0, 1.0, 18164, dtype=np.float32)
        mean_reference.publish(directory, fold=1, seed=42, case_ids=cases,
                               prediction=np.broadcast_to(delta, (len(cases), 18164)),
                               source={"fixture": True})
        return Path(directory)

    def test_published_run_is_complete_and_validated_by_the_shared_driver(self):
        cases = ["s1", "s2", "s3"]
        with tempfile.TemporaryDirectory() as tmp:
            run = self.run_dir(Path(tmp) / "fold_1/seed_42", cases)
            artifacts = validate_run(run, cases)
            self.assertEqual(artifacts.terminal["status"], "complete")
            self.assertEqual(artifacts.label, "Mean")
            self.assertEqual(artifacts.fold, 1)
            self.assertEqual(artifacts.seed, 42)
            self.assertEqual(artifacts.history, [])
            identity = mean_reference.read_json(run / "identity.json")
            self.assertTrue(identity["identity"]["deterministic_reference"])
            self.assertIn("cold_drug_semantics", identity["identity"])

    def test_republishing_the_same_reference_is_byte_deterministic(self):
        cases = ["s1", "s2"]
        delta = np.full(18164, 0.25, dtype=np.float32)
        with tempfile.TemporaryDirectory() as tmp:
            first = Path(tmp) / "first"
            second = Path(tmp) / "second"
            for directory in (first, second):
                mean_reference.publish(directory, fold=0, seed=17, case_ids=cases,
                                       prediction=np.broadcast_to(delta, (2, 18164)),
                                       source={"fixture": True})
            for name in ("prediction.npy", "case_ids.json", "terminal.json"):
                self.assertEqual((first / name).read_bytes(),
                                 (second / name).read_bytes(), name)


if __name__ == "__main__":
    unittest.main()
