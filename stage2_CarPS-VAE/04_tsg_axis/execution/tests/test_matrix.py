"""The formal matrix preflight and registration boundary."""
import sys
import unittest
import json
from pathlib import Path
from tempfile import TemporaryDirectory

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP / "execution"))

from matrix import check_registrations, match_rosters


class FormalMatrixTest(unittest.TestCase):
    def test_rejects_reordered_held_out_rows_before_registration(self):
        fold = {"train": ["a"], "valid": ["b"], "test": ["c", "d"]}
        self.assertEqual(match_rosters({"axis-1": [fold], "axis-2": [fold]}), [fold])
        changed = {**fold, "test": ["d", "c"]}
        with self.assertRaisesRegex(ValueError, "axis-2.*fold 0.*test"):
            match_rosters({"axis-1": [fold], "axis-2": [changed]})

    def test_formal_registration_requires_matched_control_protocol_and_patience(self):
        axis = "hvg-5k"
        roster = {"train": ["a"], "valid": ["b"], "test": ["c"]}
        preflight = {"training_manifest_sha256": {axis: "training"}, "folds": "folds"}
        config = {
            "arguments": {"axis": axis, "fold": 0, "smoke": False},
            "schema": "stage2-tenv-transigen-axis-run-v2",
            "protocol": "matched-plate-control-native-loss-v1",
            "scenario": "cold-drug", "seed": 42, "conditions": roster,
            "max_epochs": 1000, "patience": 30, "training_sha256": "training",
            "folds": "folds", "nonce": "fresh", "code_sha256": {"run.py": "successor"},
        }
        with TemporaryDirectory() as scratch:
            path = Path(scratch)
            identity = path / "registration.json"
            identity.write_text(json.dumps(config))
            self.assertEqual(set(check_registrations({(axis, 0): path}, [roster], preflight)), {(axis, 0)})
            for change in ({"schema": "stage2-tenv-transigen-axis-run-v1"},
                           {"protocol": "axis-wide-control"}, {"patience": 50}):
                with self.subTest(change=change):
                    identity.write_text(json.dumps(config | change))
                    with self.assertRaisesRegex(ValueError, "fresh formal registration plan"):
                        check_registrations({(axis, 0): path}, [roster], preflight)
