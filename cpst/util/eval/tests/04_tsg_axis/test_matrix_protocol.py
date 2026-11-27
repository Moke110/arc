"""Exercise the Step 04 matrix's public registration gate."""

import json
from pathlib import Path
import sys
import tempfile
import unittest

STEP = Path(__file__).resolve().parents[4] / "04_tsg_axis"
sys.path.insert(0, str(STEP / "execution"))

from matrix import check_registrations, run_paths
from tsg_axis.protocol import MAX_EPOCHS, PATIENCE, PROTOCOL, RUN_SCHEMA, shared_protocol


class MatrixProtocolTests(unittest.TestCase):
    def test_roster_and_registration_seed_must_agree(self):
        protocol = shared_protocol()
        with tempfile.TemporaryDirectory() as scratch:
            paths = run_paths(Path(scratch))
            self.assertEqual(len(paths), 150)
            key = "hvg-5k", 4, 2027
            target = paths[key]
            self.assertEqual(target.relative_to(scratch).as_posix(), "runs/hvg-5k/fold-4/seed-2027")
            target.mkdir(parents=True)
            record = {
                "arguments": {"axis": "hvg-5k", "fold": 4, "seed": 42, "smoke": False},
                "seed": 2027, "scenario": "cold-drug", "schema": RUN_SCHEMA, "protocol": PROTOCOL,
                "conditions": {}, "max_epochs": MAX_EPOCHS, "patience": PATIENCE,
                "training_sha256": "training", "folds": "legacy", "nonce": "one",
                "fold_release_sha256": protocol.release_sha256, "code_sha256": {},
            }
            config = target / "registration.json"
            evidence = {"training_manifest_sha256": {"hvg-5k": "training"}, "folds": "legacy",
                        "fold_release_sha256": protocol.release_sha256}
            config.write_text(json.dumps(record))
            with self.assertRaisesRegex(ValueError, "registration plan"):
                check_registrations({key: target}, [{}] * 5, evidence)
            record["arguments"]["seed"] = 2027
            config.write_text(json.dumps(record))
            self.assertEqual(set(check_registrations({key: target}, [{}] * 5, evidence)), {key})
