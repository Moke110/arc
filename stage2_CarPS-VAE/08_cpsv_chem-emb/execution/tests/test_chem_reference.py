"""Exercise the exact KPGT score gate used before reference fit reuse."""

from pathlib import Path
import sys
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from chem_reference import verify_reference_scores
from chem_runs import ORIGINAL_CODE, axis_run
from util.training.provenance import git_code_hashes, verify_replay_configuration
from chem_report import METRICS


class ReferenceTest(unittest.TestCase):
    def test_sealed_runtime_pin_keeps_the_scientific_registration_exact(self):
        path = "stage2_CarPS-VAE/08_cpsv_chem-emb/execution/chem_runs.py"
        pins = git_code_hashes(axis_run.REPO, ORIGINAL_CODE, (path,))
        current = {"code_sha256": {path: "new-code"}, "experiment": {"fold": 0}, "partitions": {"test": [8]}}
        sealed = {"code_sha256": pins, "partitions": {"test": [8]}}
        options = {"repo": axis_run.REPO, "revision": ORIGINAL_CODE}

        verify_replay_configuration(sealed, current, replay=True, **options)
        for invalid, replay in ((sealed, False), ({**sealed, "partitions": {"test": [9]}}, True),
                                ({**sealed, "code_sha256": {path: "wrong-code"}}, True)):
            with self.assertRaisesRegex(ValueError, "registered study"):
                verify_replay_configuration(invalid, current, replay=replay, **options)

    def test_changed_score_rejects_reuse(self):
        recorded = pd.DataFrame({
            "metric": METRICS, "value": [.3] * len(METRICS),
            "n_conditions": 7, "n_eligible": 7, "n_undefined": 0,
            "run_sha256": "original-registration",
        })
        independently_scored = recorded.drop(columns="run_sha256").copy()
        independently_scored.loc[independently_scored.metric.eq("Δ PCC"), "value"] = .30000000000000004
        with self.assertRaisesRegex(ValueError, "exactly"):
            verify_reference_scores(independently_scored, recorded, "original-registration")


if __name__ == "__main__":
    unittest.main()
