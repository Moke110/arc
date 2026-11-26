"""Exercise the exact KPGT score gate used before reference fit reuse."""

import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from chem_reference import verify_live_alpha_dependencies, verify_reference_scores
from chem_report import METRICS


class ReferenceTest(unittest.TestCase):
    def test_live_alpha_dependencies_match_the_frozen_release(self):
        stage = Path(__file__).resolve().parents[3]
        release = json.loads((stage / "07_cpsv_alpha/output/release.json").read_text())
        scratch = stage.parent / "tmp"
        scratch.mkdir(exist_ok=True)
        dependencies = (
            "05_tsg_loss/05.1_loss_metrics/execution/objective.py",
            "05_tsg_loss/05.2_mse_l1/execution/additive_objective.py",
            "05_tsg_loss/05.2_mse_l1/execution/run_sweep.py",
        )
        with TemporaryDirectory(dir=scratch) as work:
            candidate = Path(work)
            for name in dependencies:
                path = candidate / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes((stage / name).read_bytes())
            verify_live_alpha_dependencies(candidate, release)
            for name in dependencies:
                with self.subTest(dependency=name):
                    path = candidate / name
                    original = path.read_bytes()
                    path.write_bytes(original + b"\n# changed\n")
                    with self.assertRaisesRegex(ValueError, path.name):
                        verify_live_alpha_dependencies(candidate, release)
                    path.write_bytes(original)

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
