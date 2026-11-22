"""Run seam: a bounded smoke fit trains, replays and scores in fresh processes."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXECUTION = HERE.parent
for path in (str(EXECUTION), str(HERE)):
    if path not in sys.path:
        sys.path.insert(0, path)

import build  # noqa: E402
import freeze  # noqa: E402
import smoke  # noqa: E402
import tsg_view_fixtures as fixtures  # noqa: E402
import views  # noqa: E402

REPO = Path(__file__).resolve().parents[4]
#: Repo-relative (or TSG_NATIVE_SOURCE) discovery; no checkout is hard-coded.
SOURCE = fixtures.native_source(REPO)
WORKER = EXECUTION / "learner.py"
SIZE = 5000


def run(*arguments):
    return subprocess.run([sys.executable, *map(str, arguments)], capture_output=True,
                          text=True)


@unittest.skipUnless((SOURCE / "src/model.py").is_file(),
                     fixtures.dependency_message("the TranSiGen source", SOURCE))
class LearnerSmokeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        views.bind(REPO)

    @classmethod
    def setUpClassPaths(cls):
        pass

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        self.fixture = fixtures.dense_fixture(self.root)
        self.accepted, _, _ = fixtures.build_view(self.fixture, size=SIZE)
        population = fixtures.treated_population(self.accepted, "tenv")
        self.sealed = fixtures.write_sealed(self.root, "tenv", population)
        build.build_product(view="full", source="tenv", output=self.root / "views",
                            staging=self.root / "staging", repo_root_path=REPO, size=SIZE,
                            fixture=True, corpus=self.fixture["corpus"],
                            axis=self.fixture["axis"], accepted=self.accepted,
                            provenance_root=None)
        self.freeze = self.root / "freezes/full/tenv"
        freeze.freeze_product(view="full", source="tenv", views_root=self.root / "views",
                              destination=self.freeze, repo_root_path=REPO,
                              sealed=self.sealed, axis_dir=self.fixture["axis"],
                              accepted_dir=self.accepted)
        self.run_dir = self.root / "runs/tenv/full/fold-0/seed-42"

    def worker(self, phase, *extra):
        return run(WORKER, self.freeze, SOURCE, self.run_dir, phase,
                   "--sealed", self.sealed, "--repo-root", REPO, *extra)

    def terminal(self):
        return json.loads((self.run_dir / "terminal.json").read_text())

    def test_bounded_smoke_fit_trains_replays_and_scores(self):
        trained = self.worker("train", "--rows", 4, "--epochs", 1, "--device", "cpu",
                              "--fold", 0, "--seed", 42)
        self.assertEqual(trained.returncode, 0, trained.stderr)
        terminal = self.terminal()
        self.assertEqual(terminal["status"], "validation-trained-awaiting-replay")
        self.assertEqual(terminal["best_epoch"], 1)
        self.assertEqual(terminal["epochs"], 1)
        self.assertEqual(terminal["device"], "cpu")
        self.assertEqual(terminal["seed"], 42)
        self.assertGreater(terminal["peak_rss_kib"], 0)
        self.assertEqual(sorted(path.name for path in self.run_dir.iterdir()),
                         ["checkpoint.pt", "config.json", "learning_curve.csv",
                          "terminal.json", "validation_metrics.json",
                          "validation_prediction.npy"])

        # No test row may be read before the validation choice is sealed by a replay.
        premature = self.worker("test")
        self.assertNotEqual(premature.returncode, 0)
        self.assertIn("requires a replayed validation choice", premature.stderr)
        self.assertFalse((self.run_dir / "test_predictions.npy").exists())

        replayed = self.worker("replay")
        self.assertEqual(replayed.returncode, 0, replayed.stderr)
        self.assertEqual(self.terminal()["status"], "validation-complete")
        self.assertTrue((self.run_dir / "validation_manifest.json").is_file())

        scored = self.worker("test")
        self.assertEqual(scored.returncode, 0, scored.stderr)
        self.assertEqual(self.terminal()["status"], "test-awaiting-replay")
        self.assertFalse((self.run_dir / "manifest.json").exists())
        metrics = json.loads((self.run_dir / "test_metrics.json").read_text())
        self.assertEqual(metrics["n_samples"], 2)

        final = self.worker("replay")
        self.assertEqual(final.returncode, 0, final.stderr)
        terminal = self.terminal()
        self.assertEqual(terminal["status"], "complete")
        manifest = json.loads((self.run_dir / "manifest.json").read_text())
        self.assertEqual(manifest["schema"], "tsg-view-run-v1")
        self.assertEqual(manifest["origin"], {
            "source": "tenv", "view": "full", "fold": 0, "seed": 42,
            "scenario": "cold-drug-dose", "method": "no-dose"})
        self.assertIn("test_predictions.npy", manifest["outputs"])
        self.assertIn("validation_prediction.npy", manifest["outputs"])
        # A complete run re-replays bit-exactly in yet another fresh process.
        again = self.worker("replay")
        self.assertEqual(again.returncode, 0, again.stderr)
        self.assertEqual(self.terminal()["status"], "complete")

    def test_smoke_driver_bounds_the_fit_and_records_resources(self):
        receipt_path = self.root / "smoke-receipt.json"
        receipt = smoke.smoke(view="full", source="tenv", freeze=self.freeze,
                              output=self.root / "smoke", rows=4, epochs=1, device="cpu",
                              repo_root_path=REPO, sealed=self.sealed,
                              source_dir=SOURCE, record=receipt_path)
        self.assertTrue(receipt_path.is_file())
        self.assertEqual(receipt["status"], "complete")
        self.assertEqual(receipt["epochs_run"], 1)
        self.assertEqual([phase["phase"] for phase in receipt["phases"]],
                         ["train", "replay", "test", "replay"])
        self.assertGreater(receipt["children_peak_rss_kib"], 0)
        self.assertTrue(all(phase["returncode"] == 0 for phase in receipt["phases"]))
        self.assertFalse((self.root / "smoke/smoke.json").exists())

    def test_smoke_bounds_are_refused(self):
        with self.assertRaisesRegex(ValueError, "smoke rows must be in 2..64"):
            smoke.smoke(view="full", source="tenv", freeze=self.freeze,
                        output=self.root / "smoke", rows=65, epochs=1, device="cpu",
                        repo_root_path=REPO, sealed=self.sealed, source_dir=SOURCE)
        with self.assertRaisesRegex(ValueError, "smoke epochs must be in 1..3"):
            smoke.smoke(view="full", source="tenv", freeze=self.freeze,
                        output=self.root / "smoke", rows=4, epochs=4, device="cpu",
                        repo_root_path=REPO, sealed=self.sealed, source_dir=SOURCE)


if __name__ == "__main__":
    unittest.main()
