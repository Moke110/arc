"""Plan/enumeration, seal and validation-receipt seams of the formal matrix."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXECUTION = HERE.parent
for path in (str(EXECUTION), str(HERE)):
    if path not in sys.path:
        sys.path.insert(0, path)

import matrix  # noqa: E402

REPO = Path(__file__).resolve().parents[4]
SOURCE = next(candidate for candidate in (
    REPO / "pert_model/TranSiGen/source",
    Path("/home/chang/arc/pert_model/TranSiGen/source"))
    if (candidate / "src/model.py").is_file())

# The intended schedule, spelled out independently of matrix.plan(): every tenv fit
# first (the report's tenv pass must publish first), then every Tahoe fit. Within a
# source the driver recovers the schedule by (fold) group and runs all of a group's
# validations before sealing that fold's receipt, so the order here is the file order.
_FOLD_VIEWS = (
    ("full", 42), ("full", 43),
    ("merge", 42), ("merge", 43),
    ("solo", 42), ("solo", 43),
    ("full-pooled-label", 42), ("full-pooled-label", 43),
    ("solo-deepest", 42),
)
_TAHOE_VIEWS = (
    ("full", 42), ("full", 43),
    ("merge", 42), ("merge", 43),
    ("solo", 42), ("solo", 43),
)
EXPECTED_ORDER = (
    [("tenv", view, fold, seed) for fold in range(3) for view, seed in _FOLD_VIEWS]
    + [("tahoe_bulk", view, fold, seed) for fold in range(3) for view, seed in _TAHOE_VIEWS]
)


def _write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def _stub_inputs(root):
    """Contract-shaped freeze and truth manifests; only their digests are plan inputs."""
    freezes, truth = Path(root) / "freezes", Path(root) / "truth"
    for source in matrix.SOURCES:
        views = matrix.views_for(source)
        for view in views:
            directory = freezes / view / source
            _write(directory / "manifest.json",
                   {"schema": "tsg-view-freeze-v1", "outputs": {"splits.csv": "0" * 64}})
            _write(directory / "config.json",
                   {"view": view, "source": source, "n_treated": 12})
    for source in matrix.SOURCES:
        _write(truth / source / "manifest.json",
               {"schema": "tsg-view-truth-v1", "source": source})
        for fold in range(3):
            for kind in ("primary", "condition_mean"):
                _write(truth / source / kind / f"fold-{fold}" / "manifest.json",
                       {"schema": "tsg-view-truth-v1", "source": source, "fold": fold,
                        "kind": kind, "rows": 4})
    return freezes, truth


def _plan(root, **kwargs):
    freezes, truth = _stub_inputs(root)
    arguments = dict(freezes=freezes, truth=truth, source_dir=SOURCE, device="cpu",
                     repo_root_path=REPO)
    arguments.update(kwargs)
    return matrix.frozen_plan(**arguments)


def _candidate(root):
    candidate = Path(root) / "candidate"
    candidate.mkdir(parents=True, exist_ok=True)
    return candidate


class PlanEnumerationTest(unittest.TestCase):
    def test_plan_is_exactly_the_45_intended_origins_in_order(self):
        origins = matrix.plan()
        actual = [(origin["source"], origin["view"], origin["fold"], origin["seed"])
                  for origin in origins]
        self.assertEqual(actual, EXPECTED_ORDER)
        self.assertEqual(len(actual), 45)
        self.assertEqual(len(set(actual)), 45)  # run identity is unique
        self.assertEqual(sum(source == "tenv" for source, _, _, _ in actual), 27)
        self.assertEqual(sum(source == "tahoe_bulk" for source, _, _, _ in actual), 18)
        first_tahoe = next(index for index, row in enumerate(actual) if row[0] == "tahoe_bulk")
        self.assertTrue(all(row[0] == "tenv" for row in actual[:first_tahoe]))
        self.assertTrue(all(row[0] == "tahoe_bulk" for row in actual[first_tahoe:]))
        # solo-deepest is the sensitivity arm: tenv only, seed 42 only.
        deepest = [row for row in actual if row[1] == "solo-deepest"]
        self.assertEqual(deepest, [("tenv", "solo-deepest", fold, 42) for fold in range(3)])
        self.assertTrue(all(origin["scenario"] == "cold-drug-dose"
                            and origin["method"] == "no-dose" for origin in origins))

    def test_origin_paths_follow_the_identity(self):
        origin = {"source": "tenv", "view": "merge", "fold": 2, "seed": 43}
        self.assertEqual(matrix.run_relative(origin).as_posix(),
                         "runs/tenv/merge/fold-2/seed-43")
        self.assertEqual(matrix.run_path(Path("/c"), origin), Path("/c/runs/tenv/merge/fold-2/seed-43"))
        self.assertEqual(matrix.selection_path(Path("/c"), "tenv", 1),
                         Path("/c/selection/tenv/fold-1.json"))
        self.assertEqual(matrix.record_relative(origin).as_posix(),
                         "records/tenv/merge/fold-2/seed-43.json")
        self.assertEqual(matrix.freeze_path(Path("/f"), "solo-deepest", "tenv"),
                         Path("/f/solo-deepest/tenv"))

    def test_filters_select_schedule_subsets(self):
        origins = matrix.plan()
        self.assertEqual(len(matrix.filter_origins(origins, only_source="tenv")), 27)
        self.assertEqual(len(matrix.filter_origins(origins, only_source="tahoe_bulk")), 18)
        self.assertEqual(len(matrix.filter_origins(origins, only_source="tenv", only_fold=0)), 9)
        picked = matrix.filter_origins(origins, only_view="full", only_seed=43)
        self.assertEqual([(o["source"], o["fold"]) for o in picked],
                         [("tenv", 0), ("tenv", 1), ("tenv", 2),
                          ("tahoe_bulk", 0), ("tahoe_bulk", 1), ("tahoe_bulk", 2)])
        with self.assertRaisesRegex(ValueError, "no origin matches"):
            matrix.filter_origins(origins, only_source="tenv", only_view="nope")

    def test_resume_state_machine(self):
        self.assertEqual(matrix.next_action(None, "validation"), "train")
        self.assertEqual(matrix.next_action("validation-trained-awaiting-replay", "validation"),
                         "replay")
        self.assertIsNone(matrix.next_action("validation-complete", "validation"))
        self.assertEqual(matrix.next_action("validation-complete", "test"), "evaluate")
        self.assertEqual(matrix.next_action("test-awaiting-replay", "test"), "replay")
        self.assertIsNone(matrix.next_action("complete", "test"))
        for phase in ("validation", "test"):
            for broken in ("failed", "interrupted", "bogus"):
                with self.assertRaises(ValueError):
                    matrix.next_action(broken, phase)


class PlanSealTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        self.candidate = _candidate(self.root)
        self.plan = _plan(self.root)

    def test_seal_writes_plan_and_digest_manifest_before_any_fit(self):
        path = matrix.seal_plan(self.candidate, self.plan)
        self.assertEqual(path, self.candidate / "formal_plan.json")
        self.assertTrue((self.candidate / "formal_plan_manifest.json").is_file())
        self.assertFalse((self.candidate / "runs").exists())
        sealed = matrix.load_plan(self.candidate)
        self.assertEqual(sealed, self.plan)
        manifest = json.loads((self.candidate / "formal_plan_manifest.json").read_text())
        self.assertEqual(manifest["schema"], matrix.PLAN_MANIFEST_SCHEMA)
        self.assertEqual(manifest["plan_sha256"], matrix.sha256(path))
        self.assertEqual(manifest["origins"], 45)
        # Idempotent: sealing the identical plan again succeeds and does not rewrite.
        before = path.read_bytes()
        matrix.seal_plan(self.candidate, self.plan)
        self.assertEqual(path.read_bytes(), before)

    def test_plan_manifest_mismatch_fails_closed(self):
        matrix.seal_plan(self.candidate, self.plan)
        manifest_path = self.candidate / "formal_plan_manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["plan_sha256"] = "0" * 64
        manifest_path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError, "plan"):
            matrix.load_plan(self.candidate)
        # A tampered plan body with an intact manifest is equally refused.
        self.candidate = _candidate(Path(self._tmp.name) / "second")
        matrix.seal_plan(self.candidate, self.plan)
        path = self.candidate / "formal_plan.json"
        tampered = json.loads(path.read_text())
        tampered["device"] = "cpu "
        path.write_text(json.dumps(tampered))
        with self.assertRaisesRegex(ValueError, "plan"):
            matrix.load_plan(self.candidate)
        # A half-written seal (one file only) is not a plan.
        (self.candidate / "formal_plan_manifest.json").unlink()
        with self.assertRaisesRegex(ValueError, "plan"):
            matrix.load_plan(self.candidate)

    def test_plan_cannot_be_sealed_refused_or_replaced_after_fits(self):
        runs = self.candidate / "runs/tenv/full/fold-0/seed-42"
        runs.mkdir(parents=True)
        with self.assertRaisesRegex(ValueError, "missing sealed plan"):
            matrix.seal_plan(self.candidate, self.plan)
        import shutil
        shutil.rmtree(self.candidate / "runs")
        matrix.seal_plan(self.candidate, self.plan)
        other = json.loads(json.dumps(self.plan))
        other["device"] = "cpu "
        with self.assertRaisesRegex(ValueError, "plan"):
            matrix.seal_plan(self.candidate, other)

    def test_check_plan_detects_runtime_or_input_drift(self):
        matrix.seal_plan(self.candidate, self.plan)
        matrix.check_plan(self.candidate, self.plan)
        other = json.loads(json.dumps(self.plan))
        other["device"] = "cpu "
        with self.assertRaisesRegex(ValueError, "plan"):
            matrix.check_plan(self.candidate, other)


class SelectionReceiptTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        self.candidate = _candidate(self.root)
        self.plan = _plan(self.root)
        matrix.seal_plan(self.candidate, self.plan)
        self.plan_sha = matrix.plan_sha256(self.candidate)
        for origin in matrix.filter_origins(self.plan["origins"],
                                            only_source="tenv", only_fold=0):
            self._fake_validated(origin)

    def _fake_validated(self, origin, *, mse=0.5):
        run = matrix.run_path(self.candidate, origin)
        run.mkdir(parents=True, exist_ok=True)
        freeze = self.plan["freezes"][f"{origin['view']}/{origin['source']}"]
        (run / "config.json").write_text(json.dumps({
            "schema": "tsg-view-run-v1", "origin": origin,
            "plan": {"sha256": self.plan_sha},
            "freeze": {"manifest_sha256": freeze["manifest_sha256"]},
            "identity": {"mode": "formal", "device": "cpu", "smoke_rows": None}}))
        (run / "validation_prediction.npy").write_bytes(b"")
        (run / "validation_metrics.json").write_text(json.dumps({"pcc": 0.3}))
        matrix.write_json(run / "validation_manifest.json", {
            "schema": "tsg-view-validation-v1", "input_manifest_sha256":
            freeze["manifest_sha256"], "outputs": {}})
        matrix.write_json(run / "terminal.json", {
            "status": "validation-complete", "seed": origin["seed"],
            "fold": origin["fold"], "device": "cpu", "epochs": 12, "best_epoch": 7,
            "validation_mse": mse, "parameter_count": 21583000,
            "training_process_token": "t" * 64,
            "validation_replay_token": "v" * 64})
        return run

    def test_receipt_enumerates_the_fold_origins_before_any_test(self):
        path = matrix.lock_selection(self.candidate, "tenv", 0)
        self.assertEqual(path, matrix.selection_path(self.candidate, "tenv", 0))
        receipt = json.loads(path.read_text())
        self.assertEqual(receipt["schema"], matrix.SELECTION_SCHEMA)
        self.assertEqual(receipt["source"], "tenv")
        self.assertEqual(receipt["fold"], 0)
        self.assertEqual(receipt["plan_sha256"], self.plan_sha)
        self.assertEqual([(entry["view"], entry["seed"]) for entry in receipt["origins"]],
                         [("full", 42), ("full", 43), ("merge", 42), ("merge", 43),
                          ("solo", 42), ("solo", 43), ("full-pooled-label", 42),
                          ("full-pooled-label", 43), ("solo-deepest", 42)])
        first = receipt["origins"][0]
        self.assertEqual(first["best_epoch"], 7)
        self.assertEqual(first["epochs"], 12)
        self.assertEqual(first["validation_mse"], 0.5)
        self.assertEqual(first["parameter_count"], 21583000)
        self.assertIn("validation_manifest_sha256", first)
        self.assertIn("truth", receipt)
        self.assertNotIn("selected_method", receipt)

    def test_receipt_is_idempotent_and_changes_fail_closed(self):
        matrix.lock_selection(self.candidate, "tenv", 0)
        path = matrix.selection_path(self.candidate, "tenv", 0)
        before = path.read_bytes()
        # Recompute after a test evaluation: the validation seal is byte-identical.
        first = matrix.filter_origins(self.plan["origins"], only_source="tenv", only_fold=0)[0]
        terminal_path = matrix.run_path(self.candidate, first) / "terminal.json"
        terminal = json.loads(terminal_path.read_text())
        terminal.update(status="complete", test_evaluation_process_token="e" * 64)
        terminal_path.write_text(json.dumps(terminal))
        matrix.lock_selection(self.candidate, "tenv", 0)
        self.assertEqual(matrix.selection_evidence(self.candidate, "tenv", 0),
                         json.loads(path.read_text()))
        self.assertEqual(path.read_bytes(), before)
        # A changed validation MSE is a changed choice.
        second = matrix.filter_origins(self.plan["origins"], only_source="tenv", only_fold=0)[1]
        run = matrix.run_path(self.candidate, second)
        terminal = json.loads((run / "terminal.json").read_text())
        terminal["validation_mse"] = 0.25
        (run / "terminal.json").write_text(json.dumps(terminal))
        with self.assertRaisesRegex(ValueError, "selection"):
            matrix.lock_selection(self.candidate, "tenv", 0)

    def test_missing_validation_seal_fails_closed(self):
        origin = matrix.filter_origins(self.plan["origins"], only_source="tenv",
                                       only_fold=1)[0]
        self._fake_validated(origin)
        run = matrix.run_path(self.candidate, origin)
        (run / "validation_manifest.json").unlink()
        with self.assertRaisesRegex(ValueError, "validation"):
            matrix.selection_evidence(self.candidate, "tenv", 1)


class FullArmToleranceTest(unittest.TestCase):
    def test_anchors_are_the_sealed_no_dose_cold_drug_dose_values(self):
        self.assertEqual(matrix.FULL_ANCHORS["tenv"]["pcc"], 0.2654588998492864)
        self.assertEqual(matrix.FULL_ANCHORS["tenv"]["rmse"], 0.421337591468443)
        self.assertEqual(matrix.FULL_ANCHORS["tahoe_bulk"]["pcc"], 0.4226328493964055)
        self.assertEqual(matrix.FULL_ANCHORS["tahoe_bulk"]["rmse"], 0.15479625072122077)
        # The declared tolerance is explicit, not a hidden epsilon.
        self.assertGreater(matrix.FULL_TOLERANCE["origin"]["pcc_abs"], 0)
        self.assertGreater(matrix.FULL_TOLERANCE["origin"]["rmse_abs"], 0)

    def test_compare_full_arm_reports_pass_or_drift_and_never_hides_drift(self):
        anchor = matrix.FULL_ANCHORS["tenv"]
        passing = matrix.compare_full_arm("tenv", {"pcc": anchor["pcc"] + 0.001,
                                                   "rmse": anchor["rmse"] - 0.001})
        self.assertEqual(passing["status"], "PASS")
        drift = matrix.compare_full_arm("tenv", {"pcc": anchor["pcc"] + 0.2,
                                                 "rmse": anchor["rmse"]})
        self.assertEqual(drift["status"], "DRIFT")
        self.assertAlmostEqual(drift["drift"]["pcc"], 0.2)
        self.assertAlmostEqual(drift["measured"]["pcc"], anchor["pcc"] + 0.2)
        with self.assertRaisesRegex(ValueError, "anchor"):
            matrix.compare_full_arm("unknown", {"pcc": 0.1, "rmse": 0.1})


if __name__ == "__main__":
    unittest.main()
