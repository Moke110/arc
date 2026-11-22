"""Formal-driver seam: plan identity, telemetry, train/test isolation, replay, receipts."""
import json
import math
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
EXECUTION = HERE.parent
for path in (str(EXECUTION), str(HERE)):
    if path not in sys.path:
        sys.path.insert(0, path)

import build  # noqa: E402
import tsg_view_fixtures as fixtures  # noqa: E402
import formal  # noqa: E402
import freeze  # noqa: E402
import learner  # noqa: E402
import matrix  # noqa: E402
import views  # noqa: E402

REPO = Path(__file__).resolve().parents[4]
SOURCE = next(candidate for candidate in (
    REPO / "pert_model/TranSiGen/source",
    Path("/home/chang/arc/pert_model/TranSiGen/source"))
    if (candidate / "src/model.py").is_file())
SIZE = 5000
FOLD = 0


def _write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def _stub_plan_inputs(root):
    freezes, truth = Path(root) / "freezes", Path(root) / "truth"
    for source in matrix.SOURCES:
        for view in matrix.views_for(source):
            _write_json(freezes / view / source / "manifest.json",
                        {"schema": matrix.FREEZE_SCHEMA, "outputs": {}})
            _write_json(freezes / view / source / "config.json",
                        {"view": view, "source": source, "n_treated": 12})
        _write_json(truth / source / "manifest.json",
                    {"schema": matrix.TRUTH_SCHEMA, "source": source})
        for fold in range(3):
            for kind in ("primary", "condition_mean"):
                _write_json(truth / source / kind / f"fold-{fold}" / "manifest.json",
                            {"schema": matrix.TRUTH_SCHEMA, "source": source,
                             "fold": fold, "kind": kind, "rows": 2})
    return freezes, truth


def write_truth(root, freeze_dir, source="tenv", fold=FOLD):
    """Write a contract-shaped frozen truth product from a ``full`` freeze's test rows."""
    root, freeze_dir = Path(root), Path(freeze_dir)
    splits = pd.read_csv(freeze_dir / "splits.csv")
    population = pd.read_csv(freeze_dir / "population.csv")
    rows = splits.loc[(splits.scenario == matrix.SCENARIO) & (splits.fold == fold)
                      & (splits.partition == "test"), "row"].to_numpy(dtype=int)
    expression = np.load(freeze_dir / "expression.npy", mmap_mode="r")
    control = np.load(freeze_dir / "control.npy", mmap_mode="r")
    delta = (np.asarray(expression[rows], dtype=np.float32)
             - np.asarray(control[rows], dtype=np.float32)).astype(np.float32)
    primary = Path(root) / source / "primary" / f"fold-{fold}"
    primary.mkdir(parents=True, exist_ok=True)
    body = population.iloc[rows].reset_index(drop=True)
    frame = pd.DataFrame({
        "profile_id": body["profile_id"].astype(str),
        "cell_line": body["cell_line"].astype(str),
        "drug_parent": body["drug_parent_inchikey"].astype(str),
        "dose_um": body["dose_um"].astype(float),
        "partition": "test", "fold": fold, "scenario": matrix.SCENARIO})
    frame.to_csv(primary / "population.csv", index=False, lineterminator="\n")
    np.save(primary / "delta.npy", delta)
    np.save(primary / "control.npy", np.asarray(control[rows], dtype=np.float32))
    shutil.copyfile(freeze_dir / "parents.csv", primary / "parents.csv")
    shutil.copyfile(freeze_dir / "morgan.npy", primary / "morgan.npy")
    files = ("population.csv", "delta.npy", "control.npy", "parents.csv", "morgan.npy")
    _write_json(primary / "manifest.json", {
        "schema": matrix.TRUTH_SCHEMA, "source": source, "fold": fold,
        "kind": "primary", "rows": len(frame),
        "outputs": {name: views.sha256(primary / name) for name in files}})

    # Condition-mean truth: pooled Δ of every replicated test condition (merge recipe
    # shape; the fixture's exact pooling rule is not under test here).
    keys = ["cell_line", "drug_parent_inchikey", "dose_um"]
    local = pd.DataFrame({"row": np.arange(len(rows))}).join(
        body[keys].reset_index(drop=True))
    members, pooled, identities = [], [], []
    for position, group in local.groupby(keys, sort=True):
        if len(group) < 2:
            continue
        members.append(";".join(body["profile_id"].iloc[group["row"]].astype(str)))
        pooled.append(delta[group["row"]].astype(np.float64).mean(axis=0).astype(np.float32))
        identities.append({"profile_id": f"cond-{len(identities)}",
                           "cell_line": position[0], "drug_parent": position[1],
                           "dose_um": float(position[2]), "count": len(group)})
    condition = Path(root) / source / "condition_mean" / f"fold-{fold}"
    condition.mkdir(parents=True, exist_ok=True)
    condition_frame = pd.DataFrame(identities)
    condition_frame["members"] = members
    condition_frame["partition"] = "test"
    condition_frame["fold"] = fold
    condition_frame["scenario"] = matrix.SCENARIO
    condition_frame.to_csv(condition / "population.csv", index=False, lineterminator="\n")
    np.save(condition / "delta.npy", np.stack(pooled).astype(np.float32)
            if pooled else np.empty((0, SIZE), np.float32))
    files = ("population.csv", "delta.npy")
    _write_json(condition / "manifest.json", {
        "schema": matrix.TRUTH_SCHEMA, "source": source, "fold": fold,
        "kind": "condition_mean", "rows": len(condition_frame),
        "outputs": {name: views.sha256(condition / name) for name in files}})
    _write_json(Path(root) / source / "manifest.json", {
        "schema": matrix.TRUTH_SCHEMA, "source": source, "folds": [FOLD]})
    return Path(root) / source


class BlockedStream:
    """A stdout whose every write fails exactly like a full non-blocking pipe."""

    def write(self, *args, **kwargs):
        raise BlockingIOError(11, "write could not complete without blocking")

    def flush(self):
        raise BlockingIOError(11, "write could not complete without blocking")

    def isatty(self):
        return False


class DriverCodeDriftTest(unittest.TestCase):
    """A driver-only fix resumes a sealed plan; a scientific change fails closed."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        freezes, truth = _stub_plan_inputs(self.root)
        self.plan = matrix.frozen_plan(freezes=freezes, truth=truth, source_dir=SOURCE,
                                       device="cpu", mode="formal", repo_root_path=REPO)
        self.candidate = self.root / "candidate"
        matrix.seal_plan(self.candidate, self.plan)
        self.plan_path = self.candidate / "formal_plan.json"
        self.before = self.plan_path.read_bytes()

    def test_driver_only_drift_is_reported_and_the_seal_is_not_rewritten(self):
        drifted = json.loads(json.dumps(self.plan))
        drifted["code"]["formal.py"] = "f" * 64
        drift = formal.driver_code_drift(self.plan, drifted)
        self.assertEqual([entry["file"] for entry in drift], ["formal.py"])
        self.assertEqual(formal.driver_code_drift(self.plan, self.plan), [])
        self.assertEqual(matrix.load_plan(self.candidate), self.plan)
        self.assertEqual(self.plan_path.read_bytes(), self.before)

    def test_scientific_code_or_input_change_still_fails_closed(self):
        scientific = json.loads(json.dumps(self.plan))
        scientific["code"]["learner.py"] = "f" * 64
        self.assertIsNone(formal.driver_code_drift(self.plan, scientific))
        truth = json.loads(json.dumps(self.plan))
        truth["truth"]["tenv/fold-0"]["primary"]["manifest_sha256"] = "f" * 64
        self.assertIsNone(formal.driver_code_drift(self.plan, truth))


class FinalReportTest(unittest.TestCase):
    """A blocked stdout must never decide the driver's exit status."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        self.candidate = self.root / "candidate"
        self.candidate.mkdir(parents=True)
        self.report = {"schema": formal.RUN_REPORT_SCHEMA, "mode": "formal",
                       "plan_sha256": "0" * 64, "selected": 45, "complete": 45,
                       "runs": [{"origin": {"source": "tenv", "view": "full",
                                            "fold": 0, "seed": 42}}],
                       "driver_code_drift": []}

    def _run_arguments(self):
        return ["run", "--candidate", str(self.candidate), "--no-scope", "--device", "cpu"]

    def test_blocked_stdout_does_not_fail_a_completed_run(self):
        with mock.patch.object(formal, "preflight",
                               return_value={"available_ram_kib": 1, "gpu": None}), \
                mock.patch.object(formal, "run_matrix", return_value=self.report), \
                mock.patch.object(sys, "stdout", BlockedStream()):
            status = formal.main(self._run_arguments())
        self.assertEqual(status, 0)
        path = self.candidate / "formal_run_report.json"
        self.assertTrue(path.is_file())
        self.assertEqual(json.loads(path.read_text()), self.report)
        self.assertTrue(path.read_text().endswith("\n"))
        self.assertIn('"complete": 45', path.read_text())

    def test_blocked_stdout_does_not_fail_verify(self):
        report = {"schema": formal.VERIFY_SCHEMA, "status": "PASS", "origins": 45}
        with mock.patch.object(formal, "verify", return_value=report), \
                mock.patch.object(sys, "stdout", BlockedStream()):
            status = formal.main(["verify", "--candidate", str(self.candidate)])
        self.assertEqual(status, 0)

    def test_blocked_stdout_does_not_fail_plan(self):
        freezes, truth = _stub_plan_inputs(self.root)
        with mock.patch.object(sys, "stdout", BlockedStream()):
            status = formal.main([
                "plan", "--candidate", str(self.candidate / "planned"),
                "--freezes", str(freezes), "--truth", str(truth),
                "--native-source", str(SOURCE), "--repo-root", str(REPO),
                "--device", "cpu"])
        self.assertEqual(status, 0)
        self.assertEqual(len(matrix.load_plan(self.candidate / "planned")["origins"]), 45)

    def test_shutdown_flush_cannot_change_the_exit_status(self):
        with mock.patch.object(sys, "stdout", BlockedStream()):
            formal.quiet_stdout()
            self.assertIsInstance(sys.stdout, formal._NullStdout)
            self.assertEqual(sys.stdout.write("ignored"), 0)
            sys.stdout.flush()

    def test_a_genuine_phase_failure_still_propagates(self):
        failure = RuntimeError("a phase really failed")
        with mock.patch.object(formal, "preflight",
                               return_value={"available_ram_kib": 1, "gpu": None}), \
                mock.patch.object(formal, "run_matrix", side_effect=failure):
            with self.assertRaisesRegex(RuntimeError, "really failed"):
                formal.main(self._run_arguments())
        self.assertFalse((self.candidate / "formal_run_report.json").exists())


class TelemetryTest(unittest.TestCase):
    def test_meminfo_and_gpu_sample_parsers(self):
        meminfo = ("MemTotal:       24000000 kB\nMemFree: 1 kB\n"
                   "MemAvailable:   16000000 kB\nSwapFree:    8000000 kB\n")
        self.assertEqual(formal.parse_meminfo(meminfo),
                         {"mem_total_kib": 24000000, "mem_available_kib": 16000000,
                          "swap_free_kib": 8000000})
        with self.assertRaisesRegex(ValueError, "misses"):
            formal.parse_meminfo("MemTotal: 1 kB\n")
        self.assertEqual(formal.parse_gpu_line("2048, 14000, 37"),
                         {"gpu_used_mib": 2048, "gpu_free_mib": 14000,
                          "gpu_util_pct": 37})
        with self.assertRaisesRegex(ValueError, "nvidia-smi"):
            formal.parse_gpu_line("2048, nope, 37")

    def test_telemetry_summary_records_peaks_and_missing_gpu(self):
        samples = [{"mem_available_kib": 1000, "mem_total_kib": 2000, "swap_free_kib": 1,
                    "gpu_used_mib": 100, "gpu_free_mib": 900, "gpu_util_pct": 10},
                   {"mem_available_kib": 800, "mem_total_kib": 2000, "swap_free_kib": 1,
                    "gpu_used_mib": 400, "gpu_free_mib": 600, "gpu_util_pct": 30},
                   {"mem_available_kib": 900, "mem_total_kib": 2000, "swap_free_kib": 1,
                    "gpu_used_mib": None, "gpu_free_mib": None, "gpu_util_pct": None}]
        summary = formal.summarize_telemetry(samples)
        self.assertEqual(summary["samples"], 3)
        self.assertEqual(summary["min_mem_available_kib"], 800)
        self.assertEqual(summary["peak_gpu_used_mib"], 400)
        self.assertEqual(summary["mean_gpu_util_pct"], 20.0)
        empty = formal.summarize_telemetry([])
        self.assertEqual(empty["samples"], 0)
        self.assertIsNone(empty["peak_gpu_used_mib"])

    def test_sample_resources_reads_a_meminfo_file_and_injected_gpu(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "meminfo"
            path.write_text("MemTotal: 10 kB\nMemAvailable: 8 kB\nSwapFree: 2 kB\n")
            completed = SimpleNamespace(
                stdout="1234, 15000, 55\n", returncode=0)
            sample = formal.sample_resources(
                meminfo_path=path, nvidia="/bin/true",
                runner=lambda *args, **kwargs: completed)
            self.assertEqual(sample["mem_available_kib"], 8)
            self.assertEqual(sample["gpu_used_mib"], 1234)
            self.assertEqual(sample["gpu_util_pct"], 55)


class PreflightTest(unittest.TestCase):
    def test_memory_floor_is_a_hard_refusal(self):
        with tempfile.TemporaryDirectory() as tmp:
            low = Path(tmp) / "low"
            low.write_text("MemTotal: 1 kB\nMemAvailable: 1000 kB\nSwapFree: 0 kB\n")
            with self.assertRaisesRegex(RuntimeError, "below"):
                formal.memory_preflight(12 * 1024 * 1024, meminfo_path=low)
            high = Path(tmp) / "high"
            high.write_text("MemTotal: 1 kB\nMemAvailable: 20000000 kB\nSwapFree: 0 kB\n")
            self.assertEqual(formal.memory_preflight(12 * 1024 * 1024, meminfo_path=high),
                             20000000)
            value = formal.preflight("cpu", 12 * 1024 * 1024, meminfo_path=high)
            self.assertEqual(value["available_ram_kib"], 20000000)
            self.assertIsNone(value["gpu"])

    def test_scope_command_declares_the_resource_caps(self):
        command = formal.scope_command(["python", "formal.py", "run"], "unit")
        text = " ".join(command)
        self.assertIn("systemd-run --user --scope", text)
        self.assertIn("MemoryMax=10G", text)
        self.assertIn("MemorySwapMax=8G", text)


class PlanBindingTest(unittest.TestCase):
    def test_formal_plan_binds_every_freeze_truth_and_code_digest(self):
        with tempfile.TemporaryDirectory() as tmp:
            freezes, truth = _stub_plan_inputs(tmp)
            plan = matrix.frozen_plan(freezes=freezes, truth=truth, source_dir=SOURCE,
                                      device="cpu", mode="formal", repo_root_path=REPO)
            self.assertEqual(plan["mode"], "formal")
            self.assertEqual(len(plan["origins"]), 45)
            self.assertEqual(sorted(plan["freezes"]), sorted(
                f"{view}/{source}" for source in matrix.SOURCES
                for view in matrix.views_for(source)))
            self.assertEqual(len(plan["freezes"]), 8)
            self.assertEqual(len(plan["truth"]), 6)
            self.assertEqual(sorted(plan["truth_source"]), ["tahoe_bulk", "tenv"])
            self.assertEqual(set(plan["code"]), set(matrix.CODE_FILES))
            self.assertEqual(plan["caps"], matrix.CAPS)
            self.assertEqual(plan["full_arm"]["anchors"], matrix.FULL_ANCHORS)
            self.assertEqual(plan["max_epochs"], 1000)
            self.assertTrue(plan["native_source"]["model_sha256"])
            candidate = Path(tmp) / "candidate"
            matrix.seal_plan(candidate, plan)
            self.assertEqual(matrix.load_plan(candidate), plan)
            # A formal plan refuses filters: they are a scheduling device only.
            with self.assertRaisesRegex(ValueError, "filters"):
                matrix.frozen_plan(freezes=freezes, truth=truth, source_dir=SOURCE,
                                   device="cpu", mode="formal", only_source=("tenv",),
                                   repo_root_path=REPO)

    def test_probe_plan_is_the_filtered_bounded_schedule(self):
        with tempfile.TemporaryDirectory() as tmp:
            freezes, truth = _stub_plan_inputs(tmp)
            plan = matrix.frozen_plan(
                freezes=freezes, truth=truth, source_dir=SOURCE, device="cpu",
                mode="probe", probe_rows=8, probe_epochs=2,
                only_source=("tenv",), only_view=("full",), only_fold=(0,),
                only_seed=(42, 43), repo_root_path=REPO)
            self.assertEqual(plan["mode"], "probe")
            self.assertEqual(plan["probe"], {"rows": 8, "epochs": 2})
            self.assertEqual([(o["source"], o["view"], o["fold"], o["seed"])
                              for o in plan["origins"]],
                             [("tenv", "full", 0, 42), ("tenv", "full", 0, 43)])
            self.assertEqual(list(plan["freezes"]), ["full/tenv"])
            self.assertEqual(list(plan["truth"]), ["tenv/fold-0"])
            self.assertEqual(plan["max_epochs"], 2)

    def test_filters_may_not_split_a_fold_group(self):
        with tempfile.TemporaryDirectory() as tmp:
            freezes, truth = _stub_plan_inputs(tmp)
            plan = matrix.frozen_plan(freezes=freezes, truth=truth, source_dir=SOURCE,
                                      device="cpu", mode="formal", repo_root_path=REPO)
            selected = matrix.filter_origins(plan["origins"], only_seed=42)
            with self.assertRaisesRegex(ValueError, "split"):
                formal.check_group_completeness(plan, selected)
            whole = matrix.filter_origins(plan["origins"], only_source="tenv",
                                          only_fold=0)
            formal.check_group_completeness(plan, whole)  # complete group, no refusal


class ScoreDeltaTest(unittest.TestCase):
    def test_score_delta_reduction_identities(self):
        vector_metrics = learner.api(REPO)["vector_metrics"]
        obs = pd.DataFrame({"profile_id": ["a", "b", "c"], "drug_parent": ["P"] * 3,
                            "dose_um": [1.0] * 3, "cell_line": ["L"] * 3})
        control = np.zeros((3, 6), dtype=np.float32)
        truth = np.zeros((3, 6), dtype=np.float32)
        truth[1] = np.array([1, 2, 3, 4, 5, 6], dtype=np.float32)
        truth[2] = np.array([1, 2, 3, 4, 5, 6], dtype=np.float32)
        prediction = np.zeros((3, 6), dtype=np.float32)
        prediction[1] = np.array([2, 4, 6, 8, 10, 12], dtype=np.float32)
        samples = learner.score_delta(obs, truth, prediction, control, vector_metrics)
        constant_truth = samples.iloc[0]
        self.assertTrue(constant_truth["constant_truth"])
        self.assertTrue(np.isnan(constant_truth["pcc"]))
        perfect = samples.iloc[1]
        self.assertAlmostEqual(perfect["pcc"], 1.0)
        self.assertAlmostEqual(perfect["rmse"], math.sqrt(91 / 6))
        constant_prediction = samples.iloc[2]
        self.assertTrue(constant_prediction["constant_prediction"])
        self.assertEqual(constant_prediction["pcc"], 0.0)
        self.assertAlmostEqual(constant_prediction["rmse"], math.sqrt(91 / 6))
        self.assertEqual(list(samples["profile_id"]), ["a", "b", "c"])
        self.assertFalse(samples["nonfinite"].any())

    def test_condition_mean_averages_predicted_deltas_over_members(self):
        vector_metrics = learner.api(REPO)["vector_metrics"]
        condition = SimpleNamespace(
            delta=np.asarray([[2.0, 4.0, 6.0]], dtype=np.float32),
            members=[["a", "b"]],
            population=pd.DataFrame({"profile_id": ["c0"], "drug_parent": ["P"],
                                     "dose_um": [1.0], "cell_line": ["L"],
                                     "count": [2]}))
        predicted_delta = np.asarray([[1.0, 2.0, 3.0], [3.0, 6.0, 9.0]])
        samples = learner.score_condition_mean(condition, predicted_delta,
                                               {"a": 0, "b": 1}, vector_metrics)
        self.assertEqual(len(samples), 1)
        self.assertAlmostEqual(float(samples["pcc"].iloc[0]), 1.0)
        self.assertAlmostEqual(float(samples["rmse"].iloc[0]), 0.0)
        self.assertEqual(int(samples["count"].iloc[0]), 2)


@unittest.skipUnless(SOURCE.is_dir() and (SOURCE / "src/model.py").is_file(),
                     "the TranSiGen source is unavailable")
class ProbeRunTest(unittest.TestCase):
    """One bounded CPU probe through the whole formal path, then its verification."""

    @classmethod
    def setUpClass(cls):
        views.bind(REPO)
        cls._tmp = tempfile.TemporaryDirectory()
        cls.root = Path(cls._tmp.name)
        fixture = fixtures.dense_fixture(cls.root)
        accepted, _, _ = fixtures.build_view(fixture, size=SIZE)
        population = fixtures.treated_population(accepted, "tenv")
        cls.sealed = fixtures.write_sealed(cls.root, "tenv", population)
        build.build_product(view="full", source="tenv", output=cls.root / "views",
                            staging=cls.root / "staging", repo_root_path=REPO, size=SIZE,
                            fixture=True, corpus=fixture["corpus"],
                            axis=fixture["axis"], accepted=accepted,
                            provenance_root=None)
        cls.freeze = cls.root / "freezes/full/tenv"
        freeze.freeze_product(view="full", source="tenv", views_root=cls.root / "views",
                              destination=cls.freeze, repo_root_path=REPO,
                              sealed=cls.sealed, axis_dir=fixture["axis"],
                              accepted_dir=accepted)
        write_truth(cls.root / "truth", cls.freeze, "tenv", FOLD)

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def arguments(self, **overrides):
        value = dict(freezes=self.root / "freezes", truth=self.root / "truth",
                     native_source=SOURCE, device="cpu", mode="probe", probe_rows=4,
                     probe_epochs=1, only_source=("tenv",), only_view=("full",),
                     only_fold=(FOLD,), only_seed=(42,), repo_root_path=REPO)
        value.update(overrides)
        return value

    def _candidate(self, name):
        return self.root / name

    def _learner(self, candidate, phase, action, selection=None):
        context = formal.build_context(
            candidate, freezes=self.root / "freezes", truth=self.root / "truth",
            native_source=SOURCE, device="cpu", mode="probe",
            probe={"rows": 4, "epochs": 1}, repo_root_path=REPO,
            sealed=self.sealed)
        return formal.invoke(context, matrix.plan()[0], action, phase,
                             selection=selection, min_ram_kib=0,
                             telemetry_seconds=0.05)

    def test_driver_runs_receipt_test_replay_and_verifies(self):
        candidate = self._candidate("driver")
        report = formal.run_matrix(candidate, **self.arguments(), sealed=self.sealed,
                                   min_ram_kib=0, telemetry_seconds=0.05)
        self.assertEqual(report["complete"], 1)
        self.assertEqual(report["selected"], 1)
        origin = matrix.plan()[0]
        run = matrix.run_path(candidate, origin)
        self.assertEqual(json.loads((run / "terminal.json").read_text())["status"],
                         "complete")
        # The choice receipt exists before the test evaluation and is unchanged after.
        receipt = matrix.selection_path(candidate, "tenv", FOLD)
        self.assertTrue(receipt.is_file())
        self.assertEqual(json.loads(receipt.read_text()),
                         matrix.selection_evidence(candidate, "tenv", FOLD))
        metrics = json.loads((run / "test_metrics.json").read_text())
        self.assertEqual(metrics["n_primary_samples"], 2)
        self.assertGreaterEqual(metrics["n_condition_mean_samples"], 1)
        self.assertEqual(metrics["truth"], json.loads(
            (run / "config.json").read_text())["truth"])
        self.assertIn("selection_sha256", metrics)
        prediction = np.load(run / "test_predictions.npy")
        self.assertEqual(prediction.shape, (2, SIZE))
        self.assertEqual(prediction.dtype, np.float32)
        # Every phase wrote its 10-second telemetry and resource result.
        attempts = formal.attempt_records(candidate, origin)
        self.assertEqual(len(attempts), 4)
        self.assertTrue(all((candidate / "resources" / matrix.run_relative(origin).relative_to(
            "runs")) for _ in attempts))
        for attempt in (candidate / "resources" / matrix.run_relative(origin).relative_to(
                "runs")).iterdir():
            self.assertTrue((attempt / "telemetry.jsonl").is_file())
            self.assertTrue((attempt / "result.json").is_file())
            self.assertEqual(json.loads((attempt / "phase.json").read_text())["threads"], 2)
        # The per-run record carries seed, epochs, best epoch, hashes, device and cost.
        record = json.loads(matrix.record_path(candidate, origin).read_text())
        self.assertEqual(record["schema"], matrix.RECORD_SCHEMA)
        self.assertEqual(record["seed"], 42)
        self.assertEqual(record["realized_epochs"], 1)
        self.assertEqual(record["best_epoch"], 1)
        self.assertEqual(record["device"], "cpu")
        self.assertIn("test_predictions.npy", record["artifacts"])
        self.assertGreater(record["peak_rss_kib"], 0)
        self.assertGreaterEqual(record["runtime_seconds"], 0)
        self.assertEqual(record["truth"]["primary"]["manifest_sha256"],
                         metrics["truth"]["primary"]["manifest_sha256"])
        # Verify re-derives the plan, receipt, tokens and reports the full-arm drift.
        verified = formal.verify(candidate, **self.arguments())
        self.assertEqual(verified["schema"], formal.VERIFY_SCHEMA)
        self.assertEqual(verified["origins"], 1)
        self.assertEqual(len(verified["full_arm"]["origins"]), 1)
        self.assertIn(verified["full_arm"]["origins"][0]["status"], ("PASS", "DRIFT"))
        self.assertIn(verified["status"], ("PASS", "DRIFT"))

    def test_resume_after_completion_runs_no_phase_twice(self):
        candidate = self._candidate("resume")
        formal.run_matrix(candidate, **self.arguments(), sealed=self.sealed,
                          min_ram_kib=0, telemetry_seconds=0.05)
        origin = matrix.plan()[0]
        attempts = lambda: len(formal.attempt_records(candidate, origin))
        before = attempts()
        report = formal.run_matrix(candidate, **self.arguments(), sealed=self.sealed,
                                   min_ram_kib=0, telemetry_seconds=0.05)
        self.assertEqual(report["complete"], 1)
        self.assertEqual(attempts(), before)

    def test_train_reads_no_test_partition_and_no_truth_product(self):
        candidate = self._candidate("isolation")
        plan = formal.plan_for(**self.arguments())
        matrix.seal_plan(candidate, plan)
        run = matrix.run_path(candidate, matrix.plan()[0])
        partitions, original = [], learner.partition_rows

        def spy(data, scenario, fold, partition, smoke_rows):
            partitions.append(partition)
            return original(data, scenario, fold, partition, smoke_rows)

        with mock.patch.object(learner, "partition_rows", side_effect=spy), \
                mock.patch.object(learner, "load_truth",
                                  side_effect=AssertionError("train read the truth")):
            terminal = learner.train(self.freeze, SOURCE, run, rows=4, epochs=1,
                                     device="cpu", fold=FOLD, seed=42,
                                     repo_root_path=REPO, sealed=self.sealed,
                                     mode="probe", matrix_root=candidate)
        self.assertEqual(sorted(partitions), ["train", "validation"])
        self.assertNotIn("test", partitions)
        self.assertEqual(terminal["status"], "validation-trained-awaiting-replay")
        self.assertEqual(terminal["mode"], "probe")
        # The checkpoint rule restores the best epoch, and the curve agrees.
        curve = pd.read_csv(run / "learning_curve.csv")
        best = curve.loc[curve["validation_mse"].idxmin()]
        self.assertEqual(terminal["best_epoch"], int(best["epoch"]))
        self.assertEqual(terminal["validation_mse"], float(best["validation_mse"]))
        saved = __import__("torch").load(run / "checkpoint.pt", map_location="cpu",
                                         weights_only=False)
        self.assertEqual(saved["best_epoch"], terminal["best_epoch"])
        self.assertEqual(saved["validation_mse"], terminal["validation_mse"])
        prediction = np.load(run / "validation_prediction.npy")
        self.assertEqual(prediction.dtype, np.float32)
        self.assertEqual(prediction.shape[1], SIZE)
        # A fresh-process validation replay seals the choice without touching truth.
        self._learner(candidate, "validation", "replay")
        self.assertEqual(json.loads((run / "terminal.json").read_text())["status"],
                         "validation-complete")
        # Test refuses before any truth read: no receipt, then a forgeried receipt.
        with mock.patch.object(learner, "_truth_rows",
                               side_effect=AssertionError("read truth before receipt")):
            with self.assertRaisesRegex(ValueError, "requires a sealed validation"):
                learner.test(self.freeze, SOURCE, run, repo_root_path=REPO,
                             sealed=self.sealed, selection=None,
                             truth_root=self.root / "truth", matrix_root=candidate)
            bogus = matrix.selection_path(candidate, "tenv", FOLD)
            bogus.parent.mkdir(parents=True, exist_ok=True)
            bogus.write_text(json.dumps({"schema": "bogus"}))
            with self.assertRaisesRegex(ValueError, "selection receipt changed"):
                learner.test(self.freeze, SOURCE, run, repo_root_path=REPO,
                             sealed=self.sealed, selection=bogus,
                             truth_root=self.root / "truth", matrix_root=candidate)
        matrix.selection_path(candidate, "tenv", FOLD).unlink()
        receipt = matrix.lock_selection(candidate, "tenv", FOLD)
        self._learner(candidate, "test", "evaluate", selection=receipt)
        terminal = json.loads((run / "terminal.json").read_text())
        self.assertEqual(terminal["status"], "test-awaiting-replay")
        self.assertEqual(terminal["test_selection_sha256"], views.sha256(receipt))
        stored = np.load(run / "test_predictions.npy")
        self.assertEqual(stored.shape, (2, SIZE))
        # The test replay reproduces the prediction bit-exactly in a fresh process.
        self._learner(candidate, "test", "replay", selection=receipt)
        self.assertEqual(json.loads((run / "terminal.json").read_text())["status"],
                         "complete")
        self.assertTrue(np.array_equal(np.load(run / "test_predictions.npy"), stored))
        tokens = json.loads((run / "terminal.json").read_text())
        self.assertEqual(len({tokens["training_process_token"],
                              tokens["validation_replay_token"],
                              tokens["test_evaluation_process_token"],
                              tokens["test_replay_token"]}), 4)
        # Moving the truth product away makes the final replay fail, not silently pass.
        moved = self.root / "truth-aside"
        (self.root / "truth").rename(moved)
        try:
            failed = subprocess.run(
                [sys.executable, str(EXECUTION / "learner.py"), str(self.freeze),
                 str(SOURCE), str(run), "replay", "--mode", "probe", "--device", "cpu",
                 "--fold", str(FOLD), "--seed", "42", "--matrix-root", str(candidate),
                 "--truth", str(self.root / "truth"), "--selection", str(receipt),
                 "--repo-root", str(REPO), "--sealed", str(self.sealed)],
                capture_output=True, text=True)
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn("truth", failed.stderr.lower())
        finally:
            moved.rename(self.root / "truth")


class TruthBlockColumnTest(unittest.TestCase):
    """The truth block reads only the columns truth.py emits.

    The emitter publishes ``drug_parent`` (primary) and ``members`` (condition-mean). The
    block used to accept two undocumented spellings of each; a truth product written under
    a different name is not the frozen product, so it must fail closed instead of being
    silently reinterpreted.
    """

    def _block(self, root, columns, kind, value):
        import hashlib
        directory = Path(root) / kind
        directory.mkdir(parents=True)
        frame = pd.DataFrame([value], columns=columns)
        frame.to_csv(directory / "population.csv", index=False)
        rng = np.random.default_rng(0)
        delta = rng.standard_normal((1, SIZE)).astype(np.float32)
        np.save(directory / "delta.npy", delta)
        if kind == "primary":
            np.save(directory / "control.npy", delta)
            pd.DataFrame({"drug_parent_inchikey": ["PARENT"]}).to_csv(
                directory / "parents.csv", index=False)
            np.save(directory / "morgan.npy",
                    np.zeros((1, 2048), dtype=np.uint8))
        population_sha = hashlib.sha256(
            (directory / "population.csv").read_bytes()).hexdigest()
        _write_json(directory / "manifest.json", {
            "schema": matrix.TRUTH_SCHEMA, "source": "tenv", "kind": kind, "fold": 0,
            "rows": 1, "outputs": {"population.csv": {"sha256": population_sha},
                                   "delta.npy": {"sha256": hashlib.sha256(
                                       (directory / "delta.npy").read_bytes()).hexdigest()}}})
        return directory

    def _primary(self, root, parent_column):
        return self._block(
            root, ["profile_id", "cell_line", parent_column, "dose_um"], "primary",
            {"profile_id": "p0", "cell_line": "line", parent_column: "PARENT",
             "dose_um": 1.0})

    def test_canonical_primary_columns_are_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = self._primary(tmp, "drug_parent")
            block = learner.TruthBlock(directory, "tenv", 0, "primary")
            self.assertEqual(block.population["drug_parent"].tolist(), ["PARENT"])

    def test_an_alternate_parent_spelling_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = self._primary(tmp, "parent")
            with self.assertRaisesRegex(learner.ScientificFailure, "drug_parent"):
                learner.TruthBlock(directory, "tenv", 0, "primary")

    def test_an_alternate_members_spelling_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = self._block(
                tmp, ["profile_id", "cell_line", "drug_parent", "dose_um",
                      "member_profile_ids"], "condition_mean",
                {"profile_id": "c0", "cell_line": "line", "drug_parent": "PARENT",
                 "dose_um": 1.0, "member_profile_ids": "p0"})
            with self.assertRaisesRegex(learner.ScientificFailure, "members"):
                learner.TruthBlock(directory, "tenv", 0, "condition_mean")

    def test_the_canonical_members_column_is_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = self._block(
                tmp, ["profile_id", "cell_line", "drug_parent", "dose_um", "members"],
                "condition_mean",
                {"profile_id": "c0", "cell_line": "line", "drug_parent": "PARENT",
                 "dose_um": 1.0, "members": "p0"})
            block = learner.TruthBlock(directory, "tenv", 0, "condition_mean")
            self.assertEqual(block.members.tolist(), [["p0"]])


if __name__ == "__main__":
    unittest.main()
