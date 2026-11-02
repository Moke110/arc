"""Unit tests for the ticket-#6 record writer and its one-to-one / no-test-read gates.

Each test builds a synthetic 45-run output tree in a temporary directory (real
frozen `locked_config.json`, synthetic manifests/checkpoints) and then mutates it
to prove the checker notices.  No GPU, no expression matrix, no network.
"""

from __future__ import annotations

import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import _support  # noqa: E402
import write_record  # noqa: E402


def _failed_ids(report: dict) -> set[str]:
    return {gate["id"] for gate in report["checks"] if not gate["ok"]}


class RecordTestBase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def _fixture(self, **kwargs):
        return _support.prepare_out(self.tmp, **kwargs)

    def _write_and_check(self, fixture):
        write_record.write_record(data_dir=fixture["data_dir"], out_dir=fixture["out_dir"])
        return write_record.check_record(data_dir=fixture["data_dir"], out_dir=fixture["out_dir"])

    def _record_row(self, fixture, unit_id):
        spec = fixture["spec"]
        rows = fixture["out_dir"] / "record" / "candidate_record.csv"
        import csv as _csv

        with open(rows, newline="", encoding="utf-8") as handle:
            for row in _csv.DictReader(handle):
                if row["unit_id"] == unit_id:
                    return row
        raise AssertionError(f"{unit_id} missing from candidate_record.csv")


class CleanMatrixTests(RecordTestBase):
    def test_clean_matrix_passes_every_gate(self):
        fixture = self._fixture()
        passed, report = self._write_and_check(fixture)
        self.assertTrue(passed, _failed_ids(report))
        self.assertEqual(report["status_counts"], {"ok": 45, "capped": 0, "failed": 0})
        self.assertEqual(report["n_selected"], 15)
        # the default fixture prefers 1e-4 in every fold
        for selection in report["selection"]:
            self.assertEqual(selection["lr"], "0.0001")

    def test_candidate_rows_record_identity_and_selection(self):
        fixture = self._fixture()
        write_record.write_record(data_dir=fixture["data_dir"], out_dir=fixture["out_dir"])
        row = self._record_row(fixture, "cm_f0_s20261021_lr1e-4_random")
        self.assertEqual(row["status"], "ok")
        self.assertEqual(row["fold"], "0")
        self.assertEqual(row["seed"], "20261021")
        self.assertEqual(row["lr"], "0.0001")
        self.assertEqual(row["selected"], "1")
        self.assertEqual(len(row["config_sha256"]), 64)
        self.assertEqual(len(row["checkpoint_sha256"]), 64)
        self.assertEqual(len(row["data_identity"]), 64)
        self.assertTrue(row["checkpoint_worktree_relpath"].endswith("checkpoint/model.pt"))

    def test_fifteen_selected_three_per_fold_at_the_chosen_lr(self):
        fixture = self._fixture()
        write_record.write_record(data_dir=fixture["data_dir"], out_dir=fixture["out_dir"])
        import csv as _csv

        with open(fixture["out_dir"] / "record" / "selected_checkpoints.csv",
                  newline="", encoding="utf-8") as handle:
            rows = list(_csv.DictReader(handle))
        self.assertEqual(len(rows), 15)
        for fold in range(5):
            fold_rows = [row for row in rows if int(row["fold"]) == fold]
            self.assertEqual(len(fold_rows), 3)
            self.assertEqual(sorted(int(row["seed"]) for row in fold_rows),
                             [20261021, 20261022, 20261023])

    def test_record_content_hash_is_stable_across_rebuilds(self):
        fixture = self._fixture()
        first = write_record.build_record(data_dir=fixture["data_dir"], out_dir=fixture["out_dir"])
        second = write_record.build_record(data_dir=fixture["data_dir"], out_dir=fixture["out_dir"])
        self.assertEqual(first["content_sha256"], second["content_sha256"])

    def test_capped_candidate_is_recorded_but_not_a_failure(self):
        fixture = self._fixture()
        capped_unit = "cm_f4_s20261023_lr1e-3_random"
        _support.write_run(
            fixture["runs_root"],
            next(c for c in fixture["candidates"] if c.unit_id == capped_unit),
            locked=fixture["locked"], frozen=fixture["frozen"],
            auprc=_support.default_auprc(
                next(c for c in fixture["candidates"] if c.unit_id == capped_unit)
            ),
            capped=True,
            epochs=100,
        )
        passed, report = self._write_and_check(fixture)
        self.assertTrue(passed, _failed_ids(report))
        self.assertEqual(report["status_counts"]["capped"], 1)
        self.assertEqual(self._record_row(fixture, capped_unit)["status"], "capped")
        self.assertEqual(self._record_row(fixture, capped_unit)["epochs_run"], "100")


class SelectionTests(RecordTestBase):
    def test_tie_break_prefers_lower_lr(self):
        def auprc(candidate):
            base = {1e-4: 0.12, 3e-4: 0.12, 1e-3: 0.01}[candidate.lr]
            return base + 0.0 * candidate.seed_index

        fixture = self._fixture(auprc=auprc)
        passed, report = self._write_and_check(fixture)
        self.assertTrue(passed, _failed_ids(report))
        for selection in report["selection"]:
            self.assertEqual(selection["lr"], "0.0001")
        import csv as _csv

        with open(fixture["out_dir"] / "record" / "lr_selection.csv",
                  newline="", encoding="utf-8") as handle:
            rows = list(_csv.DictReader(handle))
        self.assertTrue(all(row["tie_break_applied"] == "1" for row in rows))


class FailureDetectionTests(RecordTestBase):
    def test_failed_candidate_fails_the_matrix_gate(self):
        victim = "cm_f2_s20261022_lr3e-4_random"
        fixture = self._fixture(failed=frozenset({victim}))
        passed, report = self._write_and_check(fixture)
        self.assertFalse(passed)
        self.assertIn("no_failed_candidates", _failed_ids(report))
        self.assertEqual(report["status_counts"]["failed"], 1)
        self.assertEqual(self._record_row(fixture, victim)["status"], "failed")

    def test_checkpoint_tamper_is_detected(self):
        fixture = self._fixture()
        write_record.write_record(data_dir=fixture["data_dir"], out_dir=fixture["out_dir"])
        victim = fixture["runs_root"] / "cm_f1_s20261021_lr1e-4_random" / "checkpoint" / "model.pt"
        victim.write_bytes(b"tampered")
        passed, report = write_record.check_record(
            data_dir=fixture["data_dir"], out_dir=fixture["out_dir"]
        )
        self.assertFalse(passed)
        self.assertIn("no_failed_candidates", _failed_ids(report))

    def test_data_identity_drift_is_detected(self):
        fixture = self._fixture()
        manifest_path = (
            fixture["runs_root"] / "cm_f0_s20261021_lr1e-4_random" / "manifest.json"
        )
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["data"]["artifacts"]["cm_transfer.h5ad"] = "0" * 64
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
        passed, report = self._write_and_check(fixture)
        self.assertFalse(passed)
        self.assertIn("no_failed_candidates", _failed_ids(report))

    def test_config_drift_is_detected(self):
        fixture = self._fixture()
        manifest_path = (
            fixture["runs_root"] / "cm_f0_s20261021_lr1e-4_random" / "manifest.json"
        )
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["config"]["recipe"]["max_epochs"] = 55
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
        passed, report = self._write_and_check(fixture)
        self.assertFalse(passed)
        self.assertIn("no_failed_candidates", _failed_ids(report))

    def test_prediction_artifacts_are_detected(self):
        fixture = self._fixture()
        victim = fixture["runs_root"] / "cm_f3_s20261021_lr1e-4_random"
        (victim / "evaluation" / "test").mkdir(parents=True)
        (victim / "evaluation" / "test" / "aggregate.json").write_text("{}", encoding="utf-8")
        passed, report = self._write_and_check(fixture)
        self.assertFalse(passed)
        self.assertIn("no_failed_candidates", _failed_ids(report))

    def test_record_csv_tamper_is_detected(self):
        fixture = self._fixture()
        write_record.write_record(data_dir=fixture["data_dir"], out_dir=fixture["out_dir"])
        csv_path = fixture["out_dir"] / "record" / "candidate_record.csv"
        csv_path.write_text(csv_path.read_text(encoding="utf-8") + "tampered\n", encoding="utf-8")
        passed, report = write_record.check_record(
            data_dir=fixture["data_dir"], out_dir=fixture["out_dir"]
        )
        self.assertFalse(passed)
        self.assertIn("record_artifacts_match", _failed_ids(report))

    def test_missing_matrix_run_fails_the_run_gate(self):
        fixture = self._fixture()
        write_record.write_record(data_dir=fixture["data_dir"], out_dir=fixture["out_dir"])
        (fixture["out_dir"] / "record" / "matrix_run.json").unlink()
        passed, report = write_record.check_record(
            data_dir=fixture["data_dir"], out_dir=fixture["out_dir"]
        )
        self.assertFalse(passed)
        self.assertIn("matrix_run_recorded", _failed_ids(report))

    def test_missing_resume_check_fails_the_resume_gate(self):
        fixture = self._fixture()
        write_record.write_record(data_dir=fixture["data_dir"], out_dir=fixture["out_dir"])
        (fixture["out_dir"] / "record" / "resume_check.json").unlink()
        passed, report = write_record.check_record(
            data_dir=fixture["data_dir"], out_dir=fixture["out_dir"]
        )
        self.assertFalse(passed)
        self.assertIn("resume_check_passed", _failed_ids(report))


class BuildQueueTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def test_forty_five_units_from_the_frozen_config(self):
        data_dir, _frozen = _support.prepare_input(self.tmp)
        out_dir = self.tmp / "output"
        import contextlib
        import io

        import build_queue

        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(build_queue.main([
                "--data-dir", str(data_dir), "--out-dir", str(out_dir),
                "--hash-cache", "none",
            ]), 0)
        payload = json.loads((out_dir / "units.json").read_text(encoding="utf-8"))
        self.assertEqual(payload["schema"], "cm_random.queue/v1")
        units = payload["units"]
        self.assertEqual(len(units), 45)
        self.assertEqual(len({unit["unit_id"] for unit in units}), 45)
        self.assertEqual(sorted({unit["fold"] for unit in units}), [0, 1, 2, 3, 4])
        self.assertEqual(sorted({unit["seed"] for unit in units}),
                         [20261021, 20261022, 20261023])
        self.assertEqual(sorted({unit["lr"] for unit in units}), [1e-4, 3e-4, 1e-3])
        for unit in units:
            self.assertEqual(unit["stage"], "finetune")
            self.assertEqual(unit["init"], "random")
            self.assertIn("--init", unit["args"])
            self.assertNotIn("--out", unit["args"])
            self.assertNotIn("--data-dir", unit["args"])


if __name__ == "__main__":
    unittest.main()
