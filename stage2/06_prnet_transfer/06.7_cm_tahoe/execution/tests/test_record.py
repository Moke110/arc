"""Unit tests for the ticket-#7 record writer, the warm-start gates and the queue.

Each test builds a synthetic 45-run output tree in a temporary directory (real
frozen ``locked_config.json`` with a tiny model structure, real tiny torch
checkpoints for the warm starts, synthetic run manifests) and then mutates it to
prove the checker notices.  No GPU, no expression matrix, no network.
"""

from __future__ import annotations

import contextlib
import csv as _csv
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import _support  # noqa: E402
import build_queue  # noqa: E402
import write_record  # noqa: E402

from ptlib import config as pt_config  # noqa: E402


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
        write_record.write_record(
            data_dir=fixture["data_dir"],
            out_dir=fixture["out_dir"],
            random_record_path=fixture["random_record_path"],
        )
        return write_record.check_record(
            data_dir=fixture["data_dir"],
            out_dir=fixture["out_dir"],
            random_record_path=fixture["random_record_path"],
        )

    def _record_row(self, fixture, unit_id):
        rows = fixture["out_dir"] / "record" / "candidate_record.csv"
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
        self.assertEqual(len(report["warmstart_provenance_sha256"]), 64)
        # the default fixture prefers 1e-4 in every fold
        for selection in report["selection"]:
            self.assertEqual(selection["lr"], "0.0001")
        # every gate family is present, including the ticket-#7 warm-start gates
        gate_ids = {gate["id"] for gate in report["checks"]}
        for required in (
            "warmstart_verified", "warmstart_provenance_complete", "per_run_init_provenance",
            "init_checkpoint_bytes_match", "init_checkpoints_distinct",
            "seed_checkpoint_one_to_one", "matrix_symmetry_with_ticket6",
        ):
            self.assertIn(required, gate_ids)

    def test_candidate_rows_record_init_identity(self):
        fixture = self._fixture()
        write_record.write_record(
            data_dir=fixture["data_dir"], out_dir=fixture["out_dir"],
            random_record_path=fixture["random_record_path"],
        )
        row = self._record_row(fixture, "cm_f0_s20261021_lr1e-4_tahoe")
        self.assertEqual(row["status"], "ok")
        self.assertEqual(row["init"], "tahoe")
        self.assertEqual(row["selected"], "1")
        self.assertEqual(len(row["config_sha256"]), 64)
        self.assertEqual(len(row["checkpoint_sha256"]), 64)
        self.assertEqual(len(row["init_checkpoint_sha256"]), 64)
        self.assertEqual(len(row["init_checkpoint_config_sha256"]), 64)
        self.assertEqual(row["init_checkpoint_seed"], "20261021")
        self.assertEqual(row["init_weights_identical"], "1")
        self.assertGreater(int(row["init_n_loaded_tensors"]), 0)
        self.assertIn("encoder.mean_encoder.weight", row["init_loaded_layers"])
        self.assertTrue(row["init_checkpoint_worktree_relpath"].endswith("model.pt"))

    def test_fifteen_selected_three_per_fold_at_the_chosen_lr(self):
        fixture = self._fixture()
        write_record.write_record(
            data_dir=fixture["data_dir"], out_dir=fixture["out_dir"],
            random_record_path=fixture["random_record_path"],
        )
        with open(fixture["out_dir"] / "record" / "selected_checkpoints.csv",
                  newline="", encoding="utf-8") as handle:
            rows = list(_csv.DictReader(handle))
        self.assertEqual(len(rows), 15)
        for fold in range(5):
            fold_rows = [row for row in rows if int(row["fold"]) == fold]
            self.assertEqual(len(fold_rows), 3)
            self.assertEqual(sorted(int(row["seed"]) for row in fold_rows),
                             [20261021, 20261022, 20261023])
            expected = fixture["checkpoints"][str(fold_rows[0]["seed"])]
            self.assertEqual(fold_rows[0]["init_checkpoint_sha256"], expected["sha256"])

    def test_record_content_hash_is_stable_across_rebuilds(self):
        fixture = self._fixture()
        kwargs = dict(data_dir=fixture["data_dir"], out_dir=fixture["out_dir"],
                      random_record_path=fixture["random_record_path"])
        first = write_record.build_record(**kwargs)
        second = write_record.build_record(**kwargs)
        self.assertEqual(first["content_sha256"], second["content_sha256"])

    def test_capped_candidate_is_recorded_but_not_a_failure(self):
        fixture = self._fixture()
        capped_unit = "cm_f4_s20261023_lr1e-3_tahoe"
        _support.write_run(
            fixture["runs_root"],
            next(c for c in fixture["candidates"] if c.unit_id == capped_unit),
            fixture=fixture,
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
        with open(fixture["out_dir"] / "record" / "lr_selection.csv",
                  newline="", encoding="utf-8") as handle:
            rows = list(_csv.DictReader(handle))
        self.assertTrue(all(row["tie_break_applied"] == "1" for row in rows))


class WarmStartGateTests(RecordTestBase):
    def _manifest_path(self, fixture, unit_id):
        return fixture["runs_root"] / unit_id / "manifest.json"

    def test_random_init_run_is_detected(self):
        fixture = self._fixture()
        unit_id = "cm_f0_s20261021_lr1e-4_tahoe"
        path = self._manifest_path(fixture, unit_id)
        manifest = json.loads(path.read_text(encoding="utf-8"))
        manifest["config"]["init"] = "random"
        manifest["config_sha256"] = pt_config.config_fingerprint(manifest["config"])
        path.write_text(json.dumps(manifest), encoding="utf-8")
        passed, report = self._write_and_check(fixture)
        self.assertFalse(passed)
        self.assertIn("no_failed_candidates", _failed_ids(report))
        self.assertIn("pretrained", self._record_row(fixture, unit_id)["detail"])

    def test_wrong_seed_checkpoint_is_detected(self):
        fixture = self._fixture()
        unit_id = "cm_f1_s20261022_lr3e-4_tahoe"
        path = self._manifest_path(fixture, unit_id)
        manifest = json.loads(path.read_text(encoding="utf-8"))
        other_seed = "20261021"
        manifest["init"]["source_path"] = fixture["checkpoints"][other_seed]["path"]
        manifest["init"]["checkpoint_sha256"] = fixture["checkpoints"][other_seed]["sha256"]
        path.write_text(json.dumps(manifest), encoding="utf-8")
        passed, report = self._write_and_check(fixture)
        self.assertFalse(passed)
        self.assertIn("per_run_init_provenance", _failed_ids(report))

    def test_partial_load_provenance_is_detected(self):
        fixture = self._fixture()
        unit_id = "cm_f2_s20261023_lr1e-3_tahoe"
        path = self._manifest_path(fixture, unit_id)
        manifest = json.loads(path.read_text(encoding="utf-8"))
        manifest["init"]["load_provenance"]["missing_keys"] = ["decoder.recon_decoder.weight"]
        path.write_text(json.dumps(manifest), encoding="utf-8")
        passed, report = self._write_and_check(fixture)
        self.assertFalse(passed)
        self.assertIn("per_run_init_provenance", _failed_ids(report))

    def test_missing_init_record_is_detected(self):
        fixture = self._fixture()
        unit_id = "cm_f3_s20261021_lr1e-4_tahoe"
        path = self._manifest_path(fixture, unit_id)
        manifest = json.loads(path.read_text(encoding="utf-8"))
        manifest.pop("init")
        path.write_text(json.dumps(manifest), encoding="utf-8")
        passed, report = self._write_and_check(fixture)
        self.assertFalse(passed)
        self.assertIn("per_run_init_provenance", _failed_ids(report))

    def test_changed_warmstart_checkpoint_bytes_are_detected(self):
        fixture = self._fixture()
        victim = fixture["checkpoints"]["20261021"]["path"]
        Path(victim).write_bytes(b"tampered checkpoint")
        passed, report = self._write_and_check(fixture)
        self.assertFalse(passed)
        self.assertIn("init_checkpoint_bytes_match", _failed_ids(report))

    def test_tampered_queue_warmstart_block_is_detected(self):
        fixture = self._fixture()
        queue_path = fixture["out_dir"] / "units.json"
        queue = json.loads(queue_path.read_text(encoding="utf-8"))
        queue["warmstart"]["checkpoints"]["20261021"]["sha256"] = "0" * 64
        queue_path.write_text(json.dumps(queue), encoding="utf-8")
        passed, report = self._write_and_check(fixture)
        self.assertFalse(passed)
        self.assertIn("per_run_init_provenance", _failed_ids(report))

    def test_incomplete_warmstart_block_is_detected(self):
        fixture = self._fixture()
        queue_path = fixture["out_dir"] / "units.json"
        queue = json.loads(queue_path.read_text(encoding="utf-8"))
        queue["warmstart"]["checkpoints"].pop("20261023")
        queue["warmstart"]["all_verified"] = False
        queue_path.write_text(json.dumps(queue), encoding="utf-8")
        write_record.write_record(
            data_dir=fixture["data_dir"], out_dir=fixture["out_dir"],
            random_record_path=fixture["random_record_path"],
        )
        passed, report = write_record.check_record(
            data_dir=fixture["data_dir"], out_dir=fixture["out_dir"],
            random_record_path=fixture["random_record_path"],
        )
        self.assertFalse(passed)
        failed = _failed_ids(report)
        self.assertIn("warmstart_verified", failed)
        self.assertIn("init_checkpoints_distinct", failed)


class SymmetryTests(RecordTestBase):
    def test_differing_matrix_fails_the_symmetry_gate(self):
        fixture = self._fixture()
        random_record = json.loads(
            Path(fixture["random_record_path"]).read_text(encoding="utf-8")
        )
        random_record["matrix"]["lrs"] = [1e-4, 3e-4, 1e-2]
        broken = self.tmp / "broken_random_record.json"
        broken.write_text(json.dumps(random_record), encoding="utf-8")
        write_record.write_record(
            data_dir=fixture["data_dir"], out_dir=fixture["out_dir"],
            random_record_path=broken,
        )
        passed, report = write_record.check_record(
            data_dir=fixture["data_dir"], out_dir=fixture["out_dir"],
            random_record_path=broken,
        )
        self.assertFalse(passed)
        self.assertIn("matrix_symmetry_with_ticket6", _failed_ids(report))

    def test_missing_ticket6_record_fails_the_symmetry_gate(self):
        fixture = self._fixture()
        write_record.write_record(
            data_dir=fixture["data_dir"], out_dir=fixture["out_dir"],
            random_record_path=None,
        )
        passed, report = write_record.check_record(
            data_dir=fixture["data_dir"], out_dir=fixture["out_dir"],
            random_record_path=None,
        )
        self.assertFalse(passed)
        self.assertIn("matrix_symmetry_with_ticket6", _failed_ids(report))


class FailureDetectionTests(RecordTestBase):
    def test_failed_candidate_fails_the_matrix_gate(self):
        victim = "cm_f2_s20261022_lr3e-4_tahoe"
        fixture = self._fixture(failed=frozenset({victim}))
        passed, report = self._write_and_check(fixture)
        self.assertFalse(passed)
        self.assertIn("no_failed_candidates", _failed_ids(report))
        self.assertEqual(report["status_counts"]["failed"], 1)

    def test_checkpoint_tamper_is_detected(self):
        fixture = self._fixture()
        write_record.write_record(
            data_dir=fixture["data_dir"], out_dir=fixture["out_dir"],
            random_record_path=fixture["random_record_path"],
        )
        victim = fixture["runs_root"] / "cm_f1_s20261021_lr1e-4_tahoe" / "checkpoint" / "model.pt"
        victim.write_bytes(b"tampered")
        passed, report = write_record.check_record(
            data_dir=fixture["data_dir"], out_dir=fixture["out_dir"],
            random_record_path=fixture["random_record_path"],
        )
        self.assertFalse(passed)
        self.assertIn("no_failed_candidates", _failed_ids(report))

    def test_data_identity_drift_is_detected(self):
        fixture = self._fixture()
        manifest_path = (
            fixture["runs_root"] / "cm_f0_s20261021_lr1e-4_tahoe" / "manifest.json"
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
            fixture["runs_root"] / "cm_f0_s20261021_lr1e-4_tahoe" / "manifest.json"
        )
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["config"]["recipe"]["max_epochs"] = 55
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
        passed, report = self._write_and_check(fixture)
        self.assertFalse(passed)
        self.assertIn("no_failed_candidates", _failed_ids(report))

    def test_prediction_artifacts_are_detected(self):
        fixture = self._fixture()
        victim = fixture["runs_root"] / "cm_f3_s20261021_lr1e-4_tahoe"
        (victim / "evaluation" / "test").mkdir(parents=True)
        (victim / "evaluation" / "test" / "aggregate.json").write_text("{}", encoding="utf-8")
        passed, report = self._write_and_check(fixture)
        self.assertFalse(passed)
        self.assertIn("no_failed_candidates", _failed_ids(report))

    def test_record_csv_tamper_is_detected(self):
        fixture = self._fixture()
        write_record.write_record(
            data_dir=fixture["data_dir"], out_dir=fixture["out_dir"],
            random_record_path=fixture["random_record_path"],
        )
        csv_path = fixture["out_dir"] / "record" / "candidate_record.csv"
        csv_path.write_text(csv_path.read_text(encoding="utf-8") + "tampered\n", encoding="utf-8")
        passed, report = write_record.check_record(
            data_dir=fixture["data_dir"], out_dir=fixture["out_dir"],
            random_record_path=fixture["random_record_path"],
        )
        self.assertFalse(passed)
        self.assertIn("record_artifacts_match", _failed_ids(report))

    def test_missing_matrix_run_fails_the_run_gate(self):
        fixture = self._fixture()
        write_record.write_record(
            data_dir=fixture["data_dir"], out_dir=fixture["out_dir"],
            random_record_path=fixture["random_record_path"],
        )
        (fixture["out_dir"] / "record" / "matrix_run.json").unlink()
        passed, report = write_record.check_record(
            data_dir=fixture["data_dir"], out_dir=fixture["out_dir"],
            random_record_path=fixture["random_record_path"],
        )
        self.assertFalse(passed)
        self.assertIn("matrix_run_recorded", _failed_ids(report))

    def test_missing_resume_check_fails_the_resume_gate(self):
        fixture = self._fixture()
        write_record.write_record(
            data_dir=fixture["data_dir"], out_dir=fixture["out_dir"],
            random_record_path=fixture["random_record_path"],
        )
        (fixture["out_dir"] / "record" / "resume_check.json").unlink()
        passed, report = write_record.check_record(
            data_dir=fixture["data_dir"], out_dir=fixture["out_dir"],
            random_record_path=fixture["random_record_path"],
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
        fixture = _support.build_queue_fixture(self.tmp)
        payload = fixture["queue"]
        self.assertEqual(payload["schema"], "cm_tahoe.queue/v1")
        units = payload["units"]
        self.assertEqual(len(units), 45)
        self.assertTrue(payload["complete"])
        self.assertFalse(payload["profile"])
        self.assertEqual(len({unit["unit_id"] for unit in units}), 45)
        self.assertEqual(sorted({unit["fold"] for unit in units}), [0, 1, 2, 3, 4])
        self.assertEqual(sorted({unit["seed"] for unit in units}),
                         [20261021, 20261022, 20261023])
        self.assertEqual(sorted({unit["lr"] for unit in units}), [1e-4, 3e-4, 1e-3])
        self.assertEqual(payload["warmstart"]["all_verified"], True)
        self.assertEqual(len(payload["warmstart"]["checkpoints"]), 3)
        for unit in units:
            self.assertEqual(unit["stage"], "finetune")
            self.assertEqual(unit["init"], "tahoe")
            self.assertEqual(unit["pipeline_init"], "pretrained")
            self.assertIn("--init", unit["args"])
            self.assertIn("pretrained", unit["args"])
            self.assertIn("--init-checkpoint", unit["args"])
            expected = fixture["checkpoints"][str(unit["seed"])]
            self.assertEqual(unit["init_checkpoint_sha256"], expected["sha256"])
            self.assertIn(expected["path"], unit["args"])
            self.assertNotIn("--out", unit["args"])
            self.assertNotIn("--data-dir", unit["args"])
            # a formal unit never carries the bounded profile budget
            self.assertNotIn("--epochs", unit["args"])

    def test_profile_queue_is_one_bounded_unit(self):
        fixture = _support.build_queue_fixture(self.tmp)
        profile_queue = fixture["out_dir"] / "profile" / "profile_queue.json"
        payload = build_queue.build_queue(
            data_dir=fixture["data_dir"],
            out_dir=fixture["out_dir"],
            checkpoint_root=fixture["checkpoint_root"],
            queue_path=profile_queue,
            hash_cache=None,
            profile=True,
            profile_epochs=2,
        )
        self.assertFalse(payload["complete"])
        self.assertTrue(payload["profile"])
        self.assertEqual(len(payload["units"]), 1)
        unit = payload["units"][0]
        self.assertEqual(unit["fold"], 0)
        self.assertEqual(unit["seed"], 20261021)
        self.assertEqual(unit["lr"], 1e-3)
        self.assertIn("--epochs", unit["args"])
        self.assertEqual(unit["args"][unit["args"].index("--epochs") + 1], "2")
        self.assertEqual(payload["profile_plan"]["epochs"], 2)

    def test_missing_checkpoint_refuses_to_write_a_queue(self):
        data_dir, _locked, _frozen = _support.prepare_input(
            self.tmp, structure=_support.SMALL_STRUCTURE
        )
        out_dir = self.tmp / "output_none"
        with self.assertRaises(Exception):
            build_queue.build_queue(
                data_dir=data_dir,
                out_dir=out_dir,
                checkpoint_root=self.tmp / "no_checkpoints",
                queue_path=out_dir / "units.json",
                hash_cache=None,
            )
        self.assertFalse((out_dir / "units.json").exists())


if __name__ == "__main__":
    unittest.main()
