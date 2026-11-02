"""Offline tests for the 06.5 record writer/reader/checker (stdlib unittest).

Run with the `arc` environment (no GPU, no torch, no frozen data):

    python -m unittest discover -s tests -v

The fixtures in ``_support.py`` mirror the real frozen contract (18,164 genes,
50,520 train / 8,524 valid Tahoe rows, 649 matched-control groups, plate13/14
valid), so these tests also pin the contract numbers the wrapper relies on.
"""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(TESTS_DIR.parent))

import _support  # noqa: E402
import record  # noqa: E402


class RecordTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(prefix="pt05-test-")
        self.root = Path(self._tmp.name)
        self.data_dir = _support.make_data_dir(self.root)
        self.locked = _support.locked_from(self.data_dir)

    def tearDown(self) -> None:
        self._tmp.cleanup()

    # -- helpers -----------------------------------------------------------
    def valid_record(self) -> dict:
        return _support.make_record(self.root, data_dir=self.data_dir)

    def check(self, payload: dict, **kwargs) -> tuple[list[str], list[str]]:
        kwargs.setdefault("locked", self.locked)
        kwargs.setdefault("expect_mode", payload.get("mode"))
        kwargs.setdefault("verify_frozen", True)
        kwargs.setdefault("expect_complete", True)
        return record.check_record(payload, **kwargs)

    # -- happy paths -------------------------------------------------------
    def test_formal_record_passes(self) -> None:
        payload = self.valid_record()
        errors, warnings = self.check(payload)
        self.assertEqual(errors, [])
        self.assertEqual(warnings, [])
        self.assertTrue(payload["complete"])
        self.assertEqual(payload["seeds_planned"], _support.SEEDS)
        self.assertEqual(payload["seeds_completed"], _support.SEEDS)
        self.assertTrue(payload["no_seed_cherry_picking"])
        self.assertEqual(len(payload["units"]), 3)
        # the record carries the three distinct checkpoints and config digests
        self.assertEqual(len({unit["checkpoint"]["sha256"] for unit in payload["units"]}), 3)
        self.assertEqual(len({unit["config_sha256"] for unit in payload["units"]}), 3)
        for unit in payload["units"]:
            self.assertEqual(unit["status"], "completed")
            self.assertEqual(unit["epochs_run"], 100)
            self.assertEqual(unit["stop_reason"], "max_epochs")
            self.assertEqual(unit["budget_config"]["max_epochs"], 100)
            self.assertEqual(unit["budget_config"]["min_epochs"], 20)
            self.assertEqual(unit["budget_config"]["batch_size"], 512)
            self.assertEqual(unit["budget_config"]["lr"], 0.001)
            self.assertEqual(unit["timings"]["steps_per_epoch"], 99)
            self.assertEqual(unit["resources"]["host_peak_rss_mb_pipeline"], 6000.0)
            self.assertEqual(unit["resources"]["gpu_peak_mib_process"], 2900.0)
            self.assertEqual(len(unit["validation_curve"]), 100)

    def test_record_round_trip(self) -> None:
        payload = self.valid_record()
        paths = record.write_record(payload, self.root / "record")
        self.assertTrue(paths["json"].exists())
        self.assertTrue(paths["csv"].exists())
        reloaded = record.read_record(paths["json"])
        self.assertEqual(reloaded["schema"], record.SCHEMA)
        errors, _ = self.check(reloaded)
        self.assertEqual(errors, [])
        lines = paths["csv"].read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(lines), 1 + 3 * 100)
        self.assertTrue(lines[0].startswith("unit_id,seed,budget,epoch"))

    def test_profile_record_passes_but_never_formal(self) -> None:
        payload = _support.make_profile_record(self.root, data_dir=self.data_dir)
        errors, warnings = self.check(payload, expect_mode="profile")
        self.assertEqual(errors, [])
        self.assertTrue(any("PROFILE" in warning for warning in warnings))
        errors, _ = self.check(payload, expect_mode="formal")
        self.assertTrue(any("profile record" in error for error in errors))

    def test_real_frozen_locked_config_matches_ticket(self) -> None:
        if not (_support.REAL_DATA_DIR / "locked_config.json").exists():
            self.skipTest("real frozen 06.3 input is not checked out")
        locked = record.normalize_locked_config(
            record.read_json(_support.REAL_DATA_DIR / "locked_config.json")
        )
        self.assertEqual(record.ticket_expectation_errors(locked), [])

    # -- frozen contract tampering ----------------------------------------
    def test_ticket_expectations_reject_a_changed_recipe(self) -> None:
        locked = dict(_support.LOCKED_CONFIG)
        locked["recipe"] = dict(locked["recipe"], batch_size=256)
        errors = record.ticket_expectation_errors(record.normalize_locked_config(locked))
        self.assertTrue(any("batch_size" in error for error in errors))

    def test_input_artifact_tamper_fails(self) -> None:
        payload = self.valid_record()
        payload["input"]["artifacts"]["tahoe_transfer.h5ad"] = "0" * 64
        errors, _ = self.check(payload)
        self.assertTrue(any("cache" in error or "tahoe_transfer" in error for error in errors))

    def test_manifest_change_after_record_fails(self) -> None:
        payload = self.valid_record()
        manifests = _support.record.read_json(self.data_dir / "manifest.json")
        manifests["extra"] = "changed"
        _support.write_json(self.data_dir / "manifest.json", manifests)
        errors, _ = self.check(payload)
        self.assertTrue(any("manifest changed" in error for error in errors))

    def test_cache_from_another_source_fails(self) -> None:
        payload = self.valid_record()
        payload["input"]["cache"]["source_sha256"] = "1" * 64
        errors, _ = self.check(payload)
        self.assertTrue(any("derived cache" in error for error in errors))

    # -- budget / selection tampering -------------------------------------
    def test_shortened_formal_budget_fails(self) -> None:
        payload = self.valid_record()
        payload["budget"]["max_epochs"] = 5
        errors, _ = self.check(payload)
        self.assertTrue(any("max_epochs" in error for error in errors))

    def test_recipe_override_in_formal_record_fails(self) -> None:
        payload = self.valid_record()
        payload["budget"]["overrides"] = {"max_epochs": 5}
        errors, _ = self.check(payload)
        self.assertTrue(any("overrides" in error for error in errors))

    def test_dropped_epochs_in_formal_budget_fails(self) -> None:
        payload = self.valid_record()
        for unit in payload["units"]:
            unit["budget_config"]["max_epochs"] = 5
        errors, _ = self.check(payload)
        self.assertTrue(any("budget.max_epochs" in error for error in errors))

    def test_wrong_config_digest_fails(self) -> None:
        payload = self.valid_record()
        payload["units"][0]["config_sha256"] = "2" * 64
        errors, _ = self.check(payload)
        self.assertTrue(any("config digest" in error for error in errors))

    def test_early_stop_before_min_epochs_fails(self) -> None:
        payload = self.valid_record()
        unit = payload["units"][0]
        unit["epochs_run"] = 5
        unit["stop_reason"] = "early_stop_no_auprc_improvement_20"
        unit["validation_curve"] = unit["validation_curve"][:5]
        unit["best_epoch"] = 4
        unit["best_valid_auprc"] = unit["validation_curve"][3]["valid_auprc"]
        errors, _ = self.check(payload)
        self.assertTrue(any("before min_epochs" in error for error in errors))

    def test_epochs_beyond_max_epochs_fails(self) -> None:
        payload = self.valid_record()
        payload["units"][0]["epochs_run"] = 101
        errors, _ = self.check(payload)
        self.assertTrue(any("exceeds max_epochs" in error for error in errors))

    def test_best_epoch_not_curve_argmax_fails(self) -> None:
        payload = self.valid_record()
        payload["units"][0]["best_epoch"] = 3
        errors, _ = self.check(payload)
        self.assertTrue(any("argmax" in error for error in errors))

    def test_non_finite_curve_value_fails(self) -> None:
        payload = self.valid_record()
        payload["units"][0]["validation_curve"][10]["valid_auprc"] = float("nan")
        errors, _ = self.check(payload)
        self.assertTrue(any("not finite" in error for error in errors))

    def test_checkpoint_sha_mismatch_on_disk_fails(self) -> None:
        payload = self.valid_record()
        unit = payload["units"][0]
        Path(unit["checkpoint"]["path"]).write_bytes(b"tampered checkpoint")
        errors, _ = self.check(payload)
        self.assertTrue(any("checkpoint file sha256" in error for error in errors))

    def test_checkpoint_reload_flag_fails(self) -> None:
        payload = self.valid_record()
        payload["units"][0]["checkpoint"]["reload_identical"] = False
        errors, _ = self.check(payload)
        self.assertTrue(any("strict checkpoint reload" in error for error in errors))

    def test_partial_formal_sweep_is_incomplete(self) -> None:
        payload = self.valid_record()
        payload["units"] = payload["units"][:1]
        payload["complete"] = False
        payload["seeds_completed"] = [payload["units"][0]["seed"]]
        payload["seeds_missing"] = _support.SEEDS[1:]
        errors, _ = self.check(payload)
        self.assertTrue(any("incomplete" in error for error in errors))

    def test_duplicate_completed_seed_is_rejected(self) -> None:
        payload = self.valid_record()
        duplicate = _support.make_run_dir(
            self.root, seed=_support.SEEDS[0], unit_id="tahoe_seed_dup", config_sha256="e" * 64
        )
        entry = record.run_entry(
            run_dir=duplicate["run_dir"],
            unit_id=duplicate["unit_id"],
            seed=_support.SEEDS[0],
            budget="formal",
            command=["python"],
            exit_code=0,
            expected_config_sha256="e" * 64,
            resources={"host_peak_rss_mb_pipeline": None, "gpu_peak_mib_process": 2900.0},
        )
        with self.assertRaises(record.RecordError):
            record.build_record(
                mode="formal",
                units=[*payload["units"], entry],
                budget=payload["budget"],
                expectations=payload["expected"],
                input_identity=payload["input"],
                command=["python"],
                seeds_planned=_support.SEEDS,
            )

    # -- failed units ------------------------------------------------------
    def test_failed_unit_is_recorded_not_dropped(self) -> None:
        item = _support.make_run_dir(self.root, seed=_support.SEEDS[0])
        entry = record.run_entry(
            run_dir=item["run_dir"],
            unit_id=item["unit_id"],
            seed=item["seed"],
            budget="formal",
            command=["python", "train.py"],
            exit_code=1,
            expected_config_sha256=item["config_sha256"],
            log_tail="Traceback: synthetic failure",
        )
        self.assertEqual(entry["status"], "failed")
        self.assertIn("synthetic failure", entry["log_tail"])
        payload = self.valid_record()
        payload["units"][0] = entry
        payload["complete"] = False
        payload["seeds_failed"] = [entry["seed"]]
        payload["seeds_completed"] = _support.SEEDS[1:]
        errors, _ = self.check(payload)
        self.assertTrue(any("incomplete" in error or "failed" in error for error in errors))

    def test_missing_manifest_after_exit_zero_is_a_failure(self) -> None:
        run_dir = self.root / "runs" / "broken"
        run_dir.mkdir(parents=True)
        entry = record.run_entry(
            run_dir=run_dir,
            unit_id="broken",
            seed=_support.SEEDS[0],
            budget="formal",
            command=["python"],
            exit_code=0,
            expected_config_sha256="0" * 64,
            log_tail="no manifest written",
        )
        self.assertEqual(entry["status"], "failed")
        self.assertTrue(entry["failures"])

    def test_checker_rejects_a_profile_record_with_formal_budget(self) -> None:
        payload = self.valid_record()
        payload["mode"] = "profile"
        errors, _ = self.check(payload, expect_mode="profile")
        self.assertTrue(any("formal units" in error or "budget" in error for error in errors))


if __name__ == "__main__":
    unittest.main()
