"""Unit-manifest tests: identity chain, fail-closed rules, the 30-unit shape."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd  # noqa: E402

import _support  # noqa: E402
import common  # noqa: E402
import units as units_mod  # noqa: E402


class UnitManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.root = _support.clean_dir(_support.tmp_root() / "units")
        cls.inputs = _support.make_selection_inputs(cls.root)

    def _read_record(self):
        return common.read_table(
            self.inputs["random_record"],
            numeric_columns=("fold", "seed", "seed_index", "checkpoint_bytes",
                             "lr_selected", "best_valid_auprc", "best_epoch", "epochs_run"),
        )

    def _build(self, **overrides):
        inputs = self.inputs
        args = dict(
            data_dir=inputs["data"]["dir"],
            random_record=inputs["random_record"],
            tahoe_record=inputs["tahoe_record"],
            baselines_dir=None,
        )
        args.update(overrides)
        return units_mod.build_unit_manifest(**args)

    def test_complete_manifest_has_30_verified_units(self):
        payload = self._build()
        self.assertTrue(payload["complete"])
        self.assertEqual(len(payload["units"]), 30)
        self.assertEqual(sum(1 for u in payload["units"] if u["init_name"] == "random"), 15)
        self.assertEqual(sum(1 for u in payload["units"] if u["init_name"] == "tahoe"), 15)
        folds = {int(u["fold"]) for u in payload["units"]}
        self.assertEqual(folds, {0, 1, 2, 3, 4})
        for unit in payload["units"]:
            self.assertEqual(unit["data_identity"], payload["data_identity"])
            self.assertTrue(unit["checkpoint"]["sha256_verified"])
            self.assertEqual(len(unit["checkpoint"]["sha256"]), 64)
            self.assertEqual(len(unit["config_sha256"]), 64)
        gates = payload["gates"]
        for key, value in gates.items():
            if isinstance(value, bool):
                self.assertTrue(value, f"gate {key} failed: {gates.get(key)}")

    def test_unit_ids_are_unique_and_carry_selection(self):
        payload = self._build()
        ids = [u["unit_id"] for u in payload["units"]]
        self.assertEqual(len(ids), len(set(ids)))
        for unit in payload["units"]:
            self.assertIn(f"f{unit['fold']}", unit["unit_id"])
            self.assertIn(f"s{unit['seed']}", unit["unit_id"])
            self.assertTrue(unit["unit_id"].endswith(unit["init_name"]))

    def test_tahoe_warm_start_provenance_recorded(self):
        payload = self._build()
        for unit in payload["units"]:
            if unit["init_name"] != "tahoe":
                continue
            init = unit["init"]
            self.assertEqual(init["mode"], "pretrained")
            self.assertTrue(init["source_sha256_verified"])
            self.assertTrue(init["source_file_present"])

    def test_missing_tahoe_record_fails_closed(self):
        with self.assertRaises(units_mod.UnitError):
            self._build(tahoe_record=self.root / "missing.csv", allow_incomplete=False)

    def test_allow_incomplete_builds_random_branch_only(self):
        payload = self._build(tahoe_record=self.root / "missing.csv", allow_incomplete=True)
        self.assertFalse(payload["complete"])
        self.assertEqual(len(payload["units"]), 15)
        self.assertEqual(payload["missing_inits"], ["tahoe"])

    def test_tampered_checkpoint_sha_fails(self):
        record = self._read_record()
        record.loc[0, "checkpoint_sha256"] = "0" * 64
        tampered = self.root / "tampered_sha.csv"
        record.to_csv(tampered, index=False)
        with self.assertRaises(units_mod.UnitError) as caught:
            self._build(random_record=tampered)
        self.assertIn("sha256 mismatch", str(caught.exception))

    def test_tampered_checkpoint_bytes_fails(self):
        record = self._read_record()
        record.loc[0, "checkpoint_bytes"] = int(record.loc[0, "checkpoint_bytes"]) + 1
        tampered = self.root / "tampered_bytes.csv"
        record.to_csv(tampered, index=False)
        with self.assertRaises(units_mod.UnitError) as caught:
            self._build(random_record=tampered)
        self.assertIn("size mismatch", str(caught.exception))

    def test_tampered_config_digest_fails(self):
        record = self._read_record()
        record.loc[0, "config_sha256"] = "f" * 64
        tampered = self.root / "tampered_config.csv"
        record.to_csv(tampered, index=False)
        with self.assertRaises(units_mod.UnitError) as caught:
            self._build(random_record=tampered)
        self.assertIn("config_sha256", str(caught.exception))

    def test_wrong_data_identity_fails(self):
        record = self._read_record()
        record["data_identity"] = "e" * 64
        tampered = self.root / "tampered_identity.csv"
        record.to_csv(tampered, index=False)
        with self.assertRaises(units_mod.UnitError) as caught:
            self._build(random_record=tampered)
        self.assertIn("data_identity", str(caught.exception))

    def test_duplicate_fold_seed_fails(self):
        record = self._read_record()
        duplicate = pd.concat([record, record.iloc[[0]]], ignore_index=True)
        tampered = self.root / "duplicate.csv"
        duplicate.to_csv(tampered, index=False)
        with self.assertRaises(units_mod.UnitError) as caught:
            self._build(random_record=tampered)
        self.assertIn("duplicate", str(caught.exception))

    def test_two_learning_rates_in_one_fold_fails(self):
        record = self._read_record()
        second = record.index[record["fold"] == 0][1]
        record.loc[second, "lr_selected"] = 3e-4
        tampered = self.root / "two_lrs.csv"
        record.to_csv(tampered, index=False)
        with self.assertRaises(units_mod.UnitError) as caught:
            self._build(random_record=tampered)
        self.assertIn("more than one selected LR", str(caught.exception))

    def test_missing_seed_fails(self):
        record = self._read_record()
        record = record[~((record["fold"] == 2) & (record["seed"] == _support.SEEDS[2]))]
        tampered = self.root / "missing_seed.csv"
        record.to_csv(tampered, index=False)
        with self.assertRaises(units_mod.UnitError):
            self._build(random_record=tampered)

    def test_prediction_artifact_in_candidate_run_fails(self):
        run = self.inputs["random_runs"][0]["result"]["run_dir"]
        artifact = run / "prediction"
        artifact.mkdir()
        try:
            with self.assertRaises(units_mod.UnitError) as caught:
                self._build()
            self.assertIn("prediction/evaluation artifacts", str(caught.exception))
        finally:
            artifact.rmdir()

    def test_random_unit_requires_random_manifest(self):
        import json

        run = self.inputs["random_runs"][0]["result"]["run_dir"]
        manifest_file = run / "manifest.json"
        original = json.loads(manifest_file.read_text(encoding="utf-8"))
        tampered = json.loads(manifest_file.read_text(encoding="utf-8"))
        tampered["init"]["mode"] = "pretrained"
        manifest_file.write_text(json.dumps(tampered), encoding="utf-8")
        try:
            with self.assertRaises(units_mod.UnitError):
                self._build()
        finally:
            manifest_file.write_text(json.dumps(original), encoding="utf-8")

    def test_tahoe_unit_requires_pretrained_manifest(self):
        import json

        run = self.inputs["tahoe_runs"][0]["result"]["run_dir"]
        manifest_file = run / "manifest.json"
        original = json.loads(manifest_file.read_text(encoding="utf-8"))
        tampered = json.loads(manifest_file.read_text(encoding="utf-8"))
        tampered["init"]["mode"] = "random"
        manifest_file.write_text(json.dumps(tampered), encoding="utf-8")
        try:
            with self.assertRaises(units_mod.UnitError):
                self._build()
        finally:
            manifest_file.write_text(json.dumps(original), encoding="utf-8")

    def test_frozen_identity_matches_selection_record(self):
        data = self.inputs["data"]
        identity = common.data_identity_sha256(data["frozen"])
        self.assertEqual(identity, data["data_identity"])


class ExpectedFoldCountsTests(unittest.TestCase):
    def test_reads_baseline_coverage(self):
        import pandas as pd

        root = _support.clean_dir(_support.tmp_root() / "foldref")
        baselines = root / "baselines"
        baselines.mkdir()
        pd.DataFrame(
            [
                {"fold": 0, "n_valid_rows": 188, "n_test_rows": 192, "n_test_conditions": 68, "n_test_drugs": 23},
                {"fold": 1, "n_valid_rows": 191, "n_test_rows": 188, "n_test_conditions": 67, "n_test_drugs": 23},
            ]
        ).to_csv(baselines / "fold_coverage.csv", index=False)
        reference = units_mod.expected_fold_counts(baselines)
        self.assertEqual(reference["0"]["n_test_rows"], 192)
        self.assertEqual(reference["1"]["n_valid_rows"], 191)
        self.assertIn("source", reference)


if __name__ == "__main__":
    unittest.main()
