"""Unit tests for the ticket-#7 strict warm-start loader and its provenance record."""

from __future__ import annotations

import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import _support  # noqa: E402
import warmstart  # noqa: E402

from ptlib import checkpoint as pt_checkpoint  # noqa: E402
from ptlib import model as pt_model  # noqa: E402

import torch  # noqa: E402

SEED = 20261021
STRUCT = _support.SMALL_STRUCTURE


class WarmStartTestBase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def _checkpoint(self, *, seed: int = SEED, perturb: bool = True, stage: str = "pretrain") -> Path:
        model = pt_model.build_model(STRUCT)
        pt_model.random_init(model, seed)
        if perturb:
            _support._perturb(model, seed)
        config = {
            "structure": dict(STRUCT),
            "recipe": {"lr": 1e-3},
            "metric": {"k": 200},
            "stage": stage,
            "seed": int(seed),
            "init": "random",
        }
        path = warmstart.checkpoint_path(self.tmp / "runs", seed)
        pt_checkpoint.save_checkpoint(
            path, model, config=config, seed=int(seed), n_genes=int(STRUCT["n_genes"]),
            best_valid_auprc=0.13, epoch=2,
        )
        return path


class VerifyCheckpointTests(WarmStartTestBase):
    def test_provenance_record_is_complete(self):
        path = self._checkpoint()
        record = warmstart.verify_checkpoint(SEED, path, STRUCT)
        self.assertEqual(record["seed"], SEED)
        self.assertEqual(record["sha256"], _support.common.sha256_file(path))
        self.assertEqual(record["payload_seed"], SEED)
        self.assertTrue(record["seed_matches"])
        self.assertEqual(record["checkpoint_stage"], "pretrain")
        self.assertEqual(record["gene_axis"],
                         {"checkpoint": 256, "expected": 256, "matches": True})
        self.assertTrue(record["structure_matches_run"])
        self.assertEqual(record["parameter_count"], pt_model.parameter_count(
            pt_model.build_model(STRUCT)
        ))
        provenance = record["load_provenance"]
        self.assertTrue(provenance["strict_load"])
        self.assertTrue(provenance["weights_identical"])
        self.assertEqual(provenance["missing_keys"], [])
        self.assertEqual(provenance["missing_parameter_keys"], [])
        self.assertEqual(provenance["unexpected_keys"], [])
        self.assertEqual(provenance["shape_mismatches"], [])
        self.assertEqual(provenance["dtype_mismatches"], [])
        self.assertEqual(provenance["loaded_layers"], sorted(
            pt_model.build_model(STRUCT).state_dict()
        ))
        self.assertEqual(provenance["n_loaded_tensors"], len(provenance["loaded_layers"]))
        self.assertGreater(provenance["n_loaded_parameter_tensors"], 0)
        self.assertEqual(len(record["checkpoint_config_sha256"]), 64)
        anti = record["anti_random_init"]
        self.assertFalse(anti["matches_fresh_random_init"])
        self.assertGreater(anti["n_tensors_identical_to_random_init"], 0)  # buffers start equal
        self.assertGreater(
            anti["n_tensors"] - anti["n_tensors_identical_to_random_init"], 0
        )

    def test_missing_checkpoint_fails(self):
        with self.assertRaisesRegex(warmstart.WarmStartError, "missing"):
            warmstart.verify_checkpoint(SEED, self.tmp / "nope.pt", STRUCT)

    def test_gene_axis_mismatch_fails_closed(self):
        path = self._checkpoint()
        other = {**STRUCT, "n_genes": 512}
        with self.assertRaisesRegex(warmstart.WarmStartError, "gene-axis"):
            warmstart.verify_checkpoint(SEED, path, other)

    def test_structure_mismatch_fails_closed(self):
        path = self._checkpoint()
        other = {**STRUCT, "z_dim": 16}
        with self.assertRaisesRegex(warmstart.WarmStartError, "incompatible"):
            warmstart.verify_checkpoint(SEED, path, other)

    def test_wrong_seed_in_payload_fails_closed(self):
        path = self._checkpoint(seed=SEED + 1)
        with self.assertRaisesRegex(warmstart.WarmStartError, "payload seed"):
            warmstart.verify_checkpoint(SEED, path, STRUCT)

    def test_non_pretrain_stage_fails_closed(self):
        path = self._checkpoint(stage="finetune")
        with self.assertRaisesRegex(warmstart.WarmStartError, "stage 'pretrain'"):
            warmstart.verify_checkpoint(SEED, path, STRUCT)

    def test_untrained_checkpoint_fails_the_random_init_canary(self):
        path = self._checkpoint(perturb=False)
        with self.assertRaisesRegex(warmstart.WarmStartError, "not a trained pretrain"):
            warmstart.verify_checkpoint(SEED, path, STRUCT)

    def test_renamed_key_fails_closed(self):
        path = self._checkpoint()
        payload = torch.load(path, map_location="cpu", weights_only=False)
        old = "encoder.mean_encoder.weight"
        payload["state_dict"][f"renamed.{old}"] = payload["state_dict"].pop(old)
        torch.save(payload, path)
        with self.assertRaises(warmstart.WarmStartError):
            warmstart.verify_checkpoint(SEED, path, STRUCT)

    def test_partial_state_dict_fails_closed(self):
        path = self._checkpoint()
        payload = torch.load(path, map_location="cpu", weights_only=False)
        first = next(iter(payload["state_dict"]))
        payload["state_dict"] = {first: payload["state_dict"][first]}
        torch.save(payload, path)
        with self.assertRaisesRegex(warmstart.WarmStartError, "missing keys"):
            warmstart.verify_checkpoint(SEED, path, STRUCT)


class VerifySeedsTests(WarmStartTestBase):
    def test_three_seeds_verified_with_one_digest(self):
        root = _support.make_checkpoints(self.tmp)
        first = warmstart.verify_seeds(
            (20261021, 20261022, 20261023), root, STRUCT
        )
        second = warmstart.verify_seeds((20261021, 20261022, 20261023), root, STRUCT)
        self.assertTrue(first["all_verified"])
        self.assertEqual(sorted(first["checkpoints"]), ["20261021", "20261022", "20261023"])
        self.assertEqual(first["provenance_sha256"], second["provenance_sha256"])
        shas = {record["sha256"] for record in first["checkpoints"].values()}
        self.assertEqual(len(shas), 3)

    def test_any_bad_seed_fails_the_whole_verification(self):
        root = _support.make_checkpoints(self.tmp)
        victim = warmstart.checkpoint_path(root, 20261022)
        payload = torch.load(victim, map_location="cpu", weights_only=False)
        payload["n_genes"] = 1
        torch.save(payload, victim)
        with self.assertRaises(warmstart.WarmStartError):
            warmstart.verify_seeds((20261021, 20261022, 20261023), root, STRUCT)

    def test_missing_seed_checkpoint_fails(self):
        root = _support.make_checkpoints(self.tmp, seeds=(20261021,))
        with self.assertRaises(warmstart.WarmStartError):
            warmstart.verify_seeds((20261021, 20261022, 20261023), root, STRUCT)


class VerifyManifestInitTests(WarmStartTestBase):
    def setUp(self):
        super().setUp()
        self.path = self._checkpoint()
        self.record = warmstart.verify_checkpoint(SEED, self.path, STRUCT)
        self.init = _support.init_record(self.record)

    def test_complete_record_passes(self):
        self.assertEqual(warmstart.verify_manifest_init(self.init, self.record), [])

    def test_random_init_run_is_rejected(self):
        broken = {**self.init, "mode": "random"}
        problems = warmstart.verify_manifest_init(broken, self.record)
        self.assertTrue(any("pretrained" in problem for problem in problems))

    def test_wrong_source_path_is_rejected(self):
        broken = {**self.init, "source_path": str(self.tmp / "other.pt")}
        problems = warmstart.verify_manifest_init(broken, self.record)
        self.assertTrue(any("source_path" in problem for problem in problems))

    def test_checkpoint_sha_drift_is_rejected(self):
        broken = {**self.init, "checkpoint_sha256": "0" * 64}
        problems = warmstart.verify_manifest_init(broken, self.record)
        self.assertTrue(any("sha256" in problem for problem in problems))

    def test_config_digest_drift_is_rejected(self):
        broken = {**self.init, "checkpoint_config_sha256": "0" * 64}
        problems = warmstart.verify_manifest_init(broken, self.record)
        self.assertTrue(any("config digest" in problem for problem in problems))

    def test_gene_axis_drift_is_rejected(self):
        broken = {**self.init, "gene_axis": {"checkpoint": 18164, "expected": 18164, "matches": True}}
        problems = warmstart.verify_manifest_init(broken, self.record)
        self.assertTrue(any("gene_axis" in problem for problem in problems))

    def test_partial_load_provenance_is_rejected(self):
        provenance = dict(self.init["load_provenance"])
        provenance["missing_keys"] = ["encoder.mean_encoder.weight"]
        broken = {**self.init, "load_provenance": provenance}
        problems = warmstart.verify_manifest_init(broken, self.record)
        self.assertTrue(any("missing_keys" in problem for problem in problems))

    def test_missing_provenance_is_rejected(self):
        broken = {key: value for key, value in self.init.items() if key != "load_provenance"}
        problems = warmstart.verify_manifest_init(broken, self.record)
        self.assertTrue(any("load_provenance is missing" in problem for problem in problems))

    def test_seed_drift_is_rejected(self):
        broken = {**self.init, "checkpoint_seed": SEED + 1}
        problems = warmstart.verify_manifest_init(broken, self.record)
        self.assertTrue(any("checkpoint_seed" in problem for problem in problems))


if __name__ == "__main__":
    unittest.main()
