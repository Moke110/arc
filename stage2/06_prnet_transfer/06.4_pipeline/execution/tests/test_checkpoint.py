"""Checkpoint payload, strict reload and the warm-start axis/structure guards."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import torch

import _support  # noqa: F401  (path setup)
from ptlib import checkpoint, config, model


def structure(n_genes: int = 64, hidden=(16,), z_dim: int = 8, c_dim: int = 8,
              n_dim: int = 4, drug_dim: int = 32, adaptor=(8,)) -> dict:
    return {
        "n_genes": n_genes,
        "hidden": list(hidden),
        "z_dim": z_dim,
        "c_dim": c_dim,
        "n_dim": n_dim,
        "drug_dim": drug_dim,
        "adaptor": list(adaptor),
        "dropout": 0.05,
    }


class CheckpointRoundTripTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        torch.manual_seed(0)
        self.struct = structure()
        self.model = model.build_model(self.struct)
        model.random_init(self.model, seed=11)
        self.payload_config = {
            "structure": self.struct, "recipe": {"lr": 1e-3}, "metric": {"k": 10}
        }
        self.path = self.dir / "checkpoint/model.pt"
        checkpoint.save_checkpoint(
            self.path, self.model, self.payload_config, seed=11, n_genes=64,
            best_valid_auprc=0.42, epoch=3,
        )

    def tearDown(self):
        self.tmp.cleanup()

    def test_payload_keys_are_exactly_the_contract(self):
        payload = checkpoint.read_checkpoint(self.path)
        self.assertEqual(set(payload), set(checkpoint.PAYLOAD_KEYS))
        self.assertEqual(payload["model"], "prnet")
        self.assertEqual(payload["n_genes"], 64)
        self.assertEqual(payload["seed"], 11)
        self.assertEqual(payload["epoch"], 3)
        self.assertAlmostEqual(payload["best_valid_auprc"], 0.42)

    def test_strict_reload_reproduces_parameters(self):
        payload = checkpoint.read_checkpoint(self.path)
        checkpoint.validate_checkpoint(payload, 64, self.struct)
        fresh = model.build_model(self.struct)
        checkpoint.load_into_model(fresh, payload)
        for key, value in self.model.state_dict().items():
            torch.testing.assert_close(fresh.state_dict()[key], value)

    def test_build_from_checkpoint_matches(self):
        payload = checkpoint.read_checkpoint(self.path)
        fresh, struct = checkpoint.build_from_checkpoint(payload, 64, self.struct)
        self.assertEqual(struct["z_dim"], 8)
        for key, value in self.model.state_dict().items():
            torch.testing.assert_close(fresh.state_dict()[key], value)

    def test_checkpoint_sha256(self):
        digest = checkpoint.checkpoint_sha256(self.path)
        self.assertEqual(len(digest), 64)
        self.assertEqual(digest, checkpoint.checkpoint_sha256(self.path))


class StrictnessTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.struct = structure()
        self.model = model.build_model(self.struct)
        self.path = self.dir / "model.pt"
        checkpoint.save_checkpoint(
            self.path, self.model, {"structure": self.struct}, seed=1, n_genes=64,
            best_valid_auprc=0.1, epoch=1,
        )

    def tearDown(self):
        self.tmp.cleanup()

    def _payload(self):
        return checkpoint.read_checkpoint(self.path)

    def test_missing_key_rejected(self):
        payload = self._payload()
        payload["state_dict"].pop(next(iter(payload["state_dict"])))
        with self.assertRaises(checkpoint.CheckpointError):
            checkpoint.load_into_model(model.build_model(self.struct), payload)

    def test_renamed_key_rejected(self):
        payload = self._payload()
        key = "encoder.mean_encoder.weight"
        payload["state_dict"]["encoder.mean_encoder.MEAN"] = payload["state_dict"].pop(key)
        with self.assertRaises(checkpoint.CheckpointError):
            checkpoint.load_into_model(model.build_model(self.struct), payload)

    def test_non_pipeline_checkpoint_rejected(self):
        torch.save({"state_dict": {}}, self.path)
        with self.assertRaises(checkpoint.CheckpointError):
            checkpoint.read_checkpoint(self.path)

    def test_axis_mismatch_fails_loudly(self):
        payload = self._payload()
        with self.assertRaisesRegex(checkpoint.CheckpointError, "gene-axis"):
            checkpoint.validate_checkpoint(payload, 128, self.struct)

    def test_structure_mismatch_fails_loudly(self):
        payload = self._payload()
        other = structure(z_dim=16)
        with self.assertRaisesRegex(checkpoint.CheckpointError, "incompatible"):
            checkpoint.validate_checkpoint(payload, 64, other)

    def test_missing_checkpoint_fails(self):
        with self.assertRaises(checkpoint.CheckpointError):
            checkpoint.read_checkpoint(self.dir / "nope.pt")


class NormalizeStructureTests(unittest.TestCase):
    def test_alias_structure_from_checkpoint_config(self):
        cfg = config.normalize_structure(
            {"in_dim": 100, "hidden": 32, "state_latent": 8, "drug_latent": 8,
             "noise": 4, "n_bits": 16, "adaptor_layer_sizes": [8]}
        )
        self.assertEqual(cfg["n_genes"], 100)
        self.assertEqual(cfg["hidden"], [32])
        self.assertEqual(cfg["drug_dim"], 16)


if __name__ == "__main__":
    unittest.main()
