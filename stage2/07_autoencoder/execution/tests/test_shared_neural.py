"""The shared matrix models expose the same reconstruction boundary."""

import json
from math import isfinite
from pathlib import Path
import sys
import unittest

import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "script"))

from benchmark import SHARED_NEURAL_IDS
from shared_neural import build_matrix_model
import test_benchmark


class MatrixModelTest(unittest.TestCase):
    def test_extra_queries_are_matched_and_affect_reconstruction(self):
        torch.manual_seed(17)
        paired = build_matrix_model("AAEcs_1024", 7, 16, fixture=True)
        torch.manual_seed(17)
        extended = build_matrix_model("AAEcs_1024+64", 7, 16, fixture=True)
        for key, value in paired.state_dict().items():
            other = extended.state_dict()[key]
            if key == "latent_queries":
                other = other[:, :3]
            torch.testing.assert_close(value, other, rtol=0, atol=0)

        matrix = torch.randn(2, 7, 16)
        positions = torch.randn(7, 16)
        latent = extended.encode(matrix)
        self.assertEqual(tuple(latent.shape), (2, 3, 16))
        restored = extended.decode(latent, positions)
        self.assertEqual(restored.shape, matrix.shape)
        self.assertTrue(torch.isfinite(restored).all())
        restored.square().mean().backward()
        self.assertGreater(extended.latent_queries.grad[:, 3:].abs().sum().item(), 0)
        with torch.no_grad():
            extended.latent_queries[:, 3:].zero_()
            modified = extended.decode(extended.encode(matrix), positions)
        self.assertFalse(torch.allclose(restored, modified))

    def test_partial_epoch_checkpoints_preserve_the_joint_order(self):
        from neural import training_windows

        self.assertEqual(list(training_windows(4, 7, 3)), [
            (1, 0, 3), (1, 3, 1), (2, 0, 2), (2, 2, 1),
        ])

    def test_shared_models_fit_validate_and_reload_the_common_budget(self):
        fixture = test_benchmark.BenchmarkTest()
        self.addCleanup(fixture.doCleanups)
        fixture.setUp()
        fixture.cli("freeze-split")
        fixture.cli("fit", "--output", fixture.artifact)
        interface = fixture.root / "interface"
        interface_record = fixture.cli("fit", "--model", "shared_expression", "--output", interface)
        orders = []
        for model_id in SHARED_NEURAL_IDS:
            with self.subTest(model=model_id):
                artifact = fixture.root / model_id
                record = fixture.cli("fit", "--model", model_id, "--interface", interface,
                                     "--output", artifact)
                self.assertEqual(record["metadata"]["latent_shape"], [3, 16])
                self.assertEqual(record["metadata"]["interface_manifest_sha256"],
                                 interface_record["artifact_manifest_sha256"])
                history = json.loads((artifact / "epochs.json").read_text())
                self.assertGreaterEqual(len(history), 3)
                self.assertEqual(history[-1]["optimizer_steps"], 5)
                orders.append([entry["sample_order_sha256"] for entry in history])
                score = fixture.root / f"{model_id}-validation.json"
                options = ("--artifact", artifact, "--mean-artifact", fixture.artifact,
                           "--interface", interface)
                validated = fixture.cli("score", *options, "--output", score)
                self.assertEqual(validated, fixture.cli("reload", *options, "--expected", score))
                for metrics in validated["sources"].values():
                    self.assertTrue(isfinite(metrics["mse"]))
        for order in orders[1:]:
            self.assertEqual(order, orders[0])
