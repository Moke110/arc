"""Verify the registered model, objective and measured-control prediction interfaces."""

from pathlib import Path
import os
import sys
import unittest

import numpy as np
import torch
from torch import nn

sys.path.insert(0, str((Path(__file__).resolve().parents[4] / "06_tsg_shape/06.2_round2/execution/tests" / Path(__file__).name).parents[1]))

from architecture_model import ArchitectureAdapter, ShapeAdapter, axis_run, configure_objective

ROOT = Path(os.environ.get("ARC_ROOT", (Path(__file__).resolve().parents[4] / "06_tsg_shape/06.2_round2/execution/tests" / Path(__file__).name).parents[5]))
SOURCE_ROOT = ROOT / "pert_model/TranSiGen/source"
TRANSFORM = dict(kind="training-log10-dose-population-standardization-v1", mean=0., std=1.)


class ArchitectureTests(unittest.TestCase):
    def build(self, shape, adapter=ArchitectureAdapter):
        return adapter(SOURCE_ROOT, TRANSFORM, shape).build(42, "cpu", n_genes=4)

    def test_original_preserves_seeded_step_061_predictions(self):
        model = self.build("original").eval()
        accepted = self.build("original", ShapeAdapter).eval()
        self.assertEqual(model.state_dict().keys(), accepted.state_dict().keys())
        self.assertTrue(all(torch.equal(value, accepted.state_dict()[key])
                            for key, value in model.state_dict().items()))
        control = torch.ones(2, 4)
        features = torch.ones(2, 2305)
        with torch.no_grad():
            torch.manual_seed(703)
            actual = model(control, features)
            torch.manual_seed(703)
            expected = accepted(control, features)
        self.assertTrue(all(torch.equal(left, right)
                            for left, right in zip(actual, expected, strict=True)))


    def test_no_x2_training_and_validation_only_decode_the_predicted_latent(self):
        native = ArchitectureAdapter(SOURCE_ROOT, TRANSFORM, "output-delta-nox2")
        model = native.build(42, "cpu", n_genes=4)
        configure_objective(model, "output-delta-nox2")
        rows = np.arange(2)

        class Inputs:
            control = np.ones((2, 4), dtype=np.float32)
            expression = np.arange(8, dtype=np.float32).reshape(2, 4)

            def doses(self, batch):
                return np.ones(len(batch))

            def fingerprints(self, batch):
                return np.zeros((len(batch), 2304), dtype=np.float32)

            def model_input(self, phase):
                return self

        calls = []
        with model.decoder_x2.register_forward_hook(lambda *_: calls.append("decode")):
            native.train_epoch(model, Inputs(), rows, native.optimizer(model), np.random.default_rng(42))
            values, objective, _ = axis_run.validate(
                model, Inputs(), rows,
                {"model": {"transform": TRANSFORM, "batch_size": 64}, "validation_seed": 22617},
            )
        self.assertEqual(calls, ["decode", "decode"])
        self.assertEqual(values.shape, (2, 4))
        self.assertGreater(objective, 0)
        self.assertGreater(model.mu_z2[0].weight.grad.abs().sum().item(), 0)

    def test_no_x2_objective_counts_delta_mse_twice_and_all_three_kl_terms(self):
        model = self.build("output-delta-nox2")
        configure_objective(model, "output-delta-nox2")
        control, control_rec, measured, predicted = [torch.Tensor([[value]]) for value in (4, 5, 6, 3)]
        mu1, mu2, mu_pred, logvar = [torch.Tensor([[value]]) for value in (1, 2, 3, 0)]
        loss, control_mse, treatment_mse, delta_mse, k1, k2, kp, delta_l1 = model.loss(
            control, control_rec, mu1, logvar, measured, None, mu2, logvar, predicted, mu_pred, logvar,
        )
        self.assertEqual([part.item() for part in (control_mse, treatment_mse, delta_mse, k1, k2, kp, delta_l1)],
                         [1., 0., 9., .5, 2., .5, 3.])
        self.assertAlmostEqual(loss.item(), 24.55, places=5)

    def test_direct_delta_predictions_keep_both_signs_and_the_same_control(self):
        for shape, decoder in (("output-delta", "decoder_delta"), ("output-delta-nox2", "decoder_x2")):
            with self.subTest(shape=shape), torch.no_grad():
                model = self.build(shape).eval()
                head = getattr(model, decoder)[-1]
                head.weight.zero_()
                head.bias.copy_(torch.tensor([-1., 1., -2., 2.]))
                control = torch.arange(8, dtype=torch.float32).view(2, 4)
                outputs = model(control, torch.ones(2, 2305))
                expected = torch.tensor([[-1., 1., -2., 2.]]).expand_as(control)
                self.assertTrue(torch.equal(outputs[3] - control, expected))
                measured_latent, _, _ = model.encode_x2(control)
                reconstruction = model.decode_x2(measured_latent)
                if shape == "output-delta":
                    self.assertEqual(reconstruction.shape, control.shape)
                    self.assertTrue((reconstruction >= 0).all())
                else:
                    self.assertIsNone(reconstruction)

    def test_interaction_heads_return_comparable_gaussian_parameters(self):
        for shape, parameters in (("z-interact_MLP", 106485), ("z-interact_FiLM", 106488)):
            model = self.build(shape).eval()
            self.assertEqual(sum(p.numel() for p in model.prediction_head.parameters()), parameters)
            control, drug_dose = torch.randn(2, 100), torch.randn(2, 432)
            mu, logvar = model.predictive_distribution(control, drug_dose)
            self.assertEqual(mu.shape, (2, 100))
            self.assertEqual(logvar.shape, (2, 100))
            if shape == "z-interact_FiLM":
                projected = model.prediction_head.projection(control)
                self.assertTrue(torch.equal(torch.cat((mu, logvar), dim=1), projected))
        original = self.build("original")
        heads = (original.mu_z2Fz1, original.logvar_z2Fz1)
        self.assertEqual(sum(p.numel() for head in heads for p in head.parameters()), 106600)

    def test_tied_encoders_share_batchnorm_and_residual_changes_only_the_mean(self):
        control, drug_dose = torch.randn(2, 100), torch.randn(2, 432)
        combined = torch.cat((control, drug_dose), dim=1)
        for shape in ("shared-z", "z-res", "shared-z-res"):
            with self.subTest(shape=shape):
                model = self.build(shape)
                tied = shape != "z-res"
                self.assertEqual(model.encoder_x1 is model.encoder_x2, tied)
                self.assertIsNot(model.mu_z1, model.mu_z2)
                control_bn = next(layer for layer in model.encoder_x1 if isinstance(layer, nn.BatchNorm1d))
                model.encode_x1(torch.zeros(2, 4))
                model.encode_x2(torch.ones(2, 4))
                self.assertEqual(control_bn.num_batches_tracked.item(), 2 if tied else 1)
                mu, logvar = model.predictive_distribution(control, drug_dose)
                expected = model.mu_z2Fz1(combined)
                if shape != "shared-z":
                    expected = expected + control
                self.assertTrue(torch.equal(mu, expected))
                self.assertTrue(torch.equal(logvar, model.logvar_z2Fz1(combined)))


if __name__ == "__main__":
    unittest.main()
