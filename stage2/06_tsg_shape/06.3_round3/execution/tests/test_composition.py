"""Check the composed network and objective against the registered components."""

from pathlib import Path
import sys
import unittest

import torch
from torch import nn
from torch.nn import functional as F

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from composition_model import CompositionAdapter, axis_run, configure_objective
from composition_roster import coordinates, fit_coordinates

SOURCE = axis_run.REPO / "pert_model/TranSiGen/source"
TRANSFORM = dict(kind="training-log10-dose-population-standardization-v1", mean=0., std=1.)


class CompositionTests(unittest.TestCase):
    def build(self, shape, seed=17):
        return CompositionAdapter(SOURCE, TRANSFORM, shape).build(seed, "cpu", n_genes=4)

    def test_registered_matrix(self):
        self.assertEqual(len(coordinates(False)), 72)
        self.assertEqual(len(fit_coordinates(False)), 36)
        self.assertEqual(len(coordinates(True)), 24)
        self.assertEqual({key[0] for key in fit_coordinates(False)},
                         {"3_z-interact_FiLM", "4_1+3", "5_2+3"})

    def test_film_arms_normalize_after_modulation_and_before_projection(self):
        for shape in ("3_z-interact_FiLM", "4_1+3", "5_2+3"):
            with self.subTest(shape=shape), torch.no_grad():
                model = self.build(shape).eval()
                head = model.prediction_head
                self.assertIsInstance(head.normalization, nn.LayerNorm)
                self.assertEqual(head.normalization.normalized_shape, (100,))
                self.assertFalse(head.normalization.elementwise_affine)
                self.assertEqual(list(head.normalization.parameters()), [])

                latent, features = torch.randn(2, 100), torch.randn(2, 432)
                self.assertFalse(head.conditioner(features).any())
                initial = torch.cat(model.predictive_distribution(latent, features), dim=1)
                torch.testing.assert_close(initial, head.projection(F.layer_norm(latent, (100,))))

                head.conditioner[-1].bias.copy_(torch.linspace(-.75, .75, 200))
                scale, shift = head.conditioner(features).chunk(2, dim=1)
                normalized = F.layer_norm((1 + scale) * latent + shift, (100,))
                actual = torch.cat(model.predictive_distribution(latent, features), dim=1)
                torch.testing.assert_close(actual, head.projection(normalized))

    def test_combinations_predict_signed_delta(self):
        for shape, decoder in (("4_1+3", "decoder_delta"), ("5_2+3", "decoder_x2")):
            with self.subTest(shape=shape), torch.no_grad():
                model = self.build(shape).eval()
                head = getattr(model, decoder)[-1]
                head.weight.zero_()
                signed_delta = torch.tensor([-1., 1., -2., 2.])
                head.bias.copy_(signed_delta)
                control = torch.arange(8, dtype=torch.float32).reshape(2, 4)
                predicted = model(control, torch.ones(2, 2305))[3]
                torch.testing.assert_close(predicted, control + signed_delta)

    def test_no_x2_paths_keep_the_measured_encoder_and_use_all_kl_terms(self):
        for shape in ("2_output-delta-no-x2", "5_2+3"):
            with self.subTest(shape=shape):
                model = self.build(shape).eval()
                self.assertIsNot(model.encoder_x1, model.encoder_x2)
                latent, measured_mean, _ = model.encode_x2(torch.ones(2, 4))
                self.assertIsNone(latent)
                measured_mean.sum().backward()
                self.assertGreater(model.mu_z2[0].weight.grad.abs().sum().item(), 0)
                configure_objective(model, shape)
                control, control_rec, measured, predicted = [torch.tensor([[float(v)]]) for v in (4, 5, 6, 3)]
                mu1, mu2, mu_pred, logvar = [torch.tensor([[float(v)]]) for v in (1, 2, 3, 0)]
                loss, *parts = model.loss(
                    control, control_rec, mu1, logvar, measured, None, mu2, logvar, predicted, mu_pred, logvar,
                )
                self.assertEqual([value.item() for value in parts], [1., 0., 9., .5, 2., .5, 3.])
                self.assertAlmostEqual(loss.item(), 24.55, places=5)

    def test_training_seed_changes_initial_parameters(self):
        seed_17 = self.build("4_1+3", seed=17)
        seed_42 = self.build("4_1+3", seed=42)
        self.assertFalse(torch.equal(seed_17.prediction_head.projection.weight,
                                     seed_42.prediction_head.projection.weight))


if __name__ == "__main__":
    unittest.main()
