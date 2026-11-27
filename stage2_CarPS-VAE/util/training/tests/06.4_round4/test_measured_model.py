"""Exercise the registered Round 4 build, prediction and objective seams."""

import os
from pathlib import Path
import sys
from unittest import TestCase

import torch

HERE = (Path(__file__).resolve().parents[4] / "06_tsg_shape/06.4_round4/execution/tests" / Path(__file__).name)
sys.path.insert(0, str(HERE.parents[1]))

from measured_model import MeasuredAdapter, configure_objective
from measured_roster import (CONTROL_ONLY, MEASURED_SHARED, MEASURED_SEPARATE, SOURCE_SHAPES,
                             coordinates, fit_coordinates)
from composition_model import CompositionAdapter

SOURCE = next(parent for parent in HERE.parents if (parent / "pert_model/TranSiGen/source/src/model.py").is_file())
SOURCE = SOURCE / "pert_model/TranSiGen/source"
DEVICE = os.environ.get("TSG_MEASURED_TEST_DEVICE", "cpu")
TRANSFORM = dict(kind="training-log10-dose-population-standardization-v1", mean=0., std=1.)


def build(shape, device=DEVICE):
    model = MeasuredAdapter(SOURCE, TRANSFORM, shape).build(17, device, n_genes=4)
    configure_objective(model, shape)
    return model.eval()


def objective_parts(model, control, measured, control_rec, predicted, mu1, logvar1, mu_pred, logvar_pred):
    latent, mu2, logvar2 = model.encode_x2(measured)
    return model.loss(control, control_rec, mu1, logvar1, measured, model.decode_x2(latent),
                      mu2, logvar2, predicted, mu_pred, logvar_pred)


def tensor(values):
    return torch.tensor(values, device=DEVICE)


class MeasuredModelTests(TestCase):
    def test_registered_matrix(self):
        self.assertEqual(len(coordinates(False)), 45)
        self.assertEqual(len(fit_coordinates(False)), 27)
        self.assertEqual(len(fit_coordinates(True)), 20)

    def test_control_only_omits_measured_graph_and_prior_on_prediction(self):
        model = build(CONTROL_ONLY)
        self.assertEqual([name for name, _ in model.named_children() if name.startswith("encoder")], ["encoder_x1"])
        self.assertFalse(any(name.startswith(("mu_z2.", "logvar_z2.", "decoder_x2.")) for name in model.state_dict()))
        control, measured, control_rec, predicted = [tensor([[value]]) for value in (4., 6., 5., 3.)]
        mu1, logvar1, mu_pred, logvar_pred = [tensor([[value]]) for value in (1., 0., 200., 90.)]
        loss, *parts = objective_parts(model, control, measured, control_rec, predicted, mu1, logvar1, mu_pred, logvar_pred)
        self.assertEqual([value.item() for value in parts], [1., 0., 9., .5, 0., 0., 3.])
        self.assertAlmostEqual(loss.item(), 24.30, places=5)

    def test_control_only_prediction_samples_the_gaussian(self):
        model = build(CONTROL_ONLY)
        control, features = torch.ones(2, 4, device=DEVICE), torch.ones(2, 2305, device=DEVICE)
        with torch.no_grad():
            for head in (model.mu_z2Fz1[0], model.logvar_z2Fz1[0]):
                head.weight.zero_()
                head.bias.zero_()
            torch.manual_seed(3)
            first = model(control, features)
            torch.manual_seed(19)
            second = model(control, features)
        torch.testing.assert_close(first[4], second[4])
        torch.testing.assert_close(first[5], second[5])
        self.assertFalse(torch.equal(first[6], second[6]))

    def test_measured_encoder_uses_signed_difference_only_in_objective(self):
        control = torch.full((2, 4), 3., device=DEVICE)
        measured = tensor([[1., 3., 5., 0.], [4., 5., 5., 2.]])
        features = torch.ones(2, 2305, device=DEVICE)
        for shape in (MEASURED_SHARED, MEASURED_SEPARATE):
            with self.subTest(shape=shape):
                model = build(shape)
                observed = []
                model.encoder_delta.register_forward_pre_hook(lambda _module, inputs: observed.append(inputs[0]))
                control_rec, mu1, logvar1, predicted, mu_pred, logvar_pred, _ = model(control, features)
                self.assertEqual(observed, [])
                loss = objective_parts(model, control, measured, control_rec, predicted, mu1, logvar1,
                                       mu_pred, logvar_pred)[0]
                self.assertEqual(len(observed), 1)
                torch.testing.assert_close(observed[0], tensor([[-2., 0., 2., -3.], [1., 2., 2., -1.]]))
                loss.backward()
                self.assertGreater(model.mu_delta[0].weight.grad.abs().sum().item(), 0)

    def test_signed_prediction_and_decoder_sharing(self):
        signed_delta = tensor([-2., 2., -1., 1.])
        control = torch.full((2, 4), 3., device=DEVICE)
        measured = tensor([[1., 3., 5., 0.], [4., 5., 5., 2.]])
        features = torch.ones(2, 2305, device=DEVICE)
        for shape in (CONTROL_ONLY, MEASURED_SHARED, MEASURED_SEPARATE):
            with self.subTest(shape=shape), torch.no_grad():
                model = build(shape)
                head = model.prediction_decoder[-1]
                head.weight.zero_()
                head.bias.copy_(signed_delta)
                predicted = model(control, features)[3]
                torch.testing.assert_close(predicted, control + signed_delta)
                if shape != CONTROL_ONLY:
                    self.assertIs(model.prediction_decoder is model.measurement_decoder, shape == MEASURED_SHARED)
                    if shape == MEASURED_SEPARATE:
                        reconstructed = model.measurement_decoder[-1]
                        reconstructed.weight.zero_()
                        reconstructed.bias.copy_(-signed_delta)
                    mu1, logvar1, mu_pred, logvar_pred = [tensor([[value]]) for value in (1., 0., 3., 0.)]
                    model.mu_delta[0].weight.zero_()
                    model.mu_delta[0].bias.fill_(2.)
                    model.logvar_delta[0].weight.zero_()
                    model.logvar_delta[0].bias.zero_()
                    loss, *parts = objective_parts(model, control, measured, control + 1, predicted,
                                                   mu1, logvar1, mu_pred, logvar_pred)
                    self.assertEqual([value.item() for value in parts], [
                        8., 51. if shape == MEASURED_SHARED else 43., 51., .5, 400., 100., 17.,
                    ])
                    self.assertAlmostEqual(loss.item(), (189.8 if shape == MEASURED_SHARED else 181.8), places=4)

    def test_references_preserve_seeded_predictions(self):
        control, features = torch.ones(2, 4, device=DEVICE), torch.ones(2, 2305, device=DEVICE)
        for shape, source in SOURCE_SHAPES.items():
            with self.subTest(shape=shape), torch.no_grad():
                model = build(shape)
                reference = CompositionAdapter(SOURCE, TRANSFORM, source).build(17, DEVICE, n_genes=4).eval()
                self.assertEqual(model.state_dict().keys(), reference.state_dict().keys())
                torch.manual_seed(703)
                actual = model(control, features)
                torch.manual_seed(703)
                expected = reference(control, features)
                for left, right in zip(actual, expected, strict=True):
                    self.assertTrue(torch.equal(left, right))
