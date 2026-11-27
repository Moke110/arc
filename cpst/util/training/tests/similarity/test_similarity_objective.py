"""Verify similarity through the registered matched-control loss boundary."""

from fractions import Fraction
from pathlib import Path
import sys
from unittest import TestCase

import torch

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
sys.path.insert(0, str(REPO / "cpst/05_tsg_loss/similarity"))

from similarity_objective import SimilarityObjective
from reconstruction_objective import ReconstructionObjective


def native(control, control_rec, mu1, logvar1,
           measured, treatment_rec, mu2, logvar2, predicted, mu_pred, logvar_pred):
    control_mse = (control_rec - control).square().sum()
    treatment_mse = (treatment_rec - measured).square().sum()
    delta_mse = (predicted - measured).square().sum()
    kl1, kl2, kl_prediction = (control.new_tensor(value) for value in (1., 2., 3.))
    total = control_mse + treatment_mse + delta_mse + 0.1 * kl1 + 0.1 * kl2 + 0.1 * kl_prediction
    return total, control_mse, treatment_mse, delta_mse, kl1, kl2, kl_prediction


def inputs(prediction, truth):
    measured_delta = torch.tensor(truth, dtype=torch.float64)
    control = torch.full_like(measured_delta, 5.)
    measured = control + measured_delta
    predicted = (control + torch.tensor(prediction, dtype=control.dtype)).requires_grad_()
    control_rec = (control + 1).requires_grad_()
    treatment_rec = (measured + 1).requires_grad_()
    return (control, control_rec, control, control, measured, treatment_rec,
            measured, measured, predicted, predicted, predicted)


class SimilarityObjectiveTests(TestCase):
    def test_cosine_is_bounded_and_averaged_over_conditions(self):
        truth = [[1., 0.]]
        for prediction, expected in (([[1., 0.]], 0.), ([[0., 1.]], .5),
                                     ([[-1., 0.]], 1.), ([[0., 0.]], .5)):
            with self.subTest(prediction=prediction):
                self.check_increment("cos", prediction, truth, expected)
        self.check_increment("cos", [[1., 0.], [-1., 0.], [0., 0.]], truth * 3, .5)

    def check_increment(self, metric, prediction, truth, expected):
        tensors = inputs(prediction, truth)
        base, *base_components = ReconstructionObjective(Fraction(1))(native, *tensors)
        loss, *components = SimilarityObjective(metric, Fraction(1, 2))(native, *tensors)
        self.assertAlmostEqual((loss - base).item(), .5 * expected)
        torch.testing.assert_close(components, base_components, rtol=0, atol=0)
        self.assertTrue(torch.isfinite(torch.autograd.grad(loss, tensors[8])[0]).all())

    def test_pcc_centers_each_vector_and_handles_constants(self):
        self.check_increment("pcc", [[10., 13., 16.]], [[-4., -3., -2.]], 0.)
        self.check_increment("pcc", [[16., 13., 10.]], [[-4., -3., -2.]], 1.)
        self.check_increment("pcc", [[10., 10., 10.]], [[-4., -3., -2.]], .5)
        self.check_increment("pcc", [[11., 8., 11.]], [[-4., -3., -2.]], .5)

    def test_zero_preserves_the_predecessor_loss_and_gradients(self):
        for metric in ("pcc", "cos"):
            with self.subTest(metric=metric):
                tensors = inputs([[1., 3.], [0., 0.]], [[1., 4.], [1., 0.]])
                diff = tensors[1], tensors[5], tensors[8]
                expected = ReconstructionObjective(Fraction(1))(native, *tensors)
                actual = SimilarityObjective(metric, Fraction(0))(native, *tensors)
                torch.testing.assert_close(actual, expected, rtol=0, atol=0)
                torch.testing.assert_close(torch.autograd.grad(actual[0], diff),
                                           torch.autograd.grad(expected[0], diff), rtol=0, atol=0)
