"""Verify the objective at the matched-control loss boundary."""

from fractions import Fraction
from pathlib import Path
import sys
from unittest import TestCase

import torch

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
EXECUTION = REPO / "cpst/05_tsg_loss/05.3_mse-l1_mse-x/execution"
sys.path.insert(0, str(EXECUTION))

from reconstruction_objective import ReconstructionObjective
from additive_objective import AdditiveObjective


class ReconstructionObjectiveTests(TestCase):
    def setUp(self):
        control = torch.tensor([[4., 6.], [1., 8.]])
        measured = torch.tensor([[5., 5.], [2., 10.]])
        self.control_rec = torch.tensor([[5., 8.], [3., 9.]], requires_grad=True)
        self.treatment_rec = torch.tensor([[7., 7.], [4., 12.]], requires_grad=True)
        self.predicted = torch.tensor([[6., 4.], [3., 12.]], requires_grad=True)
        self.inputs = (control, self.control_rec, control, control,
                       measured, self.treatment_rec, measured, measured,
                       self.predicted, self.predicted, self.predicted)

    @staticmethod
    def native(control, control_rec, mu1, logvar1,
               measured, treatment_rec, mu2, logvar2, predicted, mu_pred, logvar_pred):
        control_mse = (control_rec - control).square().sum()
        treatment_mse = (treatment_rec - measured).square().sum()
        delta_mse = ((predicted - control) - (measured - control)).square().sum()
        kl1, kl2, kl_prediction = (control.new_tensor(value) for value in (1., 2., 3.))
        total = control_mse + treatment_mse + delta_mse + 0.1 * kl1 + 0.1 * kl2 + 0.1 * kl_prediction
        return total, control_mse, treatment_mse, delta_mse, kl1, kl2, kl_prediction

    def test_half_weight_changes_only_the_reconstruction_gradient(self):
        total, *components = ReconstructionObjective(Fraction(1, 2))(self.native, *self.inputs)
        self.assertAlmostEqual(total.item(), 29.35, places=5)
        self.assertEqual([value.item() for value in components], [10., 16., 7., 1., 2., 3., 5.])
        gradients = torch.autograd.grad(total, (self.control_rec, self.treatment_rec, self.predicted))
        expected = (torch.tensor([[1., 2.], [2., 1.]]),
                    torch.tensor([[2., 2.], [2., 2.]]),
                    torch.tensor([[3.75, -3.75], [3.75, 5.75]]))
        torch.testing.assert_close(gradients, expected, rtol=0, atol=0)

    def test_zero_omits_both_reconstruction_terms_but_keeps_raw_components(self):
        total, *components = ReconstructionObjective(Fraction(0))(self.native, *self.inputs)
        self.assertAlmostEqual(total.item(), 16.35, places=5)
        self.assertEqual([value.item() for value in components], [10., 16., 7., 1., 2., 3., 5.])
        gradients = torch.autograd.grad(total, (self.control_rec, self.treatment_rec), allow_unused=True)
        self.assertEqual(gradients, (None, None))

    def test_one_preserves_step_05_2_loss_and_gradients_exactly(self):
        expected = AdditiveObjective(Fraction(7, 4))(self.native, *self.inputs)
        actual = ReconstructionObjective(Fraction(1))(self.native, *self.inputs)
        torch.testing.assert_close(actual, expected, rtol=0, atol=0)
        differentiable = self.control_rec, self.treatment_rec, self.predicted
        torch.testing.assert_close(torch.autograd.grad(actual[0], differentiable),
                                   torch.autograd.grad(expected[0], differentiable),
                                   rtol=0, atol=0)
