"""Exercise the additive objective with a known pair of matched-control Deltas."""

from fractions import Fraction
from pathlib import Path
import sys
from unittest import TestCase

import torch

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
sys.path.insert(0, str(REPO / "stage2_CarPS-VAE/05_tsg_loss/05.2_mse_l1/execution"))

from additive_objective import AdditiveObjective, ComponentTrace


class AdditiveObjectiveTests(TestCase):
    def setUp(self):
        control = torch.tensor([[4., 6.], [1., 8.]])
        measured = torch.tensor([[5., 5.], [2., 10.]])
        predicted = torch.tensor([[6., 4.], [3., 12.]], requires_grad=True)
        self.inputs = (control, control, control, control, measured, measured,
                       measured, measured, predicted, predicted, predicted)

        # Use known native components while retaining a differentiable native scalar.
        def native(*inputs):
            delta_mse = (inputs[8] - inputs[4]).square().sum()
            parts = tuple(torch.tensor(value) for value in (10., 20., 1., 2., 3.))
            return 30.6 + delta_mse, *parts[:2], delta_mse, *parts[2:]

        self.native = native

    def test_zero_uses_the_native_scalar_and_still_traces_l1(self):
        expected = self.native(*self.inputs)[0]
        trace = ComponentTrace(lambda *inputs: AdditiveObjective(Fraction(0))(self.native, *inputs))
        actual = trace(*self.inputs)[0]
        torch.testing.assert_close(actual, expected, rtol=0, atol=0)
        gradients = torch.autograd.grad(actual, self.inputs[8])[0]
        torch.testing.assert_close(gradients, torch.tensor([[2., -2.], [2., 4.]]), rtol=0, atol=0)
        metrics = trace.per_condition("train")
        self.assertEqual(metrics["train_delta_mse"], 3.5)
        self.assertEqual(metrics["train_delta_l1"], 2.5)

    def test_positive_weight_adds_the_full_batch_l1_sum(self):
        actual, *parts = AdditiveObjective(Fraction(3, 4))(self.native, *self.inputs)
        self.assertAlmostEqual(actual.item(), 41.35, places=5)
        self.assertEqual([value.item() for value in parts], [10., 20., 7., 1., 2., 3., 5.])
        gradients = torch.autograd.grad(actual, self.inputs[8])[0]
        torch.testing.assert_close(gradients, torch.tensor([[2.75, -2.75], [2.75, 4.75]]), rtol=0, atol=0)
