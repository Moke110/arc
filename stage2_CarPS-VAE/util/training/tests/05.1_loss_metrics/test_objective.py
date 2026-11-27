"""Check the published loss forms through the callable fit objective."""

from fractions import Fraction
from pathlib import Path
import sys
from unittest import TestCase

import torch

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
sys.path.insert(0, str(REPO / "stage2_CarPS-VAE/05_tsg_loss/05.1_loss_metrics/execution"))

from objective import delta_loss


class DeltaLossTests(TestCase):
    def test_mse_sums_both_conditions_and_all_genes(self):
        predicted = torch.tensor([[1., 0.], [2., 4.]])
        measured = torch.tensor([[0., 1.], [1., 2.]])
        self.assertEqual(delta_loss("mse", predicted, measured).item(), 7.)

    def test_alternatives_follow_the_published_reductions(self):
        predicted = torch.tensor([[1., 0.], [2., 4.]])
        measured = torch.tensor([[0., 1.], [1., 2.]])
        self.assertEqual(delta_loss("l1", predicted, measured).item(), 5.)
        self.assertAlmostEqual(delta_loss("cos", predicted, measured).item(), 2., places=6)
        self.assertAlmostEqual(delta_loss("pcc", predicted, measured).item(), 4., places=6)

    def test_angular_losses_have_finite_gradients_on_degenerate_vectors(self):
        for form in ("cos", "pcc"):
            for constant in (0., 1., 1e-10):
                with self.subTest(form=form, constant=constant):
                    predicted = torch.full((2, 3), constant, requires_grad=True)
                    measured = torch.tensor([[1., 3., 0.], [0., 0., 0.]])
                    value = delta_loss(form, predicted, measured)
                    value.backward()
                    self.assertTrue(torch.isfinite(value).item())
                    self.assertTrue(torch.isfinite(predicted.grad).all().item())


class ObjectiveTests(TestCase):
    @classmethod
    def setUpClass(cls):
        from run_loss import axis_run
        from os import environ
        source = Path(environ.get("ARC_TRANSIGEN_SOURCE", axis_run.REPO / "pert_model/TranSiGen/source"))
        if not source.is_dir():
            from unittest import SkipTest
            raise SkipTest("materialize TranSiGen source or set ARC_TRANSIGEN_SOURCE")
        cls.adapter = axis_run.kpgt.adapter(source, axis_run.kpgt.fit_transform([1, 10]))

    def test_native_anchor_exactly_preserves_scalar_and_gradients(self):
        from objective import Objective
        model = self.adapter.build(42, "cpu", n_genes=2)
        generator = torch.Generator().manual_seed(9)
        tensors = [torch.randn(3, 2, generator=generator, requires_grad=True) for _ in range(11)]
        native = model.loss(*tensors)
        actual = Objective("mse", Fraction(1))(model.loss, *tensors)
        self.assertTrue(torch.equal(native[0], actual[0]))
        expected_gradients = torch.autograd.grad(native[0], tensors)
        actual_gradients = torch.autograd.grad(actual[0], tensors)
        for expected, result in zip(expected_gradients, actual_gradients, strict=True):
            torch.testing.assert_close(result, expected, rtol=0, atol=0)

    def test_alternate_weights_and_reconstruction_paths(self):
        from objective import Objective
        model = self.adapter.build(42, "cpu", n_genes=2)
        x1 = torch.tensor([[1., 2.]])
        x2 = torch.tensor([[2., 5.]])
        x1_reconstructed = torch.tensor([[2., 0.]], requires_grad=True)
        x2_reconstructed = torch.tensor([[2., 3.]], requires_grad=True)
        predicted = torch.tensor([[4., 2.]], requires_grad=True)
        mean = torch.tensor([[1., 0.]])
        zero = torch.zeros_like(x1)
        mu_pred = torch.tensor([[2., 0.]])
        values = (x1, x1_reconstructed, mean, zero, x2, x2_reconstructed, mean, zero,
                  predicted, mu_pred, zero)
        total, *parts = Objective("l1", Fraction(3, 4))(model.loss, *values)
        self.assertEqual([part.item() for part in parts], [5., 4., 5., .5, .5, .5])
        self.assertAlmostEqual(total.item(), 12.9, places=5)
        gradients = torch.autograd.grad(total, (x1_reconstructed, x2_reconstructed, predicted))
        for actual, expected in zip(gradients, ([[2., -4.]], [[0., -4.]], [[.75, -.75]]), strict=True):
            torch.testing.assert_close(actual, torch.tensor(expected), rtol=0, atol=0)
