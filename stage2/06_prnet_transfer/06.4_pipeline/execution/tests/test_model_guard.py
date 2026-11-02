"""Regression tests for the ticket-#9 variance-head guard relaxation.

``check_architecture`` probes the model with an all-zero batch before any GPU time is
spent. Softplus is mathematically strictly positive, but float32 underflows to exactly 0
for a very negative decoder logit, and an off-distribution zero probe is not evidence that
the model is broken (four Tahoe-init checkpoints in ticket #9 hit this while their
real-input variance stayed positive and the evaluated metric uses the predicted mean
only). The guard must therefore reject non-finite or negative variance, and accept an
exact zero.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import torch

EXECUTION = Path(__file__).resolve().parents[1]
if str(EXECUTION) not in sys.path:
    sys.path.insert(0, str(EXECUTION))

from ptlib import model as model_mod  # noqa: E402

STRUCT = {
    "n_genes": 32,
    "hidden": [16],
    "z_dim": 8,
    "c_dim": 8,
    "n_dim": 4,
    "drug_dim": 8,
    "adaptor": [16],
    "dropout": 0.05,
}


def _model() -> model_mod.PGM:
    torch.manual_seed(7)
    return model_mod.build_model(STRUCT)


class VarianceHeadGuardTests(unittest.TestCase):
    def test_negative_variance_head_bias_underflows_and_passes(self):
        model = _model()
        bias = model.decoder.recon_decoder.bias
        dim = bias.numel() // 2
        with torch.no_grad():
            bias[dim:] = -200.0
        # the zero probe really does underflow to exactly zero for every gene
        model.eval()
        with torch.no_grad():
            mean, var = model(
                torch.zeros(2, STRUCT["n_genes"]),
                torch.zeros(2, STRUCT["drug_dim"]),
                torch.zeros(2, STRUCT["n_dim"]),
            )
        self.assertTrue(bool(torch.isfinite(mean).all()))
        self.assertEqual(int((var == 0).sum()), var.numel())
        model_mod.check_architecture(model, STRUCT)  # must not raise

    def test_nan_variance_head_still_raises(self):
        model = _model()
        bias = model.decoder.recon_decoder.bias
        dim = bias.numel() // 2
        with torch.no_grad():
            bias[dim] = float("nan")
        with self.assertRaises(SystemExit) as caught:
            model_mod.check_architecture(model, STRUCT)
        self.assertIn("non-finite variance", str(caught.exception))

    def test_inf_variance_head_still_raises(self):
        model = _model()
        bias = model.decoder.recon_decoder.bias
        dim = bias.numel() // 2
        with torch.no_grad():
            bias[dim] = float("inf")
        with self.assertRaises(SystemExit) as caught:
            model_mod.check_architecture(model, STRUCT)
        self.assertIn("non-finite variance", str(caught.exception))

    def test_negative_variance_still_raises(self):
        model = _model()
        model.decoder.forward = lambda z: (  # type: ignore[assignment]
            torch.zeros(z.size(0), STRUCT["n_genes"]),
            -torch.ones(z.size(0), STRUCT["n_genes"]),
        )
        with self.assertRaises(SystemExit) as caught:
            model_mod.check_architecture(model, STRUCT)
        self.assertIn("negative variance", str(caught.exception))

    def test_architecture_mismatch_still_raises(self):
        model = _model()
        wrong = dict(STRUCT)
        wrong["z_dim"] = STRUCT["z_dim"] + 1
        with self.assertRaises(SystemExit) as caught:
            model_mod.check_architecture(model, wrong)
        self.assertIn("architecture mismatch", str(caught.exception))

    def test_default_model_passes(self):
        model_mod.check_architecture(_model(), STRUCT)


if __name__ == "__main__":
    unittest.main()
