"""Checkpoint selection at the one-coordinate training seam."""
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

import numpy as np
import pandas as pd
import torch

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP / "execution"))

from run import PATIENCE, fit


class ValidationInputs:
    def __init__(self):
        self.control = np.zeros((65, 1), dtype=np.float32)
        self.expression = self.control.copy()
        self.expression[-1] = 3

    def fingerprints(self, rows):
        return self.expression[rows]

    def doses(self, rows):
        return np.ones(len(rows))

    def model_input(self, partition):
        return self


class NativeLossFixture(torch.nn.Module):
    """Deterministic stand-in for the external author model's batch sum contract."""

    def __init__(self):
        super().__init__()
        self.epoch = torch.nn.Parameter(torch.tensor(0.0), requires_grad=False)

    def forward(self, control, features):
        assert not self.training and not torch.is_grad_enabled()
        rec = control if self.epoch == 1 else control + 2
        pred = control + 1 if self.epoch == 1 else features[:, :1]
        return rec, control, control, pred, control, control, control

    def encode_x2(self, measured):
        zeros = torch.zeros_like(measured)
        return measured, zeros, zeros

    def decode_x2(self, latent):
        return latent

    def loss(self, control, rec, mu1, logvar1, measured, measured_rec, mu2, logvar2,
             pred, mu_pred, logvar_pred):
        return ((rec - control).square().sum() + (pred - measured).square().sum(),)


class TrainingFixture:
    def optimizer(self, model):
        return None

    def train_epoch(self, model, inputs, rows, optimizer, rng):
        model.train()
        model.epoch.add_(1)
        return 1.0

class ValidationTest(unittest.TestCase):
    def test_checkpoint_minimizes_per_condition_author_loss_with_30_epoch_patience(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            transform = {"kind": "training-log10-dose-population-standardization-v1", "mean": 0.0, "std": 1.0}
            config = {"partitions": {"train": [0, 1], "valid": list(range(65))},
                      "patience": PATIENCE, "max_epochs": 100, "seed": 42, "validation_seed": 22617,
                      "model": {"batch_size": 64, "transform": transform}}
            record = fit(root, TrainingFixture(), NativeLossFixture(), ValidationInputs(), config, "fixture")
            history = pd.read_csv(root / "learning_curve.csv")
            checkpoint = torch.load(root / "checkpoint.pt", weights_only=True)
            prediction = np.load(root / "validation_prediction.npy")

        self.assertEqual(record["best_epoch"], 1)
        self.assertEqual(checkpoint["epoch"], 1)
        self.assertEqual(record["epochs"], 31)
        self.assertAlmostEqual(record["validation_loss"], 68 / 65)
        self.assertAlmostEqual(history.validation_loss.iloc[0], 68 / 65)
        self.assertEqual(history.validation_mse.iloc[1], 0)
        np.testing.assert_array_equal(prediction, np.ones((65, 1), dtype=np.float32))
