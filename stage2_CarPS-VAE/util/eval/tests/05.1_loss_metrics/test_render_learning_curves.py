"""Render the complete formal learning-curve grid directly from CSV."""

import csv
from itertools import product
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from unittest import TestCase

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
RENDER = REPO / "stage2_CarPS-VAE/05_tsg_loss/05.1_loss_metrics" / "execution/render_learning_curves.py"
FIELDS = ("loss", "alpha", "fold", "seed", "run_sha256", "best_epoch", "stopping_epoch", "epoch",
          "train_loss", "validation_loss", "validation_mse")


class LearningCurveFigureTests(TestCase):
    def test_csv_requires_each_observed_epoch_and_renders_the_full_grid(self):
        with TemporaryDirectory() as scratch:
            root = Path(scratch)
            source = root / "epoch-objectives.csv"
            output = root / "figures"
            command = [sys.executable, str(RENDER), "--data", str(root), "--output", str(output)]
            coordinates = product(("mse", "cos", "pcc", "l1"), ("0.75", "1", "1.25"), range(3))
            rows = [(loss, alpha, fold, 42, f"{index:064x}", 2, 32, epoch, 2 + (epoch - 2) ** 2,
                     1 + (epoch - 2) ** 2, .5)
                    for index, (loss, alpha, fold) in enumerate(coordinates, 1)
                    for epoch in range(1, 33)]

            with source.open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(FIELDS)
                writer.writerows(rows[:-1])
            rejected = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(rejected.returncode, 0)
            self.assertIn("missing or padded epochs", rejected.stderr)
            self.assertFalse(output.exists())

            with source.open("a", newline="") as stream:
                csv.writer(stream).writerow(rows[-1])
            subprocess.run(command, capture_output=True, text=True, check=True)
            figure = (output / "learning-curves.png").read_bytes()
            self.assertEqual(figure[:8], b"\x89PNG\r\n\x1a\n")
