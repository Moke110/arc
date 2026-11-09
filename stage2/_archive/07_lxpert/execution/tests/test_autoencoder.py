"""Exercise bounded autoencoder fitting and independent reload via the CLI."""
import json
import subprocess
import sys
import unittest

import numpy as np
import torch

import test_inputs


class AutoencoderTest(unittest.TestCase):
    def setUp(self):
        self.input_fixture = test_inputs.InputReleaseTest()
        self.input_fixture.setUp()
        self.addCleanup(self.input_fixture.doCleanups)
        self.root = self.input_fixture.root
        self.release = self.input_fixture.release

    def run_autoencoder(self, action, *args, error=None):
        result = subprocess.run(
            [sys.executable, str(test_inputs.CLI), action, "--release", str(self.release),
             "--device", "cpu", *args], capture_output=True, text=True, check=False,
        )
        if error:
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn(error, result.stderr)
            return None
        self.assertEqual(result.returncode, 0, result.stderr)
        return json.loads(result.stdout)

    def test_short_fit_and_reload(self):
        self.input_fixture.prepare()
        output = self.root / "fit"
        fit = self.run_autoencoder("fit", "--output", str(output), "--fixture",
                                   "--steps", "1", "--batch-size", "2")
        self.assertEqual(fit["training_profiles"], {"tahoe": 1, "tenv": 1})
        self.assertEqual(fit["observed_genes"], {"tahoe": 2, "tenv": 3})
        self.assertEqual(fit["architecture"], {"gene_count": 3, "latent_tokens": 512,
                                              "width": 128, "heads": 8})
        reload_args = ("--checkpoint", str(output / "checkpoint.pt"), "--fixture")
        first = self.run_autoencoder("reload", *reload_args)
        self.assertEqual(first["latent_shape"], [2, 512, 128])
        self.assertEqual(first["reconstruction_shape"], [2, 3])
        self.assertTrue(first["finite"])
        self.assertLessEqual(first["max_abs_error"], first["atol"])
        self.assertEqual(first["axis_sha256"], fit["axis_sha256"])
        self.assertEqual(first["normalization"], fit["normalization"])

        checkpoint = torch.load(output / "checkpoint.pt", map_location="cpu", weights_only=True)
        for (source, row), prediction in zip(first["probe_rows"], checkpoint["probe_prediction"],
                                             strict=True):
            expression = np.load(self.release / f"{source}_expression.npy")[row]
            observed = [0, 2] if source == "tahoe" else [0, 1, 2]
            errors = [(float(prediction[gene]) - float(expression[gene])) ** 2
                      for gene in observed]
            self.assertAlmostEqual(first["probe_mse"][source], sum(errors) / len(errors), places=4)

        # A separately prepared release is still invalid if the gene order differs.
        fixture = self.input_fixture
        fixture.release = self.root / "reordered"
        axis = self.root / "tables/genes.csv"
        header, *genes = axis.read_text().splitlines()
        axis.write_text("\n".join([header, *reversed(genes)]) + "\n")
        fixture.prepare()
        self.release = fixture.release
        self.run_autoencoder("reload", *reload_args, error="checkpoint metadata mismatch")
