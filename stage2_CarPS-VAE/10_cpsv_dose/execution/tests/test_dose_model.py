"""Check the registered dose model forward seam against explicit inputs."""

import os
from pathlib import Path
import sys
import unittest

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dose_model import DoseAdapter, model_features
from dose_protocol import ARMS, DEFAULT_SOURCE

SOURCE = Path(os.getenv("CARPS_TEST_TRANSIGEN_SOURCE", DEFAULT_SOURCE))


class DoseModelTests(unittest.TestCase):
    def test_dose_transport(self):
        chemical = np.ones((3, 512), dtype=np.float32)
        doses = np.array([0.2, 1.0, 10.0])
        for arm in ARMS:
            if arm == "TranSiGen":
                continue
            with self.subTest(arm=arm):
                features = model_features(arm, chemical, doses)
                if arm == "no-dose":
                    np.testing.assert_array_equal(features, chemical)
                else:
                    self.assertEqual(features.shape, (3, 513))
                    np.testing.assert_array_equal(features[:, :512], chemical)
                    expected = ([0.07918125, 0.30103, 1.0413927] if arm.endswith("scaling")
                                else [-0.69897, 0, 1])
                    np.testing.assert_allclose(features[:, -1], expected, atol=1e-7)

    def test_forward_configuration_and_film_identity(self):
        torch.set_num_threads(1)
        control = torch.ones(2, 7)
        chemistry = torch.ones(2, 512)
        torch.manual_seed(12)
        reference = DoseAdapter(SOURCE, "no-dose").build(17, "cpu", n_genes=7).eval()
        before = torch.get_rng_state()
        with torch.no_grad():
            baseline = reference(control, chemistry)[3]
        for arm in ARMS[1:]:
            with self.subTest(arm=arm):
                model = DoseAdapter(SOURCE, arm).build(17, "cpu", n_genes=7).eval()
                self.assertTrue(torch.equal(torch.get_rng_state(), before))
                for name, value in model.state_dict().items():
                    other = reference.state_dict().get(name)
                    if other is not None and value.shape == other.shape:
                        self.assertTrue(torch.equal(value, other), name)
                self.assertEqual(model.feat_embeddings[0].in_features, 544 if arm.startswith("Xdose") else 512)
                self.assertEqual(model.mu_z2Fz1[0].in_features, 532 if arm.startswith("Zdose") else 500)
                if arm.endswith("FiLM"):
                    self.assertEqual(model.conditioner[-1].out_features, 1024 if arm.startswith("X") else 800)
                    features = torch.cat((chemistry, torch.tensor([[0.0], [1.0]])), dim=1)
                    with torch.no_grad():
                        self.assertTrue(torch.equal(model(control, features)[3], baseline))


if __name__ == "__main__":
    unittest.main()
