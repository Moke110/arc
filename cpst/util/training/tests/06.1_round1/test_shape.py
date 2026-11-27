import os
from pathlib import Path
import sys
import unittest

import torch

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
EXECUTION = REPO / "cpst/06_tsg_module/06.1_round1/execution"
sys.path.insert(0, str(EXECUTION))
from shape_model import ShapeAdapter, axis_run
from shape_roster import SHAPES


class ShapeTest(unittest.TestCase):
    def test_native_model_build_matches_registered_architectures(self):
        coordinates = {
            "original": ((1200,), (800,), 100),
            "widen_z-110": ((1200,), (800,), 110),
            "widen_z-120": ((1200,), (800,), 120),
            "widen_z-130": ((1200,), (800,), 130),
            "widen_enc-2k": ((2000,), (800,), 100),
            "widen_enc-4k": ((4000,), (800,), 100),
            "widen_enc-6k": ((6000,), (800,), 100),
            "deepen_enc": ((1200, 1200), (800,), 100),
            "deepen_dec": ((1200,), (800, 800), 100),
            "deepen_both": ((1200, 1200), (800, 800), 100),
        }
        self.assertEqual(SHAPES, coordinates)
        source = Path(os.environ.get("TSG_SOURCE", axis_run.REPO / "pert_model/TranSiGen/source"))
        transform = axis_run.kpgt.fit_transform([1., 10.])
        for key, (encoder, decoder, latent) in coordinates.items():
            with self.subTest(shape=key):
                network = ShapeAdapter(source, transform, key).build(42, "cpu").eval()
                for branch in (network.encoder_x1, network.encoder_x2):
                    self.assert_blocks(branch, (5000, *encoder))
                for branch in (network.decoder_x1, network.decoder_x2):
                    self.assert_blocks(branch[:-2], (latent, *decoder))
                    self.assertEqual((branch[-2].in_features, branch[-2].out_features), (decoder[-1], 5000))
                self.assertNotEqual(network.encoder_x1[0].weight.data_ptr(), network.encoder_x2[0].weight.data_ptr())
                self.assertNotEqual(network.decoder_x1[0].weight.data_ptr(), network.decoder_x2[0].weight.data_ptr())
                for head in (network.mu_z1, network.logvar_z1, network.mu_z2, network.logvar_z2):
                    self.assertEqual((head[0].in_features, head[0].out_features), (encoder[-1], latent))
                for head in (network.mu_z2Fz1, network.logvar_z2Fz1):
                    self.assertEqual((head[0].in_features, head[0].out_features), (latent + 432, latent))
                calls = []
                with (network.decoder_x2.register_forward_hook(lambda *args: calls.append(args[-1])),
                      torch.no_grad()):
                    treated = torch.zeros(2, 5000)
                    observed, *_ = network.encode_x2(treated)
                    self.assertIs(network.decode_x2(observed), calls[0])
                    pred = network(treated, torch.zeros(2, axis_run.kpgt.WIDTH + 1))
                    self.assertIs(pred[3], calls[1])
                    self.assertEqual(pred[-1].shape, (2, latent))
                del network
        original = axis_run.kpgt.Adapter(source, transform).build(42, "cpu")
        refit = ShapeAdapter(source, transform, "original").build(42, "cpu")
        self.assertEqual(list(original.state_dict()), list(refit.state_dict()))
        for name, value in original.state_dict().items():
            self.assertTrue(torch.equal(value, refit.state_dict()[name]), name)

    def assert_blocks(self, layers, widths):
        self.assertEqual(len(layers), 4 * (len(widths) - 1))
        for index, (input_width, output_width) in enumerate(zip(widths, widths[1:])):
            block = layers[index * 4:index * 4 + 4]
            self.assertEqual(tuple(map(type, block)),
                             (torch.nn.Linear, torch.nn.BatchNorm1d, torch.nn.ReLU, torch.nn.Dropout))
            self.assertEqual((block[0].in_features, block[0].out_features), (input_width, output_width))
            self.assertEqual(block[1].num_features, output_width)
            self.assertEqual(block[3].p, 0.1)
