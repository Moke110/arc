"""XPert compression architecture and full-axis CLI contract."""

import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

import numpy as np
import torch

from test_gene_axes import build
from test_landmark import make_fixture, run

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "script"))
from model import XPert, loss_terms


WIDTHS = (512, 1024, 2048, 4096)
CHANNELS = 256


class CompressionTest(unittest.TestCase):
    def test_all_widths_restore_the_real_gene_axis_and_share_each_map(self):
        with torch.device("meta"):
            identity = torch.empty(18164, CHANNELS)
            control = torch.empty(2, 18164)
            morgan = torch.empty(2, 2048)
            dose = torch.empty(2, dtype=torch.long)
            for width in WIDTHS:
                with self.subTest(width=width):
                    model = XPert(identity, latent_tokens=width).eval()
                    shapes = {}
                    for name in ("encoder", "decoder"):
                        def record(module, args, result, name=name):
                            shapes.setdefault(name, []).append((args[0].shape, result.shape))
                        getattr(model, name).register_forward_hook(record)
                    outputs = model(control, morgan, dose)
                    assert [value.shape for value in outputs] == [(2, 18164)] * 3
                    assert shapes["encoder"] == [((2, 18164, CHANNELS), (2, width, CHANNELS))]
                    assert shapes["decoder"] == [((2, width, CHANNELS), (2, 18164, CHANNELS))] * 2
                    assert sum(p.numel() for m in (model.encoder, model.decoder) for p in m.parameters()) == (
                        2 * 18164 * width + 2 * width**2 + 3 * width + 18164
                    )
                    terms = loss_terms(outputs, torch.empty_like(control), control)
                    weights = self.map_weights(model)
                    for i, term in enumerate(terms):
                        gradients = torch.autograd.grad(term, weights, retain_graph=i < 3)
                        assert [g.shape for g in gradients] == [p.shape for p in weights]

    @staticmethod
    def map_weights(model):
        weights = [p for layer in (model.encoder, model.decoder) for p in layer.parameters()
                   if p.ndim == 2]
        assert len(weights) == 4
        return weights

    def test_each_objective_trains_both_maps_in_a_reduced_axis_fit(self):
        generator = torch.Generator().manual_seed(23)
        identity = torch.randn(7, CHANNELS, generator=generator)
        control = torch.rand(1, 7, generator=generator)
        treated = torch.rand(1, 7, generator=generator)
        morgan = torch.ones(1, 2048)
        dose = torch.tensor([4])
        torch.manual_seed(17)
        uncompressed = XPert(identity)
        torch.manual_seed(17)
        model = XPert(identity, latent_tokens=3).eval()
        compressed = model.state_dict()
        assert all(torch.equal(value, compressed[name])
                   for name, value in uncompressed.state_dict().items())
        same_channels = control.unsqueeze(-1).expand(-1, -1, CHANNELS)
        encoded = model.encoder(same_channels)
        torch.testing.assert_close(encoded[..., 0], encoded[..., -1])
        terms = loss_terms(model(control, morgan, dose), treated, control)
        for i, term in enumerate(terms):
            gradients = torch.autograd.grad(term, self.map_weights(model), retain_graph=i < 3)
            assert all(torch.isfinite(g).all() and torch.count_nonzero(g) for g in gradients)

    def test_fit_and_replay_use_the_existing_full_axis_predictor(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            native, tenv = make_fixture(root)
            axes, _ = build(root, native, tenv)
            full = axes / "views/G18164/predictor"
            small = axes / "views/G2000/predictor"
            policy = root / "policy.json"
            policy.write_text(json.dumps({
                "kind": "qualification", "device": "cpu", "kernel": "math", "dtype": "float32",
                "microbatch": 1, "effective_batch": 1,
                "phases": {"tahoe": {"updates": 1, "validation_interval": 1, "learning_rate": .0001},
                           "tenv": {"updates": 1, "validation_interval": 1, "learning_rate": .0001}},
            }))
            fit, prediction = root / "fit", root / "prediction"
            fit_args = ("fit", "--model", "XPert_MLP512", "--policy", policy,
                        "--fold", 0, "--seed", 17, "--output", fit)
            rejected = run("script/xpert.py", *fit_args, "--predictor", small, check=False)
            assert rejected.returncode != 0
            assert "requires the G18164 predictor" in rejected.stderr
            run("script/xpert.py", *fit_args, "--predictor", full)
            run("script/xpert.py", "predict", "--predictor", full, "--checkpoint", fit / "tenv.pt",
                "--role", "validation", "--output", prediction)
            run("script/xpert.py", "replay", "--predictor", full, "--checkpoint", fit / "tenv.pt",
                "--prediction", prediction)
            report = json.loads((fit / "run.json").read_text())
            assert report["model"] == "XPert_MLP512"
            assert report["architecture"]["latent_tokens"] == 512
            assert report["architecture"]["gene_count"] == 3
            assert np.load(prediction / "delta.npy").shape[1] == 3
