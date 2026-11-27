"""Check the public alpha constructor against its selected Step 06.4 reference."""

import os
from pathlib import Path
import sys
import unittest

import torch


STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP / "output/api"))

from cpst_alpha import AlphaAdapter, objective


class AlphaModelTest(unittest.TestCase):
    def test_default_preserves_selected_checkpoint_layout_and_predictions(self):
        from measured_model import MeasuredAdapter, configure_objective
        from measured_roster import MEASURED_SHARED
        from tsg_axis.model import fit_transform

        source = Path(os.environ.get("CARPS_TEST_SOURCE", STEP.parents[1] / "pert_model/TranSiGen/source"))
        transform = fit_transform([0.1, 1, 10])
        reference = MeasuredAdapter(source, transform, MEASURED_SHARED).build(17, "cpu", n_genes=8)
        configure_objective(reference, MEASURED_SHARED)
        alpha = AlphaAdapter(source, transform).build(17, "cpu", n_genes=8)
        self.assertEqual(tuple(reference.state_dict()), tuple(alpha.state_dict()))
        for name, weight in reference.state_dict().items():
            self.assertTrue(torch.equal(weight, alpha.state_dict()[name]), name)

        x1 = torch.linspace(0, 1, 16).reshape(2, 8)
        features = torch.linspace(-1, 1, 4610).reshape(2, 2305)
        reference.eval()
        alpha.eval()
        with torch.no_grad():
            torch.manual_seed(100042)
            expected = reference(x1, features)
            torch.manual_seed(100042)
            observed = alpha(x1, features)
        for old, new in zip(expected, observed, strict=True):
            self.assertTrue(torch.equal(old, new))
        self.assertEqual(objective()["delta_l1_weight"], 1.75)


    def test_native_width_projection_preserves_downstream_initialization_and_stream(self):
        from tsg_axis.model import fit_transform

        source = Path(os.environ.get("CARPS_TEST_SOURCE", STEP.parents[1] / "pert_model/TranSiGen/source"))
        transform = fit_transform([0.1, 1, 10])
        device = os.environ.get("CARPS_TEST_DEVICE", "cpu")

        def projection(width, projected_width):
            self.assertEqual((width, projected_width), (16, 400))
            layer = torch.nn.Linear(width, projected_width)
            torch.nn.init.zeros_(layer.weight)
            return layer

        native = AlphaAdapter(source, transform).build(42, device, n_genes=8)
        cpu_rng = torch.get_rng_state()
        cuda_rng = torch.cuda.get_rng_state(device) if torch.device(device).type == "cuda" else None
        variant = AlphaAdapter(source, transform, chemical_width=16, projection_factory=projection)
        model = variant.build(42, device, n_genes=8)
        self.assertTrue(torch.equal(cpu_rng, torch.get_rng_state()))
        if cuda_rng is not None:
            self.assertTrue(torch.equal(cuda_rng, torch.cuda.get_rng_state(device)))
        native_state = native.state_dict()
        variant_state = model.state_dict()
        first_layer = {"feat_embeddings.0.weight", "feat_embeddings.0.bias"}
        self.assertEqual(set(native_state), set(variant_state))
        for name in native_state.keys() - first_layer:
            self.assertTrue(torch.equal(native_state[name], variant_state[name]), name)
        self.assertEqual(variant_state["feat_embeddings.0.weight"].shape, (400, 16))
        self.assertFalse(variant_state["feat_embeddings.0.weight"].count_nonzero())

        model.eval()
        with torch.no_grad():
            prediction = model(torch.zeros(2, 8, device=device), torch.ones(2, 17, device=device))[3]
        self.assertEqual(prediction.shape, (2, 8))
        self.assertTrue(torch.isfinite(prediction).all())


if __name__ == "__main__":
    unittest.main()
