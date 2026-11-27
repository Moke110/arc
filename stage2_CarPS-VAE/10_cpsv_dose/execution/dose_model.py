"""Apply the Step 10 dose interfaces to the released CarPS-VAE backbone."""

from hashlib import sha256
import json
import sys

import numpy as np
import torch
from torch import nn

from dose_protocol import REFERENCE, STAGE, TRANSIGEN, representation

sys.path.insert(0, str(STAGE / "08_cpsv_chem-emb/execution"))
from chem_runs import AlphaAdapter, PROJECTION_WIDTH, WIDTHS, axis_run
from carps_vae_alpha import objective as alpha_objective

DOSE_WIDTH = 32
HIDDEN_WIDTH = 64
UNIMOL_WIDTH = WIDTHS["UniMol"]


def model_features(arm, chemistry, dose_um):
    """Transport frozen chemistry and, where used, a scalar dose to the model."""
    width = WIDTHS[representation(arm)]
    if chemistry.ndim != 2 or chemistry.shape[1] != width or chemistry.dtype != np.float32:
        raise ValueError(f"{arm} requires float32 chemistry vectors of width {width}")
    if arm in (REFERENCE, TRANSIGEN):
        return chemistry
    doses = axis_run.kpgt.positive_doses(dose_um)
    if len(doses) != len(chemistry):
        raise ValueError("chemical vectors and dose values must have the same row order")
    scalar = np.log10(1.0 + doses) if arm.endswith("scaling") else np.log10(doses)
    return np.concatenate((chemistry, scalar.astype(np.float32)[:, None]), axis=1)


def replaced_linear(layer, in_features):
    """Initialize added columns and preserve weights and bias on shared inputs."""
    replacement = nn.Linear(in_features, layer.out_features, device="cpu")
    shared = min(in_features, layer.in_features)
    with torch.no_grad():
        replacement.weight[:, :shared].copy_(layer.weight.detach().cpu()[:, :shared])
        replacement.bias.copy_(layer.bias.detach().cpu())
    return replacement


def configuration(arm):
    """Describe the registered dimensions and objective before fitting."""
    chemistry = representation(arm)
    width = WIDTHS[chemistry]
    output_width = {"Xmol-FiLM": 2 * UNIMOL_WIDTH, "Zmol-FiLM": 2 * PROJECTION_WIDTH}.get(arm)
    if arm.endswith("-MLP"):
        output_width = DOSE_WIDTH
    transform = (None if arm in (REFERENCE, TRANSIGEN) else
                 "log10(1 + dose / 1 µM)" if arm.endswith("scaling") else "log10(dose / 1 µM)")
    return {
        "architecture": "XX-XX" if arm == TRANSIGEN else "XΔ-XΔ",
        "molecule_feature": chemistry, "features_dim": width, "projection_width": PROJECTION_WIDTH,
        "projection_input_width": width + (DOSE_WIDTH if arm.startswith("Xdose") else 0),
        "predictive_head_dim": 100 + PROJECTION_WIDTH + (DOSE_WIDTH if arm.startswith("Zdose") else 0),
        "conditioner": [1, HIDDEN_WIDTH, output_width] if output_width is not None else None,
        "film_final_initialization": "zeros" if arm.endswith("FiLM") else None,
        "transform": transform, "n_genes": 5000, "n_latent": 100, "encoder": [1200], "decoder": [800],
        "dropout": 0.1, "compute_dtype": "float32", "batch_size": 64, "learning_rate": 1e-4,
        "weight_decay": 1e-5, "control_input": "measured same-plate DMSO in every partition",
    }


def objective(arm):
    if arm == TRANSIGEN:
        return {"name": "author TranSiGen", "beta": 0.1, "control_mse_weight": 1.0,
                "treated_mse_weight": 1.0, "delta_mse_weight": 1.0, "delta_l1_weight": 0.0}
    return alpha_objective()


class DoseAdapter(AlphaAdapter):
    """Retain alpha's base construction and alter only the arm's dose interface."""

    def __init__(self, source, arm):
        # Alpha's constructor accepts a dict for its legacy dose view. Step 10
        # supplies its features directly at the new study's execution boundary.
        super().__init__(source, {}, chemical_width=WIDTHS[representation(arm)])
        self.arm = arm

    def network(self):
        if self.arm == TRANSIGEN:
            return self.native
        original = super().network()
        arm = self.arm

        class DoseNetwork(original):
            def condition(self, vector, dose):
                if arm.endswith("scaling"):
                    return vector * dose
                if arm.endswith("FiLM"):
                    gamma, beta = self.conditioner(dose).chunk(2, dim=1)
                    return (1 + gamma) * vector + beta
                embedding = self.conditioner(dose) if arm.endswith("-MLP") else dose.repeat(1, DOSE_WIDTH)
                return torch.cat((vector, embedding), dim=1)

            def forward(self, x1, features):
                expected = UNIMOL_WIDTH + (arm != REFERENCE)
                if features.ndim != 2 or features.shape[1] != expected:
                    raise ValueError(f"{arm} requires {expected} feature columns")
                z1, mu1, logvar1 = self.encode_x1(x1)
                x1_rec = self.decode_x1(z1)
                chemistry, dose = features[:, :UNIMOL_WIDTH], features[:, UNIMOL_WIDTH:]
                if arm.startswith("X"):
                    chemistry = self.condition(chemistry, dose)
                chemical_latent = self.feat_embeddings(chemistry)
                if arm.startswith("Z"):
                    chemical_latent = self.condition(chemical_latent, dose)
                mu, logvar = self.predictive_distribution(z1, chemical_latent)
                z_pred = self.sample_latent(mu, logvar)
                predicted = x1 + self.prediction_decoder(z_pred)
                return x1_rec, mu1, logvar1, predicted, mu, logvar, z_pred

        return DoseNetwork

    def build(self, seed, device, *, n_genes=5000):
        if self.arm == TRANSIGEN:
            return axis_run.kpgt.Adapter.build(self, seed, device, n_genes=n_genes)
        model = super().build(seed, device, n_genes=n_genes)
        dims = configuration(self.arm)
        accelerators = [self.device] if self.device.type == "cuda" else []
        with torch.random.fork_rng(devices=accelerators):
            if self.arm.startswith("Xdose"):
                model.feat_embeddings[0] = replaced_linear(
                    model.feat_embeddings[0], dims["projection_input_width"]).to(self.device)
            for head in (model.mu_z2Fz1, model.logvar_z2Fz1):
                head[0] = replaced_linear(head[0], dims["predictive_head_dim"]).to(self.device)
        with torch.random.fork_rng(devices=accelerators):
            if dims["conditioner"] is not None:
                input_width, hidden, output_width = dims["conditioner"]
                conditioner = nn.Sequential(nn.Linear(input_width, hidden), nn.ReLU(),
                                            nn.Linear(hidden, output_width))
                if self.arm.endswith("FiLM"):
                    nn.init.zeros_(conditioner[-1].weight)
                    nn.init.zeros_(conditioner[-1].bias)
                model.conditioner = conditioner.to(self.device)
        return model


def streams(device):
    """Capture training CPU and CUDA RNG states without consuming draws."""
    cuda = torch.cuda.get_rng_state(device) if device.type == "cuda" else None
    return {
        "torch_cpu_rng_sha256": sha256(torch.get_rng_state().numpy().tobytes()).hexdigest(),
        "torch_cuda_rng_sha256": sha256(cuda.numpy().tobytes()).hexdigest() if cuda is not None else None,
    }


def initialization(model, arm):
    """Record per-tensor identity to compare same-shaped parameters across arms."""
    dims = configuration(arm)
    projection = model.feat_embeddings[0]
    heads = (model.mu_z2Fz1[0], model.logvar_z2Fz1[0])
    conditioner = getattr(model, "conditioner", None)
    if ((projection.in_features, projection.out_features) !=
            (dims["projection_input_width"], dims["projection_width"])
            or any(head.in_features != dims["predictive_head_dim"] for head in heads)
            or (conditioner is None) != (dims["conditioner"] is None)):
        raise ValueError(f"{arm} model does not match its registered projection or head")
    if conditioner is not None:
        widths = [conditioner[0].in_features, conditioner[0].out_features, conditioner[-1].out_features]
        if widths != dims["conditioner"] or conditioner[-1].in_features != widths[1]:
            raise ValueError(f"{arm} has an incorrect dose conditioner")
        if arm.endswith("FiLM") and any(torch.count_nonzero(value) for value in conditioner[-1].parameters()):
            raise ValueError(f"{arm} FiLM does not start with the identity transform")
    return {
        "parameters": {
            name: {"shape": list(tensor.shape), "sha256": sha256(tensor.detach().cpu().numpy().tobytes()).hexdigest()}
            for name, tensor in sorted(model.named_parameters())
        },
        **streams(next(model.parameters()).device),
    }


def epoch_streams(device, numpy_rng):
    return {
        **streams(device),
        "numpy_rng_sha256": sha256(json.dumps(numpy_rng.bit_generator.state, sort_keys=True).encode()).hexdigest(),
    }
