"""CarPS-Transcriptome alpha model and its fixed KPGT default."""

from operator import index

import torch
from torch import nn

from .measured_model import MeasuredAdapter, configure_objective, objective_identity
from .measured_roster import MEASURED_SHARED
from tsg_axis.model import WIDTH as KPGT_WIDTH

PROJECTION_WIDTH = 400


class AlphaAdapter(MeasuredAdapter):
    """Construct the shared signed Δ model configured with the alpha objective.

    The default constructor uses the original KPGT initialization and forward.
    An experiment may supply a native chemical width and a host linear-layer
    factory accepting (input width, output width).
    """

    def __init__(self, source, transform, *, chemical_width=KPGT_WIDTH, projection_factory=None):
        super().__init__(source, transform, MEASURED_SHARED)
        self.chemical_width = index(chemical_width)
        if isinstance(chemical_width, bool) or self.chemical_width < 1:
            raise ValueError("chemical input width must be a positive integer")
        if projection_factory is None and self.chemical_width != KPGT_WIDTH:
            projection_factory = nn.Linear
        self.projection_factory = projection_factory

    def network(self):
        original = super().network()
        if self.projection_factory is None:
            return original
        width = self.chemical_width

        class ChemicalInputNetwork(original):
            def forward(self, x1, features):
                if features.ndim != 2 or features.shape[1] != width + 1:
                    raise ValueError(f"expected {width} chemical values plus one standardized dose")
                z1, mu1, logvar1 = self.encode_x1(x1)
                x1_rec = self.decode_x1(z1)
                drug_dose = torch.cat((self.feat_embeddings(features[:, :width]),
                                      features[:, width:].repeat(1, 32)), dim=1)
                mu_pred, logvar_pred = self.predictive_distribution(z1, drug_dose)
                z_pred = self.sample_latent(mu_pred, logvar_pred)
                predicted = x1 + self.prediction_decoder(z_pred)
                return x1_rec, mu1, logvar1, predicted, mu_pred, logvar_pred, z_pred

        return ChemicalInputNetwork

    def build(self, seed, device, *, n_genes=5000):
        model = super().build(seed, device, n_genes=n_genes)
        if self.projection_factory is not None:
            # Keep the RNG state inherited from the native constructor for latent sampling.
            cuda_devices = [self.device] if self.device.type == "cuda" else []
            with torch.random.fork_rng(devices=cuda_devices), torch.device("cpu"):
                projection = self.projection_factory(self.chemical_width, PROJECTION_WIDTH)
            if (not isinstance(projection, nn.Linear)
                    or projection.weight.device.type != "cpu"
                    or (projection.in_features, projection.out_features) != (self.chemical_width, PROJECTION_WIDTH)):
                raise ValueError("first chemical projection must be a host linear layer with the declared shape")
            model.feat_embeddings[0] = projection.to(device=self.device, dtype=torch.float32)
        configure_objective(model, MEASURED_SHARED)
        return model


def objective():
    """Return the inherited alpha objective weights and target identity."""
    return objective_identity(MEASURED_SHARED)
