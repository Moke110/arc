"""Compose the signed readouts with the normalized Round 3 FiLM Gaussian head."""

from pathlib import Path
import sys

from torch import nn

STEP = Path(__file__).resolve().parents[2] / "util/training/legacy_transigen/modules/composition"
ROUND_TWO = STEP.with_name("prediction_path")
sys.path.insert(0, str(ROUND_TWO / "execution"))

from .architecture_model import ArchitectureAdapter, FiLMPrediction
from .architecture_model import configure_objective as configure_round_two
from .architecture_model import objective_identity as round_two_objective
from .composition_roster import COMBINATIONS, NORMALIZED_FILM_SHAPES, ROUND_TWO_SHAPES


class NormalizedFiLMPrediction(nn.Module):
    """Apply parameter-free LayerNorm to the modulated control latent."""

    def __init__(self, head):
        super().__init__()
        self.conditioner = head.conditioner
        self.projection = head.projection
        self.normalization = nn.LayerNorm(100, elementwise_affine=False)

    def forward(self, control, drug_dose):
        scale, shift = self.conditioner(drug_dose).chunk(2, dim=1)
        return self.projection(self.normalization((1 + scale) * control + shift))


class CompositionAdapter(ArchitectureAdapter):
    def __init__(self, source, transform, shape):
        super().__init__(source, transform, ROUND_TWO_SHAPES[shape])
        self.composition = shape

    def network(self):
        network = super().network()
        if self.composition not in NORMALIZED_FILM_SHAPES:
            return network
        combined = self.composition in COMBINATIONS

        class ComposedArchitecture(network):
            def __init__(self, **kwargs):
                super().__init__(**kwargs)
                if combined:
                    del self.mu_z2Fz1, self.logvar_z2Fz1
                    self.prediction_head = FiLMPrediction()
                self.prediction_head = NormalizedFiLMPrediction(self.prediction_head)

            def predictive_distribution(self, z1, drug_dose):
                return self.prediction_head(z1, drug_dose).chunk(2, dim=1)

        return ComposedArchitecture


def configure_objective(model, shape):
    configure_round_two(model, ROUND_TWO_SHAPES[shape])


def objective_identity(shape):
    return round_two_objective(ROUND_TWO_SHAPES[shape])
