"""Compose the accepted Round 2 signed readouts with its FiLM Gaussian head."""

from pathlib import Path
import sys

STEP = Path(__file__).resolve().parents[1]
ROUND_TWO = STEP.with_name("06.2_round2")
sys.path.insert(0, str(ROUND_TWO / "execution"))

from architecture_model import ArchitectureAdapter, FiLMPrediction, axis_run
from architecture_model import configure_objective as configure_round_two
from architecture_model import objective_identity as round_two_objective
from composition_roster import COMBINATIONS, ROUND_TWO_SHAPES


class CompositionAdapter(ArchitectureAdapter):
    def __init__(self, source, transform, shape):
        super().__init__(source, transform, ROUND_TWO_SHAPES[shape])
        self.composition = shape

    def network(self):
        network = super().network()
        if self.composition not in COMBINATIONS:
            return network

        class ComposedArchitecture(network):
            def __init__(self, **kwargs):
                super().__init__(**kwargs)
                del self.mu_z2Fz1, self.logvar_z2Fz1
                self.prediction_head = FiLMPrediction()

            def predictive_distribution(self, z1, drug_dose):
                return self.prediction_head(z1, drug_dose).chunk(2, dim=1)

        return ComposedArchitecture


def configure_objective(model, shape):
    configure_round_two(model, ROUND_TWO_SHAPES[shape])


def objective_identity(shape):
    return round_two_objective(ROUND_TWO_SHAPES[shape])
