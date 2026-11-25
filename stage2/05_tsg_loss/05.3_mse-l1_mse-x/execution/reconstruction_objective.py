"""Weight the two matched-control expression reconstruction errors together."""

from pathlib import Path
import sys

PREDECESSOR_EXECUTION = Path(__file__).resolve().parents[2] / "05.2_mse_l1/execution"
sys.path.append(str(PREDECESSOR_EXECUTION))

from additive_objective import AdditiveObjective
from sweep_protocol import ALPHA_LABELS, DELTA_L1_WEIGHT


class ReconstructionObjective:
    """Keep the Step 05.2 Delta and KL terms while varying reconstruction."""

    def __init__(self, alpha):
        if alpha not in ALPHA_LABELS:
            raise ValueError(f"unsupported shared reconstruction MSE weight: {alpha}")
        self.alpha = alpha
        self.predecessor = AdditiveObjective(DELTA_L1_WEIGHT)

    def identity(self):
        return {**self.predecessor.identity(), "alpha": ALPHA_LABELS[self.alpha],
                "reconstruction_weight": float(self.alpha)}

    def __call__(self, native_loss, *inputs):
        predecessor = self.predecessor(native_loss, *inputs)
        if self.alpha == 1:
            return predecessor

        _, control_mse, treatment_mse, delta_mse, kl1, kl2, kl_prediction, delta_l1 = predecessor
        total = 0.1 * (kl1 + kl2 + kl_prediction) + delta_mse + float(DELTA_L1_WEIGHT) * delta_l1
        if self.alpha:
            total = total + float(self.alpha) * (control_mse + treatment_mse)
        return total, *predecessor[1:]
