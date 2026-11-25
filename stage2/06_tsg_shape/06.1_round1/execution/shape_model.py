"""TranSiGen expression architectures registered in ADR 0094."""

from itertools import product
from pathlib import Path
import sys

STEP = Path(__file__).resolve().parents[1]
STAGE = STEP.parents[1]
LOSS = STAGE / "05_tsg_loss"
sys.path[:0] = [str(LOSS / "05.2_mse_l1/execution"), str(LOSS / "05.1_loss_metrics/execution")]

from run_loss import axis_run

import torch

SHAPES = {
    "widen_z-100": ((1200,), (800,), 100),
    "widen_z-110": ((1200,), (800,), 110),
    "widen_z-120": ((1200,), (800,), 120),
    "widen_z-130": ((1200,), (800,), 130),
    "widen_enc": ((6000,), (800,), 100),
    "deepen_enc": ((1200, 1200), (800,), 100),
    "deepen_dec": ((1200,), (800, 800), 100),
    "deepen-both": ((1200, 1200), (800, 800), 100),
}
CONTROL = "widen_z-100"
SMOKE_SHAPES = tuple(name for name in SHAPES if name not in ("widen_z-110", "widen_z-120"))


def coordinates(smoke=False):
    return tuple(product(SMOKE_SHAPES if smoke else SHAPES, (0,) if smoke else range(3)))


def dimensions(shape):
    encoder, decoder, latent = SHAPES[shape]
    return dict(encoder=list(encoder), decoder=list(decoder), n_latent=latent,
                predictive_head_dim=latent + 400 + 32)


class ShapeAdapter(axis_run.kpgt.Adapter):
    def __init__(self, source, transform, shape):
        super().__init__(source, transform)
        self.shape = shape

    def network(self):
        dims = dimensions(self.shape)
        native = self.native
        conditioned = super().network()
        if self.shape == CONTROL:
            return conditioned

        class Shaped(conditioned):
            def __init__(self, **kwargs):
                kwargs.update(n_en_hidden=dims["encoder"], n_de_hidden=dims["decoder"],
                              n_latent=dims["n_latent"])
                # Preserve native initialization before replacing the two dose-conditioned heads.
                native.__init__(self, **kwargs)
                self.mu_z2Fz1 = torch.nn.Sequential(
                    torch.nn.Linear(dims["predictive_head_dim"], self.n_latent))
                self.logvar_z2Fz1 = torch.nn.Sequential(
                    torch.nn.Linear(dims["predictive_head_dim"], self.n_latent))

        return Shaped
