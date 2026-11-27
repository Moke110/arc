"""TranSiGen expression architectures registered in ADRs 0094 and 0095."""

from pathlib import Path
import sys

from . import native

import torch

from .shape_roster import CONTROL, dimensions


class ShapeAdapter(native.Adapter):
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
