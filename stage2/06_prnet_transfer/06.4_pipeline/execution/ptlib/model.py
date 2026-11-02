"""PRnet ``PGM`` — the perturbation-conditioned generative model used by both stages.

**Same-source copy — do not diverge.**  This is the in-file re-implementation the
frozen steps already use (the upstream ``pert_model/PRnet/source/`` tree is never
imported, per ADR 0030):

* ``stage2/05_prnet_pretrain/execution/train_prnet.py`` (``PAdaptor`` /
  ``PEncoder`` / ``PDecoder`` / ``PGM`` / ``weight_init``), itself a line-for-line
  mirror of ``pert_model/PRnet/source/models/PRnet.py``;
* ``stage2/03_models_cmbench/03.1_cold_drug_18k/execution/script/adapters/prnet/``
  for the training convention.

Only difference from 05: layer widths are read from the resolved structure, so
the #10 capacity ablation (state latent 128, hidden 256) and #11 reuse this file
without edits.  The default resolved structure is the audited recipe:
adaptor ``1024 -> 128 -> 64``, encoder ``[18164 + 64] -> 128 -> 64`` (mean only),
decoder ``[64 + 64 + 10] -> 128 -> 2 x 18164`` with ``relu`` as the mean half and
``softplus`` as the variance half.

Convention: ``PDecoder.forward`` returns ``(mean, var)`` (a tuple).  Never assume
the author trainer's concatenated-tensor convention; call ``PGM.forward`` and
unpack.
"""

from __future__ import annotations

from typing import Any, Mapping

import torch
import torch.nn as nn
from torch.nn import functional as F

__all__ = [
    "PAdaptor",
    "PEncoder",
    "PDecoder",
    "PGM",
    "weight_init",
    "build_model",
    "random_init",
    "check_architecture",
    "parameter_count",
    "gaussian_nll",
]


def _upstream_mlp(layer_sizes: list[int], dropout_rate: float) -> nn.Sequential | None:
    """The author's MLP builder: first layer bias-free, later layers with BN+dropout."""
    if len(layer_sizes) <= 1:
        return None
    layers = nn.Sequential()
    for i, (in_size, out_size) in enumerate(zip(layer_sizes[:-1], layer_sizes[1:])):
        if i == 0:
            layers.add_module(f"L{i}", nn.Linear(in_size, out_size, bias=False))
        else:
            layers.add_module(f"L{i}", nn.Linear(in_size, out_size))
            layers.add_module(f"N{i}", nn.BatchNorm1d(out_size))
            layers.add_module(f"A{i}", nn.LeakyReLU(negative_slope=0.3))
            layers.add_module(f"D{i}", nn.Dropout(p=dropout_rate))
    return layers


class PAdaptor(nn.Module):
    """Mirror of ``PRnet.PAdaptor`` — drug rFCFP4 ``drug_dim -> adaptor[-1] -> c_dim``."""

    def __init__(self, layer_sizes: list[int], comb_dimension: int, dropout_rate: float) -> None:
        super().__init__()
        self.FC = _upstream_mlp(layer_sizes, dropout_rate)
        self.comb_encoder = nn.Linear(layer_sizes[-1], comb_dimension)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if self.FC is not None:
            x = self.FC(x)
        return self.comb_encoder(x)


class PEncoder(nn.Module):
    """Mirror of ``PRnet.PEncoder`` — ``[x + c] -> hidden -> z`` (mean only)."""

    def __init__(self, layer_sizes: list[int], z_dimension: int, dropout_rate: float) -> None:
        super().__init__()
        self.FC = _upstream_mlp(layer_sizes, dropout_rate)
        self.mean_encoder = nn.Linear(layer_sizes[-1], z_dimension)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if self.FC is not None:
            x = self.FC(x)
        return self.mean_encoder(x)


class PDecoder(nn.Module):
    """Mirror of ``PRnet.PDecoder`` — ``[z + c + n] -> hidden -> 2 * x_dim``.

    ``relu(head[:, :dim])`` is the non-negative mean; ``softplus(head[:, dim:])``
    is the strictly positive variance (the author trainer applies the same
    softplus, ``PRnetTrainer.py:190``).  Returns the tuple convention.
    """

    def __init__(
        self,
        z_dimension: int,
        layer_sizes: list[int],
        x_dimension: int,
        dropout_rate: float,
    ) -> None:
        super().__init__()
        layer_sizes = [z_dimension] + layer_sizes
        self.FirstL = nn.Sequential(
            nn.Linear(layer_sizes[0], layer_sizes[1], bias=False),
            nn.BatchNorm1d(layer_sizes[1]),
            nn.LeakyReLU(negative_slope=0.3),
            nn.Dropout(p=dropout_rate),
        )
        # The author's ``HiddenL`` is empty for a one-hidden-layer config
        # (``i + 3 < len(layer_sizes)`` never fires); keep the attribute for state
        # dict parity and a no-op forward branch.
        self.HiddenL = None
        self.recon_decoder = nn.Linear(layer_sizes[-2], layer_sizes[-1])

    def forward(self, z: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        x = self.FirstL(z)
        if self.HiddenL is not None:  # pragma: no cover - mirrors the author branch
            x = self.HiddenL(x)
        head = self.recon_decoder(x)
        dim = head.size(1) // 2
        return F.relu(head[:, :dim]), F.softplus(head[:, dim:])


class PGM(nn.Module):
    """Mirror of ``PRnet.PGM`` — the perturbation-conditioned generative model."""

    def __init__(
        self,
        x_dim: int,
        c_dim: int,
        n_dim: int,
        hidden_layer_sizes: list[int],
        z_dimension: int,
        adaptor_layer_sizes: list[int],
        comb_adapt_dim: int,
        dr_rate: float,
    ) -> None:
        super().__init__()
        self.x_dim, self.c_dim, self.n_dim, self.z_dim = x_dim, c_dim, n_dim, z_dimension

        encoder_layer_sizes = [x_dim + c_dim, *hidden_layer_sizes]
        decoder_layer_sizes = [*reversed(hidden_layer_sizes), x_dim * 2]
        adaptor_layer_sizes = [comb_adapt_dim, *adaptor_layer_sizes]

        self.CombAdaptor = PAdaptor(adaptor_layer_sizes, c_dim, dr_rate)
        self.encoder = PEncoder(encoder_layer_sizes, z_dimension, dr_rate)
        self.decoder = PDecoder(
            z_dimension + c_dim + n_dim, decoder_layer_sizes, x_dim, dr_rate
        )

    def forward(
        self, x: torch.Tensor, c: torch.Tensor, n: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor]:
        c = self.CombAdaptor(c)
        z = self.encoder(torch.cat((x, c), dim=1))
        return self.decoder(torch.cat((z, c, n), dim=1))


def weight_init(module: nn.Module) -> None:
    """The author's ``weight_init`` (``PRnetTrainer.py:491-504``)."""
    if isinstance(module, nn.Linear):
        module.weight.data.normal_(0, 0.02)
        if module.bias is not None:
            module.bias.data.zero_()
    elif isinstance(module, nn.BatchNorm1d):
        module.weight.data.normal_(1, 0.02)
        if module.bias is not None:
            module.bias.data.zero_()


def build_model(struct: Mapping[str, Any]) -> PGM:
    """Construct an un-initialised :class:`PGM` from a canonical structure."""
    return PGM(
        x_dim=int(struct["n_genes"]),
        c_dim=int(struct["c_dim"]),
        n_dim=int(struct["n_dim"]),
        hidden_layer_sizes=[int(w) for w in struct["hidden"]],
        z_dimension=int(struct["z_dim"]),
        adaptor_layer_sizes=[int(w) for w in struct["adaptor"]],
        comb_adapt_dim=int(struct["drug_dim"]),
        dr_rate=float(struct["dropout"]),
    )


def random_init(model: PGM, seed: int) -> PGM:
    """The documented author random initialisation, seeded for reproducibility.

    ``normal_(0, 0.02)`` weights, ``normal_(1, 0.02)`` BatchNorm scales, zero
    biases.  The run seed makes the initialisation a recorded property of the run
    (so a random-init candidate can be reproduced exactly).
    """
    torch.manual_seed(int(seed))
    model.apply(weight_init)
    return model


def parameter_count(model: PGM) -> int:
    return int(sum(p.numel() for p in model.parameters()))


def check_architecture(model: PGM, struct: Mapping[str, Any]) -> None:
    """Assert the resolved PRnet shapes before any GPU time is spent.

    Shape arithmetic, not a hard-coded 18,164: the same check guards the capacity
    ablations and the synthetic-fixture tests.
    """
    n_genes = int(struct["n_genes"])
    hidden = [int(w) for w in struct["hidden"]]
    z_dim = int(struct["z_dim"])
    c_dim = int(struct["c_dim"])
    n_dim = int(struct["n_dim"])
    drug_dim = int(struct["drug_dim"])
    adaptor = [int(w) for w in struct["adaptor"]]
    expected = {
        "adaptor in": (model.CombAdaptor.FC[0].weight.shape, (adaptor[0], drug_dim)),
        "adaptor out": (model.CombAdaptor.comb_encoder.weight.shape, (c_dim, adaptor[-1])),
        "encoder in": (model.encoder.FC[0].weight.shape, (hidden[0], n_genes + c_dim)),
        "encoder out": (model.encoder.mean_encoder.weight.shape, (z_dim, hidden[0])),
        "decoder in": (
            model.decoder.FirstL[0].weight.shape,
            (hidden[0], z_dim + c_dim + n_dim),
        ),
        "decoder out": (model.decoder.recon_decoder.weight.shape, (2 * n_genes, hidden[0])),
    }
    for name, (got, want) in expected.items():
        if tuple(got) != tuple(want):
            raise SystemExit(
                f"architecture mismatch at {name}: got {tuple(got)}, want {tuple(want)}"
            )
    model.eval()
    device = next(model.parameters()).device
    with torch.no_grad():
        mean, var = model(
            torch.zeros(2, n_genes, device=device),
            torch.zeros(2, drug_dim, device=device),
            torch.zeros(2, n_dim, device=device),
        )
    if tuple(mean.shape) != (2, n_genes) or tuple(var.shape) != (2, n_genes):
        raise SystemExit(f"forward shape mismatch: mean {tuple(mean.shape)} var {tuple(var.shape)}")
    if not bool(torch.all(var > 0)):
        raise SystemExit("decoder variance head is not strictly positive")
    model.train()


def gaussian_nll(
    model: PGM,
    control: torch.Tensor,
    label: torch.Tensor,
    noise: torch.Tensor,
    target: torch.Tensor,
) -> torch.Tensor:
    """``GUSS`` — Gaussian negative log-likelihood; **may be negative, never clamped**."""
    mean, var = model(control, label, noise)
    return F.gaussian_nll_loss(mean, target, var, reduction="mean")
