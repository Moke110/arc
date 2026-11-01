#!/usr/bin/env python
"""Stage 2 · 05 — PRnet pretrain bring-up on the Tahoe-100M pseudobulk view (GPU).

This is the *pipeline bring-up* the user asked for: train PRnet's perturbation-
conditioned generative model on the model-facing Tahoe view until a loss-epoch
curve exists.  No CM finetune, no cold-drug evaluation, no benchmark comparison.

SAME-SOURCE COPY — DO NOT DIVERGE
---------------------------------
The model and the training loop are re-implemented here (the frozen upstream tree
is **not** imported) and are deliberately line-for-line equivalents of:

* ``pert_model/PRnet/source/models/PRnet.py``      — ``PGM`` (l.70-141),
  ``PEncoder`` (l.144-181), ``PDecoder`` (l.183-233), ``PAdaptor`` (l.238-271).
  Resolved sizes: adaptor ``1024 -> 128 -> 64``; encoder ``[18164+64] -> 128 -> 64``
  (mean only); decoder ``[64+64+10] -> 128 -> 2*18164`` with ``relu(head[:, :dim])``
  as the mean and ``softplus(head[:, dim:])`` as the variance (the author's trainer
  applies that softplus itself, ``PRnetTrainer.py:190``; folding it into the decoder
  is numerically identical and is what spec §2 states).
* ``pert_model/PRnet/source/trainer/PRnetTrainer.py:156-258`` — ``Adam(lr=1e-3,
  weight_decay=1e-8)``, ``ReduceLROnPlateau(mode='min', factor=0.5, patience=10,
  min_lr=1e-8)`` stepped on the **validation** loss, ``GaussianNLLLoss`` ("GUSS"),
  early stop after 20 epochs without improvement, ``noise_dim=10``.
* ``pert_model/PRnet/source/data/_utils.py:44-68`` — ``Drug_dose_encoder``:
  ``GetMorganFingerprintAsBitVect(mol, 2, useFeatures=True, nBits=1024)`` (rFCFP4)
  multiplied by ``log10(dose_uM + 1)``.
* ``stage2/03_models_cmbench/03.1_cold_drug_18k/execution/script/adapters/prnet/run.py:139-279``
  — the frozen 18k benchmark loop; this step keeps its optimiser / scheduler / batch /
  epoch / early-stop / seed values so the two runs are comparable.

Contract: ``.scratch/stage2-05-prnet-pretrain/spec.md`` §2.

Loss note: ``Gaussian NLL`` is a negative log-likelihood and **can be negative**.
It is never clamped, ``abs()``-ed or otherwise "fixed" here.

Usage (from the repository root)::

    PY=/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python
    $PY stage2/05_prnet_pretrain/execution/train_prnet.py --check
    $PY stage2/05_prnet_pretrain/execution/train_prnet.py          # full 100-epoch run
    $PY stage2/05_prnet_pretrain/execution/train_prnet.py --smoke   # 2 epochs, fixture
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import platform
import resource
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
import torch.nn as nn
from torch.nn import functional as F

# --------------------------------------------------------------------------- #
# Frozen contract constants (spec §2; identical to the 03.1 PRnet adapter)
# --------------------------------------------------------------------------- #
N_GENES = 18164
DRUG_DIM = 1024
COMB_DIM = 64
Z_DIM = 64
NOISE_DIM = 10
HIDDEN = (128,)
ADAPTOR_HIDDEN = (128,)
DROP_RATE = 0.05

BATCH_SIZE = 512
MAX_EPOCHS = 100
LR = 1e-3
WEIGHT_DECAY = 1e-8
SCHED_FACTOR = 0.5
SCHED_PATIENCE = 10
MIN_LR = 1e-8
EARLY_STOP_PATIENCE = 20
MIN_EPOCHS = 20
SEED = 20261021

VALID_PLATES = ("plate13", "plate14")
SMOKE_EPOCHS = 2
SMOKE_TRAIN_ROWS = 3000
SMOKE_VALID_ROWS = 1000

VIEW_DEFAULT = "stage2/04_tahoe_data/output/pretrain_view/tahoe_cm_schema.h5ad"
OUTPUT_DEFAULT = "stage2/05_prnet_pretrain/output"

#: candidate ``uns`` keys holding the ``control_group`` labels when
#: ``uns['control_mean']`` is stored as a bare array.
CONTROL_LABEL_KEYS = (
    "control_mean_index",
    "control_groups",
    "control_group_index",
    "control_mean_groups",
)

log = logging.getLogger("prnet_pretrain")


# --------------------------------------------------------------------------- #
# Paths / logging
# --------------------------------------------------------------------------- #
def find_repo_root(start: Path) -> Path:
    """Walk up from *start* to the repository root (the directory holding ``.git``).

    Never counts directory levels: the Stage 2 layout has already moved once
    (ADR 0030) and a hard-coded depth fails silently when it moves again.
    """
    for candidate in (start, *start.parents):
        if (candidate / ".git").exists():
            return candidate
    raise SystemExit(f"cannot find the repository root: no .git at or above {start}")


def resolve(path: Path, root: Path) -> Path:
    return path if path.is_absolute() else (root / path)


def setup_logging(run_log: Path | None) -> None:
    """Log to stdout, and to ``run.log`` in **overwrite** mode (reruns never append)."""
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)-5s %(message)s", datefmt="%H:%M:%S")
    stream = logging.StreamHandler(sys.stdout)
    stream.setFormatter(fmt)
    log.addHandler(stream)
    if run_log is not None:
        file_handler = logging.FileHandler(run_log, mode="w", encoding="utf-8")
        file_handler.setFormatter(fmt)
        log.addHandler(file_handler)


def peak_rss_mb() -> float:
    # ru_maxrss is KiB on Linux.
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


# --------------------------------------------------------------------------- #
# Model — mirrors pert_model/PRnet/source/models/PRnet.py (see module docstring)
# --------------------------------------------------------------------------- #
def _upstream_mlp(layer_sizes: list[int], dropout_rate: float) -> nn.Sequential | None:
    """The author's ``for i, (in_size, out_size) in enumerate(zip(...))`` builder.

    Mirrors ``PAdaptor.__init__``/``PEncoder.__init__``: the first layer has no bias,
    every later layer is ``Linear -> BatchNorm1d -> LeakyReLU(0.3) -> Dropout``.
    """
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
    """Mirror of ``PRnet.PAdaptor`` — drug rFCFP4 ``1024 -> 128 -> 64``."""

    def __init__(
        self,
        layer_sizes: list[int],
        comb_dimension: int = COMB_DIM,
        dropout_rate: float = DROP_RATE,
    ) -> None:
        super().__init__()
        self.FC = _upstream_mlp(layer_sizes, dropout_rate)
        self.comb_encoder = nn.Linear(layer_sizes[-1], comb_dimension)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if self.FC is not None:
            x = self.FC(x)
        return self.comb_encoder(x)


class PEncoder(nn.Module):
    """Mirror of ``PRnet.PEncoder`` — ``[x + comb] -> 128 -> z`` (mean only)."""

    def __init__(self, layer_sizes: list[int], z_dimension: int, dropout_rate: float) -> None:
        super().__init__()
        self.FC = _upstream_mlp(layer_sizes, dropout_rate)
        self.mean_encoder = nn.Linear(layer_sizes[-1], z_dimension)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if self.FC is not None:
            x = self.FC(x)
        return self.mean_encoder(x)


class PDecoder(nn.Module):
    """Mirror of ``PRnet.PDecoder`` — ``[z+comb+noise] -> 128 -> 2*x_dim``.

    Returns ``(mean, var)``: ``relu(head[:, :dim])`` is the mean and
    ``softplus(head[:, dim:])`` the variance (spec §2; the author's trainer applies
    the same softplus to this head at ``PRnetTrainer.py:190``).
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
        # The author builds an (empty) ``HiddenL`` here because
        # ``i + 3 < len(layer_sizes)`` never fires for [138, 128, 36288]; an empty
        # Sequential is a no-op, so None is equivalent and keeps the state_dict equal.
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
        x_dim: int = N_GENES,
        c_dim: int = COMB_DIM,
        n_dim: int = NOISE_DIM,
        hidden_layer_sizes: tuple[int, ...] = HIDDEN,
        z_dimension: int = Z_DIM,
        adaptor_layer_sizes: tuple[int, ...] = ADAPTOR_HIDDEN,
        comb_adapt_dim: int = DRUG_DIM,
        dr_rate: float = DROP_RATE,
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


def check_architecture(pgm: PGM) -> None:
    """Assert the audited PRnet dimensions of spec §2 before wasting GPU time."""
    expected = {
        "adaptor in": (pgm.CombAdaptor.FC[0].weight.shape, (HIDDEN[0], DRUG_DIM)),
        "adaptor out": (pgm.CombAdaptor.comb_encoder.weight.shape, (COMB_DIM, HIDDEN[0])),
        "encoder in": (pgm.encoder.FC[0].weight.shape, (HIDDEN[0], N_GENES + COMB_DIM)),
        "encoder out": (pgm.encoder.mean_encoder.weight.shape, (Z_DIM, HIDDEN[0])),
        "decoder in": (pgm.decoder.FirstL[0].weight.shape, (HIDDEN[0], Z_DIM + COMB_DIM + NOISE_DIM)),
        "decoder out": (pgm.decoder.recon_decoder.weight.shape, (2 * N_GENES, HIDDEN[0])),
    }
    for name, (got, want) in expected.items():
        if tuple(got) != tuple(want):
            raise SystemExit(f"architecture mismatch at {name}: got {tuple(got)}, want {tuple(want)}")
    pgm.eval()
    device = next(pgm.parameters()).device
    with torch.no_grad():
        mean, var = pgm(
            torch.zeros(2, N_GENES, device=device),
            torch.zeros(2, DRUG_DIM, device=device),
            torch.zeros(2, NOISE_DIM, device=device),
        )
    if tuple(mean.shape) != (2, N_GENES) or tuple(var.shape) != (2, N_GENES):
        raise SystemExit(f"forward shape mismatch: mean {tuple(mean.shape)} var {tuple(var.shape)}")
    if not bool(torch.all(var > 0)):
        raise SystemExit("decoder variance head is not strictly positive")
    pgm.train()


# --------------------------------------------------------------------------- #
# View loading (CSR layer only — never densify; spec §2)
# --------------------------------------------------------------------------- #
@dataclass
class ViewMeta:
    n_obs: int
    n_vars: int
    var_names: np.ndarray
    obs: pd.DataFrame
    control_values: np.ndarray  # float32 (n_groups, n_vars)
    control_labels: list[str]
    control_mean_encoding: str


def read_sparse_layer(path: Path, layer: str, n_obs: int, n_vars: int) -> sp.csr_matrix:
    """Read one CSR layer from an ``.h5ad`` without touching ``X`` or densifying."""
    import h5py

    with h5py.File(path, "r") as handle:
        if "layers" not in handle or layer not in handle["layers"]:
            raise SystemExit(
                f"{path} has no layers['{layer}']; this is not the 04 pretrain view "
                "(expected `X` raw counts + `layers['log1p_cp10k']`)"
            )
        node = handle["layers"][layer]
        if isinstance(node, h5py.Group) and {"data", "indices", "indptr"} <= set(node.keys()):
            data = node["data"][:]
            indices = node["indices"][:]
            indptr = node["indptr"][:]
            shape = tuple(int(x) for x in node.attrs.get("shape", (n_obs, n_vars)))
            matrix = sp.csr_matrix((data, indices, indptr), shape=shape)
        else:  # dense fallback: tolerated for tiny fixtures only
            dense = np.asarray(node[:], dtype=np.float32)
            if dense.nbytes > 2 * 1024**3:
                raise SystemExit(
                    f"layers['{layer}'] is dense ({dense.nbytes / 1024**3:.1f} GiB); "
                    "the contract requires CSR (deviation D0)"
                )
            matrix = sp.csr_matrix(dense)

    if matrix.shape != (n_obs, n_vars):
        raise SystemExit(f"layers['{layer}'] shape {matrix.shape} != view shape {(n_obs, n_vars)}")
    if int(matrix.indptr[0]) != 0 or int(matrix.indptr[-1]) != matrix.data.size:
        raise SystemExit(f"layers['{layer}'] indptr is inconsistent")
    if matrix.data.dtype != np.float32:
        matrix = matrix.astype(np.float32)
    return matrix


def resolve_control_mean(
    uns: dict, var_names: np.ndarray, obs_groups: np.ndarray, n_groups_expected: int
) -> tuple[np.ndarray, list[str], str]:
    """Materialise ``uns['control_mean']`` as a float32 ``(groups, genes)`` matrix.

    The production encoding is spec §1.4's ``(700, 18164)`` ``DataFrame`` indexed by
    ``control_group``.  It is written by 04 with anndata 0.12.16 (``dataframe/0.2.0``):
    anndata 0.10.8's *writer* overflows the HDF5 attribute that holds the 18,164 column
    names, but its *reader* handles the 0.2.0 encoding, which is all this step needs.
    Accepted forms, in preference order:

    1. ``DataFrame`` ``(groups, genes)`` — the production path,
    2. ``DataFrame`` ``(genes, groups)``, transposed (700 column names do fit),
    3. ``ndarray``/sparse ``(groups, genes)`` + a label array under one of
       :data:`CONTROL_LABEL_KEYS` (robustness fallback, used by the synthetic fixture
       when it is generated by an anndata too old to write the DataFrame),
    4. as (3) with the labels omitted — the index is then taken to be the sorted
       ``control_group`` values, which is what spec §1.4 prescribes.

    Every case is cross-checked against the view's ``var`` order and ``obs`` groups.
    The resolved encoding is returned so ``config.json`` / ``model_card.md`` record it.
    """
    n_vars = len(var_names)
    raw = uns.get("control_mean")
    if raw is None:
        raise SystemExit("uns['control_mean'] is missing; this is not the 04 pretrain view")

    labels: list[str] | None = None
    encoding = "unknown"
    if isinstance(raw, pd.DataFrame):
        frame = raw
        if frame.shape == (n_groups_expected, n_vars):
            encoding = "dataframe"
            if [str(c) for c in frame.columns] == [str(v) for v in var_names]:
                values = frame.to_numpy(dtype=np.float64)
            else:
                frame = frame.reindex(columns=[str(v) for v in var_names])
                if frame.shape != (n_groups_expected, n_vars):
                    raise SystemExit("uns['control_mean'] columns do not match the view var axis")
                values = frame.to_numpy(dtype=np.float64)
            labels = [str(i) for i in frame.index]
        elif frame.shape == (n_vars, n_groups_expected):
            encoding = "dataframe_transposed"
            if [str(i) for i in frame.index] != [str(v) for v in var_names]:
                raise SystemExit("transposed uns['control_mean'] rows do not match the view var axis")
            values = frame.to_numpy(dtype=np.float64).T
            labels = [str(c) for c in frame.columns]
        else:
            raise SystemExit(f"uns['control_mean'] DataFrame shape {frame.shape} is not recognised")
    else:
        values = raw.toarray() if sp.issparse(raw) else np.asarray(raw)
        values = np.asarray(values, dtype=np.float64)
        if values.ndim != 2:
            raise SystemExit(f"uns['control_mean'] is {values.ndim}-D, expected 2-D")
        if values.shape == (n_vars, n_groups_expected):
            values = values.T
        elif values.shape != (n_groups_expected, n_vars):
            raise SystemExit(
                f"uns['control_mean'] shape {values.shape} is neither "
                f"{(n_groups_expected, n_vars)} nor {(n_vars, n_groups_expected)}"
            )
        for key in CONTROL_LABEL_KEYS:
            if key in uns:
                labels = [str(x) for x in np.asarray(uns[key]).reshape(-1)]
                encoding = "ndarray+labels"
                break

    if labels is None:
        labels = sorted({str(g) for g in obs_groups})
        encoding = "ndarray+sorted_obs_order"
        log.warning(
            "uns['control_mean'] carries no label array; assuming the sorted "
            "control_group order prescribed by spec §1.4"
        )

    if len(labels) != values.shape[0]:
        raise SystemExit(
            f"uns['control_mean'] has {values.shape[0]} rows but {len(labels)} labels"
        )
    positions = pd.Index(labels).get_indexer(pd.Index([str(g) for g in obs_groups]))
    if (positions < 0).any():
        missing = sorted({str(g) for g, p in zip(obs_groups, positions) if p < 0})[:3]
        raise SystemExit(f"rows reference control_group(s) absent from control_mean: {missing}")
    log.info(
        "control_mean: %d groups x %d genes, encoding=%s, labels from %s",
        values.shape[0],
        values.shape[1],
        encoding,
        "the stored index" if "dataframe" in encoding or "labels" in encoding else "sorted obs order",
    )
    return values.astype(np.float32), labels, encoding


def load_view(view: Path, with_layer: bool) -> tuple[ViewMeta, sp.csr_matrix | None]:
    import anndata as ad

    started = time.perf_counter()
    backed = ad.read_h5ad(view, backed="r")
    try:
        obs = backed.obs.copy()
        var_names = np.asarray(backed.var_names.astype(str))
        n_obs, n_vars = int(backed.n_obs), int(backed.n_vars)
        uns = dict(backed.uns)
    finally:
        file_handle = getattr(backed, "file", None)
        if file_handle is not None:
            file_handle.close()
    log.info(
        "read view metadata: %d x %d, obs %d cols, uns %s (%.1fs, rss %.0f MiB)",
        n_obs,
        n_vars,
        obs.shape[1],
        sorted(uns)[:6],
        time.perf_counter() - started,
        peak_rss_mb(),
    )

    if n_vars != N_GENES:
        raise SystemExit(f"view has {n_vars} genes, spec §2 requires {N_GENES}")
    if "control_group" not in obs.columns:
        raise SystemExit("obs['control_group'] is missing")
    obs_groups = obs["control_group"].astype(str).to_numpy()
    n_groups_expected = int(pd.unique(obs_groups).size)
    control_values, control_labels, control_mean_encoding = resolve_control_mean(
        uns, var_names, obs_groups, n_groups_expected
    )

    layer = None
    if with_layer:
        started = time.perf_counter()
        layer = read_sparse_layer(view, "log1p_cp10k", n_obs, n_vars)
        log.info(
            "read layers['log1p_cp10k']: CSR %s nnz=%d (%.1fs, rss %.0f MiB)",
            layer.shape,
            layer.nnz,
            time.perf_counter() - started,
            peak_rss_mb(),
        )
    return ViewMeta(
        n_obs, n_vars, var_names, obs, control_values, control_labels, control_mean_encoding
    ), layer


# --------------------------------------------------------------------------- #
# Selection / split / features
# --------------------------------------------------------------------------- #
@dataclass
class Selection:
    rows: np.ndarray  # view row positions of the training samples
    fp_ix: np.ndarray  # index into the fingerprint cache
    dose_scale: np.ndarray  # log10(dose_uM + 1)
    cg_ix: np.ndarray  # index into control_values
    fingerprints: np.ndarray  # (unique SMILES, 1024) rFCFP4 cache
    train: np.ndarray  # boolean mask over `rows`
    valid: np.ndarray
    counts: dict[str, int]


def build_selection(meta: ViewMeta) -> Selection:
    obs = meta.obs
    for column in ("is_control", "has_structure", "smiles", "dose_um", "plate", "control_group"):
        if column not in obs.columns:
            raise SystemExit(f"obs['{column}'] is missing; not the 04 pretrain view")

    is_control = obs["is_control"].to_numpy(dtype=bool)
    has_structure = obs["has_structure"].to_numpy(dtype=bool)
    smiles = obs["smiles"].fillna("").astype(str).to_numpy()
    dose = obs["dose_um"].to_numpy(dtype=np.float64)
    plate = obs["plate"].astype(str).to_numpy()

    counts: dict[str, int] = {
        "n_obs": int(len(obs)),
        "n_controls": int(is_control.sum()),
        "n_perturbed": int((~is_control).sum()),
        "dropped_no_structure": int(((~is_control) & ~has_structure).sum()),
        "dropped_empty_smiles": 0,
        "dropped_bad_dose": 0,
    }

    keep = (~is_control) & has_structure
    empty = keep & (smiles == "")
    counts["dropped_empty_smiles"] = int(empty.sum())
    if counts["dropped_empty_smiles"]:
        log.warning(
            "dropping %d non-control row(s) with has_structure but an empty SMILES",
            counts["dropped_empty_smiles"],
        )
    keep &= smiles != ""
    bad_dose = ~np.isfinite(dose) | (dose < 0)
    bad = keep & bad_dose
    counts["dropped_bad_dose"] = int(bad.sum())
    if counts["dropped_bad_dose"]:
        log.warning("dropping %d non-control row(s) with a non-finite/negative dose", counts["dropped_bad_dose"])
    keep &= ~bad_dose

    rows = np.flatnonzero(keep).astype(np.int64)
    valid = np.isin(plate[rows], VALID_PLATES)
    counts["n_samples"] = int(rows.size)
    counts["n_train"] = int((~valid).sum())
    counts["n_valid"] = int(valid.sum())
    if rows.size == 0:
        raise SystemExit("no training samples survived the §2 row filter")
    if counts["n_valid"] == 0:
        raise SystemExit(
            f"no validation rows: none of the training samples sits on {VALID_PLATES}"
        )
    if counts["n_train"] == 0:
        raise SystemExit(f"no training rows: every sample sits on {VALID_PLATES}")

    uniq_smiles, fp_ix = np.unique(smiles[rows], return_inverse=True)
    fingerprints = morgan_fingerprints(uniq_smiles)
    counts["n_unique_smiles"] = int(uniq_smiles.size)
    counts["n_fingerprints_all_zero"] = int((fingerprints.sum(axis=1) == 0).sum())
    if counts["n_fingerprints_all_zero"]:
        log.warning("%d SMILES produced an all-zero rFCFP4", counts["n_fingerprints_all_zero"])

    cg_positions = pd.Index(meta.control_labels).get_indexer(
        pd.Index(obs["control_group"].astype(str).to_numpy()[rows])
    )
    if (cg_positions < 0).any():  # resolve_control_mean already checked; keep a guard
        raise SystemExit("control_group resolution failed for selected rows")

    dose_scale = np.log10(dose[rows] + 1.0).astype(np.float32)
    log.info(
        "samples: %d (%d train / %d valid) | controls %d | perturbed %d | "
        "dropped: %d no-structure, %d empty-smiles, %d bad-dose | %d unique SMILES",
        counts["n_samples"],
        counts["n_train"],
        counts["n_valid"],
        counts["n_controls"],
        counts["n_perturbed"],
        counts["dropped_no_structure"],
        counts["dropped_empty_smiles"],
        counts["dropped_bad_dose"],
        counts["n_unique_smiles"],
    )
    return Selection(
        rows=rows,
        fp_ix=fp_ix.astype(np.int32),
        dose_scale=dose_scale,
        cg_ix=cg_positions.astype(np.int32),
        fingerprints=fingerprints,
        train=~valid,
        valid=valid,
        counts=counts,
    )


def morgan_fingerprints(unique_smiles: np.ndarray, n_bits: int = DRUG_DIM) -> np.ndarray:
    """rFCFP4: ``GetMorganFingerprintAsBitVect(mol, 2, useFeatures=True, nBits=1024)``.

    One fingerprint per **unique** SMILES (≈377 in the real view), cached; the dose
    scaling is applied per row (spec §2).  Mirrors ``_utils.Drug_dose_encoder``.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem

    RDLogger.DisableLog("rdApp.*")
    table = np.zeros((len(unique_smiles), n_bits), dtype=np.float32)
    failures = 0
    for i, smiles in enumerate(unique_smiles):
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            failures += 1
            continue
        bit_string = AllChem.GetMorganFingerprintAsBitVect(
            mol, 2, useFeatures=True, nBits=n_bits
        ).ToBitString()
        if len(bit_string) != n_bits:
            raise SystemExit(f"RDKit returned {len(bit_string)} bits, expected {n_bits}")
        table[i] = np.array(list(bit_string), dtype=np.float32)
    if failures:
        log.warning("%d/%d SMILES failed to parse in RDKit (zero fingerprint)", failures, len(unique_smiles))
    return table


class BatchGatherer:
    """Builds ``(control, fingerprint, target)`` batches from CSR + cached features."""

    def __init__(
        self,
        fingerprints: np.ndarray,
        control_values: np.ndarray,
        layer: sp.csr_matrix,
        device: torch.device,
    ) -> None:
        self._fingerprints = fingerprints
        self._control_values = control_values
        self._layer = layer
        self._device = device

    def batch(self, split: SplitData, idx: np.ndarray) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        label = self._fingerprints[split.fp_ix[idx]]  # fancy indexing copies
        label *= split.dose_scale[idx, None]
        control = self._control_values[split.cg_ix[idx]]
        target = self._layer[split.rows[idx]].toarray()
        return (
            torch.from_numpy(control).to(self._device),
            torch.from_numpy(np.ascontiguousarray(label)).to(self._device),
            torch.from_numpy(np.ascontiguousarray(target)).to(self._device),
        )


@dataclass
class SplitData:
    name: str
    rows: np.ndarray
    fp_ix: np.ndarray
    dose_scale: np.ndarray
    cg_ix: np.ndarray

    def __len__(self) -> int:
        return int(self.rows.size)

    def subset(self, keep: np.ndarray) -> "SplitData":
        return SplitData(
            name=self.name,
            rows=self.rows[keep],
            fp_ix=self.fp_ix[keep],
            dose_scale=self.dose_scale[keep],
            cg_ix=self.cg_ix[keep],
        )


# --------------------------------------------------------------------------- #
# Training
# --------------------------------------------------------------------------- #
def gaussian_nll(model: PGM, control, label, noise, target) -> torch.Tensor:
    """``GUSS`` — Gaussian negative log-likelihood; **may be negative, never clamped**."""
    mean, var = model(control, label, noise)
    return F.gaussian_nll_loss(mean, target, var, reduction="mean")


def save_checkpoint(path: Path, model: PGM, config: dict, seed: int, best_valid: float, epoch: int) -> None:
    state = {k: v.detach().cpu() for k, v in model.state_dict().items()}
    torch.save(
        {
            "model": "prnet",
            "state_dict": state,
            "config": config,
            "seed": seed,
            "n_genes": N_GENES,
            "best_valid_loss": best_valid,
            "epoch": epoch,
        },
        path,
    )


def epoch_batches(n_rows: int, batch_size: int = BATCH_SIZE) -> list[tuple[int, int]]:
    """Contiguous ``[start, stop)`` slices over one shuffled epoch.

    A trailing batch of a **single** row is folded into the previous slice instead of
    being dropped: ``BatchNorm1d`` needs more than one sample in train mode, and
    dropping the row would silently bias the epoch mean.  Every row is trained exactly
    once per epoch.
    """
    if n_rows < 2:
        raise SystemExit(f"need at least 2 training rows to train BatchNorm1d, got {n_rows}")
    starts = list(range(0, n_rows, batch_size))
    if len(starts) > 1 and n_rows - starts[-1] == 1:
        starts.pop()  # fold the singleton into the previous batch
    stops = starts[1:] + [n_rows]
    return list(zip(starts, stops))


def train(
    model: PGM,
    optimizer: torch.optim.Optimizer,
    scheduler: torch.optim.lr_scheduler.ReduceLROnPlateau,
    gatherer: BatchGatherer,
    train_split: SplitData,
    valid_split: SplitData,
    max_epochs: int,
    checkpoint_path: Path,
    config: dict,
    smoke: bool,
) -> dict:
    device = next(model.parameters()).device
    cpu_generator = torch.Generator(device="cpu").manual_seed(SEED)
    history: list[dict] = []
    best_valid = math.inf
    best_epoch = -1
    stale = 0
    stopped_early = False

    for epoch in range(max_epochs):
        started = time.perf_counter()
        model.train()
        order = torch.randperm(len(train_split), generator=cpu_generator).numpy()
        total, seen = 0.0, 0
        # Batch-size-weighted epoch loss: Σ(NLL_batch × n_batch) / Σn_batch.  A mean of
        # per-batch means would over-weight the final (small) batch.
        for start, stop in epoch_batches(len(train_split)):
            idx = order[start:stop]
            control, label, target = gatherer.batch(train_split, idx)
            noise = torch.randn(control.shape[0], NOISE_DIM, generator=cpu_generator).to(device)
            loss = gaussian_nll(model, control, label, noise, target)
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()
            total += float(loss.item()) * idx.size
            seen += int(idx.size)
        train_loss = total / max(seen, 1)

        model.eval()
        valid_total, valid_seen = 0.0, 0
        with torch.no_grad():
            for start in range(0, len(valid_split), BATCH_SIZE):
                stop = min(start + BATCH_SIZE, len(valid_split))
                idx = np.arange(start, stop)
                control, label, target = gatherer.batch(valid_split, idx)
                noise = torch.zeros(control.shape[0], NOISE_DIM, device=device)
                loss = gaussian_nll(model, control, label, noise, target)
                valid_total += float(loss.item()) * idx.size
                valid_seen += int(idx.size)
        valid_loss = valid_total / max(valid_seen, 1)
        scheduler.step(valid_loss)
        lr = float(optimizer.param_groups[0]["lr"])
        seconds = time.perf_counter() - started

        history.append(
            {
                "epoch": epoch + 1,
                "train_loss": train_loss,
                "valid_loss": valid_loss,
                "lr": lr,
                "seconds": seconds,
            }
        )

        improved = valid_loss < best_valid
        if improved:
            best_valid = valid_loss
            best_epoch = epoch + 1
            stale = 0
            save_checkpoint(checkpoint_path, model, config, SEED, best_valid, epoch)
        else:
            stale += 1
        min_epochs_met = (epoch + 1) >= MIN_EPOCHS
        log.info(
            "epoch %d/%d | train=%.6f | valid=%.6f | lr=%.2e | %.1fs | "
            "patience %d/%d | min_epochs %s%s",
            epoch + 1,
            max_epochs,
            train_loss,
            valid_loss,
            lr,
            seconds,
            stale,
            EARLY_STOP_PATIENCE,
            "met" if min_epochs_met else "not met",
            " (smoke: early stop disabled)" if smoke else "",
        )
        if improved:
            log.info("new best valid loss %.6f -> %s", best_valid, checkpoint_path.name)
        if not math.isfinite(train_loss) or not math.isfinite(valid_loss):
            log.warning("epoch %d produced a non-finite loss (not clamped)", epoch + 1)

        if not smoke and min_epochs_met and stale > EARLY_STOP_PATIENCE:
            stopped_early = True
            log.info(
                "early stop at epoch %d: valid loss has not improved for %d epochs",
                epoch + 1,
                stale,
            )
            break

    if best_epoch < 0:  # no finite improvement: still leave a usable checkpoint
        save_checkpoint(checkpoint_path, model, config, SEED, best_valid, len(history) - 1)
        log.warning("no epoch improved the valid loss; saved the final weights instead")

    return {
        "history": history,
        "best_valid_loss": best_valid,
        "best_epoch": best_epoch,
        "epochs_run": len(history),
        "stopped_early": stopped_early,
        "stop_epoch": len(history),
    }


# --------------------------------------------------------------------------- #
# Outputs
# --------------------------------------------------------------------------- #
def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")


def write_summary_table(path: Path, history: list[dict]) -> None:
    frame = pd.DataFrame(history, columns=["epoch", "train_loss", "valid_loss", "lr", "seconds"])
    frame.to_csv(path, index=False, float_format="%.10g")


def write_summary_png(
    path: Path,
    history: list[dict],
    best_epoch: int,
    best_valid: float,
    stopped_early: bool,
    stop_epoch: int,
    smoke: bool,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    epochs = [h["epoch"] for h in history]
    train_losses = [h["train_loss"] for h in history]
    valid_losses = [h["valid_loss"] for h in history]

    train_color, valid_color, best_color = "#3f6fa8", "#e08a4a", "#c0392b"
    fig, axis = plt.subplots(figsize=(8.0, 4.8))
    axis.plot(
        epochs, train_losses, color=train_color, marker="o", markersize=4, lw=1.8,
        label="train", clip_on=False,
    )
    axis.plot(
        epochs, valid_losses, color=valid_color, marker="s", markersize=4, lw=1.8,
        label="valid", clip_on=False,
    )
    if best_epoch > 0:
        axis.plot(
            [best_epoch],
            [best_valid],
            marker="*",
            markersize=16,
            color=best_color,
            linestyle="none",
            label=f"best epoch {best_epoch} ({best_valid:.5f})",
            zorder=5,
            clip_on=False,
        )
        axis.annotate(
            f"{best_valid:.5f}",
            (best_epoch, best_valid),
            textcoords="offset points",
            xytext=(9, 9),
            fontsize=7,
            color=best_color,
        )
    elif best_epoch < 0:
        axis.axhline(0.0, color="#9aa0a6", lw=0.8, ls="--")
    if stopped_early:
        axis.axvline(
            stop_epoch,
            color="#6b7280",
            lw=1.2,
            ls=":",
            label=f"early stop (epoch {stop_epoch})",
        )

    axis.set_xlabel("epoch")
    axis.set_ylabel("Gaussian NLL (may be negative)")
    title = "PRnet pretrain — Gaussian NLL per epoch"
    if smoke:
        title += "  [smoke: synthetic fixture, 2 epochs]"
    axis.set_title(title)
    axis.grid(True, color="#e3e3e3", lw=0.6)
    axis.legend(frameon=False, fontsize=8)
    # Epoch is a discrete count: integer ticks only, and the axis spans exactly
    # [1, epochs_run] so the curve fills the panel instead of floating inside it.
    # A single-epoch history cannot use [1, 1] as a limit, so pad that one case.
    if max(epochs) > 1:
        axis.set_xlim(1, max(epochs))
    else:
        axis.set_xlim(0.5, 1.5)
    axis.xaxis.set_major_locator(MaxNLocator(integer=True))
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_model_card(
    path: Path,
    args: argparse.Namespace,
    view: Path,
    meta: ViewMeta,
    selection: Selection,
    config: dict,
    result: dict,
    device: torch.device,
    device_name: str,
    timing: dict,
) -> None:
    history = result["history"]
    first, last = history[0], history[-1]
    best_epoch = result["best_epoch"]
    best_valid = result["best_valid_loss"]
    lines = [
        "# PRnet pretrain — model card",
        "",
        "Auto-generated by `stage2/05_prnet_pretrain/execution/train_prnet.py`"
        " (spec: `.scratch/stage2-05-prnet-pretrain/spec.md` §2).",
        "",
        f"- **Mode:** {'SMOKE (synthetic fixture, 2 epochs)' if args.smoke else 'full pretrain run'}",
        f"- **Started:** {config['started']}",
        f"- **Finished:** {timing['finished']}",
        f"- **Wall time:** {timing['wall_time_s']:.1f} s end to end "
        f"({timing['train_loop_s']:.1f} s inside the {result['epochs_run']} training epochs)",
        f"- **Peak memory:** {timing['peak_rss_mb']:.0f} MiB host RSS "
        f"(the view is read as CSR, never densified)",
        f"- **Command:** `{' '.join(config['command'])}`",
        f"- **Device:** `{device}` — {device_name}",
        f"- **Env:** torch {torch.__version__}, CUDA available "
        f"{torch.cuda.is_available()}, python {platform.python_version()}",
        "",
        "## What was trained",
        "",
        "PRnet `PGM` (perturbation-conditioned generative model) on the Tahoe-100M",
        "pseudobulk pretrain view, from random initialisation.",
        "",
        f"- view: `{view}`",
        f"- shape: {meta.n_obs:,} rows x {meta.n_vars:,} genes; `X` is raw UMI counts (unused here),",
        "  the target is `layers['log1p_cp10k']`",
        f"- samples: {selection.counts['n_samples']:,} "
        f"({selection.counts['n_train']:,} train / {selection.counts['n_valid']:,} valid)",
        f"- controls (excluded): {selection.counts['n_controls']:,}, "
        f"perturbed: {selection.counts['n_perturbed']:,}",
        f"- dropped `has_structure == False`: {selection.counts['dropped_no_structure']:,} rows",
        f"- split: by plate, valid = {', '.join(VALID_PLATES)}",
        f"- control input: `uns['control_mean']` — {len(meta.control_labels)} `(plate, cell_line)` groups",
        f"  (on-disk encoding read: `{meta.control_mean_encoding}`)",
        "",
        "## Architecture",
        "",
        "Mirrors `pert_model/PRnet/source/models/PRnet.py` (re-implemented in-file, upstream tree not imported):",
        "",
        "- `PAdaptor` (drug adaptor): 1024 rFCFP4 -> 128 -> 64",
        "- `PEncoder`: [18164 + 64] -> 128 -> 64 (mean only)",
        "- `PDecoder`: [64 + 64 + 10] -> 128 -> 2 x 18164; `relu(head[:, :dim])` = mean,",
        "  `softplus(head[:, dim:])` = variance",
        "- dropout 0.05; weights `normal(0, 0.02)`, BatchNorm `normal(1, 0.02)`",
        "",
        "## Training recipe (same-source copy — do not diverge)",
        "",
        "Identical to `pert_model/PRnet/source/trainer/PRnetTrainer.py:156-258` and",
        "`stage2/03_models_cmbench/03.1_cold_drug_18k/execution/script/adapters/prnet/run.py:139-279`:",
        "",
        f"- loss: `F.gaussian_nll_loss(mean, target, var, reduction='mean')` (author `GUSS`);",
        "  it **can be negative** and is never clamped, `abs()`-ed or otherwise altered",
        f"- Adam(lr={LR:g}, weight_decay={WEIGHT_DECAY:g})",
        f"- ReduceLROnPlateau(mode='min', factor={SCHED_FACTOR:g}, patience={SCHED_PATIENCE}, "
        f"min_lr={MIN_LR:g}) stepped on **valid** loss",
        f"- batch {BATCH_SIZE}, max_epochs {MAX_EPOCHS}, early stop after {EARLY_STOP_PATIENCE} "
        f"stale epochs, min_epochs {MIN_EPOCHS}, seed {SEED}",
        f"- fingerprint: RDKit `GetMorganFingerprintAsBitVect(mol, 2, useFeatures=True, nBits=1024)`",
        "  x `log10(dose_uM + 1)`",
        "",
        "**Epoch loss is batch-size-weighted.** `fit_history.json` / `summary.png` report",
        "`sum(NLL_batch x n_batch) / sum(n_batch)` for both splits — not a mean of per-batch",
        "means — so the final partial batch cannot bias an epoch. A trailing batch of a",
        "single row is folded into the previous batch (BatchNorm1d needs >1 sample) and no",
        "row is dropped.",
        "",
        "## Result",
        "",
        f"- epochs run: {result['epochs_run']}",
        (
            f"- early stop: **fired at epoch {result['stop_epoch']}** — the valid loss had not "
            f"improved for more than {EARLY_STOP_PATIENCE} epochs and the {MIN_EPOCHS}-epoch "
            "floor was met"
            if result["stopped_early"]
            else f"- early stop: not fired — the run used its full {result['epochs_run']}-epoch "
            f"allocation (patience {EARLY_STOP_PATIENCE}, min_epochs {MIN_EPOCHS})"
        ),
        f"- first epoch: train {first['train_loss']:.6f} / valid {first['valid_loss']:.6f}",
        f"- best epoch {best_epoch}: valid {best_valid:.6f} (best validation loss)",
        f"- last epoch: train {last['train_loss']:.6f} / valid {last['valid_loss']:.6f}",
        "",
        "Artifacts (all in `output/`, overwritten on every rerun): `run.log`,",
        "`fit_history.json`, `summary.png`, `summary_table.csv`, `checkpoint/model.pt`,",
        "`config.json`, `model_card.md`.",
        "",
        "## Deviations inherited from the 04 view (D0–D3)",
        "",
        "- **D0 layout**: `X`/`log1p_cp10k` are CSR instead of the CM views' dense `int32`.",
        "  Dense would be 4.87 GB (9.74 GB with the layer), not 4.87 TB — the plan's units error —",
        "  and at 83% density this 16.32 GB CSR artifact is the larger container; CSR mirrors the",
        "  corpus's CSR parts and the slice-write assembly.",
        "- **D1 control keying**: `control_group` = `tahoe::<plate>::<cell_line_id>` (700 groups)",
        "  instead of the corpus's 14 plate-level groups, so the control mean never mixes cell lines.",
        "- **D2 gene axis**: 18,164 slots, `FAM153B` all-zero (no Tahoe measurement).",
        "- **D3 provenance gaps**: `time_h`/`cell_type`/`donor` NaN, `atc`/`moa`/`mesh_classes`",
        "  empty, `pharm_*` False — recorded rather than invented.",
        "",
        "## Caveats",
        "",
        "- Pipeline bring-up only: no CM finetune, no cold-drug evaluation, no comparison",
        "  against `03_models_cmbench`. A negative Gaussian NLL is expected and correct.",
        "- The control input is the group mean of raw counts, **not** log1p_cp10k; this",
        "  mirrors the author's dataset (control row and target row share the same space).",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "PRnet pretrain bring-up on the Tahoe-100M pseudobulk view (spec "
            "§2). Trains the author PGM with Gaussian NLL and writes the loss curve."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--view", type=Path, default=Path(VIEW_DEFAULT),
                        help="pretrain view .h5ad (04 output)")
    parser.add_argument("--output", type=Path, default=Path(OUTPUT_DEFAULT),
                        help="run directory (overwritten on rerun)")
    parser.add_argument("--device", default="cuda", choices=("cuda", "cpu"),
                        help="torch device; falls back to cpu if cuda is unavailable")
    parser.add_argument("--seed", type=int, default=SEED, help="RNG seed")
    parser.add_argument("--smoke", action="store_true",
                        help=f"pipeline check: {SMOKE_EPOCHS} epochs on "
                             f"{SMOKE_TRAIN_ROWS}/{SMOKE_VALID_ROWS} sampled rows")
    parser.add_argument("--check", action="store_true",
                        help="print the resolved paths, view shape and split counts, then exit")
    return parser.parse_args(argv)


def resolve_device(requested: str) -> tuple[torch.device, str]:
    available = torch.cuda.is_available()
    if requested.startswith("cuda") and available:
        return torch.device("cuda"), torch.cuda.get_device_name(0)
    if requested.startswith("cuda"):
        log.warning("--device %s requested but CUDA is unavailable; using CPU", requested)
        return torch.device("cpu"), platform.processor() or "cpu"
    return torch.device("cpu"), platform.processor() or "cpu"


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    root = find_repo_root(Path(__file__).resolve().parent)
    view = resolve(args.view, root)
    output = resolve(args.output, root)

    if not view.exists():
        raise SystemExit(
            f"view not found: {view}\n"
            "Build it first: `python stage2/04_tahoe_data/execution/build_pretrain_view.py`\n"
            "or smoke-test this step on the synthetic fixture:\n"
            "  python stage2/05_prnet_pretrain/execution/_smoke/make_fixture.py\n"
            "  python stage2/05_prnet_pretrain/execution/train_prnet.py --smoke "
            "--view <fixture>.h5ad"
        )

    output.mkdir(parents=True, exist_ok=True)
    setup_logging(None if args.check else output / "run.log")
    started = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    run_t0 = time.perf_counter()
    command = [sys.executable, *sys.argv] if argv is None else ["train_prnet.py", *argv]
    log.info("PRnet pretrain bring-up — spec §2 (stage2/05_prnet_pretrain)")
    log.info("repo root %s", root)
    log.info("view %s", view)
    log.info("output %s", output)
    log.info("smoke=%s check=%s", args.smoke, args.check)

    device, device_name = resolve_device(args.device)
    log.info(
        "device=%s cuda_available=%s torch=%s gpu=%s",
        device,
        torch.cuda.is_available(),
        torch.__version__,
        device_name,
    )

    torch.manual_seed(args.seed)
    np.random.seed(args.seed % 2**32)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(args.seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

    meta, layer = load_view(view, with_layer=not args.check)
    selection = build_selection(meta)
    if args.check:
        log.info(
            "check OK: view %s has %d x %d, %d control groups, %d samples",
            view,
            meta.n_obs,
            meta.n_vars,
            len(meta.control_labels),
            selection.counts["n_samples"],
        )
        return 0

    # The fingerprint cache built by `build_selection` is kept as-is: its row order is
    # exactly the `fp_ix` encoding, so nothing is recomputed here.
    train_split = SplitData(
        "train",
        rows=selection.rows[selection.train],
        fp_ix=selection.fp_ix[selection.train],
        dose_scale=selection.dose_scale[selection.train],
        cg_ix=selection.cg_ix[selection.train],
    )
    valid_split = SplitData(
        "valid",
        rows=selection.rows[selection.valid],
        fp_ix=selection.fp_ix[selection.valid],
        dose_scale=selection.dose_scale[selection.valid],
        cg_ix=selection.cg_ix[selection.valid],
    )
    if args.smoke:
        train_split = train_split.subset(np.arange(min(SMOKE_TRAIN_ROWS, len(train_split))))
        valid_split = valid_split.subset(np.arange(min(SMOKE_VALID_ROWS, len(valid_split))))
        log.info(
            "smoke: train=%d valid=%d rows, %d epochs",
            len(train_split),
            len(valid_split),
            SMOKE_EPOCHS,
        )

    assert layer is not None
    gatherer = BatchGatherer(
        fingerprints=selection.fingerprints,
        control_values=meta.control_values,
        layer=layer,
        device=device,
    )

    model = PGM().to(device)
    model.apply(weight_init)
    check_architecture(model)
    n_params = sum(p.numel() for p in model.parameters())
    log.info("PGM built: %d parameters; architecture matches spec §2", n_params)

    optimizer = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode="min",
        factor=SCHED_FACTOR,
        patience=SCHED_PATIENCE,
        min_lr=MIN_LR,
    )

    max_epochs = SMOKE_EPOCHS if args.smoke else MAX_EPOCHS
    config = {
        "program": "stage2/05_prnet_pretrain/execution/train_prnet.py",
        "spec": ".scratch/stage2-05-prnet-pretrain/spec.md §2",
        "same_source_copy": [
            "pert_model/PRnet/source/models/PRnet.py",
            "pert_model/PRnet/source/trainer/PRnetTrainer.py",
            "pert_model/PRnet/source/data/_utils.py",
            "stage2/03_models_cmbench/03.1_cold_drug_18k/execution/script/adapters/prnet/run.py",
        ],
        "started": started,
        "command": command,
        "smoke": bool(args.smoke),
        "seed": args.seed,
        "device": str(device),
        "device_name": device_name,
        "torch": torch.__version__,
        "cuda_available": bool(torch.cuda.is_available()),
        "view": str(view.relative_to(root)) if view.is_relative_to(root) else str(view),
        "view_shape": [meta.n_obs, meta.n_vars],
        "n_control_groups": len(meta.control_labels),
        "control_mean": {
            "n_groups": len(meta.control_labels),
            "n_genes": meta.n_vars,
            "encoding": meta.control_mean_encoding,
        },
        "split": {"valid_plates": list(VALID_PLATES), "n_train": len(train_split), "n_valid": len(valid_split)},
        "rows": selection.counts,
        "model": {
            "class": "PGM",
            "x_dimension": N_GENES,
            "drug_dimension": DRUG_DIM,
            "comb_dimension": COMB_DIM,
            "noise_dim": NOISE_DIM,
            "hidden_layer_sizes": list(HIDDEN),
            "adaptor_layer_sizes": list(ADAPTOR_HIDDEN),
            "z_dimension": Z_DIM,
            "dr_rate": DROP_RATE,
            "n_parameters": n_params,
        },
        "loss": "F.gaussian_nll_loss(mean, target, var, reduction='mean')",
        "loss_reduction": {
            "per_epoch": "sum(NLL_batch * n_batch) / sum(n_batch) — batch-size weighted",
            "trailing_singleton_batch": "folded into the previous batch (BatchNorm1d needs >1 sample)",
        },
        "optimizer": {"name": "adam", "lr": LR, "weight_decay": WEIGHT_DECAY},
        "scheduler": {
            "name": "reduce_on_plateau",
            "mode": "min",
            "factor": SCHED_FACTOR,
            "patience": SCHED_PATIENCE,
            "min_lr": MIN_LR,
        },
        "batch_size": BATCH_SIZE,
        "max_epochs": max_epochs,
        "early_stop_patience": EARLY_STOP_PATIENCE,
        "min_epochs": MIN_EPOCHS,
        "peak_rss_mb_after_load": round(peak_rss_mb(), 1),
    }
    write_json(output / "config.json", config)

    train_batches = epoch_batches(len(train_split))
    valid_batches = math.ceil(len(valid_split) / BATCH_SIZE)
    log.info(
        "training: %d epochs | split: train %d rows (%d batch%s, last %d rows), "
        "valid %d rows (%d batch%s) | per-epoch loss = sum(NLL_batch * n_batch) / sum(n_batch)",
        max_epochs,
        len(train_split),
        len(train_batches),
        "" if len(train_batches) == 1 else "es",
        train_batches[-1][1] - train_batches[-1][0],
        len(valid_split),
        valid_batches,
        "" if valid_batches == 1 else "es",
    )
    checkpoint_path = output / "checkpoint" / "model.pt"
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    result = train(
        model=model,
        optimizer=optimizer,
        scheduler=scheduler,
        gatherer=gatherer,
        train_split=train_split,
        valid_split=valid_split,
        max_epochs=max_epochs,
        checkpoint_path=checkpoint_path,
        config=config,
        smoke=args.smoke,
    )

    finished = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    wall_time_s = time.perf_counter() - run_t0
    train_loop_s = sum(float(h["seconds"]) for h in result["history"])
    peak_rss = peak_rss_mb()
    timing = {
        "finished": finished,
        "wall_time_s": round(wall_time_s, 1),
        "train_loop_s": round(train_loop_s, 1),
        "peak_rss_mb": round(peak_rss, 1),
    }
    fit_history = {
        "run": {
            "smoke": bool(args.smoke),
            "device": str(device),
            "device_name": device_name,
            "view": config["view"],
            "seed": args.seed,
            "best_epoch": result["best_epoch"],
            "best_valid_loss": result["best_valid_loss"],
            "epochs_run": result["epochs_run"],
            "stopped_early": result["stopped_early"],
            "stop_epoch": result["stop_epoch"],
            "started": config["started"],
            "finished": finished,
            "wall_time_s": timing["wall_time_s"],
            "train_loop_s": timing["train_loop_s"],
            "peak_rss_mb": timing["peak_rss_mb"],
        },
        "history": result["history"],
    }
    write_json(output / "fit_history.json", fit_history)
    write_summary_table(output / "summary_table.csv", result["history"])
    write_summary_png(
        output / "summary.png",
        result["history"],
        result["best_epoch"],
        result["best_valid_loss"],
        result["stopped_early"],
        result["stop_epoch"],
        args.smoke,
    )
    write_model_card(
        output / "model_card.md",
        args,
        view,
        meta,
        selection,
        config,
        result,
        device,
        device_name,
        timing,
    )

    log.info(
        "done: %d epochs, best valid %.6f at epoch %d%s | wall %.1f s | peak rss %.0f MiB",
        result["epochs_run"],
        result["best_valid_loss"],
        result["best_epoch"],
        " (early stop)" if result["stopped_early"] else "",
        timing["wall_time_s"],
        peak_rss,
    )
    for name in ("run.log", "config.json", "fit_history.json", "summary_table.csv",
                 "summary.png", "checkpoint/model.pt", "model_card.md"):
        path = output / name
        log.info("wrote %s (%s)", path, f"{path.stat().st_size:,} bytes" if path.exists() else "MISSING")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
