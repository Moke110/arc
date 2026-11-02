"""Ticket #7 (06.7_cm_tahoe): the Tahoe-init candidate matrix and LR selection.

Pure logic, no torch and no data access.  ``build_queue.py``, ``run_matrix.py``,
``write_record.py`` and ``tests/`` all import this module, so the matrix and the
selection rule have exactly one definition.

The matrix is the ticket-#6 matrix — 5 folds x 3 seeds x 3 learning rates = 45
**full-parameter CM finetunes** — with exactly one change: every candidate warm
starts from that seed's ticket-#5 formal Tahoe pretrain checkpoint.  Folds, seeds,
learning rates, budget, early-stopping rule and per-fold LR selection are derived
from the same frozen ``06.3_input/output/locked_config.json`` (``n_folds``,
``seeds``, ``lrs.cm``) as ticket #6, so the two branches stay symmetric and paired.

Selection rule (locked config ``lr_selection``): for each fold, choose the common
LR that maximises the mean over the three seeds of the best validation direction
AUPRC@200 at k=200; on a tie choose the lower LR.  Test outcomes never enter the
selection — the pipeline's trainer only sees the fold's train/valid roles and the
test split is postponed to ticket #9.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

__all__ = [
    "MatrixError",
    "SelectionError",
    "Candidate",
    "UNIT_ID_RE",
    "INIT_MODE",
    "TRAIN_INIT_FLAG",
    "LR_TIE_ABS_TOL",
    "COMPLETE_STATUSES",
    "lr_label",
    "unit_id_for",
    "parse_unit_id",
    "matrix_spec",
    "candidates",
    "checkpoint_path",
    "seed_checkpoint_map",
    "train_args",
    "derive_status",
    "select_lr",
    "fold_selections",
]

#: Branch label of this substep (unit-id suffix, record column, ticket #9 pairing).
INIT_MODE = "tahoe"
#: The ticket-#4 trainer flag that actually performs the warm start.
TRAIN_INIT_FLAG = "pretrained"

#: Two LR means inside this absolute tolerance count as a tie ("lower LR on tie").
LR_TIE_ABS_TOL = 1e-12

#: Candidate statuses that count as a usable validation score.
COMPLETE_STATUSES = ("ok", "capped")

UNIT_ID_RE = re.compile(r"^cm_f(?P<fold>\d+)_s(?P<seed>\d+)_lr(?P<lr>[0-9.eE+-]+)_tahoe$")

#: checkpoint location inside ticket #5's formal run tree
CHECKPOINT_RELPATH = Path("checkpoint") / "model.pt"


class MatrixError(ValueError):
    """The frozen locked config cannot define the ticket-#7 matrix."""


class SelectionError(RuntimeError):
    """A fold's LR selection is incomplete or ambiguous."""


@dataclass(frozen=True)
class Candidate:
    """One ``(fold, seed, lr)`` training unit of the Tahoe-init matrix."""

    fold: int
    seed: int
    seed_index: int
    lr: float
    lr_label: str
    order: int
    unit_id: str

    def as_row(self) -> dict[str, Any]:
        return {
            "unit_id": self.unit_id,
            "fold": self.fold,
            "seed": self.seed,
            "seed_index": self.seed_index,
            "lr": self.lr,
            "lr_label": self.lr_label,
            "order": self.order,
        }


def lr_label(lr: float) -> str:
    """Stable short label for a learning rate: ``1e-4`` / ``3e-4`` / ``1e-3``."""
    value = float(lr)
    if not math.isfinite(value) or value <= 0:
        raise MatrixError(f"learning rate must be a positive finite number, got {lr!r}")
    exponent = math.floor(math.log10(value))
    mantissa = round(value / (10.0 ** exponent), 6)
    # 10.0 / 0.1 style rounding can land on 10; normalise to the next decade.
    if abs(mantissa - 10.0) < 1e-9:
        mantissa, exponent = 1.0, exponent + 1
    return f"{mantissa:g}e{exponent}"


def unit_id_for(fold: int, seed: int, lr: float) -> str:
    return f"cm_f{int(fold)}_s{int(seed)}_lr{lr_label(lr)}_{INIT_MODE}"


def parse_unit_id(unit_id: str) -> dict[str, Any]:
    match = UNIT_ID_RE.match(str(unit_id))
    if not match:
        raise MatrixError(f"not a ticket-#7 unit id: {unit_id!r}")
    return {
        "fold": int(match.group("fold")),
        "seed": int(match.group("seed")),
        "lr_label": match.group("lr"),
    }


def matrix_spec(locked: Mapping[str, Any]) -> dict[str, Any]:
    """Derive the matrix from the frozen locked config; fail closed on drift."""
    try:
        n_folds = int(locked["n_folds"])
        seeds = [int(seed) for seed in locked["seeds"]]
        lrs = [float(lr) for lr in locked["lrs"]["cm"]]
    except (KeyError, TypeError, ValueError) as exc:
        raise MatrixError(f"locked config has no usable folds/seeds/lrs.cm block: {exc}") from exc
    if n_folds < 1:
        raise MatrixError(f"n_folds must be positive, got {n_folds}")
    if len(seeds) != 3:
        raise MatrixError(f"the spec fixes three seeds, locked config has {seeds!r}")
    if len(set(seeds)) != len(seeds):
        raise MatrixError(f"locked config seeds are not unique: {seeds!r}")
    if len(lrs) != 3:
        raise MatrixError(f"the spec fixes three CM learning rates, locked config has {lrs!r}")
    if any(not math.isfinite(lr) or lr <= 0 for lr in lrs):
        raise MatrixError(f"locked config lrs.cm are not positive finite numbers: {lrs!r}")
    if len(set(lrs)) != len(lrs):
        raise MatrixError(f"locked config lrs.cm are not unique: {lrs!r}")
    if any(a >= b for a, b in zip(lrs, lrs[1:])):
        # the tie-break is "lower LR first", which needs the ascending candidate order
        raise MatrixError(f"locked config lrs.cm must be strictly ascending: {lrs!r}")
    labels = [lr_label(lr) for lr in lrs]
    if len(set(labels)) != len(labels):
        raise MatrixError(f"learning-rate labels collide: {labels!r}")
    return {
        "init": INIT_MODE,
        "pipeline_init": TRAIN_INIT_FLAG,
        "n_folds": n_folds,
        "folds": list(range(n_folds)),
        "seeds": seeds,
        "lrs": lrs,
        "lr_labels": labels,
        "n_candidates": n_folds * len(seeds) * len(lrs),
        "candidate_order": "ascending LR (the tie-break order)",
    }


def candidates(spec: Mapping[str, Any]) -> list[Candidate]:
    """All matrix units in a deterministic order: fold, then ascending LR, then seed."""
    out: list[Candidate] = []
    order = 0
    lrs: Sequence[float] = spec["lrs"]
    for fold in spec["folds"]:
        for lr in lrs:
            for seed_index, seed in enumerate(spec["seeds"]):
                out.append(
                    Candidate(
                        fold=int(fold),
                        seed=int(seed),
                        seed_index=int(seed_index),
                        lr=float(lr),
                        lr_label=lr_label(lr),
                        order=order,
                        unit_id=unit_id_for(fold, seed, lr),
                    )
                )
                order += 1
    return out


def checkpoint_path(checkpoint_root: Path | str, seed: int) -> Path:
    """``<checkpoint_root>/tahoe_seed<seed>/checkpoint/model.pt`` (ticket #5 layout)."""
    return Path(checkpoint_root) / f"tahoe_seed{int(seed)}" / CHECKPOINT_RELPATH


def seed_checkpoint_map(
    seeds: Iterable[int], checkpoint_root: Path | str
) -> dict[int, Path]:
    """One checkpoint path per matrix seed; every path must exist (fail closed)."""
    mapping: dict[int, Path] = {}
    for seed in seeds:
        path = checkpoint_path(checkpoint_root, int(seed))
        if not path.is_file():
            raise MatrixError(
                f"seed {int(seed)} has no ticket-#5 formal checkpoint at {path}; "
                "the formal pretrain sweep must finish before the ticket-#7 matrix can run"
            )
        mapping[int(seed)] = path
    return mapping


def train_args(
    candidate: Candidate,
    init_checkpoint: Path | str,
    hash_cache: str | None = None,
) -> list[str]:
    """The exact ``train.py`` arguments for one candidate (queue owns out/data-dir)."""
    args = [
        "--fold", str(candidate.fold),
        "--seed", str(candidate.seed),
        "--lr", repr(float(candidate.lr)),
        "--init", TRAIN_INIT_FLAG,
        "--init-checkpoint", str(init_checkpoint),
        "--device", "cuda",
    ]
    if hash_cache:
        # Memoise the frozen-artifact sha256 by (size, mtime) so 45 runs do not
        # each re-hash the 15.3 GB frozen input; the queue forwards unit args.
        args.extend(["--hash-cache", str(hash_cache)])
    return args


# --------------------------------------------------------------------------- #
# status derivation from one run directory
# --------------------------------------------------------------------------- #
def derive_status(
    manifest: Mapping[str, Any] | None,
    fit_history: Mapping[str, Any] | None,
    checkpoint_exists: bool,
    checkpoint_sha256_matches: bool,
) -> tuple[str, str]:
    """Classify a candidate as ``ok`` / ``capped`` / ``failed`` plus a detail string.

    ``ok``     — completed, checkpoint reloaded identical, stopped early on
                 validation AUPRC (the expected outcome).
    ``capped`` — completed identically but consumed the whole 100-epoch budget
                 without early stopping (the budget bound, not a failure).
    ``failed`` — anything else (missing/incomplete manifest, missing or tampered
                 checkpoint, degenerate no-finite-AUPRC run).
    """
    if manifest is None:
        return "failed", "run directory holds no manifest.json (not completed)"
    if manifest.get("status") != "completed":
        return "failed", f"manifest status is {manifest.get('status')!r}"
    checkpoint = manifest.get("checkpoint") or {}
    if not checkpoint_exists:
        return "failed", "checkpoint file is missing"
    if not checkpoint_sha256_matches:
        return "failed", "checkpoint sha256 does not match the manifest"
    if checkpoint.get("reload_identical") is not True:
        return "failed", "checkpoint did not strict-reload identically"
    auprc = manifest.get("best_valid_auprc")
    if auprc is None or not math.isfinite(float(auprc)):
        return "failed", "no finite validation AUPRC"
    if fit_history is not None:
        history = fit_history.get("history") or []
        epochs_run = manifest.get("epochs_run")
        if epochs_run is not None and int(epochs_run) != len(history):
            return "failed", (
                f"manifest epochs_run={epochs_run} disagrees with fit_history rows={len(history)}"
            )
    stop_reason = str(manifest.get("stop_reason") or "")
    if stop_reason == "max_epochs":
        return "capped", f"consumed the full budget ({manifest.get('epochs_run')} epochs)"
    if stop_reason.startswith("early_stop"):
        return "ok", stop_reason
    return "failed", f"unexpected stop reason {stop_reason!r}"


# --------------------------------------------------------------------------- #
# per-fold LR selection
# --------------------------------------------------------------------------- #
def _mean(values: Sequence[float]) -> float:
    return math.fsum(values) / float(len(values))


def select_lr(
    cands: Iterable[Candidate],
    records: Mapping[str, Mapping[str, Any]],
    *,
    tol: float = LR_TIE_ABS_TOL,
) -> list[dict[str, Any]]:
    """Choose each fold's common LR from the three seeds' validation AUPRC.

    ``records`` maps ``unit_id`` to a row with at least ``status`` and
    ``best_valid_auprc``.  A fold's LR is eligible only when all three seeds of
    that LR are complete (``ok``/``capped``) with a finite AUPRC: the spec's
    rule is a mean over three seeds, never a mean over whatever survived.  When
    no LR is complete the fold is returned with ``selected_lr=None`` and
    ``complete=False`` — callers must treat that as a failed matrix, not invent
    a pick from partial seeds.
    """
    cands = list(cands)
    if not cands:
        raise SelectionError("no candidates to select from")
    by_fold: dict[int, list[Candidate]] = {}
    for candidate in cands:
        by_fold.setdefault(candidate.fold, []).append(candidate)

    folds: list[dict[str, Any]] = []
    for fold in sorted(by_fold):
        fold_candidates = by_fold[fold]
        lrs = sorted({candidate.lr for candidate in fold_candidates})
        n_seeds = len({candidate.seed for candidate in fold_candidates})
        score_rows: list[dict[str, Any]] = []
        for lr in lrs:  # ascending: candidate order == lower LR first
            lr_candidates = sorted(
                (c for c in fold_candidates if c.lr == lr), key=lambda c: c.seed_index
            )
            values: list[float] = []
            used: list[int] = []
            missing: list[int] = []
            for candidate in lr_candidates:
                row = records.get(candidate.unit_id) or {}
                status = str(row.get("status") or "")
                value = row.get("best_valid_auprc")
                if status in COMPLETE_STATUSES and value is not None and math.isfinite(float(value)):
                    values.append(float(value))
                    used.append(int(candidate.seed))
                else:
                    missing.append(int(candidate.seed))
            score_rows.append(
                {
                    "fold": fold,
                    "lr": lr,
                    "lr_label": lr_label(lr),
                    "n_seeds": n_seeds,
                    "n_seeds_ok": len(values),
                    # a partial mean is never reported as a candidate score
                    "mean_auprc": _mean(values) if len(values) == n_seeds else None,
                    "seeds_used": used,
                    "seeds_missing": missing,
                    "complete": len(values) == n_seeds,
                }
            )

        eligible = [row for row in score_rows if row["complete"] and row["mean_auprc"] is not None]
        if not eligible:
            folds.append(
                {
                    "fold": fold,
                    "selected_lr": None,
                    "selected_lr_label": None,
                    "selected_mean_auprc": None,
                    "tie_break_applied": False,
                    "complete": False,
                    "n_seeds": n_seeds,
                    "note": "no learning rate has a complete three-seed validation mean",
                    "lr_scores": score_rows,
                }
            )
            continue

        best = max(row["mean_auprc"] for row in eligible)
        tied = [row for row in eligible if abs(row["mean_auprc"] - best) <= tol]
        chosen = tied[0]  # ascending LR order => lower LR first
        folds.append(
            {
                "fold": fold,
                "selected_lr": chosen["lr"],
                "selected_lr_label": chosen["lr_label"],
                "selected_mean_auprc": chosen["mean_auprc"],
                "tie_break_applied": len(tied) > 1,
                "complete": True,
                "n_seeds": n_seeds,
                "note": (
                    "lower LR on tie" if len(tied) > 1 else "single highest three-seed mean"
                ),
                "lr_scores": score_rows,
            }
        )
    return folds


def fold_selections(
    cands: Iterable[Candidate],
    records: Mapping[str, Mapping[str, Any]],
) -> dict[int, dict[str, Any]]:
    """Convenience: ``select_lr`` keyed by fold."""
    return {row["fold"]: row for row in select_lr(cands, records)}
