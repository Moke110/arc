"""The evaluator: the only code path that reads truth Δ for scored rows.

Two entry points:

* :func:`score_validation` — called by the trainer every epoch.  It scores the
  in-memory model on the validation rows with zero noise, aggregates replicate
  predictions **within a condition first**, then computes direction AUPRC@k over
  the condition matrix and macro-averages.  It also returns the per-epoch
  validation cost split (gather / forward / AUPRC) and auxiliary Δ RMSE /
  Pearson / Spearman.
* :func:`run_evaluation` — the offline ``evaluate.py`` path.  It reads a
  ``predict.py`` output directory (predicted expression + Δ + row ids), verifies
  every row id against the frozen split, computes the **experimental truth Δ**
  from matched DMSO controls (not from the model input), scores it with the same
  frozen metric, and writes ``metrics.csv`` / ``per_condition.csv`` /
  ``aggregate.json`` / ``deltas.npz``.

Truth control rules: CM truth Δ is ``treated log1p - per-plate matched DMSO
mean`` (:class:`ptlib.data.CMControlPool`); Tahoe truth Δ is ``treated log1p -
matched (plate, cell_line) control mean``.  Neither is ever a model input at
validation/test time — the model sees only c* (CM) or the matched control of a
non-scored split (Tahoe).
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd
import torch

from . import PIPELINE_VERSION
from . import data as data_mod
from . import h5io
from . import manifest as manifest_mod
from . import metrics as metrics_mod

__all__ = [
    "EvaluateOptions",
    "score_validation",
    "run_evaluation",
    "expected_eval_rows",
    "load_predictions",
]

log = logging.getLogger("ptlib.evaluate")


@dataclass
class EvaluateOptions:
    stage: str
    pred_dir: Path
    out_dir: Path
    data_dir: Path
    role: str
    fold: int | None = None
    hash_cache: Path | None = None
    table_cache: Path | None = None
    command: list[str] | None = None


def _as_device_batch(array: np.ndarray, device: torch.device) -> torch.Tensor:
    return torch.from_numpy(np.ascontiguousarray(array)).to(device)


# --------------------------------------------------------------------------- #
# Per-epoch validation scoring
# --------------------------------------------------------------------------- #
def score_validation(
    model: torch.nn.Module,
    dataset: Any,
    rows: np.ndarray,
    condition_ids: np.ndarray,
    k: int,
    batch_size: int,
    device: torch.device,
    n_dim: int,
    truth_baseline: np.ndarray | None = None,
) -> dict[str, Any]:
    """Deterministic validation scoring (zero noise, fixed order).

    ``condition_ids`` maps each validation row to its condition (``drug x dose``
    for CM, ``drug x dose x cell_line`` for Tahoe).  ``truth_baseline`` is the
    matched-control matrix used only for the truth Δ (CM); when ``None`` the
    matched control equals the model input (Tahoe).
    """
    rows = np.asarray(rows, dtype=np.int64)
    condition_ids = np.asarray(condition_ids, dtype=np.int64)
    if rows.size != condition_ids.size:
        raise ValueError("rows and condition_ids must be the same length")
    n_cond = int(condition_ids.max()) + 1 if condition_ids.size else 0
    n_genes = int(dataset.frozen.n_genes)
    sum_pred = np.zeros((n_cond, n_genes), dtype=np.float64)
    sum_true = np.zeros((n_cond, n_genes), dtype=np.float64)
    counts = np.zeros(n_cond, dtype=np.int64)
    nll_total, nll_seen = 0.0, 0
    gather_s = forward_s = 0.0
    zero_noise = torch.zeros(1, n_dim, device=device)

    with torch.no_grad():
        for start in range(0, rows.size, int(batch_size)):
            stop = min(start + int(batch_size), rows.size)
            batch_rows = rows[start:stop]
            batch_cond = condition_ids[start:stop]

            t0 = time.perf_counter()
            control, label, target = dataset.batch(batch_rows)
            input_baseline = control  # the model input baseline
            truth_base = (
                np.asarray(truth_baseline[start:stop], dtype=np.float32)
                if truth_baseline is not None
                else input_baseline
            )
            gather_s += time.perf_counter() - t0

            t0 = time.perf_counter()
            noise = zero_noise.expand(batch_rows.size, n_dim)
            mean, var = model(
                _as_device_batch(control, device),
                _as_device_batch(label, device),
                noise,
            )
            nll = torch.nn.functional.gaussian_nll_loss(
                mean, _as_device_batch(target, device), var, reduction="mean"
            )
            forward_s += time.perf_counter() - t0
            nll_total += float(nll.item()) * int(batch_rows.size)
            nll_seen += int(batch_rows.size)

            pred = mean.detach().cpu().numpy().astype(np.float64)
            true_delta = np.asarray(target, dtype=np.float64) - np.asarray(
                truth_base, dtype=np.float64
            )
            pred_delta = pred - np.asarray(input_baseline, dtype=np.float64)
            np.add.at(sum_pred, batch_cond, pred_delta)
            np.add.at(sum_true, batch_cond, true_delta)
            np.add.at(counts, batch_cond, 1)

    if np.any(counts == 0):
        raise ValueError("a validation condition has no scored row")
    cond_pred = (sum_pred / counts[:, None]).astype(np.float64)
    cond_true = (sum_true / counts[:, None]).astype(np.float64)

    t0 = time.perf_counter()
    direction = metrics_mod.direction_auprc_matrix(cond_true, cond_pred, k)
    regression = metrics_mod.regression_metrics_matrix(cond_true, cond_pred)
    metrics_s = time.perf_counter() - t0

    macro_auprc = metrics_mod.macro(direction["auprc"])
    total_s = gather_s + forward_s + metrics_s
    return {
        "auprc": macro_auprc["mean"],
        "auprc_sd": macro_auprc["sd"],
        "auprc_up": metrics_mod.macro(direction["auprc_up"])["mean"],
        "auprc_down": metrics_mod.macro(direction["auprc_down"])["mean"],
        "rmse": metrics_mod.macro(regression["rmse"])["mean"],
        "pcc": metrics_mod.macro(regression["pcc"])["mean"],
        "spearman": metrics_mod.macro(regression["spearman"])["mean"],
        "valid_nll": nll_total / max(nll_seen, 1),
        "n_conditions": int(n_cond),
        "n_rows": int(rows.size),
        "timings": {
            "gather_s": gather_s,
            "forward_s": forward_s,
            "metrics_s": metrics_s,
            "total_s": total_s,
        },
    }


# --------------------------------------------------------------------------- #
# Offline evaluation of a prediction directory
# --------------------------------------------------------------------------- #
def _condition_index(condition_keys: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    unique, inverse = np.unique(np.asarray(condition_keys, dtype=object), return_inverse=True)
    return unique, inverse.astype(np.int64)


def expected_eval_rows(
    frozen: data_mod.FrozenInput, stage: str, role: str, fold: int | None
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """Return ``(rows, condition_ids, condition_names, info)`` for a scored split.

    Uses only metadata/control reads (plus the truth read, done by the caller).
    """
    structures = frozen.structures()
    if stage == "finetune":
        if fold is None:
            raise ValueError("--fold is required for --stage finetune")
        folds = data_mod.load_folds(frozen)
        source = data_mod.CMPredictSource(frozen, folds, int(fold), structures)
        rows = source.rows_for(role)
        keys = source.condition_key(rows)
        names, cond = _condition_index(keys)
        return rows, cond, names, {"fold": int(fold), "fold_table": folds.report}
    if stage == "pretrain":
        if role not in ("valid", "train"):
            raise ValueError("Tahoe has no test split; role must be 'train' or 'valid'")
        source = data_mod.TahoePredictSource(frozen, structures)
        rows = source.rows_for(role)
        dataset = source.dataset
        keys = np.asarray(
            [dataset._condition[dataset._row_pos[int(r)]] for r in rows], dtype=object
        )
        names, cond = _condition_index(keys)
        return rows, cond, names, {"pairs": dataset.counts.get("pairs_check")}
    raise ValueError(f"unknown stage {stage!r}")


def load_predictions(pred_dir: Path, frozen: data_mod.FrozenInput) -> dict[str, Any]:
    """Load and sanity-check a ``predict.py`` output directory."""
    pred_dir = Path(pred_dir)
    for name in ("pred_expr.npy", "delta.npy", "meta.csv", "manifest.json"):
        if not (pred_dir / name).exists():
            raise data_mod.DataError(f"prediction directory is missing {name}: {pred_dir}")
    pred_expr = np.load(pred_dir / "pred_expr.npy")
    delta = np.load(pred_dir / "delta.npy")
    meta = pd.read_csv(pred_dir / "meta.csv")
    if "row" not in meta.columns:
        raise data_mod.DataError("prediction meta.csv has no 'row' column")
    if pred_expr.shape != delta.shape:
        raise data_mod.DataError(
            f"pred_expr {pred_expr.shape} and delta {delta.shape} differ"
        )
    if pred_expr.shape[0] != len(meta):
        raise data_mod.DataError("prediction matrix rows do not match meta.csv")
    if pred_expr.shape[1] != frozen.n_genes:
        raise data_mod.DataError(
            f"predictions have {pred_expr.shape[1]} genes, frozen axis has {frozen.n_genes}"
        )
    if not np.all(np.isfinite(pred_expr)):
        raise data_mod.DataError("predictions contain non-finite values")
    return {
        "pred_expr": pred_expr,
        "delta": delta,
        "meta": meta,
        "manifest": manifest_mod.read_json(pred_dir / "manifest.json"),
    }


def _truth_for_rows(
    frozen: data_mod.FrozenInput,
    stage: str,
    rows: np.ndarray,
    control_pool: data_mod.CMControlPool | None,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Return ``(treated_log1p, matched_control_baseline, info)`` for the scored rows."""
    layer, sorted_rows = h5io.read_layer_rows(frozen.cm_h5ad if stage == "finetune" else frozen.tahoe_h5ad,
                                              rows, frozen.n_genes)
    if not np.array_equal(sorted_rows, np.asarray(rows, dtype=np.int64)):
        raise data_mod.DataError("truth row read returned an unexpected order")
    if stage == "finetune":
        if control_pool is None:
            control_pool = data_mod.CMControlPool(frozen)
        baseline = control_pool.mean_for_rows(rows)
        return layer, baseline, {
            "source": "CM per-plate matched DMSO mean (evaluation-only)",
            "n_control_wells": int(control_pool.control_rows.size),
            "n_plates": int(len(control_pool.plate_means)),
        }
    structures = frozen.structures()
    source = data_mod.TahoePredictSource(frozen, structures)
    baseline = source.baseline(rows)
    return layer, baseline, {
        "source": "Tahoe matched (plate, cell_line) control mean",
    }


def run_evaluation(options: EvaluateOptions) -> dict[str, Any]:
    """Score a prediction directory against the frozen truth; write the record."""
    started = time.perf_counter()
    started_utc = manifest_mod.now_utc()
    options.out_dir = Path(options.out_dir)
    options.out_dir.mkdir(parents=True, exist_ok=True)
    frozen = data_mod.open_frozen_input(
        options.data_dir, hash_cache=options.hash_cache, table_cache=options.table_cache
    )
    predictions = load_predictions(options.pred_dir, frozen)
    rows, condition_ids, condition_names, split_info = expected_eval_rows(
        frozen, options.stage, options.role, options.fold
    )
    meta = predictions["meta"]
    if len(meta) != rows.size or not np.array_equal(
        np.asarray(meta["row"], dtype=np.int64), rows
    ):
        raise data_mod.DataError(
            "prediction rows do not match the frozen split exactly (order included); "
            f"got {len(meta)} rows, expected {rows.size} — a prediction directory from a "
            "different fold/role/split cannot be scored"
        )

    control_pool = data_mod.CMControlPool(frozen) if options.stage == "finetune" else None
    treated, truth_baseline, truth_info = _truth_for_rows(
        frozen, options.stage, rows, control_pool
    )
    truth_delta = treated.astype(np.float64) - truth_baseline.astype(np.float64)
    pred_delta = np.asarray(predictions["delta"], dtype=np.float64)
    # The prediction's Δ must be exactly prediction - the input baseline the model saw.
    input_baseline = predictions["pred_expr"].astype(np.float64) - pred_delta
    if options.stage == "finetune":
        source = data_mod.CMPredictSource(
            frozen,
            data_mod.load_folds(frozen),
            int(options.fold),
            frozen.structures(),
        )
        expected_baseline = np.broadcast_to(source.c_star, (rows.size, frozen.n_genes))
        baseline_delta = float(np.max(np.abs(input_baseline - expected_baseline)))
        if baseline_delta > 1e-4:
            raise data_mod.DataError(
                f"predictions were not computed against this fold's c* (max abs diff "
                f"{baseline_delta:.3g}); refusing to score a mismatched baseline"
            )

    n_cond = len(condition_names)
    cond_pred = metrics_mod.aggregate_rows(pred_delta, condition_ids, n_cond)
    cond_true = metrics_mod.aggregate_rows(truth_delta, condition_ids, n_cond)

    k_values = [k for k in metrics_mod.K_LIST if k <= frozen.n_genes]
    if int(frozen.config["metric"]["k"]) not in k_values:
        k_values.append(int(frozen.config["metric"]["k"]))
    k_values = sorted(set(k_values))
    by_k: dict[int, dict[str, np.ndarray]] = {
        k: metrics_mod.direction_auprc_matrix(cond_true, cond_pred, k) for k in k_values
    }
    regression = metrics_mod.regression_metrics_matrix(cond_true, cond_pred)

    # per-condition table
    frames = {
        "condition": condition_names.astype(str),
        "n_replicates": np.zeros(n_cond, dtype=np.int64),
    }
    np.add.at(frames["n_replicates"], condition_ids, 1)
    per_condition = pd.DataFrame(frames)
    for k in k_values:
        per_condition[f"auprc_{k}"] = by_k[k]["auprc"]
        per_condition[f"auprc_up_{k}"] = by_k[k]["auprc_up"]
        per_condition[f"auprc_down_{k}"] = by_k[k]["auprc_down"]
        per_condition[f"n_up_{k}"] = by_k[k]["n_up"]
    per_condition["rmse"] = regression["rmse"]
    per_condition["pcc"] = regression["pcc"]
    per_condition["spearman"] = regression["spearman"]
    per_condition["pred_constant"] = regression["pred_constant"]
    per_condition["true_constant"] = regression["true_constant"]
    per_condition["regression_nan_reason"] = [
        "constant_prediction" if pc else ("constant_true" if tc else "")
        for pc, tc in zip(regression["pred_constant"], regression["true_constant"])
    ]
    per_condition.to_csv(options.out_dir / "per_condition.csv", index=False, float_format="%.10g")

    # long metrics.csv + aggregate.json
    metric_rows: list[dict[str, Any]] = []
    aggregates: dict[str, Any] = {
        "stage": options.stage,
        "role": options.role,
        "fold": options.fold,
        "n_conditions": int(n_cond),
        "n_rows": int(rows.size),
        "k_list": k_values,
        "headline_k": int(frozen.config["metric"]["k"]),
        "by_k": {},
        "regression": {},
    }
    for k in k_values:
        aggregates["by_k"][str(k)] = {}
        for direction_key, column in (
            ("macro", "auprc"),
            ("up", "auprc_up"),
            ("down", "auprc_down"),
        ):
            values = by_k[k][column]
            stats = metrics_mod.macro(values)
            aggregates["by_k"][str(k)][direction_key] = stats
            metric_rows.append(
                {
                    "level": "macro", "condition": "ALL", "k": k,
                    "direction": direction_key, "metric": "auprc",
                    "value": stats["mean"], "sd": stats["sd"], "nan": stats["nan"],
                }
            )
    for metric_name in ("rmse", "pcc", "spearman"):
        stats = metrics_mod.macro(regression[metric_name])
        aggregates["regression"][metric_name] = stats
        metric_rows.append(
            {
                "level": "macro", "condition": "ALL", "k": 0, "direction": "macro",
                "metric": metric_name, "value": stats["mean"], "sd": stats["sd"],
                "nan": stats["nan"],
            }
        )
    pd.DataFrame(metric_rows).to_csv(
        options.out_dir / "metrics.csv", index=False, float_format="%.10g"
    )
    manifest_mod.write_json(options.out_dir / "aggregate.json", aggregates)
    np.savez_compressed(
        options.out_dir / "deltas.npz",
        conditions=condition_names.astype(str),
        true_delta=cond_true.astype(np.float32),
        pred_delta=cond_pred.astype(np.float32),
        condition_ids=condition_ids.astype(np.int32),
        rows=np.asarray(rows, dtype=np.int64),
        k=np.asarray(k_values, dtype=np.int32),
    )

    pred_manifest = predictions["manifest"]
    manifest = manifest_mod.build_manifest(
        "evaluate",
        stage=options.stage,
        status="completed",
        pipeline_version=PIPELINE_VERSION,
        command=list(options.command or []),
        git=manifest_mod.git_info(),
        data=frozen.record(),
        fold=options.fold,
        role=options.role,
        predictions={
            "dir": str(Path(options.pred_dir)),
            "rows": int(rows.size),
            "pred_expr_sha256": manifest_mod.sha256_file(Path(options.pred_dir) / "pred_expr.npy"),
            "delta_sha256": manifest_mod.sha256_file(Path(options.pred_dir) / "delta.npy"),
            "checkpoint_sha256": (pred_manifest.get("checkpoint") or {}).get("sha256"),
        },
        metric={
            "name": "direction_auprc",
            "k_list": k_values,
            "headline_k": int(frozen.config["metric"]["k"]),
            "aggregation": "replicates averaged within condition, then direction AP, macro over conditions",
            "plateau_mode": "max",
        },
        truth=truth_info,
        timings={
            "started_utc": started_utc,
            "finished_utc": manifest_mod.now_utc(),
            "wall_s": time.perf_counter() - started,
        },
        peak_rss_mb=manifest_mod.peak_rss_mb(),
        noise_policy={"inference": "zero noise, model output read from predictions", "validation": "n/a"},
        warnings=[],
    )
    manifest_mod.write_json(options.out_dir / "manifest.json", manifest)
    log.info(
        "evaluated %s/%s fold=%s: %d rows, %d conditions, headline AUPRC@%d=%.6f",
        options.stage, options.role, options.fold, rows.size, n_cond,
        int(frozen.config["metric"]["k"]),
        aggregates["by_k"][str(int(frozen.config["metric"]["k"]))]["macro"]["mean"],
    )
    return {"aggregate": aggregates, "out_dir": str(options.out_dir)}
