"""Aggregation for ticket #9 — every number derives from the evaluator's frozen output.

The evaluator (`ptlib.evaluate`) writes, per unit, a ``per_condition.csv`` (condition =
``drug_key|dose_µM``, direction AUPRC@200 and auxiliary Δ metrics), an ``aggregate.json``
(macro statistics over conditions) and ``deltas.npz`` (condition-level true/predicted Δ).
This module only *aggregates* those files:

* per unit, per fold×seed, per drug×dose;
* the paired Tahoe-vs-random gain **by drug** — doses, seeds and candidate runs are never
  treated as independent drugs;
* the paired comparison with #8's primary ``delta_mean`` baseline on the identical test
  population;
* dose / Tahoe-seen / low-response / toxicity strata (labels that do not exist in this
  round are reported as ``unknown``, never inferred from transcriptomic magnitude);
* paired bootstrap uncertainty resampling **drugs** (the independent sampling unit).

No metric is recomputed from expression here: the direction AUPRC and auxiliary Δ metrics
come from the frozen evaluator's own output. The only recomputation is a consistency check
of the macro mean against the evaluator's ``aggregate.json``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

import common

__all__ = [
    "K",
    "BOOTSTRAP_SEED",
    "N_BOOTSTRAP",
    "UnitEval",
    "load_unit_eval",
    "unit_results_frame",
    "fold_seed_frame",
    "per_drug_dose_frame",
    "paired_by_drug_frame",
    "drug_gains_frame",
    "bootstrap_mean",
    "paired_summary",
    "strata_frame",
    "load_baseline_conditions",
    "load_baselines_summary",
    "read_tahoe_seen_drug_keys",
    "condition_parts",
]

K = 200
N_BOOTSTRAP = 10000
BOOTSTRAP_SEED = 20261102
MODEL_INITS = ("random", "tahoe")
PRIMARY_BASELINE = "delta_mean"
BASELINE_ORDER = ("no_change", "delta_mean", "delta_mean_dose", "ridge_fp")


class AggregateError(common.EvaluationError):
    """A unit's evaluation output cannot be aggregated."""


def condition_parts(condition: str) -> tuple[str, float]:
    """``drug_key|dose`` -> ``(drug_key, dose)`` (the evaluator's CM condition key)."""
    text = str(condition)
    if "|" not in text:
        raise AggregateError(f"not a CM condition key: {condition!r}")
    drug, dose = text.rsplit("|", 1)
    return drug, float(dose)


@dataclass
class UnitEval:
    """One unit's frozen evaluation output."""

    unit: dict[str, Any]
    aggregate: dict[str, Any]
    per_condition: pd.DataFrame
    condition_keys: np.ndarray
    rows: np.ndarray
    true_delta: np.ndarray
    pred_delta: np.ndarray
    eval_manifest: dict[str, Any]
    eval_dir: Path

    @property
    def unit_id(self) -> str:
        return str(self.unit["unit_id"])

    @property
    def fold(self) -> int:
        return int(self.unit["fold"])

    @property
    def seed(self) -> int:
        return int(self.unit["seed"])

    @property
    def init(self) -> str:
        return str(self.unit["init_name"])


def load_unit_eval(unit: Mapping[str, Any], out: Path | str) -> UnitEval:
    """Load + sanity-check one unit's evaluator output (never re-derives a metric)."""
    out = Path(out)
    unit_id = str(unit["unit_id"])
    eval_dir = out / "units" / unit_id / "eval"
    for name in ("aggregate.json", "per_condition.csv", "deltas.npz", "manifest.json"):
        if not (eval_dir / name).exists():
            raise AggregateError(
                f"unit {unit_id} has no evaluation output ({eval_dir / name} missing); "
                "run score.py first"
            )
    aggregate = common.load_json(eval_dir / "aggregate.json")
    manifest = common.load_json(eval_dir / "manifest.json")
    raw_frame = common.read_table(
        eval_dir / "per_condition.csv",
        string_columns=("condition", "regression_nan_reason"),
    )
    numeric = [
        column for column in raw_frame.columns
        if column not in ("condition", "regression_nan_reason")
    ]
    frame = raw_frame.copy()
    for column in numeric:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    payload = np.load(eval_dir / "deltas.npz")
    condition_keys = np.asarray(payload["conditions"], dtype=str)
    rows = np.asarray(payload["rows"], dtype=np.int64)
    true_delta = np.asarray(payload["true_delta"], dtype=np.float64)
    pred_delta = np.asarray(payload["pred_delta"], dtype=np.float64)

    if str(aggregate.get("stage")) != "finetune":
        raise AggregateError(f"unit {unit_id}: aggregate stage is not finetune")
    if int(aggregate.get("fold", -1)) != int(unit["fold"]):
        raise AggregateError(f"unit {unit_id}: aggregate fold does not match the unit")
    if int(aggregate.get("headline_k", -1)) != K:
        raise AggregateError(f"unit {unit_id}: headline k is not {K}")
    if int(aggregate.get("n_rows", -1)) != rows.size:
        raise AggregateError(f"unit {unit_id}: aggregate row count disagrees with deltas.npz")
    if len(frame) != condition_keys.size:
        raise AggregateError(f"unit {unit_id}: per_condition rows disagree with deltas.npz")
    if not np.array_equal(
        np.asarray(frame["condition"].astype(str), dtype=str), condition_keys
    ):
        raise AggregateError(f"unit {unit_id}: per_condition condition order disagrees")
    k_column = f"auprc_{K}"
    if k_column not in frame.columns:
        raise AggregateError(f"unit {unit_id}: per_condition.csv has no {k_column}")
    # macro consistency: the evaluator's aggregate == mean over its own per-condition column
    macro = aggregate["by_k"][str(K)]["macro"]
    observed = float(np.nanmean(frame[k_column].to_numpy(dtype=np.float64)))
    if not math.isclose(float(macro["mean"]), observed, rel_tol=1e-6, abs_tol=1e-9):
        raise AggregateError(
            f"unit {unit_id}: aggregate macro AUPRC@{K} {macro['mean']} != "
            f"per-condition mean {observed}"
        )
    # direction-set sizes recorded by #8's schema; recompute the down set from frozen truth
    n_pos = (true_delta > 0).sum(axis=1)
    n_neg = (true_delta < 0).sum(axis=1)
    expected_up = np.minimum(n_pos, K)
    expected_down = np.minimum(n_neg, K)
    if int(frame[f"n_up_{K}"].to_numpy(dtype=np.int64).sum()) != int(expected_up.sum()):
        raise AggregateError(f"unit {unit_id}: per_condition up-set sizes disagree with truth Δ")
    frame[f"n_down_{K}"] = expected_down.astype(np.int64)
    return UnitEval(
        unit=dict(unit),
        aggregate=aggregate,
        per_condition=frame,
        condition_keys=condition_keys,
        rows=rows,
        true_delta=true_delta,
        pred_delta=pred_delta,
        eval_manifest=manifest,
        eval_dir=eval_dir,
    )


# --------------------------------------------------------------------------- #
# Unit / fold / seed tables
# --------------------------------------------------------------------------- #
def _macro(aggregate: Mapping[str, Any], direction: str) -> dict[str, Any]:
    return dict(aggregate["by_k"][str(K)][direction])


def unit_results_frame(evals: Sequence[UnitEval], score_records: Mapping[str, Mapping[str, Any]]) -> pd.DataFrame:
    """One row per evaluated unit: header metrics, auxiliary Δ metrics, training budget."""
    rows: list[dict[str, Any]] = []
    for unit_eval in evals:
        unit = unit_eval.unit
        aggregate = unit_eval.aggregate
        macro = _macro(aggregate, "macro")
        up = _macro(aggregate, "up")
        down = _macro(aggregate, "down")
        regression = aggregate.get("regression") or {}
        name = lambda key: (regression.get(key) or {}).get("mean")  # noqa: E731
        training = unit.get("training") or {}
        score = score_records.get(unit_eval.unit_id) or {}
        steps = score.get("steps") or {}
        rows.append(
            {
                "unit_id": unit_eval.unit_id,
                "init": unit_eval.init,
                "fold": unit_eval.fold,
                "seed": unit_eval.seed,
                "seed_index": unit.get("seed_index"),
                "lr": float(unit["lr"]),
                "checkpoint_sha256": unit["checkpoint"]["sha256"],
                "config_sha256": unit["config_sha256"],
                "data_identity": unit["data_identity"],
                "n_rows": int(aggregate["n_rows"]),
                "n_conditions": int(aggregate["n_conditions"]),
                "auprc_200": float(macro["mean"]),
                "auprc_sd_200": macro.get("sd"),
                "auprc_up_200": float(up["mean"]),
                "auprc_down_200": float(down["mean"]),
                "auprc_nan_conditions": int(macro.get("nan") or 0),
                "rmse": name("rmse"),
                "pcc": name("pcc"),
                "spearman": name("spearman"),
                "best_valid_auprc": unit.get("best_valid_auprc"),
                "best_epoch": training.get("best_epoch"),
                "epochs_run": training.get("epochs_run"),
                "stop_reason": training.get("stop_reason"),
                "train_wall_s": training.get("train_wall_s"),
                "train_peak_rss_mb": training.get("train_peak_rss_mb"),
                "eval_wall_s": (
                    float(steps.get("predict", {}).get("wall_s", 0.0))
                    + float(steps.get("evaluate", {}).get("wall_s", 0.0))
                )
                if steps
                else None,
                "eval_host_peak_rss_mb": max(
                    (
                        float(steps[step].get("peak_rss_mb", 0.0))
                        for step in ("predict", "evaluate")
                        if steps.get(step)
                    ),
                    default=None,
                ),
                "eval_gpu_peak_mb": (score.get("gpu") or {}).get("peak_mb"),
                "status": score.get("status", "unknown"),
            }
        )
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(["init", "fold", "seed"]).reset_index(drop=True)


def fold_seed_frame(unit_results: pd.DataFrame) -> pd.DataFrame:
    """Per fold × init: three-seed mean/sd for the headline and auxiliary metrics."""
    if unit_results.empty:
        return unit_results
    metric_columns = [
        "auprc_200", "auprc_up_200", "auprc_down_200", "rmse", "pcc", "spearman",
        "best_valid_auprc",
    ]
    grouped = unit_results.groupby(["fold", "init"], dropna=False)
    rows: list[dict[str, Any]] = []
    for (fold, init), group in grouped:
        row: dict[str, Any] = {
            "fold": int(fold),
            "init": str(init),
            "n_seeds": int(group["seed"].nunique()),
            "seeds": ",".join(str(int(s)) for s in sorted(group["seed"].unique())),
            "n_rows": int(group["n_rows"].iloc[0]),
            "n_conditions": int(group["n_conditions"].iloc[0]),
        }
        for column in metric_columns:
            values = pd.to_numeric(group[column], errors="coerce").to_numpy(dtype=np.float64)
            row[f"{column}_mean"] = float(np.nanmean(values)) if np.any(np.isfinite(values)) else np.nan
            row[f"{column}_sd"] = (
                float(np.nanstd(values, ddof=1)) if np.sum(np.isfinite(values)) > 1 else np.nan
            )
        rows.append(row)
    frame = pd.DataFrame(rows)
    return frame.sort_values(["fold", "init"]).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Baselines (#8's tracked record)
# --------------------------------------------------------------------------- #
def load_baseline_conditions(
    baselines_dir: Path | str, *, required_folds: Iterable[int] | None = None
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Load #8's per-condition record and verify it covers every fold."""
    baselines_dir = Path(baselines_dir)
    path = baselines_dir / "per_drug_dose.csv"
    if not path.exists():
        raise AggregateError(f"#8 per-condition record is missing: {path}")
    frame = common.read_table(
        path,
        numeric_columns=(
            "fold", "dose_um", "n_replicates", f"auprc_{K}", f"auprc_up_{K}",
            f"auprc_down_{K}", f"n_up_{K}", f"n_down_{K}", "rmse", "pcc", "spearman",
        ),
        string_columns=("baseline", "drug_key", "condition", "regression_nan_reason"),
    )
    required = ("fold", "baseline", "drug_key", "dose_um", "condition", f"auprc_{K}")
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise AggregateError(f"{path} is missing columns {missing}")
    folds = sorted(int(value) for value in frame["fold"].unique())
    if required_folds is not None and set(folds) != {int(v) for v in required_folds}:
        raise AggregateError(
            f"#8's record covers folds {folds}, expected {sorted(int(v) for v in required_folds)}"
        )
    record = {
        "path": str(path),
        "sha256": common.sha256_file(path),
        "n_rows": int(len(frame)),
        "folds": folds,
        "baselines": sorted(frame["baseline"].astype(str).unique()),
        "n_conditions_per_fold": {
            str(fold): int(len(frame[(frame["fold"] == fold) & (frame["baseline"] == PRIMARY_BASELINE)]))
            for fold in folds
        },
    }
    return frame, record


def load_baselines_summary(baselines_dir: Path | str) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Load #8's summary table (macro-over-folds rows are the published bar)."""
    path = Path(baselines_dir) / "summary_table.csv"
    if not path.exists():
        raise AggregateError(f"#8 summary table is missing: {path}")
    frame = common.read_table(
        path,
        numeric_columns=(
            "auprc", "auprc_up", "auprc_down", "rmse", "pcc", "spearman",
        ),
        string_columns=("scope", "fold", "baseline", "role"),
    )
    record = {"path": str(path), "sha256": common.sha256_file(path)}
    return frame, record


def _baseline_pivot(baselines: pd.DataFrame, fold: int) -> pd.DataFrame:
    """Wide per-condition AUPRC@200 for one fold (one column per baseline)."""
    subset = baselines[baselines["fold"] == int(fold)]
    if subset.empty:
        raise AggregateError(f"#8's record has no conditions for fold {fold}")
    pivot = subset.pivot_table(
        index="condition", columns="baseline", values=f"auprc_{K}", aggfunc="first"
    )
    pivot.columns = [f"{str(column)}_auprc_200" for column in pivot.columns]
    return pivot.reset_index()


# --------------------------------------------------------------------------- #
# Per drug x dose
# --------------------------------------------------------------------------- #
def per_drug_dose_frame(
    evals: Sequence[UnitEval], baselines: pd.DataFrame | None = None
) -> pd.DataFrame:
    """Per unit × condition predictions/metrics, joined to #8's identical-population bars."""
    frames: list[pd.DataFrame] = []
    for unit_eval in evals:
        frame = unit_eval.per_condition.copy()
        frame.insert(0, "unit_id", unit_eval.unit_id)
        frame.insert(1, "init", unit_eval.init)
        frame.insert(2, "fold", unit_eval.fold)
        frame.insert(3, "seed", unit_eval.seed)
        frame.insert(4, "lr", float(unit_eval.unit["lr"]))
        drugs = []
        doses = []
        for condition in frame["condition"].astype(str):
            drug, dose = condition_parts(condition)
            drugs.append(drug)
            doses.append(dose)
        frame.insert(5, "drug_key", drugs)
        frame.insert(6, "dose_um", doses)
        if baselines is not None:
            pivot = _baseline_pivot(baselines, unit_eval.fold)
            frame = frame.merge(pivot, on="condition", how="left", validate="one_to_one")
            if frame[f"{PRIMARY_BASELINE}_auprc_200"].isna().any():
                missing = frame.loc[
                    frame[f"{PRIMARY_BASELINE}_auprc_200"].isna(), "condition"
                ].tolist()[:5]
                raise AggregateError(
                    f"unit {unit_eval.unit_id}: #8 has no delta_mean row for conditions "
                    f"{missing} (test populations differ)"
                )
            reference_conditions = set(pivot["condition"].astype(str))
            model_conditions = set(frame["condition"].astype(str))
            if model_conditions != reference_conditions:
                raise AggregateError(
                    f"fold {unit_eval.fold}: the model's scored conditions differ from "
                    f"#8's identical test population "
                    f"(+{sorted(model_conditions - reference_conditions)[:3]}, "
                    f"-{sorted(reference_conditions - model_conditions)[:3]})"
                )
        frames.append(frame)
    if not frames:
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    columns = [
        "unit_id", "init", "fold", "seed", "lr", "drug_key", "dose_um", "condition",
        "n_replicates", f"auprc_{K}", f"auprc_up_{K}", f"auprc_down_{K}",
        f"n_up_{K}", f"n_down_{K}", "rmse", "pcc", "spearman", "pred_constant",
        "true_constant", "regression_nan_reason",
    ]
    columns += [column for column in combined.columns if column.endswith("_auprc_200")
                and column not in columns]
    return combined[columns].reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Paired by drug (the resampling unit is the drug)
# --------------------------------------------------------------------------- #
def paired_by_drug_frame(per_drug_dose: pd.DataFrame) -> pd.DataFrame:
    """Per (fold, seed, drug): init arms and #8 baselines, all on the same conditions."""
    if per_drug_dose.empty:
        return per_drug_dose
    value = f"auprc_{K}"
    baseline_columns = [
        column for column in per_drug_dose.columns if column.endswith("_auprc_200")
    ]
    # drug-level = equal-weight mean over the drug's dose conditions, within each unit
    grouped = per_drug_dose.groupby(
        ["fold", "seed", "init", "drug_key"], dropna=False
    )
    drug_level = grouped[[value, *baseline_columns]].mean().reset_index()
    counts = grouped.size().reset_index(name="n_doses")
    drug_level = drug_level.merge(
        counts, on=["fold", "seed", "init", "drug_key"], how="left", validate="one_to_one"
    )

    random = drug_level[drug_level["init"] == "random"].rename(
        columns={value: "random_auprc_200", "n_doses": "random_n_doses"}
    )
    tahoe = drug_level[drug_level["init"] == "tahoe"].rename(
        columns={value: "tahoe_auprc_200", "n_doses": "tahoe_n_doses"}
    )
    for frame in (random, tahoe):
        for column in baseline_columns:
            frame.drop(columns=[column], inplace=True, errors="ignore")
    paired = random.merge(
        tahoe, on=["fold", "seed", "drug_key"], how="outer", suffixes=("", "_tahoe")
    )
    paired["n_doses"] = paired["random_n_doses"].fillna(paired["tahoe_n_doses"])
    if paired[["random_auprc_200", "tahoe_auprc_200"]].isna().any().any():
        incomplete = paired[
            paired[["random_auprc_200", "tahoe_auprc_200"]].isna().any(axis=1)
        ]
        raise AggregateError(
            f"paired init comparison incomplete for {len(incomplete)} (fold, seed, drug) "
            f"rows (e.g. {incomplete.head(2).to_dict('records')})"
        )
    # baselines are init-independent: re-join the per-fold values from the same frame
    baseline_level = (
        drug_level.groupby(["fold", "drug_key"], dropna=False)[baseline_columns]
        .mean()
        .reset_index()
    )
    paired = paired.drop(
        columns=[c for c in baseline_columns if c in paired.columns], errors="ignore"
    ).merge(baseline_level, on=["fold", "drug_key"], how="left", validate="many_to_one")
    paired["gain_tahoe_minus_random"] = (
        paired["tahoe_auprc_200"] - paired["random_auprc_200"]
    )
    paired[f"gain_tahoe_minus_{PRIMARY_BASELINE}"] = (
        paired["tahoe_auprc_200"] - paired[f"{PRIMARY_BASELINE}_auprc_200"]
    )
    paired[f"gain_random_minus_{PRIMARY_BASELINE}"] = (
        paired["random_auprc_200"] - paired[f"{PRIMARY_BASELINE}_auprc_200"]
    )
    paired["gain_tahoe_minus_no_change"] = (
        paired["tahoe_auprc_200"] - paired.get("no_change_auprc_200", np.nan)
    )
    paired["gain_random_minus_no_change"] = (
        paired["random_auprc_200"] - paired.get("no_change_auprc_200", np.nan)
    )
    columns = [
        "fold", "seed", "drug_key", "n_doses",
        "random_auprc_200", "tahoe_auprc_200", "gain_tahoe_minus_random",
        f"{PRIMARY_BASELINE}_auprc_200", f"gain_tahoe_minus_{PRIMARY_BASELINE}",
        f"gain_random_minus_{PRIMARY_BASELINE}",
        "no_change_auprc_200", "gain_tahoe_minus_no_change", "gain_random_minus_no_change",
        "delta_mean_dose_auprc_200", "ridge_fp_auprc_200",
    ]
    columns = [column for column in columns if column in paired.columns]
    return paired[columns].sort_values(["fold", "seed", "drug_key"]).reset_index(drop=True)


def drug_gains_frame(paired: pd.DataFrame) -> pd.DataFrame:
    """One row per drug: seed-mean arms/gains and their seed spread (bootstrap input)."""
    if paired.empty:
        return paired
    value_columns = [
        "random_auprc_200", "tahoe_auprc_200", "gain_tahoe_minus_random",
        f"{PRIMARY_BASELINE}_auprc_200", f"gain_tahoe_minus_{PRIMARY_BASELINE}",
        f"gain_random_minus_{PRIMARY_BASELINE}",
        "no_change_auprc_200", "gain_tahoe_minus_no_change", "gain_random_minus_no_change",
    ]
    value_columns = [column for column in value_columns if column in paired.columns]
    rows: list[dict[str, Any]] = []
    for (fold, drug), group in paired.groupby(["fold", "drug_key"], dropna=False):
        row: dict[str, Any] = {
            "fold": int(fold),
            "drug_key": str(drug),
            "n_seeds": int(group["seed"].nunique()),
        }
        for column in value_columns:
            values = pd.to_numeric(group[column], errors="coerce").to_numpy(dtype=np.float64)
            row[f"{column}_mean"] = float(np.nanmean(values))
            row[f"{column}_seed_sd"] = (
                float(np.nanstd(values, ddof=1)) if np.sum(np.isfinite(values)) > 1 else np.nan
            )
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["fold", "drug_key"]).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Bootstrap (drugs are the independent sampling unit)
# --------------------------------------------------------------------------- #
def bootstrap_mean(
    values: Sequence[float],
    *,
    n_bootstrap: int = N_BOOTSTRAP,
    seed: int = BOOTSTRAP_SEED,
) -> dict[str, Any]:
    """Mean with a 95 % percentile bootstrap CI, resampling the given (drug) values."""
    array = np.asarray([value for value in values if _finite(value)], dtype=np.float64)
    result: dict[str, Any] = {
        "mean": None, "ci_low": None, "ci_high": None,
        "n": int(array.size), "n_bootstrap": int(n_bootstrap), "seed": int(seed),
        "sd": None, "min": None, "max": None,
    }
    if array.size == 0:
        return result
    result["mean"] = float(array.mean())
    result["sd"] = float(array.std(ddof=1)) if array.size > 1 else 0.0
    result["min"] = float(array.min())
    result["max"] = float(array.max())
    if array.size == 1 or n_bootstrap <= 0:
        result["ci_low"] = result["ci_high"] = float(array.mean())
        return result
    rng = np.random.default_rng(int(seed))
    indices = rng.integers(0, array.size, size=(int(n_bootstrap), array.size))
    means = array[indices].mean(axis=1)
    result["ci_low"] = float(np.percentile(means, 2.5))
    result["ci_high"] = float(np.percentile(means, 97.5))
    return result


def _finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _verdict(stats: Mapping[str, Any]) -> str:
    """Pre-registered reading of the paired evidence (recorded in the report)."""
    mean, low = stats.get("mean"), stats.get("ci_low")
    if mean is None:
        return "not_computable"
    if low is not None and float(low) > 0:
        return "positive"
    if float(mean) <= 0:
        return "negative"
    return "inconclusive"


def paired_summary(
    drug_gains: pd.DataFrame,
    gain_key: str,
    *,
    n_bootstrap: int = N_BOOTSTRAP,
    seed: int = BOOTSTRAP_SEED,
) -> dict[str, Any]:
    """Paired gain summary: overall CI + per-fold + seed stability (all by drug)."""
    column = f"{gain_key}_mean"
    if drug_gains.empty or column not in drug_gains.columns:
        return {"gain": gain_key, "n_drugs": 0, "overall": bootstrap_mean([], n_bootstrap=n_bootstrap, seed=seed)}
    stats = bootstrap_mean(drug_gains[column], n_bootstrap=n_bootstrap, seed=seed)
    per_fold: dict[str, Any] = {}
    for fold, group in drug_gains.groupby("fold"):
        values = group[column].to_numpy(dtype=np.float64)
        per_fold[str(int(fold))] = {
            "mean": float(np.nanmean(values)),
            "n_drugs": int(values.size),
        }
    # seed stability: mean over drugs of each seed's paired value, then mean/sd over seeds
    seed_means = (
        drug_gains.groupby("fold")[column].mean().to_dict() if not drug_gains.empty else {}
    )
    # per-seed structure lives in paired_by_drug; recompute a seed-level mean across folds
    stability = None
    if "n_seeds" in drug_gains.columns:
        stability = {
            "n_seeds": int(drug_gains["n_seeds"].max()),
            "per_drug_seed_sd_mean": float(
                np.nanmean(drug_gains.get(f"{gain_key}_seed_sd", pd.Series(dtype=float)))
            )
            if f"{gain_key}_seed_sd" in drug_gains.columns
            else None,
        }
    return {
        "gain": gain_key,
        "metric": f"direction AUPRC@{K} (macro over conditions), drug-level",
        "n_drugs": int(len(drug_gains)),
        "n_bootstrap": int(n_bootstrap),
        "seed": int(seed),
        "overall": stats,
        f"per_fold_{gain_key}": per_fold,
        "fold_mean_of_drugs": seed_means,
        "seed_stability": stability,
        "verdict": _verdict(stats),
    }


# --------------------------------------------------------------------------- #
# Strata
# --------------------------------------------------------------------------- #
def read_tahoe_seen_drug_keys(data_dir: Path | str) -> tuple[set[str], dict[str, Any]]:
    """Drug keys treated in Tahoe (a *seen structure*, never a CM test label)."""
    data_dir = Path(data_dir)
    tahoe = data_dir / "tahoe_transfer.h5ad"
    if not tahoe.exists():
        raise AggregateError(f"frozen Tahoe view is missing: {tahoe}")
    common.pipeline_modules()
    from ptlib import h5io  # noqa: WPS433 - pipeline reader, no second parser

    meta = h5io.open_obs_var(tahoe)
    obs = meta.obs
    treated = ~np.asarray(obs["is_control"]).astype(bool)
    if "has_structure" in obs.columns:
        treated &= np.asarray(obs["has_structure"]).astype(bool)
    keys = obs.loc[treated, "drug_key"].astype(str)
    seen = sorted({key for key in keys.tolist() if key and key.lower() != "nan"})
    record = {
        "path": str(tahoe),
        "sha256": common.sha256_file(tahoe),
        "n_tahoe_observed_treated_rows": int(treated.sum()),
        "n_tahoe_seen_drug_keys": len(seen),
        "drug_keys": seen,
    }
    return set(seen), record


def strata_frame(
    per_drug_dose: pd.DataFrame,
    drug_gains: pd.DataFrame,
    *,
    tahoe_seen: set[str],
    n_bootstrap: int = N_BOOTSTRAP,
    seed: int = BOOTSTRAP_SEED,
) -> pd.DataFrame:
    """Dose / Tahoe-seen / unknown-label strata, all paired by drug.

    ``response_label`` and ``toxicity_label`` are ``unknown`` on every row of this round
    (06.2 §6: no per-well viability/response readout exists and it must never be inferred
    from transcriptomic magnitude). They are reported as an explicit ``unknown`` stratum
    with the pre-registered criterion, never silently dropped.
    """
    if per_drug_dose.empty:
        return pd.DataFrame()
    model = per_drug_dose.copy()
    model["tahoe_seen"] = model["drug_key"].astype(str).isin(tahoe_seen)
    value = f"auprc_{K}"
    rows: list[dict[str, Any]] = []

    def _add(
        stratum_type: str,
        stratum_value: str,
        subset: pd.DataFrame,
        drug_subset: pd.DataFrame,
        *,
        label_known: bool,
        note: str,
        gain_key: str = "gain_tahoe_minus_random",
        ci_key: str = "gain_tahoe_minus_random_mean",
    ) -> None:
        arm_means = (
            subset.groupby(["fold", "init"])[value].mean().unstack("init")
            if not subset.empty
            else pd.DataFrame()
        )
        random_auprc = (
            float(np.nanmean(arm_means["random"])) if "random" in arm_means else None
        )
        tahoe_auprc = (
            float(np.nanmean(arm_means["tahoe"])) if "tahoe" in arm_means else None
        )
        if ci_key in drug_subset.columns:
            stats = bootstrap_mean(
                drug_subset[ci_key], n_bootstrap=n_bootstrap, seed=seed
            )
        else:
            stats = bootstrap_mean([])
        delta_mean = None
        if f"{PRIMARY_BASELINE}_auprc_200_mean" in drug_subset.columns:
            delta_mean = float(np.nanmean(drug_subset[f"{PRIMARY_BASELINE}_auprc_200_mean"]))
        rows.append(
            {
                "stratum_type": stratum_type,
                "stratum_value": stratum_value,
                "label_known": bool(label_known),
                "n_conditions": int(subset["condition"].nunique()) if not subset.empty else 0,
                "n_drugs": int(drug_subset["drug_key"].nunique()) if not drug_subset.empty else 0,
                "n_seeds": int(subset["seed"].nunique()) if not subset.empty else 0,
                "random_auprc_200": random_auprc,
                "tahoe_auprc_200": tahoe_auprc,
                "delta_mean_auprc_200": delta_mean,
                "gain_tahoe_minus_random": stats.get("mean"),
                "gain_tahoe_minus_random_ci_low": stats.get("ci_low"),
                "gain_tahoe_minus_random_ci_high": stats.get("ci_high"),
                "gain_tahoe_minus_delta_mean": (
                    float(
                        np.nanmean(
                            drug_subset[f"gain_tahoe_minus_{PRIMARY_BASELINE}_mean"]
                        )
                    )
                    if f"gain_tahoe_minus_{PRIMARY_BASELINE}_mean" in drug_subset.columns
                    else None
                ),
                "n_bootstrap": int(n_bootstrap) if stats.get("n") else 0,
                "note": note,
            }
        )

    for dose in sorted(model["dose_um"].dropna().unique()):
        subset = model[model["dose_um"] == dose]
        drug_subset = drug_gains.merge(
            subset[["fold", "drug_key"]].drop_duplicates(), on=["fold", "drug_key"], how="inner"
        )
        _add(
            "dose",
            f"{float(dose):g}",
            subset,
            drug_subset,
            label_known=True,
            note="pre-registered dose stratum (measured temp_env doses)",
        )
    for seen_value, note in (
        (True, "drug_key observed in Tahoe treated rows (structure seen; no CM label)"),
        (False, "drug_key never observed in Tahoe treated rows"),
    ):
        subset = model[model["tahoe_seen"] == seen_value]
        drug_subset = drug_gains.merge(
            subset[["fold", "drug_key"]].drop_duplicates(), on=["fold", "drug_key"], how="inner"
        )
        _add(
            "tahoe_seen",
            str(seen_value),
            subset,
            drug_subset,
            label_known=True,
            note=note,
        )
    for label_type in ("response_label", "toxicity_label"):
        rows.append(
            {
                "stratum_type": label_type,
                "stratum_value": "unknown",
                "label_known": False,
                "n_conditions": int(model["condition"].nunique()),
                "n_drugs": int(model["drug_key"].nunique()),
                "n_seeds": int(model["seed"].nunique()),
                "random_auprc_200": None,
                "tahoe_auprc_200": None,
                "delta_mean_auprc_200": None,
                "gain_tahoe_minus_random": None,
                "gain_tahoe_minus_random_ci_low": None,
                "gain_tahoe_minus_random_ci_high": None,
                "gain_tahoe_minus_delta_mean": None,
                "n_bootstrap": 0,
                "note": (
                    "UNKNOWN by pre-registered criterion: no per-well response/toxicity "
                    "readout exists (06.2 §6) and it is never inferred from transcriptomic "
                    "magnitude; no labelled stratum is claimed this round"
                ),
            }
        )
    return pd.DataFrame(rows)
