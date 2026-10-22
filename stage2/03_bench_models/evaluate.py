"""Evaluation path and leaderboard for the Stage 2 benchmark (ticket 07/12).

The evaluator reads per-fold `pred_test.npy` / `group_ids.json`, averages the
per-well predictions within `group_id`, scores the group mean, and writes:

    results/metrics.csv
    results/per_drug.csv
    results/leaderboard.csv
    results/bootstrap.npz
    results/README.md

The metric definitions live in `metrics.py`; this module owns aggregation and
artifact layout.  It never touches model internals.
"""

from __future__ import annotations

import argparse
import json
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from cases import BenchData, save_npz
from metrics import group_metrics

K_LIST = (100, 200, 500)
BOOTSTRAP_REPS = 2_000
BOOTSTRAP_SEED = 20261022

METRIC_COLUMNS = [
    "model",
    "fold",
    "level",
    "subset",
    "group_id",
    "drug_key",
    "dose_um",
    "k",
    "direction",
    "metric",
    "value",
    "nan_reason",
]
PER_DRUG_COLUMNS = [
    "model",
    "fold",
    "drug_key",
    "n_groups",
    "auprc",
    "pcc",
    "spearman",
    "r2",
    "rmse",
]
LEADERBOARD_COLUMNS = [
    "rank",
    "model",
    "headline_auprc",
    "headline_sd",
    "ci_low",
    "ci_high",
    "fold0",
    "fold1",
    "fold2",
    "fold3",
    "fold4",
    "mean_k100",
    "mean_k500",
    "pcc",
    "spearman",
    "r2",
    "rmse",
    "pooled_rmse",
    "nrmse",
    "nan_counts",
    "statistically_tied",
    "ideal_fit",
    "n_folds",
    "winner_eligible",
    "status",
]


def _mean(values) -> float:
    vals = [float(v) for v in values if v is not None and not math.isnan(float(v))]
    return float(np.mean(vals)) if vals else math.nan


def _sd(values) -> float:
    vals = [float(v) for v in values if v is not None and not math.isnan(float(v))]
    if len(vals) < 2:
        return math.nan
    return float(np.std(vals, ddof=1))


def _min(values) -> float:
    vals = [float(v) for v in values if v is not None and not math.isnan(float(v))]
    return float(np.min(vals)) if vals else math.nan


def _nan_count(values) -> int:
    return int(sum(1 for v in values if v is None or math.isnan(float(v))))


def _ideal_fit(value) -> str:
    """Normalise the registry's ideal-form label to `full` / `partial` / `none`."""

    if isinstance(value, str):
        return value
    return "full" if bool(value) else "none"


def group_mean_predictions(pred: np.ndarray, group_ids) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Average per-well predictions within `group_id`.

    Returns `(group_ids_sorted, group_pred float64, counts)`.
    """

    group_ids = np.asarray(group_ids).astype(str)
    if pred.shape[0] != len(group_ids):
        raise ValueError(f"pred rows {pred.shape[0]} != group ids {len(group_ids)}")
    unique, inverse = np.unique(group_ids, return_inverse=True)
    sums = np.zeros((len(unique), pred.shape[1]), dtype=np.float64)
    np.add.at(sums, inverse, pred.astype(np.float64))
    counts = np.bincount(inverse, minlength=len(unique)).astype(np.int64)
    if (counts == 0).any():
        raise ValueError("empty group in prediction averaging")
    return unique, sums / counts[:, None], counts


def _nan_reason(values, pred: np.ndarray, true: np.ndarray, metric: str) -> str:
    if not math.isnan(float(values)):
        return ""
    if metric in {"pcc", "spearman"} and float(np.std(pred)) == 0.0:
        return "constant_prediction"
    if metric in {"pcc", "spearman", "r2", "nrmse"} and float(np.std(true)) == 0.0:
        return "constant_true"
    return "undefined"


def score_groups(
    group_ids: np.ndarray,
    true: np.ndarray,
    pred: np.ndarray,
    meta: pd.DataFrame,
) -> list[dict]:
    """Score every test group at k=100/200/500 and the regression metrics."""

    records = []
    for i, group in enumerate(group_ids):
        t = np.asarray(true[i], dtype=np.float64)
        p = np.asarray(pred[i], dtype=np.float64)
        by_k = {k: group_metrics(t, p, k) for k in K_LIST}
        row = meta.loc[str(group)]
        records.append(
            {
                "group_id": str(group),
                "drug_key": str(row["drug_key"]),
                "dose_um": float(row["dose_um"]),
                "atc_moa": bool(row["pharm_atc_moa"]),
                "metrics": by_k,
                "regression": by_k[200],
            }
        )
    return records


def aggregate(records: list[dict]) -> dict:
    """Macro aggregate over a set of group records."""

    out: dict = {"n_groups": len(records), "k": {}, "regression": {}}
    for k in K_LIST:
        out["k"][k] = {}
        for direction, key in (("macro", "auprc"), ("up", "auprc_up"), ("down", "auprc_down")):
            vals = [r["metrics"][k][key] for r in records]
            out["k"][k][direction] = {
                "mean": _mean(vals),
                "sd": _sd(vals),
                "min": _min(vals),
                "nan": _nan_count(vals),
            }
    for metric in ("pcc", "spearman", "r2", "rmse", "nrmse"):
        vals = [r["regression"][metric] for r in records]
        out["regression"][metric] = {
            "mean": _mean(vals),
            "sd": _sd(vals),
            "min": _min(vals),
            "nan": _nan_count(vals),
        }
    return out


def _metric_row(
    model: str,
    fold: int,
    level: str,
    subset: str,
    group_id: str,
    drug_key: str,
    dose_um: float,
    k: int,
    direction: str,
    metric: str,
    value: float,
    nan_reason: str,
) -> dict:
    return {
        "model": model,
        "fold": fold,
        "level": level,
        "subset": subset,
        "group_id": group_id,
        "drug_key": drug_key,
        "dose_um": dose_um,
        "k": int(k),
        "direction": direction,
        "metric": metric,
        "value": float(value),
        "nan_reason": nan_reason,
    }


def _rows_from_aggregate(
    model: str,
    fold: int,
    level: str,
    subset: str,
    agg: dict,
    group_id: str = "",
    drug_key: str = "",
    dose_um: float = math.nan,
) -> list[dict]:
    rows = []
    no_groups = agg["n_groups"] == 0
    for k in K_LIST:
        for direction in ("macro", "up", "down"):
            stats = agg["k"][k][direction]
            reason = "no_groups" if no_groups else ("all_nan" if math.isnan(stats["mean"]) else "")
            sd_reason = (
                "no_groups"
                if no_groups
                else ("single_group" if math.isnan(stats["sd"]) and agg["n_groups"] < 2 else ("all_nan" if math.isnan(stats["sd"]) else ""))
            )
            rows.append(
                _metric_row(
                    model, fold, level, subset, group_id, drug_key, dose_um,
                    k, direction, "auprc", stats["mean"], reason,
                )
            )
            rows.append(
                _metric_row(
                    model, fold, level, subset, group_id, drug_key, dose_um,
                    k, direction, "auprc_sd", stats["sd"], sd_reason,
                )
            )
    for metric in ("pcc", "spearman", "r2", "rmse", "nrmse"):
        stats = agg["regression"][metric]
        reason = "no_groups" if no_groups else ("all_nan" if math.isnan(stats["mean"]) else "")
        sd_reason = (
            "no_groups"
            if no_groups
            else ("single_group" if math.isnan(stats["sd"]) and agg["n_groups"] < 2 else ("all_nan" if math.isnan(stats["sd"]) else ""))
        )
        rows.append(
            _metric_row(
                model, fold, level, subset, group_id, drug_key, dose_um,
                0, "macro", metric, stats["mean"], reason,
            )
        )
        rows.append(
            _metric_row(
                model, fold, level, subset, group_id, drug_key, dose_um,
                0, "macro", f"{metric}_sd", stats["sd"], sd_reason,
            )
        )
    return rows


def _headline_rows(
    model: str, fold_aggs: dict, subset: str
) -> tuple[list[dict], dict]:
    """Aggregate fold-level macro metrics across folds."""

    rows = []
    summary = {"k": {}, "regression": {}}
    for k in K_LIST:
        summary["k"][k] = {}
        for direction in ("macro", "up", "down"):
            vals = [
                fold_aggs[f]["k"][k][direction]["mean"]
                for f in sorted(fold_aggs)
            ]
            stats = {"mean": _mean(vals), "sd": _sd(vals), "min": _min(vals)}
            summary["k"][k][direction] = stats
            rows.append(
                _metric_row(
                    model, -1, "headline", subset, "", "", math.nan,
                    k, direction, "auprc", stats["mean"],
                    "all_nan" if math.isnan(stats["mean"]) else "",
                )
            )
            rows.append(
                _metric_row(
                    model, -1, "headline", subset, "", "", math.nan,
                    k, direction, "auprc_sd", stats["sd"],
                    "all_nan" if math.isnan(stats["sd"]) else "",
                )
            )
            rows.append(
                _metric_row(
                    model, -1, "headline", subset, "", "", math.nan,
                    k, direction, "auprc_min", stats["min"],
                    "all_nan" if math.isnan(stats["min"]) else "",
                )
            )
    for metric in ("pcc", "spearman", "r2", "rmse", "nrmse"):
        vals = [fold_aggs[f]["regression"][metric]["mean"] for f in sorted(fold_aggs)]
        stats = {"mean": _mean(vals), "sd": _sd(vals), "min": _min(vals)}
        summary["regression"][metric] = stats
        rows.append(
            _metric_row(
                model, -1, "headline", subset, "", "", math.nan,
                0, "macro", metric, stats["mean"],
                "all_nan" if math.isnan(stats["mean"]) else "",
            )
        )
        rows.append(
            _metric_row(
                model, -1, "headline", subset, "", "", math.nan,
                0, "macro", f"{metric}_sd", stats["sd"],
                "all_nan" if math.isnan(stats["sd"]) else "",
            )
        )
        rows.append(
            _metric_row(
                model, -1, "headline", subset, "", "", math.nan,
                0, "macro", f"{metric}_min", stats["min"],
                "all_nan" if math.isnan(stats["min"]) else "",
            )
        )
    return rows, summary


def _load_fold(
    output_dir: Path,
    model: str,
    fold: int,
    data: BenchData,
) -> dict | None:
    run_dir = output_dir / model / f"fold_{fold}"
    status_path = run_dir / "status.json"
    pred_path = run_dir / "pred_test.npy"
    group_path = run_dir / "group_ids.json"
    if not (status_path.exists() and pred_path.exists() and group_path.exists()):
        return None
    status = json.loads(status_path.read_text())
    if status.get("status") not in {"completed", "timed_out"}:
        return None
    if not pred_path.exists():
        return None
    pred = np.load(pred_path, mmap_mode="r")
    group_ids = np.asarray(json.loads(group_path.read_text())["test"]).astype(str)
    if pred.shape != (len(group_ids), 18164):
        raise ValueError(f"{model} fold {fold}: prediction shape {pred.shape}")
    unique, group_pred, counts = group_mean_predictions(np.asarray(pred), group_ids)
    missing = [g for g in unique if g not in data.group_row]
    if missing:
        raise ValueError(f"{model} fold {fold}: unknown groups {missing[:3]}")
    true = np.stack([data.group_delta[data.group_row[g]] for g in unique]).astype(np.float64)
    meta = data.group_meta.loc[list(unique)]
    records = score_groups(unique, true, group_pred, meta)
    status["n_groups"] = int(len(unique))
    status["n_wells"] = int(len(group_ids))
    return {
        "status": status,
        "group_ids": unique,
        "group_pred": group_pred.astype(np.float32),
        "true": true.astype(np.float32),
        "records": records,
        "counts": counts,
    }


def _bootstrap(
    models: list[str],
    fold_records: dict[str, dict[int, list[dict]]],
    leaderboard: pd.DataFrame,
) -> dict:
    """Stratified cluster bootstrap over test drugs, seed 20261022."""

    models_sorted = sorted(models)
    model_index = {m: i for i, m in enumerate(models_sorted)}
    available_folds = sorted({f for m in models for f in fold_records[m]})
    if not available_folds:
        return {"models": np.asarray([], dtype=str)}

    fold_drugs: dict[int, list[str]] = {}
    fold_drug_groups: dict[int, dict[str, list[str]]] = {}
    for fold in available_folds:
        any_records = next(iter(fold_records[m][fold] for m in models if fold in fold_records[m]))
        drugs = sorted({r["drug_key"] for r in any_records})
        fold_drugs[fold] = drugs
        mapping: dict[str, list[str]] = {d: [] for d in drugs}
        for r in any_records:
            mapping[r["drug_key"]].append(r["group_id"])
        for d in mapping:
            mapping[d].sort()
        fold_drug_groups[fold] = mapping

    # drug_ap[fold] shape (n_drugs, max_groups, n_models)
    max_groups = max(len(v) for mapping in fold_drug_groups.values() for v in mapping.values())
    drug_ap: dict[int, np.ndarray] = {}
    for fold in available_folds:
        drugs = fold_drugs[fold]
        arr = np.full((len(drugs), max_groups, len(models_sorted)), np.nan, dtype=np.float64)
        for m in models_sorted:
            if fold not in fold_records[m]:
                continue
            by_group = {r["group_id"]: r["metrics"][200]["auprc"] for r in fold_records[m][fold]}
            for di, drug in enumerate(drugs):
                for gi, group in enumerate(fold_drug_groups[fold][drug]):
                    if group in by_group:
                        arr[di, gi, model_index[m]] = by_group[group]
        drug_ap[fold] = arr

    rng = np.random.default_rng(BOOTSTRAP_SEED)
    n_reps = BOOTSTRAP_REPS
    drug_idx = np.full((n_reps, 5, 92), -1, dtype=np.int32)
    for fold in available_folds:
        n_drugs = len(fold_drugs[fold])
        if n_drugs != 92:
            raise ValueError(f"bootstrap expects 92 test drugs, fold {fold} has {n_drugs}")
        drug_idx[:, fold, :] = rng.integers(0, n_drugs, size=(n_reps, n_drugs))

    fold_auprc = np.full((n_reps, 5, len(models_sorted)), np.nan, dtype=np.float64)
    for fold in available_folds:
        sampled = drug_ap[fold][drug_idx[:, fold, :]]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            fold_auprc[:, fold, :] = np.nanmean(sampled, axis=(1, 2))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        headline = np.nanmean(fold_auprc, axis=1)

    ci_low = np.full(len(models_sorted), np.nan)
    ci_high = np.full(len(models_sorted), np.nan)
    for i in range(len(models_sorted)):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            lo, hi = np.nanpercentile(headline[:, i], [2.5, 97.5])
        ci_low[i] = float(lo)
        ci_high[i] = float(hi)

    # Paired adjacent comparisons in leaderboard order, same resamples.
    pair_a: list[str] = []
    pair_b: list[str] = []
    pair_low: list[float] = []
    pair_high: list[float] = []
    pair_tied: list[bool] = []
    ordered = leaderboard.sort_values("rank", kind="mergesort")["model"].tolist()
    for a, b in zip(ordered, ordered[1:]):
        ia, ib = model_index[a], model_index[b]
        diff = headline[:, ia] - headline[:, ib]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            lo, hi = np.nanpercentile(diff, [2.5, 97.5])
        pair_a.append(a)
        pair_b.append(b)
        pair_low.append(float(lo))
        pair_high.append(float(hi))
        pair_tied.append(bool(lo <= 0.0 <= hi))

    return {
        "seed": np.asarray(BOOTSTRAP_SEED, dtype=np.int64),
        "drug_idx": drug_idx,
        "fold_auprc": fold_auprc,
        "headline": headline,
        "models": np.asarray(models_sorted, dtype=str),
        "ci_low": ci_low,
        "ci_high": ci_high,
        "pair_a": np.asarray(pair_a, dtype=str),
        "pair_b": np.asarray(pair_b, dtype=str),
        "pair_low": np.asarray(pair_low, dtype=np.float64),
        "pair_high": np.asarray(pair_high, dtype=np.float64),
        "pair_tied": np.asarray(pair_tied, dtype=bool),
        "ordered": np.asarray(ordered, dtype=str),
    }


def _write_readme(path: Path) -> None:
    path.write_text(
        "# Stage 2 benchmark results\n\n"
        "Produced by `stage2/03_bench_models/evaluate.py` (called by `run.sh`).\n\n"
        "## Files\n\n"
        "| file | content |\n|---|---|\n"
        "| `metrics.csv` | long metrics: group/drug/fold/headline, all and atc-moa subsets; `k=0` rows are regression metrics |\n"
        "| `per_drug.csv` | per model/fold/drug mean over doses; AUPRC is k=200 macro |\n"
        "| `leaderboard.csv` | one row per model; rank by headline k=200 macro AUPRC |\n"
        "| `bootstrap.npz` | stratified 92-drug bootstrap, seed 20261022, 2,000 replicates |\n"
        "| `run_cards.csv` | per model/fold wall-clock, peak GPU, final/best epoch |\n"
        "| `fidelity.csv` | one row per declared deviation from the per-model ledgers |\n"
        "| `../REPORT.md` | the executed benchmark report and its limitations |\n\n"
        "## Conventions\n\n"
        "- Headline = mean of the available fold-level macro AUPRCs at k=200.  A\n"
        "  model needs all five folds to be `winner_eligible`; a dry run with one\n"
        "  fold is reported but not eligible.\n"
        "- AUPRC is always direction-macro and group-macro; there is no micro\n"
        "  AUPRC pooling.  `auprc_up`/`auprc_down` are audit rows only.\n"
        "- Regression rows use `k=0` and `direction=macro`; R² uses the true mean\n"
        "  and RMSE is in raw `log1p_cp10k` Δ units.\n"
        "- A constant prediction has NaN PCC/Spearman, counted in `nan_counts` and\n"
        "  labelled `constant_prediction`; it is not imputed as 0.\n"
        "- `ideal_fit` is `full`, `partial`, or `none` (baselines); it is a report\n"
        "  label, not a second leaderboard.\n"
        "- `subset=atc-moa` is descriptive; it never creates a second leaderboard.\n",
        encoding="utf-8",
    )


def evaluate_run(
    output_dir: Path,
    results_dir: Path,
    models: list[str] | None = None,
    folds: list[int] | None = None,
    registry: dict | None = None,
    data_dir: Path | None = None,
) -> dict:
    output_dir = Path(output_dir)
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    data_dir = Path(data_dir) if data_dir else output_dir.parent / "data"
    data = BenchData(data_dir)
    registry = registry or {}

    if models is None:
        models = sorted(
            p.name for p in output_dir.iterdir() if p.is_dir() and p.name != "_cases"
        )
    models = list(models)
    if folds is None:
        folds = list(range(5))

    all_rows: list[dict] = []
    per_drug_rows: list[dict] = []
    fold_records: dict[str, dict[int, list[dict]]] = {}
    fold_aggs: dict[str, dict[int, dict]] = {}
    fold_aggs_atc: dict[str, dict[int, dict]] = {}
    fold_status: dict[str, dict[int, dict]] = {}

    for model in models:
        fold_records[model] = {}
        fold_aggs[model] = {}
        fold_aggs_atc[model] = {}
        fold_status[model] = {}
        for fold in folds:
            loaded = _load_fold(output_dir, model, fold, data)
            if loaded is None:
                continue
            status = loaded["status"]
            fold_status[model][fold] = status
            records = loaded["records"]
            fold_records[model][fold] = records

            agg = aggregate(records)
            agg_atc = aggregate([r for r in records if r["atc_moa"]])
            fold_aggs[model][fold] = agg
            fold_aggs_atc[model][fold] = agg_atc
            all_rows.extend(_rows_from_aggregate(model, fold, "fold", "all", agg))
            all_rows.extend(_rows_from_aggregate(model, fold, "fold", "atc-moa", agg_atc))

            # Group-level rows: every group appears under `subset=all`; atc-moa
            # groups also get an atc-moa view row.
            for i, r in enumerate(records):
                subsets = ["all"] + (["atc-moa"] if r["atc_moa"] else [])
                for subset in subsets:
                    for k in K_LIST:
                        m = r["metrics"][k]
                        for direction, key in (("macro", "auprc"), ("up", "auprc_up"), ("down", "auprc_down")):
                            val = m[key]
                            reason = _nan_reason(
                                val, loaded["group_pred"][i], loaded["true"][i], "auprc"
                            )
                            all_rows.append(
                                _metric_row(
                                    model, fold, "group", subset, r["group_id"], r["drug_key"],
                                    r["dose_um"], k, direction, "auprc", val, reason,
                                )
                            )
                    for metric in ("pcc", "spearman", "r2", "rmse", "nrmse"):
                        val = r["regression"][metric]
                        reason = _nan_reason(
                            val,
                            loaded["group_pred"][i],
                            loaded["true"][i],
                            metric,
                        )
                        all_rows.append(
                            _metric_row(
                                model, fold, "group", subset, r["group_id"], r["drug_key"],
                                r["dose_um"], 0, "macro", metric, val, reason,
                            )
                        )

            # Per-drug rows and drug-level metrics.
            for subset, subset_records in (("all", records), ("atc-moa", [r for r in records if r["atc_moa"]])):
                by_drug: dict[str, list[dict]] = {}
                for r in subset_records:
                    by_drug.setdefault(r["drug_key"], []).append(r)
                for drug, drug_records in sorted(by_drug.items()):
                    agg_drug = aggregate(drug_records)
                    all_rows.extend(
                        _rows_from_aggregate(
                            model, fold, "drug", subset, agg_drug, drug_key=drug
                        )
                    )
                    if subset == "all":
                        per_drug_rows.append(
                            {
                                "model": model,
                                "fold": fold,
                                "drug_key": drug,
                                "n_groups": len(drug_records),
                                "auprc": agg_drug["k"][200]["macro"]["mean"],
                                "pcc": agg_drug["regression"]["pcc"]["mean"],
                                "spearman": agg_drug["regression"]["spearman"]["mean"],
                                "r2": agg_drug["regression"]["r2"]["mean"],
                                "rmse": agg_drug["regression"]["rmse"]["mean"],
                            }
                        )

        if not fold_aggs[model]:
            continue

        # Headline rows per subset.
        for subset, by_fold in (("all", fold_aggs[model]), ("atc-moa", fold_aggs_atc[model])):
            rows, _ = _headline_rows(model, by_fold, subset)
            all_rows.extend(rows)

    # Leaderboard.
    leaderboard_rows = []
    for model in models:
        by_fold = fold_aggs.get(model, {})
        if not by_fold:
            continue
        folds_present = sorted(by_fold)
        fold_values = {
            f: by_fold[f]["k"][200]["macro"]["mean"] for f in folds_present
        }
        headline = _mean(fold_values.values())
        headline_sd = _sd(fold_values.values())
        mean_k100 = _mean(by_fold[f]["k"][100]["macro"]["mean"] for f in folds_present)
        mean_k500 = _mean(by_fold[f]["k"][500]["macro"]["mean"] for f in folds_present)
        pcc = _mean(by_fold[f]["regression"]["pcc"]["mean"] for f in folds_present)
        spearman = _mean(by_fold[f]["regression"]["spearman"]["mean"] for f in folds_present)
        r2 = _mean(by_fold[f]["regression"]["r2"]["mean"] for f in folds_present)
        rmse = _mean(by_fold[f]["regression"]["rmse"]["mean"] for f in folds_present)
        nrmse = _mean(by_fold[f]["regression"]["nrmse"]["mean"] for f in folds_present)
        # Pooled MSE is the mean of the per-group MSEs; RMSE is its square root.
        sq = [
            r["regression"]["rmse"] ** 2
            for f in folds_present
            for r in fold_records[model][f]
            if not math.isnan(r["regression"]["rmse"])
        ]
        pooled_rmse = float(np.sqrt(np.mean(sq))) if sq else math.nan
        nan_counts = 0
        for f in folds_present:
            for r in fold_records[model][f]:
                for k in K_LIST:
                    nan_counts += _nan_count([r["metrics"][k]["auprc"], r["metrics"][k]["auprc_up"], r["metrics"][k]["auprc_down"]])
                nan_counts += _nan_count([r["regression"][m] for m in ("pcc", "spearman", "r2", "rmse", "nrmse")])
        status = "completed" if len(folds_present) == 5 and all(
            fold_status[model].get(f, {}).get("status") == "completed" for f in folds_present
        ) else "partial"
        winner_eligible = bool(
            len(folds_present) == 5
            and all(fold_status[model].get(f, {}).get("status") == "completed" for f in folds_present)
        )
        leaderboard_rows.append(
            {
                "model": model,
                "headline_auprc": headline,
                "headline_sd": headline_sd,
                "ci_low": math.nan,
                "ci_high": math.nan,
                "fold0": fold_values.get(0, math.nan),
                "fold1": fold_values.get(1, math.nan),
                "fold2": fold_values.get(2, math.nan),
                "fold3": fold_values.get(3, math.nan),
                "fold4": fold_values.get(4, math.nan),
                "mean_k100": mean_k100,
                "mean_k500": mean_k500,
                "pcc": pcc,
                "spearman": spearman,
                "r2": r2,
                "rmse": rmse,
                "pooled_rmse": pooled_rmse,
                "nrmse": nrmse,
                "nan_counts": int(nan_counts),
                "statistically_tied": False,
                "ideal_fit": _ideal_fit(registry.get(model, {}).get("ideal_fit", False)),
                "n_folds": len(folds_present),
                "winner_eligible": winner_eligible,
                "status": status,
            }
        )

    leaderboard = pd.DataFrame(leaderboard_rows, columns=LEADERBOARD_COLUMNS[1:])
    if len(leaderboard):
        leaderboard = leaderboard.sort_values(
            ["headline_auprc", "mean_k100", "mean_k500", "pcc", "model"],
            ascending=[False, False, False, False, True],
            na_position="last",
            kind="mergesort",
        ).reset_index(drop=True)
        leaderboard.insert(0, "rank", np.arange(1, len(leaderboard) + 1))
        # Bootstrap (only for models with at least one fold).
        boot = _bootstrap([m for m in models if m in fold_aggs and fold_aggs[m]], fold_records, leaderboard)
        if "models" in boot and len(boot["models"]):
            boot_models = list(map(str, boot["models"]))
            ci_low = {m: float(boot["ci_low"][i]) for i, m in enumerate(boot_models)}
            ci_high = {m: float(boot["ci_high"][i]) for i, m in enumerate(boot_models)}
            leaderboard["ci_low"] = leaderboard["model"].map(ci_low)
            leaderboard["ci_high"] = leaderboard["model"].map(ci_high)
            ties: set[str] = set()
            for a, b, tied in zip(boot["pair_a"], boot["pair_b"], boot["pair_tied"]):
                if bool(tied):
                    ties.add(str(a))
                    ties.add(str(b))
            leaderboard["statistically_tied"] = leaderboard["model"].isin(ties)
        save_npz(results_dir / "bootstrap.npz", {k: v for k, v in boot.items() if k != "ordered"})
    else:
        # No completed model: write empty but schema-valid artifacts.
        leaderboard = pd.DataFrame(columns=LEADERBOARD_COLUMNS)
        leaderboard.to_csv(results_dir / "leaderboard.csv", index=False, lineterminator="\n")

    metrics = pd.DataFrame(all_rows, columns=METRIC_COLUMNS)
    if len(metrics):
        metrics = metrics.sort_values(
            ["model", "fold", "level", "subset", "group_id", "k", "direction", "metric"],
            kind="mergesort",
            na_position="last",
        ).reset_index(drop=True)
    metrics.to_csv(
        results_dir / "metrics.csv",
        index=False,
        float_format="%.17g",
        lineterminator="\n",
    )
    per_drug = pd.DataFrame(per_drug_rows, columns=PER_DRUG_COLUMNS)
    if len(per_drug):
        per_drug = per_drug.sort_values(
            ["model", "fold", "drug_key"], kind="mergesort"
        ).reset_index(drop=True)
    per_drug.to_csv(
        results_dir / "per_drug.csv",
        index=False,
        float_format="%.17g",
        lineterminator="\n",
    )
    if len(leaderboard):
        leaderboard = leaderboard[LEADERBOARD_COLUMNS]
        leaderboard.to_csv(
            results_dir / "leaderboard.csv",
            index=False,
            float_format="%.17g",
            lineterminator="\n",
        )
    _write_readme(results_dir / "README.md")
    return {
        "leaderboard": leaderboard,
        "metrics_rows": len(metrics),
        "per_drug_rows": len(per_drug),
        "models": models,
        "folds": folds,
    }


def summarize_predictions(
    true: np.ndarray, pred: np.ndarray, group_ids: np.ndarray
) -> dict:
    """Per-fold summary used by the runner's `metrics.json`."""

    unique, group_pred, counts = group_mean_predictions(np.asarray(pred), group_ids)
    true = np.asarray(true, dtype=np.float64)
    records = []
    for i, group in enumerate(unique):
        by_k = {k: group_metrics(true[i], group_pred[i], k) for k in K_LIST}
        records.append({"metrics": by_k, "regression": by_k[200]})
    agg = aggregate(records)
    return {
        "n_groups": int(len(unique)),
        "n_wells": int(len(group_ids)),
        "auprc": {str(k): agg["k"][k]["macro"]["mean"] for k in K_LIST},
        "auprc_up": {str(k): agg["k"][k]["up"]["mean"] for k in K_LIST},
        "auprc_down": {str(k): agg["k"][k]["down"]["mean"] for k in K_LIST},
        "pcc": agg["regression"]["pcc"]["mean"],
        "spearman": agg["regression"]["spearman"]["mean"],
        "r2": agg["regression"]["r2"]["mean"],
        "rmse": agg["regression"]["rmse"]["mean"],
        "nrmse": agg["regression"]["nrmse"]["mean"],
        "nan_counts": int(
            sum(agg["k"][k][d]["nan"] for k in K_LIST for d in ("macro", "up", "down"))
            + sum(agg["regression"][m]["nan"] for m in ("pcc", "spearman", "r2", "rmse", "nrmse"))
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--model", action="append", default=None)
    parser.add_argument("--fold", action="append", type=int, default=None)
    parser.add_argument("--data-dir", type=Path, default=None)
    parser.add_argument("--registry", type=Path, default=None)
    args = parser.parse_args(argv)
    registry = json.loads(args.registry.read_text()) if args.registry else {}
    result = evaluate_run(
        args.output_dir,
        args.results_dir,
        models=args.model,
        folds=sorted(set(args.fold)) if args.fold else None,
        registry=registry,
        data_dir=args.data_dir,
    )
    print(f"wrote {result['metrics_rows']} metric rows, {result['per_drug_rows']} per-drug rows")
    if len(result["leaderboard"]):
        print(result["leaderboard"][["rank", "model", "headline_auprc", "n_folds"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
