#!/usr/bin/env python3
"""Reduce the sealed `tsg_view` run matrix to deterministic report tables.

The report layer is a *read-only* projection: it consumes the completed run
directories written by :mod:`learner` (through :mod:`matrix`), the frozen truth
product (``output/truth``), its diagnostics, and the accepted Stage 04
``Mean_cell`` rows.  It never rebuilds a delta, never reads a checkpoint, and never
scores a model.  Every table is sorted explicitly and written with
``lineterminator="\\n"``; values that feed a plot are parsed with
``float_precision="round_trip"`` so re-rendering cannot move the last bits.

Estimands (plan C5):

* a run's score is sample-equal over its profiles (``protocol.summary``);
* the headline is the **equal-fold mean**: per seed, the mean over the three folds;
  then the mean over seeds.  Fold sample SD (ddof=1) and seed sample SD (ddof=1) are
  separate columns, never pooled into one error bar;
* ``sample_equal_pooled`` weights every profile equally across the whole view and is
  recorded for comparison only -- it is never the headline;
* ``solo-deepest`` is a one-seed sensitivity fit: it is tagged ``arm_type`` and its
  paired comparison is written to its own file, never mixed with the main arms.

The per-dose table is keyed to the fold's held dose.  Each fold holds exactly one
dose, so dose and fold are confounded by design; that is stated in every row.

Mean_cell provenance: the accepted Stage 04 cold-drug-dose ``Mean_cell`` rows come
from the passed round-1 baselines report
(``stage2/04_tsg_dose/execution/historical/report_axis_local/baselines/data/round1_runs.csv``,
digest pinned by ``log/round1-baselines-acceptance.json``).  Its test rows are the
same frozen truth rows used here (sample counts are asserted equal), so the baseline
is reused rather than recomputed.  The frozen truth product publishes only
``delta.npy`` -- not the treated expression and vehicle rows a recomputation would
need -- so a missing/changed provenance file is a hard error, never a silent
substitute.

CLI::

    report.py --runs <output/runs> --truth <output/truth>
              --report <output/report> --figure-data <figure/data>
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import statistics
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import matrix  # noqa: E402
from views import SOURCES, sha256, write_json  # noqa: E402

SCHEMA = "tsg-view-report-v1"
METRICS = ("pcc", "spearman", "rmse", "r2")
METRIC_COLUMNS = tuple(f"delta_{metric}" for metric in METRICS)
SOURCE_ORDER = ("tenv", "tahoe_bulk")
VIEW_ORDER = matrix.VIEW_ORDER
ALL_VIEWS = ("full", "merge", "solo", "full-pooled-label", "solo-deepest")
FOLDS = (0, 1, 2)
MAIN_SEEDS = (42, 43)
ARM_MAIN = "main"
ARM_SOLO_DEEPEST = "solo-deepest sensitivity"
ERROR_LEVEL = 1.96

RUN_ESTIMAND = "per-profile sample-equal Delta metrics of one run (equal weight per profile)"
FOLD_ESTIMAND = ("equal-fold mean and fold sample SD (ddof=1); sample_equal_pooled is "
                 "recorded for comparison and is never the headline")
SEED_ESTIMAND = ("headline = mean of per-seed equal-fold means; seed sample SD (ddof=1); "
                 "fold and seed spread are reported separately, never pooled")
CONDITION_CAVEAT = ("the replicated-condition population is not a random subset of the source: "
                    "repeats concentrate on multi-plate, multi-article conditions")
DOSE_DESIGN_NOTE = ("dose and fold are confounded by design: each fold holds exactly one dose, "
                    "so dose and fold effects cannot be disentangled")
PAIR_ESTIMAND = ("equal-fold mean of paired (view - full) differences per seed, then equal-seed "
                 "mean; half-width = 1.96*sqrt((fold_sd/sqrt(n_folds))**2 + "
                 "(seed_sd/sqrt(n_seeds))**2), fold_sd = equal-seed mean of the per-seed fold "
                 "sample SDs (ddof=1), seed_sd = sample SD of the per-seed equal-fold means "
                 "(ddof=1); fold and seed spread are never pooled into one bar")
SD_ESTIMAND = ("fold_sd = equal-seed mean of the per-seed fold sample SDs (ddof=1); "
               "seed_sd = sample SD of the per-seed equal-fold means (ddof=1); drawn as two "
               "separate error bars, never pooled")
PAIR_ESTIMAND_SINGLE_SEED = ("equal-fold mean of paired (view - full) differences; half-width = "
                             "1.96*fold_sd/sqrt(n_folds) with fold_sd the fold sample SD (ddof=1); "
                             "single-seed sensitivity, no seed spread, reported separately from "
                             "the main arms")

MEAN_CELL_PATH = ("stage2/04_tsg_dose/execution/historical/report_axis_local/baselines/"
                  "data/round1_runs.csv")
MEAN_CELL_RECEIPT = "stage2/04_tsg_dose/log/round1-baselines-acceptance.json"
MEAN_CELL_SHA256 = "52d209006d0885f92fd792d0071573295abba8f834ef8430de59019c04a6baeb"
MEAN_CELL_DEFINITION = ("Step 04 axis-local Mean_cell: per cell line, the pooled mean of the "
                        "fitted treated expression plus that cell line's unique original "
                        "vehicles, scored Delta-space on the same frozen cold-drug-dose test rows")
NOISE_FLOOR_ESTIMAND = ("median over replicated conditions of sqrt(mean over genes of the "
                        "within-condition Delta variance, ddof=0) (population SD; headline)")
PCC_CEILING_ESTIMAND = ("profile-weighted mean of sqrt(rho_c) over the source's replicated "
                        "conditions; heuristic diagnostic, not a measurement-error model")
CEILING_LABEL = "heuristic diagnostic"

TABLE_COLUMNS = {
    "metrics_primary": (
        "source", "view", "arm_type", "fold", "seed", "n_samples", "n_defined_pcc",
        "n_defined_spearman", "n_defined_rmse", "n_defined_r2", "delta_pcc",
        "delta_spearman", "delta_rmse", "delta_r2", "constant_truth_count",
        "constant_prediction_count", "nonfinite_count", "truth", "estimand",
        "run_manifest_sha256", "run_metrics_sha256"),
    "fold_spread": (
        "truth", "source", "view", "arm_type", "seed", "metric", "equal_fold_mean",
        "fold_sd", "n_defined_folds", "sample_equal_pooled", "estimand"),
    "seed_spread": (
        "truth", "source", "view", "arm_type", "metric", "headline", "seed_sd", "n_seeds",
        "n_folds_total", "sample_equal_pooled", "estimand"),
    "metrics_condition_mean": (
        "source", "view", "arm_type", "fold", "seed", "n_conditions", "n_profiles",
        "n_defined_pcc", "n_defined_spearman", "n_defined_rmse", "n_defined_r2",
        "delta_pcc", "delta_spearman", "delta_rmse", "delta_r2",
        "truth_replicated_condition_count", "caveat_flag", "caveat", "truth", "estimand",
        "run_metrics_sha256"),
    "metrics_per_dose": (
        "source", "view", "arm_type", "fold", "seed", "held_dose_um", "held_dose_label",
        "dose_fold_confounded", "metric", "value", "n_samples", "n_defined",
        "n_profiles_total", "estimand", "design_note"),
    "paired_vs_full": (
        "source", "view", "arm_type", "metric", "delta_mean", "low", "high", "half_width",
        "fold_sd", "seed_sd", "n_folds", "n_seeds", "n_pairs", "missing_pairs",
        "interval_estimator", "status"),
    "baseline_mean_cell": (
        "source", "view", "view_scope", "fold", "metric", "value", "n_samples", "baseline",
        "scenario", "seed", "estimand", "provenance_path", "provenance_sha256",
        "acceptance_receipt", "acceptance_receipt_sha256", "truth_rows_match"),
    "diagnostics_beside_scores": (
        "source", "view", "arm_type", "n_seeds", "n_folds", "delta_pcc", "delta_spearman",
        "delta_rmse", "delta_r2", "noise_floor", "noise_floor_mean", "noise_floor_ddof",
        "noise_floor_estimand", "pcc_ceiling", "pcc_ceiling_label", "pcc_ceiling_estimand",
        "pcc_ceiling_input_conditions", "pcc_ceiling_input_profiles",
        "rmse_minus_noise_floor", "pcc_minus_ceiling", "truth", "status",
        "source_manifest_sha256", "diagnostics_sha256"),
    "run_status": (
        "source", "view", "arm_type", "fold", "seed", "status", "detail", "run_relative"),
    "learning_curves": (
        "source", "view", "arm_type", "fold", "seed", "epoch", "train_loss",
        "validation_mse", "best_validation"),
}

FIGURE_COLUMNS = {
    "learning_curves": TABLE_COLUMNS["learning_curves"],
    "learning_curves_mean": (
        "source", "view", "arm_type", "epoch", "mean_train_loss", "mean_validation_mse",
        "n_runs"),
    "metrics_bars": (
        "source", "view", "arm_type", "metric", "headline", "fold_sd", "seed_sd", "n_folds",
        "n_seeds", "sd_estimand", "mean_cell", "mean_cell_fold_sd"),
    "metrics_points": (
        "source", "view", "arm_type", "fold", "seed", "metric", "value", "x_offset"),
    "metrics_baseline": (
        "source", "metric", "mean_cell", "mean_cell_fold_sd", "n_folds", "provenance_sha256"),
    "metrics_per_dose_bars": (
        "source", "fold", "held_dose_um", "held_dose_label", "view", "arm_type", "metric",
        "value", "seed_sd", "n_seeds", "n_samples", "dose_fold_confounded", "design_note"),
    "metrics_per_dose_points": (
        "source", "fold", "held_dose_um", "view", "arm_type", "metric", "seed", "value",
        "x_offset"),
    "views_effect": (
        "source", "view", "arm_type", "metric", "delta_mean", "low", "high", "fold_sd",
        "seed_sd", "n_pairs", "status", "interval_estimator"),
}

TABLE_COLUMNS["paired_vs_full_solo_deepest"] = TABLE_COLUMNS["paired_vs_full"]

FOLD_OFFSETS = {0: -0.16, 1: 0.0, 2: 0.16}
SEED_OFFSETS = {42: -0.05, 43: 0.05}


class ReportError(ValueError):
    """A report cannot be assembled without hiding a scientific fact."""


# --------------------------------------------------------------------------- #
# small deterministic helpers
# --------------------------------------------------------------------------- #
def read_json(path):
    return json.loads(Path(path).read_text())


def read_csv(path):
    return pd.read_csv(Path(path), float_precision="round_trip")


def write_csv(frame, path, columns):
    """Write a table with a fixed header, LF endings and an explicit row order."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = frame.reindex(columns=list(columns))
    frame.to_csv(path, index=False, lineterminator="\n")
    return path


def arm_type(view):
    return ARM_SOLO_DEEPEST if view == "solo-deepest" else ARM_MAIN


def metric_label(metric):
    return f"delta_{metric}"


def _mean(values):
    values = [float(value) for value in values if value is not None
              and isinstance(value, (int, float)) and math.isfinite(float(value))]
    return float(statistics.fmean(values)) if values else None


def _stdev(values):
    values = [float(value) for value in values if value is not None
              and isinstance(value, (int, float)) and math.isfinite(float(value))]
    return float(statistics.stdev(values)) if len(values) > 1 else None


def _weighted(values, weights):
    pairs = [(float(value), float(weight)) for value, weight in zip(values, weights)
             if value is not None and weight]
    if not pairs:
        return None
    total = sum(weight for _, weight in pairs)
    return float(sum(value * weight for value, weight in pairs) / total) if total else None


def _is_finite(value):
    return isinstance(value, (int, float)) and math.isfinite(float(value))


def _sorted(frame, *, by_source=True, by_view=True, by_metric=False, by_fold=False,
            by_seed=False, by_epoch=False, truth_order=False):
    """Sort on categorical columns in declared, not lexical, order."""
    frame = frame.copy()
    keys = []
    if truth_order:
        frame["_truth"] = pd.Categorical(frame["truth"], ("primary", "condition_mean"),
                                         ordered=True)
        keys.append("_truth")
    if by_source:
        frame["_source"] = pd.Categorical(frame["source"], SOURCE_ORDER, ordered=True)
        keys.append("_source")
    if by_view:
        frame["_view"] = pd.Categorical(frame["view"], ALL_VIEWS, ordered=True)
        keys.append("_view")
    if by_metric:
        frame["_metric"] = pd.Categorical(frame["metric"], METRIC_COLUMNS, ordered=True)
        keys.append("_metric")
    if by_fold:
        keys.append("fold")
    if by_seed:
        keys.append("seed")
    if by_epoch:
        keys.append("epoch")
    frame = frame.sort_values(keys, kind="mergesort")
    return frame.drop(columns=[key for key in keys if key.startswith("_")])


def _summary(samples):
    """The sample-equal run summary, recomputed from the per-profile table."""
    result = {"n_samples": int(len(samples))}
    for metric in METRICS:
        values = pd.to_numeric(samples[metric], errors="coerce").to_numpy(dtype=float)
        finite = values[np.isfinite(values)]
        result[metric] = float(finite.mean()) if len(finite) else None
        result[f"n_defined_{metric}"] = int(len(finite))
    result["constant_truth_count"] = int(samples["constant_truth"].sum())
    result["constant_prediction_count"] = int(samples["constant_prediction"].sum())
    result["nonfinite_count"] = int(samples["nonfinite"].sum())
    return result


# --------------------------------------------------------------------------- #
# run inventory and loading
# --------------------------------------------------------------------------- #
def plan_origins():
    return matrix.plan()


def discover_run_dirs(runs_root):
    """Map ``(source, view, fold, seed)`` to the on-disk run directory."""
    runs_root = Path(runs_root)
    found = {}
    if not runs_root.is_dir():
        return found
    for run in sorted(runs_root.glob("*/*/fold-*/seed-*")):
        if not run.is_dir():
            continue
        try:
            seed = int(run.name.split("-", 1)[1])
            fold = int(run.parent.name.split("-", 1)[1])
        except (IndexError, ValueError) as error:
            raise ReportError(f"malformed run path: {run}") from error
        source = run.parent.parent.parent.name
        view = run.parent.parent.name
        found[(source, view, fold, seed)] = run
    return found


def _status(run):
    if run is None or not Path(run).is_dir():
        return "missing", "missing run directory"
    run = Path(run)
    terminal_path = run / "terminal.json"
    if not terminal_path.is_file():
        return "incomplete", "missing terminal.json"
    try:
        terminal = read_json(terminal_path)
    except json.JSONDecodeError:
        return "incomplete", "unreadable terminal.json"
    if terminal.get("status") != "complete":
        return "incomplete", f"terminal status {terminal.get('status')!r} (requires 'complete')"
    for name in ("config.json", "manifest.json", "test_metrics.json",
                 "test_per_sample.csv", "learning_curve.csv"):
        if not (run / name).is_file():
            return "incomplete", f"missing {name}"
    return "complete", ""


def _verified_outputs(run, manifest, names):
    outputs = manifest.get("outputs")
    if not isinstance(outputs, dict):
        raise ReportError(f"run manifest has no outputs mapping: {run}")
    digests = {}
    for name in names:
        expected = outputs.get(name)
        if expected is None:
            raise ReportError(f"run manifest does not register {name}: {run}")
        actual = sha256(Path(run) / name)
        if actual != expected:
            raise ReportError(f"run artifact digest mismatch ({name}): {run}")
        digests[name] = actual
    return digests


def load_complete_run(run, origin, truth_counts):
    """Read one complete run; fail closed on identity, digest or scoring drift."""
    run = Path(run)
    config = read_json(run / "config.json")
    stored_origin = config.get("origin")
    if stored_origin != origin:
        raise ReportError(f"run origin does not match its directory: {run}")
    if config.get("mode") != "formal":
        raise ReportError(f"report consumes formal runs only (mode={config.get('mode')!r}): {run}")
    terminal = read_json(run / "terminal.json")
    if int(terminal.get("fold")) != int(origin["fold"]) or int(terminal.get("seed")) != int(
            origin["seed"]):
        raise ReportError(f"terminal fold/seed do not match the origin: {run}")
    manifest = read_json(run / "manifest.json")
    if manifest.get("schema") != matrix.RUN_SCHEMA:
        raise ReportError(f"unexpected run manifest schema: {run}")
    digests = _verified_outputs(run, manifest, ("test_metrics.json", "test_per_sample.csv",
                                                "learning_curve.csv"))
    metrics = read_json(run / "test_metrics.json")
    samples = read_csv(run / "test_per_sample.csv")
    source = origin["source"]
    fold = int(origin["fold"])
    held_dose = truth_counts[source][fold]["held_dose_um"]
    if len(samples) != truth_counts[source][fold]["primary"]:
        raise ReportError(
            f"run covers {len(samples)} primary profiles but the frozen truth has "
            f"{truth_counts[source][fold]['primary']}: {run}")
    doses = sorted({float(value) for value in samples["dose_um"].to_numpy()})
    if doses != [held_dose]:
        raise ReportError(f"run test doses {doses} do not match the held dose {held_dose}: {run}")
    primary = metrics.get("primary")
    if not isinstance(primary, dict):
        raise ReportError(f"test_metrics.json carries no primary summary: {run}")
    computed = _summary(samples)
    for metric in METRICS:
        stored = primary.get(metric)
        measured = computed[metric]
        if stored is None and measured is None:
            continue
        if stored is None or measured is None or not math.isclose(
                float(stored), float(measured), rel_tol=1e-12, abs_tol=1e-12):
            raise ReportError(f"run primary {metric} summary does not reduce from its "
                              f"per-profile rows: {run}")
        if int(primary.get(f"undefined_{metric}_count", -1)) != (
                computed["n_samples"] - computed[f"n_defined_{metric}"]):
            raise ReportError(f"run primary {metric} undefined count is inconsistent: {run}")
    condition = metrics.get("condition_mean")
    if not isinstance(condition, dict):
        raise ReportError(f"test_metrics.json carries no condition_mean summary: {run}")
    n_conditions = int(metrics.get("n_condition_mean_samples", -1))
    if n_conditions != truth_counts[source][fold]["condition_mean"]:
        raise ReportError(
            f"run covers {n_conditions} replicated conditions but the frozen truth has "
            f"{truth_counts[source][fold]['condition_mean']}: {run}")
    values = {metric_label(metric): computed[metric] for metric in METRICS}
    condition_values = {metric_label(metric): condition.get(metric) for metric in METRICS}
    n_defined = {metric_label(metric): computed[f"n_defined_{metric}"] for metric in METRICS}
    n_defined_condition = {
        metric_label(metric): int(n_conditions) - int(condition.get(f"undefined_{metric}_count", 0))
        for metric in METRICS}
    return {
        "source": source, "view": origin["view"], "arm_type": arm_type(origin["view"]),
        "fold": fold, "seed": int(origin["seed"]), "n_samples": int(len(samples)),
        "n_conditions": n_conditions, "held_dose_um": held_dose,
        "held_dose_label": f"{held_dose:g} \u00b5M",
        "values": values, "condition_values": condition_values,
        "n_defined": n_defined, "n_defined_condition": n_defined_condition,
        "counts": {key: computed[key] for key in ("constant_truth_count",
                                                  "constant_prediction_count",
                                                  "nonfinite_count")},
        "run_manifest_sha256": sha256(run / "manifest.json"),
        "run_metrics_sha256": digests["test_metrics.json"],
        "learning_curve": read_csv(run / "learning_curve.csv"),
    }


# --------------------------------------------------------------------------- #
# truth, diagnostics and the accepted Mean_cell rows
# --------------------------------------------------------------------------- #
def load_truth(truth_root, source):
    base = Path(truth_root) / source
    manifest_path = base / "manifest.json"
    if not manifest_path.is_file():
        raise ReportError(f"missing frozen truth source manifest: {manifest_path}")
    manifest = read_json(manifest_path)
    if manifest.get("schema") != matrix.TRUTH_SCHEMA or manifest.get("source") != source:
        raise ReportError(f"unexpected frozen truth manifest: {manifest_path}")
    folds = {}
    for fold in FOLDS:
        primary = read_csv(base / "primary" / f"fold-{fold}" / "population.csv")
        condition = read_csv(base / "condition_mean" / f"fold-{fold}" / "population.csv")
        primary_test = primary[primary["partition"] == "test"]
        condition_test = condition[condition["partition"] == "test"]
        if not len(primary_test) or not len(condition_test):
            raise ReportError(f"empty frozen truth test partition: {source}/fold-{fold}")
        doses = sorted({float(value) for value in primary_test["dose_um"].to_numpy()})
        if len(doses) != 1:
            raise ReportError(f"fold {fold} of {source} does not hold exactly one dose: {doses}")
        members = set()
        for row in condition_test.to_dict("records"):
            ids = [member for member in str(row["members"]).split(";") if member]
            if len(ids) != int(row["count"]):
                raise ReportError(
                    f"condition-mean member count does not match its list "
                    f"({source}/fold-{fold})")
            members.update(ids)
        unknown = members - set(primary_test["profile_id"].astype(str))
        if unknown:
            raise ReportError(
                f"condition-mean truth references profiles absent from the primary truth "
                f"({source}/fold-{fold}): {sorted(unknown)[:3]}")
        if len(members) > len(primary_test):
            raise ReportError(
                f"condition-mean truth covers more profiles than the primary test rows "
                f"({source}/fold-{fold})")
        folds[fold] = {
            "primary": int(len(primary_test)),
            "condition_mean": int(len(condition_test)),
            "held_dose_um": doses[0],
        }
    return {"manifest_sha256": sha256(manifest_path), "folds": folds, "manifest": manifest}


def load_diagnostics(truth_root, source, truth_manifest_sha256):
    path = Path(truth_root) / source / "diagnostics.json"
    if not path.is_file():
        raise ReportError(f"missing frozen truth diagnostics: {path}")
    value = read_json(path)
    if value.get("schema") != "tsg-view-truth-diagnostics-v1" or value.get("source") != source:
        raise ReportError(f"unexpected diagnostics product: {path}")
    noise = value.get("noise_floor") or {}
    ceiling = value.get("pcc_ceiling") or {}
    if noise.get("median") is None or ceiling.get("value") is None:
        raise ReportError(f"diagnostics product misses the headline values: {path}")
    if ceiling.get("label") != CEILING_LABEL:
        raise ReportError(f"single-sample ceiling is not labelled {CEILING_LABEL!r}: {path}")
    return {
        "source": source,
        "path": path,
        "sha256": sha256(path),
        "truth_manifest_sha256": truth_manifest_sha256,
        "noise_floor": float(noise["median"]),
        "noise_floor_mean": float(noise["mean"]) if noise.get("mean") is not None else None,
        "noise_floor_ddof": 0,
        "pcc_ceiling": float(ceiling["value"]),
        "pcc_ceiling_label": str(ceiling["label"]),
        "pcc_ceiling_input_conditions": int(ceiling.get("input_conditions", -1)),
        "pcc_ceiling_input_profiles": int(ceiling.get("input_profiles", -1)),
    }


def initialize_mean_cell(repo_root, source_path=None, source_sha256=None, receipt_path=None):
    """Resolve and verify the accepted Stage 04 cold-drug-dose Mean_cell rows."""
    repo_root = Path(repo_root)
    path = Path(source_path) if source_path is not None else repo_root / MEAN_CELL_PATH
    receipt = Path(receipt_path) if receipt_path is not None else repo_root / MEAN_CELL_RECEIPT
    expected = source_sha256 if source_sha256 is not None else MEAN_CELL_SHA256
    if not path.is_file():
        raise ReportError(
            f"the accepted Step 04 cold-drug-dose Mean_cell rows are missing: {path}; the "
            "frozen truth product publishes Delta only, so the baseline cannot be "
            "recomputed from it without the treated expression and vehicle rows")
    digest = sha256(path)
    if digest != expected:
        raise ReportError(f"Mean_cell provenance digest mismatch: {path} is {digest}, "
                          f"expected {expected}")
    receipt_digest = sha256(receipt) if receipt.is_file() else None
    if source_path is None:
        if not receipt.is_file():
            raise ReportError(f"missing accepted Mean_cell receipt: {receipt}")
        accepted = (read_json(receipt).get("files") or {}).get("data/round1_runs.csv")
        if accepted != digest:
            raise ReportError(
                f"Mean_cell rows {digest} do not match the accepted baselines receipt "
                f"({accepted})")
    def recorded(item):
        if item is None:
            return ""
        try:
            return str(Path(item).resolve().relative_to(repo_root.resolve()))
        except ValueError:
            return str(item)

    return {"path": path, "recorded_path": recorded(path), "sha256": digest,
            "receipt": receipt if receipt.is_file() else None,
            "receipt_recorded": recorded(receipt) if receipt.is_file() else "",
            "receipt_sha256": receipt_digest, "definition": MEAN_CELL_DEFINITION}


def load_mean_cell(mean_cell, truth, sources):
    """Select the accepted Delta-space Mean_cell rows and pin them to our truth rows."""
    table = read_csv(mean_cell["path"])
    required = {"source", "scenario", "seed", "method", "fold", "n_samples", "pcc",
                "spearman", "rmse", "r2"}
    if not required.issubset(table.columns):
        raise ReportError("Mean_cell provenance table misses required columns")
    selected = table[(table["scenario"] == "cold-drug-dose")
                     & (table["method"] == "Mean_cell") & (table["seed"] == 42)]
    rows = []
    for source in sources:
        for fold in FOLDS:
            block = selected[(selected["source"] == source) & (selected["fold"] == fold)]
            if len(block) != 1:
                raise ReportError(f"expected one accepted Mean_cell row per fold: {source}/fold-{fold}")
            block = block.iloc[0]
            if int(block["n_samples"]) != truth[source]["folds"][fold]["primary"]:
                raise ReportError(
                    f"Mean_cell test rows ({int(block['n_samples'])}) are not the frozen truth "
                    f"rows ({truth[source]['folds'][fold]['primary']}): {source}/fold-{fold}")
            for metric in METRICS:
                rows.append({
                    "source": source, "view": "all",
                    "view_scope": "shared frozen truth; identical for every view",
                    "fold": fold, "metric": metric_label(metric),
                    "value": float(block[metric]), "n_samples": int(block["n_samples"]),
                    "baseline": "Mean_cell", "scenario": "cold-drug-dose", "seed": 42,
                    "estimand": MEAN_CELL_DEFINITION,
                    "provenance_path": mean_cell["recorded_path"],
                    "provenance_sha256": mean_cell["sha256"],
                    "acceptance_receipt": mean_cell["receipt_recorded"],
                    "acceptance_receipt_sha256": mean_cell["receipt_sha256"] or "",
                    "truth_rows_match": True,
                })
    return pd.DataFrame(rows, columns=TABLE_COLUMNS["baseline_mean_cell"])


# --------------------------------------------------------------------------- #
# tables
# --------------------------------------------------------------------------- #
def primary_table(runs):
    rows = []
    for run in runs:
        row = {"source": run["source"], "view": run["view"], "arm_type": run["arm_type"],
               "fold": run["fold"], "seed": run["seed"], "n_samples": run["n_samples"],
               "truth": "primary", "estimand": RUN_ESTIMAND,
               "run_manifest_sha256": run["run_manifest_sha256"],
               "run_metrics_sha256": run["run_metrics_sha256"]}
        row.update(run["counts"])
        for metric in METRICS:
            row[metric_label(metric)] = run["values"][metric_label(metric)]
            row[f"n_defined_{metric}"] = run["n_defined"][metric_label(metric)]
        rows.append(row)
    return pd.DataFrame(rows, columns=TABLE_COLUMNS["metrics_primary"])


def condition_mean_table(runs):
    rows = []
    for run in runs:
        row = {"source": run["source"], "view": run["view"], "arm_type": run["arm_type"],
               "fold": run["fold"], "seed": run["seed"],
               "n_conditions": run["n_conditions"], "n_profiles": run["n_samples"],
               "truth_replicated_condition_count": run["n_conditions"],
               "caveat_flag": True, "caveat": CONDITION_CAVEAT,
               "truth": "condition_mean", "estimand": RUN_ESTIMAND,
               "run_metrics_sha256": run["run_metrics_sha256"]}
        for metric in METRICS:
            row[metric_label(metric)] = run["condition_values"][metric_label(metric)]
            row[f"n_defined_{metric}"] = run["n_defined_condition"][metric_label(metric)]
        rows.append(row)
    return pd.DataFrame(rows, columns=TABLE_COLUMNS["metrics_condition_mean"])


def _metric_frame(runs, truth, value_key, defined_key):
    rows = []
    for run in runs:
        for metric in METRICS:
            rows.append({
                "truth": truth, "source": run["source"], "view": run["view"],
                "arm_type": run["arm_type"], "fold": run["fold"], "seed": run["seed"],
                "metric": metric_label(metric), "value": run[value_key][metric_label(metric)],
                "n_defined": run[defined_key][metric_label(metric)]})
    return pd.DataFrame(rows, columns=("truth", "source", "view", "arm_type", "fold", "seed",
                                       "metric", "value", "n_defined"))


def spread_tables(runs):
    """Equal-fold means with separate fold and seed spreads for both truths."""
    frames = [_metric_frame(runs, "primary", "values", "n_defined"),
              _metric_frame(runs, "condition_mean", "condition_values", "n_defined_condition")]
    frame = pd.concat(frames, ignore_index=True)
    fold_rows, seed_rows = [], []
    for key, block in frame.groupby(["truth", "source", "view", "arm_type", "seed", "metric"],
                                    sort=True):
        truth, source, view, kind, seed, metric = key
        values = [row["value"] for _, row in block.iterrows()]
        equal = _mean(values)
        fold_rows.append({
            "truth": truth, "source": source, "view": view, "arm_type": kind, "seed": int(seed),
            "metric": metric, "equal_fold_mean": equal, "fold_sd": _stdev(values),
            "n_defined_folds": int(len([value for value in values if value is not None
                                        and _is_finite(value)])),
            "sample_equal_pooled": _weighted(
                values, [row["n_defined"] for _, row in block.iterrows()]),
            "estimand": FOLD_ESTIMAND})
    for key, block in frame.groupby(["truth", "source", "view", "arm_type", "metric"],
                                    sort=True):
        truth, source, view, kind, metric = key
        block = block.sort_values(["seed", "fold"])
        seed_means, folds = [], set()
        for seed, seed_block in block.groupby("seed", sort=True):
            values = [row["value"] for _, row in seed_block.iterrows()]
            mean = _mean(values)
            if mean is not None:
                seed_means.append(mean)
            folds.update(int(fold) for fold in seed_block["fold"])
        seed_rows.append({
            "truth": truth, "source": source, "view": view, "arm_type": kind, "metric": metric,
            "headline": _mean(seed_means), "seed_sd": _stdev(seed_means),
            "n_seeds": len(seed_means), "n_folds_total": len(folds),
            "sample_equal_pooled": _weighted(
                [row["value"] for _, row in block.iterrows()],
                [row["n_defined"] for _, row in block.iterrows()]),
            "estimand": SEED_ESTIMAND})
    folds = pd.DataFrame(fold_rows, columns=TABLE_COLUMNS["fold_spread"])
    seeds = pd.DataFrame(seed_rows, columns=TABLE_COLUMNS["seed_spread"])
    return (_sorted(folds, by_metric=True, by_seed=True, truth_order=True),
            _sorted(seeds, by_metric=True, truth_order=True))


def per_dose_table(runs):
    rows = []
    for run in runs:
        for metric in METRICS:
            name = metric_label(metric)
            rows.append({
                "source": run["source"], "view": run["view"], "arm_type": run["arm_type"],
                "fold": run["fold"], "seed": run["seed"],
                "held_dose_um": run["held_dose_um"], "held_dose_label": run["held_dose_label"],
                "dose_fold_confounded": True, "metric": name,
                "value": run["values"][name], "n_samples": run["n_samples"],
                "n_defined": run["n_defined"][name], "n_profiles_total": run["n_samples"],
                "estimand": ("sample-equal Delta metric over the fold's held-dose test profiles "
                             "(each fold holds exactly one dose)"),
                "design_note": DOSE_DESIGN_NOTE})
    return _sorted(pd.DataFrame(rows, columns=TABLE_COLUMNS["metrics_per_dose"]),
                   by_metric=True, by_fold=True, by_seed=True)


def _pair_interval(pairs, *, single_seed):
    """Summarize paired differences with fold and seed spread kept separate."""
    by_seed = {}
    for pair in pairs:
        by_seed.setdefault(pair["seed"], []).append(pair)
    seed_means = [_mean([pair["delta"] for pair in block]) for _, block in sorted(by_seed.items())]
    seed_means = [value for value in seed_means if value is not None]
    fold_sds = [_stdev([pair["delta"] for pair in block]) for _, block in sorted(by_seed.items())]
    fold_sds = [value for value in fold_sds if value is not None]
    headline = _mean(seed_means)
    fold_sd = _mean(fold_sds)
    seed_sd = _stdev(seed_means) if len(seed_means) > 1 else None
    n_folds = len({pair["fold"] for pair in pairs})
    n_seeds = len(seed_means)
    if not pairs:
        return {"delta_mean": None, "low": None, "high": None, "half_width": None,
                "fold_sd": None, "seed_sd": None, "n_folds": 0, "n_seeds": 0,
                "status": "no-pairs"}
    if headline is None or fold_sd is None or n_folds < 2:
        return {"delta_mean": headline, "low": None, "high": None, "half_width": None,
                "fold_sd": fold_sd, "seed_sd": seed_sd, "n_folds": n_folds,
                "n_seeds": n_seeds, "status": "insufficient-folds"}
    if single_seed or seed_sd is None:
        half = ERROR_LEVEL * fold_sd / math.sqrt(n_folds)
        return {"delta_mean": headline, "low": headline - half, "high": headline + half,
                "half_width": half, "fold_sd": fold_sd, "seed_sd": seed_sd,
                "n_folds": n_folds, "n_seeds": n_seeds, "status": "ok"}
    half = ERROR_LEVEL * math.sqrt((fold_sd / math.sqrt(n_folds)) ** 2
                                   + (seed_sd / math.sqrt(n_seeds)) ** 2)
    return {"delta_mean": headline, "low": headline - half, "high": headline + half,
            "half_width": half, "fold_sd": fold_sd, "seed_sd": seed_sd,
            "n_folds": n_folds, "n_seeds": n_seeds, "status": "ok"}


def paired_tables(runs):
    """Paired view-vs-`full` differences; `solo-deepest` goes to its own table."""
    primary = [run for run in runs]
    index = {}
    for run in primary:
        for metric in METRIC_COLUMNS:
            index[(run["source"], run["view"], run["fold"], run["seed"], metric)] = run["values"][metric]
    main_rows, solo_rows = [], []
    for source in SOURCE_ORDER:
        views = VIEW_ORDER[source]
        for view in views:
            if view == "full":
                continue
            target = solo_rows if view == "solo-deepest" else main_rows
            for metric in METRIC_COLUMNS:
                pairs, missing = [], 0
                for fold in FOLDS:
                    for seed in matrix.seeds_for(view):
                        mine = index.get((source, view, fold, seed, metric))
                        full = index.get((source, "full", fold, seed, metric))
                        if mine is None or full is None:
                            missing += 1
                            continue
                        pairs.append({"fold": fold, "seed": seed, "delta": float(mine) - float(full)})
                summary = _pair_interval(pairs, single_seed=(view == "solo-deepest"))
                summary.update({
                    "source": source, "view": view, "arm_type": arm_type(view),
                    "metric": metric, "n_pairs": len(pairs), "missing_pairs": missing,
                    "interval_estimator": (PAIR_ESTIMAND_SINGLE_SEED if view == "solo-deepest"
                                           else PAIR_ESTIMAND)})
                target.append(summary)
    main = _sorted(pd.DataFrame(main_rows, columns=TABLE_COLUMNS["paired_vs_full"]),
                   by_metric=True)
    solo = _sorted(pd.DataFrame(solo_rows, columns=TABLE_COLUMNS["paired_vs_full"]),
                   by_metric=True)
    return main, solo


def diagnostics_table(diagnostics, seeds):
    primary_seeds = [row for row in seeds.to_dict("records") if row["truth"] == "primary"]
    lookup = {(row["source"], row["view"], row["metric"]): row for row in primary_seeds}
    rows = []
    for source in SOURCE_ORDER:
        diag = diagnostics.get(source)
        if diag is None:
            continue
        for view in VIEW_ORDER[source]:
            row = {"source": source, "view": view, "arm_type": arm_type(view),
                   "n_seeds": 0, "n_folds": 0, "truth": "primary"}
            observed = False
            for metric in METRIC_COLUMNS:
                seed_row = lookup.get((source, view, metric))
                value = seed_row["headline"] if seed_row is not None else None
                row[metric] = value
                if seed_row is not None:
                    observed = True
                    row["n_seeds"] = max(row["n_seeds"], int(seed_row["n_seeds"]))
                    row["n_folds"] = max(row["n_folds"], int(seed_row["n_folds_total"]))
            row["status"] = "ok" if observed else "pending"
            row.update({
                "noise_floor": diag["noise_floor"], "noise_floor_mean": diag["noise_floor_mean"],
                "noise_floor_ddof": diag["noise_floor_ddof"],
                "noise_floor_estimand": NOISE_FLOOR_ESTIMAND,
                "pcc_ceiling": diag["pcc_ceiling"], "pcc_ceiling_label": diag["pcc_ceiling_label"],
                "pcc_ceiling_estimand": PCC_CEILING_ESTIMAND,
                "pcc_ceiling_input_conditions": diag["pcc_ceiling_input_conditions"],
                "pcc_ceiling_input_profiles": diag["pcc_ceiling_input_profiles"],
                "rmse_minus_noise_floor": (None if row["delta_rmse"] is None
                                           else row["delta_rmse"] - diag["noise_floor"]),
                "pcc_minus_ceiling": (None if row["delta_pcc"] is None
                                      else row["delta_pcc"] - diag["pcc_ceiling"]),
                "source_manifest_sha256": diag["truth_manifest_sha256"],
                "diagnostics_sha256": diag["sha256"]})
            rows.append(row)
    frame = pd.DataFrame(rows, columns=TABLE_COLUMNS["diagnostics_beside_scores"])
    return _sorted(frame)


def status_table(origins, found):
    rows = []
    for origin in origins:
        key = matrix.origin_key(origin)
        status, detail = _status(found.get(key))
        rows.append({
            "source": origin["source"], "view": origin["view"],
            "arm_type": arm_type(origin["view"]), "fold": int(origin["fold"]),
            "seed": int(origin["seed"]), "status": status, "detail": detail,
            "run_relative": str(matrix.run_relative(origin))})
    order = {key: index for index, key in enumerate(matrix.origin_key(origin)
                                                    for origin in origins)}
    frame = pd.DataFrame(rows, columns=TABLE_COLUMNS["run_status"])
    frame["_order"] = [order[(row["source"], row["view"], row["fold"], row["seed"])]
                       for row in frame.to_dict("records")]
    return frame.sort_values("_order", kind="mergesort").drop(columns="_order")


# --------------------------------------------------------------------------- #
# learning curves and figure-data projections
# --------------------------------------------------------------------------- #
def learning_curves(runs):
    rows = []
    for run in runs:
        curve = run["learning_curve"].copy()
        required = {"epoch", "train_loss", "validation_mse"}
        if not required.issubset(curve.columns):
            raise ReportError(f"learning curve misses required columns: "
                              f"{run['source']}/{run['view']}/fold-{run['fold']}/seed-{run['seed']}")
        curve = curve.sort_values("epoch")
        epochs = [int(value) for value in curve["epoch"]]
        if epochs != list(range(1, len(epochs) + 1)):
            raise ReportError("learning curve epochs are not contiguous 1..N: "
                              f"{run['source']}/{run['view']}/fold-{run['fold']}/seed-{run['seed']}")
        losses = pd.to_numeric(curve["validation_mse"], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(losses).all():
            raise ReportError("learning curve carries a non-finite validation loss")
        best = int(np.argmin(losses))
        for position, (_, row) in enumerate(curve.iterrows()):
            rows.append({
                "source": run["source"], "view": run["view"], "arm_type": run["arm_type"],
                "fold": run["fold"], "seed": run["seed"], "epoch": int(row["epoch"]),
                "train_loss": float(row["train_loss"]),
                "validation_mse": float(row["validation_mse"]),
                "best_validation": int(position == best)})
    return _sorted(pd.DataFrame(rows, columns=TABLE_COLUMNS["learning_curves"]),
                   by_fold=True, by_seed=True)


def learning_curve_means(curves):
    rows = []
    for key, block in curves.groupby(["source", "view", "arm_type", "epoch"], sort=True):
        source, view, kind, epoch = key
        rows.append({
            "source": source, "view": view, "arm_type": kind, "epoch": int(epoch),
            "mean_train_loss": _mean(block["train_loss"]),
            "mean_validation_mse": _mean(block["validation_mse"]),
            "n_runs": int(len(block))})
    return _sorted(pd.DataFrame(rows, columns=FIGURE_COLUMNS["learning_curves_mean"]),
                   by_epoch=True)


def figure_data(runs, seeds, fold_spread, per_dose, paired, solo_paired, baseline,
                sources_present):
    """Project the report tables onto the exact values the renderer draws."""
    primary = seeds[seeds["truth"] == "primary"]
    bars, points = [], []
    for row in primary.to_dict("records"):
        bars.append({
            "source": row["source"], "view": row["view"], "arm_type": row["arm_type"],
            "metric": row["metric"], "headline": row["headline"], "fold_sd": None,
            "seed_sd": row["seed_sd"], "n_folds": row["n_folds_total"],
            "n_seeds": row["n_seeds"], "sd_estimand": SD_ESTIMAND,
            "mean_cell": None, "mean_cell_fold_sd": None})
    for bar in bars:
        fold_rows = [row for row in fold_spread.to_dict("records")
                     if row["truth"] == "primary" and row["source"] == bar["source"]
                     and row["view"] == bar["view"] and row["metric"] == bar["metric"]]
        bar["fold_sd"] = _mean([row["fold_sd"] for row in fold_rows])
        block = baseline[(baseline["source"] == bar["source"])
                         & (baseline["metric"] == bar["metric"])]
        values = [float(value) for value in block["value"]]
        bar["mean_cell"] = _mean(values)
        bar["mean_cell_fold_sd"] = _stdev(values)
    for run in runs:
        for metric in METRIC_COLUMNS:
            points.append({
                "source": run["source"], "view": run["view"], "arm_type": run["arm_type"],
                "fold": run["fold"], "seed": run["seed"], "metric": metric,
                "value": run["values"][metric],
                "x_offset": round(FOLD_OFFSETS[run["fold"]] + SEED_OFFSETS.get(run["seed"], 0.0), 6)})
    baseline_rows = []
    for key, block in baseline.groupby(["source", "metric"], sort=True):
        source, metric = key
        values = [float(value) for value in block["value"]]
        baseline_rows.append({"source": source, "metric": metric, "mean_cell": _mean(values),
                              "mean_cell_fold_sd": _stdev(values), "n_folds": int(len(values)),
                              "provenance_sha256": block["provenance_sha256"].iloc[0]})
    dose_bars, dose_points = [], []
    for key, block in per_dose.groupby(["source", "fold", "view", "arm_type", "metric"],
                                       sort=True):
        source, fold, view, kind, metric = key
        values = [row["value"] for row in block.to_dict("records")]
        first = block.iloc[0]
        dose_bars.append({
            "source": source, "fold": int(fold), "held_dose_um": float(first["held_dose_um"]),
            "held_dose_label": first["held_dose_label"], "view": view, "arm_type": kind,
            "metric": metric, "value": _mean(values), "seed_sd": _stdev(values),
            "n_seeds": int(len(block)), "n_samples": int(first["n_samples"]),
            "dose_fold_confounded": True, "design_note": DOSE_DESIGN_NOTE})
        for row in block.to_dict("records"):
            dose_points.append({
                "source": source, "fold": int(fold), "held_dose_um": float(row["held_dose_um"]),
                "view": view, "arm_type": kind, "metric": metric, "seed": int(row["seed"]),
                "value": row["value"],
                "x_offset": SEED_OFFSETS.get(int(row["seed"]), 0.0)})
    effects = []
    for table in (paired, solo_paired):
        for row in table.to_dict("records"):
            if row["source"] not in sources_present:
                continue
            effects.append({key: row[key] for key in FIGURE_COLUMNS["views_effect"]})
    return {
        "learning_curves": learning_curves(runs),
        "learning_curves_mean": learning_curve_means(learning_curves(runs)),
        "metrics_bars": pd.DataFrame(bars, columns=FIGURE_COLUMNS["metrics_bars"]),
        "metrics_points": pd.DataFrame(points, columns=FIGURE_COLUMNS["metrics_points"]),
        "metrics_baseline": pd.DataFrame(baseline_rows,
                                         columns=FIGURE_COLUMNS["metrics_baseline"]),
        "metrics_per_dose_bars": pd.DataFrame(dose_bars,
                                              columns=FIGURE_COLUMNS["metrics_per_dose_bars"]),
        "metrics_per_dose_points": pd.DataFrame(dose_points,
                                                columns=FIGURE_COLUMNS["metrics_per_dose_points"]),
        "views_effect": pd.DataFrame(effects, columns=FIGURE_COLUMNS["views_effect"]),
    }


# --------------------------------------------------------------------------- #
# orchestration
# --------------------------------------------------------------------------- #
def build_report(*, runs_root, truth_root, repo_root, mean_cell_source=None,
                 mean_cell_sha256=None, mean_cell_receipt=None):
    runs_root = Path(runs_root)
    truth_root = Path(truth_root)
    repo_root = Path(repo_root)
    origins = plan_origins()
    found = discover_run_dirs(runs_root)
    planned = {matrix.origin_key(origin) for origin in origins}
    unexpected = sorted(set(found) - planned)
    if unexpected:
        raise ReportError(f"run root carries origins outside the sealed plan: {unexpected}")
    truth = {source: load_truth(truth_root, source) for source in SOURCES}
    statuses = status_table(origins, found)
    runs = []
    for origin in origins:
        status, _ = _status(found.get(matrix.origin_key(origin)))
        if status == "complete":
            runs.append(load_complete_run(found[matrix.origin_key(origin)], origin,
                                          {source: value["folds"] for source, value in truth.items()}))
    sources = [source for source in SOURCE_ORDER
               if any(run["source"] == source for run in runs)]
    diagnostics = {source: load_diagnostics(truth_root, source,
                                            truth[source]["manifest_sha256"])
                   for source in sources}
    mean_cell = initialize_mean_cell(repo_root, mean_cell_source, mean_cell_sha256,
                                     mean_cell_receipt)
    baseline = load_mean_cell(mean_cell, truth, sources)
    primary = primary_table(runs)
    condition = condition_mean_table(runs)
    folds, seeds = spread_tables(runs)
    per_dose = per_dose_table(runs)
    paired, solo_paired = paired_tables(runs)
    diagnostics_beside = diagnostics_table(diagnostics, seeds)
    curves = learning_curves(runs)
    figures = figure_data(runs, seeds, folds, per_dose, paired, solo_paired, baseline,
                          sources)
    report = {
        "metrics_primary": primary,
        "fold_spread": folds,
        "seed_spread": seeds,
        "metrics_condition_mean": condition,
        "metrics_per_dose": per_dose,
        "paired_vs_full": paired,
        "paired_vs_full_solo_deepest": solo_paired,
        "baseline_mean_cell": baseline,
        "diagnostics_beside_scores": diagnostics_beside,
        "run_status": statuses,
        "learning_curves": curves,
    }
    return {"tables": report, "figure_data": figures, "truth": truth,
            "diagnostics": diagnostics, "mean_cell": mean_cell, "runs": runs,
            "statuses": statuses}


def write_report(result, report_dir, figure_dir):
    report_dir = Path(report_dir)
    figure_dir = Path(figure_dir)
    report_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    files = {}
    for name, frame in result["tables"].items():
        path = write_csv(frame, report_dir / f"{name}.csv", TABLE_COLUMNS[name])
        files[f"{name}.csv"] = sha256(path)
    for name, frame in result["figure_data"].items():
        path = write_csv(frame, figure_dir / f"{name}.csv", FIGURE_COLUMNS[name])
        files[f"figure/{name}.csv"] = sha256(path)
    return files


def report_manifest(result, files):
    truth = result["truth"]
    diagnostics = result["diagnostics"]
    statuses = result["statuses"]
    complete = statuses[statuses["status"] == "complete"]
    missing = statuses[statuses["status"] != "complete"]
    return {
        "schema": SCHEMA,
        "complete_origins": int(len(complete)),
        "pending_origins": int(len(missing)),
        "sources": list(result["diagnostics"]),
        "truth_manifest_sha256": {source: value["manifest_sha256"]
                                  for source, value in truth.items()},
        "diagnostics_sha256": {source: value["sha256"]
                               for source, value in diagnostics.items()},
        "mean_cell": {"path": result["mean_cell"]["recorded_path"],
                      "sha256": result["mean_cell"]["sha256"],
                      "acceptance_receipt": (result["mean_cell"]["receipt_recorded"] or None),
                      "acceptance_receipt_sha256": result["mean_cell"]["receipt_sha256"],
                      "definition": result["mean_cell"]["definition"]},
        "estimands": {"run": RUN_ESTIMAND, "fold": FOLD_ESTIMAND, "seed": SEED_ESTIMAND,
                      "condition_mean_caveat": CONDITION_CAVEAT,
                      "dose_design_note": DOSE_DESIGN_NOTE,
                      "paired": PAIR_ESTIMAND,
                      "paired_single_seed": PAIR_ESTIMAND_SINGLE_SEED,
                      "noise_floor": NOISE_FLOOR_ESTIMAND,
                      "pcc_ceiling": PCC_CEILING_ESTIMAND},
        "tables": {name: digest for name, digest in sorted(files.items())
                   if not name.startswith("figure/")},
        "figure_data": {name: digest for name, digest in sorted(files.items())
                        if name.startswith("figure/")},
    }


def repo_root():
    return Path(__file__).resolve().parents[3]


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--runs", type=Path, required=True,
                        help="run root carrying <source>/<view>/fold-<f>/seed-<s>/")
    parser.add_argument("--truth", type=Path, required=True, help="frozen truth product root")
    parser.add_argument("--report", type=Path, required=True, help="report-table output directory")
    parser.add_argument("--figure-data", type=Path, required=True,
                        help="figure/data directory for the plotted CSVs")
    parser.add_argument("--root", type=Path, default=None,
                        help="checkout that carries the accepted Stage 04 Mean_cell rows")
    parser.add_argument("--mean-cell-source", type=Path, default=None,
                        help="override the accepted Mean_cell rows (tests)")
    parser.add_argument("--mean-cell-sha256", default=None,
                        help="expected digest of --mean-cell-source")
    parser.add_argument("--mean-cell-receipt", type=Path, default=None)
    args = parser.parse_args(argv)
    root = args.root if args.root is not None else repo_root()
    result = build_report(runs_root=args.runs, truth_root=args.truth, repo_root=root,
                          mean_cell_source=args.mean_cell_source,
                          mean_cell_sha256=args.mean_cell_sha256,
                          mean_cell_receipt=args.mean_cell_receipt)
    files = write_report(result, args.report, args.figure_data)
    manifest = report_manifest(result, files)
    write_json(Path(args.report) / "report_manifest.json", manifest)
    summary = {
        "schema": SCHEMA,
        "complete_origins": manifest["complete_origins"],
        "pending_origins": manifest["pending_origins"],
        "sources": manifest["sources"],
        "tables": sorted(name for name in files if not name.startswith("figure/")),
        "figure_data": sorted(name for name in files if name.startswith("figure/")),
    }
    sys.stdout.write(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
