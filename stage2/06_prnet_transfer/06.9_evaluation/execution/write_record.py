#!/usr/bin/env python
"""Aggregate the frozen evaluations into the ticket-#9 record and report inputs.

    python execution/write_record.py --phase test          # the tracked record
    python execution/write_record.py --phase valid --allow-subset   # harness dry run
    python execution/write_record.py --phase test --check   # re-derive and byte-compare

Everything is recomputed from the evaluator's frozen per-unit output plus #6/#7/#8's
tracked records. The full test record lands in ``output/record/``; a valid/subset dry run
lands in ``output/dry_run/`` and is labelled as such, so dry-run numbers can never be
mistaken for the first-round test result.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

import aggregate as aggregate_mod  # noqa: E402
import checks as checks_mod  # noqa: E402
import common  # noqa: E402
import figures as figures_mod  # noqa: E402
from common import (  # noqa: E402
    DEFAULT_BASELINES_DIR,
    DEFAULT_DATA_DIR,
    DEFAULT_OUT,
    load_json,
    now_utc,
    peak_rss_mb,
    sha256_file,
    write_json,
    write_text_atomic,
)
from units import load_unit_manifest  # noqa: E402

SCHEMA = "stage2/06.9_evaluation/record/v1"
SUMMARY_SCHEMA = "stage2/06.9_evaluation/summary/v1"
TICKET = "#9 (Moke110/arc)"

RECORD_FILES = (
    "unit_results.csv",
    "fold_seed_table.csv",
    "per_drug_dose.csv",
    "paired_by_drug.csv",
    "drug_gains.csv",
    "strata.csv",
    "summary.json",
    "manifest.json",
    "record_check.json",
    "tahoe_seen_drugs.json",
)

NEGATIVE_RESULT_CLAUSE = (
    "A negative result is a valid first-round outcome. If the paired Tahoe-vs-random gain "
    "or the advantage over the simple baselines is not positive (95% bootstrap CI over "
    "drugs excluding zero), the study reports that no transfer gain was demonstrated and "
    "that PRnet on this frozen input should not be described as a drug-screening tool. No "
    "tuning, re-selection or seed/fold/LR change is performed after seeing the test "
    "result; the 30 units are the frozen #6/#7 selection."
)

SCOPE_STATEMENT = (
    "This round predicts general 24 h iPSC-CM transcriptomic responses to drugs never "
    "measured in the temp_env CM training set, on the single fixed training-consensus "
    "background c*. It evaluates the 30 checkpoints selected before the test split was "
    "scored. It is not a hypoxia/reoxygenation efficacy claim, not a ScreenSeq or DToxS "
    "result, not a dose-extrapolation claim (the three measured temp_env doses 0.2/1/10 µM "
    "are reported as strata), and not a donor- or batch-specific claim. Low-response and "
    "toxicity labels are unknown in this round and are reported as unknown."
)

DECISION_RULE = (
    "Pre-registered reading: the primary evidence is the drug-level paired mean "
    "difference in direction AUPRC@200 with a 95% percentile bootstrap CI resampling "
    "drugs (B=10000, seed=20261102). 'positive' = point estimate > 0 and CI low > 0; "
    "'negative' = point estimate <= 0; otherwise 'inconclusive'. Both acceptance "
    "questions (vs random init; vs the simple baselines) are reported regardless of sign, "
    "and all three seeds and five folds are shown — no best-seed selection."
)


DEVIATIONS = [
    {
        "id": "variance_head_guard_relaxed",
        "component": (
            "stage2/06_prnet_transfer/06.4_pipeline/execution/ptlib/model.py::"
            "check_architecture"
        ),
        "fix_commit": "813c0bb",
        "what": (
            "The all-zero architecture probe required the decoder variance to be strictly "
            "positive. Four selected Tahoe-init checkpoints underflow softplus to exactly 0 "
            "on that off-distribution input; the guard now rejects only non-finite or "
            "negative variance and accepts an exact zero."
        ),
        "cause": (
            "softplus(raw variance logit) is mathematically > 0 but float32 underflows to 0 "
            "for very negative logits; the fast-learning-rate Tahoe checkpoints retain a "
            "more negative variance head than the random-init ones."
        ),
        "evidence": {
            "weights_non_finite": 0,
            "zero_probe_var_zeros": {
                "e9_f0_s20261023_lr1e-4_tahoe": "464 of 36,328 (2 x 18,164)",
                "e9_f1_s20261023_lr1e-4_tahoe": "300 of 36,328",
                "e9_f4_s20261022_lr1e-4_tahoe": "2 of 36,328",
                "e9_f4_s20261023_lr1e-4_tahoe": "482 of 36,328",
            },
            "real_input_var_min": {
                "e9_f0_s20261023_lr1e-4_tahoe": 2.186267226284011e-12,
                "e9_f1_s20261023_lr1e-4_tahoe": 3.797551162421087e-11,
                "e9_f4_s20261022_lr1e-4_tahoe": 8.789152131788436e-13,
                "e9_f4_s20261023_lr1e-4_tahoe": 1.1615077788462713e-08,
            },
            "real_input_var_zeros": 0,
            "real_input_rows": 32,
            "metric_uses": (
                "predicted mean only: ptlib.predict.predict_rows returns the mean; the "
                "variance is unused by Δ or the direction-AUPRC@200 metric"
            ),
        },
        "unit_ids": [
            "e9_f0_s20261023_lr1e-4_tahoe",
            "e9_f1_s20261023_lr1e-4_tahoe",
            "e9_f4_s20261022_lr1e-4_tahoe",
            "e9_f4_s20261023_lr1e-4_tahoe",
        ],
        "effect_on_results": (
            "none on any predicted value: guard-only change, no model, loss, checkpoint or "
            "selection changed; the 26 already-scored units are byte-identical (eval "
            "manifest + unit record hashes) and only the 4 failed units were re-scored"
        ),
        "recovery": "4/4 re-scored successfully; final 30/30 units completed",
    }
]


class RecordError(common.EvaluationError):
    """The record cannot be derived from the frozen artifacts."""


# --------------------------------------------------------------------------- #
# Derived tables
# --------------------------------------------------------------------------- #
def _git_info() -> dict[str, Any]:
    info: dict[str, Any] = {"sha": None, "branch": None, "dirty": None}
    try:
        for key, command in (
            ("sha", ["git", "-C", str(common.ROOT), "rev-parse", "HEAD"]),
            ("branch", ["git", "-C", str(common.ROOT), "rev-parse", "--abbrev-ref", "HEAD"]),
            ("dirty", ["git", "-C", str(common.ROOT), "status", "--porcelain"]),
        ):
            done = subprocess.run(command, capture_output=True, text=True, timeout=15, check=False)
            if done.returncode == 0:
                if key == "dirty":
                    info[key] = bool(done.stdout.strip())
                else:
                    info[key] = done.stdout.strip()
    except (OSError, subprocess.SubprocessError):  # pragma: no cover - odd git setups
        pass
    return info


def _score_records(out: Path) -> dict[str, dict[str, Any]]:
    record: dict[str, dict[str, Any]] = {}
    units_root = out / "units"
    if not units_root.is_dir():
        return record
    for path in sorted(units_root.glob("*/record.json")):
        try:
            record[path.parent.name] = load_json(path)
        except (OSError, json.JSONDecodeError):
            continue
    return record


def _scored_evals(units_payload: Mapping[str, Any], out: Path) -> tuple[list[Any], list[str]]:
    evals: list[Any] = []
    missing: list[str] = []
    for unit in units_payload["units"]:
        eval_dir = out / "units" / unit["unit_id"] / "eval"
        if not (eval_dir / "aggregate.json").exists():
            missing.append(unit["unit_id"])
            continue
        evals.append(aggregate_mod.load_unit_eval(unit, out))
    return evals, missing


def _models_block(
    unit_results: pd.DataFrame, fold_seed: pd.DataFrame
) -> dict[str, Any]:
    models: dict[str, Any] = {}
    for init in aggregate_mod.MODEL_INITS:
        unit_rows = unit_results[unit_results["init"] == init]
        if unit_rows.empty:
            continue
        per_fold: dict[str, Any] = {}
        for _, row in fold_seed[fold_seed["init"] == init].iterrows():
            per_fold[str(int(row["fold"]))] = {
                "auprc_200": float(row["auprc_200_mean"]),
                "auprc_200_sd_over_seeds": _nan_or_float(row["auprc_200_sd"]),
                "auprc_up_200": float(row["auprc_up_200_mean"]),
                "auprc_down_200": float(row["auprc_down_200_mean"]),
                "rmse": _nan_or_float(row["rmse_mean"]),
                "pcc": _nan_or_float(row["pcc_mean"]),
                "spearman": _nan_or_float(row["spearman_mean"]),
                "n_seeds": int(row["n_seeds"]),
            }
        fold_means = [entry["auprc_200"] for entry in per_fold.values()]
        fold_ups = [entry["auprc_up_200"] for entry in per_fold.values()]
        fold_downs = [entry["auprc_down_200"] for entry in per_fold.values()]
        seed_means = (
            unit_rows.groupby("seed")["auprc_200"].mean().sort_index().to_dict()
        )
        seed_values = np.asarray(list(seed_means.values()), dtype=np.float64)
        models[init] = {
            "n_units": int(len(unit_rows)),
            "macro_auprc_200_over_folds": float(np.mean(fold_means)) if fold_means else None,
            "macro_auprc_up_200_over_folds": float(np.mean(fold_ups)) if fold_ups else None,
            "macro_auprc_down_200_over_folds": (
                float(np.mean(fold_downs)) if fold_downs else None
            ),
            "per_fold": per_fold,
            "seed_stability": {
                "per_seed_mean_over_folds": {
                    str(int(seed)): float(value) for seed, value in seed_means.items()
                },
                "mean": float(np.mean(seed_values)) if seed_values.size else None,
                "sd": float(np.std(seed_values, ddof=1)) if seed_values.size > 1 else 0.0,
            },
            "auxiliary": {
                metric: _nan_or_float(unit_rows[metric].mean())
                for metric in ("rmse", "pcc", "spearman")
            },
            "seed_mean_auprc_200": _nan_or_float(unit_rows["auprc_200"].mean()),
        }
    return models


def _nan_or_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def _primary_bar(
    baselines_summary: pd.DataFrame | None, baselines_record: Mapping[str, Any] | None
) -> dict[str, Any] | None:
    if baselines_summary is None or baselines_summary.empty:
        return None
    macro = baselines_summary[baselines_summary["scope"] == "macro_over_folds"]
    if macro.empty:
        return None
    entry: dict[str, Any] = {"primary": aggregate_mod.PRIMARY_BASELINE, "bars": {}}
    for _, row in macro.iterrows():
        name = str(row["baseline"])
        entry["bars"][name] = {
            "auprc_200": _nan_or_float(row.get("auprc")),
            "auprc_up_200": _nan_or_float(row.get("auprc_up")),
            "auprc_down_200": _nan_or_float(row.get("auprc_down")),
            "rmse": _nan_or_float(row.get("rmse")),
            "pcc": _nan_or_float(row.get("pcc")),
            "spearman": _nan_or_float(row.get("spearman")),
            "role": str(row.get("role", "")),
        }
    primary = entry["bars"].get(aggregate_mod.PRIMARY_BASELINE)
    entry["auprc_200"] = primary.get("auprc_200") if primary else None
    entry["source"] = (baselines_record or {}).get("path")
    entry["source_sha256"] = (baselines_record or {}).get("sha256")
    entry["note"] = (
        "#8's published bar: macro over five folds, identical 951-row / 347-condition "
        "test population"
    )
    return entry


def _test_population(
    evals: Sequence[Any], paired: pd.DataFrame, units_payload: Mapping[str, Any]
) -> dict[str, Any]:
    per_fold: dict[str, Any] = {}
    for unit_eval in evals:
        fold = str(int(unit_eval.fold))
        if fold not in per_fold:
            per_fold[fold] = {
                "n_rows": int(unit_eval.rows.size),
                "n_conditions": int(unit_eval.condition_keys.size),
                "n_drugs": int(
                    len({aggregate_mod.condition_parts(key)[0] for key in unit_eval.condition_keys})
                ),
                "row_sha256": common.sha256_bytes(np.asarray(unit_eval.rows, dtype=np.int64).tobytes()),
                "conditions_sha256": common.sha256_bytes(
                    np.asarray(unit_eval.condition_keys, dtype=str).tobytes()
                ),
            }
    return {
        "n_rows": int(sum(entry["n_rows"] for entry in per_fold.values())),
        "n_conditions": int(sum(entry["n_conditions"] for entry in per_fold.values())),
        "n_drugs": int(len(paired[["fold", "drug_key"]].drop_duplicates()))
        if not paired.empty
        else int(sum(entry["n_drugs"] for entry in per_fold.values())),
        "n_folds": len(per_fold),
        "per_fold": per_fold,
        "source": "frozen evaluator output; identical to #8's test population",
    }


def _budgets(
    unit_results: pd.DataFrame, score_summary: Mapping[str, Any] | None
) -> dict[str, Any]:
    def _stats(values: pd.Series) -> dict[str, Any]:
        numeric = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=np.float64)
        if numeric.size == 0:
            return {"n": 0}
        return {
            "n": int(numeric.size),
            "min": float(numeric.min()),
            "mean": float(numeric.mean()),
            "max": float(numeric.max()),
            "sum": float(numeric.sum()),
        }

    training: dict[str, Any] = {}
    for init in aggregate_mod.MODEL_INITS:
        subset = unit_results[unit_results["init"] == init]
        if subset.empty:
            continue
        stop_reasons = subset["stop_reason"].astype(str).value_counts().to_dict()
        training[init] = {
            "epochs_run": _stats(subset["epochs_run"]),
            "train_wall_s": _stats(subset["train_wall_s"]),
            "train_peak_rss_mb": _stats(subset["train_peak_rss_mb"]),
            "stop_reasons": {str(key): int(value) for key, value in stop_reasons.items()},
            "capped_units": int(
                sum(1 for reason in subset["stop_reason"].astype(str) if "max_epochs" in reason)
            ),
        }
    summary = dict(score_summary or {})
    failures = list(summary.get("failures") or [])
    return {
        "evaluation": {
            "units_planned": int(summary.get("units_planned", len(unit_results))),
            "units_completed": int(summary.get("units_completed", len(unit_results))),
            "units_failed": int(summary.get("units_failed", len(failures))),
            "failures": failures,
            "wall_s_total": summary.get("wall_s"),
            "wall_s_per_unit": _stats(unit_results["eval_wall_s"]),
            "host_peak_rss_mb": summary.get("host_peak_rss_mb")
            or _stats(unit_results["eval_host_peak_rss_mb"]).get("max"),
            "gpu_peak_mb": summary.get("gpu_peak_mb"),
            "predict_repeats": 2,
            "preflight": summary.get("preflight"),
        },
        "training": training,
    }


def _strata_block(
    strata: pd.DataFrame, tahoe_seen_record: Mapping[str, Any] | None
) -> dict[str, Any]:
    has_table = not strata.empty and "stratum_type" in strata.columns
    return {
        "n_rows": int(len(strata)) if has_table else 0,
        "response_label": "unknown",
        "toxicity_label": "unknown",
        "unknown_labels": {
            "response_label": {
                "value": "unknown",
                "n_conditions": int(
                    strata.loc[
                        strata["stratum_type"] == "response_label", "n_conditions"
                    ].iloc[0]
                ) if has_table and (strata["stratum_type"] == "response_label").any() else 0,
                "criterion": (
                    "pre-registered: no per-well response readout exists in this round "
                    "(06.2 §6); response is never inferred from transcriptomic magnitude"
                ),
            },
            "toxicity_label": {
                "value": "unknown",
                "n_conditions": int(
                    strata.loc[
                        strata["stratum_type"] == "toxicity_label", "n_conditions"
                    ].iloc[0]
                ) if has_table and (strata["stratum_type"] == "toxicity_label").any() else 0,
                "criterion": (
                    "pre-registered: no per-well viability/toxicity readout exists in this "
                    "round (06.2 §6); toxicity is never inferred from transcriptomic magnitude"
                ),
            },
        },
        "tahoe_seen": {
            "n_tahoe_seen_drug_keys": (tahoe_seen_record or {}).get("n_tahoe_seen_drug_keys"),
            "definition": "drug_key observed in Tahoe treated rows (structure seen, no CM label)",
            "sha256": (tahoe_seen_record or {}).get("sha256"),
        },
        "dose_strata": sorted(
            strata.loc[strata["stratum_type"] == "dose", "stratum_value"].astype(str).tolist()
        ) if has_table else [],
        "table": "strata.csv",
    }


def _questions(
    drug_gains: pd.DataFrame,
    paired: pd.DataFrame,
    models: Mapping[str, Any],
    primary_bar: Mapping[str, Any] | None,
) -> dict[str, Any]:
    questions: dict[str, Any] = {}
    random_key = "gain_tahoe_minus_random"
    questions["q_a_tahoe_vs_random"] = aggregate_mod.paired_summary(
        drug_gains, random_key
    )
    if f"gain_tahoe_minus_{aggregate_mod.PRIMARY_BASELINE}_mean" in drug_gains.columns:
        questions["q_b_tahoe_vs_primary_baseline"] = aggregate_mod.paired_summary(
            drug_gains, f"gain_tahoe_minus_{aggregate_mod.PRIMARY_BASELINE}"
        )
        questions["q_b_random_vs_primary_baseline"] = aggregate_mod.paired_summary(
            drug_gains, f"gain_random_minus_{aggregate_mod.PRIMARY_BASELINE}"
        )
    if "gain_tahoe_minus_no_change_mean" in drug_gains.columns:
        questions["q_constant_no_change"] = {
            "tahoe_vs_no_change": aggregate_mod.paired_summary(
                drug_gains, "gain_tahoe_minus_no_change"
            ),
            "random_vs_no_change": aggregate_mod.paired_summary(
                drug_gains, "gain_random_minus_no_change"
            ),
            "no_change_auprc_200": (
                (primary_bar or {}).get("bars", {}).get("no_change", {}).get("auprc_200")
            ),
        }
    # seed-level stability of the headline paired gain
    if not paired.empty:
        per_seed = (
            paired.groupby("seed")["gain_tahoe_minus_random"].mean().sort_index().to_dict()
        )
        values = np.asarray(list(per_seed.values()), dtype=np.float64)
        questions["q_a_tahoe_vs_random"]["per_seed_mean_gain"] = {
            str(int(seed)): float(value) for seed, value in per_seed.items()
        }
        questions["q_a_tahoe_vs_random"]["per_seed_gain_sd"] = (
            float(np.std(values, ddof=1)) if values.size > 1 else 0.0
        )
        questions["q_a_tahoe_vs_random"]["model_means"] = {
            init: (models.get(init) or {}).get("macro_auprc_200_over_folds")
            for init in aggregate_mod.MODEL_INITS
        }
        questions["q_a_tahoe_vs_random"]["reading"] = (
            "Does Tahoe pretraining beat random init on unseen-drug CM responses? "
            "Paired by fold+seed, aggregated by drug (doses/seeds/runs are not "
            "independent drugs)."
        )
    return questions


# --------------------------------------------------------------------------- #
# Record construction
# --------------------------------------------------------------------------- #
def build_record(
    *,
    phase: str,
    units_path: Path,
    out: Path,
    data_dir: Path,
    baselines_dir: Path | None,
    allow_subset: bool,
    target: Path,
    make_figures: bool = True,
    tahoe_seen_keys: set[str] | None = None,
) -> dict[str, Any]:
    started = time.perf_counter()
    units_payload = load_unit_manifest(
        units_path, require_complete=(phase == "test" and not allow_subset)
    )
    units = units_payload["units"]
    evals, missing = _scored_evals(units_payload, out)
    if missing and not allow_subset:
        raise RecordError(
            f"{len(missing)} unit(s) have no evaluation output (e.g. {missing[:3]}); "
            "run score.py first, or pass --allow-subset for a labelled dry run"
        )
    if not evals:
        raise RecordError("no scored unit found; run score.py first")
    subset = len(evals) != len(units)
    dry_run = phase != "test" or subset or not units_payload.get("complete")
    if dry_run and target == out / "record":
        raise RecordError(
            "a dry-run/subset record must not be written into output/record/; "
            "it would be indistinguishable from the first-round test record"
        )

    score_summary_path = out / "record" / "score_summary.json"
    score_summary = (
        load_json(score_summary_path) if score_summary_path.exists() else None
    )
    if score_summary and score_summary.get("phase") != phase:
        raise RecordError(
            f"score_summary.json is for phase {score_summary.get('phase')!r}, not {phase!r}"
        )
    score_records = _score_records(out)

    unit_results = aggregate_mod.unit_results_frame(evals, score_records)
    fold_seed = aggregate_mod.fold_seed_frame(unit_results)
    baselines_summary: pd.DataFrame | None = None
    baselines_conditions: pd.DataFrame | None = None
    baselines_record: dict[str, Any] | None = None
    if phase == "test" and baselines_dir is not None:
        baselines_summary, _ = aggregate_mod.load_baselines_summary(baselines_dir)
        baselines_conditions, baselines_record = aggregate_mod.load_baseline_conditions(
            baselines_dir,
            required_folds=[entry["fold"] for entry in units if entry["init_name"] == "random"],
        )
    per_drug_dose = aggregate_mod.per_drug_dose_frame(evals, baselines_conditions)

    paired = pd.DataFrame()
    drug_gains = pd.DataFrame()
    strata = pd.DataFrame()
    tahoe_seen_record: dict[str, Any] | None = None
    if phase == "test" and {"random", "tahoe"} <= {unit_eval.init for unit_eval in evals}:
        try:
            paired = aggregate_mod.paired_by_drug_frame(per_drug_dose)
            drug_gains = aggregate_mod.drug_gains_frame(paired)
            if tahoe_seen_keys is None:
                seen_keys, tahoe_seen_record = aggregate_mod.read_tahoe_seen_drug_keys(data_dir)
            else:
                seen_keys = set(tahoe_seen_keys)
                tahoe_seen_record = {
                    "injected_by_test": True,
                    "n_tahoe_seen_drug_keys": len(seen_keys),
                    "drug_keys": sorted(seen_keys),
                }
            strata = aggregate_mod.strata_frame(
                per_drug_dose, drug_gains, tahoe_seen=seen_keys
            )
        except aggregate_mod.AggregateError as exc:
            if not allow_subset:
                raise
            print(f"warning: paired/strata aggregation skipped ({exc})", file=sys.stderr)

    models = _models_block(unit_results, fold_seed)
    primary_bar = _primary_bar(baselines_summary, baselines_record)
    questions = _questions(drug_gains, paired, models, primary_bar)
    population = _test_population(evals, paired, units_payload)
    budgets = _budgets(unit_results, score_summary)
    summary: dict[str, Any] = {
        "schema": SUMMARY_SCHEMA,
        "kind": "first_round_evaluation",
        "ticket": TICKET,
        "phase": phase,
        "dry_run": bool(dry_run),
        "subset": bool(subset),
        "k": aggregate_mod.K,
        "created_utc": now_utc(),
        "units_manifest": {
            "path": str(units_path),
            "sha256": sha256_file(units_path),
            "created_utc": units_payload.get("created_utc"),
            "complete": bool(units_payload.get("complete")),
            "n_units": len(units),
            "n_scored": len(evals),
        },
        "data_identity": units_payload.get("data_identity"),
        "frozen_manifest_sha256": (units_payload.get("frozen") or {}).get("manifest_sha256"),
        "test_population": population,
        "models": models,
        "primary_bar": primary_bar,
        "questions": questions,
        "strata": _strata_block(strata, tahoe_seen_record),
        "budgets": budgets,
        "deviations": DEVIATIONS,
        "negative_result_clause": NEGATIVE_RESULT_CLAUSE,
        "scope_statement": SCOPE_STATEMENT,
        "decision_rule": DECISION_RULE,
        "bootstrap": {
            "n_bootstrap": aggregate_mod.N_BOOTSTRAP,
            "seed": aggregate_mod.BOOTSTRAP_SEED,
            "unit": "drug (fold pooled; doses/seeds/runs are not independent drugs)",
            "interval": "percentile 95% (2.5, 97.5)",
        },
    }

    target.mkdir(parents=True, exist_ok=True)
    table_hashes: dict[str, str] = {}

    def _write_csv(name: str, frame: pd.DataFrame) -> None:
        path = target / name
        frame.to_csv(path, index=False, float_format="%.10g")
        table_hashes[name] = sha256_file(path)

    _write_csv("unit_results.csv", unit_results)
    _write_csv("fold_seed_table.csv", fold_seed)
    if not per_drug_dose.empty:
        _write_csv("per_drug_dose.csv", per_drug_dose)
    if not paired.empty:
        _write_csv("paired_by_drug.csv", paired)
        _write_csv("drug_gains.csv", drug_gains)
    if not strata.empty:
        _write_csv("strata.csv", strata)
    if tahoe_seen_record is not None:
        write_json(target / "tahoe_seen_drugs.json", tahoe_seen_record)
        table_hashes["tahoe_seen_drugs.json"] = sha256_file(target / "tahoe_seen_drugs.json")

    written_tables = list(table_hashes)
    figures: list[dict[str, Any]] = []
    if make_figures and not dry_run and not fold_seed.empty:
        figures = figures_mod.make_figures(
            fold_seed=fold_seed,
            unit_results=unit_results,
            drug_gains=drug_gains,
            strata=strata,
            summary=summary,
            out_dir=common.DEFAULT_OUT / "figures",
        )

    summary["figures"] = figures
    summary["tables"] = table_hashes
    write_json(target / "summary.json", summary)
    table_hashes["summary.json"] = sha256_file(target / "summary.json")

    expected_files = list(written_tables)
    check = checks_mod.run_gates(
        phase=phase,
        units_payload=units_payload,
        evals=evals,
        unit_results=unit_results,
        per_drug_dose=per_drug_dose,
        paired=paired,
        drug_gains=drug_gains,
        summary=summary,
        baselines_record=baselines_record,
        out=out,
        expected_files=expected_files,
        units_manifest_sha256=sha256_file(units_path),
        record_dir=target,
    )
    write_json(target / "record_check.json", check)
    manifest = {
        "schema": SCHEMA,
        "kind": "evaluation_record",
        "ticket": TICKET,
        "phase": phase,
        "dry_run": bool(dry_run),
        "subset": bool(subset),
        "created_utc": summary["created_utc"],
        "command": ["write_record.py", "--phase", phase],
        "git": _git_info(),
        "units_manifest": summary["units_manifest"],
        "data_identity": summary["data_identity"],
        "frozen": units_payload.get("frozen"),
        "sources": {
            "random_record": units_payload.get("sources", {}).get("random_record"),
            "tahoe_record": units_payload.get("sources", {}).get("tahoe_record"),
            "baselines": baselines_record,
        },
        "metric": {
            "name": "direction_auprc",
            "headline_k": aggregate_mod.K,
            "aggregation": (
                "replicates averaged within condition, direction AP, macro over "
                "conditions; then macro over folds of the three-seed mean"
            ),
            "selection_untouched": "the 30 units are #6/#7's frozen selection",
        },
        "baselines": {
            "path": str(baselines_dir) if baselines_dir else None,
            "record": baselines_record,
            "primary": aggregate_mod.PRIMARY_BASELINE,
        },
        "bootstrap": summary["bootstrap"],
        "deviations": DEVIATIONS,
        "gates": checks_mod.gate_summary(check),
        "record_sha256": table_hashes,
        "figures": figures,
        "environment": _environment(),
        "timings": {"wall_s": time.perf_counter() - started},
        "peak_rss_mb": peak_rss_mb(),
    }
    # manifest.json cannot hash itself: its own sha256 is deliberately absent, and the
    # re-derivation check compares the other tracked files byte for byte.
    write_json(target / "manifest.json", manifest)
    return {
        "summary": summary,
        "check": check,
        "manifest": manifest,
        "unit_results": unit_results,
        "fold_seed": fold_seed,
        "per_drug_dose": per_drug_dose,
        "paired": paired,
        "drug_gains": drug_gains,
        "strata": strata,
        "target": target,
        "dry_run": dry_run,
    }


def _environment() -> dict[str, Any]:
    import platform

    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "python_executable": sys.executable,
    }


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Aggregate the 30 frozen unit evaluations into the ticket-#9 record.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--units", type=Path, default=DEFAULT_OUT / "record" / "units.json"
    )
    parser.add_argument("--phase", choices=("valid", "test"), required=True)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="score output root")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--baselines-dir", type=Path, default=DEFAULT_BASELINES_DIR)
    parser.add_argument(
        "--target", type=Path, default=None,
        help="record directory (default output/record for the full test run, "
             "output/dry_run otherwise)",
    )
    parser.add_argument("--allow-subset", action="store_true",
                        help="allow a partial/valid harness dry run (labelled)")
    parser.add_argument("--no-figures", action="store_true")
    parser.add_argument("--check", action="store_true",
                        help="re-derive and byte-compare against the existing record, write nothing")
    return parser


def main(argv: list[str] | None = None) -> int:
    ns = build_parser().parse_args(argv)
    out = Path(ns.out)
    target = Path(ns.target) if ns.target is not None else None
    if ns.check:
        if target is None:
            target = out / "record"
        if not (target / "summary.json").exists():
            print(f"error: no record to check at {target}", file=sys.stderr)
            return 2
        reference = load_json(target / "summary.json")
        reference_manifest = load_json(target / "manifest.json")
        existing = dict(reference_manifest.get("record_sha256") or {})
        rebuilt = build_record(
            phase=ns.phase,
            units_path=ns.units,
            out=out,
            data_dir=ns.data_dir,
            baselines_dir=ns.baselines_dir,
            allow_subset=ns.allow_subset,
            target=out / "_check_record",
            make_figures=False,
        )
        problems: list[str] = []
        for name, digest in existing.items():
            if name == "summary.json":
                # carries created_utc/plot hashes; compared field-wise below instead
                continue
            path = target / name
            check_path = Path(rebuilt["target"]) / name
            if not path.exists() or not check_path.exists():
                problems.append(f"{name}: missing")
            elif sha256_file(path) != sha256_file(check_path):
                problems.append(f"{name}: content changed on re-derivation")
        # summary carries created_utc/figures; compare its content fields explicitly
        # volatile run evidence: `created_utc` is a timestamp, `figures`/`tables` are
        # output hashes, and `budgets` describes the run (wall time, RSS, GPU, preflight
        # memory) rather than being derivable from the frozen predictions. Everything
        # scientific (identity, population, models, questions, strata, scope) must match.
        volatile = ("created_utc", "figures", "tables", "budgets")
        a = {k: v for k, v in reference.items() if k not in volatile}
        b = {
            k: v
            for k, v in load_json(Path(rebuilt["target"]) / "summary.json").items()
            if k not in volatile
        }
        if a != b:
            problems.append("summary.json: content fields differ on re-derivation")
        if problems:
            print("record check FAILED:", file=sys.stderr)
            for problem in problems:
                print(f"  - {problem}", file=sys.stderr)
            return 1
        print(f"record check passed: {len(existing)} file(s) byte-identical on re-derivation")
        return 0

    if target is None:
        target = out / ("record" if ns.phase == "test" and not ns.allow_subset else "dry_run")
    try:
        result = build_record(
            phase=ns.phase,
            units_path=ns.units,
            out=out,
            data_dir=ns.data_dir,
            baselines_dir=ns.baselines_dir if ns.phase == "test" else None,
            allow_subset=ns.allow_subset,
            target=target,
            make_figures=not ns.no_figures,
        )
    except common.EvaluationError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    summary = result["summary"]
    print(json.dumps(
        {
            "target": str(target),
            "phase": summary["phase"],
            "dry_run": summary["dry_run"],
            "n_scored": summary["units_manifest"]["n_scored"],
            "test_population": {
                key: summary["test_population"][key]
                for key in ("n_rows", "n_conditions", "n_drugs", "n_folds")
            },
            "models": {
                init: block.get("macro_auprc_200_over_folds")
                for init, block in summary["models"].items()
            },
            "primary_bar": (summary.get("primary_bar") or {}).get("auprc_200"),
            "q_a_tahoe_vs_random": {
                key: ((summary["questions"].get("q_a_tahoe_vs_random") or {}).get(key))
                for key in ("overall", "verdict")
            },
            "gates": checks_mod.gate_summary(result["check"]),
        },
        indent=2,
    ))
    if not result["check"]["passed"] and not result["dry_run"]:
        print("error: record gates failed; see record_check.json", file=sys.stderr)
        return 1
    if result["dry_run"]:
        print(
            f"NOTE: dry-run record written to {target} (NOT the first-round test record)",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
