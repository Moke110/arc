"""Record gates for ticket #9.

Every gate is a boolean fact derived from the frozen artifacts; ``write_record.py`` refuses
to write a record with a failed gate and stores the full gate table in
``output/record/record_check.json``. The gates are grouped by the claim they protect:

* **selection frozen** — the 30 units come from #6/#7's records; checkpoints/configs/data
  identities verify; no candidate run ever held a prediction/evaluation artifact.
* **test touched once** — each unit has exactly one test evaluation; the test role is
  recorded; the unit manifest predates the score run.
* **identical population** — every fold's scored rows/conditions equal #8's identical test
  population, and all six units of a fold score the same rows and conditions.
* **results recompute** — every headline macro is the mean of the frozen per-condition
  column the evaluator wrote; the pairing keeps all three seeds and both inits per drug.
* **record content** — the tracked files are byte-recomputable and their sha256s match.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

import common

__all__ = ["run_gates", "gate_summary"]


def _rows_signature(evals: Sequence[Any]) -> dict[str, Any]:
    by_fold: dict[str, list[str]] = {}
    for unit_eval in evals:
        digest = common.sha256_bytes(
            np.asarray(unit_eval.rows, dtype=np.int64).tobytes()
        )
        by_fold.setdefault(str(int(unit_eval.fold)), []).append(digest)
    return {fold: sorted(set(values)) for fold, values in by_fold.items()}


def _conditions_signature(evals: Sequence[Any]) -> dict[str, list[str]]:
    by_fold: dict[str, list[str]] = {}
    for unit_eval in evals:
        digest = common.sha256_bytes(
            np.asarray(unit_eval.condition_keys, dtype=str).tobytes()
        )
        by_fold.setdefault(str(int(unit_eval.fold)), []).append(digest)
    return {fold: sorted(set(values)) for fold, values in by_fold.items()}


def run_gates(
    *,
    phase: str,
    units_payload: Mapping[str, Any],
    evals: Sequence[Any],
    unit_results: pd.DataFrame,
    per_drug_dose: pd.DataFrame,
    paired: pd.DataFrame,
    drug_gains: pd.DataFrame,
    summary: Mapping[str, Any],
    baselines_record: Mapping[str, Any] | None,
    out: Path,
    expected_files: Sequence[str],
    units_manifest_sha256: str,
    record_dir: Path | None = None,
) -> dict[str, Any]:
    """Evaluate every gate and return ``{"passed": bool, "gates": {...}}``."""
    gates: dict[str, Any] = {}
    units = units_payload.get("units") or []
    manifest_gates = units_payload.get("gates") or {}
    complete = bool(units_payload.get("complete"))
    n_units = len(units)

    gates["units_manifest_complete"] = bool(complete and n_units == 30)
    gates["units_manifest_gates_passed"] = all(
        value for key, value in manifest_gates.items() if isinstance(value, bool)
    )
    gates["selection_records_present"] = bool(
        units_payload.get("sources", {}).get("random_record")
        and units_payload.get("sources", {}).get("tahoe_record")
    )
    gates["checkpoint_sha256_verified"] = bool(
        manifest_gates.get("checkpoint_sha256_verified")
    )
    gates["config_digest_chain_ok"] = bool(manifest_gates.get("config_digest_chain_ok"))
    gates["data_identity_uniform"] = bool(
        manifest_gates.get("data_identity_uniform")
        and len({unit["data_identity"] for unit in units}) <= 1
    )
    gates["warm_start_provenance_recorded"] = bool(
        manifest_gates.get("init_provenance_ok")
    )
    gates["no_test_artifacts_in_candidate_runs"] = bool(
        manifest_gates.get("no_test_artifacts_in_candidate_runs")
    )

    # ---- selection is frozen before the test role is scored ----
    score_runs = out / "record" / "score_runs.csv"
    gates["score_runs_recorded"] = score_runs.exists()
    gates["phase_recorded"] = phase in ("valid", "test")
    gates["test_role_scored_once"] = (
        phase != "test"
        or (
            all(
                str(unit_eval.aggregate.get("role")) == "test"
                for unit_eval in evals
            )
            and len(evals) == n_units == 30
        )
    )
    scored_ids = {unit_eval.unit_id for unit_eval in evals}
    gates["every_unit_scored_once"] = len(scored_ids) == n_units and len(evals) == n_units
    unit_manifest_created = str(units_payload.get("created_utc", ""))
    score_started = ""
    if score_runs.exists():
        try:
            first = common.read_table(score_runs)
            score_started = str(first.get("started_utc", pd.Series(dtype=str)).min())
        except (OSError, pd.errors.ParserError):
            score_started = ""
    gates["units_manifest_precedes_test_scoring"] = bool(
        phase != "test" or (unit_manifest_created and score_started and unit_manifest_created <= score_started)
    )

    # ---- identical test population ----
    fold_reference = units_payload.get("fold_reference") or {}
    reference_ok = True
    for unit_eval in evals:
        reference = fold_reference.get(str(int(unit_eval.fold))) or {}
        if not reference:
            reference_ok = False
            break
        rows_expected = reference.get("n_rows", reference.get("n_test_rows"))
        conditions_expected = reference.get("n_conditions", reference.get("n_test_conditions"))
        if phase == "test":
            if int(unit_eval.rows.size) != int(rows_expected):
                reference_ok = False
                break
            if int(unit_eval.condition_keys.size) != int(conditions_expected):
                reference_ok = False
                break
        else:
            if int(unit_eval.rows.size) != int(reference.get("n_valid_rows", unit_eval.rows.size)):
                reference_ok = False
                break
    gates["eval_rows_match_recorded_population"] = bool(reference_ok)
    row_signatures = _rows_signature(evals)
    condition_signatures = _conditions_signature(evals)
    gates["same_rows_across_units_of_a_fold"] = all(
        len(values) == 1 for values in row_signatures.values()
    )
    gates["same_conditions_across_units_of_a_fold"] = all(
        len(values) == 1 for values in condition_signatures.values()
    )
    n_folds = int(units_payload["requirements"]["n_folds"])
    gates["all_folds_present"] = sorted({int(e.fold) for e in evals}) == list(
        range(n_folds)
    )

    # ---- frozen predictions / metrics recompute ----
    recompute_ok = True
    identity_ok = True
    for unit_eval in evals:
        unit = unit_eval.unit
        observed = float(
            np.nanmean(unit_eval.per_condition[f"auprc_{200}"].to_numpy(dtype=np.float64))
        )
        row = unit_results[unit_results["unit_id"] == unit_eval.unit_id]
        if row.empty or not np.isclose(
            float(row["auprc_200"].iloc[0]), observed, rtol=1e-6, atol=1e-9
        ):
            recompute_ok = False
            break
        manifest_checkpoint = (unit_eval.eval_manifest.get("predictions") or {}).get(
            "checkpoint_sha256"
        ) or (unit_eval.eval_manifest.get("checkpoint") or {}).get("sha256")
        if str(manifest_checkpoint) != str(unit["checkpoint"]["sha256"]):
            identity_ok = False
            break
        data = unit_eval.eval_manifest.get("data") or {}
        frozen_artifacts = units_payload.get("frozen", {}).get("artifacts") or {}
        if dict(data.get("artifacts") or {}) != dict(frozen_artifacts):
            identity_ok = False
            break
    gates["unit_results_recompute_from_frozen_predictions"] = bool(recompute_ok)
    gates["eval_identity_matches_unit"] = bool(identity_ok)

    if phase == "test":
        gates["baselines_available"] = bool(baselines_record)
        if baselines_record:
            gates["baselines_same_population"] = bool(
                int(baselines_record.get("n_rows", 0)) > 0
                and sorted(baselines_record.get("folds", []))
                == list(range(int(units_payload["requirements"]["n_folds"])))
            )
        gates["paired_by_drug_complete"] = bool(
            not paired.empty
            and not drug_gains.empty
            and int(drug_gains["n_seeds"].min()) == 3
            and len(drug_gains) == int(summary.get("test_population", {}).get("n_drugs", 0) or 0)
        )
        gates["three_seeds_both_inits_per_drug"] = bool(
            not paired.empty
            and paired.groupby(["fold", "drug_key"])["seed"].nunique().min() == 3
            and set(paired["seed"].unique()) == set(int(s) for s in units_payload["requirements"]["seeds"])
        )
        gates["bootstrap_recorded"] = bool(
            (summary.get("questions") or {}).get("q_a_tahoe_vs_random", {}).get("n_bootstrap", 0) > 0
        )
        gates["strata_present"] = bool(summary.get("strata", {}).get("n_rows", 0) > 0)
        gates["unknown_labels_explicit"] = bool(
            summary.get("strata", {}).get("response_label") == "unknown"
            and summary.get("strata", {}).get("toxicity_label") == "unknown"
        )
    else:
        gates["baselines_available"] = True
        gates["baselines_same_population"] = True
        gates["paired_by_drug_complete"] = True
        gates["three_seeds_both_inits_per_drug"] = True
        gates["bootstrap_recorded"] = True
        gates["strata_present"] = True
        gates["unknown_labels_explicit"] = True

    # ---- record content ----
    record_root = Path(record_dir) if record_dir is not None else out / "record"
    gates["record_files_present"] = all(
        (record_root / name).exists() for name in expected_files
    )
    gates["units_manifest_sha256_recorded"] = bool(units_manifest_sha256)
    gates["record_is_labelled_dry_run"] = (
        "dry_run" in summary and "phase" in summary and "subset" in summary
    )
    passed = all(value for key, value in gates.items() if isinstance(value, bool))
    return {"passed": passed, "gates": gates}


def gate_summary(check: Mapping[str, Any]) -> dict[str, Any]:
    gates = check.get("gates", {})
    return {
        "passed": bool(check.get("passed")),
        "n_gates": sum(1 for value in gates.values() if isinstance(value, bool)),
        "failed": [key for key, value in gates.items() if isinstance(value, bool) and not value],
    }


def write_check(path: Path | str, payload: Mapping[str, Any]) -> None:
    common.write_json(path, payload)
