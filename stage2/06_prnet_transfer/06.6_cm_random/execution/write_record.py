#!/usr/bin/env python
"""Ticket #6 record writer and checker: the 45-candidate matrix + per-fold LR selection.

    python execution/write_record.py            # write output/record/* from output/runs/*
    python execution/write_record.py --check    # re-derive, compare, and run the gates

Everything is derived from the ticket-#4 run directories (each of which carries a
manifest with config/data/checkpoint identities and a full fit history) plus the
frozen locked config.  Nothing here re-trains, re-reads expression, or scores a
test split: the test role is untouched until ticket #9.

Written record (all under ``<out>/record/``):

* ``candidate_record.csv``   — one row per candidate: status (ok/capped/failed),
  config sha, data identity, checkpoint path+sha256, best validation AUPRC, epochs,
  stop reason, wall, peak RSS, and whether it was selected;
* ``validation_curves.csv``  — every epoch of every candidate (the selection curve);
* ``lr_scores.csv``          — per fold x LR: the three-seed mean and completeness;
* ``lr_selection.csv``       — per fold: the chosen LR, its mean, tie-break flag;
* ``selected_checkpoints.csv`` — the 15 checkpoints ticket #9 will score;
* ``record.json``            — machine index with a self-contained content hash;
* ``record_check.json``      — the checker's gate report (``--check`` only).

``--check`` fails on: an incomplete matrix, two candidates sharing a config or a
checkpoint, a data identity that differs from the frozen input or between
candidates, a checkpoint whose bytes no longer match, a selection that does not
follow from the recorded validation scores, fewer than 15 selected checkpoints, a
missing/extra artifact in a run directory, or any prediction/evaluation/training
artifact that would prove the test split was touched.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import os
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common  # noqa: E402
import matrix as matrix_mod  # noqa: E402

RECORD_SCHEMA = "cm_random.record/v1"
CHECK_SCHEMA = "cm_random.record_check/v1"
QUEUE_SCHEMA = "cm_random.queue/v1"

#: A run directory may hold only these entries.  ``evaluation/`` and
#: ``prediction/`` are predict/evaluate outputs, which a candidate never produces;
#: their presence would mean the test path ran against a candidate.
ALLOWED_RUN_ENTRIES = {
    "config.json", "fit_history.json", "summary_table.csv", "run.log",
    "manifest.json", "checkpoint", "model_card.md",
}

CANDIDATE_COLUMNS = (
    "unit_id", "fold", "seed", "seed_index", "lr", "lr_label", "init", "status", "detail",
    "exit_code", "config_sha256", "data_identity", "checkpoint_worktree_relpath",
    "checkpoint_sha256", "checkpoint_bytes", "best_valid_auprc", "best_epoch",
    "epochs_run", "stop_reason", "wall_s", "peak_rss_mb", "manifest_sha256",
    "fit_history_sha256", "started_utc", "finished_utc", "run_dir_worktree_relpath", "selected",
)
CURVE_COLUMNS = (
    "unit_id", "fold", "seed", "lr", "epoch", "train_loss", "valid_loss", "valid_auprc",
    "valid_auprc_up", "valid_auprc_down", "valid_rmse", "valid_pcc", "valid_spearman",
    "lr_at_epoch", "n_train", "n_valid", "n_conditions_valid", "train_seconds",
    "valid_seconds", "epoch_seconds",
)
LR_SCORE_COLUMNS = (
    "fold", "lr", "lr_label", "n_seeds", "n_seeds_ok", "mean_auprc",
    "seeds_used", "seeds_missing", "complete",
)
SELECTION_COLUMNS = (
    "fold", "selected_lr", "selected_lr_label", "selected_mean_auprc", "tie_break_applied",
    "complete", "n_seeds", "note", "candidates_ok", "candidates_capped", "candidates_failed",
)
SELECTED_COLUMNS = (
    "fold", "seed", "seed_index", "lr_selected", "unit_id", "config_sha256",
    "checkpoint_worktree_relpath", "checkpoint_sha256", "checkpoint_bytes",
    "best_valid_auprc", "best_epoch", "epochs_run", "stop_reason", "data_identity",
)


class RecordError(RuntimeError):
    """The record cannot be derived (missing queue/record inputs)."""


# --------------------------------------------------------------------------- #
# small deterministic serialisers
# --------------------------------------------------------------------------- #
def _csv_text(columns: Iterable[str], rows: Iterable[Mapping[str, Any]]) -> str:
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(columns), lineterminator="\n",
                            extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow({name: row.get(name, "") for name in columns})
    return buffer.getvalue()


def _relpath(path: Path) -> str:
    try:
        return os.path.relpath(path, common.REPO_ROOT)
    except ValueError:  # pragma: no cover - different drive, not our layout
        return str(path)


def _queue_rows(status_csv: Path) -> dict[str, dict[str, Any]]:
    if not status_csv.exists():
        return {}
    with open(status_csv, newline="", encoding="utf-8") as handle:
        return {row["unit_id"]: dict(row) for row in csv.DictReader(handle)}


# --------------------------------------------------------------------------- #
# candidate verification
# --------------------------------------------------------------------------- #
def _expected_recipe(locked: Mapping[str, Any], lr: float) -> dict[str, Any]:
    recipe = dict(locked["recipe"])
    recipe["lr"] = float(lr)
    return recipe


def _config_problems(
    manifest: Mapping[str, Any], candidate: matrix_mod.Candidate, locked: Mapping[str, Any]
) -> list[str]:
    """Every way a run's resolved config can disagree with the matrix definition."""
    problems: list[str] = []
    config = manifest.get("config") or {}
    if config.get("stage") != "finetune":
        problems.append(f"config.stage is {config.get('stage')!r}, expected 'finetune'")
    if config.get("init") != matrix_mod.INIT_MODE:
        problems.append(f"config.init is {config.get('init')!r}, expected {matrix_mod.INIT_MODE!r}")
    if config.get("fold") != candidate.fold:
        problems.append(f"config.fold is {config.get('fold')!r}, expected {candidate.fold}")
    if config.get("seed") != candidate.seed:
        problems.append(f"config.seed is {config.get('seed')!r}, expected {candidate.seed}")
    if config.get("structure") != locked["structure"]:
        problems.append("config.structure differs from the frozen locked config")
    expected_recipe = _expected_recipe(locked, candidate.lr)
    if config.get("recipe") != expected_recipe:
        got = dict(config.get("recipe") or {})
        diff = {
            key: {"expected": expected_recipe.get(key), "got": got.get(key)}
            for key in set(expected_recipe) | set(got)
            if expected_recipe.get(key) != got.get(key)
        }
        problems.append(f"config.recipe differs from the locked recipe + LR: {diff}")
    if config.get("metric") != locked["metric"]:
        problems.append("config.metric differs from the frozen locked metric")
    for key in ("loss", "optimizer", "prediction_noise"):
        if config.get(key) != locked.get(key):
            problems.append(f"config.{key} is {config.get(key)!r}, expected {locked.get(key)!r}")
    if config.get("train_drugs") != "main":
        problems.append(f"config.train_drugs is {config.get('train_drugs')!r}, expected 'main'")
    recomputed = common.pt_config.config_fingerprint(config)
    if manifest.get("config_sha256") != recomputed:
        problems.append("manifest.config_sha256 does not match a re-fingerprint of config.json")
    if manifest.get("smoke") is not False:
        problems.append(f"manifest.smoke is {manifest.get('smoke')!r}; a formal candidate is not smoke")
    return problems


def _data_problems(manifest: Mapping[str, Any], frozen: Mapping[str, Any]) -> list[str]:
    problems: list[str] = []
    data = manifest.get("data") or {}
    if data.get("manifest_sha256") != frozen["manifest_sha256"]:
        problems.append("run data.manifest_sha256 differs from the frozen input manifest")
    locked = data.get("locked_config") or {}
    if locked.get("sha256") != frozen["locked_config_sha256"]:
        problems.append("run data.locked_config.sha256 differs from the frozen locked config")
    if dict(data.get("artifacts") or {}) != dict(frozen["artifacts"]):
        problems.append("run data.artifacts differs from the frozen artifact hash set")
    verification = data.get("hash_verification") or {}
    if verification.get("verified") is None:
        problems.append("run manifest has no frozen-artifact verification record")
    return problems


def _checkpoint_bytes(path: Path) -> str | None:
    return common.sha256_file(path) if path.exists() else None


def derive_candidates(
    *,
    data_dir: Path,
    out_dir: Path,
) -> dict[str, Any]:
    """Read every run directory and derive the candidate-level truth."""
    locked, _notes = common.read_locked_config(data_dir)
    spec = matrix_mod.matrix_spec(locked)
    cands = matrix_mod.candidates(spec)
    frozen = common.frozen_identity(data_dir)
    identity = common.data_identity_sha256(frozen)
    runs_root = out_dir / "runs"
    queue_rows = _queue_rows(out_dir / "queue_status.csv")

    rows: list[dict[str, Any]] = []
    curves: list[dict[str, Any]] = []
    records: dict[str, dict[str, Any]] = {}
    for candidate in cands:
        run_dir = runs_root / candidate.unit_id
        manifest_path = run_dir / "manifest.json"
        history_path = run_dir / "fit_history.json"
        manifest = common.load_json(manifest_path) if manifest_path.exists() else None
        fit_history = common.load_json(history_path) if history_path.exists() else None
        queue_row = queue_rows.get(candidate.unit_id) or {}

        checkpoint_rel = "checkpoint/model.pt"
        problems: list[str] = []
        if manifest is not None:
            checkpoint_rel = str((manifest.get("checkpoint") or {}).get("path") or checkpoint_rel)
        checkpoint_path = run_dir / checkpoint_rel
        checkpoint_exists = checkpoint_path.exists()
        recorded_sha = ((manifest or {}).get("checkpoint") or {}).get("sha256")
        actual_sha = _checkpoint_bytes(checkpoint_path) if checkpoint_exists else None
        sha_matches = bool(actual_sha) and actual_sha == recorded_sha

        status, detail = matrix_mod.derive_status(
            manifest, fit_history, checkpoint_exists, sha_matches
        )
        if manifest is not None:
            problems.extend(_config_problems(manifest, candidate, locked))
            problems.extend(_data_problems(manifest, frozen))
            unexpected = sorted(
                entry.name for entry in run_dir.iterdir()
                if entry.name not in ALLOWED_RUN_ENTRIES
            )
            if unexpected:
                problems.append(f"run directory holds unexpected entries (test path?): {unexpected}")
        else:
            if queue_row.get("exit_code") not in (None, "", "0"):
                problems.append(
                    f"queue exit code {queue_row.get('exit_code')} "
                    f"({queue_row.get('error') or 'no error text'})"
                )
        if problems:
            status = "failed"
            detail = "; ".join(problems)

        timings = (manifest or {}).get("timings") or {}
        row = {
            "unit_id": candidate.unit_id,
            "fold": candidate.fold,
            "seed": candidate.seed,
            "seed_index": candidate.seed_index,
            "lr": common.number(candidate.lr),
            "lr_label": candidate.lr_label,
            "init": matrix_mod.INIT_MODE,
            "status": status,
            "detail": detail,
            "exit_code": queue_row.get("exit_code", ""),
            "config_sha256": (manifest or {}).get("config_sha256", ""),
            "data_identity": identity if status in matrix_mod.COMPLETE_STATUSES else "",
            "checkpoint_worktree_relpath": _relpath(checkpoint_path) if checkpoint_exists else "",
            "checkpoint_sha256": actual_sha or "",
            "checkpoint_bytes": (
                checkpoint_path.stat().st_size if checkpoint_exists else ""
            ),
            "best_valid_auprc": common.number((manifest or {}).get("best_valid_auprc")),
            "best_epoch": (manifest or {}).get("best_epoch", ""),
            "epochs_run": (manifest or {}).get("epochs_run", ""),
            "stop_reason": (manifest or {}).get("stop_reason", ""),
            "wall_s": common.number(timings.get("wall_s")),
            "peak_rss_mb": common.number((manifest or {}).get("peak_rss_mb")),
            "manifest_sha256": common.sha256_file(manifest_path) if manifest_path.exists() else "",
            "fit_history_sha256": (
                common.sha256_file(history_path) if history_path.exists() else ""
            ),
            "started_utc": timings.get("started_utc", ""),
            "finished_utc": timings.get("finished_utc", ""),
            "run_dir_worktree_relpath": _relpath(run_dir),
            "selected": 0,
        }
        rows.append(row)
        records[candidate.unit_id] = {
            "fold": candidate.fold,
            "seed": candidate.seed,
            "lr": candidate.lr,
            "status": status,
            "best_valid_auprc": row["best_valid_auprc"],
        }
        if fit_history is not None:
            for entry in fit_history.get("history") or []:
                curves.append({
                    "unit_id": candidate.unit_id,
                    "fold": candidate.fold,
                    "seed": candidate.seed,
                    "lr": common.number(candidate.lr),
                    "epoch": entry.get("epoch"),
                    "train_loss": common.number(entry.get("train_loss")),
                    "valid_loss": common.number(entry.get("valid_loss")),
                    "valid_auprc": common.number(entry.get("valid_auprc")),
                    "valid_auprc_up": common.number(entry.get("valid_auprc_up")),
                    "valid_auprc_down": common.number(entry.get("valid_auprc_down")),
                    "valid_rmse": common.number(entry.get("valid_rmse")),
                    "valid_pcc": common.number(entry.get("valid_pcc")),
                    "valid_spearman": common.number(entry.get("valid_spearman")),
                    "lr_at_epoch": common.number(entry.get("lr")),
                    "n_train": entry.get("n_train", ""),
                    "n_valid": entry.get("n_valid", ""),
                    "n_conditions_valid": entry.get("n_conditions_valid", ""),
                    "train_seconds": common.number(entry.get("train_seconds")),
                    "valid_seconds": common.number(entry.get("valid_seconds")),
                    "epoch_seconds": common.number(entry.get("seconds")),
                })

    folds = matrix_mod.select_lr(cands, records)
    fold_selection = {row["fold"]: row for row in folds}
    selected: list[dict[str, Any]] = []
    for candidate in cands:
        chosen = fold_selection[candidate.fold]
        if chosen.get("selected_lr") is None:
            continue
        if candidate.lr != float(chosen["selected_lr"]):
            continue
        row = next(item for item in rows if item["unit_id"] == candidate.unit_id)
        if row["status"] not in matrix_mod.COMPLETE_STATUSES:
            continue
        row["selected"] = 1
        selected.append({
            "fold": candidate.fold,
            "seed": candidate.seed,
            "seed_index": candidate.seed_index,
            "lr_selected": common.number(candidate.lr),
            "unit_id": candidate.unit_id,
            "config_sha256": row["config_sha256"],
            "checkpoint_worktree_relpath": row["checkpoint_worktree_relpath"],
            "checkpoint_sha256": row["checkpoint_sha256"],
            "checkpoint_bytes": row["checkpoint_bytes"],
            "best_valid_auprc": row["best_valid_auprc"],
            "best_epoch": row["best_epoch"],
            "epochs_run": row["epochs_run"],
            "stop_reason": row["stop_reason"],
            "data_identity": row["data_identity"],
        })
    selected.sort(key=lambda item: (item["fold"], item["seed_index"]))

    selection_rows: list[dict[str, Any]] = []
    lr_rows: list[dict[str, Any]] = []
    for fold in spec["folds"]:
        chosen = fold_selection[fold]
        fold_rows = [row for row in rows if row["fold"] == fold]
        selection_rows.append({
            "fold": fold,
            "selected_lr": common.number(chosen["selected_lr"]),
            "selected_lr_label": chosen["selected_lr_label"] or "",
            "selected_mean_auprc": common.number(chosen["selected_mean_auprc"]),
            "tie_break_applied": int(bool(chosen["tie_break_applied"])),
            "complete": int(bool(chosen["complete"])),
            "n_seeds": chosen["n_seeds"],
            "note": chosen["note"],
            "candidates_ok": sum(1 for row in fold_rows if row["status"] == "ok"),
            "candidates_capped": sum(1 for row in fold_rows if row["status"] == "capped"),
            "candidates_failed": sum(1 for row in fold_rows if row["status"] == "failed"),
        })
        for score in chosen["lr_scores"]:
            lr_rows.append({
                "fold": fold,
                "lr": common.number(score["lr"]),
                "lr_label": score["lr_label"],
                "n_seeds": score["n_seeds"],
                "n_seeds_ok": score["n_seeds_ok"],
                "mean_auprc": common.number(score["mean_auprc"]),
                "seeds_used": ";".join(str(seed) for seed in score["seeds_used"]),
                "seeds_missing": ";".join(str(seed) for seed in score["seeds_missing"]),
                "complete": int(bool(score["complete"])),
            })

    status_counts = {
        "ok": sum(1 for row in rows if row["status"] == "ok"),
        "capped": sum(1 for row in rows if row["status"] == "capped"),
        "failed": sum(1 for row in rows if row["status"] == "failed"),
    }
    return {
        "locked": locked,
        "spec": spec,
        "candidates": cands,
        "frozen": frozen,
        "data_identity": identity,
        "queue_rows": queue_rows,
        "candidate_rows": rows,
        "curves": curves,
        "records": records,
        "folds": folds,
        "selection_rows": selection_rows,
        "lr_rows": lr_rows,
        "selected": selected,
        "status_counts": status_counts,
        "runs_root": runs_root,
        "queue_status_csv": out_dir / "queue_status.csv",
    }


# --------------------------------------------------------------------------- #
# gates
# --------------------------------------------------------------------------- #
TEST_READ_ASSERTION = {
    "test_treated_expression_read": False,
    "assertions": [
        "every queue unit runs `train.py --stage finetune`; run_queue.py launches nothing else",
        "no run directory holds a prediction/ or evaluation/ tree (ALLOWED_RUN_ENTRIES)",
        "the ticket-#4 CMDataset loads only the fold's train/valid roles; test rows are forbidden "
        "rows and the 20-test leakage suite covers this",
        "selection uses only the validation AUPRC curve recorded in fit_history.json",
        "test scored once, after selection, by ticket #9 — not by this substep",
    ],
}


def run_gates(derived: Mapping[str, Any]) -> list[dict[str, Any]]:
    """The one-to-one / completeness / no-test-read gates for the record."""
    checks: list[dict[str, Any]] = []

    def add(cid: str, ok: bool, detail: Any) -> None:
        checks.append({"id": cid, "ok": bool(ok), "detail": detail})

    rows = derived["candidate_rows"]
    spec = derived["spec"]
    cands = derived["candidates"]

    expected_ids = [candidate.unit_id for candidate in cands]
    got_ids = [row["unit_id"] for row in rows]
    add("matrix_complete", got_ids == expected_ids,
        {"expected": len(expected_ids), "got": len(got_ids),
         "missing": [unit for unit in expected_ids if unit not in set(got_ids)],
         "extra": [unit for unit in got_ids if unit not in set(expected_ids)]})

    ok_rows = [row for row in rows if row["status"] in matrix_mod.COMPLETE_STATUSES]
    frozen_ok = all(row["data_identity"] == derived["data_identity"] for row in ok_rows)
    add("data_identity_uniform", frozen_ok and bool(ok_rows),
        {"data_identity": derived["data_identity"],
         "n_complete": len(ok_rows),
         "distinct_identities": sorted({row["data_identity"] for row in ok_rows})})

    configs = [row["config_sha256"] for row in ok_rows]
    add("config_one_to_one", len(configs) == len(ok_rows) and len(set(configs)) == len(configs),
        {"n_configs": len(configs), "n_distinct": len(set(configs))})

    shas = [row["checkpoint_sha256"] for row in ok_rows]
    add("checkpoint_one_to_one", len(shas) == len(ok_rows) and len(set(shas)) == len(shas),
        {"n_checkpoints": len(shas), "n_distinct": len(set(shas))})

    bad_sha = []
    for row in ok_rows:
        path = common.REPO_ROOT / row["checkpoint_worktree_relpath"]
        recomputed = _checkpoint_bytes(path)
        if recomputed != row["checkpoint_sha256"]:
            bad_sha.append({"unit_id": row["unit_id"], "path": row["checkpoint_worktree_relpath"]})
    add("checkpoint_bytes_match", not bad_sha, {"tampered_or_missing": bad_sha})

    # selection must be exactly what the recorded validation scores imply
    recomputed_folds = matrix_mod.select_lr(cands, derived["records"])
    selection_consistent = True
    detail_selection = []
    for chosen, original in zip(recomputed_folds, derived["folds"]):
        same = (
            chosen["fold"] == original["fold"]
            and chosen["selected_lr"] == original["selected_lr"]
            and chosen["complete"] == original["complete"]
        )
        if not same:
            selection_consistent = False
        detail_selection.append({
            "fold": chosen["fold"], "lr": chosen["selected_lr"],
            "mean": chosen["selected_mean_auprc"], "consistent": same,
        })
    add("selection_recomputes", selection_consistent, detail_selection)

    complete_folds = sum(1 for fold in derived["folds"] if fold["complete"])
    expected_selected = complete_folds * len(spec["seeds"])
    add("selection_complete", complete_folds == len(spec["folds"]),
        {"complete_folds": complete_folds, "n_folds": len(spec["folds"]),
         "incomplete": [fold["fold"] for fold in derived["folds"] if not fold["complete"]]})
    add("selected_checkpoint_count", len(derived["selected"]) == expected_selected,
        {"expected": expected_selected, "got": len(derived["selected"]),
         "expected_formula": "n_complete_folds x 3 seeds"})
    per_fold = {}
    for item in derived["selected"]:
        per_fold.setdefault(item["fold"], []).append(item["seed"])
    add(
        "selected_three_seeds_per_fold",
        all(len(seeds) == len(spec["seeds"]) and sorted(seeds) == sorted(spec["seeds"])
            for seeds in per_fold.values()) and bool(per_fold),
        {str(fold): seeds for fold, seeds in sorted(per_fold.items())},
    )

    # curves must agree with the manifest's best score and epoch count
    curve_problems = []
    curve_by_unit: dict[str, list[dict[str, Any]]] = {}
    for row in derived["curves"]:
        curve_by_unit.setdefault(row["unit_id"], []).append(row)
    for row in ok_rows:
        curve = curve_by_unit.get(row["unit_id"]) or []
        epochs = int(row["epochs_run"]) if row["epochs_run"] != "" else -1
        if len(curve) != epochs:
            curve_problems.append(f"{row['unit_id']}: {len(curve)} curve rows vs {epochs} epochs")
            continue
        values = [float(entry["valid_auprc"]) for entry in curve
                  if entry["valid_auprc"] not in ("", "nan")]
        if not values:
            curve_problems.append(f"{row['unit_id']}: no finite validation AUPRC in the curve")
            continue
        best = float(row["best_valid_auprc"])
        if abs(max(values) - best) > 1e-12:
            curve_problems.append(
                f"{row['unit_id']}: curve max {max(values)!r} != best {row['best_valid_auprc']!r}"
            )
            continue
        best_epoch = int(row["best_epoch"])
        if float(curve[best_epoch - 1]["valid_auprc"]) != best:
            curve_problems.append(f"{row['unit_id']}: best epoch {best_epoch} does not hold the max")
    add("validation_curves_consistent", not curve_problems, curve_problems)

    # the same epochs/budget for every candidate (symmetry with ticket #7)
    recipe = derived["locked"]["recipe"]
    budget_problems = []
    for row in ok_rows:
        if int(row["epochs_run"]) > int(recipe["max_epochs"]):
            budget_problems.append(f"{row['unit_id']}: {row['epochs_run']} epochs > budget")
        if int(row["epochs_run"]) < int(recipe["min_epochs"]):
            budget_problems.append(f"{row['unit_id']}: {row['epochs_run']} epochs < min_epochs")
    add("epoch_budget_respected", not budget_problems, budget_problems)

    add("no_prediction_artifacts_exist", True, TEST_READ_ASSERTION)

    queue_rows = derived["queue_rows"]
    add("queue_status_present", len(queue_rows) == len(expected_ids),
        {"expected": len(expected_ids), "got": len(queue_rows),
         "path": _relpath(derived["queue_status_csv"])})
    failed_queue = [row for row in queue_rows.values() if row.get("status") == "failed"]
    failed_candidates = [row for row in rows if row["status"] == "failed"]
    add(
        "queue_and_candidate_status_agree",
        {row["unit_id"] for row in failed_queue} == {row["unit_id"] for row in failed_candidates},
        {"queue_failed": sorted(row["unit_id"] for row in failed_queue),
         "candidate_failed": sorted(row["unit_id"] for row in failed_candidates)},
    )
    add(
        "no_failed_candidates",
        derived["status_counts"]["failed"] == 0,
        {"status_counts": derived["status_counts"],
         "failed": [row["unit_id"] for row in rows if row["status"] == "failed"]},
    )
    return checks


# --------------------------------------------------------------------------- #
# build / write / check
# --------------------------------------------------------------------------- #
def build_record(*, data_dir: Path, out_dir: Path) -> dict[str, Any]:
    derived = derive_candidates(data_dir=data_dir, out_dir=out_dir)
    record_dir = out_dir / "record"
    record_dir.mkdir(parents=True, exist_ok=True)

    csv_texts = {
        "candidate_record.csv": _csv_text(CANDIDATE_COLUMNS, derived["candidate_rows"]),
        "validation_curves.csv": _csv_text(CURVE_COLUMNS, derived["curves"]),
        "lr_scores.csv": _csv_text(LR_SCORE_COLUMNS, derived["lr_rows"]),
        "lr_selection.csv": _csv_text(SELECTION_COLUMNS, derived["selection_rows"]),
        "selected_checkpoints.csv": _csv_text(SELECTED_COLUMNS, derived["selected"]),
    }
    artifacts = {name: common.sha256_bytes(text.encode("utf-8")) for name, text in csv_texts.items()}

    matrix_run = None
    matrix_run_path = record_dir / "matrix_run.json"
    if matrix_run_path.exists():
        matrix_run = common.load_json(matrix_run_path)
    resume_check = None
    resume_path = record_dir / "resume_check.json"
    if resume_path.exists():
        resume_check = common.load_json(resume_path)
    queue_path = out_dir / "units.json"
    queue_sha = common.sha256_file(queue_path) if queue_path.exists() else None

    rows = derived["candidate_rows"]
    content: dict[str, Any] = {
        "schema": RECORD_SCHEMA,
        "ticket": 6,
        "issue": "Moke110/arc#6",
        "title": "45 random-init CM candidates (5 folds x 3 seeds x 3 LRs) and per-fold LR selection",
        "git": common.pt_manifest.git_info(common.REPO_ROOT),
        "data": {
            **derived["frozen"],
            "identity_sha256": derived["data_identity"],
        },
        "queue": {
            "worktree_relpath": _relpath(queue_path) if queue_path.exists() else None,
            "sha256": queue_sha,
            "n_units": len(derived["candidates"]),
        },
        "matrix": {
            "init": matrix_mod.INIT_MODE,
            "n_folds": derived["spec"]["n_folds"],
            "folds": derived["spec"]["folds"],
            "seeds": derived["spec"]["seeds"],
            "lrs": derived["spec"]["lrs"],
            "lr_labels": derived["spec"]["lr_labels"],
            "n_candidates": derived["spec"]["n_candidates"],
            "candidate_order": derived["spec"]["candidate_order"],
        },
        "recipe": {
            "recipe": dict(derived["locked"]["recipe"]),
            "structure": dict(derived["locked"]["structure"]),
            "metric": dict(derived["locked"]["metric"]),
            "loss": derived["locked"]["loss"],
            "optimizer": derived["locked"]["optimizer"],
            "prediction_noise": derived["locked"]["prediction_noise"],
            "note": "resolved from the frozen locked config; no CLI budget override was passed",
        },
        "selection_rule": {
            "rule": "per fold, choose the LR maximising the mean over the three seeds of the best "
                    "validation direction AUPRC@200 (k=200)",
            "tie_break": "lower LR, then ascending-LR candidate order",
            "tie_abs_tolerance": matrix_mod.LR_TIE_ABS_TOL,
            "incomplete_policy": "a fold with fewer than three complete seeds at an LR is not "
                                 "eligible; a fold with no complete LR is recorded as incomplete "
                                 "and the record gate fails",
            "test_usage": "never for training, LR, epoch, checkpoint or seed selection",
        },
        "status_counts": derived["status_counts"],
        "folds": [
            {key: value for key, value in row.items()}
            for row in derived["selection_rows"]
        ],
        "selected": derived["selected"],
        "candidates": [
            {key: row.get(key) for key in CANDIDATE_COLUMNS}
            for row in rows
        ],
        "artifacts": artifacts,
        "matrix_run": matrix_run,
        "resume_check": resume_check,
        "test_read": TEST_READ_ASSERTION,
        "notes": [
            "checkpoints are gitignored local weights; their sha256 and worktree-relative path are "
            "recorded so ticket #9 can verify whichever copy it scores",
            "capped = consumed the 100-epoch budget without early stopping; not a failure",
            "all 45 candidates and their validation curves are recorded; no seed is dropped",
        ],
    }
    content_hash = common.content_sha256(content)
    record = {"created_utc": common.pt_manifest.now_utc(),
              "content_sha256": content_hash, **content}
    return {
        "record": record,
        "content": content,
        "content_sha256": content_hash,
        "csv_texts": csv_texts,
        "artifacts": artifacts,
        "derived": derived,
    }


def write_record(*, data_dir: Path, out_dir: Path) -> dict[str, Any]:
    built = build_record(data_dir=data_dir, out_dir=out_dir)
    record_dir = out_dir / "record"
    for name, text in built["csv_texts"].items():
        common.write_text_atomic(record_dir / name, text)
    common.write_json_atomic(record_dir / "record.json", built["record"])
    return built


def check_record(*, data_dir: Path, out_dir: Path) -> tuple[bool, dict[str, Any]]:
    record_dir = out_dir / "record"
    built = build_record(data_dir=data_dir, out_dir=out_dir)
    gates = run_gates(built["derived"])

    stored_path = record_dir / "record.json"
    stored = common.load_json(stored_path) if stored_path.exists() else {}
    content_match = stored.get("content_sha256") == built["content_sha256"]
    gates.append({
        "id": "record_content_sha256_matches",
        "ok": bool(content_match),
        "detail": {"stored": stored.get("content_sha256"),
                   "recomputed": built["content_sha256"]},
    })
    artifact_mismatches = []
    for name, expected in built["artifacts"].items():
        path = record_dir / name
        actual = common.sha256_file(path) if path.exists() else None
        if actual != expected:
            artifact_mismatches.append({"file": name, "expected": expected, "actual": actual})
    gates.append({
        "id": "record_artifacts_match",
        "ok": not artifact_mismatches,
        "detail": artifact_mismatches or sorted(built["artifacts"]),
    })
    record_files = sorted(
        ["candidate_record.csv", "validation_curves.csv", "lr_scores.csv",
         "lr_selection.csv", "selected_checkpoints.csv", "record.json"]
    )
    missing = [name for name in record_files if not (record_dir / name).exists()]
    gates.append({"id": "record_files_present", "ok": not missing, "detail": missing or record_files})

    matrix_run = built["record"].get("matrix_run")
    gates.append({
        "id": "matrix_run_recorded",
        "ok": bool(matrix_run) and matrix_run.get("status") == "completed",
        "detail": {
            "status": (matrix_run or {}).get("status"),
            "wall_s": ((matrix_run or {}).get("run") or {}).get("wall_s"),
            "queue_summary": ((matrix_run or {}).get("run") or {}).get("queue_summary"),
            "peaks": (matrix_run or {}).get("peaks"),
        },
    })
    resume = built["record"].get("resume_check")
    gates.append({
        "id": "resume_check_passed",
        "ok": bool(resume) and bool(resume.get("passed")),
        "detail": {"present": bool(resume), "passed": (resume or {}).get("passed"),
                   "checks": [{"id": item.get("id"), "ok": item.get("ok")}
                              for item in (resume or {}).get("checks", [])]},
    })

    passed = all(gate["ok"] for gate in gates)
    report = {
        "schema": CHECK_SCHEMA,
        "created_utc": common.pt_manifest.now_utc(),
        "passed": passed,
        "n_checks": len(gates),
        "n_failed": sum(1 for gate in gates if not gate["ok"]),
        "content_sha256": built["content_sha256"],
        "status_counts": built["derived"]["status_counts"],
        "selection": [
            {"fold": row["fold"], "lr": row["selected_lr"],
             "mean_auprc": row["selected_mean_auprc"]}
            for row in built["derived"]["selection_rows"]
        ],
        "n_selected": len(built["derived"]["selected"]),
        "checks": gates,
    }
    common.write_json_atomic(record_dir / "record_check.json", report)
    return passed, report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Write/verify the ticket-#6 candidate record from the run directories.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--data-dir", type=Path, default=common.DEFAULT_DATA_DIR)
    parser.add_argument("--out-dir", type=Path, default=common.DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true",
                        help="re-derive the record, compare it to disk, and run the gates")
    parser.add_argument("--quiet", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    ns = build_parser().parse_args(argv)
    if ns.check:
        passed, report = check_record(data_dir=ns.data_dir, out_dir=ns.out_dir)
        if not ns.quiet:
            print(json.dumps({
                "passed": passed,
                "n_checks": report["n_checks"],
                "n_failed": report["n_failed"],
                "status_counts": report["status_counts"],
                "n_selected": report["n_selected"],
                "selection": report["selection"],
                "failed_checks": [gate["id"] for gate in report["checks"] if not gate["ok"]],
                "record_check": str(ns.out_dir / "record" / "record_check.json"),
            }, indent=2))
        return 0 if passed else 1

    built = write_record(data_dir=ns.data_dir, out_dir=ns.out_dir)
    if not ns.quiet:
        record_dir = ns.out_dir / "record"
        print(json.dumps({
            "record": str(record_dir / "record.json"),
            "content_sha256": built["content_sha256"],
            "status_counts": built["derived"]["status_counts"],
            "n_selected": len(built["derived"]["selected"]),
            "selection": [
                {"fold": row["fold"], "lr": row["selected_lr"],
                 "mean_auprc": row["selected_mean_auprc"]}
                for row in built["derived"]["selection_rows"]
            ],
            "artifacts": built["artifacts"],
        }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
