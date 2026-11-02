"""Shared synthetic-run fixtures for the 06.6 unit tests.

The tests never touch the real frozen expression matrices: they copy the
**tracked** frozen `locked_config.json` (so the matrix definition and the config
fingerprint are the real ones), write a synthetic `manifest.json` for the frozen
artifacts, and synthesise one run directory per candidate in the exact shape the
ticket-#4 trainer writes.  That is enough to exercise the record writer and every
gate without a GPU or the 16 GB input.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path
from typing import Any, Callable

EXECUTION = Path(__file__).resolve().parent.parent
if str(EXECUTION) not in sys.path:
    sys.path.insert(0, str(EXECUTION))

import common  # noqa: E402
import matrix as matrix_mod  # noqa: E402

from ptlib import config as pt_config  # noqa: E402
from ptlib import manifest as pt_manifest  # noqa: E402

REAL_LOCKED_CONFIG = (
    common.REPO_ROOT / "stage2" / "06_prnet_transfer" / "06.3_input" / "output" / "locked_config.json"
)


def fake_hash(label: str) -> str:
    return common.sha256_bytes(f"06.6-fixture:{label}".encode("utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def prepare_input(tmp: Path) -> tuple[Path, dict[str, Any]]:
    """A synthetic frozen-input directory + its identity (real locked config)."""
    data_dir = tmp / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(REAL_LOCKED_CONFIG, data_dir / "locked_config.json")
    artifacts = {name: fake_hash(name) for name in pt_manifest.REQUIRED_ARTIFACTS}
    write_json(
        data_dir / "manifest.json",
        {
            "schema": "06.3_input/frozen_manifest/1",
            "artifacts": artifacts,
        },
    )
    return data_dir, common.frozen_identity(data_dir)


def default_auprc(candidate: matrix_mod.Candidate) -> float:
    """Decreasing in LR within every fold: the 1e-4 candidate always wins."""
    lr_index = {1e-4: 0, 3e-4: 1, 1e-3: 2}[float(candidate.lr)]
    return 0.10 + 0.01 * candidate.fold - 0.005 * lr_index - 0.001 * candidate.seed_index


def _config_for(candidate: matrix_mod.Candidate, locked: dict[str, Any]) -> dict[str, Any]:
    return {
        "structure": dict(locked["structure"]),
        "recipe": {**locked["recipe"], "lr": float(candidate.lr)},
        "metric": dict(locked["metric"]),
        "loss": locked["loss"],
        "optimizer": locked["optimizer"],
        "dose_encoding": dict(locked.get("dose_encoding") or {}),
        "prediction_noise": locked["prediction_noise"],
        "train_drugs": "main",
        "stage": "finetune",
        "seed": candidate.seed,
        "fold": candidate.fold,
        "init": matrix_mod.INIT_MODE,
    }


def _history(best_auprc: float, epochs: int, *, peak_at: int | None = None) -> list[dict[str, Any]]:
    peak = int(peak_at if peak_at is not None else epochs)
    entries = []
    for epoch in range(1, epochs + 1):
        # linear ramp to the peak epoch, then a small decay
        if epoch <= peak:
            auprc = best_auprc - 0.002 * (peak - epoch)
        else:
            auprc = best_auprc - 0.001 * (epoch - peak)
        entries.append({
            "epoch": epoch,
            "train_loss": 1.5 - 0.01 * epoch,
            "valid_loss": 1.4 - 0.01 * epoch,
            "valid_auprc": auprc,
            "valid_auprc_up": auprc * 1.05,
            "valid_auprc_down": auprc * 0.95,
            "valid_rmse": 0.4,
            "valid_pcc": 0.05,
            "valid_spearman": 0.03,
            "lr": 3e-4,
            "train_seconds": 0.7,
            "valid_seconds": 0.4,
            "n_train": 571,
            "n_valid": 188,
            "n_conditions_valid": 63,
            "seconds": 1.1,
        })
    return entries


def write_run(
    runs_root: Path,
    candidate: matrix_mod.Candidate,
    *,
    locked: dict[str, Any],
    frozen: dict[str, Any],
    auprc: float,
    epochs: int = 22,
    capped: bool = False,
    peak_at: int | None = None,
    config_mutator: Callable[[dict[str, Any]], None] | None = None,
    data_mutator: Callable[[dict[str, Any]], None] | None = None,
    omit_manifest: bool = False,
    extra_entries: tuple[str, ...] = (),
) -> Path:
    run_dir = runs_root / candidate.unit_id
    (run_dir / "checkpoint").mkdir(parents=True, exist_ok=True)
    checkpoint_bytes = f"fake-weights:{candidate.unit_id}".encode("utf-8")
    checkpoint_path = run_dir / "checkpoint" / "model.pt"
    checkpoint_path.write_bytes(checkpoint_bytes)
    checkpoint_sha = common.sha256_file(checkpoint_path)

    config = _config_for(candidate, locked)
    if config_mutator is not None:
        config_mutator(config)
    config_sha = pt_config.config_fingerprint(config)
    history = _history(auprc, epochs, peak_at=peak_at)
    best_epoch = peak_at if peak_at is not None else epochs
    stop_reason = "max_epochs" if capped else "early_stop_no_auprc_improvement_20"

    write_json(run_dir / "config.json", config)
    write_json(run_dir / "fit_history.json", {
        "run": {"stop_reason": stop_reason, "best_valid_auprc": auprc},
        "history": history,
    })
    (run_dir / "summary_table.csv").write_text("epoch,valid_auprc\n", encoding="utf-8")
    (run_dir / "run.log").write_text("fixture run\n", encoding="utf-8")
    (run_dir / "model_card.md").write_text("# fixture\n", encoding="utf-8")

    data = {
        "manifest_sha256": frozen["manifest_sha256"],
        "locked_config": {"sha256": frozen["locked_config_sha256"]},
        "artifacts": dict(frozen["artifacts"]),
        "hash_verification": {"verified": len(frozen["artifacts"]), "cached": []},
    }
    if data_mutator is not None:
        data_mutator(data)
    manifest = {
        "schema": pt_manifest.SCHEMA,
        "kind": "train",
        "stage": "finetune",
        "status": "completed",
        "smoke": False,
        "config": config,
        "config_sha256": config_sha,
        "data": data,
        "fold": candidate.fold,
        "seed": candidate.seed,
        "init": {"mode": matrix_mod.INIT_MODE},
        "checkpoint": {
            "path": "checkpoint/model.pt",
            "sha256": checkpoint_sha,
            "reload_identical": True,
            "best_valid_auprc": auprc,
            "best_epoch": best_epoch,
        },
        "best_valid_auprc": auprc,
        "best_epoch": best_epoch,
        "epochs_run": len(history),
        "stop_reason": stop_reason,
        "peak_rss_mb": 1330.0,
        "timings": {
            "wall_s": 42.5,
            "started_utc": "2026-11-02T07:00:00+00:00",
            "finished_utc": "2026-11-02T07:00:42+00:00",
        },
    }
    if not omit_manifest:
        write_json(run_dir / "manifest.json", manifest)
    for name in extra_entries:
        target = run_dir / name
        target.mkdir(parents=True, exist_ok=True)
        (target / "manifest.json").write_text("{}\n", encoding="utf-8")
    return run_dir


def prepare_out(
    tmp: Path,
    *,
    auprc: Callable[[matrix_mod.Candidate], float] = default_auprc,
    capped: frozenset[str] = frozenset(),
    failed: frozenset[str] = frozenset(),
    **run_kwargs: Any,
) -> dict[str, Any]:
    """A complete 45-run output tree + queue status, ready for the record writer."""
    data_dir, frozen = prepare_input(tmp)
    out_dir = tmp / "output"
    locked, _notes = common.read_locked_config(data_dir)
    spec = matrix_mod.matrix_spec(locked)
    cands = matrix_mod.candidates(spec)
    runs_root = out_dir / "runs"
    runs_root.mkdir(parents=True, exist_ok=True)

    status_rows = []
    for candidate in cands:
        if candidate.unit_id in failed:
            # an interrupted/failed unit: run dir present but no manifest
            run_dir = runs_root / candidate.unit_id
            (run_dir / "checkpoint").mkdir(parents=True, exist_ok=True)
            status_rows.append({
                "unit_id": candidate.unit_id, "stage": "finetune", "status": "failed",
                "exit_code": "1", "started_utc": "", "finished_utc": "", "wall_s": "0",
                "run_dir": str(run_dir), "error": "exit code 1; see logs",
            })
            continue
        write_run(
            runs_root,
            candidate,
            locked=locked,
            frozen=frozen,
            auprc=float(auprc(candidate)),
            capped=candidate.unit_id in capped,
            **run_kwargs,
        )
        status_rows.append({
            "unit_id": candidate.unit_id, "stage": "finetune", "status": "completed",
            "exit_code": "0", "started_utc": "", "finished_utc": "",
            "wall_s": "42.5", "run_dir": str(runs_root / candidate.unit_id), "error": "",
        })

    # units.json: the real builder, pointed at the synthetic data dir
    import contextlib
    import io

    import build_queue  # noqa: E402  (local import: uses ptlib.queue at module load)

    with contextlib.redirect_stdout(io.StringIO()):
        build_queue.main([
            "--data-dir", str(data_dir),
            "--out-dir", str(out_dir),
            "--hash-cache", "none",
        ])
    with open(out_dir / "queue_status.csv", "w", newline="", encoding="utf-8") as handle:
        import csv as _csv

        writer = _csv.DictWriter(handle, fieldnames=list(status_rows[0]))
        writer.writeheader()
        writer.writerows(status_rows)

    # run-level records the gates require
    write_json(out_dir / "record" / "matrix_run.json", {
        "schema": "cm_random.matrix_run/v1",
        "status": "completed",
        "run": {"wall_s": 1234.5, "exit_code": 0, "queue_summary": {"completed": len(cands)},
                "units_failed": []},
        "peaks": {"host_process_tree_rss_mib": 2800.0, "gpu_used_mib": 4100.0},
    })
    write_json(out_dir / "record" / "resume_check.json", {
        "schema": "cm_random.resume_check/v1",
        "passed": True,
        "checks": [{"id": "finished_matrix_all_skipped", "ok": True}],
    })
    return {
        "data_dir": data_dir,
        "out_dir": out_dir,
        "frozen": frozen,
        "locked": locked,
        "spec": spec,
        "candidates": cands,
        "runs_root": runs_root,
        "status_rows": status_rows,
    }
