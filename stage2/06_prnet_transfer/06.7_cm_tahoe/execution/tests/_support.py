"""Shared synthetic fixtures for the 06.7 (ticket #7) unit tests.

The tests never touch the real 16.2 GB expression matrices.  They copy the
**tracked** frozen `locked_config.json` (so the matrix definition, the recipe and
the config fingerprints are the real ones), optionally shrink only the model
*structure* to a tiny one so the warm-start checkpoint can be a real torch model,
and synthesise one run directory per candidate in the exact shape the ticket-#4
trainer writes — including the full ``init`` warm-start provenance block.  That is
enough to exercise the matrix, the strict loader and the record gates without a GPU
or the frozen input.
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
import warmstart  # noqa: E402

from ptlib import checkpoint as pt_checkpoint  # noqa: E402
from ptlib import config as pt_config  # noqa: E402
from ptlib import manifest as pt_manifest  # noqa: E402
from ptlib import model as pt_model  # noqa: E402

import torch  # noqa: E402

REAL_LOCKED_CONFIG = (
    common.REPO_ROOT
    / "stage2"
    / "06_prnet_transfer"
    / "06.3_input"
    / "output"
    / "locked_config.json"
)
REAL_RANDOM_RECORD = (
    common.REPO_ROOT
    / "stage2"
    / "06_prnet_transfer"
    / "06.6_cm_random"
    / "output"
    / "record"
    / "record.json"
)

#: tiny PRnet structure — same builder, same key names, ~40k parameters
SMALL_STRUCTURE: dict[str, Any] = {
    "n_genes": 256,
    "hidden": [16],
    "z_dim": 8,
    "c_dim": 8,
    "n_dim": 4,
    "drug_dim": 32,
    "adaptor": [8],
    "dropout": 0.05,
}


def fake_hash(label: str) -> str:
    return common.sha256_bytes(f"06.7-fixture:{label}".encode("utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def prepare_input(tmp: Path, *, structure: dict[str, Any] | None = None) -> tuple[Path, dict[str, Any], dict[str, Any]]:
    """A synthetic frozen-input directory (real locked config, optional tiny structure)."""
    data_dir = tmp / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    locked = json.loads(REAL_LOCKED_CONFIG.read_text(encoding="utf-8"))
    if structure is not None:
        locked["structure"] = dict(structure)
    write_json(data_dir / "locked_config.json", locked)
    artifacts = {name: fake_hash(name) for name in pt_manifest.REQUIRED_ARTIFACTS}
    write_json(
        data_dir / "manifest.json",
        {"schema": "06.3_input/frozen_manifest/1", "artifacts": artifacts},
    )
    normalized, _notes = common.read_locked_config(data_dir)
    return data_dir, normalized, common.frozen_identity(data_dir)


def _perturb(model: "pt_model.PGM", seed: int) -> None:
    """Make the saved weights differ from ``random_init(seed)`` (a trained checkpoint)."""
    generator = torch.Generator().manual_seed(int(seed) + 7919)
    with torch.no_grad():
        for parameter in model.parameters():
            parameter.add_(torch.randn(parameter.shape, generator=generator) * 0.05)
        for buffer in model.buffers():
            if buffer.dtype.is_floating_point:
                buffer.add_(torch.randn(buffer.shape, generator=generator) * 0.01)


def make_checkpoints(
    tmp: Path,
    *,
    structure: dict[str, Any] = SMALL_STRUCTURE,
    seeds: tuple[int, ...] = (20261021, 20261022, 20261023),
) -> Path:
    """Write one real tiny torch checkpoint per seed in ticket #5's layout."""
    root = tmp / "checkpoints"
    for seed in seeds:
        model = pt_model.build_model(structure)
        pt_model.random_init(model, seed)
        _perturb(model, seed)
        path = warmstart.checkpoint_path(root, seed)
        model_config = {
            "structure": dict(structure),
            "recipe": {"lr": 1e-3},
            "metric": {"name": "direction_auprc", "plateau_mode": "max", "k": 200},
            "stage": "pretrain",
            "train_drugs": "main",
            "seed": int(seed),
            "init": "random",
        }
        pt_checkpoint.save_checkpoint(
            path,
            model,
            config=model_config,
            seed=int(seed),
            n_genes=int(structure["n_genes"]),
            best_valid_auprc=0.12 + 0.001 * (int(seed) % 10),
            epoch=2,
        )
    return root


def default_auprc(candidate: matrix_mod.Candidate) -> float:
    """Decreasing in LR within every fold: the 1e-4 candidate always wins."""
    lr_index = {1e-4: 0, 3e-4: 1, 1e-3: 2}[float(candidate.lr)]
    return 0.10 + 0.01 * candidate.fold - 0.005 * lr_index - 0.001 * candidate.seed_index


def config_for(candidate: matrix_mod.Candidate, locked: dict[str, Any]) -> dict[str, Any]:
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
        "init": matrix_mod.TRAIN_INIT_FLAG,
    }


def init_record(checkpoint_record: dict[str, Any]) -> dict[str, Any]:
    """The ``manifest.init`` block the ticket-#4 trainer writes for a warm start."""
    return {
        "mode": "pretrained",
        "source_path": checkpoint_record["path"],
        "checkpoint_sha256": checkpoint_record["sha256"],
        "checkpoint_epoch": checkpoint_record["checkpoint_epoch"],
        "checkpoint_best_valid_auprc": checkpoint_record["checkpoint_best_valid_auprc"],
        "checkpoint_seed": checkpoint_record["payload_seed"],
        "checkpoint_config_sha256": checkpoint_record["checkpoint_config_sha256"],
        "checkpoint_structure": checkpoint_record["checkpoint_structure"],
        "gene_axis": checkpoint_record["gene_axis"],
        "structure_matches_run": True,
        "load_provenance": checkpoint_record["load_provenance"],
        "parameter_count": checkpoint_record["parameter_count"],
    }


def build_queue_fixture(
    tmp: Path,
    *,
    structure: dict[str, Any] | None = SMALL_STRUCTURE,
) -> dict[str, Any]:
    """Synthetic input + checkpoints + the real queue builder -> queue payload."""
    import contextlib
    import io

    import build_queue  # noqa: E402

    data_dir, locked, frozen = prepare_input(tmp, structure=structure)
    checkpoint_root = make_checkpoints(tmp, structure=locked["structure"])
    out_dir = tmp / "output"
    with contextlib.redirect_stdout(io.StringIO()):
        payload = build_queue.build_queue(
            data_dir=data_dir,
            out_dir=out_dir,
            checkpoint_root=checkpoint_root,
            queue_path=out_dir / "units.json",
            hash_cache=None,
        )
    spec = matrix_mod.matrix_spec(locked)
    return {
        "data_dir": data_dir,
        "out_dir": out_dir,
        "locked": locked,
        "frozen": frozen,
        "spec": spec,
        "candidates": matrix_mod.candidates(spec),
        "queue": payload,
        "checkpoints": payload["warmstart"]["checkpoints"],
        "checkpoint_root": checkpoint_root,
    }


def _history(best_auprc: float, epochs: int, *, peak_at: int | None = None) -> list[dict[str, Any]]:
    peak = int(peak_at if peak_at is not None else epochs)
    entries = []
    for epoch in range(1, epochs + 1):
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
    fixture: dict[str, Any],
    auprc: float,
    epochs: int = 22,
    capped: bool = False,
    peak_at: int | None = None,
    config_mutator: Callable[[dict[str, Any]], None] | None = None,
    data_mutator: Callable[[dict[str, Any]], None] | None = None,
    init_mutator: Callable[[dict[str, Any]], None] | None = None,
    omit_manifest: bool = False,
    extra_entries: tuple[str, ...] = (),
) -> Path:
    """Write one synthetic ticket-#4 run directory (config/history/manifest/ckpt)."""
    locked = fixture["locked"]
    frozen = fixture["frozen"]
    checkpoint_record = fixture["checkpoints"][str(candidate.seed)]
    run_dir = runs_root / candidate.unit_id
    (run_dir / "checkpoint").mkdir(parents=True, exist_ok=True)
    checkpoint_bytes = f"fake-weights:{candidate.unit_id}".encode("utf-8")
    checkpoint_path = run_dir / "checkpoint" / "model.pt"
    checkpoint_path.write_bytes(checkpoint_bytes)
    checkpoint_sha = common.sha256_file(checkpoint_path)

    config = config_for(candidate, locked)
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
    (run_dir / "run.log").write_text(
        "warm start: loaded the fixture checkpoint\n", encoding="utf-8"
    )
    (run_dir / "model_card.md").write_text("# fixture\n", encoding="utf-8")

    data = {
        "manifest_sha256": frozen["manifest_sha256"],
        "locked_config": {"sha256": frozen["locked_config_sha256"]},
        "artifacts": dict(frozen["artifacts"]),
        "hash_verification": {"verified": len(frozen["artifacts"]), "cached": []},
    }
    if data_mutator is not None:
        data_mutator(data)
    init = init_record(checkpoint_record)
    if init_mutator is not None:
        init_mutator(init)
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
        "init": init,
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


def write_random_record(tmp: Path, fixture: dict[str, Any]) -> Path:
    """A ticket-#6 record stub with the same matrix/recipe and ``init='random'``.

    Only the fields the ticket-#7 symmetry gate reads are present; the real
    ticket-#6 record is not needed (and would carry the full 18,164-gene structure
    that these synthetic fixtures intentionally shrink).
    """
    locked = fixture["locked"]
    spec = fixture["spec"]
    path = tmp / "random_record.json"
    write_json(path, {
        "schema": "cm_random.record/v1",
        "matrix": {
            "init": "random",
            "n_folds": spec["n_folds"],
            "folds": spec["folds"],
            "seeds": spec["seeds"],
            "lrs": spec["lrs"],
            "n_candidates": spec["n_candidates"],
            "candidate_order": spec["candidate_order"],
        },
        "recipe": {
            "recipe": dict(locked["recipe"]),
            "structure": dict(locked["structure"]),
            "metric": dict(locked["metric"]),
            "loss": locked["loss"],
            "optimizer": locked["optimizer"],
            "prediction_noise": locked["prediction_noise"],
        },
        "selection_rule": {
            "rule": "per fold, choose the LR maximising the mean over the three seeds of the "
                    "best validation direction AUPRC@200 (k=200)",
            "tie_break": "lower LR, then ascending-LR candidate order",
        },
    })
    return path


def prepare_out(
    tmp: Path,
    *,
    auprc: Callable[[matrix_mod.Candidate], float] = default_auprc,
    capped: frozenset[str] = frozenset(),
    failed: frozenset[str] = frozenset(),
    **run_kwargs: Any,
) -> dict[str, Any]:
    """A complete 45-run output tree + queue status, ready for the record writer."""
    fixture = build_queue_fixture(tmp)
    locked = fixture["locked"]
    frozen = fixture["frozen"]
    cands = fixture["candidates"]
    out_dir = fixture["out_dir"]
    runs_root = out_dir / "runs"
    runs_root.mkdir(parents=True, exist_ok=True)

    status_rows = []
    for candidate in cands:
        if candidate.unit_id in failed:
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
            fixture=fixture,
            auprc=float(auprc(candidate)),
            capped=candidate.unit_id in capped,
            **run_kwargs,
        )
        status_rows.append({
            "unit_id": candidate.unit_id, "stage": "finetune", "status": "completed",
            "exit_code": "0", "started_utc": "", "finished_utc": "",
            "wall_s": "42.5", "run_dir": str(runs_root / candidate.unit_id), "error": "",
        })

    with open(out_dir / "queue_status.csv", "w", newline="", encoding="utf-8") as handle:
        import csv as _csv

        writer = _csv.DictWriter(handle, fieldnames=list(status_rows[0]))
        writer.writeheader()
        writer.writerows(status_rows)

    write_json(out_dir / "record" / "matrix_run.json", {
        "schema": "cm_tahoe.matrix_run/v1",
        "status": "completed",
        "run": {"wall_s": 1234.5, "exit_code": 0, "queue_summary": {"completed": len(cands)},
                "units_failed": []},
        "peaks": {"host_process_tree_rss_mib": 2800.0, "gpu_used_mib": 4100.0},
        "warmstart": {"provenance_sha256": fixture["queue"]["warmstart"]["provenance_sha256"]},
    })
    write_json(out_dir / "record" / "resume_check.json", {
        "schema": "cm_tahoe.resume_check/v1",
        "passed": True,
        "checks": [{"id": "finished_matrix_all_skipped", "ok": True}],
    })
    fixture.update({
        "runs_root": runs_root,
        "status_rows": status_rows,
        "random_record_path": write_random_record(tmp, fixture),
    })
    return fixture
