"""Synthetic fixtures for the 06.5 record/checker tests (no GPU, no frozen data).

Every fixture is a miniature of the real thing: a frozen 06.3 input directory
with a manifest/locked-config/contract, and one or three fake 06.4 run
directories with a manifest/config/fit_history/model_card/checkpoint.  The
numbers mirror the real frozen contract (18,164 genes, 50,520 train / 8,524
valid Tahoe rows, 649 control groups, plate13+plate14 valid) so a test failure
also flags a contract drift.
"""

from __future__ import annotations

import copy
import json
import sys
from pathlib import Path
from typing import Any, Mapping

EXECUTION = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(EXECUTION))

import record  # noqa: E402

STEP = EXECUTION.parent
REAL_DATA_DIR = STEP.parent / "06.3_input" / "output"
SEEDS = [20261021, 20261022, 20261023]
N_GENES = 18164
N_TRAIN = 50520
N_VALID = 8524
N_GROUPS = 649
N_ROWS = 67018

LOCKED_CONFIG: dict[str, Any] = {
    "schema": "06.3_input/locked_config/1",
    "version": "1.0.0",
    "contract_id": "stage2/06_prnet_transfer/06.3_input/contract@1.0.0",
    "structure": {
        "n_genes": N_GENES,
        "hidden": [128],
        "z_dim": 64,
        "c_dim": 64,
        "n_dim": 10,
        "drug_dim": 1024,
        "adaptor": [128],
        "dropout": 0.05,
    },
    "recipe": {
        "batch_size": 512,
        "max_epochs": 100,
        "min_epochs": 20,
        "early_stop_patience": 20,
        "sched_factor": 0.5,
        "sched_patience": 10,
        "min_lr": 1e-08,
        "weight_decay": 1e-08,
    },
    "metric": {"name": "direction_auprc", "k": 200, "plateau_mode": "max"},
    "loss": "gaussian_nll",
    "optimizer": "adam",
    "prediction_noise": "zero",
    "n_folds": 5,
    "seeds": SEEDS,
    "lrs": {"pretrain": 0.001, "cm": [0.0001, 0.0003, 0.001]},
    "dose_encoding": {"fingerprint": "FCFP4", "n_bits": 1024},
}


def sha256_bytes(payload: bytes) -> str:
    import hashlib

    return hashlib.sha256(payload).hexdigest()


def write_json(path: Path, payload: Mapping) -> None:  # type: ignore[name-defined]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def make_data_dir(root: Path, *, locked: dict[str, Any] | None = None) -> Path:
    """A miniature frozen input directory with real files (so sha256 works)."""
    data_dir = Path(root) / "06.3_input" / "output"
    data_dir.mkdir(parents=True, exist_ok=True)
    tahoe_bytes = b"synthetic tahoe transfer view"
    (data_dir / "tahoe_transfer.h5ad").write_bytes(tahoe_bytes)
    artifacts = {
        name: sha256_bytes(name.encode("utf-8"))
        for name in (
            "cm_transfer.h5ad",
            "folds.parquet",
            "consensus_controls.npz",
            "consensus_membership.parquet",
            "tahoe_control_groups.parquet",
            "tahoe_pairs.npz",
            "locked_config.json",
        )
    }
    artifacts["tahoe_transfer.h5ad"] = sha256_bytes(tahoe_bytes)
    manifest = {
        "artifacts": [
            {"path": name, "role": "synthetic", "sha256": digest, "bytes": 10,
             "tracked": True}
            for name, digest in artifacts.items()
        ],
        "build": {"tahoe_valid_plates": ["plate13", "plate14"], "fold_seed": 20261102,
                  "n_folds": 5},
        "config": {"digest": "f" * 64, "path": "locked_config.json",
                   "sha256": "e" * 64},
        "counts": {
            "tahoe": {
                "rows": N_ROWS,
                "usable_train": N_TRAIN,
                "usable_valid": N_VALID,
                "split_train_rows": 57422,
                "split_valid_rows": 9596,
                "control_groups": N_GROUPS,
                "control_members": 1294,
                "pairs": 60861,
                "controls": 1400,
                "treated": 65618,
            },
            "cm": {},
        },
    }
    contract = {
        "schema": "06.3_input/contract/1",
        "gene_axis": {"n_genes": N_GENES, "identity_columns": ["hgnc_id"],
                      "measurement_mask": "var.gene_measured", "order": "genes.csv"},
        "control_semantics": {
            "tahoe_train_validation": "mean of the same (plate, cell_line) group's QC-pass "
                                      "control rows on the log layer; never across cell lines"
        },
    }
    write_json(data_dir / "locked_config.json", locked or LOCKED_CONFIG)
    write_json(data_dir / "manifest.json", manifest)
    write_json(data_dir / "contract.json", contract)
    return data_dir


def make_run_dir(
    root: Path,
    *,
    seed: int,
    unit_id: str | None = None,
    budget: str = "formal",
    config_sha256: str = "d" * 64,
    epochs: int = 100,
    stop_reason: str = "max_epochs",
    best_epoch: int = 87,
    best_valid_auprc: float = 0.0612,
    n_train: int = N_TRAIN,
    n_valid: int = N_VALID,
    control_groups: int = N_GROUPS,
    peak_rss_mb: float = 6000.0,
    wall_s: float = 20000.0,
    recipe: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Write a fake finished run directory and return the (seed, unit_id, dir)."""
    unit_id = unit_id or f"tahoe_seed{seed}"
    run_dir = Path(root) / "runs" / unit_id
    run_dir.mkdir(parents=True, exist_ok=True)
    effective_recipe = dict(LOCKED_CONFIG["recipe"]) | {"lr": 0.001}
    if recipe:
        effective_recipe.update(recipe)
    checkpoint_dir = run_dir / "checkpoint"
    checkpoint_dir.mkdir(exist_ok=True)
    checkpoint_bytes = f"checkpoint for seed {seed} {budget}".encode("utf-8")
    (checkpoint_dir / "model.pt").write_bytes(checkpoint_bytes)
    checkpoint_sha = sha256_bytes(checkpoint_bytes)

    # a plausible curve: rise, oscillate, then stale after the best epoch
    curve = []
    for epoch in range(1, epochs + 1):
        if epoch == best_epoch:
            auprc = best_valid_auprc
        else:
            auprc = min(
                0.019 + 0.0004 * min(epoch, best_epoch - 1) + (0.0002 if epoch % 3 == 0 else 0.0),
                best_valid_auprc - 0.001,
            )
        curve.append(
            {
                "epoch": epoch,
                "train_loss": -2.0 - 0.001 * epoch,
                "valid_loss": -1.5 - 0.0005 * epoch,
                "valid_auprc": auprc,
                "valid_auprc_up": auprc,
                "valid_auprc_down": auprc,
                "valid_rmse": 0.4,
                "valid_pcc": 0.2,
                "valid_spearman": 0.1,
                "lr": 0.001,
                "train_seconds": 120.0,
                "valid_seconds": 45.0,
                "seconds": 165.0,
                "n_train": n_train,
                "n_valid": n_valid,
            }
        )
    config_payload = {
        "structure": dict(LOCKED_CONFIG["structure"]),
        "recipe": effective_recipe,
        "metric": dict(LOCKED_CONFIG["metric"]),
        "loss": "gaussian_nll",
        "optimizer": "adam",
        "prediction_noise": "zero",
        "train_drugs": "main",
        "stage": "pretrain",
        "seed": seed,
        "fold": None,
        "init": "random",
        "config_sha256": config_sha256,
        "run": {"stage": "pretrain", "seed": seed, "smoke": budget == "smoke"},
    }
    manifest = {
        "kind": "train",
        "stage": "pretrain",
        "status": "completed",
        "created_utc": "2026-11-02T00:00:00+00:00",
        "command": ["train.py", "--stage", "pretrain", "--seed", str(seed)],
        "git": {"sha": "deadbeef", "branch": "feat/x", "dirty": False},
        "data": {
            "artifacts": {},
            "gene_axis": {"n_genes": N_GENES, "measured": N_GENES - 1},
            "hash_verification": {"verified": 8, "cached": [], "bytes": 15_300_000_000},
            "manifest_sha256": "a" * 64,
        },
        "config": config_payload,
        "config_sha256": config_sha256,
        "structure": dict(LOCKED_CONFIG["structure"]),
        "seed": seed,
        "fold": None,
        "device": {"requested": "cuda", "resolved": "cuda", "name": "synthetic GPU"},
        "init": {"mode": "random", "parameter_count": 7185000},
        "checkpoint": {
            "path": "checkpoint/model.pt",
            "sha256": checkpoint_sha,
            "reload_identical": True,
            "best_valid_auprc": best_valid_auprc,
            "best_epoch": best_epoch,
        },
        "noise_policy": {"train": "seeded", "validation": "zero", "inference": "zero"},
        "timings": {"wall_s": wall_s, "train_loop_s": wall_s - 500,
                    "per_epoch_train_s": [120.0] * epochs},
        "peak_rss_mb": peak_rss_mb,
        "epochs_run": epochs,
        "stop_reason": stop_reason,
        "best_valid_auprc": best_valid_auprc,
        "best_epoch": best_epoch,
        "metric": dict(LOCKED_CONFIG["metric"]),
        "data_counts": {
            "n_train": n_train,
            "n_valid": n_valid,
            "control_groups": control_groups,
            "smoke_subset": budget == "smoke",
        },
        "smoke": budget == "smoke",
    }
    write_json(run_dir / "manifest.json", manifest)
    write_json(run_dir / "config.json", config_payload)
    write_json(run_dir / "fit_history.json", {"run": {"seed": seed}, "history": curve})
    (run_dir / "model_card.md").write_text(
        f"# PRnet pretrain — model card\n\n- seed {seed}\n- checkpoint sha256 {checkpoint_sha}\n",
        encoding="utf-8",
    )
    return {
        "unit_id": unit_id,
        "seed": seed,
        "budget": budget,
        "run_dir": run_dir,
        "config_sha256": config_sha256,
        "checkpoint_sha256": checkpoint_sha,
    }


def input_identity(data_dir: Path) -> dict[str, Any]:
    identity = record.frozen_input_identity(data_dir)
    tahoe_sha = identity["artifacts"]["tahoe_transfer.h5ad"]
    identity["cache"] = {
        "path": "/synthetic/cache.npy",
        "source_sha256": tahoe_sha,
        "layer_sha256": "c" * 64,
        "file_sha256": "b" * 64,
        "shape": [N_ROWS, N_GENES],
        "dtype": "float32",
        "bytes": 4869259936,
        "sample_max_abs_diff": 0.0,
        "schema": "ptlib.tahoe_log1p_cache/v1",
    }
    identity["table_cache"] = {"path": "/synthetic/input_cache", "manifest_sha256": "9" * 64,
                               "tables": {}}
    return identity


def locked_from(data_dir: Path) -> dict[str, Any]:
    return record.normalize_locked_config(record.read_json(data_dir / "locked_config.json"))


def make_record(
    root: Path,
    *,
    data_dir: Path,
    seeds: list[int] | None = None,
    mode: str = "formal",
    locked_sha: str | None = None,
) -> dict[str, Any]:
    """A valid record with three (or one) completed units."""
    seeds = seeds or SEEDS
    identities = [
        make_run_dir(root, seed=seed, config_sha256=("d%d" % (index + 1)).ljust(64, "0"))
        for index, seed in enumerate(seeds)
    ]
    units = []
    for item in identities:
        entry = record.run_entry(
            run_dir=item["run_dir"],
            unit_id=item["unit_id"],
            seed=item["seed"],
            budget="formal",
            command=["python", "train.py", "--stage", "pretrain", "--seed", str(item["seed"])],
            exit_code=0,
            expected_config_sha256=item["config_sha256"],
            wall_s_outer=20100.0,
            resources={
                "host_peak_rss_mb_pipeline": None,
                "host_peak_rss_mb_sampler": 6050.0,
                "gpu_peak_mib_device": 3050.0,
                "gpu_peak_mib_process": 2900.0,
                "gpu_baseline_mib": 1050.0,
                "gpu_peak_mib_device_minus_baseline": 2000.0,
                "gpu_after_mib": 1050.0,
                "samples": 20000,
                "interval_s": 1.0,
            },
        )
        units.append(entry)
    expected = {
        "locked_config_sha256": locked_sha or record.sha256_file(data_dir / "locked_config.json"),
        "structure": dict(LOCKED_CONFIG["structure"]),
        "budget": dict(LOCKED_CONFIG["recipe"]) | {"lr": 0.001},
        "metric": dict(LOCKED_CONFIG["metric"]),
        "loss": "gaussian_nll",
        "optimizer": "adam",
        "prediction_noise": "zero",
        "n_folds": 5,
        "seeds": seeds,
        "lrs": dict(LOCKED_CONFIG["lrs"]),
        "config_sha256": {item["unit_id"]: item["config_sha256"] for item in identities},
        "config_sha256_by_seed": {str(item["seed"]): item["config_sha256"] for item in identities},
    }
    budget = {
        "kind": "locked",
        "smoke": False,
        "row_subset": False,
        "profile_budgets": None,
        "overrides": {},
        "overrides_by_unit": {},
        **(dict(LOCKED_CONFIG["recipe"]) | {"lr": 0.001}),
    }
    return record.build_record(
        mode=mode,
        units=units,
        budget=budget,
        expectations=expected,
        input_identity=input_identity(data_dir),
        command=["python", "run_pretrain.py", "--mode", mode],
        git={"sha": "deadbeef", "branch": "feat/x", "dirty": False},
        seeds_planned=seeds,
        seeds_requested=seeds,
        host={"python": "synthetic", "torch": "2.0.1+cu118"},
    )


def copy_record(payload: dict[str, Any]) -> dict[str, Any]:
    return copy.deepcopy(payload)


def refresh_input_shas(payload: dict[str, Any]) -> dict[str, Any]:
    """Recompute the manifest/locked hashes after a fixture edit (for tamper tests)."""
    data_dir = Path(payload["input"]["data_dir"])
    payload["input"]["manifest_sha256"] = record.sha256_file(data_dir / "manifest.json")
    payload["input"]["locked_config_sha256"] = record.sha256_file(data_dir / "locked_config.json")
    return payload


def make_profile_record(
    root: Path,
    *,
    data_dir: Path,
    budget: str = "smoke",
) -> dict[str, Any]:
    """A single-unit profile record (bounded smoke or one full-width epoch)."""
    seed = SEEDS[0]
    unit_id = f"tahoe_seed{seed}_{budget}"
    if budget == "smoke":
        overrides = {"max_epochs": 2, "min_epochs": 1}
        epochs, n_train, n_valid, stop = 2, 3000, 1000, "max_epochs"
        best_epoch, best_auprc, wall = 2, 0.0131, 61.0
    else:
        overrides = {"max_epochs": 1, "min_epochs": 1}
        epochs, n_train, n_valid, stop = 1, N_TRAIN, N_VALID, "max_epochs"
        best_epoch, best_auprc, wall = 1, 0.0187, 268.0
    item = make_run_dir(
        root,
        seed=seed,
        unit_id=unit_id,
        budget=budget,
        config_sha256=("c%d" % (1 if budget == "smoke" else 2)).ljust(64, "0"),
        epochs=epochs,
        stop_reason=stop,
        best_epoch=best_epoch,
        best_valid_auprc=best_auprc,
        n_train=n_train,
        n_valid=n_valid,
        peak_rss_mb=2500.0 if budget == "smoke" else 5200.0,
        wall_s=wall,
        recipe=overrides,
    )
    entry = record.run_entry(
        run_dir=item["run_dir"],
        unit_id=item["unit_id"],
        seed=item["seed"],
        budget=budget,
        command=["python", "train.py", "--stage", "pretrain", "--smoke"],
        exit_code=0,
        expected_config_sha256=item["config_sha256"],
        wall_s_outer=wall + 1.0,
        resources={
            "host_peak_rss_mb_pipeline": None,
            "host_peak_rss_mb_sampler": 2550.0 if budget == "smoke" else 5250.0,
            "gpu_peak_mib_device": 3000.0,
            "gpu_peak_mib_process": 2900.0,
            "gpu_baseline_mib": 1050.0,
            "gpu_peak_mib_device_minus_baseline": 1950.0,
            "gpu_after_mib": 1050.0,
            "samples": int(wall),
            "interval_s": 1.0,
        },
    )
    expected = {
        "locked_config_sha256": record.sha256_file(data_dir / "locked_config.json"),
        "structure": dict(LOCKED_CONFIG["structure"]),
        "budget": dict(LOCKED_CONFIG["recipe"]) | {"lr": 0.001},
        "metric": dict(LOCKED_CONFIG["metric"]),
        "loss": "gaussian_nll",
        "optimizer": "adam",
        "prediction_noise": "zero",
        "n_folds": 5,
        "seeds": list(SEEDS),
        "lrs": dict(LOCKED_CONFIG["lrs"]),
        "config_sha256": {unit_id: item["config_sha256"]},
        "config_sha256_by_seed": {},
    }
    record_budget = {
        "kind": "profile",
        "smoke": budget == "smoke",
        "row_subset": True,
        "profile_budgets": [budget],
        "overrides": {},
        "overrides_by_unit": {unit_id: dict(overrides)},
        **(dict(LOCKED_CONFIG["recipe"]) | {"lr": 0.001}),
    }
    return record.build_record(
        mode="profile",
        units=[entry],
        budget=record_budget,
        expectations=expected,
        input_identity=input_identity(data_dir),
        command=["python", "run_pretrain.py", "--mode", "profile"],
        git={"sha": "deadbeef", "branch": "feat/x", "dirty": False},
        seeds_planned=[seed],
        seeds_requested=[seed],
        host={"python": "synthetic", "torch": "2.0.1+cu118"},
    )
