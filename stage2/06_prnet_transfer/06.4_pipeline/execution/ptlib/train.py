"""The one trainer for both stages.

``--stage pretrain`` trains on the frozen Tahoe view (matched ``(plate,
cell_line)`` control input, treated log1p target).  ``--stage finetune`` trains
one CM fold (same-plate matched control input for the training wells, the fold's
fixed consensus c* as the validation/inference input).

Both stages:

* optimise the author Gaussian NLL (``GUSS``; may be negative, never clamped);
* select the checkpoint and early-stop on **validation direction AUPRC@k (max)**;
* step ``ReduceLROnPlateau`` with **max** semantics on that AUPRC;
* save the same checkpoint payload and reload it strictly;
* record per-epoch validation cost (gather / forward / AUPRC split), wall time,
  peak RSS, stop reason, git SHA, data hashes and checkpoint identity.

Noise policy: training noise is drawn per batch from a CPU ``torch.Generator``
seeded with the run seed; validation and inference use zero noise, fixed row
order, ``model.eval()`` — deterministic and repeatable.
"""

from __future__ import annotations

import logging
import math
import platform
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd
import torch

from . import PIPELINE_VERSION, checkpoint as checkpoint_mod
from . import config as config_mod
from . import data as data_mod
from . import evaluate as evaluate_mod
from . import manifest as manifest_mod
from . import metrics as metrics_mod
from . import model as model_mod
from . import modelcard as modelcard_mod
from . import paths
from . import runtime

__all__ = [
    "TrainOptions",
    "TrainResult",
    "seed_everything",
    "resolve_device",
    "resolve_run_config",
    "run_training",
    "epoch_batches",
]

log = logging.getLogger("ptlib.train")

SMOKE_EPOCHS = 2
SMOKE_TAHOE_TRAIN_ROWS = 3000
SMOKE_TAHOE_VALID_ROWS = 1000


@dataclass
class TrainOptions:
    stage: str
    out_dir: Path
    data_dir: Path = paths.DEFAULT_FROZEN_INPUT
    seed: int = 20261021
    device: str = "cuda"
    structure_spec: str | None = None
    init: str = "random"
    init_checkpoint: Path | None = None
    train_drugs: str = "main"
    fold: int | None = None
    epochs: int | None = None
    lr: float | None = None
    batch_size: int | None = None
    min_epochs: int | None = None
    early_stop_patience: int | None = None
    sched_patience: int | None = None
    max_train_rows: int | None = None
    max_valid_rows: int | None = None
    tahoe_cache: Path | None = None
    cache_mode: str = "auto"
    verify_cache: str = "metadata"
    hash_cache: Path | None = None
    table_cache: Path | None = None
    smoke: bool = False
    check: bool = False
    model_card: bool = True
    command: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.out_dir = Path(self.out_dir)
        self.data_dir = Path(self.data_dir)
        if self.stage not in ("pretrain", "finetune"):
            raise SystemExit(f"--stage must be pretrain or finetune, got {self.stage!r}")
        if self.stage == "finetune" and self.fold is None:
            raise SystemExit("--fold is required for --stage finetune")
        if self.init not in ("random", "pretrained"):
            raise SystemExit(f"--init must be random or pretrained, got {self.init!r}")
        if self.init == "pretrained" and not self.init_checkpoint:
            raise SystemExit("--init pretrained requires --init-checkpoint")


@dataclass
class TrainResult:
    run_dir: Path
    checkpoint_path: Path
    best_valid_auprc: float
    best_epoch: int
    epochs_run: int
    stop_reason: str
    history: list[dict[str, Any]]
    manifest: dict[str, Any]

    def summary(self) -> dict[str, Any]:
        return {
            "run_dir": str(self.run_dir),
            "checkpoint": str(self.checkpoint_path),
            "best_valid_auprc": self.best_valid_auprc,
            "best_epoch": self.best_epoch,
            "epochs_run": self.epochs_run,
            "stop_reason": self.stop_reason,
        }


seed_everything = runtime.seed_everything
resolve_device = runtime.resolve_device


def epoch_batches(n_rows: int, batch_size: int) -> list[tuple[int, int]]:
    """Contiguous ``[start, stop)`` slices; a trailing singleton folds into the previous."""
    if n_rows < 2:
        raise SystemExit(f"need at least 2 training rows to train BatchNorm1d, got {n_rows}")
    starts = list(range(0, n_rows, int(batch_size)))
    if len(starts) > 1 and n_rows - starts[-1] == 1:
        starts.pop()
    stops = starts[1:] + [n_rows]
    return list(zip(starts, stops))


def _apply_smoke(options: TrainOptions) -> None:
    if not options.smoke:
        return
    if options.epochs is None:
        options.epochs = SMOKE_EPOCHS
    if options.min_epochs is None:
        options.min_epochs = 1
    if options.stage == "pretrain":
        if options.max_train_rows is None:
            options.max_train_rows = SMOKE_TAHOE_TRAIN_ROWS
        if options.max_valid_rows is None:
            options.max_valid_rows = SMOKE_TAHOE_VALID_ROWS


def resolve_run_config(
    frozen_config: Mapping[str, Any],
    options: TrainOptions,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Resolve the frozen config + explicit overrides into the run config.

    Pure (no data access): the queue uses it to predict a unit's config
    fingerprint for the resume check.
    """
    cfg = {
        "structure": dict(frozen_config["structure"]),
        "recipe": dict(frozen_config["recipe"]),
        "metric": dict(frozen_config["metric"]),
        "loss": frozen_config.get("loss", "gaussian_nll"),
        "optimizer": frozen_config.get("optimizer", "adam"),
        "dose_encoding": dict(frozen_config.get("dose_encoding", {})),
        "prediction_noise": frozen_config.get("prediction_noise", "zero"),
        "train_drugs": options.train_drugs,
        "stage": options.stage,
        "seed": int(options.seed),
        "fold": options.fold,
        "init": options.init,
    }
    struct, struct_notes = config_mod.structure_from_spec(options.structure_spec, cfg["structure"])
    cfg["structure"] = struct
    overrides = {
        "max_epochs": options.epochs,
        "batch_size": options.batch_size,
        "min_epochs": options.min_epochs,
        "early_stop_patience": options.early_stop_patience,
        "sched_patience": options.sched_patience,
    }
    recipe, applied = config_mod.apply_recipe_overrides(cfg["recipe"], overrides)
    cfg["recipe"] = recipe
    # learning rate is a stage-level knob (pretrain LR vs the CM candidates)
    if options.lr is not None:
        cfg["recipe"]["lr"] = float(options.lr)
    else:
        lrs = frozen_config.get("lrs", {}) or {}
        if options.stage == "pretrain":
            cfg["recipe"]["lr"] = float(lrs.get("pretrain", 1e-3))
        else:
            candidates = lrs.get("cm") or [1e-3]
            cfg["recipe"]["lr"] = float(candidates[-1])  # deterministic default; #6 sets it
    config_mod.validate_config(cfg)
    notes = {
        "structure": struct_notes,
        "recipe_overrides": applied,
        "locked_config_mapping": [],
    }
    return cfg, notes


def _resolve_tahoe_cache(
    frozen: data_mod.FrozenInput,
    options: TrainOptions,
) -> tuple[Any, dict[str, Any]]:
    """Resolve the derived memmap cache per ``--cache-mode``."""
    cache_path = Path(options.tahoe_cache) if options.tahoe_cache else (
        paths.DEFAULT_OUTPUT / "cache" / data_mod_cache_name()
    )
    source_sha = frozen.artifact_hashes.get("tahoe_transfer.h5ad")
    mode = options.cache_mode
    exists = cache_path.exists()
    if mode in ("off", None):
        return None, {"mode": "off", "path": str(cache_path)}
    if mode == "require" and not exists:
        raise SystemExit(
            f"--cache-mode require but the derived Tahoe cache is missing at {cache_path}; "
            "build it with --cache-mode build (one-off, ~5 GB)."
        )
    if mode == "auto" and not exists:
        return None, {"mode": "absent", "path": str(cache_path)}
    if mode == "build":
        record = cache_mod_build(frozen, cache_path, source_sha, log)
        cache = cache_mod_open(cache_path, source_sha, frozen.tahoe_h5ad, options.verify_cache, log)
        return cache, {"mode": "built", "path": str(cache_path), **record}
    cache = cache_mod_open(cache_path, source_sha, frozen.tahoe_h5ad, options.verify_cache, log)
    return cache, {
        "mode": "reused",
        "path": str(cache_path),
        "verify": options.verify_cache,
        "sample_max_abs_diff": cache.sample_max_abs_diff,
        "layer_sha256": cache.record.get("layer_sha256"),
        "file_sha256": cache.record.get("file_sha256"),
        "source_sha256": cache.record.get("source_sha256"),
        "build_seconds": cache.record.get("build_seconds"),
    }


def data_mod_cache_name() -> str:
    from . import cache as cache_mod

    return cache_mod.DEFAULT_CACHE_NAME


def cache_mod_build(frozen, cache_path, source_sha, logger):
    from . import cache as cache_mod

    return cache_mod.build_cache(
        frozen.tahoe_h5ad, cache_path, source_sha, log=logger
    )


def cache_mod_open(cache_path, source_sha, tahoe_h5ad, verify, logger):
    from . import cache as cache_mod

    return cache_mod.open_cache(
        cache_path, source_sha, source_h5ad=tahoe_h5ad, verify=verify
    )


def _build_dataset(stage: str, frozen, options, logger):
    structures = frozen.structures()
    info: dict[str, Any] = {}
    if stage == "pretrain":
        cache, cache_info = _resolve_tahoe_cache(frozen, options)
        dataset = data_mod.TahoeDataset(
            frozen=frozen,
            structures=structures,
            max_train_rows=options.max_train_rows,
            max_valid_rows=options.max_valid_rows,
            seed=options.seed,
            cache=cache,
            full=(options.max_train_rows is None and cache is None),
        ).prepare()
        info["cache"] = cache_info
        info["counts"] = dataset.counts
        info["structures"] = structures.sources
        info["rows"] = {"train": int(dataset.rows_train.size), "valid": int(dataset.rows_valid.size)}
        return dataset, info
    folds = data_mod.load_folds(frozen)
    dataset = data_mod.CMDataset(
        frozen=frozen,
        fold=int(options.fold),
        fold_table=folds,
        structures=structures,
        train_drugs=options.train_drugs,
    ).prepare()
    info["counts"] = dataset.counts
    info["structures"] = structures.sources
    info["folds"] = folds.report
    info["consensus"] = dataset.consensus_check
    return dataset, info


def _build_model_for_run(frozen, cfg, options, device, logger):
    struct = cfg["structure"]
    model = model_mod.build_model(struct)
    init_record: dict[str, Any] = {"mode": options.init}
    if options.init == "random":
        model_mod.random_init(model, options.seed)
        init_record["source_path"] = None
        init_record["checkpoint_sha256"] = None
    else:
        path = Path(options.init_checkpoint)
        payload = checkpoint_mod.read_checkpoint(path)
        checkpoint_mod.validate_checkpoint(payload, int(struct["n_genes"]), struct)
        checkpoint_mod.load_into_model(model, payload)
        init_record["source_path"] = str(path)
        init_record["checkpoint_sha256"] = checkpoint_mod.checkpoint_sha256(path)
        init_record["checkpoint_epoch"] = int(payload.get("epoch", -1))
        init_record["checkpoint_best_valid_auprc"] = float(payload.get("best_valid_auprc", math.nan))
        init_record["checkpoint_seed"] = int(payload.get("seed", -1))
        journal_init(logger, init_record)
    model.to(device)
    model_mod.check_architecture(model, struct)
    init_record["parameter_count"] = model_mod.parameter_count(model)
    return model, init_record


def journal_init(logger, record: Mapping[str, Any]) -> None:
    logger.info(
        "warm start: loaded %s (sha256 %s…, epoch %s, run seed %s)",
        record.get("source_path"),
        str(record.get("checkpoint_sha256"))[:12],
        record.get("checkpoint_epoch"),
        record.get("checkpoint_seed"),
    )


def _tensor(batch: np.ndarray, device: torch.device) -> torch.Tensor:
    return torch.from_numpy(np.ascontiguousarray(batch)).to(device)


def run_training(options: TrainOptions) -> TrainResult:
    """Run one training unit; returns its result and writes all run artifacts."""
    _apply_smoke(options)
    if options.smoke:
        log.info("SMOKE mode: budgets are pipeline checks, not scientific results")
    seed_everything(options.seed)
    started_wall = time.perf_counter()
    started_utc = manifest_mod.now_utc()
    options.out_dir.mkdir(parents=True, exist_ok=True)

    frozen = data_mod.open_frozen_input(
        options.data_dir, hash_cache=options.hash_cache, table_cache=options.table_cache
    )
    cfg, cfg_notes = resolve_run_config(frozen.config, options)
    cfg_notes["locked_config_mapping"] = frozen.config_notes.get("mapping", [])
    cfg_notes["locked_config_defaulted"] = frozen.config_notes.get("defaulted", [])
    device, device_name = resolve_device(options.device)

    dataset, data_info = _build_dataset(options.stage, frozen, options, log)
    model, init_record = _build_model_for_run(frozen, cfg, options, device, log)

    recipe = cfg["recipe"]
    k = int(cfg["metric"]["k"])
    n_dim = int(cfg["structure"]["n_dim"])
    lr = float(recipe["lr"])
    optimizer = torch.optim.Adam(
        model.parameters(), lr=lr, weight_decay=float(recipe["weight_decay"])
    )
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode="max",
        factor=float(recipe["sched_factor"]),
        patience=int(recipe["sched_patience"]),
        min_lr=float(recipe["min_lr"]),
    )

    if options.stage == "pretrain":
        train_rows = dataset.rows_train
        valid_rows = dataset.rows_valid
    else:
        train_rows = dataset.train_rows
        valid_rows = dataset.valid_rows
    valid_conditions = dataset.condition_ids(valid_rows)
    truth = None
    if options.stage == "finetune":
        truth = data_mod.CMControlPool(frozen)
        matched = truth.mean_for_rows(valid_rows)
        valid_truth_baseline = matched
    else:
        valid_truth_baseline = None  # matched control == the model input for Tahoe

    logger = log
    epoch_generator = torch.Generator(device="cpu").manual_seed(int(options.seed))
    max_epochs = int(recipe["max_epochs"])
    min_epochs = int(recipe["min_epochs"])
    patience = int(recipe["early_stop_patience"])
    checkpoint_path = options.out_dir / checkpoint_mod.CHECKPOINT_RELPATH
    history: list[dict[str, Any]] = []
    best_auprc = -math.inf
    best_epoch = 0
    stale = 0
    stop_reason = "max_epochs"

    for epoch in range(max_epochs):
        epoch_started = time.perf_counter()
        model.train()
        order = torch.randperm(int(train_rows.size), generator=epoch_generator).numpy()
        total_loss, seen = 0.0, 0
        train_started = time.perf_counter()
        for start, stop in epoch_batches(int(train_rows.size), int(recipe["batch_size"])):
            idx = order[start:stop]
            rows = train_rows[idx]
            control, label, target = dataset.batch(rows)
            noise = torch.randn(
                int(len(idx)), n_dim, generator=epoch_generator
            ).to(device)
            loss = model_mod.gaussian_nll(
                model,
                _tensor(control, device),
                _tensor(label, device),
                noise,
                _tensor(target, device),
            )
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()
            total_loss += float(loss.item()) * int(idx.size)
            seen += int(idx.size)
        train_loss = total_loss / max(seen, 1)
        train_seconds = time.perf_counter() - train_started

        model.eval()
        scored = evaluate_mod.score_validation(
            model=model,
            dataset=dataset,
            rows=valid_rows,
            condition_ids=valid_conditions,
            k=k,
            batch_size=int(recipe["batch_size"]),
            device=device,
            n_dim=n_dim,
            truth_baseline=valid_truth_baseline,
        )
        scheduler.step(scored["auprc"])
        lr_now = float(optimizer.param_groups[0]["lr"])
        seconds = time.perf_counter() - epoch_started

        improved = bool(
            np.isfinite(scored["auprc"]) and scored["auprc"] > best_auprc
        )
        if improved:
            best_auprc = float(scored["auprc"])
            best_epoch = epoch + 1
            stale = 0
            checkpoint_mod.save_checkpoint(
                checkpoint_path,
                model,
                config={"structure": cfg["structure"], "recipe": cfg["recipe"],
                        "metric": cfg["metric"], "stage": options.stage,
                        "train_drugs": options.train_drugs, "seed": options.seed,
                        "init": options.init},
                seed=int(options.seed),
                n_genes=int(cfg["structure"]["n_genes"]),
                best_valid_auprc=best_auprc,
                epoch=epoch + 1,
            )
        else:
            stale += 1

        entry = {
            "epoch": epoch + 1,
            "train_loss": train_loss,
            "valid_loss": scored["valid_nll"],
            "valid_auprc": float(scored["auprc"]),
            "valid_auprc_up": float(scored["auprc_up"]),
            "valid_auprc_down": float(scored["auprc_down"]),
            "valid_rmse": float(scored["rmse"]),
            "valid_pcc": float(scored["pcc"]),
            "valid_spearman": float(scored["spearman"]),
            "lr": lr_now,
            "train_seconds": train_seconds,
            "valid_seconds": scored["timings"]["total_s"],
            "valid_gather_seconds": scored["timings"]["gather_s"],
            "valid_forward_seconds": scored["timings"]["forward_s"],
            "valid_metrics_seconds": scored["timings"]["metrics_s"],
            "n_train": int(train_rows.size),
            "n_valid": int(valid_rows.size),
            "n_conditions_valid": int(scored["n_conditions"]),
            "seconds": seconds,
        }
        history.append(entry)
        logger.info(
            "epoch %d/%d | train_nll=%.6f | valid_nll=%.6f | valid_AUPRC@%d=%.6f | "
            "lr=%.2e | train %.1fs | valid %.1fs (gather %.1fs/forward %.1fs/ap %.1fs) | "
            "best %.6f @%d | stale %d/%d%s",
            epoch + 1, max_epochs, train_loss, scored["valid_nll"], k, scored["auprc"],
            lr_now, train_seconds, scored["timings"]["total_s"],
            scored["timings"]["gather_s"], scored["timings"]["forward_s"],
            scored["timings"]["metrics_s"], best_auprc, best_epoch, stale, patience,
            " (smoke: early stop disabled)" if options.smoke else "",
        )
        if not math.isfinite(train_loss) or not math.isfinite(scored["valid_nll"]):
            logger.warning(
                "epoch %d produced a non-finite NLL (not clamped): train=%r valid=%r",
                epoch + 1, train_loss, scored["valid_nll"],
            )

        if not options.smoke and (epoch + 1) >= min_epochs and stale > patience:
            stop_reason = f"early_stop_no_auprc_improvement_{patience}"
            logger.info("early stop at epoch %d: no AUPRC improvement for %d epochs", epoch + 1, stale)
            break

    if best_epoch == 0:  # no finite improvement: keep a usable last-epoch checkpoint
        logger.warning("no epoch improved validation AUPRC; saving the final weights")
        checkpoint_mod.save_checkpoint(
            checkpoint_path, model, config={"structure": cfg["structure"], "recipe": cfg["recipe"],
                                            "metric": cfg["metric"], "stage": options.stage,
                                            "train_drugs": options.train_drugs,
                                            "seed": options.seed, "init": options.init},
            seed=int(options.seed), n_genes=int(cfg["structure"]["n_genes"]),
            best_valid_auprc=math.nan, epoch=len(history),
        )
        stop_reason = "no_finite_auprc"
        best_epoch = len(history)

    wall_s = time.perf_counter() - started_wall
    finished_utc = manifest_mod.now_utc()

    # Write the run record.
    config_payload = dict(cfg)
    config_payload["config_sha256"] = config_mod.config_fingerprint(cfg)
    config_payload["run"] = {
        "stage": options.stage,
        "fold": options.fold,
        "seed": options.seed,
        "device": device_name,
        "command": list(options.command),
        "smoke": bool(options.smoke),
        "structure_notes": cfg_notes.get("structure"),
        "recipe_overrides": cfg_notes.get("recipe_overrides"),
        "locked_config_mapping": cfg_notes.get("locked_config_mapping"),
        "locked_config_defaulted": cfg_notes.get("locked_config_defaulted"),
        "tahoe_cache": data_info.get("cache"),
        "data_counts": data_info.get("counts"),
        "structures": data_info.get("structures"),
        "consensus": data_info.get("consensus"),
        "folds": data_info.get("folds"),
        "started_utc": started_utc,
        "finished_utc": finished_utc,
        "wall_s": wall_s,
    }
    manifest_mod.write_json(options.out_dir / "config.json", config_payload)
    fit_history = {
        "run": {
            "stage": options.stage,
            "fold": options.fold,
            "seed": options.seed,
            "started_utc": started_utc,
            "finished_utc": finished_utc,
            "wall_s": wall_s,
            "peak_rss_mb": manifest_mod.peak_rss_mb(),
            "stop_reason": stop_reason,
            "best_valid_auprc": best_auprc if math.isfinite(best_auprc) else None,
            "best_epoch": best_epoch,
            "metric": cfg["metric"],
            "smoke": bool(options.smoke),
        },
        "history": history,
    }
    manifest_mod.write_json(options.out_dir / "fit_history.json", fit_history)
    pd.DataFrame(history).to_csv(
        options.out_dir / "summary_table.csv", index=False, float_format="%.10g"
    )

    checkpoint_sha = checkpoint_mod.checkpoint_sha256(checkpoint_path)
    # Reload strictly and confirm the state round-trips (a warm-start contract check).
    payload = checkpoint_mod.read_checkpoint(checkpoint_path)
    reload_model = model_mod.build_model(payload["config"]["structure"])
    checkpoint_mod.load_into_model(reload_model, payload)
    checkpoint_hash_match = all(
        torch.equal(reload_model.state_dict()[key].cpu(), value)
        for key, value in payload["state_dict"].items()
    )

    timings = {
        "started_utc": started_utc,
        "finished_utc": finished_utc,
        "wall_s": wall_s,
        "train_loop_s": float(sum(entry["train_seconds"] for entry in history)),
        "per_epoch_valid_s": [float(entry["valid_seconds"]) for entry in history],
        "per_epoch_valid_gather_s": [float(entry["valid_gather_seconds"]) for entry in history],
        "per_epoch_valid_forward_s": [float(entry["valid_forward_seconds"]) for entry in history],
        "per_epoch_valid_metrics_s": [float(entry["valid_metrics_seconds"]) for entry in history],
        "per_epoch_valid_total_s": float(sum(entry["valid_seconds"] for entry in history)),
        "per_epoch_train_s": [float(entry["train_seconds"]) for entry in history],
    }
    manifest = manifest_mod.build_manifest(
        "train",
        stage=options.stage,
        status="completed",
        pipeline_version=PIPELINE_VERSION,
        command=list(options.command),
        git=manifest_mod.git_info(paths.ROOT),
        data=frozen.record(),
        locked_config={
            "path": str(frozen.data_dir / "locked_config.json"),
            "sha256": manifest_mod.sha256_file(frozen.data_dir / "locked_config.json"),
            "mapping": frozen.config_notes.get("mapping", []),
            "defaulted": frozen.config_notes.get("defaulted", []),
        },
        config=cfg,
        config_sha256=config_mod.config_fingerprint(cfg),
        structure=cfg["structure"],
        seed=int(options.seed),
        fold=options.fold,
        device=manifest_mod.device_record(options.device, str(device), device_name),
        init=init_record,
        checkpoint={
            "path": checkpoint_mod.CHECKPOINT_RELPATH,
            "sha256": checkpoint_sha,
            "reload_identical": bool(checkpoint_hash_match),
            "best_valid_auprc": best_auprc if math.isfinite(best_auprc) else None,
            "best_epoch": best_epoch,
        },
        noise_policy={
            "train": "per-batch torch.randn from the seeded CPU generator",
            "validation": "zero noise, model.eval(), fixed row order",
            "inference": "zero noise, model.eval(), fixed row order",
        },
        timings=timings,
        peak_rss_mb=manifest_mod.peak_rss_mb(),
        epochs_run=len(history),
        stop_reason=stop_reason,
        best_valid_auprc=best_auprc if math.isfinite(best_auprc) else None,
        best_epoch=best_epoch,
        metric={"name": "direction_auprc", "k": k, "plateau_mode": "max"},
        smoke=bool(options.smoke),
        data_counts=data_info.get("counts"),
        tahoe_cache=data_info.get("cache"),
        warnings=[],
    )
    manifest_mod.write_json(options.out_dir / "manifest.json", manifest)

    result = TrainResult(
        run_dir=options.out_dir,
        checkpoint_path=checkpoint_path,
        best_valid_auprc=best_auprc,
        best_epoch=best_epoch,
        epochs_run=len(history),
        stop_reason=stop_reason,
        history=history,
        manifest=manifest,
    )
    if options.model_card:
        modelcard_mod.write_model_card(options.out_dir)
    log.info(
        "finished %s: epochs=%d stop=%s best_valid_AUPRC@%d=%.6f wall=%.1fs peak_rss=%.0f MiB",
        options.stage, len(history), stop_reason, k,
        best_auprc if math.isfinite(best_auprc) else float("nan"),
        wall_s, manifest_mod.peak_rss_mb(),
    )
    return result


# --------------------------------------------------------------------------- #
# CLI argument surface (shared by ``train.py`` and ``run_queue.py``)
# --------------------------------------------------------------------------- #
def _cache_default() -> Path:
    return paths.DEFAULT_OUTPUT / "cache" / data_mod_cache_name()


def build_parser() -> "argparse.ArgumentParser":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "PRnet formal-transfer trainer (ticket #4). One invocation = one unit: "
            "--stage pretrain (Tahoe, matched controls) or --stage finetune (one CM "
            "fold, fixed consensus c* at validation). Selection/early stop use "
            "validation direction AUPRC@k (max)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--stage", choices=("pretrain", "finetune"), required=True)
    parser.add_argument("--data-dir", type=Path, default=paths.DEFAULT_FROZEN_INPUT,
                        help="frozen 06.3 input directory")
    parser.add_argument("--out", type=Path, required=True, help="run output directory")
    parser.add_argument("--fold", type=int, default=None,
                        help="CM fold (required for --stage finetune)")
    parser.add_argument("--seed", type=int, default=20261021)
    parser.add_argument("--device", default="cuda", choices=("cpu", "cuda", "cuda:0"),
                        help="torch device; cuda is a hard failure when unavailable")
    parser.add_argument("--structure", dest="structure_spec", default=None,
                        help="structure preset (hidden128|latent128|hidden256) or a JSON file")
    parser.add_argument("--init", choices=("random", "pretrained"), default="random")
    parser.add_argument("--init-checkpoint", type=Path, default=None,
                        help="strict warm-start checkpoint (required for --init pretrained)")
    parser.add_argument("--train-drugs", choices=("main", "expanded"), default="main",
                        help="#11 extension: 'expanded' adds eligibility-excluded QC-pass wells")
    parser.add_argument("--epochs", type=int, default=None)
    parser.add_argument("--lr", type=float, default=None)
    parser.add_argument("--batch-size", type=int, default=None)
    parser.add_argument("--min-epochs", type=int, default=None)
    parser.add_argument("--early-stop-patience", type=int, default=None)
    parser.add_argument("--sched-patience", type=int, default=None)
    parser.add_argument("--max-train-rows", type=int, default=None,
                        help="deterministic Tahoe train subset (smoke / bring-up)")
    parser.add_argument("--max-valid-rows", type=int, default=None)
    parser.add_argument("--tahoe-cache", type=Path, default=None,
                        help=f"derived log1p memmap (default {_cache_default()})")
    parser.add_argument("--cache-mode", choices=("auto", "require", "build", "off"),
                        default="auto",
                        help="derived Tahoe cache policy (build once, then it is reused)")
    parser.add_argument("--verify-cache", choices=("metadata", "deep"), default="metadata",
                        help="deep also recomputes the cache full-array checksum")
    parser.add_argument("--hash-cache", type=Path, default=None,
                        help="memoise frozen-artifact sha256 by (size, mtime)")
    parser.add_argument("--table-cache", type=Path, default=None,
                        help="derived CSV table cache (pinned env has no parquet engine)")
    parser.add_argument("--smoke", action="store_true",
                        help="2 epochs / Tahoe subset; pipeline check, not a result")
    parser.add_argument("--check", action="store_true",
                        help="verify data/config/architecture and print the resolved run, then exit")
    parser.add_argument("--no-model-card", action="store_true", help="skip model_card.md")
    return parser


def options_from_args(argv: list[str]) -> TrainOptions:
    parser = build_parser()
    ns = parser.parse_args(argv)
    return TrainOptions(
        stage=ns.stage,
        out_dir=ns.out,
        data_dir=ns.data_dir,
        seed=ns.seed,
        device=ns.device,
        structure_spec=ns.structure_spec,
        init=ns.init,
        init_checkpoint=ns.init_checkpoint,
        train_drugs=ns.train_drugs,
        fold=ns.fold,
        epochs=ns.epochs,
        lr=ns.lr,
        batch_size=ns.batch_size,
        min_epochs=ns.min_epochs,
        early_stop_patience=ns.early_stop_patience,
        sched_patience=ns.sched_patience,
        max_train_rows=ns.max_train_rows,
        max_valid_rows=ns.max_valid_rows,
        tahoe_cache=ns.tahoe_cache,
        cache_mode=ns.cache_mode,
        verify_cache=ns.verify_cache,
        hash_cache=ns.hash_cache,
        table_cache=ns.table_cache,
        smoke=ns.smoke,
        check=ns.check,
        model_card=not ns.no_model_card,
        command=["train.py", *argv],
    )


def described_run(options: TrainOptions) -> dict[str, Any]:
    """Resolve and describe a run without reading any expression matrix."""
    _apply_smoke(options)
    frozen = data_mod.open_frozen_input(
        options.data_dir, hash_cache=options.hash_cache, table_cache=options.table_cache
    )
    cfg, notes = resolve_run_config(frozen.config, options)
    return {
        "stage": options.stage,
        "fold": options.fold,
        "seed": options.seed,
        "init": options.init,
        "structure": cfg["structure"],
        "recipe": cfg["recipe"],
        "metric": cfg["metric"],
        "config_sha256": config_mod.config_fingerprint(cfg),
        "frozen_input": frozen.record(),
        "notes": notes,
    }


def check_run(options: TrainOptions) -> dict[str, Any]:
    """Verify frozen input, resolved config, dataset counts and architecture; no training."""
    _apply_smoke(options)
    frozen = data_mod.open_frozen_input(
        options.data_dir, hash_cache=options.hash_cache, table_cache=options.table_cache
    )
    cfg, notes = resolve_run_config(frozen.config, options)
    structures = frozen.structures()
    info: dict[str, Any] = {"stage": options.stage, "fold": options.fold}
    if options.stage == "pretrain":
        cache = None
        cache_info: dict[str, Any] = {"mode": "off"}
        if options.cache_mode != "off":
            cache, cache_info = _resolve_tahoe_cache(frozen, options)
        dataset = data_mod.TahoeDataset(
            frozen=frozen,
            structures=structures,
            max_train_rows=options.max_train_rows,
            max_valid_rows=options.max_valid_rows,
            seed=options.seed,
            cache=cache,
            load_targets=False,
        ).prepare()
        info["counts"] = dataset.counts
        info["rows"] = {"train": int(dataset.rows_train.size), "valid": int(dataset.rows_valid.size)}
        info["cache"] = cache_info
    else:
        folds = data_mod.load_folds(frozen)
        dataset = data_mod.CMDataset(
            frozen=frozen,
            fold=int(options.fold),
            fold_table=folds,
            structures=structures,
            train_drugs=options.train_drugs,
        ).prepare()
        info["counts"] = dataset.counts
        info["consensus"] = dataset.consensus_check
        info["folds"] = folds.report
    model = model_mod.build_model(cfg["structure"])
    model_mod.check_architecture(model, cfg["structure"])
    info["model"] = {
        "structure": cfg["structure"],
        "parameter_count": model_mod.parameter_count(model),
    }
    info["recipe"] = cfg["recipe"]
    info["metric"] = cfg["metric"]
    info["config_sha256"] = config_mod.config_fingerprint(cfg)
    info["frozen_input"] = frozen.record()
    info["notes"] = notes
    info["smoke"] = bool(options.smoke)
    info["ok"] = True
    return info
