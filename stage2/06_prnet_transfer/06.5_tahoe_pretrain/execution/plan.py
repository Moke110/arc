"""Resolve the ticket-#5 run plan through the 06.4 pipeline's own config code.

Why this module exists: the wrapper must invoke ``train.py`` with the **locked
config** and must not re-derive the recipe itself.  The 06.4 pipeline already
normalises the frozen ``locked_config.json`` and resolves a run config
(``ptlib.config`` + ``ptlib.train.resolve_run_config``); this module calls those
exact functions to compute the config digest each planned run must record, and
to state the expected budget before a single step is run.

Import contract: ``run_pretrain.py`` puts
``<repo>/stage2/06_prnet_transfer/06.4_pipeline/execution`` on ``sys.path``
**before** importing this module; the import below then binds the pinned GPU
environment's ``ptlib``.  A wrong interpreter (no torch) fails at import — which
is the intended fail-closed behaviour, because the training runs need it anyway.

The only budget this module *adds* is the profile budget (2-epoch `--smoke`
subset, or one full-width epoch).  A formal unit applies **no** overrides: the
resolved recipe is the frozen recipe and the learning rate is ``lrs.pretrain``.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

_EXECUTION = Path(__file__).resolve().parent
_PIPELINE_EXECUTION = _EXECUTION.parent.parent / "06.4_pipeline" / "execution"
if str(_PIPELINE_EXECUTION) not in sys.path:
    sys.path.insert(0, str(_PIPELINE_EXECUTION))

from ptlib import config as config_mod
from ptlib import paths as ptlib_paths
from ptlib import train as trainlib

__all__ = [
    "PROFILE_BUDGETS",
    "load_locked",
    "expectation_block",
    "plan_unit",
    "plan_units",
    "pipeline_paths",
]

#: The profile budgets the wrapper can run.  Both are pipeline/throughput
#: checks; neither is training evidence.
PROFILE_BUDGETS = ("smoke", "one-epoch")


def pipeline_paths() -> dict[str, Path]:
    """The 06.4 pipeline locations this substep reuses (never copies)."""
    return {
        "root": Path(ptlib_paths.ROOT),
        "execution": Path(ptlib_paths.EXECUTION),
        "pipeline": Path(ptlib_paths.PIPELINE),
        "frozen_input": Path(ptlib_paths.DEFAULT_FROZEN_INPUT),
        "output": Path(ptlib_paths.DEFAULT_OUTPUT),
    }


def load_locked(data_dir: Path | str) -> tuple[dict[str, Any], Mapping[str, Any], dict[str, Any]]:
    """Read + normalise the frozen locked config exactly as the pipeline does.

    Returns ``(raw, frozen_config, notes)``.  ``raw`` is kept for the independent
    ticket-expectation check; ``frozen_config`` is what ``resolve_run_config``
    consumes.
    """
    path = Path(data_dir) / "locked_config.json"
    raw = config_mod.load_json(path)
    frozen_config, notes = config_mod.normalize_locked_config(raw)
    return raw, frozen_config, notes


def _apply_profile_budget(options: "trainlib.TrainOptions", budget: str) -> None:
    if budget == "formal":
        return
    if budget == "smoke":
        # Use the pipeline's own smoke mutator so the resolved config is exactly
        # what ``train.py --smoke`` will compute.  Fail closed if it disappears.
        apply_smoke = getattr(trainlib, "_apply_smoke", None)
        if apply_smoke is None:
            raise SystemExit(
                "06.4 ptlib.train no longer exposes _apply_smoke; the 06.5 profile plan "
                "cannot guarantee the smoke budget. Update 06.5/execution/plan.py."
            )
        options.smoke = True
        apply_smoke(options)
        return
    if budget == "one-epoch":
        options.epochs = 1
        options.min_epochs = 1
        return
    raise SystemExit(f"unknown profile budget {budget!r} (known: {PROFILE_BUDGETS})")


def plan_unit(
    frozen_config: Mapping[str, Any],
    seed: int,
    budget: str,
    *,
    data_dir: Path | str,
    device: str = "cuda",
) -> dict[str, Any]:
    """Resolve one unit without touching data: budget, config digest, overrides."""
    options = trainlib.TrainOptions(
        stage="pretrain",
        out_dir=Path("/nonexistent/06.5-plan"),
        data_dir=Path(data_dir),
        seed=int(seed),
        device=device,
    )
    _apply_profile_budget(options, budget)
    resolved, notes = trainlib.resolve_run_config(frozen_config, options)
    fingerprint = config_mod.config_fingerprint(resolved)
    if budget == "formal":
        if notes.get("recipe_overrides"):
            raise SystemExit(
                "a formal unit resolved recipe overrides "
                f"{notes['recipe_overrides']!r}; the formal budget must be the locked recipe"
            )
        if int(resolved["recipe"]["max_epochs"]) != int(frozen_config["recipe"]["max_epochs"]):
            raise SystemExit("a formal unit did not resolve the locked max_epochs")
        if float(resolved["recipe"]["lr"]) != float((frozen_config.get("lrs") or {})["pretrain"]):
            raise SystemExit("a formal unit did not resolve lrs.pretrain from the locked config")
    return {
        "seed": int(seed),
        "budget": budget,
        "config": resolved,
        "config_sha256": fingerprint,
        "overrides": dict(notes.get("recipe_overrides") or {}),
        "budget_config": {
            **{key: resolved["recipe"][key] for key in (
                "batch_size", "max_epochs", "min_epochs", "early_stop_patience",
                "sched_factor", "sched_patience", "min_lr", "weight_decay", "lr",
            )},
            "metric": dict(resolved["metric"]),
        },
    }


def plan_units(
    frozen_config: Mapping[str, Any],
    seeds: Sequence[int],
    budgets: Sequence[str],
    *,
    data_dir: Path | str,
    device: str = "cuda",
) -> list[dict[str, Any]]:
    """One planned unit per (seed, budget), in run order (seed-major)."""
    units: list[dict[str, Any]] = []
    for seed in seeds:
        for budget in budgets:
            units.append(
                plan_unit(frozen_config, int(seed), budget, data_dir=data_dir, device=device)
            )
    return units


def expectation_block(
    frozen_config: Mapping[str, Any],
    locked_config_sha256: str,
    planned: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """The record's ``expected`` block: the locked reference + every config digest.

    Digests are keyed by ``unit_id`` (a profile record can run two budgets of the
    same seed) and, for formal units, also by seed.
    """
    recipe = dict(frozen_config["recipe"])
    recipe["lr"] = float((frozen_config.get("lrs") or {})["pretrain"])
    fingerprints: dict[str, str] = {}
    by_seed: dict[str, str] = {}
    for unit in planned:
        fingerprints[str(unit["unit_id"])] = unit["config_sha256"]
        if unit["budget"] == "formal":
            by_seed[str(unit["seed"])] = unit["config_sha256"]
    return {
        "locked_config_sha256": locked_config_sha256,
        "structure": dict(frozen_config["structure"]),
        "budget": recipe,
        "metric": dict(frozen_config["metric"]),
        "loss": frozen_config.get("loss"),
        "optimizer": frozen_config.get("optimizer"),
        "prediction_noise": frozen_config.get("prediction_noise"),
        "n_folds": int(frozen_config.get("n_folds", 0)),
        "seeds": [int(seed) for seed in frozen_config.get("seeds", [])],
        "lrs": dict(frozen_config.get("lrs") or {}),
        "config_sha256": fingerprints,
        "config_sha256_by_seed": by_seed,
    }
