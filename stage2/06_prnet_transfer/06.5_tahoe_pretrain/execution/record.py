"""Record writer, reader and checker for the three seeded Tahoe pretrains (#5).

The training itself lives in ``06.4_pipeline`` (``train.py --stage pretrain``);
this module owns the *evidence* around it: one machine-readable record that
proves what was trained, on which frozen input, with which resolved config, for
how long, at what cost, and how the run stopped.  It is deliberately
**stdlib-only** so the checker runs in any environment and the tests need no
GPU, no torch and no frozen data.

Two files are written per record:

* ``<record>.json`` — the record itself (schema below); the authoritative
  artifact.  It embeds the input identity (frozen artifact sha256s, gene axis,
  matched-control groups, derived-cache identity), the locked-config digest
  contract, the resolved per-run config digest, the full validation curve and
  the measured resources.
* ``<record>_curve.csv`` — one row per (unit, epoch) so a reader can plot the
  validation AUPRC curves without parsing JSON.

Schema ``06.5_tahoe_pretrain/record/1``::

    {
      "schema", "kind", "mode" ("profile"|"formal"), "complete",
      "created_utc", "git", "command", "budget", "selection", "expected",
      "seeds_planned", "seeds_requested", "seeds_completed", "seeds_failed",
      "no_seed_cherry_picking",
      "input": {data_dir, manifest_sha256, locked_config_sha256, contract_sha256,
                artifacts{}, gene_axis{}, control_groups{}, tahoe_split{},
                cache{}, table_cache{}},
      "units": [ unit_entry ... ],
      "summary": {...}
    }

Each unit entry is one training invocation: ``unit_id``, ``seed``, ``budget``
(formal / smoke / one-epoch), the exact ``command``, ``exit_code``,
``status``, ``stop_reason``, ``epochs_run``, ``best_epoch``,
``best_valid_auprc``, ``config``/``config_sha256``/``expected_config_sha256``,
``checkpoint`` (path + sha256 + strict-reload flag), ``model_card``,
``validation_curve``, ``timings``, ``resources`` and ``failures``.

The checker never trusts the record: it re-reads each run directory's own
``manifest.json`` / ``fit_history.json`` when present, recomputes checkpoint
hashes, compares the recorded budget against the frozen locked config and the
ticket's locked recipe, and enforces that a formal record contains exactly the
three frozen seeds — no seed is dropped after the fact.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

__all__ = [
    "SCHEMA",
    "KIND",
    "CURVE_COLUMNS",
    "STOP_REASONS",
    "TICKET_EXPECTATIONS",
    "RecordError",
    "now_utc",
    "sha256_file",
    "read_json",
    "write_json",
    "find_repo_root",
    "git_info",
    "normalize_locked_config",
    "ticket_expectation_errors",
    "frozen_input_identity",
    "run_entry",
    "build_record",
    "write_record",
    "read_record",
    "check_record",
    "resource_report",
]

SCHEMA = "06.5_tahoe_pretrain/record/1"
KIND = "tahoe_pretrain_seeds"

#: The curve columns copied out of the 06.4 ``fit_history.json`` per epoch.
CURVE_COLUMNS = (
    "epoch",
    "train_loss",
    "valid_loss",
    "valid_auprc",
    "valid_auprc_up",
    "valid_auprc_down",
    "valid_rmse",
    "valid_pcc",
    "valid_spearman",
    "lr",
    "train_seconds",
    "valid_seconds",
)

#: The only stop reasons the 06.4 trainer can emit.
STOP_REASONS = (
    "max_epochs",
    "early_stop_no_auprc_improvement_20",
    "no_finite_auprc",
)

#: The ticket's locked recipe.  ``06.3_input/output/locked_config.json`` is
#: frozen evidence, but the wrapper promises *this* budget: if the frozen file
#: ever disagrees, the run is refused rather than silently re-scoped.  This is
#: an independent acceptance table, not a second copy of pipeline logic.
TICKET_EXPECTATIONS: dict[str, Any] = {
    "structure": {
        "n_genes": 18164,
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
        "min_lr": 1e-8,
        "weight_decay": 1e-8,
    },
    "lrs": {"pretrain": 1e-3},
    "metric": {"name": "direction_auprc", "k": 200, "plateau_mode": "max"},
    "loss": "gaussian_nll",
    "optimizer": "adam",
    "prediction_noise": "zero",
    "n_folds": 5,
    "seeds": [20261021, 20261022, 20261023],
}

#: Floating-point tolerance for comparing frozen recipe values.
_REL_TOL = 1e-12


class RecordError(RuntimeError):
    """A record, a run directory or the frozen input cannot be trusted."""


# --------------------------------------------------------------------------- #
# small stdlib helpers
# --------------------------------------------------------------------------- #
def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: Path | str, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(chunk_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path | str) -> dict[str, Any]:
    p = Path(path)
    try:
        payload = json.loads(p.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise RecordError(f"JSON file not found: {p}") from exc
    except json.JSONDecodeError as exc:
        raise RecordError(f"JSON file is not valid JSON: {p} ({exc})") from exc
    if not isinstance(payload, dict):
        raise RecordError(f"JSON file must hold an object: {p}")
    return payload


def write_json(path: Path | str, payload: Mapping[str, Any]) -> Path:
    """Write JSON atomically (temp file + rename) so a killed run leaves no half file."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    os.replace(tmp, p)
    return p


def find_repo_root(start: Path | str | None = None) -> Path:
    """ADR 0030 §4: walk up to ``.git``; never count ``../..`` levels."""
    path = Path(start).resolve() if start is not None else Path(__file__).resolve()
    if path.is_file():
        path = path.parent
    for candidate in (path, *path.parents):
        if (candidate / ".git").exists():
            return candidate
    raise RecordError(f"no .git directory at or above {path}")


def git_info(root: Path | str) -> dict[str, Any]:
    """Record the git SHA/branch/dirty flag of the tree that produced the record."""
    root = Path(root)

    def _git(*args: str) -> str:
        try:
            out = subprocess.run(
                ["git", "-C", str(root), *args],
                capture_output=True, text=True, timeout=30, check=False,
            )
        except (OSError, subprocess.SubprocessError):
            return ""
        return out.stdout.strip()

    return {
        "sha": _git("rev-parse", "HEAD"),
        "branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "dirty": bool(_git("status", "--porcelain")),
        "root": str(root),
    }


# --------------------------------------------------------------------------- #
# frozen-config expectations
# --------------------------------------------------------------------------- #
def _as_float(value: Any) -> float | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    return None


def _same_number(a: Any, b: Any) -> bool:
    fa, fb = _as_float(a), _as_float(b)
    if fa is None or fb is None:
        return fa is None and fb is None
    return math.isclose(fa, fb, rel_tol=_REL_TOL, abs_tol=0.0)


def normalize_locked_config(raw: Mapping[str, Any]) -> dict[str, Any]:
    """Extract the frozen blocks this substep consumes from the locked config.

    The ``06.3`` file stores the canonical blocks at the top level
    (``structure``/``recipe``/``metric``/``lrs``/``seeds``/...).  Missing blocks
    fail closed: silently defaulting here would hide a changed contract.
    """
    if not isinstance(raw, Mapping):
        raise RecordError("locked config must be a JSON object")
    try:
        structure = dict(raw["structure"])
        recipe = dict(raw["recipe"])
        metric = dict(raw["metric"])
        lrs = dict(raw["lrs"])
        seeds = [int(s) for s in raw["seeds"]]
    except (KeyError, TypeError, ValueError) as exc:
        raise RecordError(f"locked config is missing a canonical block: {exc}") from exc
    return {
        "structure": structure,
        "recipe": recipe,
        "metric": metric,
        "lrs": lrs,
        "seeds": seeds,
        "loss": raw.get("loss"),
        "optimizer": raw.get("optimizer"),
        "prediction_noise": raw.get("prediction_noise"),
        "n_folds": raw.get("n_folds"),
        "contract_id": raw.get("contract_id"),
        "version": raw.get("version"),
    }


def ticket_expectation_errors(locked: Mapping[str, Any]) -> list[str]:
    """Compare the frozen locked config against the ticket's locked recipe.

    A non-empty list means the frozen input and this ticket disagree; the
    wrapper refuses to train on it (that is a contract change, not a knob).
    """
    errors: list[str] = []
    struct = locked.get("structure") or {}
    for key, expected in TICKET_EXPECTATIONS["structure"].items():
        actual = struct.get(key)
        if key in ("hidden", "adaptor"):
            actual = [int(x) for x in actual] if actual is not None else None
        if actual != expected:
            errors.append(f"locked structure.{key}={actual!r}, ticket expects {expected!r}")
    recipe = locked.get("recipe") or {}
    for key, expected in TICKET_EXPECTATIONS["recipe"].items():
        actual = recipe.get(key)
        if not _same_number(actual, expected):
            errors.append(f"locked recipe.{key}={actual!r}, ticket expects {expected!r}")
    metric = locked.get("metric") or {}
    for key, expected in TICKET_EXPECTATIONS["metric"].items():
        actual = metric.get(key)
        if key == "k":
            if _as_float(actual) is None or int(actual) != int(expected):
                errors.append(f"locked metric.k={actual!r}, ticket expects {expected!r}")
        elif actual != expected:
            errors.append(f"locked metric.{key}={actual!r}, ticket expects {expected!r}")
    for key in ("loss", "optimizer", "prediction_noise"):
        actual = locked.get(key)
        expected = TICKET_EXPECTATIONS[key]
        if actual != expected:
            errors.append(f"locked {key}={actual!r}, ticket expects {expected!r}")
    if not _same_number(locked.get("n_folds"), TICKET_EXPECTATIONS["n_folds"]):
        errors.append(f"locked n_folds={locked.get('n_folds')!r}, ticket expects 5")
    lrs = locked.get("lrs") or {}
    if not _same_number(lrs.get("pretrain"), TICKET_EXPECTATIONS["lrs"]["pretrain"]):
        errors.append(f"locked lrs.pretrain={lrs.get('pretrain')!r}, ticket expects 1e-3")
    seeds = [int(s) for s in (locked.get("seeds") or [])]
    if seeds != TICKET_EXPECTATIONS["seeds"]:
        errors.append(f"locked seeds={seeds!r}, ticket expects {TICKET_EXPECTATIONS['seeds']!r}")
    return errors


def _budget_from_recipe(recipe: Mapping[str, Any], metric: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "batch_size": int(recipe["batch_size"]),
        "max_epochs": int(recipe["max_epochs"]),
        "min_epochs": int(recipe["min_epochs"]),
        "early_stop_patience": int(recipe["early_stop_patience"]),
        "sched_factor": float(recipe["sched_factor"]),
        "sched_patience": int(recipe["sched_patience"]),
        "min_lr": float(recipe["min_lr"]),
        "weight_decay": float(recipe["weight_decay"]),
        "lr": float(recipe["lr"]),
        "metric": {"name": metric.get("name"), "k": int(metric.get("k", 0)),
                   "plateau_mode": metric.get("plateau_mode")},
    }


def frozen_input_identity(
    data_dir: Path | str,
    contract: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Read the frozen 06.3 manifest + contract into the record's identity block.

    Hashes are the frozen manifest's own values (proof of *which* artifacts the
    run must verify); the run manifest then records that verification happened.
    """
    data_dir = Path(data_dir)
    manifest_path = data_dir / "manifest.json"
    locked_path = data_dir / "locked_config.json"
    if not manifest_path.exists():
        raise RecordError(f"frozen input manifest is missing: {manifest_path}")
    if not locked_path.exists():
        raise RecordError(f"frozen locked config is missing: {locked_path}")
    manifest = read_json(manifest_path)
    locked_raw = read_json(locked_path)

    artifacts: dict[str, str] = {}
    sizes: dict[str, int] = {}
    block = manifest.get("artifacts")
    if isinstance(block, Mapping):
        items = block.items()
    elif isinstance(block, list):
        items = []
        for item in block:
            if not isinstance(item, Mapping) or "path" not in item:
                raise RecordError(f"frozen manifest artifact entry is malformed: {item!r}")
            items.append((item["path"], item))
    else:
        raise RecordError("frozen manifest has no artifacts block")
    for name, payload in items:
        if isinstance(payload, Mapping):
            digest = payload.get("sha256")
            if payload.get("bytes") is not None:
                sizes[str(name)] = int(payload["bytes"])
        else:
            digest = payload
        if not isinstance(digest, str) or len(digest) != 64:
            raise RecordError(f"frozen manifest artifact {name!r} has no sha256")
        artifacts[str(name)] = digest
    if "tahoe_transfer.h5ad" not in artifacts:
        raise RecordError("frozen manifest does not hash tahoe_transfer.h5ad")

    counts = manifest.get("counts") or {}
    tahoe_counts = counts.get("tahoe") or {}
    build = manifest.get("build") or {}
    contract_raw = dict(contract or {})
    contract_path = data_dir / "contract.json"
    if not contract_raw and contract_path.exists():
        contract_raw = read_json(contract_path)
    gene_axis = dict(contract_raw.get("gene_axis") or {})
    control_semantics = dict(contract_raw.get("control_semantics") or {})
    if not gene_axis.get("n_genes"):
        gene_axis["n_genes"] = int(locked_raw.get("structure", {}).get("n_genes", 0))

    return {
        "data_dir": str(data_dir.resolve()),
        "manifest_sha256": sha256_file(manifest_path),
        "locked_config_sha256": sha256_file(locked_path),
        "contract_sha256": sha256_file(contract_path) if contract_path.exists() else None,
        "artifacts": artifacts,
        "artifact_bytes": sizes,
        "gene_axis": {
            "n_genes": int(gene_axis.get("n_genes", 0)),
            "identity_columns": gene_axis.get("identity_columns"),
            "measurement_mask": gene_axis.get("measurement_mask"),
            "order": gene_axis.get("order"),
        },
        "control_groups": {
            "n_groups": tahoe_counts.get("control_groups"),
            "n_control_members": tahoe_counts.get("control_members"),
            "n_pairs": tahoe_counts.get("pairs"),
            "rule": control_semantics.get("tahoe_train_validation"),
        },
        "tahoe_split": {
            "valid_plates": build.get("tahoe_valid_plates"),
            "usable_train": tahoe_counts.get("usable_train"),
            "usable_valid": tahoe_counts.get("usable_valid"),
            "split_train_rows": tahoe_counts.get("split_train_rows"),
            "split_valid_rows": tahoe_counts.get("split_valid_rows"),
            "n_controls": tahoe_counts.get("controls"),
            "n_treated": tahoe_counts.get("treated"),
        },
        "counts_csv": counts,
    }


# --------------------------------------------------------------------------- #
# per-run evidence
# --------------------------------------------------------------------------- #
def _curve_from_history(history: Mapping[str, Any]) -> tuple[list[dict[str, Any]], list[str]]:
    errors: list[str] = []
    entries = history.get("history")
    if not isinstance(entries, list) or not entries:
        return [], ["fit_history.json has no per-epoch history"]
    curve: list[dict[str, Any]] = []
    for expected_epoch, entry in enumerate(entries, start=1):
        if not isinstance(entry, Mapping):
            errors.append(f"history entry {expected_epoch} is not an object")
            continue
        if int(entry.get("epoch", -1)) != expected_epoch:
            errors.append(
                f"history epoch sequence is broken at entry {expected_epoch} "
                f"(got {entry.get('epoch')!r})"
            )
        row: dict[str, Any] = {}
        for column in CURVE_COLUMNS:
            value = entry.get(column)
            if column == "epoch":
                row[column] = int(value) if value is not None else None
            else:
                row[column] = None if value is None else float(value)
        curve.append(row)
    return curve, errors


def run_entry(
    *,
    run_dir: Path | str,
    unit_id: str,
    seed: int,
    budget: str,
    command: Sequence[str],
    exit_code: int,
    expected_config_sha256: str | None,
    wall_s_outer: float | None = None,
    resources: Mapping[str, Any] | None = None,
    log_tail: str | None = None,
) -> dict[str, Any]:
    """Read one finished (or failed) run directory into a record unit entry."""
    run_dir = Path(run_dir)
    entry: dict[str, Any] = {
        "unit_id": unit_id,
        "seed": int(seed),
        "budget": budget,
        "command": [str(part) for part in command],
        "exit_code": int(exit_code),
        "status": "failed",
        "run_dir": str(run_dir),
        "expected_config_sha256": expected_config_sha256,
        "failures": [],
        "resources": dict(resources or {}),
        "timings": {"wall_s_outer": wall_s_outer},
    }
    manifest_path = run_dir / "manifest.json"
    if exit_code != 0:
        entry["failures"].append(f"trainer exited {exit_code}")
        if log_tail:
            entry["log_tail"] = log_tail[-4000:]
        return entry
    if not manifest_path.exists():
        entry["failures"].append(f"trainer exited 0 but wrote no manifest: {manifest_path}")
        if log_tail:
            entry["log_tail"] = log_tail[-4000:]
        return entry
    try:
        manifest = read_json(manifest_path)
        config = read_json(run_dir / "config.json") if (run_dir / "config.json").exists() else {}
        history = read_json(run_dir / "fit_history.json") if (run_dir / "fit_history.json").exists() else {}
    except RecordError as exc:
        entry["failures"].append(str(exc))
        return entry

    curve, curve_errors = _curve_from_history(history)
    entry["failures"].extend(curve_errors)
    checkpoint = manifest.get("checkpoint") or {}
    timings = manifest.get("timings") or {}
    timings_record = dict(timings)
    timings_record.update(entry["timings"])
    timings_record["wall_s_pipeline"] = timings.get("wall_s")
    timings_record["train_loop_s"] = timings.get("train_loop_s")
    if curve:
        train_seconds = [row["train_seconds"] for row in curve if row.get("train_seconds") is not None]
        valid_seconds = [row["valid_seconds"] for row in curve if row.get("valid_seconds") is not None]
        epoch_seconds = [
            (row.get("train_seconds") or 0.0) + (row.get("valid_seconds") or 0.0)
            for row in curve
        ]
        n_steps = (manifest.get("data_counts") or {}).get("n_train")
        timings_record["epoch_seconds"] = {
            "mean": sum(epoch_seconds) / len(epoch_seconds),
            "median": sorted(epoch_seconds)[len(epoch_seconds) // 2],
            "min": min(epoch_seconds),
            "max": max(epoch_seconds),
            "train_mean": sum(train_seconds) / len(train_seconds) if train_seconds else None,
            "valid_mean": sum(valid_seconds) / len(valid_seconds) if valid_seconds else None,
            "n_epochs": len(curve),
        }
        if n_steps and train_seconds:
            steps = math.ceil(int(n_steps) / int((config.get("recipe") or {}).get("batch_size", 512)))
            timings_record["steps_per_epoch"] = steps
            timings_record["s_per_step"] = sum(train_seconds) / (steps * len(train_seconds))

    entry.update(
        {
            "status": "completed" if not curve_errors else "invalid",
            "stop_reason": manifest.get("stop_reason"),
            "epochs_run": manifest.get("epochs_run"),
            "best_epoch": manifest.get("best_epoch"),
            "best_valid_auprc": manifest.get("best_valid_auprc"),
            "config": config,
            "config_sha256": manifest.get("config_sha256"),
            "config_matches_expected": (
                expected_config_sha256 is not None
                and manifest.get("config_sha256") == expected_config_sha256
            ),
            "manifest_sha256": sha256_file(manifest_path),
            "input": {
                "artifacts": ((manifest.get("data") or {}).get("artifacts") or {}),
                "gene_axis": ((manifest.get("data") or {}).get("gene_axis") or {}),
                "hash_verification": ((manifest.get("data") or {}).get("hash_verification") or {}),
                "data_counts": manifest.get("data_counts") or {},
                "cache": manifest.get("tahoe_cache") or {},
                "manifest_sha256": ((manifest.get("data") or {}).get("manifest_sha256")),
            },
            "checkpoint": {
                "path": str((run_dir / str(checkpoint.get("path", ""))).resolve())
                if checkpoint.get("path") else None,
                "path_rel": str(checkpoint.get("path")),
                "sha256": checkpoint.get("sha256"),
                "reload_identical": checkpoint.get("reload_identical"),
                "best_epoch": checkpoint.get("best_epoch"),
                "best_valid_auprc": checkpoint.get("best_valid_auprc"),
            },
            "model_card": {"path": str(run_dir / "model_card.md")},
            "budget_config": _budget_from_recipe(
                config.get("recipe") or {}, config.get("metric") or {}
            ),
            "validation": _curve_summary(curve),
            "validation_curve": curve,
            "timings": timings_record,
            "resources": {
                **dict(resources or {}),
                "host_peak_rss_mb_pipeline": manifest.get("peak_rss_mb"),
            },
            "device": manifest.get("device"),
            "git": manifest.get("git"),
            "metric": manifest.get("metric"),
            "noise_policy": manifest.get("noise_policy"),
            "init": manifest.get("init"),
        }
    )
    card = run_dir / "model_card.md"
    if card.exists():
        entry["model_card"]["sha256"] = sha256_file(card)
        entry["model_card"]["exists"] = True
    else:
        entry["model_card"]["exists"] = False
    return entry


def _curve_summary(curve: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    if not curve:
        return {"n_epochs": 0}
    auprc = [(row.get("epoch"), row.get("valid_auprc")) for row in curve]
    finite = [(epoch, value) for epoch, value in auprc if value is not None and math.isfinite(value)]
    best = max(finite, key=lambda pair: pair[1]) if finite else (None, None)
    first_argmax = next(
        (epoch for epoch, value in finite if value == max(v for _, v in finite)), None
    )
    return {
        "n_epochs": len(curve),
        "first": curve[0],
        "last": curve[-1],
        "best": {"epoch": best[0], "valid_auprc": best[1]},
        "best_epoch_is_first_argmax": best[0] == first_argmax,
    }


# --------------------------------------------------------------------------- #
# record assembly
# --------------------------------------------------------------------------- #
def build_record(
    *,
    mode: str,
    units: Sequence[Mapping[str, Any]],
    budget: Mapping[str, Any],
    expectations: Mapping[str, Any],
    input_identity: Mapping[str, Any],
    command: Sequence[str],
    git: Mapping[str, Any] | None = None,
    seeds_planned: Sequence[int] | None = None,
    seeds_requested: Sequence[int] | None = None,
    notes: Sequence[str] | None = None,
    host: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Assemble a record from finished unit entries (no I/O)."""
    if mode not in ("profile", "formal"):
        raise RecordError(f"record mode must be profile or formal, got {mode!r}")
    seeds_planned = [int(s) for s in (seeds_planned or expectations.get("seeds") or [])]
    seeds_requested = [int(s) for s in (seeds_requested or seeds_planned)]
    unit_ids = [str(unit.get("unit_id")) for unit in units]
    if len(set(unit_ids)) != len(unit_ids):
        raise RecordError(f"duplicate unit ids in the record: {unit_ids}")
    completed = [unit for unit in units if unit.get("status") == "completed"]
    failed = [unit for unit in units if unit.get("status") != "completed"]
    seeds_completed = sorted({int(unit["seed"]) for unit in completed})
    seeds_failed = sorted({int(unit["seed"]) for unit in failed})
    missing = [seed for seed in seeds_planned if seed not in seeds_completed]
    if mode == "formal":
        # one formal unit per frozen seed; any duplication is a cherry-pick smell
        formal_seeds = [int(unit["seed"]) for unit in units if unit.get("status") == "completed"]
        duplicated = sorted({seed for seed in formal_seeds if formal_seeds.count(seed) > 1})
        if duplicated:
            raise RecordError(f"duplicate completed seeds in a formal record: {duplicated}")
    complete = bool(units) and not failed and not missing
    by_seed: dict[str, Any] = {}
    for unit in units:
        if unit.get("status") != "completed":
            by_seed[str(unit["seed"])] = None
            continue
        value = unit.get("best_valid_auprc")
        if value is not None:
            previous = by_seed.get(str(unit["seed"]))
            by_seed[str(unit["seed"])] = max(previous, float(value)) if previous is not None else float(value)
    record = {
        "schema": SCHEMA,
        "kind": KIND,
        "mode": mode,
        "complete": complete,
        "created_utc": now_utc(),
        "git": dict(git or {}),
        "host": dict(host or {}),
        "command": [str(part) for part in command],
        "budget": dict(budget),
        "selection": {
            "metric": "direction_auprc",
            "k": 200,
            "plateau_mode": "max",
            "checkpoint_selection": "validation direction AUPRC@200 (max)",
            "validation_split": "frozen tahoe_split valid = plate13 + plate14 (matched "
                                "(plate, cell_line) DMSO control mean, log1p)",
            "test_usage": "the Tahoe test/valid split never informs LR, epoch, seed or capacity",
        },
        "expected": dict(expectations),
        "seeds_planned": seeds_planned,
        "seeds_requested": seeds_requested,
        "seeds_completed": seeds_completed,
        "seeds_failed": seeds_failed,
        "seeds_missing": missing,
        "no_seed_cherry_picking": True,
        "input": dict(input_identity),
        "units": [dict(unit) for unit in units],
        "notes": list(notes or []),
        "summary": {
            "n_units": len(units),
            "n_completed": len(completed),
            "n_failed": len(failed),
            "stop_reasons": {str(unit["unit_id"]): unit.get("stop_reason") for unit in units},
            "epochs_run": {str(unit["unit_id"]): unit.get("epochs_run") for unit in units},
            "best_valid_auprc": by_seed,
            "best_valid_auprc_by_unit": {
                str(unit["unit_id"]): unit.get("best_valid_auprc") for unit in units
            },
            "host_peak_rss_mb_max": max(
                (float((unit.get("resources") or {}).get("host_peak_rss_mb_pipeline") or 0.0)
                 for unit in units),
                default=0.0,
            ),
            "gpu_peak_mib_max": max(
                (
                    float(
                        (unit.get("resources") or {}).get("gpu_peak_mib_attributed")
                        or (unit.get("resources") or {}).get("gpu_peak_mib_process")
                        or (unit.get("resources") or {}).get("gpu_peak_mib_device")
                        or 0.0
                    )
                    for unit in units
                ),
                default=0.0,
            ),
            "gpu_peak_source": {
                str(unit.get("unit_id")): (
                    "process" if (unit.get("resources") or {}).get("gpu_peak_mib_process")
                    else "device-minus-baseline"
                    if (unit.get("resources") or {}).get("gpu_peak_mib_attributed")
                    else "device"
                    if (unit.get("resources") or {}).get("gpu_peak_mib_device")
                    else None
                )
                for unit in units
            },
            "wall_s_total": sum(
                float((unit.get("timings") or {}).get("wall_s_pipeline") or 0.0)
                for unit in units
            ),
        },
    }
    return record


def write_record(
    record: Mapping[str, Any],
    out_dir: Path | str,
    stem: str | None = None,
) -> dict[str, Path]:
    """Write ``<stem>.json`` + ``<stem>_curve.csv`` into *out_dir* (atomic JSON)."""
    out_dir = Path(out_dir)
    stem = stem or ("tahoe_pretrain_record" if record.get("mode") == "formal" else "profile_record")
    json_path = write_json(out_dir / f"{stem}.json", record)
    csv_path = out_dir / f"{stem}_curve.csv"
    lines = ["unit_id,seed,budget," + ",".join(CURVE_COLUMNS)]
    for unit in record.get("units", []):
        for row in unit.get("validation_curve") or []:
            values = [str(unit.get("unit_id")), str(unit.get("seed")), str(unit.get("budget"))]
            for column in CURVE_COLUMNS:
                value = row.get(column)
                values.append("" if value is None else repr(value) if isinstance(value, float) else str(value))
            lines.append(",".join(values))
    tmp = csv_path.with_suffix(csv_path.suffix + ".tmp")
    tmp.write_text("\n".join(lines) + "\n", encoding="utf-8")
    os.replace(tmp, csv_path)
    return {"json": json_path, "csv": csv_path}


def read_record(path: Path | str) -> dict[str, Any]:
    record = read_json(path)
    if record.get("schema") != SCHEMA:
        raise RecordError(f"record schema {record.get('schema')!r}, expected {SCHEMA!r}")
    if record.get("kind") != KIND:
        raise RecordError(f"record kind {record.get('kind')!r}, expected {KIND!r}")
    return record


# --------------------------------------------------------------------------- #
# checker
# --------------------------------------------------------------------------- #
def _check_locked(record: Mapping[str, Any], locked: Mapping[str, Any]) -> list[str]:
    """Re-derive the record's ``expected`` block from the frozen locked config."""
    errors: list[str] = []
    expected = record.get("expected") or {}
    for message in ticket_expectation_errors(locked):
        errors.append(f"frozen locked config: {message}")
    if expected.get("locked_config_sha256") != (record.get("input") or {}).get("locked_config_sha256"):
        errors.append("expected.locked_config_sha256 disagrees with input.locked_config_sha256")

    def _compare(block: str, got: Any, want: Any) -> None:
        if got != want:
            errors.append(f"expected.{block}={got!r}, locked config says {want!r}")

    # compare the canonical keys this substep consumes, not every decision-record
    # field the frozen file carries (e.g. metric.condition, metric.true_set)
    expected_structure = expected.get("structure") or {}
    locked_structure = locked.get("structure") or {}
    for key, want in expected_structure.items():
        if locked_structure.get(key) != want:
            errors.append(
                f"expected.structure.{key}={want!r}, locked config says "
                f"{locked_structure.get(key)!r}"
            )
    expected_metric = expected.get("metric") or {}
    locked_metric = locked.get("metric") or {}
    for key, want in expected_metric.items():
        got = locked_metric.get(key)
        if key == "k":
            if not _same_number(got, want):
                errors.append(f"expected.metric.k={want!r}, locked config says {got!r}")
        elif got != want:
            errors.append(f"expected.metric.{key}={want!r}, locked config says {got!r}")
    for key in ("loss", "optimizer", "prediction_noise"):
        _compare(key, expected.get(key), locked.get(key))
    if not _same_number(expected.get("n_folds"), locked.get("n_folds")):
        errors.append(f"expected.n_folds={expected.get('n_folds')!r}, locked {locked.get('n_folds')!r}")
    if [int(s) for s in (expected.get("seeds") or [])] != [int(s) for s in locked["seeds"]]:
        errors.append("expected.seeds disagree with the locked seeds")
    locked_budget = dict(locked.get("recipe") or {})
    locked_budget["lr"] = float((locked.get("lrs") or {}).get("pretrain", math.nan))
    recorded_budget = expected.get("budget") or {}
    for key, want in locked_budget.items():
        got = recorded_budget.get(key)
        if not _same_number(got, want):
            errors.append(f"expected.budget.{key}={got!r}, locked config says {want!r}")
    return errors


def _check_input(record: Mapping[str, Any], verify_frozen: bool) -> list[str]:
    errors: list[str] = []
    identity = record.get("input") or {}
    artifacts = identity.get("artifacts") or {}
    tahoe_sha = artifacts.get("tahoe_transfer.h5ad")
    if not isinstance(tahoe_sha, str) or len(tahoe_sha) != 64:
        errors.append("record input.artifacts has no tahoe_transfer.h5ad sha256")
    gene_axis = identity.get("gene_axis") or {}
    if int(gene_axis.get("n_genes", -1)) != int(TICKET_EXPECTATIONS["structure"]["n_genes"]):
        errors.append(f"record gene axis has {gene_axis.get('n_genes')!r} genes, expected 18,164")
    groups = identity.get("control_groups") or {}
    if not groups.get("n_groups"):
        errors.append("record input.control_groups is empty (matched (plate, cell_line) controls)")
    split = identity.get("tahoe_split") or {}
    plates = split.get("valid_plates")
    if plates is not None and list(plates) != ["plate13", "plate14"]:
        errors.append(f"record valid split plates {plates!r}, expected ['plate13', 'plate14']")
    cache = identity.get("cache") or {}
    if cache.get("source_sha256") != tahoe_sha:
        errors.append("derived cache was not built from the recorded frozen tahoe_transfer.h5ad")
    if verify_frozen:
        try:
            frozen = frozen_input_identity(str(identity.get("data_dir")))
        except RecordError as exc:
            errors.append(f"cannot re-read the frozen input: {exc}")
        else:
            if frozen["manifest_sha256"] != identity.get("manifest_sha256"):
                errors.append(
                    "frozen input manifest changed since the record was written "
                    f"({identity.get('manifest_sha256')} -> {frozen['manifest_sha256']})"
                )
            for name, digest in artifacts.items():
                if name in frozen["artifacts"] and frozen["artifacts"][name] != digest:
                    errors.append(f"frozen artifact {name} sha256 changed since the record")
            if frozen["tahoe_split"].get("valid_plates") and split.get("valid_plates"):
                if list(frozen["tahoe_split"]["valid_plates"]) != list(split["valid_plates"]):
                    errors.append("record valid split plates disagree with the frozen manifest")
    return errors


def _check_unit(unit: Mapping[str, Any], record: Mapping[str, Any], *, strict: bool) -> list[str]:
    errors: list[str] = []
    expected = record.get("expected") or {}
    mode = record.get("mode")
    uid = unit.get("unit_id")
    if unit.get("status") not in ("completed", "invalid", "failed"):
        errors.append(f"{uid}: unknown status {unit.get('status')!r}")
    if unit.get("status") != "completed":
        return errors
    # the unit's own resolved budget (a profile record may mix smoke + full-width)
    unit_budget = unit.get("budget_config") or record.get("budget") or {}
    reference_budget = dict(expected.get("budget") or {})
    fingerprint = expected.get("config_sha256") or {}
    want_fp = (
        fingerprint.get(str(uid))
        or fingerprint.get(str(unit.get("seed")))
        or fingerprint.get(unit.get("seed"))
    )
    got_fp = unit.get("config_sha256")
    if want_fp is None:
        errors.append(f"{uid}: no expected config digest recorded for seed {unit.get('seed')}")
    elif got_fp != want_fp:
        errors.append(f"{uid}: config digest {got_fp} != expected {want_fp}")
    if unit.get("config_matches_expected") is not True:
        errors.append(f"{uid}: trainer manifest config digest does not match the plan")

    stop = unit.get("stop_reason")
    patience = int(unit_budget.get("early_stop_patience", -1))
    allowed = {f"early_stop_no_auprc_improvement_{patience}", "max_epochs", "no_finite_auprc"}
    if stop not in allowed:
        errors.append(f"{uid}: unexpected stop reason {stop!r} (allowed {sorted(allowed)})")
    epochs = unit.get("epochs_run")
    if not isinstance(epochs, int) or epochs < 1:
        errors.append(f"{uid}: epochs_run {epochs!r} is not a positive integer")
    else:
        max_epochs = int(unit_budget.get("max_epochs", -1))
        min_epochs = int(unit_budget.get("min_epochs", -1))
        if max_epochs > 0 and epochs > max_epochs:
            errors.append(f"{uid}: epochs_run {epochs} exceeds max_epochs {max_epochs}")
        if stop == "max_epochs" and max_epochs > 0 and epochs != max_epochs:
            errors.append(f"{uid}: stop_reason max_epochs but epochs_run {epochs} != {max_epochs}")
        if isinstance(stop, str) and stop.startswith("early_stop") and epochs < min_epochs:
            errors.append(f"{uid}: early stop at epoch {epochs} before min_epochs {min_epochs}")

    curve = unit.get("validation_curve") or []
    if len(curve) != epochs:
        errors.append(f"{uid}: curve has {len(curve)} rows but epochs_run is {epochs}")
    finite_auprc = []
    for index, row in enumerate(curve, start=1):
        if int(row.get("epoch", -1)) != index:
            errors.append(f"{uid}: curve epoch sequence broken at row {index}")
        for column in ("train_loss", "valid_loss", "valid_auprc", "lr"):
            value = row.get(column)
            if value is None or not math.isfinite(float(value)):
                errors.append(f"{uid}: curve epoch {index} {column} is not finite ({value!r})")
        auprc = row.get("valid_auprc")
        if auprc is not None and math.isfinite(float(auprc)):
            finite_auprc.append((index, float(auprc)))
    if finite_auprc:
        best_epoch, best_value = max(finite_auprc, key=lambda pair: pair[1])
        if unit.get("best_epoch") != best_epoch:
            errors.append(
                f"{uid}: best_epoch {unit.get('best_epoch')!r} is not the curve argmax {best_epoch}"
            )
        recorded_best = unit.get("best_valid_auprc")
        if recorded_best is None or not math.isclose(
            float(recorded_best), best_value, rel_tol=1e-9, abs_tol=1e-12
        ):
            errors.append(
                f"{uid}: best_valid_auprc {recorded_best!r} != curve max {best_value!r}"
            )
        if not unit.get("validation", {}).get("best_epoch_is_first_argmax", True):
            errors.append(f"{uid}: best epoch is not the first occurrence of the maximum")

    # the resolved recipe must be the locked reference budget, except for the
    # documented profile overrides (recorded explicitly on the unit)
    if mode == "formal":
        for key, want in reference_budget.items():
            got = unit_budget.get(key)
            if not _same_number(got, want):
                errors.append(f"{uid}: budget.{key}={got!r}, locked reference says {want!r}")
    else:
        overrides = (record.get("budget") or {}).get("overrides_by_unit", {}).get(str(uid), {})
        for key in ("batch_size", "early_stop_patience", "sched_factor", "sched_patience",
                    "min_lr", "weight_decay", "lr"):
            if not _same_number(unit_budget.get(key), reference_budget.get(key)):
                if str(key) not in {str(k) for k in overrides}:
                    errors.append(
                        f"{uid}: budget.{key} differs from the locked reference but is not a "
                        "recorded profile override"
                    )

    checkpoint = unit.get("checkpoint") or {}
    sha = checkpoint.get("sha256")
    if not isinstance(sha, str) or len(sha) != 64:
        errors.append(f"{uid}: checkpoint sha256 is missing or malformed")
    if checkpoint.get("reload_identical") is not True:
        errors.append(f"{uid}: the trainer did not confirm a strict checkpoint reload")
    if strict:
        path = checkpoint.get("path")
        if not path or not Path(path).exists():
            errors.append(f"{uid}: checkpoint file is missing: {path}")
        elif isinstance(sha, str) and sha256_file(path) != sha:
            errors.append(f"{uid}: checkpoint file sha256 does not match the record")
        card = unit.get("model_card") or {}
        if not card.get("exists"):
            errors.append(f"{uid}: run model_card.md is missing")

    resources = unit.get("resources") or {}
    pipeline_peak = resources.get("host_peak_rss_mb_pipeline")
    if not isinstance(pipeline_peak, (int, float)) or float(pipeline_peak) <= 0:
        errors.append(f"{uid}: host peak RSS is not recorded")
    if strict:
        gpu_peak = (
            resources.get("gpu_peak_mib_attributed")
            or resources.get("gpu_peak_mib_process")
            or resources.get("gpu_peak_mib_device")
        )
        if not isinstance(gpu_peak, (int, float)) or float(gpu_peak) <= 0:
            errors.append(f"{uid}: GPU peak is not recorded")
        timings = unit.get("timings") or {}
        wall = timings.get("wall_s_pipeline")
        if not isinstance(wall, (int, float)) or float(wall) <= 0:
            errors.append(f"{uid}: pipeline wall time is not recorded")
        if mode == "formal" and not (timings.get("steps_per_epoch") or 0) > 0:
            errors.append(f"{uid}: steps per epoch is not recorded")

    # input counts: a formal full-width run must use the whole frozen valid split
    data_counts = ((unit.get("input") or {}).get("data_counts") or {})
    split = (record.get("input") or {}).get("tahoe_split") or {}
    if mode == "formal":
        expected_pairs = {
            "n_train": split.get("usable_train"),
            "n_valid": split.get("usable_valid"),
            "control_groups": (record.get("input") or {}).get("control_groups", {}).get("n_groups"),
        }
        for key, want in expected_pairs.items():
            got = data_counts.get(key)
            if want is not None and got != want:
                errors.append(f"{uid}: data_counts.{key}={got!r}, frozen contract says {want!r}")
        if data_counts.get("smoke_subset"):
            errors.append(f"{uid}: a formal run used a row subset")
    return errors


def check_record(
    record: Mapping[str, Any],
    *,
    verify_frozen: bool = True,
    locked: Mapping[str, Any] | None = None,
    expect_mode: str | None = None,
    strict: bool = True,
    expect_complete: bool = True,
) -> tuple[list[str], list[str]]:
    """Validate a record.  Returns ``(errors, warnings)``; empty errors = pass."""
    errors: list[str] = []
    warnings: list[str] = []
    if record.get("schema") != SCHEMA:
        errors.append(f"schema {record.get('schema')!r}, expected {SCHEMA!r}")
        return errors, warnings
    if record.get("kind") != KIND:
        errors.append(f"kind {record.get('kind')!r}, expected {KIND!r}")
    mode = record.get("mode")
    if mode not in ("profile", "formal"):
        errors.append(f"mode {mode!r} is neither profile nor formal")
        return errors, warnings
    if expect_mode is not None and mode != expect_mode:
        errors.append(f"record mode is {mode!r}, expected {expect_mode!r}")
    if mode == "profile":
        warnings.append("PROFILE record: pipeline/throughput evidence, not a formal pretrain")
        if expect_mode == "formal":
            errors.append("a profile record can never stand in for the formal pretrain record")
        if record.get("no_seed_cherry_picking") is not True:
            warnings.append("profile record does not assert no_seed_cherry_picking")
        if expect_complete and not record.get("complete"):
            warnings.append("profile record is incomplete (a profile unit failed)")
    else:
        if record.get("complete") is not True and expect_complete:
            errors.append(
                "formal record is incomplete "
                f"(seeds_missing={record.get('seeds_missing')}, "
                f"seeds_failed={record.get('seeds_failed')})"
            )
        if record.get("no_seed_cherry_picking") is not True:
            errors.append("formal record must assert no_seed_cherry_picking")
        formal = [unit for unit in record.get("units", []) if unit.get("budget") == "formal"]
        if expect_complete and len(formal) != 3:
            errors.append(f"formal record has {len(formal)} formal units, expected 3")
        planned = [int(s) for s in record.get("seeds_planned") or []]
        locked_seeds = [int(s) for s in (locked or {}).get("seeds", TICKET_EXPECTATIONS["seeds"])]
        if planned != locked_seeds:
            errors.append(f"formal record plans seeds {planned}, frozen seeds are {locked_seeds}")
        unit_seeds = [int(unit.get("seed")) for unit in formal if unit.get("status") == "completed"]
        if expect_complete and sorted(unit_seeds) != sorted(planned):
            errors.append(
                f"formal record completed seeds {sorted(set(unit_seeds))}, planned {sorted(planned)}"
            )

    budget = record.get("budget") or {}
    budget_kind = budget.get("kind")
    if mode == "formal" and budget_kind not in (None, "locked"):
        errors.append(f"formal record has budget.kind={budget_kind!r}, expected 'locked'")
    if mode == "profile" and budget_kind == "locked":
        errors.append("profile record carries the locked formal budget kind")
    if mode == "formal":
        if budget.get("smoke"):
            errors.append("formal record has smoke=True")
        if budget.get("row_subset"):
            errors.append("formal record used a row subset")
        if budget.get("max_epochs") != TICKET_EXPECTATIONS["recipe"]["max_epochs"]:
            errors.append(
                f"formal budget max_epochs={budget.get('max_epochs')!r}, "
                f"locked is {TICKET_EXPECTATIONS['recipe']['max_epochs']}"
            )
        if budget.get("min_epochs") != TICKET_EXPECTATIONS["recipe"]["min_epochs"]:
            errors.append(
                f"formal budget min_epochs={budget.get('min_epochs')!r}, "
                f"locked is {TICKET_EXPECTATIONS['recipe']['min_epochs']}"
            )
        for key, expected in TICKET_EXPECTATIONS["recipe"].items():
            if key in ("max_epochs", "min_epochs"):
                continue
            if not _same_number(budget.get(key), expected):
                errors.append(f"formal budget {key}={budget.get(key)!r}, ticket expects {expected!r}")
        if not _same_number(budget.get("lr"), TICKET_EXPECTATIONS["lrs"]["pretrain"]):
            errors.append(
                f"formal budget lr={budget.get('lr')!r}, locked lrs.pretrain is 1e-3"
            )
        overrides = budget.get("overrides") or {}
        if overrides:
            errors.append(f"formal run applied recipe overrides: {sorted(overrides)}")

    if not record.get("units"):
        errors.append("record has no units")
    for unit in record.get("units", []):
        errors.extend(_check_unit(unit, record, strict=strict))

    if locked is not None:
        errors.extend(_check_locked(record, locked))
    errors.extend(_check_input(record, verify_frozen))

    # cross-run evidence: same input identity, distinct checkpoints
    checkpoint_shas = [
        (unit.get("checkpoint") or {}).get("sha256")
        for unit in record.get("units", [])
        if (unit.get("checkpoint") or {}).get("sha256")
    ]
    if len(set(checkpoint_shas)) != len(checkpoint_shas):
        errors.append("two units share a checkpoint sha256 (the record is not three distinct runs)")
    return errors, warnings


def resource_report(
    *,
    host_peak_rss_mb_pipeline: float | None,
    host_peak_rss_mb_sampler: float | None,
    gpu_peak_mib_device: float | None,
    gpu_peak_mib_process: float | None,
    gpu_baseline_mib: float | None,
    gpu_after_mib: float | None,
    samples: int,
    interval_s: float,
    gpu_process_attribution: str | None = None,
) -> dict[str, Any]:
    """Normalise one run's measured resources into the record's shape.

    ``gpu_peak_mib_attributed`` is the per-process peak when the driver can
    report it, else the device peak minus the baseline (only meaningful when the
    run is alone on the GPU — which formal runs must be).  The raw device and
    baseline numbers are kept so the reader can judge the attribution.
    """
    device_attributed = None
    if gpu_peak_mib_device is not None and gpu_baseline_mib is not None:
        device_attributed = max(0.0, float(gpu_peak_mib_device) - float(gpu_baseline_mib))
    attributed = gpu_peak_mib_process if gpu_peak_mib_process else device_attributed
    return {
        "host_peak_rss_mb_pipeline": host_peak_rss_mb_pipeline,
        "host_peak_rss_mb_sampler": host_peak_rss_mb_sampler,
        "gpu_peak_mib_attributed": attributed,
        "gpu_peak_mib_device": gpu_peak_mib_device,
        "gpu_peak_mib_process": gpu_peak_mib_process,
        "gpu_peak_mib_device_minus_baseline": device_attributed,
        "gpu_baseline_mib": gpu_baseline_mib,
        "gpu_after_mib": gpu_after_mib,
        "gpu_process_attribution": gpu_process_attribution,
        "samples": int(samples),
        "interval_s": float(interval_s),
    }
