"""Pipeline configuration: locked-config resolution and validation.

Two layers exist and must not be confused:

* the **frozen locked config** produced by ticket #3 (``06.3_input/output/
  locked_config.json``) — structure, dose encoding, recipe budgets, metric
  definition, seeds and LR candidates.  It is read-only evidence; this module
  normalises it (field aliases are accepted and recorded) and validates it.
* the **resolved run config** — the frozen config plus explicit CLI overrides
  (smoke budgets, ``--structure`` presets for the #10 capacity ablations,
  ``--train-drugs`` for #11).  Every resolved run records both.

Nothing here guesses silently: missing required values fail closed, and every
alias/default used is returned in ``notes`` so the run manifest can record how
the frozen config was interpreted.

Canonical structure (the parameter shapes depend on all of these)::

    n_genes   int    fixed gene axis length (18,164 in the study)
    hidden    [int]  encoder/decoder hidden widths (one layer: [128])
    z_dim     int    state latent (64)
    c_dim     int    drug latent (64)
    n_dim     int    noise dimension (10)
    drug_dim  int    dose-fingerprint dimension (1024)
    adaptor   [int]  drug-adaptor hidden widths ([128])
    dropout   float  dropout rate (0.05)

Canonical recipe: ``batch_size``, ``max_epochs``, ``min_epochs``,
``early_stop_patience``, ``sched_factor``, ``sched_patience``, ``min_lr``,
``weight_decay``.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

__all__ = [
    "ConfigError",
    "CANONICAL_STRUCTURE_DEFAULTS",
    "CANONICAL_RECIPE_DEFAULTS",
    "STRUCTURE_PRESETS",
    "normalize_locked_config",
    "normalize_structure",
    "read_locked_config",
    "load_json",
    "validate_structure",
    "validate_recipe",
    "validate_config",
    "structure_from_spec",
    "structure_compatible",
    "apply_recipe_overrides",
    "canonical_json",
    "sha256_payload",
    "config_fingerprint",
]


class ConfigError(ValueError):
    """A configuration is missing, inconsistent or outside its contract."""


CANONICAL_STRUCTURE_DEFAULTS: dict[str, Any] = {
    "n_genes": 18164,
    "hidden": [128],
    "z_dim": 64,
    "c_dim": 64,
    "n_dim": 10,
    "drug_dim": 1024,
    "adaptor": [128],
    "dropout": 0.05,
}

CANONICAL_RECIPE_DEFAULTS: dict[str, Any] = {
    "batch_size": 512,
    "max_epochs": 100,
    "min_epochs": 20,
    "early_stop_patience": 20,
    "sched_factor": 0.5,
    "sched_patience": 10,
    "min_lr": 1e-8,
    "weight_decay": 1e-8,
}

#: Named structures for the capacity ablation (#10).  ``hidden128`` is the
#: first-round baseline; the other two are the audited post-first-round
#: ablations (state latent 128, hidden width 256).
STRUCTURE_PRESETS: dict[str, dict[str, Any]] = {
    "hidden128": {"hidden": [128]},
    "latent128": {"z_dim": 128},
    "hidden256": {"hidden": [256]},
}

_STRUCTURE_ALIASES: dict[str, tuple[str, ...]] = {
    "n_genes": ("n_genes", "in_dim", "x_dimension", "x_dim", "n_genes_axis"),
    "hidden": ("hidden", "hidden_layer_sizes", "hidden_sizes", "hidden_width"),
    "z_dim": ("z_dim", "z_dimension", "state_latent", "state_latent_dim"),
    "c_dim": ("c_dim", "c_dimension", "comb_dimension", "drug_latent", "drug_latent_dim"),
    "n_dim": ("n_dim", "noise", "noise_dim", "noise_dimension"),
    "drug_dim": ("drug_dim", "drug_dimension", "n_bits", "fingerprint_bits"),
    "adaptor": ("adaptor", "adaptor_layer_sizes", "adaptor_sizes"),
    "dropout": ("dropout", "dr_rate", "drop_rate", "dropout_rate"),
}

_RECIPE_ALIASES: dict[str, tuple[str, ...]] = {
    "batch_size": ("batch_size", "batch"),
    "max_epochs": ("max_epochs", "epochs", "max_epoch"),
    "min_epochs": ("min_epochs", "min_epoch"),
    "early_stop_patience": ("early_stop_patience", "patience", "early_stopping_patience"),
    "sched_factor": ("sched_factor", "plateau_factor", "lr_factor"),
    "sched_patience": ("sched_patience", "plateau_patience", "lr_patience"),
    "min_lr": ("min_lr", "scheduler_min_lr"),
    "weight_decay": ("weight_decay", "wd"),
}

#: spellings that are accepted without being reported as an alias mapping
_ACCEPTED: dict[str, tuple[str, ...]] = {
    "metric_k": ("metric_k", "k"),
}

_METRIC_ALIASES: dict[str, tuple[str, ...]] = {
    "metric_k": ("metric_k", "k", "auprc_k", "top_k"),
    "metric_name": ("metric", "metric_name", "selection_metric"),
}

_TOP_ALIASES: dict[str, tuple[str, ...]] = {
    "seeds": ("seeds", "seed_list", "fixed_seeds"),
    "n_folds": ("n_folds", "folds", "num_folds"),
}


def load_json(path: Path | str) -> dict[str, Any]:
    """Read a JSON object, failing closed on anything else."""
    p = Path(path)
    try:
        payload = json.loads(p.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ConfigError(f"config file not found: {p}") from exc
    except json.JSONDecodeError as exc:
        raise ConfigError(f"config file is not valid JSON: {p} ({exc})") from exc
    if not isinstance(payload, dict):
        raise ConfigError(f"config file must hold a JSON object: {p}")
    return payload


def _get_alias(mapping: Mapping[str, Any], aliases: Iterable[str]) -> tuple[Any, str | None]:
    for alias in aliases:
        if alias in mapping and mapping[alias] is not None:
            return mapping[alias], alias
    return None, None


def _as_int_list(value: Any, name: str) -> list[int]:
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        value = [value]
    if not isinstance(value, (list, tuple)):
        raise ConfigError(f"{name} must be an int or a list of ints, got {value!r}")
    out: list[int] = []
    for item in value:
        if isinstance(item, bool) or int(item) != item:
            raise ConfigError(f"{name} must contain integers, got {item!r}")
        out.append(int(item))
    return out


def _normalize_structure(raw: Mapping[str, Any], notes: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise ConfigError(f"structure block must be an object, got {type(raw).__name__}")
    struct: dict[str, Any] = {}
    for canonical, aliases in _STRUCTURE_ALIASES.items():
        value, used = _get_alias(raw, aliases)
        if value is None:
            continue
        if used != canonical:
            notes["mapping"].append(f"structure.{used} -> structure.{canonical}")
        struct[canonical] = value
    # Normalise list-valued fields.
    if "hidden" in struct:
        struct["hidden"] = _as_int_list(struct["hidden"], "hidden")
    if "adaptor" in struct:
        struct["adaptor"] = _as_int_list(struct["adaptor"], "adaptor")
    for key in ("n_genes", "z_dim", "c_dim", "n_dim", "drug_dim"):
        if key in struct:
            value = struct[key]
            if isinstance(value, bool) or int(value) != value:
                raise ConfigError(f"structure.{key} must be an integer, got {value!r}")
            struct[key] = int(value)
    if "dropout" in struct:
        struct["dropout"] = float(struct["dropout"])
    return struct


def _normalize_recipe(raw: Mapping[str, Any], notes: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise ConfigError(f"recipe block must be an object, got {type(raw).__name__}")
    recipe: dict[str, Any] = {}
    for canonical, aliases in _RECIPE_ALIASES.items():
        value, used = _get_alias(raw, aliases)
        if value is None:
            continue
        if used != canonical:
            notes["mapping"].append(f"recipe.{used} -> recipe.{canonical}")
        if canonical in {"max_epochs", "min_epochs", "batch_size", "early_stop_patience",
                         "sched_patience"}:
            if isinstance(value, bool) or int(value) != value:
                raise ConfigError(f"recipe.{canonical} must be an integer, got {value!r}")
            recipe[canonical] = int(value)
        else:
            recipe[canonical] = float(value)
    return recipe


def _parse_metric_name(name: Any) -> int | None:
    """Accept ``direction_auprc@200`` / ``direction AUPRC@200`` and return 200."""
    if not isinstance(name, str):
        return None
    lowered = name.lower().replace(" ", "_").replace("-", "_")
    if "auprc" not in lowered:
        raise ConfigError(
            f"selection metric must be direction AUPRC (the frozen study metric), got {name!r}"
        )
    if "@" in lowered:
        tail = lowered.split("@", 1)[1]
        digits = "".join(ch for ch in tail if ch.isdigit())
        if digits:
            return int(digits)
    digits = "".join(ch for ch in lowered if ch.isdigit())
    return int(digits) if digits else None


def normalize_locked_config(raw: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Normalise a frozen locked config into the canonical pipeline config.

    Returns ``(canonical, notes)`` where ``notes`` records every field alias and
    every documented default that was used, so the run manifest can state how
    the frozen artifact was interpreted.
    """
    if not isinstance(raw, Mapping):
        raise ConfigError("locked config must be a JSON object")
    notes: dict[str, Any] = {"mapping": [], "defaulted": []}

    structure_raw, _ = _get_alias(raw, ("structure", "model", "model_structure"))
    if structure_raw is None:
        # Some builds inline the structure keys at the top level.
        structure_raw = raw
    struct = _normalize_structure(structure_raw, notes)

    recipe_raw, _ = _get_alias(raw, ("recipe", "optim", "optimizer_config", "optimization"))
    if recipe_raw is None:
        recipe_raw = raw
    recipe = _normalize_recipe(recipe_raw, notes)

    for key, value in CANONICAL_STRUCTURE_DEFAULTS.items():
        if key not in struct:
            struct[key] = list(value) if isinstance(value, list) else value
            notes["defaulted"].append(f"structure.{key}={value}")

    metric_raw, _ = _get_alias(raw, ("metric", "metric_definition", "selection"))
    metric_name: Any = None
    metric_k: Any = None
    if isinstance(metric_raw, Mapping):
        metric_k, used_k = _get_alias(metric_raw, _METRIC_ALIASES["metric_k"])
        metric_name, used_n = _get_alias(metric_raw, _METRIC_ALIASES["metric_name"])
        for canonical, used in (("metric_k", used_k), ("metric_name", used_n)):
            if used is not None and used != canonical and used not in _ACCEPTED.get(canonical, ()):
                notes["mapping"].append(f"metric.{used} -> metric.{canonical}")
    else:
        metric_name = metric_raw
    k_from_name = _parse_metric_name(metric_name) if metric_name is not None else None

    cfg: dict[str, Any] = dict(raw)
    cfg["structure"] = struct
    cfg["recipe"] = recipe

    metric_out: dict[str, Any] = {"name": "direction_auprc", "plateau_mode": "max"}
    if metric_k is not None:
        if isinstance(metric_k, bool) or int(metric_k) != metric_k:
            raise ConfigError(f"metric k must be an integer, got {metric_k!r}")
        metric_out["k"] = int(metric_k)
    elif k_from_name is not None:
        metric_out["k"] = k_from_name
        notes["defaulted"].append("metric.k from the metric name")
    cfg["metric"] = metric_out

    for canonical, aliases in _TOP_ALIASES.items():
        value, used = _get_alias(raw, aliases)
        if value is None:
            continue
        if used != canonical:
            notes["mapping"].append(f"{used} -> {canonical}")
        cfg[canonical] = value

    # Dose encoding: the frozen contract fixes it; accept a nested block too.
    dose_raw, _ = _get_alias(raw, ("dose_encoding", "dose", "drug_encoding"))
    if isinstance(dose_raw, Mapping):
        cfg["dose_encoding"] = dict(dose_raw)
    else:
        cfg["dose_encoding"] = {
            "fingerprint": "FCFP4",
            "n_bits": struct.get("drug_dim", CANONICAL_STRUCTURE_DEFAULTS["drug_dim"]),
            "radius": 2,
            "use_features": True,
            "dose_scale": "log10(dose_uM+1)",
        }
        notes["defaulted"].append("dose_encoding (frozen contract: FCFP4 r2 1024 x log10(uM+1))")

    # Fill missing recipe values with the frozen spec defaults, recording each.
    for key, value in CANONICAL_RECIPE_DEFAULTS.items():
        if key not in cfg["recipe"]:
            cfg["recipe"][key] = value
            notes["defaulted"].append(f"recipe.{key}={value}")

    for key, value in (
        ("loss", "gaussian_nll"),
        ("optimizer", "adam"),
        ("prediction_noise", "zero"),
    ):
        if key not in cfg:
            cfg[key] = value
            notes["defaulted"].append(f"{key}={value}")
    if "n_folds" not in cfg:
        cfg["n_folds"] = 5
        notes["defaulted"].append("n_folds=5")
    if "seeds" not in cfg:
        cfg["seeds"] = [20261021, 20261022, 20261023]
        notes["defaulted"].append("seeds=[20261021, 20261022, 20261023]")
    lrs_raw, _ = _get_alias(raw, ("lrs", "learning_rates", "lr"))
    if not isinstance(lrs_raw, Mapping):
        lrs_raw = {}
    if "pretrain" not in lrs_raw:
        lrs_raw = {**lrs_raw, "pretrain": 1e-3}
        notes["defaulted"].append("lrs.pretrain=1e-3")
    if "cm" not in lrs_raw:
        lrs_raw = {**lrs_raw, "cm": [1e-4, 3e-4, 1e-3]}
        notes["defaulted"].append("lrs.cm=[1e-4,3e-4,1e-3]")
    cfg["lrs"] = {k: v for k, v in lrs_raw.items()}

    validate_config(cfg)
    return cfg, notes


def normalize_structure(raw: Mapping[str, Any]) -> dict[str, Any]:
    """Normalise a structure block (or a config that inlines one) to canonical keys."""
    if not isinstance(raw, Mapping):
        raise ConfigError(f"structure block must be an object, got {type(raw).__name__}")
    notes: dict[str, Any] = {"mapping": [], "defaulted": []}
    return _normalize_structure(raw, notes)


def read_locked_config(path: Path | str) -> tuple[dict[str, Any], dict[str, Any]]:
    return normalize_locked_config(load_json(path))


def validate_structure(struct: Mapping[str, Any]) -> None:
    """Fail loudly on a structure that cannot build the audited PRnet ``PGM``."""
    for key in ("n_genes", "hidden", "z_dim", "c_dim", "n_dim", "drug_dim", "adaptor", "dropout"):
        if key not in struct:
            raise ConfigError(f"structure.{key} is missing")
    if int(struct["n_genes"]) < 1:
        raise ConfigError(f"structure.n_genes must be positive, got {struct['n_genes']}")
    for key in ("z_dim", "c_dim", "n_dim", "drug_dim"):
        if int(struct[key]) < 1:
            raise ConfigError(f"structure.{key} must be positive, got {struct[key]}")
    for key in ("hidden", "adaptor"):
        widths = struct[key]
        if not widths:
            raise ConfigError(f"structure.{key} must be a non-empty list")
        if any(int(w) < 1 for w in widths):
            raise ConfigError(f"structure.{key} must hold positive widths, got {widths}")
    if not 0.0 <= float(struct["dropout"]) < 1.0:
        raise ConfigError(f"structure.dropout must be in [0, 1), got {struct['dropout']}")


def validate_recipe(recipe: Mapping[str, Any]) -> None:
    required = (
        "batch_size", "max_epochs", "min_epochs", "early_stop_patience",
        "sched_factor", "sched_patience", "min_lr", "weight_decay",
    )
    for key in required:
        if key not in recipe:
            raise ConfigError(f"recipe.{key} is missing")
    if int(recipe["batch_size"]) < 2:
        raise ConfigError(
            f"recipe.batch_size must be >= 2 (BatchNorm1d in train mode), got {recipe['batch_size']}"
        )
    if int(recipe["max_epochs"]) < 1:
        raise ConfigError("recipe.max_epochs must be >= 1")
    if not 1 <= int(recipe["min_epochs"]) <= int(recipe["max_epochs"]):
        raise ConfigError(
            f"recipe.min_epochs ({recipe['min_epochs']}) must be in 1..max_epochs "
            f"({recipe['max_epochs']})"
        )
    if int(recipe["early_stop_patience"]) < 0 or int(recipe["sched_patience"]) < 0:
        raise ConfigError("patience values must be non-negative")
    if not 0.0 < float(recipe["sched_factor"]) < 1.0:
        raise ConfigError(
            f"recipe.sched_factor must be in (0, 1), got {recipe['sched_factor']}"
        )
    if float(recipe["min_lr"]) < 0.0:
        raise ConfigError("recipe.min_lr must be non-negative")
    if float(recipe["weight_decay"]) < 0.0:
        raise ConfigError("recipe.weight_decay must be non-negative")


def validate_config(cfg: Mapping[str, Any]) -> None:
    """Cross-block validation of the canonical config."""
    struct = cfg.get("structure")
    if not isinstance(struct, Mapping):
        raise ConfigError("config.structure is missing")
    validate_structure(struct)
    recipe = cfg.get("recipe")
    if not isinstance(recipe, Mapping):
        raise ConfigError("config.recipe is missing")
    validate_recipe(recipe)

    metric = cfg.get("metric")
    if not isinstance(metric, Mapping):
        raise ConfigError("config.metric is missing")
    if "k" not in metric:
        raise ConfigError("config.metric.k (the direction AUPRC k) is missing")
    k = metric["k"]
    if isinstance(k, bool) or int(k) != k or int(k) < 1:
        raise ConfigError(f"config.metric.k must be a positive integer, got {k!r}")
    if int(k) > int(struct["n_genes"]):
        raise ConfigError(
            f"config.metric.k={k} exceeds the gene axis ({struct['n_genes']} genes)"
        )
    if metric.get("name") != "direction_auprc":
        raise ConfigError(
            "config.metric.name must be 'direction_auprc' (the frozen study metric)"
        )
    if metric.get("plateau_mode") != "max":
        raise ConfigError("config.metric.plateau_mode must be 'max' (AUPRC is maximised)")

    if cfg.get("loss") != "gaussian_nll":
        raise ConfigError("config.loss must be 'gaussian_nll' (the author GUSS loss)")
    if cfg.get("optimizer") != "adam":
        raise ConfigError("config.optimizer must be 'adam'")
    if cfg.get("prediction_noise") not in (None, "zero"):
        raise ConfigError(
            "config.prediction_noise must be 'zero': validation/inference are deterministic"
        )
    n_folds = cfg.get("n_folds", 5)
    if isinstance(n_folds, bool) or int(n_folds) != n_folds or int(n_folds) < 2:
        raise ConfigError(f"config.n_folds must be an integer >= 2, got {n_folds!r}")
    seeds = cfg.get("seeds")
    if seeds is not None:
        if not isinstance(seeds, (list, tuple)) or not seeds:
            raise ConfigError("config.seeds must be a non-empty list")
        for seed in seeds:
            if isinstance(seed, bool) or int(seed) != seed:
                raise ConfigError(f"config.seeds must be integers, got {seed!r}")
    lrs = cfg.get("lrs") or {}
    if "pretrain" in lrs and float(lrs["pretrain"]) <= 0:
        raise ConfigError("config.lrs.pretrain must be positive")
    if "cm" in lrs:
        if not isinstance(lrs["cm"], (list, tuple)) or not lrs["cm"]:
            raise ConfigError("config.lrs.cm must be a non-empty list")
        if any(float(lr) <= 0 for lr in lrs["cm"]):
            raise ConfigError("config.lrs.cm must hold positive learning rates")


def structure_from_spec(spec: str | None, base: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Resolve ``--structure`` against the frozen structure.

    * ``None`` → the frozen structure unchanged.
    * a preset name (``hidden128`` / ``latent128`` / ``hidden256``) → that
      capacity override on top of the frozen structure.
    * a path to a JSON object → a partial structure override.
    """
    notes: dict[str, Any] = {"source": "locked", "requested": spec}
    struct = {key: (list(value) if isinstance(value, list) else value)
              for key, value in base.items()}
    if spec is None:
        validate_structure(struct)
        return struct, notes
    if spec in STRUCTURE_PRESETS:
        struct.update(STRUCTURE_PRESETS[spec])
        notes["source"] = "preset"
        notes["preset"] = spec
        validate_structure(struct)
        return struct, notes
    path = Path(spec)
    if path.exists():
        raw = load_json(path)
        for key, aliases in _STRUCTURE_ALIASES.items():
            value, used = _get_alias(raw, aliases)
            if value is None:
                continue
            if key in {"hidden", "adaptor"}:
                value = _as_int_list(value, key)
            elif key in {"n_genes", "z_dim", "c_dim", "n_dim", "drug_dim"}:
                value = int(value)
            elif key == "dropout":
                value = float(value)
            struct[key] = value
        notes["source"] = "json"
        notes["path"] = str(path.resolve())
        validate_structure(struct)
        return struct, notes
    raise ConfigError(
        f"--structure {spec!r} is neither a preset ({', '.join(sorted(STRUCTURE_PRESETS))}) "
        "nor an existing JSON file"
    )


def structure_compatible(
    a: Mapping[str, Any], b: Mapping[str, Any]
) -> tuple[bool, list[str]]:
    """Return whether two structures produce identically shaped parameters."""
    mismatches: list[str] = []
    for key in ("n_genes", "z_dim", "c_dim", "n_dim", "drug_dim"):
        if int(a.get(key, -1)) != int(b.get(key, -2)):
            mismatches.append(f"{key}: {a.get(key)!r} != {b.get(key)!r}")
    for key in ("hidden", "adaptor"):
        if [int(x) for x in a.get(key, [])] != [int(x) for x in b.get(key, [])]:
            mismatches.append(f"{key}: {a.get(key)!r} != {b.get(key)!r}")
    if float(a.get("dropout", -1.0)) != float(b.get("dropout", -2.0)):
        mismatches.append(f"dropout: {a.get('dropout')!r} != {b.get('dropout')!r}")
    return (not mismatches), mismatches


def apply_recipe_overrides(
    recipe: Mapping[str, Any], overrides: Mapping[str, Any]
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Apply explicit CLI overrides and record them; never invent one."""
    resolved = dict(recipe)
    applied: dict[str, Any] = {}
    for key, value in overrides.items():
        if value is None:
            continue
        resolved[key] = value
        applied[key] = value
    validate_recipe(resolved)
    return resolved, applied


def canonical_json(payload: Any) -> str:
    """Stable JSON encoding used for content hashes."""
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)


def sha256_payload(payload: Any) -> str:
    return hashlib.sha256(canonical_json(payload).encode("utf-8")).hexdigest()


def config_fingerprint(cfg: Mapping[str, Any]) -> str:
    """Content hash of the config blocks that define the run."""
    subset = {
        "structure": cfg.get("structure"),
        "recipe": cfg.get("recipe"),
        "metric": cfg.get("metric"),
        "loss": cfg.get("loss"),
        "optimizer": cfg.get("optimizer"),
        "dose_encoding": cfg.get("dose_encoding"),
        "prediction_noise": cfg.get("prediction_noise"),
        "train_drugs": cfg.get("train_drugs"),
        "stage": cfg.get("stage"),
        "seed": cfg.get("seed"),
        "init": cfg.get("init"),
        "fold": cfg.get("fold"),
    }
    return sha256_payload(subset)
