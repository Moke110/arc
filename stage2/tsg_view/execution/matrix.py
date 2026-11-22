#!/usr/bin/env python3
"""The sealed 45-fit formal matrix: plan, run identity, paths and choice receipt.

The plan is a deterministic data structure, written once and digest-sealed **before**
the first fit. Every run is identified by ``(source, view, fold, seed)``; the schedule
keeps every tenv fit before every Tahoe fit so the report's tenv pass can publish first.

This module owns:

* :func:`plan` -- the exact 45 intended origins, in schedule order.
* :func:`frozen_plan` -- the plan plus every freeze/truth/code digest it binds, so a
  resumed driver can prove the runtime, inputs and schedule did not move.
* :func:`seal_plan`/:func:`load_plan`/:func:`check_plan` -- the plan file and its
  digest manifest, with fail-closed behaviour on any mismatch.
* :func:`selection_evidence`/:func:`lock_selection` -- the validation-locked choice
  receipt for one ``(source, fold)`` group, sealed before any test row is read.
* run-directory, registration-free run-record and resource-attempt paths.
"""
import platform
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from common import (SEEDS, SOURCE_ORDER, VIEW_ORDER, VIEW_SEEDS,  # noqa: E402
                    dumps, read_json, relative, seeds_for)
from common import views_for as _views_for  # noqa: E402
from views import SOURCES, VIEWS, repo_root, sha256, write_json  # noqa: E402

PLAN_SCHEMA = "tsg-view-formal-plan-v1"
PLAN_MANIFEST_SCHEMA = "tsg-view-formal-plan-manifest-v1"
SELECTION_SCHEMA = "tsg-view-validation-selection-v1"
RECORD_SCHEMA = "tsg-view-run-record-v1"
RUN_SCHEMA = "tsg-view-run-v1"
VALIDATION_SCHEMA = "tsg-view-validation-v1"
FREEZE_SCHEMA = "tsg-view-freeze-v1"
TRUTH_SCHEMA = "tsg-view-truth-v1"
TRUTH_SCHEMAS = (TRUTH_SCHEMA, "tsg-view-truth-fold-v1")
METHOD = "no-dose"
SCENARIO = "cold-drug-dose"
#: The schedule has one source of truth (``common``): the plan, the figure grid and the
#: report tables all order views and seeds the same way.
PATIENCE = 50
MAX_EPOCHS = 1000
VALIDATION_SEED = 22617
INFERENCE_OFFSET = 100000
CAPS = {"memory_max": "10G", "memory_swap_max": "8G", "threads": 2,
        "systemd_scope": "systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G"}

#: The run directory artifact set of plan C4 (``manifest.json`` is the digest envelope).
EXPECTED_OUTPUTS = frozenset({
    "config.json", "checkpoint.pt", "learning_curve.csv", "validation_prediction.npy",
    "validation_metrics.json", "validation_manifest.json", "terminal.json",
    "test_predictions.npy", "test_per_sample.csv", "test_metrics.json",
})

#: Sealed no-dose cold-drug-dose anchors (Step 04 axis-local round-1 summary CSV,
#: tenv lines 58-61 / tahoe lines 86-89). These are equal-fold means at seed 42; the
#: generated ``full`` arm must reproduce them within the tolerance below, and every
#: measured drift is recorded in the verify report instead of being hidden.
FULL_ANCHORS = {
    "tenv": {"pcc": 0.2654588998492864, "pcc_sd": 0.01749819417069059,
             "rmse": 0.421337591468443, "rmse_sd": 0.011933315642387388,
             "spearman": 0.2860951512686545, "r2": -0.08243900295460854},
    "tahoe_bulk": {"pcc": 0.4226328493964055, "pcc_sd": 0.025321961689764784,
                   "rmse": 0.15479625072122077, "rmse_sd": 0.011476647255689428,
                   "spearman": 0.3584354787379628, "r2": -0.28114642671533124},
}
#: Declared replication tolerance. ``origin`` admits the sealed fold-to-fold spread
#: (about 2 SD of the anchor), ``fold_mean`` is the strict matrix-level check.
FULL_TOLERANCE = {
    "origin": {"pcc_abs": 0.05, "rmse_abs": 0.02},
    "fold_mean": {"pcc_abs": 0.01, "rmse_abs": 0.005},
}
CODE_FILES = ("matrix.py", "formal.py", "learner.py", "inputs.py", "freeze.py",
              "views/__init__.py")


# --------------------------------------------------------------------------- #
# schedule
# --------------------------------------------------------------------------- #
def views_for(source):
    return _views_for(source)


def plan():
    """The exact 45 intended origins: 42 main fits plus 3 ``solo-deepest`` fits.

    Order is part of the contract: all tenv runs first (``solo-deepest`` included, so
    the tenv pass is complete before Tahoe starts), then all Tahoe runs. Within a
    source the order is ``fold``, then view, then seed so the driver's fold groups are
    contiguous and the schedule is deterministic.
    """
    origins = []
    for source in SOURCE_ORDER:
        for fold in range(3):
            for view in VIEW_ORDER[source]:
                for seed in seeds_for(view):
                    origins.append({"source": source, "view": view, "fold": fold,
                                    "seed": seed, "scenario": SCENARIO, "method": METHOD})
    return origins


def origin_key(origin):
    return (origin["source"], origin["view"], int(origin["fold"]), int(origin["seed"]))


def filter_origins(origins, *, only_source=(), only_view=(), only_fold=(), only_seed=()):
    """Select a schedule subset; an empty selection is an error, never a no-op run."""
    sources = {only_source} if isinstance(only_source, str) else set(only_source)
    views = {only_view} if isinstance(only_view, str) else set(only_view)
    folds = {only_fold} if isinstance(only_fold, int) else set(only_fold)
    seeds = {only_seed} if isinstance(only_seed, int) else set(only_seed)
    chosen = [origin for origin in origins
              if (not sources or origin["source"] in sources)
              and (not views or origin["view"] in views)
              and (not folds or int(origin["fold"]) in folds)
              and (not seeds or int(origin["seed"]) in seeds)]
    if not chosen:
        raise ValueError("no origin matches the requested filters")
    return chosen


# --------------------------------------------------------------------------- #
# paths
# --------------------------------------------------------------------------- #
def run_relative(origin):
    return Path("runs") / origin["source"] / origin["view"] / (
        f"fold-{int(origin['fold'])}") / f"seed-{int(origin['seed'])}"


def run_path(root, origin):
    return Path(root) / run_relative(origin)


def record_relative(origin):
    return Path("records") / origin["source"] / origin["view"] / (
        f"fold-{int(origin['fold'])}") / f"seed-{int(origin['seed'])}.json"


def record_path(root, origin):
    return Path(root) / record_relative(origin)


def selection_relative(source, fold):
    return Path("selection") / source / f"fold-{int(fold)}.json"


def selection_path(root, source, fold):
    return Path(root) / selection_relative(source, fold)


def freeze_path(freezes, view, source):
    return Path(freezes) / view / source


def truth_path(truth, source, fold, kind=("primary", "condition_mean")):
    return Path(truth) / source / kind / f"fold-{int(fold)}"


def resources_path(root, origin, phase, action, stamp):
    attempt = Path("resources") / run_relative(origin).relative_to("runs")
    return Path(root) / attempt / f"{phase}-{action}-{stamp}"


# --------------------------------------------------------------------------- #
# sealing
# --------------------------------------------------------------------------- #
def _bind_freezes(origins, freezes, root):
    bindings = {}
    for origin in origins:
        key = f"{origin['view']}/{origin['source']}"
        if key in bindings:
            continue
        directory = freeze_path(freezes, origin["view"], origin["source"])
        manifest_path, config_path = directory / "manifest.json", directory / "config.json"
        for path in (manifest_path, config_path):
            if not path.is_file():
                raise ValueError(f"freeze manifest is missing: {path}")
        manifest = read_json(manifest_path)
        if manifest.get("schema") != FREEZE_SCHEMA:
            raise ValueError(f"freeze schema mismatch: {manifest_path}")
        config = read_json(config_path)
        if (config.get("view"), config.get("source")) != (origin["view"],
                                                          origin["source"]):
            raise ValueError(f"freeze identity mismatch: {directory}")
        bindings[key] = {"view": origin["view"], "source": origin["source"],
                         "n_treated": int(config["n_treated"]),
                         "manifest_sha256": sha256(manifest_path),
                         "config_sha256": sha256(config_path),
                         "path": relative(directory, root)}
    return bindings


def _bind_truth(origins, truth, root):
    bindings, source_bindings = {}, {}
    for origin in origins:
        source, fold = origin["source"], int(origin["fold"])
        if source not in source_bindings:
            source_manifest = Path(truth) / source / "manifest.json"
            if not source_manifest.is_file():
                raise ValueError(f"truth source manifest is missing: {source_manifest}")
            if read_json(source_manifest).get("schema") not in TRUTH_SCHEMAS:
                raise ValueError(f"truth schema mismatch: {source_manifest}")
            source_bindings[source] = {"manifest_sha256": sha256(source_manifest),
                                       "path": relative(Path(truth) / source, root)}
        key = f"{source}/fold-{fold}"
        if key in bindings:
            continue
        entry = {}
        for kind in ("primary", "condition_mean"):
            directory = truth_path(truth, source, fold, kind)
            manifest_path = directory / "manifest.json"
            if not manifest_path.is_file():
                raise ValueError(f"{kind} truth manifest is missing: {manifest_path}")
            manifest = read_json(manifest_path)
            if manifest.get("schema") not in TRUTH_SCHEMAS:
                raise ValueError(f"{kind} truth schema mismatch: {manifest_path}")
            if (manifest.get("source"), int(manifest.get("fold", -1)), manifest.get("kind")
                    ) != (source, fold, kind):
                raise ValueError(f"{kind} truth identity mismatch: {manifest_path}")
            entry[kind] = {"manifest_sha256": sha256(manifest_path),
                           "rows": manifest.get("rows"),
                           "path": relative(directory, root)}
        bindings[key] = entry
    return bindings, source_bindings


def _code_bindings():
    bindings = {}
    for name in CODE_FILES:
        path = HERE / name
        if not path.is_file():
            raise ValueError(f"formal code file is missing: {path}")
        bindings[name] = sha256(path)
    return bindings


def _environment(device):
    import torch
    return {"python": platform.python_version(), "torch": str(torch.__version__),
            "numpy": np.__version__, "cuda": torch.version.cuda, "device": str(device)}


def frozen_plan(*, freezes, truth, source_dir, device="cuda", mode="formal",
                probe_rows=None, probe_epochs=None, only_source=(), only_view=(),
                only_fold=(), only_seed=(), repo_root_path=None):
    """Bind the schedule to every input digest it will consume.

    A formal plan always contains exactly the 45 intended origins; filters are a run
    time scheduling device and are refused here. A bounded probe plan is the filtered
    schedule itself, so its fold-group selection receipts stay well defined.
    """
    if mode not in ("formal", "probe"):
        raise ValueError("formal plan mode must be 'formal' or 'probe'")
    root = repo_root(repo_root_path)
    if mode == "formal":
        if any((only_source, only_view, only_fold, only_seed)):
            raise ValueError("the formal plan is fixed at 45 origins; filters apply to run")
        origins = plan()
        max_epochs = MAX_EPOCHS
        probe = None
    else:
        if probe_rows is None or probe_epochs is None:
            raise ValueError("a bounded probe plan requires rows and epochs")
        if not 2 <= int(probe_rows) <= 64:
            raise ValueError("probe rows must be in 2..64")
        if not 1 <= int(probe_epochs) <= 3:
            raise ValueError("probe epochs must be in 1..3")
        origins = filter_origins(plan(), only_source=only_source, only_view=only_view,
                                 only_fold=only_fold, only_seed=only_seed)
        max_epochs = int(probe_epochs)
        probe = {"rows": int(probe_rows), "epochs": int(probe_epochs)}
    source_dir = Path(source_dir)
    model_path = source_dir / "src/model.py"
    if not model_path.is_file():
        raise ValueError(f"the native TranSiGen source is unavailable: {source_dir}")
    import learner
    native = learner.api(root)["models"]
    bindings = _bind_freezes(origins, freezes, root)
    truth_bindings, truth_source = _bind_truth(origins, truth, root)
    return {
        "schema": PLAN_SCHEMA,
        "mode": mode,
        "probe": probe,
        "origins": origins,
        "scenario": SCENARIO,
        "method": METHOD,
        "seeds": list(SEEDS),
        "device": str(device),
        "max_epochs": max_epochs,
        "patience": PATIENCE,
        "validation_seed": VALIDATION_SEED,
        "inference_offset": INFERENCE_OFFSET,
        "caps": CAPS,
        "freezes": bindings,
        "truth": truth_bindings,
        "truth_source": truth_source,
        "code": _code_bindings(),
        "learner_runtime": learner.runtime_identity(root),
        "native_source": {"path": relative(source_dir, root),
                          "model_sha256": sha256(model_path),
                          "source_digest": native.SOURCE_DIGEST},
        "environment": _environment(device),
        "full_arm": {"anchors": FULL_ANCHORS, "tolerance": FULL_TOLERANCE},
    }


def _validate_plan(value):
    if not isinstance(value, dict) or value.get("schema") != PLAN_SCHEMA:
        raise ValueError("formal plan schema mismatch")
    mode = value.get("mode")
    if mode not in ("formal", "probe"):
        raise ValueError("formal plan mode is unknown")
    origins = value.get("origins")
    if not isinstance(origins, list) or not origins:
        raise ValueError("formal plan carries no origins")
    if mode == "formal" and origins != plan():
        raise ValueError("formal plan is not exactly the 45 intended origins")
    keys, seen = [], set()
    for origin in origins:
        try:
            key = origin_key(origin)
        except (KeyError, TypeError) as error:
            raise ValueError("formal plan origin is malformed") from error
        if key in seen:
            raise ValueError("formal plan repeats a run identity")
        seen.add(key)
        if (origin.get("source") not in SOURCES or origin.get("view") not in VIEWS
                or int(origin.get("fold", -1)) not in range(3)
                or int(origin.get("seed", -1)) not in SEEDS
                or origin.get("scenario") != SCENARIO or origin.get("method") != METHOD):
            raise ValueError(f"formal plan origin is outside the registered matrix: {origin}")
        keys.append(key)
    if mode == "formal" and len(keys) != 45:
        raise ValueError("formal plan does not enumerate 45 fits")
    for section in ("scenario", "method", "seeds", "device", "max_epochs", "patience",
                    "validation_seed", "inference_offset", "caps", "freezes", "truth",
                    "truth_source", "code", "learner_runtime", "native_source",
                    "environment", "full_arm"):
        if section not in value:
            raise ValueError(f"formal plan is missing its {section} identity")
    if mode == "probe" and value.get("probe") != {"rows": int(value["probe"]["rows"]),
                                                  "epochs": int(value["max_epochs"])}:
        raise ValueError("bounded probe plan identity mismatch")
    return value


def seal_plan(candidate, plan_value):
    """Write the plan and its digest seal before the first fit, or verify the seal."""
    candidate = Path(candidate)
    _validate_plan(plan_value)
    plan_path = candidate / "formal_plan.json"
    manifest_path = candidate / "formal_plan_manifest.json"
    if not plan_path.exists() and not manifest_path.exists():
        for tree in ("runs", "selection"):
            base = candidate / tree
            if base.exists() and any(base.iterdir()):
                raise ValueError(
                    "missing sealed plan for an existing formal candidate; refuse to "
                    "retro-seal a run directory")
        candidate.mkdir(parents=True, exist_ok=True)
        write_json(plan_path, plan_value)
        write_json(manifest_path, {
            "schema": PLAN_MANIFEST_SCHEMA,
            "plan_sha256": sha256(plan_path),
            "origins": len(plan_value["origins"]),
            "mode": plan_value["mode"]})
        return plan_path
    load_plan(candidate)
    if read_json(plan_path) != plan_value:
        raise ValueError("sealed formal plan differs from the frozen plan")
    return plan_path


def load_plan(candidate):
    """Read and digest-verify the sealed plan; a missing half is a hard failure."""
    candidate = Path(candidate)
    plan_path = candidate / "formal_plan.json"
    manifest_path = candidate / "formal_plan_manifest.json"
    if not plan_path.is_file() or not manifest_path.is_file():
        raise ValueError("the sealed formal plan is missing")
    manifest = read_json(manifest_path)
    if manifest.get("schema") != PLAN_MANIFEST_SCHEMA:
        raise ValueError("formal plan manifest schema mismatch")
    if manifest.get("plan_sha256") != sha256(plan_path):
        raise ValueError("formal plan digest does not match its manifest")
    value = read_json(plan_path)
    _validate_plan(value)
    if int(manifest.get("origins", -1)) != len(value["origins"]):
        raise ValueError("formal plan manifest origin count mismatch")
    return value


def plan_sha256(candidate):
    return sha256(Path(candidate) / "formal_plan.json")


def check_plan(candidate, expected):
    sealed = load_plan(candidate)
    if sealed != expected:
        raise ValueError("sealed formal plan, runtime or input identity changed")
    return sealed


# --------------------------------------------------------------------------- #
# resume machine
# --------------------------------------------------------------------------- #
def next_action(status, phase):
    transitions = {
        "validation": {None: "train",
                       "validation-trained-awaiting-replay": "replay",
                       "validation-complete": None,
                       "test-awaiting-replay": None,
                       "complete": None},
        "test": {"validation-complete": "evaluate",
                 "test-awaiting-replay": "replay",
                 "complete": None},
    }
    if phase not in transitions:
        raise ValueError(f"unknown run phase: {phase}")
    if status not in transitions[phase]:
        raise ValueError(f"cannot resume {phase} from {status}")
    return transitions[phase][status]


# --------------------------------------------------------------------------- #
# validation-locked choice receipt
# --------------------------------------------------------------------------- #
def selection_evidence(candidate, source, fold):
    """Recompute the fold's validation evidence; deterministic and test-independent."""
    candidate = Path(candidate)
    sealed = load_plan(candidate)
    if source not in SOURCES or int(fold) not in range(3):
        raise ValueError(f"unknown selection group: {source}/fold-{fold}")
    origins = [origin for origin in sealed["origins"]
               if origin["source"] == source and int(origin["fold"]) == int(fold)]
    if not origins:
        raise ValueError(f"the plan has no origin in {source}/fold-{fold}")
    entries = []
    for origin in origins:
        run = run_path(candidate, origin)
        terminal_path = run / "terminal.json"
        if not terminal_path.is_file():
            raise ValueError(f"validation is not trained: {run}")
        terminal = read_json(terminal_path)
        if terminal.get("status") not in ("validation-complete", "test-awaiting-replay",
                                          "complete"):
            raise ValueError(f"validation is not sealed for {run}: {terminal.get('status')}")
        validation_path = run / "validation_manifest.json"
        if not validation_path.is_file():
            raise ValueError(f"the fold's validation seal is missing: {validation_path}")
        validation = read_json(validation_path)
        if validation.get("schema") != VALIDATION_SCHEMA:
            raise ValueError(f"validation seal schema mismatch: {validation_path}")
        config = read_json(run / "config.json")
        if (config.get("origin") or {}).get("view") != origin["view"] or (
                config.get("origin") or {}).get("seed") != origin["seed"]:
            raise ValueError(f"run config identity mismatch: {run}")
        if (config.get("plan") or {}).get("sha256") != plan_sha256(candidate):
            raise ValueError(f"run is bound to another formal plan: {run}")
        binding = sealed["freezes"].get(f"{origin['view']}/{source}")
        if binding is None or (config.get("freeze") or {}).get(
                "manifest_sha256") != binding["manifest_sha256"]:
            raise ValueError(f"run is bound to another freeze: {run}")
        entries.append({
            "view": origin["view"],
            "seed": int(origin["seed"]),
            "validation_manifest_sha256": sha256(validation_path),
            "validation_mse": terminal.get("validation_mse"),
            "best_epoch": terminal.get("best_epoch"),
            "epochs": terminal.get("epochs"),
            "parameter_count": terminal.get("parameter_count"),
        })
    return {"schema": SELECTION_SCHEMA, "source": source, "fold": int(fold),
            "plan_sha256": plan_sha256(candidate),
            "truth": sealed["truth"][f"{source}/fold-{int(fold)}"],
            "origins": entries}


def write_selection(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(dumps(value) + "\n")
    temporary.replace(path)
    return path


def lock_selection(candidate, source, fold):
    """Seal one fold's choice receipt before any test row is read; idempotent-or-fail."""
    path = selection_path(candidate, source, fold)
    evidence = selection_evidence(candidate, source, fold)
    if path.is_file():
        existing = read_json(path)
        if existing != evidence:
            raise ValueError("sealed validation selection receipt changed; refuse to "
                             "re-open the test partition")
        return path
    return write_selection(path, evidence)


# --------------------------------------------------------------------------- #
# per-run record
# --------------------------------------------------------------------------- #
def run_record(candidate, origin, phases=()):
    """The C4 per-run record: seed, epochs, best epoch, hashes, device, cost."""
    candidate = Path(candidate)
    run = run_path(candidate, origin)
    config = read_json(run / "config.json")
    terminal = read_json(run / "terminal.json")
    manifest = read_json(run / "manifest.json")
    metrics = read_json(run / "test_metrics.json")
    peak_rss = [int(terminal.get("peak_rss_kib") or 0)]
    peak_gpu = [0]
    runtime = float(terminal.get("elapsed_seconds") or 0.0)
    for phase in phases:
        runtime += float(phase.get("elapsed_seconds") or 0.0)
        peak_rss.append(int(phase.get("peak_rss_kib") or 0))
        peak_gpu.append(int(phase.get("peak_gpu_used_mib") or 0))
    return {
        "schema": RECORD_SCHEMA,
        "origin": config["origin"],
        "plan_sha256": (config.get("plan") or {}).get("sha256"),
        "mode": (config.get("identity") or {}).get("mode"),
        "seed": terminal.get("seed"),
        "fold": terminal.get("fold"),
        "device": terminal.get("device"),
        "realized_epochs": terminal.get("epochs"),
        "best_epoch": terminal.get("best_epoch"),
        "validation_mse": terminal.get("validation_mse"),
        "runtime_seconds": round(runtime, 3),
        "training_seconds": terminal.get("elapsed_seconds"),
        "peak_rss_kib": max(peak_rss),
        "peak_gpu_used_mib": max(peak_gpu),
        "truth": terminal.get("truth"),
        "selection_sha256": terminal.get("test_selection_sha256"),
        "artifacts": manifest.get("outputs"),
        "metrics": {"primary": metrics.get("primary"),
                    "condition_mean": metrics.get("condition_mean")},
        "phases": list(phases),
    }


def write_record(candidate, origin, record=None):
    path = record_path(candidate, origin)
    value = record if record is not None else run_record(candidate, origin)
    path.parent.mkdir(parents=True, exist_ok=True)
    write_json(path, value)
    return path


# --------------------------------------------------------------------------- #
# full-arm tolerance
# --------------------------------------------------------------------------- #
def compare_full_arm(source, metrics, *, level="origin"):
    """Compare measured ``full`` scores with the sealed anchor and record the drift."""
    if source not in FULL_ANCHORS:
        raise ValueError(f"no sealed no-dose anchor for source: {source}")
    anchor = FULL_ANCHORS[source]
    tolerance = FULL_TOLERANCE.get(level)
    if tolerance is None:
        raise ValueError(f"unknown full-arm tolerance level: {level}")
    measured = {"pcc": float(metrics["pcc"]), "rmse": float(metrics["rmse"])}
    drift = {"pcc": measured["pcc"] - anchor["pcc"],
             "rmse": measured["rmse"] - anchor["rmse"]}
    within = (abs(drift["pcc"]) <= tolerance["pcc_abs"]
              and abs(drift["rmse"]) <= tolerance["rmse_abs"])
    return {"status": "PASS" if within else "DRIFT", "level": level,
            "source": source, "anchor": anchor, "measured": measured,
            "drift": drift, "tolerance": tolerance}
