#!/usr/bin/env python3
"""One `tsg_view` run: validation-first native TranSiGen fit, replay, then scoring.

The worker implements the sealed Step 04 recipe with the view-specific identity
``(source, view, fold, seed)`` and the plan C4 artifact set:

``config.json``, ``checkpoint.pt``, ``learning_curve.csv``,
``validation_prediction.npy``, ``validation_metrics.json``, ``terminal.json``,
``test_predictions.npy``, ``test_per_sample.csv``, ``test_metrics.json``,
``manifest.json`` (``tsg-view-run-v1``) plus ``validation_manifest.json``, the
validation-choice seal.

Phase discipline (one process per phase, exact fresh-process replay):

1. ``train``  -- validation treated-expression MSE, ``EarlyStopping(patience=50)``,
   restore the best epoch; never reads a test row.
2. ``replay`` -- a fresh process re-derives the best checkpoint and reproduces the
   validation predictions bit-exactly; this seals the choice.
3. ``test``   -- only after ``validation-complete``; scores the held partition.
4. ``replay`` -- reproduces the test predictions bit-exactly and writes the manifest.
"""
import argparse
import json
import os
import platform
import resource
import secrets
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from inputs import Inputs  # noqa: E402
from views import SOURCES, VIEWS, bind, repo_root, sha256  # noqa: E402

RUN_SCHEMA = "tsg-view-run-v1"
VALIDATION_SCHEMA = "tsg-view-validation-v1"
METHOD = "no-dose"
METHOD_IDENTITY = "tsg-dose-no-dose-native-f32-v1"
SCENARIO = "cold-drug-dose"
TRAINING_SEEDS = (42, 43)
VALIDATION_SEED = 22617
INFERENCE_OFFSET = 100000
PATIENCE = 50
MAX_EPOCHS = 1000
METRICS = ("pcc", "spearman", "rmse", "r2")
PROCESS_TOKEN = secrets.token_hex(32)


class ScientificFailure(ValueError):
    """A numerical or exact-replay failure that cannot be resumed."""


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def read_json(path):
    return json.loads(Path(path).read_text())


def api(root):
    """Bind the shared implementations after the Step 02/04 API roots are on sys.path."""
    modules = bind(root)
    from tenv5k.learning import EarlyStopping, fixed_rng
    from tsg_dose import model as models
    from tsg_round1.protocol import score_prediction, summary
    return {"root": modules["root"], "root_modules": modules, "EarlyStopping": EarlyStopping,
            "fixed_rng": fixed_rng, "models": models, "score_prediction": score_prediction,
            "summary": summary}


def runtime_identity(root):
    """Digests of this worker and every reused shared module (replay refuses drift)."""
    root = Path(root)
    owned = [HERE / "learner.py", HERE / "inputs.py", HERE / "views/__init__.py"]
    shared = [
        root / "stage2/02_models_tenv-5k/output/api/tenv5k/model_transigen.py",
        root / "stage2/02_models_tenv-5k/output/api/tenv5k/learning.py",
        root / "stage2/02_models_tenv-5k/output/api/tenv5k/protocol.py",
        root / "stage2/02_models_tenv-5k/output/api/tenv5k/inputs.py",
        root / "stage2/04_tsg_dose/output/api/tsg_dose/model.py",
        root / "stage2/04_tsg_dose/output/api/tsg_round1/protocol.py",
    ]
    identity = {}
    for path in owned + shared:
        if not path.is_file():
            raise ScientificFailure(f"reused runtime module is missing: {path}")
        key = f"shared/{path.name}" if path in shared else path.name
        identity[key] = sha256(path)
    return identity


def environment(device):
    import torch
    return {"python": platform.python_version(), "torch": str(torch.__version__),
            "numpy": np.__version__, "device": str(device),
            "cuda": torch.version.cuda, "threads": torch.get_num_threads()}


def smoke_subset(data, rows, limit):
    """Round-robin over dose buckets, then sort: the Step 04 bounded-smoke rule."""
    if limit is None:
        return rows
    buckets = [list(rows[data.doses(rows) == dose]) for dose in sorted(set(data.doses(rows)))]
    chosen = []
    for level in range(max(map(len, buckets))):
        for bucket in buckets:
            if level < len(bucket):
                chosen.append(bucket[level])
            if len(chosen) == limit:
                return np.array(sorted(chosen), dtype=int)
    return np.array(sorted(chosen), dtype=int)


def predict(native, model, data, rows, seed, fixed_rng):
    with fixed_rng(seed, next(model.parameters()).device):
        value = np.asarray(native.predict(model, data, rows), dtype=np.float32)
    if not np.isfinite(value).all():
        raise ScientificFailure("non-finite prediction")
    return value


def validation_mse(prediction, data, rows):
    """The sealed checkpoint metric: treated-expression MSE over every validation gene."""
    return float(np.mean(np.square(
        np.asarray(prediction, dtype=np.float32)
        - np.asarray(data.expression[rows], dtype=np.float32))))


def _prepare_run(run):
    run = Path(run)
    if run.exists() and any(run.iterdir()):
        raise ValueError("a view run requires a new empty directory")
    run.mkdir(parents=True, exist_ok=True)
    return run


def train(freeze, source_path, run, *, rows=None, epochs=MAX_EPOCHS, device="cpu",
          fold=0, seed=42, repo_root_path=None, sealed=None):
    run = _prepare_run(run)
    started = time.monotonic()
    if not 1 <= int(epochs) <= MAX_EPOCHS:
        raise ValueError(f"epochs must be in 1..{MAX_EPOCHS}")
    if rows is not None and not 2 <= int(rows) <= 64:
        raise ValueError("smoke rows must be in 2..64")
    if rows is not None and int(epochs) > 3:
        raise ValueError("a bounded smoke fit is limited to 3 epochs")
    if rows is None and int(epochs) != MAX_EPOCHS:
        raise ValueError("a formal run requires the sealed 1000-epoch budget")
    if int(fold) not in range(3):
        raise ValueError("fold must be 0, 1 or 2")
    if int(seed) not in TRAINING_SEEDS:
        raise ValueError(f"training seed must be one of {TRAINING_SEEDS}")
    root = repo_root(repo_root_path)
    helper = api(root)
    data = Inputs(freeze, sealed=sealed, repo_root=root)
    config = data.config
    if config["view"] not in VIEWS or config["source"] not in SOURCES:
        raise ValueError("unknown view freeze identity")
    origin = {"source": config["source"], "view": config["view"], "fold": int(fold),
              "seed": int(seed), "scenario": SCENARIO, "method": METHOD}
    train_rows = data.rows(SCENARIO, origin["fold"], "train")
    valid_rows = data.rows(SCENARIO, origin["fold"], "validation")
    if rows is not None:
        train_rows = smoke_subset(data, train_rows, rows)
        valid_rows = smoke_subset(data, valid_rows, rows)
    if len(train_rows) < 2 or not len(valid_rows):
        raise ValueError("empty fitting partition")
    models = helper["models"]
    native = models.adapter(source_path, METHOD, None)
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    import torch
    torch.use_deterministic_algorithms(True)
    torch.set_num_threads(2)
    backend = models.configure_backend(device)
    config = {
        "schema": RUN_SCHEMA,
        "origin": origin,
        "freeze": {"manifest_sha256": sha256(Path(freeze) / "manifest.json"),
                   "view_product_sha256": data.manifest["view_product"]["sha256"],
                   "sealed_splits_sha256": data.config["sealed"]["splits_sha256"]},
        "identity": {"method": METHOD, "method_identity": METHOD_IDENTITY,
                     "native_source_sha256": models.SOURCE_DIGEST,
                     "max_epochs": int(epochs), "patience": PATIENCE,
                     "validation_seed": VALIDATION_SEED,
                     "inference_seed": origin["seed"] + INFERENCE_OFFSET,
                     "smoke_rows": None if rows is None else int(rows),
                     "device": str(device), "backend": backend,
                     "runtime_identity": runtime_identity(root),
                     "environment": environment(device)},
        "train_order": train_rows.tolist(),
        "validation_order": valid_rows.tolist(),
    }
    write_json(run / "config.json", config)
    model = native.build(origin["seed"], device)
    optimizer = native.optimizer(model)
    stopper = helper["EarlyStopping"](patience=PATIENCE, max_epochs=int(epochs))
    shuffle = np.random.default_rng(origin["seed"])
    history = []
    try:
        for epoch in range(1, int(epochs) + 1):
            try:
                loss = native.train_epoch(model, data, train_rows, optimizer, shuffle)
            except ValueError as error:
                if "non-finite native TranSiGen loss" not in str(error):
                    raise
                raise ScientificFailure(f"nonfinite native loss at epoch {epoch}") from error
            prediction = predict(native, model, data, valid_rows, VALIDATION_SEED,
                                 helper["fixed_rng"])
            mse = validation_mse(prediction, data, valid_rows)
            if not np.isfinite([loss, mse]).all():
                raise ScientificFailure(f"nonfinite training evidence at epoch {epoch}")
            improved = stopper.observe(mse, epoch)
            history.append({"epoch": epoch, "train_loss": loss, "validation_mse": mse,
                            "best_mse": stopper.best, "stale": stopper.stale,
                            "lr": optimizer.param_groups[0]["lr"]})
            if improved:
                torch.save({"state_dict": model.state_dict(), "config": config,
                            "best_epoch": epoch, "validation_mse": mse,
                            "input_manifest_sha256": config["freeze"]["manifest_sha256"]},
                           run / "checkpoint.pt")
            if stopper.stopped:
                break
        saved = torch.load(run / "checkpoint.pt", map_location=device, weights_only=False)
        model.load_state_dict(saved["state_dict"])
        validation = predict(native, model, data, valid_rows, VALIDATION_SEED,
                             helper["fixed_rng"])
        restored = validation_mse(validation, data, valid_rows)
        if restored != stopper.best:
            raise ScientificFailure(
                "restored checkpoint validation MSE differs from the best epoch")
        np.save(run / "validation_prediction.npy", validation)
        samples = helper["score_prediction"](
            data.population.iloc[valid_rows], data.expression[valid_rows].astype(np.float64),
            validation, data.control[valid_rows],
            prediction_ids=data.population["profile_id"].iloc[valid_rows].astype(str))
        write_json(run / "validation_metrics.json", helper["summary"](samples))
        pd.DataFrame(history).to_csv(run / "learning_curve.csv", index=False)
        terminal = {"status": "validation-trained-awaiting-replay",
                    "training_process_token": PROCESS_TOKEN, "training_pid": os.getpid(),
                    "smoke": rows is not None, "seed": origin["seed"],
                    "fold": origin["fold"],
                    "view": origin["view"], "source": origin["source"],
                    "device": str(device), "epochs": len(history),
                    "best_epoch": stopper.best_epoch, "validation_mse": stopper.best,
                    "elapsed_seconds": round(time.monotonic() - started, 3),
                    "peak_rss_kib": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
                    "parameter_count": sum(p.numel() for p in model.parameters()),
                    "weights_updated": True}
        write_json(run / "terminal.json", terminal)
        return terminal
    except Exception as error:
        write_json(run / "terminal.json", {
            "status": "failed", "failure_type": type(error).__name__, "error": str(error),
            "attempt_origin": "train"})
        raise


def load(freeze, run, source_path, repo_root_path=None, sealed=None):
    """Rebuild and re-verify a checkpointed run from its config alone."""
    root = repo_root(repo_root_path)
    helper = api(root)
    data = Inputs(freeze, sealed=sealed, repo_root=root)
    run = Path(run)
    config = read_json(run / "config.json")
    if config.get("schema") != RUN_SCHEMA:
        raise ValueError("run config schema mismatch")
    if sha256(Path(freeze) / "manifest.json") != config["freeze"]["manifest_sha256"]:
        raise ValueError("view freeze changed after training")
    identity = config["identity"]
    if identity["runtime_identity"] != runtime_identity(root):
        raise ValueError("the shared learner runtime changed after training")
    import torch
    saved = torch.load(run / "checkpoint.pt", map_location="cpu", weights_only=False)
    if saved["config"] != config or saved["input_manifest_sha256"] != config["freeze"][
            "manifest_sha256"]:
        raise ValueError("checkpoint configuration or input digest mismatch")
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    torch.use_deterministic_algorithms(True)
    torch.set_num_threads(2)
    if environment(identity["device"]) != identity["environment"]:
        raise ValueError("replay environment mismatch")
    for partition, key in (("train", "train_order"), ("validation", "validation_order")):
        rows = data.rows(SCENARIO, config["origin"]["fold"], partition)
        if identity["smoke_rows"] is not None:
            rows = smoke_subset(data, rows, identity["smoke_rows"])
        if rows.tolist() != config[key]:
            raise ValueError("checkpoint partition identity mismatch")
    models = helper["models"]
    if models.configure_backend(identity["device"]) != identity["backend"]:
        raise ValueError("backend policy mismatch")
    if any(value.is_floating_point() and value.dtype != torch.float32
           for value in saved["state_dict"].values()):
        raise ValueError("checkpoint precision mismatch")
    native = models.adapter(source_path, METHOD, None)
    model = native.build(config["origin"]["seed"], identity["device"])
    model.load_state_dict(saved["state_dict"])
    return data, native, model, config, helper


def replay(freeze, source_path, run, *, repo_root_path=None, sealed=None):
    """Fresh-process exact replay of every registered prediction, then the manifest."""
    run = Path(run)
    terminal = read_json(run / "terminal.json")
    if (terminal.get("training_pid") == os.getpid()
            or terminal.get("training_process_token") == PROCESS_TOKEN):
        raise ValueError("replay must run in a fresh process")
    if terminal["status"] not in ("validation-trained-awaiting-replay", "validation-complete",
                                  "test-awaiting-replay", "complete"):
        raise ValueError("cannot replay a failed run")
    if (run / "manifest.json").exists():
        verify(run)
    data, native, model, config, helper = load(freeze, run, source_path,
                                               repo_root_path=repo_root_path, sealed=sealed)
    identity = config["identity"]
    validation_path = run / "validation_prediction.npy"
    test_path = run / "test_predictions.npy"
    if not validation_path.is_file():
        raise ValueError("missing registered validation prediction")
    if terminal["status"] in ("test-awaiting-replay", "complete") and not test_path.is_file():
        raise ValueError("missing registered test prediction")
    if test_path.exists() and terminal.get("test_evaluation_process_token") == PROCESS_TOKEN:
        raise ValueError("test replay must run in a process other than the evaluation")
    replayed = {
        "validation": predict(native, model, data,
                              np.asarray(config["validation_order"], dtype=int),
                              identity["validation_seed"], helper["fixed_rng"])}
    rows = data.rows(SCENARIO, config["origin"]["fold"], "test")
    if identity["smoke_rows"] is not None:
        rows = smoke_subset(data, rows, identity["smoke_rows"])
    replayed["test"] = predict(native, model, data, rows, identity["inference_seed"],
                               helper["fixed_rng"])
    for partition, path in (("validation", validation_path), ("test", test_path)):
        if not path.exists():
            continue
        stored = np.load(path, mmap_mode="r")
        if stored.dtype != np.float32 or not np.array_equal(stored, replayed[partition]):
            raise ScientificFailure(f"exact fresh-process {partition} replay differs")
    if not (run / "validation_manifest.json").exists():
        write_json(run / "validation_manifest.json", {
            "schema": VALIDATION_SCHEMA,
            "input_manifest_sha256": config["freeze"]["manifest_sha256"],
            "parameter_count": terminal["parameter_count"],
            "validation_order": config["validation_order"],
            "outputs": {path.name: sha256(path) for path in sorted(run.iterdir())
                        if path.is_file() and path.name not in (
                            "manifest.json", "terminal.json", "validation_manifest.json")}})
    complete = test_path.exists()
    terminal.update({"status": "complete" if complete else "validation-complete",
                     "replay_pid": os.getpid(), "replay_process_token": PROCESS_TOKEN,
                     "replay": "exact fresh-process native sampled prediction"})
    write_json(run / "terminal.json", terminal)
    manifest(run)
    return terminal


def test(freeze, source_path, run, *, repo_root_path=None, sealed=None, selection=None):
    """Score the held partition only after a replayed validation choice exists."""
    run = Path(run)
    terminal = read_json(run / "terminal.json")
    if terminal["status"] != "validation-complete":
        raise ValueError("test scoring requires a replayed validation choice")
    verify(run)
    data, native, model, config, helper = load(freeze, run, source_path,
                                               repo_root_path=repo_root_path, sealed=sealed)
    identity = config["identity"]
    if identity["smoke_rows"] is None:
        if selection is None or not Path(selection).is_file():
            raise ValueError("a formal run requires a sealed validation selection receipt")
    elif selection is not None:
        raise ValueError("a bounded smoke run does not consume a selection receipt")
    rows = data.rows(SCENARIO, config["origin"]["fold"], "test")
    if identity["smoke_rows"] is not None:
        rows = smoke_subset(data, rows, identity["smoke_rows"])
    if not len(rows):
        raise ValueError("empty held partition")
    prediction = predict(native, model, data, rows, identity["inference_seed"],
                         helper["fixed_rng"])
    np.save(run / "test_predictions.npy", prediction)
    samples = helper["score_prediction"](
        data.population.iloc[rows], data.expression[rows].astype(np.float64), prediction,
        data.control[rows],
        prediction_ids=data.population["profile_id"].iloc[rows].astype(str))
    samples.to_csv(run / "test_per_sample.csv", index=False)
    write_json(run / "test_metrics.json", helper["summary"](samples))
    terminal.update({"status": "test-awaiting-replay",
                     "test_evaluation_process_token": PROCESS_TOKEN,
                     "test_selection_sha256": None if selection is None else sha256(selection),
                     "test_samples": int(len(samples))})
    write_json(run / "terminal.json", terminal)
    # A run without a manifest is expected between evaluation and test replay.
    (run / "manifest.json").unlink(missing_ok=True)
    return terminal


def manifest(run):
    run = Path(run)
    files = {path.name: sha256(path) for path in sorted(run.iterdir())
             if path.is_file() and path.name != "manifest.json"}
    write_json(run / "manifest.json", {
        "schema": RUN_SCHEMA,
        "origin": read_json(run / "config.json")["origin"],
        "freeze": read_json(run / "config.json")["freeze"],
        "outputs": files})
    return files


def verify(run):
    run = Path(run)
    registered = read_json(run / "manifest.json")
    if registered.get("schema") != RUN_SCHEMA:
        raise ValueError("run manifest schema mismatch")
    actual = {path.name for path in run.iterdir()
              if path.is_file() and path.name != "manifest.json"}
    if actual != set(registered["outputs"]):
        raise ValueError("run manifest coverage mismatch")
    for name, digest in registered["outputs"].items():
        if Path(name).name != name or sha256(run / name) != digest:
            raise ValueError(f"run artifact digest mismatch: {name}")
    return registered


def execute(freeze, source_path, run, phase, **kwargs):
    run = Path(run)
    if phase == "train":
        return train(freeze, source_path, run, **{key: kwargs[key] for key in (
            "rows", "epochs", "device", "fold", "seed", "repo_root_path", "sealed")})
    if phase == "replay":
        return replay(freeze, source_path, run, repo_root_path=kwargs["repo_root_path"],
                      sealed=kwargs["sealed"])
    if phase == "test":
        return test(freeze, source_path, run, repo_root_path=kwargs["repo_root_path"],
                    sealed=kwargs["sealed"], selection=kwargs.get("selection"))
    raise ValueError(f"unknown run phase: {phase}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("freeze", type=Path)
    parser.add_argument("source_path", type=Path, help="TranSiGen source directory")
    parser.add_argument("run", type=Path)
    parser.add_argument("phase", choices=("train", "replay", "test"))
    parser.add_argument("--rows", type=int, metavar="2..64")
    parser.add_argument("--epochs", type=int, default=MAX_EPOCHS)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--fold", type=int, default=0, choices=range(3))
    parser.add_argument("--seed", type=int, default=42, choices=TRAINING_SEEDS)
    parser.add_argument("--selection", type=Path)
    parser.add_argument("--sealed", type=Path)
    parser.add_argument("--repo-root", type=Path)
    args = parser.parse_args()
    result = execute(args.freeze, args.source_path, args.run, args.phase, rows=args.rows,
                     epochs=args.epochs, device=args.device, fold=args.fold,
                     seed=args.seed, selection=args.selection,
                     sealed=args.sealed, repo_root_path=args.repo_root)
    sys.stdout.write(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
