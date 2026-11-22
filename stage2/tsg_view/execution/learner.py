#!/usr/bin/env python3
"""One `tsg_view` run: validation-first native TranSiGen fit, replay, then scoring.

The worker implements the sealed Step 04 recipe with the view-specific identity
``(source, view, fold, seed)`` and the plan C4 artifact set:

``config.json``, ``checkpoint.pt``, ``learning_curve.csv``,
``validation_prediction.npy``, ``validation_metrics.json``, ``terminal.json``,
``test_predictions.npy``, ``test_per_sample.csv``, ``test_metrics.json``,
``manifest.json`` (``tsg-view-run-v1``) plus ``validation_manifest.json``, the
validation-choice per-run seal.

Three modes:

* ``smoke``  -- bounded rows/epochs; evaluates on the view's own test partition (the
  ticket 01-04 smoke contract). No plan and no truth product are involved.
* ``probe``  -- bounded rows/epochs **through the formal path**: bound to a sealed probe
  plan, validation-locked by a fold choice receipt, and evaluated only against the
  frozen truth product. The two bounded GPU probes use this mode.
* ``formal`` -- the sealed 1000-epoch budget, bound to the 45-fit formal plan, evaluated
  only against the frozen truth product.

Phase discipline (one process per phase, exact fresh-process replay):

1. ``train``  -- validation treated-expression MSE, ``EarlyStopping(patience=50)``,
   restore the best epoch; never reads a test row (in particular never reads the truth
   product).
2. ``replay`` -- a fresh process re-derives the best checkpoint and reproduces the
   validation predictions bit-exactly; this seals the choice (and re-reproduces test
   predictions once they exist).
3. ``test``   -- only after ``validation-complete`` and after the fold's choice receipt
   exists and matches a recomputation; scores the frozen primary and condition-mean
   truth, never the view's test partition.
4. ``replay`` -- reproduces the test predictions bit-exactly and writes the manifest.

Evaluation is delta-space (plan C7.10): the model reads only ``truth.control`` and the
truth product's copied Morgan fingerprints, predicts expression, and scoring derives
``Δ = prediction - control``. Expression is never reconstructed as ``Δ + control``.
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
TRUTH_SCHEMA = "tsg-view-truth-v1"
TRUTH_SCHEMAS = (TRUTH_SCHEMA, "tsg-view-truth-fold-v1")
METHOD = "no-dose"
METHOD_IDENTITY = "tsg-dose-no-dose-native-f32-v1"
SCENARIO = "cold-drug-dose"
TRAINING_SEEDS = (42, 43)
VALIDATION_SEED = 22617
INFERENCE_OFFSET = 100000
PATIENCE = 50
MAX_EPOCHS = 1000
MODES = ("smoke", "probe", "formal")
METRICS = ("pcc", "spearman", "rmse", "r2")
PROCESS_TOKEN = secrets.token_hex(32)
TRUTH_KINDS = ("primary", "condition_mean")


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
    from tenv5k.protocol import vector_metrics
    from tsg_dose import model as models
    from tsg_round1.protocol import score_prediction, summary
    return {"root": modules["root"], "root_modules": modules, "EarlyStopping": EarlyStopping,
            "fixed_rng": fixed_rng, "models": models, "score_prediction": score_prediction,
            "summary": summary, "vector_metrics": vector_metrics}


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


# --------------------------------------------------------------------------- #
# frozen truth product (view-independent, plan C3) and its inference shim
# --------------------------------------------------------------------------- #
class TruthBlock:
    """One verified frozen-truth fold directory (primary or condition-mean)."""

    def __init__(self, directory, source, fold, kind):
        self.directory = Path(directory)
        self.source, self.fold, self.kind = source, int(fold), kind
        manifest_path = self.directory / "manifest.json"
        if not manifest_path.is_file():
            raise ScientificFailure(f"truth manifest is missing: {manifest_path}")
        self.manifest = read_json(manifest_path)
        if self.manifest.get("schema") not in TRUTH_SCHEMAS:
            raise ScientificFailure(f"truth manifest schema mismatch: {manifest_path}")
        if (self.manifest.get("source"), self.manifest.get("kind"),
                int(self.manifest.get("fold", -1))) != (source, kind, self.fold):
            raise ScientificFailure(f"truth identity mismatch: {manifest_path}")
        self.manifest_sha256 = sha256(manifest_path)
        for name, declared in (self.manifest.get("outputs") or {}).items():
            digest = declared.get("sha256") if isinstance(declared, dict) else declared
            if sha256(self.directory / name) != digest:
                raise ScientificFailure(f"truth artifact digest mismatch: {name}")
        self.population = pd.read_csv(self.directory / "population.csv")
        if "profile_id" not in self.population.columns:
            raise ScientificFailure(f"truth population misses profile_id: {self.directory}")
        parent_column = "drug_parent" if "drug_parent" in self.population.columns else "parent"
        if parent_column not in self.population.columns:
            raise ScientificFailure(f"truth population misses its parent column: {self.directory}")
        self.population = self.population.rename(columns={parent_column: "drug_parent"})
        if self.population["profile_id"].astype(str).duplicated().any():
            raise ScientificFailure(f"truth population repeats a profile_id: {self.directory}")
        for name in ("cell_line", "dose_um"):
            if name not in self.population.columns:
                raise ScientificFailure(f"truth population misses {name}: {self.directory}")
        for name, expected in (("fold", self.fold), ("scenario", SCENARIO)):
            if name in self.population.columns and not self.population[name].eq(
                    expected).all():
                raise ScientificFailure(f"truth population {name} mismatch: {self.directory}")
        if "partition" in self.population.columns and not self.population["partition"].eq(
                "test").all():
            raise ScientificFailure(f"truth population is not the test partition: {self.directory}")
        self.delta = np.load(self.directory / "delta.npy", mmap_mode="r")
        shape = (len(self.population), 5000)
        if self.delta.shape != shape or self.delta.dtype != np.float32:
            raise ScientificFailure(f"truth delta shape or dtype mismatch: {self.directory}")
        self.control = None
        self.parents = None
        self.morgan = None
        if kind == "primary":
            control_path = self.directory / "control.npy"
            if not control_path.is_file():
                raise ScientificFailure(f"primary truth control is missing: {control_path}")
            self.control = np.load(control_path, mmap_mode="r")
            if self.control.shape != shape or self.control.dtype != np.float32:
                raise ScientificFailure(f"truth control shape or dtype mismatch: {self.directory}")
            parents_path, morgan_path = (self.directory / "parents.csv",
                                         self.directory / "morgan.npy")
            if not parents_path.is_file() or not morgan_path.is_file():
                raise ScientificFailure(
                    f"primary truth fingerprints are missing: {self.directory}")
            self.parents = pd.read_csv(parents_path)
            self.morgan = np.load(morgan_path, mmap_mode="r")
            if (self.morgan.shape != (len(self.parents), 2048)
                    or self.morgan.dtype != np.uint8
                    or "drug_parent_inchikey" not in self.parents.columns):
                raise ScientificFailure(f"truth Morgan identity mismatch: {self.directory}")
            if not set(self.population["drug_parent"]) <= set(
                    self.parents["drug_parent_inchikey"]):
                raise ScientificFailure(f"truth population parent is not in parents.csv: "
                                        f"{self.directory}")
            self._positions = self.population["drug_parent"].map(
                dict(zip(self.parents["drug_parent_inchikey"],
                         range(len(self.parents))))).to_numpy(dtype=int)
        elif "members" in self.population.columns:
            self.members = self.population["members"].astype(str).map(_split_members)
        else:
            for alternate in ("member_profile_ids", "member_ids"):
                if alternate in self.population.columns:
                    self.members = self.population[alternate].astype(str).map(_split_members)
                    break
            else:
                raise ScientificFailure(
                    f"condition-mean truth misses its members column: {self.directory}")
        if kind == "condition_mean" and "count" in self.population.columns and not (
                self.population["count"].astype(int).eq(
                    [len(items) for items in self.members]).all()):
            raise ScientificFailure(f"condition-mean truth count mismatch: {self.directory}")

    @property
    def indices(self):
        return np.arange(len(self.population))

    def identities(self):
        columns = ["profile_id", "drug_parent", "dose_um", "cell_line"]
        return self.population[columns].copy()

    def inference_rows(self):
        if self.kind != "primary":
            raise ScientificFailure("only the primary truth carries control inputs")
        return InferenceRows(self.control, self.morgan, self._positions)


class InferenceRows:
    """The frozen-truth `InferenceRows` shim: ``control`` plus truth fingerprints.

    ``no-dose`` inference reads only these two inputs, so any view's model can be
    scored against the view-independent truth rows without touching the view's
    expression or test partition.
    """

    def __init__(self, control, morgan, positions):
        self.control = np.asarray(control, dtype=np.float32)
        self._morgan = morgan
        self._positions = np.asarray(positions, dtype=int)

    def __len__(self):
        return len(self.control)

    def fingerprints(self, rows):
        rows = np.asarray(rows, dtype=int)
        return np.asarray(self._morgan[self._positions[rows]], dtype=np.float32)


def _split_members(value):
    return [member for member in str(value).split(";") if member]


def truth_directory(truth_root, source, fold, kind):
    return Path(truth_root) / source / kind / f"fold-{int(fold)}"


def load_truth(truth_root, source, fold, kind):
    return TruthBlock(truth_directory(truth_root, source, fold, kind), source, fold, kind)


def check_truth_source(truth_root, source, expected_sha256):
    """Re-derive the source-wide truth manifest digest locked into the run config."""
    path = Path(truth_root) / source / "manifest.json"
    if not path.is_file():
        raise ScientificFailure(f"truth source manifest is missing: {path}")
    if sha256(path) != expected_sha256:
        raise ScientificFailure("truth source manifest changed after training")


def _samples_frame(identities, rows):
    """Join identity columns to metric rows, keeping summary's columns when empty."""
    columns = list(METRICS) + ["constant_truth", "constant_prediction", "nonfinite"]
    values = pd.DataFrame(rows, columns=columns) if rows else pd.DataFrame(
        {column: [] for column in columns})
    return identities.reset_index(drop=True).join(values)


def score_delta(obs, truth_delta, prediction, control, metrics):
    """Per-sample Δ metrics with Δ read directly from the frozen truth (never rebuilt)."""
    truth = np.asarray(truth_delta, dtype=np.float64)
    predicted = np.asarray(prediction, dtype=np.float64)
    matched = np.asarray(control, dtype=np.float64)
    if truth.shape != predicted.shape or truth.shape != matched.shape or len(obs) != len(truth):
        raise ScientificFailure("prediction identity/shape mismatch against the truth")
    predicted_delta = predicted - matched
    rows = []
    for measured, guess in zip(truth, predicted_delta):
        finite = bool(np.isfinite(measured).all() and np.isfinite(guess).all())
        if finite:
            values = metrics(measured, guess)
        else:
            values = dict.fromkeys(METRICS, float("nan"))
            values.update(constant_truth=False, constant_prediction=False)
        rows.append(dict(values, nonfinite=not finite))
    identities = pd.DataFrame({
        "profile_id": obs["profile_id"].astype(str).to_numpy(),
        "drug_parent": obs["drug_parent"].to_numpy(),
        "dose_um": np.asarray(obs["dose_um"], dtype=float),
        "cell_line": obs["cell_line"].astype(str).to_numpy()})
    return _samples_frame(identities, rows)


def score_condition_mean(condition_truth, predicted_delta, member_positions, metrics):
    """Per-condition metrics after averaging predicted Δ over the condition's members."""
    truth = np.asarray(condition_truth.delta, dtype=np.float64)
    rows = []
    for position, members in enumerate(condition_truth.members):
        try:
            average = np.asarray(
                predicted_delta[[member_positions[str(member)] for member in members]],
                dtype=np.float64).mean(axis=0)
        except KeyError as error:
            raise ScientificFailure(
                f"condition-mean truth member is absent from the primary truth: {error}") from error
        finite = bool(np.isfinite(truth[position]).all() and np.isfinite(average).all())
        if finite:
            values = metrics(truth[position], average)
        else:
            values = dict.fromkeys(METRICS, float("nan"))
            values.update(constant_truth=False, constant_prediction=False)
        rows.append(dict(values, nonfinite=not finite))
    identities = condition_truth.population[
        ["profile_id", "drug_parent", "dose_um", "cell_line", "count"]].copy()
    identities["profile_id"] = identities["profile_id"].astype(str)
    return _samples_frame(identities, rows)


def partition_rows(data, scenario, fold, partition, smoke_rows):
    rows = data.rows(scenario, fold, partition)
    return smoke_subset(data, rows, smoke_rows) if smoke_rows is not None else rows


# --------------------------------------------------------------------------- #
# training
# --------------------------------------------------------------------------- #
def _prepare_run(run):
    run = Path(run)
    if run.exists() and any(run.iterdir()):
        raise ValueError("a view run requires a new empty directory")
    run.mkdir(parents=True, exist_ok=True)
    return run


def _mode_bounds(mode, rows, epochs):
    if mode not in MODES:
        raise ValueError(f"run mode must be one of {MODES}")
    if mode == "formal":
        if rows is not None:
            raise ValueError("a formal run has no bounded row subset")
        if int(epochs) != MAX_EPOCHS:
            raise ValueError("a formal run requires the sealed 1000-epoch budget")
        return
    if rows is None or not 2 <= int(rows) <= 64:
        raise ValueError("bounded rows must be in 2..64")
    if not 1 <= int(epochs) <= 3:
        raise ValueError("a bounded run is limited to 3 epochs")
    if mode == "smoke" and int(epochs) > 3:
        raise ValueError("a bounded smoke fit is limited to 3 epochs")


def sealed_plan(matrix_root, origin, mode, epochs, rows, freeze_sha256):
    """Prove the run identity is in the sealed plan; returns the plan document."""
    import matrix
    if matrix_root is None:
        raise ValueError(f"a {mode} run requires --matrix-root")
    plan = matrix.load_plan(matrix_root)
    if plan.get("mode") != mode:
        raise ValueError(f"the sealed plan is a {plan.get('mode')} plan, not {mode}")
    matching = [item for item in plan["origins"] if matrix.origin_key(item) == (
        origin["source"], origin["view"], origin["fold"], origin["seed"])]
    if len(matching) != 1 or matching[0] != origin:
        raise ValueError("the run identity is not in the sealed plan")
    binding = plan["freezes"].get(f"{origin['view']}/{origin['source']}")
    if binding is None or binding["manifest_sha256"] != freeze_sha256:
        raise ValueError("the view freeze is not the sealed plan's freeze")
    if mode == "formal":
        if int(epochs) != int(plan["max_epochs"]):
            raise ValueError("formal epochs differ from the sealed plan")
    elif plan.get("probe") != {"rows": int(rows), "epochs": int(epochs)}:
        raise ValueError("bounded probe identity differs from the sealed plan")
    return plan


def train(freeze, source_path, run, *, rows=None, epochs=MAX_EPOCHS, device="cpu",
          fold=0, seed=42, repo_root_path=None, sealed=None, mode=None,
          matrix_root=None):
    run = _prepare_run(run)
    started = time.monotonic()
    if mode is None:
        mode = "smoke" if rows is not None else "formal"
    _mode_bounds(mode, rows, epochs)
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
    freeze_sha256 = sha256(Path(freeze) / "manifest.json")
    plan_identity = None
    if mode != "smoke":
        plan_document = sealed_plan(matrix_root, origin, mode, epochs, rows, freeze_sha256)
        plan_identity = {
            "sha256": sha256(Path(matrix_root) / "formal_plan.json"),
            "mode": mode,
            "truth": plan_document["truth"][f"{origin['source']}/fold-{origin['fold']}"],
            "truth_source_sha256": plan_document["truth_source"][origin["source"]][
                "manifest_sha256"],
            "full_arm": plan_document["full_arm"],
        }
    train_rows = partition_rows(data, SCENARIO, origin["fold"], "train", rows)
    valid_rows = partition_rows(data, SCENARIO, origin["fold"], "validation", rows)
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
        "mode": mode,
        "plan": plan_identity,
        "freeze": {"manifest_sha256": freeze_sha256,
                   "view_product_sha256": data.manifest["view_product"]["sha256"],
                   "sealed_splits_sha256": data.config["sealed"]["splits_sha256"]},
        "truth": plan_identity["truth"] if plan_identity else None,
        "identity": {"method": METHOD, "method_identity": METHOD_IDENTITY,
                     "native_source_sha256": models.SOURCE_DIGEST,
                     "max_epochs": int(epochs), "patience": PATIENCE,
                     "validation_seed": VALIDATION_SEED,
                     "inference_seed": origin["seed"] + INFERENCE_OFFSET,
                     "smoke_rows": None if rows is None else int(rows),
                     "mode": mode, "plan_sha256": plan_identity["sha256"]
                     if plan_identity else None,
                     "device": str(device), "backend": backend,
                     "threads": torch.get_num_threads(),
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
                            "input_manifest_sha256": freeze_sha256},
                           run / "checkpoint.pt")
            if stopper.stopped:
                break
        saved = torch.load(run / "checkpoint.pt", map_location=device, weights_only=False)
        model.load_state_dict(saved["state_dict"])
        if saved["best_epoch"] != stopper.best_epoch:
            raise ScientificFailure("checkpoint best epoch differs from the stopping record")
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
                    "smoke": mode == "smoke", "mode": mode,
                    "seed": origin["seed"], "fold": origin["fold"],
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


# --------------------------------------------------------------------------- #
# verification / replay / scoring
# --------------------------------------------------------------------------- #
def load(freeze, run, source_path, repo_root_path=None, sealed=None, matrix_root=None):
    """Rebuild and re-verify a checkpointed run from its config alone."""
    root = repo_root(repo_root_path)
    helper = api(root)
    data = Inputs(freeze, sealed=sealed, repo_root=root)
    run = Path(run)
    config = read_json(run / "config.json")
    if config.get("schema") != RUN_SCHEMA:
        raise ValueError("run config schema mismatch")
    mode = config.get("mode")
    if mode not in MODES:
        raise ValueError("run mode mismatch")
    if sha256(Path(freeze) / "manifest.json") != config["freeze"]["manifest_sha256"]:
        raise ValueError("view freeze changed after training")
    if mode != "smoke":
        sealed_plan(matrix_root, config["origin"], mode,
                    config["identity"]["max_epochs"], config["identity"]["smoke_rows"],
                    config["freeze"]["manifest_sha256"])
        if (config.get("plan") or {}).get("sha256") != sha256(
                Path(matrix_root) / "formal_plan.json"):
            raise ValueError("run is bound to another sealed plan")
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
        rows = partition_rows(data, SCENARIO, config["origin"]["fold"], partition,
                              identity["smoke_rows"])
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


def _truth_rows(config, truth_root):
    """Load and digest-verify the frozen truth rows the run was sealed against."""
    truth = config.get("truth") or {}
    if truth_root is None:
        raise ValueError(f"a {config.get('mode')} run requires the frozen truth product")
    check_truth_source(truth_root, config["origin"]["source"],
                       config["plan"]["truth_source_sha256"])
    primary = load_truth(truth_root, config["origin"]["source"], config["origin"]["fold"],
                         "primary")
    if primary.manifest_sha256 != truth["primary"]["manifest_sha256"]:
        raise ScientificFailure("primary truth changed after training")
    condition = load_truth(truth_root, config["origin"]["source"],
                           config["origin"]["fold"], "condition_mean")
    if condition.manifest_sha256 != truth["condition_mean"]["manifest_sha256"]:
        raise ScientificFailure("condition-mean truth changed after training")
    return primary, condition


def _check_selection(config, selection, matrix_root):
    """Require, before any test read, that the fold receipt is present and unchanged."""
    import matrix
    if selection is None or not Path(selection).is_file():
        raise ValueError("a formal run requires a sealed validation selection receipt")
    receipt = read_json(selection)
    expected = matrix.selection_evidence(matrix_root, config["origin"]["source"],
                                         config["origin"]["fold"])
    if receipt != expected:
        raise ValueError("sealed validation selection receipt changed; refuse to read the "
                         "frozen test truth")
    return sha256(Path(selection))


def _evaluate_truth(config, truth_root, prediction, helper, primary=None, condition=None):
    """Primary and condition-mean metrics in delta space against the frozen truth."""
    if primary is None or condition is None:
        primary, condition = _truth_rows(config, truth_root)
    predicted = np.asarray(prediction, dtype=np.float64)
    matched = np.asarray(primary.control, dtype=np.float64)
    predicted_delta = predicted - matched
    samples = score_delta(primary.identities(), primary.delta, prediction,
                          primary.control, helper["vector_metrics"])
    positions = {str(profile): index for index, profile in
                 enumerate(primary.population["profile_id"].astype(str))}
    condition_samples = score_condition_mean(condition, predicted_delta, positions,
                                             helper["vector_metrics"])
    return samples, condition_samples, primary, condition


def replay(freeze, source_path, run, *, repo_root_path=None, sealed=None, truth_root=None,
           matrix_root=None, selection=None):
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
                                               repo_root_path=repo_root_path,
                                               sealed=sealed, matrix_root=matrix_root)
    identity = config["identity"]
    mode = config["mode"]
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
    if test_path.exists():
        if mode != "smoke":
            _check_selection(config, selection, matrix_root)
            primary, _ = _truth_rows(config, truth_root)
            test_rows = primary.indices
            inference = primary.inference_rows()
        else:
            inference = data
            test_rows = partition_rows(data, SCENARIO, config["origin"]["fold"], "test",
                                       identity["smoke_rows"])
        replayed["test"] = predict(native, model, inference, test_rows,
                                   identity["inference_seed"], helper["fixed_rng"])
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
    if complete:
        terminal["test_replay_token"] = PROCESS_TOKEN
    else:
        terminal["validation_replay_token"] = PROCESS_TOKEN
    terminal.update({"status": "complete" if complete else "validation-complete",
                     "replay_pid": os.getpid(), "replay_process_token": PROCESS_TOKEN,
                     "replay": "exact fresh-process native sampled prediction"})
    write_json(run / "terminal.json", terminal)
    manifest(run)
    return terminal


def test(freeze, source_path, run, *, repo_root_path=None, sealed=None, selection=None,
         truth_root=None, matrix_root=None):
    """Score the frozen truth only after a replayed validation choice exists."""
    run = Path(run)
    terminal = read_json(run / "terminal.json")
    if terminal["status"] != "validation-complete":
        raise ValueError("test scoring requires a replayed validation choice")
    config_preview = read_json(run / "config.json")
    mode = config_preview.get("mode")
    verify(run)
    data, native, model, config, helper = load(freeze, run, source_path,
                                               repo_root_path=repo_root_path,
                                               sealed=sealed, matrix_root=matrix_root)
    identity = config["identity"]
    if mode == "smoke":
        if selection is not None:
            raise ValueError("a bounded smoke run does not consume a selection receipt")
        rows = partition_rows(data, SCENARIO, config["origin"]["fold"], "test",
                              identity["smoke_rows"])
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
                         "test_selection_sha256": None,
                         "test_samples": int(len(samples))})
    else:
        selection_sha256 = _check_selection(config, selection, matrix_root)
        primary, condition = _truth_rows(config, truth_root)
        if not len(primary.indices):
            raise ValueError("empty frozen truth partition")
        inference = primary.inference_rows()
        prediction = predict(native, model, inference, primary.indices,
                             identity["inference_seed"], helper["fixed_rng"])
        np.save(run / "test_predictions.npy", prediction)
        samples, condition_samples, primary, condition = _evaluate_truth(
            config, truth_root, prediction, helper, primary, condition)
        samples.to_csv(run / "test_per_sample.csv", index=False)
        primary_summary = helper["summary"](samples)
        metrics = dict(primary_summary)
        metrics.update({
            "primary": primary_summary,
            "condition_mean": helper["summary"](condition_samples),
            "truth": config["truth"],
            "selection_sha256": selection_sha256,
            "delta_space": "prediction - truth.control; truth delta read from the "
                           "frozen product (never rebuilt)",
            "n_primary_samples": int(len(samples)),
            "n_condition_mean_samples": int(len(condition_samples))})
        write_json(run / "test_metrics.json", metrics)
        terminal.update({"status": "test-awaiting-replay",
                         "test_evaluation_process_token": PROCESS_TOKEN,
                         "test_selection_sha256": selection_sha256,
                         "truth": config["truth"],
                         "test_samples": int(len(samples)),
                         "n_primary_samples": int(len(samples)),
                         "n_condition_mean_samples": int(len(condition_samples))})
    write_json(run / "terminal.json", terminal)
    # A run without a manifest is expected between evaluation and test replay.
    (run / "manifest.json").unlink(missing_ok=True)
    return terminal


def manifest(run):
    run = Path(run)
    files = {path.name: sha256(path) for path in sorted(run.iterdir())
             if path.is_file() and path.name != "manifest.json"}
    config = read_json(run / "config.json")
    value = {
        "schema": RUN_SCHEMA,
        "origin": config["origin"],
        "mode": config.get("mode"),
        "plan": config.get("plan") and config["plan"]["sha256"],
        "freeze": config["freeze"],
        "truth": config.get("truth"),
        "outputs": files}
    write_json(run / "manifest.json", value)
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
            "rows", "epochs", "device", "fold", "seed", "repo_root_path", "sealed",
            "mode", "matrix_root")})
    if phase == "replay":
        return replay(freeze, source_path, run, repo_root_path=kwargs["repo_root_path"],
                      sealed=kwargs["sealed"], truth_root=kwargs.get("truth_root"),
                      matrix_root=kwargs.get("matrix_root"),
                      selection=kwargs.get("selection"))
    if phase == "test":
        return test(freeze, source_path, run, repo_root_path=kwargs["repo_root_path"],
                    sealed=kwargs["sealed"], selection=kwargs.get("selection"),
                    truth_root=kwargs.get("truth_root"),
                    matrix_root=kwargs.get("matrix_root"))
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
    parser.add_argument("--mode", choices=MODES,
                        help="smoke (bounded, view test rows), probe (bounded, frozen "
                             "truth + choice receipt) or formal (sealed 1000 epochs)")
    parser.add_argument("--matrix-root", type=Path,
                        help="the sealed formal/probe candidate root")
    parser.add_argument("--truth", type=Path,
                        help="the frozen truth product root (output/truth)")
    parser.add_argument("--selection", type=Path)
    parser.add_argument("--sealed", type=Path)
    parser.add_argument("--repo-root", type=Path)
    args = parser.parse_args()
    mode = args.mode or ("smoke" if args.rows is not None else "formal")
    result = execute(args.freeze, args.source_path, args.run, args.phase, rows=args.rows,
                     epochs=args.epochs, device=args.device, fold=args.fold,
                     seed=args.seed, selection=args.selection,
                     sealed=args.sealed, repo_root_path=args.repo_root, mode=mode,
                     matrix_root=args.matrix_root, truth_root=args.truth)
    sys.stdout.write(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
