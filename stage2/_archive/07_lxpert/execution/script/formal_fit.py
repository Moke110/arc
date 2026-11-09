"""Freeze a Step 07 training recipe before formal development and full-data refit."""
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile
import time

import numpy as np
import torch

from autoencoder import TranscriptomeAutoencoder, reconstruction_mse
from inputs import SOURCES, sha256
from short_fit import (
    CHECKPOINT_SCHEMA, ExpressionRelease, attention_runtime, measure, probe, runtime, train_batch,
)


RECORD_SCHEMA = 1
OPTIMIZER = {
    "name": "AdamW", "lr": 0.001, "betas": [0.9, 0.999], "eps": 1e-8,
    "weight_decay": 0.01, "amsgrad": False, "maximize": False, "foreach": False,
    "capturable": False, "differentiable": False, "fused": False,
}
SCRIPT_FILES = ("lxpert.py", "inputs.py", "autoencoder.py", "short_fit.py", "formal_fit.py")


@contextmanager
def publish_directory(output):
    if output.exists():
        raise ValueError(f"output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent, prefix=f".{output.name}-") as staging:
        directory = Path(staging)
        yield directory
        directory.rename(output)


def write_record(path, record):
    path.write_text(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    digest = sha256(path)
    path.with_suffix(".sha256").write_text(digest + "\n")
    return digest


def read_record(path):
    content = path.read_bytes()
    digest = hashlib.sha256(content).hexdigest()
    if digest != path.with_suffix(".sha256").read_text().strip():
        raise ValueError(f"record digest mismatch: {path.name}")
    record = json.loads(content)
    if record["schema_version"] != RECORD_SCHEMA:
        raise ValueError(f"unsupported record schema: {path.name}")
    return record, digest


def contract(release, execution):
    return {
        "input": release.contract,
        "source_sha256": release.manifest["source_sha256"],
        "source_membership": release.manifest["sources"],
        "split_seed": release.manifest["split_seed"],
        "fixture": release.manifest["fixture"],
        "runtime": {
            **execution,
            "attention_backend": "flash_attention" if execution["device"] == "cuda" else "torch_cpu",
            "torch_version": str(torch.__version__),
            "numpy_version": np.__version__,
            "cuda_version": torch.version.cuda,
        },
        "script_sha256": {name: sha256(Path(__file__).with_name(name)) for name in SCRIPT_FILES},
    }


def freeze(args, source_pins):
    if min(args.batch_size, args.max_epochs, args.max_seconds) < 1 or args.seed < 0:
        raise ValueError("batch size, epoch and seconds ceilings must be positive and seed non-negative")
    if not args.fixture and args.batch_size != 32:
        raise ValueError("full-axis formal fitting requires the admitted batch size of 32")
    with publish_directory(args.run_dir) as directory:
        execution = runtime(args.device, args.fixture)
        release = ExpressionRelease(args.release, args.fixture, source_pins)
        plan = {
            "schema_version": RECORD_SCHEMA,
            "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
            **contract(release, execution),
            "recipe": {
                "seed": args.seed, "batch_size": args.batch_size,
                "max_epochs": args.max_epochs, "max_seconds": args.max_seconds,
                "optimizer": OPTIMIZER,
                "sampler": "NumPy PCG64 permutation of combined rows, keyed by seed and one-based epoch",
                "objective": "mean of per-profile observed-target reconstruction MSE",
                "selection": "minimum tenv development validation MSE, earliest epoch on ties",
            },
        }
        digest = write_record(directory / "plan.json", plan)
    return {**plan, "plan_sha256": digest}


def load_plan(args, release, execution):
    plan, digest = read_record(args.run_dir / "plan.json")
    for key, value in contract(release, execution).items():
        if plan[key] != value:
            raise ValueError(f"frozen plan mismatch: {key}")
    if plan["recipe"]["optimizer"] != OPTIMIZER:
        raise ValueError("frozen plan mismatch: optimizer")
    return plan, digest


def new_training(release, recipe, device):
    torch.manual_seed(recipe["seed"])
    model = TranscriptomeAutoencoder(release.manifest["gene_count"])
    initial_state = hashlib.sha256()
    for name, tensor in model.state_dict().items():
        initial_state.update(name.encode() + b"\0")
        initial_state.update(tensor.numpy().tobytes())
    model.to(device)
    options = {key: value for key, value in recipe["optimizer"].items() if key != "name"}
    options["betas"] = tuple(options["betas"])
    return model, torch.optim.AdamW(model.parameters(), **options), initial_state.hexdigest()


def rows_for_role(release, role):
    return {source: release.rows(source, role) for source in SOURCES}


def batches(rows, size):
    for start in range(0, len(rows), size):
        yield [(SOURCES[int(source)], int(row)) for source, row in rows[start:start + size]]


def combine_rows(by_source):
    return np.concatenate([
        np.column_stack((np.full(len(by_source[source]), index), by_source[source]))
        for index, source in enumerate(SOURCES)
    ])


def check_deadline(deadline):
    if time.perf_counter() >= deadline:
        raise ValueError("formal fit exceeded its frozen wall time ceiling")


def train_epoch(release, model, optimizer, by_source, recipe, epoch, device, deadline):
    rows = combine_rows(by_source)
    rng = np.random.Generator(np.random.PCG64([recipe["seed"], epoch]))
    rng.shuffle(rows)
    model.train()
    total, steps = 0.0, 0
    visits = dict.fromkeys(SOURCES, 0)
    for members in batches(rows, recipe["batch_size"]):
        check_deadline(deadline)
        expression, observed = release.batch(members, device)
        value = train_batch(model, optimizer, expression, observed, device)
        total += value * len(members)
        steps += 1
        for source, _ in members:
            visits[source] += 1
    if device == "cuda":
        torch.cuda.synchronize()
    check_deadline(deadline)
    return {
        "epoch": epoch,
        "steps": steps,
        "training_profiles": visits,
        "training_mse": total / len(rows),
        "sample_order_sha256": hashlib.sha256(rows.astype("<i8", copy=False).tobytes()).hexdigest(),
    }


@torch.inference_mode()
def validate(release, model, by_source, size, device, deadline):
    model.eval()
    scores = {}
    for source, rows in by_source.items():
        total = 0.0
        for start in range(0, len(rows), size):
            check_deadline(deadline)
            members = [(source, int(row)) for row in rows[start:start + size]]
            expression, observed = release.batch(members, device)
            with attention_runtime(device):
                losses = reconstruction_mse(model(expression), expression, observed)
            total += float(losses.sum(dtype=torch.float64))
        if not np.isfinite(total):
            raise ValueError(f"{source}: non-finite development validation MSE")
        scores[source] = total / len(rows)
    return scores


def progress(stage, record):
    print(json.dumps({"stage": stage, **record}, sort_keys=True), file=sys.stderr, flush=True)


def selected_epochs(epochs):
    return min(epochs, key=lambda row: (row["validation_mse"]["tenv"], row["epoch"]))["epoch"]


def develop(args, source_pins):
    start = time.perf_counter()
    with publish_directory(args.run_dir / "development") as directory:
        execution = runtime(args.device, args.fixture)
        release = ExpressionRelease(args.release, args.fixture, source_pins)
        plan, plan_sha = load_plan(args, release, execution)
        recipe = plan["recipe"]
        deadline = start + recipe["max_seconds"]
        model, optimizer, initial_sha = new_training(release, recipe, args.device)
        training = rows_for_role(release, "train")
        validation = rows_for_role(release, "validation")
        if any(len(rows) == 0 for roles in (training, validation) for rows in roles.values()):
            raise ValueError("each source requires development training and validation profiles")
        epochs = []
        for epoch in range(1, recipe["max_epochs"] + 1):
            epoch_start = time.perf_counter()
            record = train_epoch(release, model, optimizer, training, recipe, epoch,
                                 args.device, deadline)
            record["validation_mse"] = validate(release, model, validation, recipe["batch_size"],
                                                args.device, deadline)
            record["validation_profiles"] = {source: len(rows) for source, rows in validation.items()}
            check_deadline(deadline)
            record["elapsed_seconds"] = time.perf_counter() - epoch_start
            epochs.append(record)
            progress("development", record)
        report = {
            "schema_version": RECORD_SCHEMA, "stage": "development",
            "plan_sha256": plan_sha, "initial_model_sha256": initial_sha,
            "epochs": epochs, "selected_epochs": selected_epochs(epochs), **measure(start, args.device),
        }
        digest = write_record(directory / "record.json", report)
        check_deadline(deadline)
    return {**report, "development_sha256": digest}


def check_development(record, plan, plan_sha):
    epochs = record["epochs"]
    if (record["stage"] != "development" or record["plan_sha256"] != plan_sha
            or [row["epoch"] for row in epochs] != list(range(1, plan["recipe"]["max_epochs"] + 1))
            or record["selected_epochs"] != selected_epochs(epochs)):
        raise ValueError("development record does not match frozen selection recipe")
    for epoch in epochs:
        for role, key in (("train", "training_profiles"), ("validation", "validation_profiles")):
            counts = {source: plan["source_membership"][source]["development"][role]["profiles"]
                      for source in SOURCES}
            if epoch[key] != counts:
                raise ValueError(f"development epoch coverage mismatch: {role}")
        mse = epoch["validation_mse"]
        if set(mse) != set(SOURCES) or not all(np.isfinite(value) and value >= 0 for value in mse.values()):
            raise ValueError("invalid source-specific development MSE")


def refit(args, source_pins):
    start = time.perf_counter()
    with publish_directory(args.output) as directory:
        execution = runtime(args.device, args.fixture)
        release = ExpressionRelease(args.release, args.fixture, source_pins)
        plan, plan_sha = load_plan(args, release, execution)
        development_path = args.run_dir / "development/record.json"
        development, development_sha = read_record(development_path)
        check_development(development, plan, plan_sha)
        recipe = plan["recipe"]
        deadline = start + recipe["max_seconds"]
        model, optimizer, initial_sha = new_training(release, recipe, args.device)
        if initial_sha != development["initial_model_sha256"]:
            raise ValueError("fresh initialization differs from the frozen development seed")
        training = {source: np.arange(len(release.profiles[source]), dtype=np.int64)
                    for source in SOURCES}
        epochs = []
        for epoch in range(1, development["selected_epochs"] + 1):
            epoch_start = time.perf_counter()
            record = train_epoch(release, model, optimizer, training, recipe, epoch,
                                 args.device, deadline)
            record["elapsed_seconds"] = time.perf_counter() - epoch_start
            epochs.append(record)
            progress("full_data_refit", record)

        probe_rows = [(source, int(release.rows(source, "validation")[0])) for source in SOURCES]
        prediction, probe_report = probe(release, model, probe_rows, args.device)
        metadata = {
            "schema_version": CHECKPOINT_SCHEMA, **release.contract, **execution,
            "fit_stage": "full_data_refit",
            "plan_sha256": plan_sha, "development_sha256": development_sha,
            "selected_epochs": development["selected_epochs"], "initial_model_sha256": initial_sha,
            "training_profiles": {source: len(rows) for source, rows in training.items()},
        }
        checkpoint = directory / "checkpoint.pt"
        torch.save({"metadata": metadata, "model": model.state_dict(), "probe_rows": probe_rows,
                    "probe_prediction": prediction}, checkpoint)
        for path, name, expected in (
            (args.run_dir / "plan.json", "plan.json", plan_sha),
            (development_path, "development.json", development_sha),
        ):
            copied = directory / name
            shutil.copyfile(path, copied)
            if sha256(copied) != expected:
                raise ValueError(f"record changed during refit: {name}")
            copied.with_suffix(".sha256").write_text(expected + "\n")
        report = {
            "schema_version": RECORD_SCHEMA, "checkpoint_metadata": metadata,
            "checkpoint_sha256": sha256(checkpoint), "epochs": epochs, "probe": probe_report,
            "earlier_tenv_test_expression_in_refit": not args.fixture,
            **measure(start, args.device),
        }
        digest = write_record(directory / "refit.json", report)
        check_deadline(deadline)
    return {**metadata, **report, "refit_sha256": digest}


def verify_formal_release(directory, metadata, checkpoint_sha, release):
    plan, plan_sha = read_record(directory / "plan.json")
    development, development_sha = read_record(directory / "development.json")
    refit, refit_sha = read_record(directory / "refit.json")
    check_development(development, plan, plan_sha)
    counts = {source: len(release.profiles[source]) for source in SOURCES}
    fixture = release.manifest["fixture"]
    if (plan_sha != metadata["plan_sha256"] or development_sha != metadata["development_sha256"]
            or plan["input"] != release.contract
            or plan["fixture"] is not fixture
            or refit["earlier_tenv_test_expression_in_refit"] is not (not fixture)
            or development["selected_epochs"] != metadata["selected_epochs"]
            or development["initial_model_sha256"] != metadata["initial_model_sha256"]
            or metadata["training_profiles"] != counts or refit["checkpoint_metadata"] != metadata
            or refit["checkpoint_sha256"] != checkpoint_sha
            or [row["epoch"] for row in refit["epochs"]] != list(range(1, metadata["selected_epochs"] + 1))
            or any(row["training_profiles"] != counts for row in refit["epochs"])):
        raise ValueError("formal checkpoint provenance mismatch")
    return {"fit_stage": metadata["fit_stage"], "plan_sha256": plan_sha,
            "development_sha256": development_sha, "refit_sha256": refit_sha}
