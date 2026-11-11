"""Fit, predict and replay Step 09 XPert from predictor-only inputs."""

import argparse
from contextlib import contextmanager
import json
import math
from pathlib import Path
import resource
import tempfile
import time

import numpy as np
import torch
from torch.nn.attention import SDPBackend, sdpa_kernel

from model import LOSS_WEIGHTS, MLP_WIDTHS, XPert, loss_terms
from predictor import AXIS_MODELS, Predictor, Profiles, sha256


PHASES = ("tahoe", "tenv")
KERNELS = {"flash": SDPBackend.FLASH_ATTENTION, "efficient": SDPBackend.EFFICIENT_ATTENTION,
           "math": SDPBackend.MATH}
RUN_SCHEMA = "step09-xpert-v1"
RERUN_CONTRACT = "step09-round1-rerun-v1"


def write_json(path, payload):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def save_checkpoint(path, payload):
    temporary = path.with_suffix(".tmp")
    torch.save(payload, temporary)
    temporary.replace(path)


def load_policy(path, fixture):
    policy = json.loads(path.read_text())
    rerun = policy.get("training_contract") == RERUN_CONTRACT
    fields = {"kind", "device", "kernel", "dtype", "microbatch", "effective_batch", "phases"}
    if rerun:
        fields.add("training_contract")
    if set(policy) != fields:
        raise ValueError("invalid policy fields")
    if policy["kind"] not in ("formal", "qualification") or (fixture and policy["kind"] == "formal"):
        raise ValueError("invalid fit purpose for this input")
    if policy["device"] not in ("cpu", "cuda") or policy["kernel"] not in KERNELS:
        raise ValueError("unsupported device or attention kernel")
    if policy["dtype"] not in ("float32", "bfloat16"):
        raise ValueError("unsupported precision")
    if policy["device"] == "cuda" and policy["kernel"] == "math":
        raise ValueError("CUDA runs must pin a fused attention kernel")
    if policy["device"] == "cpu" and (policy["kernel"], policy["dtype"]) != ("math", "float32"):
        raise ValueError("CPU qualification requires math attention and float32")
    if policy["kind"] == "formal" and policy["device"] != "cuda":
        raise ValueError("formal policy requires CUDA")
    microbatch, effective_batch = policy["microbatch"], policy["effective_batch"]
    if (type(microbatch) is not int or type(effective_batch) is not int
            or not 0 < microbatch <= effective_batch or effective_batch % microbatch):
        raise ValueError("effective batch must be a positive multiple of microbatch")
    if set(policy["phases"]) != set(PHASES):
        raise ValueError("policy requires Tahoe and tenv phases")
    if rerun:
        if effective_batch != 16 or microbatch not in (2, 4):
            raise ValueError("rerun requires effective batch 16 and microbatch 2 or 4")
        for source, phase in policy["phases"].items():
            fields = {"min_passes", "max_passes", "validation_interval", "learning_rate"}
            if policy["kind"] == "qualification" and "max_updates" in phase:
                fields.add("max_updates")
            if set(phase) != fields:
                raise ValueError("invalid rerun phase policy")
            if (phase["min_passes"] != 1 or phase["max_passes"] != (40 if source == "tahoe" else 5)
                    or phase["validation_interval"] != (50 if source == "tahoe" else 10)
                    or phase["learning_rate"] != 1e-4):
                raise ValueError("rerun phase differs from ADR 0060")
            if "max_updates" in phase and (type(phase["max_updates"]) is not int or phase["max_updates"] <= 0):
                raise ValueError("max_updates must be a positive integer")
        return policy
    for phase in policy["phases"].values():
        fields = {"updates", "validation_interval", "learning_rate"}
        if policy["kind"] == "qualification" and "validation_limit" in phase:
            fields.add("validation_limit")
        if set(phase) != fields:
            raise ValueError("invalid phase policy")
        for field in fields - {"learning_rate"}:
            if type(phase[field]) is not int or phase[field] <= 0:
                raise ValueError(f"{field} must be a positive integer")
        if not math.isfinite(phase["learning_rate"]) or phase["learning_rate"] <= 0:
            raise ValueError("learning rate must be finite and positive")
    return policy


@contextmanager
def precision(policy):
    with sdpa_kernel(KERNELS[policy["kernel"]]), torch.autocast(
        device_type=policy["device"], dtype=torch.bfloat16, enabled=policy["dtype"] == "bfloat16",
    ):
        yield


def batch_indices(size, batch_size, seed):
    """Stream shuffled rows without losing the tail at an epoch boundary."""
    rng = np.random.default_rng(seed)
    shuffled = rng.permutation(size)
    offset = 0
    while True:
        chunks, remaining = [], batch_size
        while remaining:
            n = min(remaining, size - offset)
            chunks.append(shuffled[offset:offset+n])
            remaining -= n
            offset += n
            if offset == size:
                shuffled = rng.permutation(size)
                offset = 0
        yield np.concatenate(chunks)


def batches(rows, microbatch):
    for start in range(0, len(rows), microbatch):
        yield rows[start:start+microbatch]


def assert_finite(terms, phase):
    if not torch.isfinite(terms).all().item():
        raise FloatingPointError(f"{phase}: non-finite objective")


def validate(model, cohort, policy, phase):
    model.eval()
    total = torch.zeros(len(LOSS_WEIGHTS), device=policy["device"])
    n = min(len(cohort), phase.get("validation_limit", len(cohort)))
    with torch.inference_mode(), precision(policy):
        for rows in batches(np.arange(n), policy["microbatch"]):
            inputs, target = cohort.batch(rows, policy["device"])
            total += loss_terms(model(*inputs), target, inputs[0]) * len(rows)
    assert_finite(total, "validation")
    return total / n


def train_phase(model, source, train, valid, output, run_key, policy, resume):
    device = policy["device"]
    phase = policy["phases"][source]
    checkpoint, last_path = output / f"{source}.pt", output / f"{source}_last.pt"

    def candidate(update):
        return output / f".{source}_best_{update}.pt"

    optimizer = torch.optim.AdamW(model.parameters(), lr=phase["learning_rate"], weight_decay=1e-5)
    batch_size, microbatch = policy["effective_batch"], policy["microbatch"]
    stream = batch_indices(len(train), batch_size, [run_key["seed"], PHASES.index(source)])
    weights = torch.tensor(LOSS_WEIGHTS, device=device)
    previous, elapsed, history = 0, 0.0, []
    best_update, best_loss = 0, math.inf
    if resume and last_path.exists():
        saved = torch.load(last_path, map_location="cpu", weights_only=True)
        if saved["run_key"] != run_key or saved["source"] != source:
            raise ValueError("resume state belongs to a different fit")
        previous, elapsed, history = saved["update"], saved["wall_seconds"], saved["history"]
        best_update, best_loss = saved["best_update"], saved["best_loss"]
        model.load_state_dict(saved["model"])
        optimizer.load_state_dict(saved["optimizer"])
        torch.set_rng_state(saved["cpu_rng"])
        if device == "cuda":
            torch.cuda.set_rng_state(saved["cuda_rng"])
        for _ in range(previous):
            next(stream)

    start = time.monotonic()
    window_start = start
    training_total = torch.zeros(len(LOSS_WEIGHTS), device=device)
    for update in range(previous + 1, phase["updates"] + 1):
        model.train()
        optimizer.zero_grad(set_to_none=True)
        for rows in batches(next(stream), microbatch):
            inputs, target = train.batch(rows, device)
            with precision(policy):
                outputs = model(*inputs)
            terms = loss_terms(outputs, target, inputs[0])
            weighted = (terms @ weights) * (len(rows) / batch_size)
            weighted.backward()
            training_total += terms.detach() * (len(rows) / batch_size)
        optimizer.step()

        if update % phase["validation_interval"] and update != phase["updates"]:
            continue
        assert_finite(training_total, source)
        validation_start = time.monotonic()
        validation = validate(model, valid, policy, phase)
        value = (validation @ weights).item()
        train_steps = update - (history[-1]["update"] if history else 0)
        row = {"update": update, "validation_loss": value, "validation_terms": validation.tolist(),
               "train_terms": (training_total / train_steps).tolist(),
               "train_seconds": validation_start - window_start,
               "validation_seconds": time.monotonic() - validation_start}
        history.append(row)
        training_total.zero_()
        previous_best = best_update
        if value < best_loss:
            best_update, best_loss = update, value
            save_checkpoint(candidate(update), {
                "schema": RUN_SCHEMA, "run_key": run_key, "source": source,
                "update": update, "validation_loss": value, "model": model.state_dict(),
            })
        save_checkpoint(last_path, {
            "run_key": run_key, "source": source, "update": update,
            "model": model.state_dict(), "optimizer": optimizer.state_dict(),
            "cpu_rng": torch.get_rng_state(),
            "cuda_rng": torch.cuda.get_rng_state() if device == "cuda" else None,
            "history": history, "best_update": best_update, "best_loss": best_loss,
            "wall_seconds": elapsed + time.monotonic() - start,
        })
        # The previous selection remains valid until the new last state commits.
        if previous_best and previous_best != best_update:
            candidate(previous_best).unlink(missing_ok=True)
        print(json.dumps({"source": source, **row}), flush=True)
        window_start = time.monotonic()

    selected = candidate(best_update)
    best = torch.load(selected, map_location="cpu", weights_only=True)
    if best["run_key"] != run_key or best["source"] != source or best["update"] != best_update:
        raise ValueError("selected checkpoint does not match the fit")
    model.load_state_dict(best["model"])
    staging = checkpoint.with_suffix(".link")
    staging.unlink(missing_ok=True)
    staging.hardlink_to(selected)
    staging.replace(checkpoint)
    return {"history": history, "best_update": best_update, "best_loss": best_loss,
            "updates": phase["updates"], "checkpoint_sha256": sha256(checkpoint),
            "wall_minutes": (elapsed + time.monotonic() - start) / 60}



def stopping_state(losses):
    """Replay the two strict phase-local stop counters from finite measurements."""
    increases = above_best = 0
    best = math.inf
    previous = None
    for value in losses:
        if not math.isfinite(value):
            raise FloatingPointError("non-finite validation objective")
        increases = increases + 1 if previous is not None and value > previous else 0
        above_best = above_best + 1 if value > best else 0
        best = min(best, value)
        previous = value
    return increases, above_best


def train_rerun_phase(model, source, train, valid, output, run_key, policy, resume=False):
    """ADR 0060 phase: a single atomic best slot and replayable exposure records."""
    phase, device = policy["phases"][source], policy["device"]
    batch_size = policy["effective_batch"]
    minimum = math.ceil(len(train) / batch_size)
    maximum = math.ceil(phase["max_passes"] * len(train) / batch_size)
    limit = min(maximum, phase.get("max_updates", maximum))
    checkpoint = output / "best.pt"
    initial_hash = sha256(checkpoint) if checkpoint.exists() else None
    optimizer = torch.optim.AdamW(model.parameters(), lr=phase["learning_rate"], weight_decay=1e-5)
    stream = batch_indices(len(train), batch_size, [run_key["seed"], PHASES.index(source)])
    weights = torch.tensor(LOSS_WEIGHTS, device=device)
    history, best_loss, best_update = [], math.inf, 0
    total = torch.zeros(len(LOSS_WEIGHTS), device=device)
    start = window_start = time.monotonic()
    stop_reason = None
    for update in range(1, limit + 1):
        model.train()
        optimizer.zero_grad(set_to_none=True)
        for rows in batches(next(stream), policy["microbatch"]):
            inputs, target = train.batch(rows, device)
            with precision(policy):
                outputs = model(*inputs)
            terms = loss_terms(outputs, target, inputs[0])
            assert_finite(terms, source)
            ((terms @ weights) * (len(rows) / batch_size)).backward()
            total += terms.detach() * (len(rows) / batch_size)
        optimizer.step()
        if update % phase["validation_interval"] and update != limit:
            continue
        validation_start = time.monotonic()
        validation = validate(model, valid, policy, phase)
        value = (validation @ weights).item()
        if not math.isfinite(value):
            raise FloatingPointError("non-finite validation objective")
        window_updates = update - (history[-1]["update"] if history else 0)
        row = {"update": update, "training_well_presentations": update * batch_size,
               "passes": update * batch_size / len(train), "validation_wells": len(valid),
               "validation_loss": value, "validation_terms": validation.tolist(),
               "train_terms": (total / window_updates).tolist(),
               "train_loss": ((total / window_updates) @ weights).item(),
               "train_seconds": validation_start - window_start,
               "validation_seconds": time.monotonic() - validation_start}
        history.append(row)
        total.zero_()
        increases, above_best = stopping_state([point["validation_loss"] for point in history])
        row.update(strict_increases=increases, strict_above_best=above_best,
                   stopping_eligible=update >= minimum)
        if value < best_loss:
            best_loss, best_update = value, update
            save_checkpoint(checkpoint, {
                "schema": RUN_SCHEMA, "training_contract": RERUN_CONTRACT,
                "run_key": run_key, "source": source, "update": update,
                "training_well_presentations": update * batch_size,
                "validation_loss": value, "model": model.state_dict(),
            })
        if update == limit:
            stop_reason = "qualification_limit" if limit < maximum else "maximum_passes"
        elif update >= minimum and increases >= 3:
            stop_reason = "three_strict_increases"
        elif update >= minimum and above_best >= 10:
            stop_reason = "ten_strict_above_best"
        write_json(output / f"{source}_history.json", {
            "source": source, "run_key": run_key, "history": history,
            "status": "phase_complete" if stop_reason else "running",
            "stop_reason": stop_reason, "best_update": best_update, "best_loss": best_loss,
        })
        print(json.dumps({"source": source, **row, "stop_reason": stop_reason}), flush=True)
        window_start = time.monotonic()
        if stop_reason:
            break
    best = torch.load(checkpoint, map_location="cpu", weights_only=True)
    if best["run_key"] != run_key or best["source"] != source or best["update"] != best_update:
        raise ValueError("selected checkpoint does not match the phase")
    model.load_state_dict(best["model"])
    return {"history": history, "best_update": best_update, "best_loss": best_loss,
            "updates": update, "minimum_updates": minimum, "maximum_updates": maximum,
            "training_wells": len(train), "validation_wells": len(valid),
            "training_well_presentations": update * batch_size,
            "passes": update * batch_size / len(train),
            "maximum_budget_overshoot_wells": maximum * batch_size - phase["max_passes"] * len(train),
            "best_training_well_presentations": best_update * batch_size,
            "best_passes": best_update * batch_size / len(train),
            "initial_checkpoint_sha256": initial_hash,
            "checkpoint_sha256": sha256(checkpoint), "stop_reason": stop_reason,
            "wall_minutes": (time.monotonic() - start) / 60}


def model_recipe(predictor, name):
    name = name or predictor.model
    if name in MLP_WIDTHS:
        if predictor.model != "XPert_G18164":
            raise ValueError(f"{name} requires the G18164 predictor")
    elif name != predictor.model:
        raise ValueError("selected model does not match its predictor gene axis")
    return name, MLP_WIDTHS.get(name)


def fit(predictor_path, output, policy_path, fold, seed, resume=False, model_name=None):
    predictor = Predictor(predictor_path, fold)
    name, latent_tokens = model_recipe(predictor, model_name)
    if seed not in predictor.manifest["training_seeds"]:
        raise ValueError("seed is outside the frozen run matrix")
    policy = load_policy(policy_path, predictor.manifest["fixture"])
    cohorts = {source: [predictor.dataset(source, role, targets=True)
                        for role in ("train", "validation")]
               for source in PHASES}
    run_key = {"predictor_sha256": predictor.digest, "fold": fold, "seed": seed,
               "policy": policy, "policy_sha256": sha256(policy_path), "model": name}
    rerun = policy.get("training_contract") == RERUN_CONTRACT
    if rerun and resume:
        raise ValueError("rerun interruption requires a fresh run directory")
    if resume:
        request = json.loads((output / "request.json").read_text())
        if request != run_key:
            raise ValueError("resume request differs from original fit")
        if (output / "run.json").exists():
            raise ValueError("fit already finished")
    else:
        output.mkdir(parents=True, exist_ok=False)
        write_json(output / "request.json", run_key)
    device = policy["device"]
    torch.manual_seed(seed)
    model = XPert(predictor.identity(), latent_tokens).to(device)
    if device == "cuda":
        torch.cuda.reset_peak_memory_stats()
    phases = {}
    for source in PHASES:
        train, valid = cohorts[source]
        trainer = train_rerun_phase if rerun else train_phase
        try:
            phases[source] = trainer(model, source, train, valid, output, run_key, policy, resume)
        except Exception as error:
            if rerun:
                write_json(output / "failure.json", {
                    "schema": RUN_SCHEMA, "status": "failed", "source": source,
                    "run_key": run_key, "error_type": type(error).__name__, "error": str(error),
                    "completed_phases": phases,
                })
            raise
    report = {
        "schema": RUN_SCHEMA, "status": "complete", **run_key,
        "phases": phases, "training_minutes": sum(phase["wall_minutes"] for phase in phases.values()),
        "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "peak_vram_allocated_bytes": torch.cuda.max_memory_allocated() if device == "cuda" else 0,
        "peak_vram_reserved_bytes": torch.cuda.max_memory_reserved() if device == "cuda" else 0,
        "torch_version": str(torch.__version__), "numpy_version": np.__version__,
        "architecture": model.architecture(),
        "architecture_sha256": sha256(Path(__file__).with_name("model.py")),
    }
    if rerun:
        report["training_contract"] = RERUN_CONTRACT
        report["purpose"] = "formal" if policy["kind"] == "formal" else "smoke"
        report["checkpoint"] = "best.pt"
        report["checkpoint_sha256"] = sha256(output / "best.pt")
    write_json(output / "run.json", report)
    return report



def require_completed_rerun(checkpoint, saved):
    if saved.get("training_contract") != RERUN_CONTRACT:
        return
    report_path = checkpoint.parent / "run.json"
    if not report_path.exists():
        raise ValueError("rerun checkpoint belongs to an incomplete run")
    report = json.loads(report_path.read_text())
    if (report.get("status") != "complete" or report.get("checkpoint_sha256") != sha256(checkpoint)
            or any(report.get(key) != saved["run_key"].get(key)
                   for key in saved["run_key"])):
        raise ValueError("rerun checkpoint has no matching completed run record")


def predict(predictor_path, checkpoint, role, output, limit=None):
    saved = torch.load(checkpoint, map_location="cpu", weights_only=True)
    require_completed_rerun(checkpoint, saved)
    if saved["schema"] != RUN_SCHEMA or saved["source"] != "tenv":
        raise ValueError("prediction requires a selected tenv checkpoint")
    run_key = saved["run_key"]
    predictor = Predictor(predictor_path, run_key["fold"])
    if predictor.digest != run_key["predictor_sha256"]:
        raise ValueError("checkpoint and predictor release do not match")
    _, latent_tokens = model_recipe(predictor, run_key["model"])
    policy = run_key["policy"]
    if limit is not None and (policy["kind"] != "qualification" or type(limit) is not int or limit <= 0):
        raise ValueError("a positive prediction limit is available only for qualification")
    cohort = predictor.dataset("tenv", role, targets=False)
    n_cases = len(cohort) if limit is None else min(len(cohort), limit)
    model = XPert(predictor.identity(), latent_tokens).to(policy["device"])
    model.load_state_dict(saved["model"])
    model.eval()
    if output.exists():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent, prefix=".predict-") as temporary:
        directory = Path(temporary) / "prediction"
        directory.mkdir()
        values = np.lib.format.open_memmap(directory / "delta.npy", mode="w+", dtype=np.float32,
                                           shape=(n_cases, len(predictor.genes)))
        with torch.inference_mode(), precision(policy):
            for rows in batches(np.arange(n_cases), policy["microbatch"]):
                inputs, _ = cohort.batch(rows, policy["device"])
                delta = model(*inputs)[2].float().cpu().numpy()
                if not np.isfinite(delta).all():
                    raise FloatingPointError("prediction contains non-finite delta")
                values[rows] = delta
        values.flush()
        del values
        cohort.cases.iloc[:n_cases][["case_id", "group_id", "drug_key", "dose_um", "control_index"]].to_csv(
            directory / "cases.csv", index=False
        )
        predictor.genes[["symbol", "hgnc_id"]].to_csv(directory / "genes.csv", index=False)
        manifest = {
            "schema": RUN_SCHEMA, "role": role, "run_key": run_key,
            "checkpoint_sha256": sha256(checkpoint),
            "files": {name: sha256(directory / name) for name in ("delta.npy", "genes.csv", "cases.csv")},
        }
        if limit is not None:
            manifest["qualification_limit"] = limit
        write_json(directory / "manifest.json", manifest)
        directory.rename(output)


def replay(predictor, checkpoint, prediction):
    manifest = json.loads((prediction / "manifest.json").read_text())
    with tempfile.TemporaryDirectory(dir=prediction.parent, prefix=".replay-") as temporary:
        regenerated = Path(temporary) / "prediction"
        predict(predictor, checkpoint, manifest["role"], regenerated, manifest.get("qualification_limit"))
        for name in ("manifest.json", "delta.npy", "cases.csv", "genes.csv"):
            if sha256(prediction / name) != sha256(regenerated / name):
                raise ValueError(f"checkpoint replay differs: {name}")
    print(json.dumps({"status": "verified", "prediction_manifest_sha256":
                      sha256(prediction / "manifest.json")}), flush=True)



def evaluate_loss(predictor_path, checkpoint, evaluator, output):
    """Evaluate selected weights on hash-verified native-axis held-out wells."""
    saved = torch.load(checkpoint, map_location="cpu", weights_only=True)
    require_completed_rerun(checkpoint, saved)
    if saved["schema"] != RUN_SCHEMA or saved["source"] != "tenv":
        raise ValueError("test objective requires selected tenv weights")
    key = saved["run_key"]
    predictor = Predictor(predictor_path, key["fold"])
    if predictor.digest != key["predictor_sha256"]:
        raise ValueError("checkpoint and predictor differ")
    evaluator = evaluator.absolute()
    if evaluator.resolve() != evaluator:
        raise ValueError("evaluator cannot be a symlink")
    manifest = json.loads((evaluator / "manifest.json").read_text())
    if (manifest["schema"] != "step09-evaluator-v1"
            or manifest["predictor_manifest_sha256"] != predictor.digest):
        raise ValueError("evaluator and predictor differ")
    relative = f"tenv/fold_{key['fold']}/test_expression.npy"
    target_path = evaluator / relative
    if target_path.resolve() != target_path or sha256(target_path) != manifest["files"][relative]:
        raise ValueError("test expression identity mismatch")
    treated = np.load(target_path, mmap_mode="r", allow_pickle=False)
    inputs = predictor.dataset("tenv", "test", targets=False)
    cohort = Profiles(inputs.cases, inputs.controls, inputs.morgan, treated, len(predictor.genes))
    _, latent = model_recipe(predictor, key["model"])
    policy = key["policy"]
    model = XPert(predictor.identity(), latent).to(policy["device"])
    model.load_state_dict(saved["model"])
    terms = validate(model, cohort, policy, {})
    payload = {"schema": "step09-native-test-loss-v1", "run_key": key,
               "checkpoint_sha256": sha256(checkpoint), "predictor_sha256": predictor.digest,
               "evaluator_sha256": sha256(evaluator / "manifest.json"),
               "test_expression_sha256": manifest["files"][relative],
               "test_wells": len(cohort), "gene_count": len(predictor.genes),
               "loss_terms": terms.tolist(), "loss_weights": list(LOSS_WEIGHTS),
               "test_loss": (terms @ torch.tensor(LOSS_WEIGHTS, device=policy["device"])).item()}
    output.parent.mkdir(parents=True, exist_ok=True)
    write_json(output, payload)
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    fit_parser = commands.add_parser("fit")
    fit_parser.add_argument("--policy", type=Path, required=True)
    fit_parser.add_argument("--fold", type=int, required=True)
    fit_parser.add_argument("--seed", type=int, required=True)
    fit_parser.add_argument("--output", type=Path, required=True)
    fit_parser.add_argument("--model", choices=(*AXIS_MODELS.values(), *MLP_WIDTHS),
                            help="select a model on its required predictor axis")
    fit_parser.add_argument("--resume", action="store_true")
    predict_parser = commands.add_parser("predict")
    predict_parser.add_argument("--role", choices=("validation", "test"), required=True)
    predict_parser.add_argument("--output", type=Path, required=True)
    predict_parser.add_argument("--limit", type=int, help="limit prediction rows in qualification runs")
    loss_parser = commands.add_parser("evaluate-loss")
    loss_parser.add_argument("--evaluator", type=Path, required=True)
    loss_parser.add_argument("--output", type=Path, required=True)
    replay_parser = commands.add_parser("replay")
    replay_parser.add_argument("--prediction", type=Path, required=True)
    for subparser in (fit_parser, predict_parser, replay_parser, loss_parser):
        subparser.add_argument("--predictor", type=Path, required=True)
    for subparser in (predict_parser, replay_parser, loss_parser):
        subparser.add_argument("--checkpoint", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "fit":
        fit(args.predictor, args.output, args.policy, args.fold, args.seed, args.resume, args.model)
    elif args.command == "predict":
        predict(args.predictor, args.checkpoint, args.role, args.output, args.limit)
    elif args.command == "evaluate-loss":
        evaluate_loss(args.predictor, args.checkpoint, args.evaluator, args.output)
    else:
        replay(args.predictor, args.checkpoint, args.prediction)


if __name__ == "__main__":
    main()
