"""Run the frozen Step 09 Round 1 model queue and verify its checkpoint predictions."""

import argparse
import fcntl
import json
from pathlib import Path
import subprocess
import sys

from round_contract import read_json, sha256, verify_committed, verify_committed_promotion


SCRIPT = Path(__file__).resolve().parent
STEP = SCRIPT.parents[1]
RECORD = STEP / "execution/round1_record"


def save_status(path, status):
    pending = path.with_suffix(".tmp")
    pending.write_text(json.dumps(status, indent=2) + "\n")
    pending.replace(path)


def check_digest(path, expected):
    if sha256(path) != expected:
        raise ValueError(f"frozen artifact has changed: {path}")


def job(output, model, stage, status, status_path, args):
    resources = output / "resources" / model
    attempt = status["attempts"].get(stage, 0) + 1
    resource_path = resources / f"{stage}-{attempt}"
    while resource_path.exists():
        attempt += 1
        resource_path = resources / f"{stage}-{attempt}"
    status["attempts"][stage] = attempt
    status.update(status="active", stage=stage, exit_code=None)
    save_status(status_path, status)
    return_code = subprocess.call([str(SCRIPT / "run_capped.sh"), str(resource_path),
                                   sys.executable, str(SCRIPT / "xpert.py"), *map(str, args)])
    if return_code:
        status.update(status="failed", stage=stage, exit_code=return_code)
        save_status(status_path, status)
        raise SystemExit(f"{model}: {stage} failed. See {resource_path}")
    return resource_path


def run_model(output, predictor, model, plan, policy, role):
    fit = output / "fits" / model
    checkpoint = fit / "tenv.pt"
    prediction = output / "predictions" / role / model
    status_path = output / "status" / f"{model}.json"
    status_path.parent.mkdir(parents=True, exist_ok=True)
    status = read_json(status_path) if status_path.exists() else {"model": model, "attempts": {}}
    if role == "test" and not status.get("validation_replay_sha256"):
        raise ValueError(f"{model}: test prediction requires a verified validation prediction")
    if role == "validation" and not (fit / "run.json").is_file():
        args = ["fit", "--model", model, "--predictor", predictor, "--policy", policy,
                "--fold", plan["fold"], "--seed", plan["seed"], "--output", fit]
        if (fit / "request.json").is_file():
            args += ["--resume"]
        job(output, model, "fit", status, status_path, args)
    run = read_json(fit / "run.json")
    if (run["status"], run["model"], run["policy_sha256"], run["predictor_sha256"], run["fold"], run["seed"]) != (
        "complete", model, plan["policy_sha256"], sha256(predictor / "manifest.json"), plan["fold"], plan["seed"]
    ):
        raise ValueError(f"{model}: fit does not match the frozen Round 1 plan")
    for source in ("tahoe", "tenv"):
        check_digest(fit / f"{source}.pt", run["phases"][source]["checkpoint_sha256"])
    if not prediction.exists():
        job(output, model, f"predict-{role}", status, status_path,
            ["predict", "--predictor", predictor, "--checkpoint", checkpoint,
             "--role", role, "--output", prediction])
    digest = sha256(prediction / "manifest.json")
    for name, expected in read_json(prediction / "manifest.json")["files"].items():
        check_digest(prediction / name, expected)
    replay_key = f"{role}_replay_sha256"
    if status.get(replay_key) != digest:
        replay = job(output, model, f"replay-{role}", status, status_path,
                     ["replay", "--predictor", predictor, "--checkpoint", checkpoint,
                      "--prediction", prediction])
        verification = read_json(replay / "stdout.log")
        if verification != {"status": "verified", "prediction_manifest_sha256": digest}:
            raise ValueError(f"{model}: invalid checkpoint replay proof")
        status[replay_key] = digest
    status.update(status="complete", stage=f"replay-{role}", exit_code=0)
    save_status(status_path, status)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--role", required=True, choices=("validation", "test"))
    parser.add_argument("--models", nargs="+")
    parser.add_argument("--output", type=Path, default=STEP / "output/round1")
    parser.add_argument("--corpus", type=Path, default=STEP / "output/corpus")
    parser.add_argument("--plan", type=Path, default=RECORD / "plan.json")
    parser.add_argument("--policy", type=Path, default=RECORD / "policy.json")
    parser.add_argument("--promotion", type=Path, default=RECORD / "promotion.json")
    args = parser.parse_args()

    verify_committed(args.plan)
    verify_committed(args.policy)
    plan = read_json(args.plan)
    if (plan["schema"], plan["fold"], plan["seed"]) != ("step09-round1-plan-v1", 0, 17):
        raise ValueError("expected the frozen Round 1 plan")
    check_digest(args.policy, plan["policy_sha256"])
    check_digest(args.corpus / "gene_axes/manifest.json", plan["gene_axes_manifest_sha256"])
    for name, digest in plan["fit_sources"].items():
        check_digest(SCRIPT / name, digest)
    if args.role == "test":
        promotion = verify_committed_promotion(args.promotion)
        if promotion["plan_sha256"] != sha256(args.plan):
            raise ValueError("committed promotion belongs to a different Round 1 plan")
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / ".queue.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        models = {entry["name"]: entry for entry in plan["models"]}
        for name in args.models or models:
            entry = models[name]
            predictor = args.corpus / entry["view"]
            check_digest(predictor / "manifest.json", entry["predictor_sha256"])
            run_model(args.output, predictor, name, plan, args.policy, args.role)


if __name__ == "__main__":
    main()
