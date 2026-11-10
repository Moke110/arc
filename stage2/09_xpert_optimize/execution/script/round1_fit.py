"""Run the frozen Step 09 Round 1 model queue and verify its checkpoint predictions."""

import argparse
import fcntl
import json
from pathlib import Path
import subprocess
import sys

from round_contract import DEFAULT_WORK, read_json, role_result, sha256, verify_committed, verify_committed_promotion


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


class JobFailed(RuntimeError):
    def __init__(self, return_code, resource_path):
        super().__init__(f"worker exited with {return_code}. See {resource_path}")
        self.return_code = return_code


def job(work, model, stage, status, status_path, args):
    resources = work / "resources" / model
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
        raise JobFailed(return_code, resource_path)
    return resource_path


def run_model(output, work, predictor, model, entry, plan, policy, role, status, status_path):
    fit = work / "fits" / model
    checkpoint = fit / "tenv.pt"
    prediction = output / "predictions" / role / model
    status.pop(f"{role}_failure", None)
    status.update(status="active", active_role=role, stage=f"verify-input-{role}", exit_code=None)
    save_status(status_path, status)
    check_digest(predictor / "manifest.json", entry["predictor_sha256"])
    if role == "validation" and not (fit / "run.json").is_file():
        args = ["fit", "--model", model, "--predictor", predictor, "--policy", policy,
                "--fold", plan["fold"], "--seed", plan["seed"], "--output", fit]
        if (fit / "request.json").is_file():
            args += ["--resume"]
        job(work, model, "fit", status, status_path, args)
    status["stage"] = f"verify-fit-{role}"
    run = read_json(fit / "run.json")
    if (run["status"], run["model"], run["policy_sha256"], run["predictor_sha256"], run["fold"], run["seed"]) != (
        "complete", model, plan["policy_sha256"], sha256(predictor / "manifest.json"), plan["fold"], plan["seed"]
    ):
        raise ValueError(f"{model}: fit does not match the frozen Round 1 plan")
    for source in ("tahoe", "tenv"):
        check_digest(fit / f"{source}.pt", run["phases"][source]["checkpoint_sha256"])
    if not prediction.exists():
        job(work, model, f"predict-{role}", status, status_path,
            ["predict", "--predictor", predictor, "--checkpoint", checkpoint,
             "--role", role, "--output", prediction])
    status["stage"] = f"verify-prediction-{role}"
    digest = sha256(prediction / "manifest.json")
    for name, expected in read_json(prediction / "manifest.json")["files"].items():
        check_digest(prediction / name, expected)
    replay_key = f"{role}_replay_sha256"
    if status.get(replay_key) != digest:
        replay = job(work, model, f"replay-{role}", status, status_path,
                     ["replay", "--predictor", predictor, "--checkpoint", checkpoint,
                      "--prediction", prediction])
        status["stage"] = f"verify-replay-{role}"
        verification = read_json(replay / "stdout.log")
        if verification != {"status": "verified", "prediction_manifest_sha256": digest}:
            raise ValueError(f"{model}: invalid checkpoint replay proof")
        status[replay_key] = digest
    status.pop("active_role", None)
    status.update(status="complete", stage=f"replay-{role}", exit_code=0)
    save_status(status_path, status)


def run_queue(output, work, corpus, plan, policy, role, names=None):
    models = {entry["name"]: entry for entry in plan["models"]}
    status_dir = work / "status"
    status_dir.mkdir(parents=True, exist_ok=True)
    failures = []
    for name in names or models:
        entry = models[name]
        status_path = status_dir / f"{name}.json"
        status = read_json(status_path) if status_path.exists() else {"model": name, "attempts": {}}
        if role == "test":
            validation, _ = role_result(status, "validation")
            if validation == "failed" and names is None:
                print(f"{name}: skipped test after terminal validation failure", flush=True)
                continue
            if validation != "complete":
                failures.append(name)
                print(f"{name}: test requires a verified validation replay", file=sys.stderr)
                continue
        try:
            run_model(output, work, corpus / entry["view"], name, entry, plan, policy, role,
                      status, status_path)
        except (JobFailed, OSError, ValueError, TypeError, KeyError, KeyboardInterrupt) as error:
            interrupted = isinstance(error, KeyboardInterrupt)
            exit_code = error.return_code if isinstance(error, JobFailed) else 1
            if interrupted:
                exit_code = 130
            failure = {"stage": status["stage"], "exit_code": exit_code, "message": str(error)}
            status.pop(f"{role}_replay_sha256", None)
            status.pop("active_role", None)
            status[f"{role}_failure"] = failure
            status.update(status="failed", exit_code=exit_code)
            save_status(status_path, status)
            if interrupted:
                raise
            failures.append(name)
            print(f"{name}: {failure['stage']} failed. {error}", file=sys.stderr)
    return failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--role", required=True, choices=("validation", "test"))
    parser.add_argument("--models", nargs="+")
    parser.add_argument("--output", type=Path, default=STEP / "output/round1")
    parser.add_argument("--work-root", type=Path, default=DEFAULT_WORK)
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
    args.work_root.mkdir(parents=True, exist_ok=True)
    with (args.work_root / ".queue.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        failures = run_queue(args.output, args.work_root, args.corpus, plan, args.policy, args.role, args.models)
        if failures:
            raise SystemExit(f"Round 1 {args.role} failures: {', '.join(failures)}")


if __name__ == "__main__":
    main()
