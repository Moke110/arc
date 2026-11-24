"""Preregister, train and exactly replay all 18 formal TranSiGen axis fits."""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

STEP = Path(__file__).resolve().parents[1]
STAGE = STEP.parent
REPO = STAGE.parent
sys.path.insert(0, str(STEP / "output/api"))

from tsg_axis.inputs import AxisInputs
from training_data import TENV_AXIS_STEPS, open_tenv_axis_training, tenv_axis_output, sha256 as file_sha256
from stage2_util.derivation.release import write_json
from stage2_util.publication import publication_lock

SCHEMA = "stage2-tenv-transigen-axis-matrix-v1"
PARTITIONS = ("train", "valid", "test")
PHASES = {"train": "training.json", "predict": "predictions.json", "replay": "manifest.json"}


def match_rosters(axes):
    """Require the same condition order in every partition of every fold."""
    first, *rest = axes.items()
    reference = first[1]
    for axis, folds in rest:
        for fold, (current, expected) in enumerate(zip(folds, reference, strict=True)):
            for partition in PARTITIONS:
                if current[partition] != expected[partition]:
                    raise ValueError(f"{axis} fold {fold} {partition} condition order differs")
    return reference


def read(path):
    return json.loads(path.read_bytes())


def roster_digest(values):
    return sha256(json.dumps(values, separators=(",", ":")).encode()).hexdigest()


def preflight(step03):
    rosters, releases = {}, {}
    common = None
    for axis in TENV_AXIS_STEPS:
        with open_tenv_axis_training(axis, step03) as release:
            data = AxisInputs(release)
            rosters[axis] = [{part: list(data.condition_ids(data.rows(fold, part))) for part in PARTITIONS}
                             for fold in range(3)]
            inputs = release.manifest["inputs"]
            identity = {"raw_sha256": inputs["raw_sha256"], "folds": inputs["folds"],
                        "members_sha256": release.manifest["outputs"]["members.csv"],
                        "exclusions_sha256": release.manifest["outputs"]["exclusions.csv"],
                        "conditions_sha256": roster_digest(data.condition_ids(np.arange(len(data.treatments))))}
            if common is None:
                common = identity
            if common != identity:
                raise ValueError(f"{axis} differs from common raw provenance or complete membership")
            releases[axis] = release.manifest_sha256
    matched = match_rosters(rosters)
    counts = [{part: len(values) for part, values in fold.items()} for fold in matched]
    if any(not count for fold in counts for count in fold.values()):
        raise ValueError("formal folds must contain all three partitions")
    evidence = {"training_manifest_sha256": releases, **common, "fold_counts": counts,
                "ordered_fold_condition_sha256": [{part: roster_digest(values)
                                                   for part, values in fold.items()} for fold in matched]}
    return matched, evidence


def run_paths(root):
    return {(axis, fold): root / "runs" / axis / f"fold-{fold}"
            for axis, fold in product(TENV_AXIS_STEPS, range(3))}


def invoke(log, phase, args):
    log.mkdir(parents=True, exist_ok=True)
    attempt = log / f"{phase}-{time.time_ns()}"
    with attempt.with_suffix(".log").open("xb") as output:
        subprocess.run(["/usr/bin/time", "-v", "-o", str(attempt.with_suffix(".resources")),
                        sys.executable, str(STEP / "execution/run.py"), phase, *args],
                       stdout=output, stderr=subprocess.STDOUT, check=True)


def check_registrations(paths, rosters, preflight_record):
    registrations = {}
    nonces = set()
    code = None
    for (axis, fold), path in paths.items():
        config = read(path / "registration.json")
        args = config["arguments"]
        if ((args["axis"], args["fold"], args["smoke"], config["scenario"], config["seed"])
                != (axis, fold, False, "cold-drug", 42)
                or config["conditions"] != rosters[fold] or config["max_epochs"] != 1000
                or config["patience"] != 50 or config["training_sha256"] !=
                preflight_record["training_manifest_sha256"][axis]
                or config["folds"] != preflight_record["folds"] or config["nonce"] in nonces):
            raise ValueError(f"{axis} fold {fold} differs from the fresh formal registration plan")
        if code is None:
            code = config["code_sha256"]
        if code != config["code_sha256"]:
            raise ValueError("formal registrations use different runtime source")
        nonces.add(config["nonce"])
        registrations[axis, fold] = file_sha256(path / "registration.json")
    return registrations


def verify_run(path, axis, fold, registration):
    manifest = read(path / "manifest.json")
    training, prediction = (read(path / PHASES[phase]) for phase in ("train", "predict"))
    if ((manifest["axis"], manifest["fold"], manifest["smoke"], manifest["seed"],
         manifest["validation_replay"], manifest["test_replay"]) != (axis, fold, False, 42, "exact", "exact")
        or any(record["registration_sha256"] != registration for record in (training, prediction, manifest))
        or len({record["process_token"] for record in (training, prediction, manifest)}) != 3):
        raise ValueError(f"{axis} fold {fold} lacks distinct, exact fresh-process evidence")
    if any(file_sha256(path / name) != digest for name, digest in manifest["outputs"].items()):
        raise ValueError(f"{axis} fold {fold} run payload changed")
    return {"axis": axis, "fold": fold, "registration_sha256": registration,
            "run_manifest_sha256": file_sha256(path / "manifest.json"),
            "epochs": training["epochs"], "best_epoch": training["best_epoch"],
            "validation_mse": training["validation_mse"], "test_conditions": manifest["test_conditions"],
            "train_wall_seconds": training["wall_seconds"], "train_peak_rss_kib": training["peak_rss_kib"]}


def execute(root, step03, representations, source, device):
    root.mkdir(parents=True, exist_ok=True)
    with publication_lock(root) as output:
        paths = run_paths(output)
        rosters, evidence = preflight(step03)
        preflight_path = output / "preflight.json"
        plan = {"schema": SCHEMA, "device": device, "seed": 42, "scenario": "cold-drug",
                "driver_sha256": file_sha256(Path(__file__)), **evidence}
        if preflight_path.exists() and read(preflight_path) != plan:
            raise ValueError("formal matrix preflight changed")
        if not preflight_path.exists():
            write_json(preflight_path, plan)
        print(json.dumps({"phase": "preflight", "matched_axes": len(TENV_AXIS_STEPS),
                          "fold_counts": evidence["fold_counts"]}), flush=True)
        for (axis, fold), path in paths.items():
            if (path / "registration.json").is_file():
                continue
            args = ["--axis", axis, "--fold", str(fold), "--device", device,
                    "--training", str(tenv_axis_output(axis, step03)), "--representations", str(representations),
                    "--source", str(source), "--output", str(path)]
            invoke(output / "logs" / axis / f"fold-{fold}", "register", args)
        registrations = check_registrations(paths, rosters, evidence)
        print(json.dumps({"phase": "registered", "coordinates": len(registrations)}), flush=True)
        completed = []
        for (axis, fold), path in paths.items():
            print(json.dumps({"phase": "fit", "axis": axis, "fold": fold,
                              "completed": len(completed)}), flush=True)
            for phase, receipt in PHASES.items():
                if not (path / receipt).is_file():
                    invoke(output / "logs" / axis / f"fold-{fold}", phase, ["--output", str(path)])
            result = verify_run(path, axis, fold, registrations[axis, fold])
            completed.append(result)
            print(json.dumps({"phase": "complete", "completed": len(completed), **result}), flush=True)
        receipt = {"schema": SCHEMA, "preflight_sha256": file_sha256(preflight_path), "fits": completed,
                   "resource_logs": {path.relative_to(output).as_posix(): file_sha256(path)
                                     for path in sorted((output / "logs").rglob("*")) if path.is_file()}}
        write_json(output / "manifest.json", receipt)
        return {"result": "PASS", "fits": len(completed), "manifest_sha256": file_sha256(output / "manifest.json")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step03", type=Path, default=STAGE / "03_data_train")
    parser.add_argument("--representations", type=Path, default=STAGE / "util/output/representations")
    parser.add_argument("--source", type=Path, default=REPO / "pert_model/TranSiGen/source")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cuda")
    args = parser.parse_args()
    print(json.dumps(execute(args.output.resolve(), args.step03.resolve(), args.representations.resolve(),
                             args.source.resolve(), args.device), sort_keys=True))


if __name__ == "__main__":
    main()
