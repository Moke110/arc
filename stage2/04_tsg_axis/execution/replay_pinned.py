"""Replay sealed axis or downstream study checkpoints with their original Git runtime."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

STEP = Path(__file__).resolve().parents[1]
REPO = next(parent for parent in STEP.parents if (parent / ".git").exists())
sys.path.insert(0, str(REPO / "stage2/util/output/api"))

from stage2_util.derivation.release import check_outputs, read_json, sha256, write_json

AXIS_COMMIT = "986f8eacb2cda5b9e0793c03e12a27761aaa1c36"
LEGACY_COMMIT = "e4587a87d246f315d13389d6218d26320a5e9aba"


def require_runtime(runtime, commit):
    head = subprocess.check_output(["git", "-C", str(runtime), "rev-parse", "HEAD"], text=True).strip()
    status = subprocess.check_output(["git", "-C", str(runtime), "status", "--porcelain",
                                      "--untracked-files=all", "--", "stage2"], text=True)
    if head != commit or status:
        raise ValueError(f"replay requires a clean runtime at {commit}")


def replay_axes(args):
    require_runtime(args.runtime, AXIS_COMMIT)
    script = args.runtime / "stage2/04_tsg_axis/execution/run.py"
    sys.path.insert(0, str(REPO / "stage2/03_data_train/output/api"))
    from training_data import tenv_axis_output

    path = args.matrix / "manifest.json"
    matrix_sha = sha256(path)
    if matrix_sha != read_json(STEP / "log/316-formal.json")["matrix"]["manifest_sha256"]:
        raise ValueError("matrix differs from the accepted formal replay receipt")
    matrix = read_json(path)
    if len(matrix["fits"]) != 18 or len({(row["axis"], row["fold"]) for row in matrix["fits"]}) != 18:
        raise ValueError("sealed comparison requires 18 distinct coordinates")
    if args.scratch.resolve().is_relative_to(args.matrix.resolve()):
        raise ValueError("replay scratch must be outside the sealed matrix")
    args.scratch.mkdir(parents=True, exist_ok=True)
    for record in matrix["fits"]:
        axis, fold = record["axis"], record["fold"]
        path = args.matrix / "runs" / axis / f"fold-{fold}"
        if sha256(path / "manifest.json") != record["run_manifest_sha256"]:
            raise ValueError(f"sealed prediction identity differs: {axis}, fold {fold}")
        with TemporaryDirectory(prefix="axis-replay-", dir=args.scratch) as staging:
            copy = Path(staging) / "run"
            shutil.copytree(path, copy, copy_function=os.link)
            command = [sys.executable, str(script), "replay", "--output", str(copy),
                       "--training", str(tenv_axis_output(axis, args.step03)),
                       "--representations", str(args.representations), "--source", str(args.source)]
            subprocess.run(command, check=True, stdout=subprocess.DEVNULL)
        print(json.dumps({"axis": axis, "fold": fold, "replay": "exact"}), flush=True)
    return {"runtime": AXIS_COMMIT, "fits": 18, "matrix_sha256": matrix_sha}


def replay_dose(study, args):
    import numpy as np
    from tsg_kpgt.learner import load, predict, verify_candidate

    verify_candidate(args.run)
    data, native, model, config = load(args.inputs, args.run, args.source)
    rows = {"validation": np.asarray(config["validation_order"]),
            "test": data.rows(config["scenario"], config["fold"], "test")}
    seeds = {"validation": config["validation_seed"], "test": config["inference_seed"]}
    for phase, positions in rows.items():
        original = np.load(args.run / f"{phase}_prediction.npy", mmap_mode="r", allow_pickle=False)
        if not np.array_equal(original, predict(native, model, data, positions, seeds[phase])):
            raise ValueError(f"sealed dose {phase} predictions differ")
    return {"study": study, "validation": "exact", "test_prediction": "exact"}


def replay_width(study, args):
    runner = args.runtime / "stage2" / "05_tsg_width_dose" / "execution" / "learner.py"
    config = read_json(args.run / "config.json")
    manifest = read_json(args.run / "manifest.json")
    if manifest["identity"] != config["identity"] or manifest["origin"] != config["origin"]:
        raise ValueError("sealed width configuration differs")
    check_outputs(args.run, {"outputs": manifest["files"]})
    for name, digest in config["identity"]["code"].items():
        if Path(name).name != name or sha256(runner.parent / name) != digest:
            raise ValueError(f"sealed width code differs: {name}")
    for name, digest in config["identity"]["source"].items():
        if Path(name).name != name or sha256(args.source / "src" / name) != digest:
            raise ValueError(f"sealed width native code differs: {name}")
    if args.scratch.resolve().is_relative_to(args.run.resolve()):
        raise ValueError("replay scratch must be outside the sealed width run")
    args.scratch.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix="width-replay-", dir=args.scratch) as staging:
        copy = Path(staging) / "run"
        shutil.copytree(args.run, copy)
        terminal = read_json(copy / "terminal.json")
        if terminal["status"] != "complete":
            raise ValueError("width checkpoint must come from a complete fit")
        write_json(copy / "terminal.json", {**terminal, "status": "trained"})
        subprocess.run([sys.executable, str(runner), str(args.inputs), str(args.source), str(copy), "test"],
                       check=True, stdout=subprocess.DEVNULL)
        for name in ("test_metrics.json", "test_per_sample.csv"):
            if sha256(copy / name) != manifest["files"][name]:
                raise ValueError(f"sealed width test rows differ: {name}")
    return {"study": study, "validation": "exact", "test_scores": "exact"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    axes = commands.add_parser("axes")
    axes.add_argument("--matrix", type=Path, required=True)
    axes.add_argument("--step03", type=Path, required=True)
    axes.add_argument("--representations", type=Path, required=True)
    axes.add_argument("--scratch", type=Path, default=REPO / "tmp/axis-replay")
    for study in ("dose", "width"):
        command = commands.add_parser(study)
        command.add_argument("--inputs", type=Path, required=True)
        command.add_argument("--run", type=Path, required=True)
    commands.choices["width"].add_argument("--scratch", type=Path, required=True)
    for command in commands.choices.values():
        command.add_argument("--runtime", type=Path, required=True)
        command.add_argument("--source", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "axes":
        result = replay_axes(args)
    else:
        require_runtime(args.runtime, LEGACY_COMMIT)
        api = args.runtime / "stage2" / "04_tsg_dose" / "output/api"
        sys.path.insert(0, str(api))
        replay = {"dose": replay_dose, "width": replay_width}[args.command]
        result = replay(args.command, args) | {"runtime": LEGACY_COMMIT, "run_sha256": sha256(args.run / "manifest.json")}
    print(json.dumps({"result": "PASS", **result}, sort_keys=True))


if __name__ == "__main__":
    main()
