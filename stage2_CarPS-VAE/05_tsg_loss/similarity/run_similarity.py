"""Register, fit and exactly replay one bounded TranSiGen Delta similarity coordinate."""

import argparse
from fractions import Fraction
from functools import partial
import json
import os
from pathlib import Path
import resource
from tempfile import TemporaryDirectory
import time

import torch

from similarity_objective import ALPHAS, METRICS, SimilarityObjective
from run_reconstruction import bind as bind_reconstruction

from run_sweep import fit
from run_loss import (AXIS, DEFAULT_REPRESENTATIONS, DEFAULT_SOURCE, DEFAULT_TRAINING,
                      INPUT_PATHS, PHASES, SEALED_FILES, axis_run)

EXECUTION = Path(__file__).resolve().parent
SCHEMA = "stage2-tenv-transigen-bounded-delta-similarity-run-v1"
PROTOCOL = "matched-plate-control-bounded-delta-similarity-v1"


def bind(args, data, vectors):
    objective = SimilarityObjective(args["metric"], Fraction(args["alpha"]))
    config = bind_reconstruction({**args, "alpha": "1"}, data, vectors)
    config.update(schema=SCHEMA, protocol=PROTOCOL, objective=objective.identity())
    config["model"]["validation_objective"] = "registered bounded Delta similarity objective per condition"
    for name in ("run_similarity.py", "similarity_objective.py"):
        path = EXECUTION / name
        config["code_sha256"][path.relative_to(axis_run.REPO).as_posix()] = axis_run.digest(path)
    return config


def register(output, args):
    with axis_run.publication_lock(output) as root:
        if root.exists():
            raise ValueError("run registration requires an unused output directory")
        with axis_run.open_inputs(args["training"], args["representations"]) as (data, vectors):
            config = bind(args, data, vectors)
        root.mkdir()
        record = root / "registration.json"
        axis_run.write_json(record, {**config, "arguments": args, "nonce": axis_run.PROCESS})
        return axis_run.digest(record)


def execute(output, phase, paths, metric):
    start = time.monotonic()
    with axis_run.publication_lock(output) as root:
        registration = root / "registration.json"
        config = axis_run.read_json(registration)
        registration_sha = axis_run.digest(registration)
        if set(paths) - set(INPUT_PATHS):
            raise ValueError("only verified input locations can be overridden")
        args = {**config["arguments"], **paths}
        if args["metric"] != metric:
            raise ValueError("run registration belongs to a different similarity study")
        receipt = PHASES[phase]
        sealed = None
        if (root / receipt).exists():
            if phase != "replay":
                raise ValueError(f"run phase already recorded: {phase}")
            sealed = axis_run.verified_receipt(root, receipt, registration_sha, SEALED_FILES)
        if phase != "train":
            training = axis_run.verified_receipt(root, PHASES["train"], registration_sha, axis_run.TRAINING_FILES)
        if phase == "replay":
            axis_run.verified_receipt(root, PHASES["predict"], registration_sha, axis_run.PREDICTION_FILES)
        os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
        with axis_run.open_inputs(args["training"], args["representations"]) as (data, vectors):
            identity = {key: value for key, value in config.items() if key not in {"arguments", "nonce"}}
            if bind(args, data, vectors) != identity:
                raise ValueError("registered experiment, inputs, code or runtime changed")
            torch.use_deterministic_algorithms(True)
            axis_run.kpgt.configure_backend(args["device"])
            native = axis_run.kpgt.adapter(args["source"], config["model"]["transform"])
            model = native.build(config["seed"], args["device"], n_genes=len(data.genes))
            objective = SimilarityObjective(metric, Fraction(config["objective"]["alpha"]))
            model.loss = partial(objective, model.loss)
            with TemporaryDirectory(prefix=".phase-", dir=root) as scratch_path:
                scratch = Path(scratch_path)
                if phase == "train":
                    details = fit(scratch, native, model, data, config, registration_sha)
                    files = axis_run.TRAINING_FILES
                else:
                    values, test = axis_run.predict_test(root, native, model, data, config, registration_sha, training)
                    details = axis_run.record_predictions(root, scratch, values, test, data, config, replay=phase == "replay")
                    files = axis_run.PREDICTION_FILES
                if sealed is not None:
                    return sealed
                if phase == "replay":
                    details.update(test_replay="exact", axis=AXIS, fold=args["fold"],
                                   smoke=args["smoke"], seed=config["seed"], schema=SCHEMA,
                                   protocol=PROTOCOL, objective=config["objective"])
                    hashes = axis_run.file_hashes(root, SEALED_FILES)
                else:
                    hashes = axis_run.file_hashes(scratch, files)
                    for name in files:
                        (scratch / name).replace(root / name)
                record = {
                    **details, "outputs": hashes, "registration_sha256": registration_sha,
                    "process_token": axis_run.PROCESS, "pid": os.getpid(),
                    "wall_seconds": round(time.monotonic() - start, 3),
                    "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                }
                axis_run.write_json(scratch / receipt, record)
                (scratch / receipt).replace(root / receipt)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="phase", required=True)
    registration = commands.add_parser("register")
    registration.add_argument("--alpha", choices=ALPHAS, required=True)
    registration.add_argument("--fold", choices=range(3), type=int, required=True)
    registration.add_argument("--training", type=Path, default=DEFAULT_TRAINING)
    registration.add_argument("--representations", type=Path, default=DEFAULT_REPRESENTATIONS)
    registration.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    registration.add_argument("--device", choices=("cpu", "cuda"), default="cuda")
    registration.add_argument("--smoke", action="store_true")
    for phase in PHASES:
        command = commands.add_parser(phase)
        for field in INPUT_PATHS:
            command.add_argument(f"--{field}", type=Path, help="same pinned input at a new location")
    for command in commands.choices.values():
        command.add_argument("--output", type=Path, required=True)
        command.add_argument("--metric", choices=METRICS, required=True)
    args = vars(parser.parse_args())
    phase, output = args.pop("phase"), args.pop("output").resolve()
    options = {key: str(value.resolve()) if isinstance(value, Path) else value
               for key, value in args.items() if value is not None}
    if phase == "register":
        result = {"registration_sha256": register(output, {"axis": AXIS, **options})}
    else:
        metric = options.pop("metric")
        result = execute(output, phase, options, metric)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
