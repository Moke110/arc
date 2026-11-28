"""Register and execute one independent CarPS-Transcriptome axis fit and replay."""

import argparse
from copy import deepcopy
from pathlib import Path
import json
import os

import torch

from axis_protocol import ARM, AXES, REFERENCE_AXIS, SEEDS, coordinates
from axis_model import AxisAdapter, paired_initialization
import chem_runs
from chem_runs import ChemicalInputs, PHASES, PHASE_FILES, SEALED_FILES, axis_run
from evaluate_loss import axis_evaluate, check_trace
from evaluate_sweep import TRACE_COLUMNS
from run_measured import DEFAULT_REPRESENTATIONS, DEFAULT_SOURCE
from util.training.transigen import fit_with_diagnostics, kl_diagnostics
import tenv_axes
from util.training import registered
from util.training.provenance import verify_replay_configuration, verify_sealed_code

SCHEMA = "cpst-unimol-axis-run-v1"
PROTOCOL = "tenv-matched-control-cpst-unimol-axis-v1"
ORIGINAL_CODE = "d3395ca0e0a055dc28e8e431740f5ec640ee4247"
INHERITED_IDENTITIES = ("raw_sha256", "folds", "native_sha256", "partitions", "conditions", "environment",
                        "scenario", "validation_seed", "test_seed", "max_epochs", "patience")


def axis_code_paths():
    return (Path(__file__), Path(__file__).with_name("axis_model.py"),
            Path(__file__).with_name("axis_protocol.py"), Path(tenv_axes.__file__))


def historical_code_paths():
    """Use the accepted alpha ancestry and Step 09's four axis dependencies."""
    axis_paths = tuple(path.relative_to(axis_run.REPO).as_posix() for path in axis_code_paths())
    return (*chem_runs.historical_code_paths(), *axis_paths)


def arguments(key, bounded):
    axis, fold, seed = key
    if type(fold) is not int or type(seed) is not int or key not in coordinates(bounded):
        raise ValueError(f"unregistered Step 09 coordinate: {key}")
    return {"axis": axis, "fold": fold, "seed": seed, "smoke": bounded, "device": "cuda"}


def step08_config(fold, seed, step03, representations, source, bounded):
    """Recreate the pinned native reference configuration using the verified releases."""
    training = tenv_axes.tenv_axis_output(REFERENCE_AXIS, step03)
    with chem_runs.open_inputs(training, representations) as (data, vectors):
        return chem_runs.bind((ARM, fold, seed), data, vectors, source, bounded)


def bind(key, data, vectors, source, bounded, reference):
    """Bind just the axis and expression widths to the inherited UniMol run."""
    args = arguments(key, bounded)
    current = axis_run.bind({**args, "source": str(source)}, data, vectors)
    if (any(current[field] != reference[field] for field in INHERITED_IDENTITIES)
            or current["model"]["transform"] != reference["model"]["transform"]):
        raise ValueError(f"axis {key[0]} differs from the inherited reference population or fit protocol")
    for field in ("kpgt_manifest_sha256", "kpgt_vectors_sha256"):
        if current[field] != reference[field]:
            raise ValueError("frozen chemical representation release changed")
    config = deepcopy(reference)
    config.update(schema=SCHEMA, protocol=PROTOCOL)
    for field in ("training_sha256", "axis_sha256", "axis_family_manifest_sha256"):
        config[field] = current[field]
    width = len(data.genes)
    model = config["model"]
    model["n_genes"] = width
    model["adaptations"] = [f"{width} gene interface" if item == "5000 gene interface" else item
                            for item in model["adaptations"]]
    for path in axis_code_paths():
        config["code_sha256"][path.relative_to(axis_run.REPO).as_posix()] = axis_run.digest(path)
    return axis_run.transigen.bind_experiment(args, config)


def configure_cuda():
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    torch.use_deterministic_algorithms(True)
    axis_run.kpgt.configure_backend("cuda")


def register(output, key, bounded, step03, representations, source):
    def definition():
        args = arguments(key, bounded)
        reference = step08_config(key[1], key[2], step03, representations, source, bounded)
        training = tenv_axes.tenv_axis_output(key[0], step03)
        with chem_runs.open_inputs(training, representations) as (data, vectors):
            config = bind(key, data, vectors, source, bounded, reference)
        return {**config, "arguments": args}

    return registered.register(output, definition)


def execute(output, phase, step03, representations, source):
    configure_cuda()
    with registered.phase(output, phase, PHASE_FILES) as run:
        config = run.config
        args = config["arguments"]
        key = args["axis"], args["fold"], args["seed"]
        if args != arguments(key, args["smoke"]):
            raise ValueError("run arguments differ from its coordinate")
        reference = step08_config(key[1], key[2], step03, representations, source, args["smoke"])
        training = tenv_axes.tenv_axis_output(key[0], step03)
        with chem_runs.open_inputs(training, representations) as (data, vectors):
            identity = {name: value for name, value in config.items() if name not in ("arguments", "nonce")}
            current = bind(key, data, vectors, source, args["smoke"], reference)
            verify_replay_configuration(identity, current, replay=run.sealed is not None,
                                        repo=axis_run.REPO, revision=ORIGINAL_CODE,
                                        historical_paths=historical_code_paths())
            adapter = AxisAdapter(source, config["model"]["transform"])
            model = adapter.build(key[2], "cuda", n_genes=len(data.genes))
            model_data = ChemicalInputs(data, vectors, ARM)
            if phase == "train":
                initial = paired_initialization(model)
                details = fit_with_diagnostics(run.scratch, adapter, model, model_data, config, run.registration_sha)
                details.update(initialization=initial, epoch_streams=adapter.epoch_streams)
            else:
                values, rows = axis_run.predict_test(run.root, adapter, model, model_data, config,
                                                     run.registration_sha, run.training)
                details = axis_run.record_predictions(run.root, run.scratch, values, rows, model_data, config,
                                                      replay=phase == "replay")
        return run.publish(details, {**args, "schema": SCHEMA, "protocol": PROTOCOL, "objective": config["objective"]})


def read_run(path, bounded):
    """Read and validate all independently generated run evidence."""
    registration = path / "registration.json"
    config = axis_run.read_json(registration)
    verify_sealed_code(config, axis_run.REPO, ORIGINAL_CODE, historical_paths=historical_code_paths())
    args = config["arguments"]
    key = args["axis"], args["fold"], args["seed"]
    registration_sha = axis_run.digest(registration)
    if (args != arguments(key, bounded) or config["schema"] != SCHEMA or config["protocol"] != PROTOCOL
            or config["seed"] != key[2] or config["alpha_release_sha256"] != chem_runs.alpha_release_digest()
            or config["model"]["n_genes"] != int(key[0].split("-")[1][0]) * 1000):
        raise ValueError(f"invalid sealed Step 09 registration: {path}")
    training, prediction, replay = (
        axis_run.verified_receipt(path, PHASES[phase], registration_sha, PHASE_FILES[phase]) for phase in PHASES
    )
    if (len({config["nonce"], training["process_token"], prediction["process_token"], replay["process_token"]}) != 4
            or any(record.get(field) != expected for record in (training, prediction, replay)
                   for field, expected in {**args, "schema": SCHEMA, "protocol": PROTOCOL,
                                           "objective": config["objective"]}.items())
            or any(record["cuda_peak_allocated_bytes"] <= 0 for record in (training, prediction, replay))
            or any(replay.get(field) != "exact" for field in
                   ("validation_replay", "validation_objective_replay", "test_replay"))
            or len(training["epoch_streams"]) != training["epochs"]
            or any(replay["outputs"][name] != digest for record in (training, prediction)
                   for name, digest in record["outputs"].items())):
        raise ValueError(f"invalid Step 09 independent phase evidence: {path}")
    trace = axis_evaluate._csv(path / "learning_curve.csv")
    if tuple(trace.columns) != TRACE_COLUMNS:
        raise ValueError(f"invalid Step 09 objective trace: {path}")
    check_trace(trace, training, config)
    if not axis_evaluate._csv(path / "kl-diagnostics.csv").equals(kl_diagnostics(trace, config["model"]["n_latent"])):
        raise ValueError(f"invalid Step 09 KL trace: {path}")
    return key, (path, config, replay, training, trace)


def add_inputs(parser):
    for name, default in (("step03", axis_run.STAGE / "03_data_train"),
                          ("representations", DEFAULT_REPRESENTATIONS), ("source", DEFAULT_SOURCE)):
        parser.add_argument(f"--{name}", type=Path, default=default)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="phase", required=True)
    registration = commands.add_parser("register")
    registration.add_argument("--axis", choices=AXES, required=True)
    registration.add_argument("--fold", choices=range(3), type=int, required=True)
    registration.add_argument("--seed", choices=SEEDS, type=int, required=True)
    registration.add_argument("--bounded", action="store_true")
    for phase in PHASES:
        commands.add_parser(phase)
    for command in commands.choices.values():
        add_inputs(command)
        command.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    inputs = tuple(path.resolve() for path in (args.step03, args.representations, args.source))
    if args.phase == "register":
        result = register(args.output, (args.axis, args.fold, args.seed), args.bounded, *inputs)
    else:
        result = execute(args.output, args.phase, *inputs)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
