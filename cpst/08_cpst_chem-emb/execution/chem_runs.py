"""Step 08 registration, fitting and independent prediction replay."""

import argparse
from contextlib import contextmanager
from dataclasses import dataclass
from hashlib import sha256
from itertools import product
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import torch

STEP = Path(__file__).resolve().parents[1]
STAGE = STEP.parent
ROUND_FOUR = STAGE / "06_tsg_module/06.4_round4"
sys.path[:0] = [str(STAGE), str(ROUND_FOUR / "execution")]

from util.training.cpst_model import CPSTAdapter, PROJECTION_WIDTH
from run_measured import AXIS, DEFAULT_REPRESENTATIONS, DEFAULT_SOURCE, DEFAULT_TRAINING
from run_measured import PHASES, SEALED_FILES, PHASE_FILES, axis_run, bind as bind_measured
from evaluate_loss import axis_evaluate, check_trace
from evaluate_sweep import TRACE_COLUMNS
from util.training.transigen import fit_with_diagnostics, kl_diagnostics
from measured_roster import MEASURED_SHARED, SEEDS
from cpst_util.representations import ARMS, WIDTHS
from util.training.cpst_reference import alpha_release_digest, verify_release
from util.training import registered
from util.training.chemical import aligned_features
from util.training.provenance import verify_replay_configuration, verify_sealed_code

SCHEMA = "cpst-alpha-chemical-run-v1"
PROTOCOL = "tenv-matched-control-cpst-alpha-chemical-v1"
REFERENCE = "KPGT"
ORIGINAL_CODE = "ac62a177a54aa99e6ac626728af12394f3e0d940"


def coordinates(smoke=False):
    """Return the preregistered representation × cold-drug fold × seed roster."""
    return tuple(product(ARMS, (0,) if smoke else range(3), (17,) if smoke else SEEDS))


@dataclass(frozen=True)
class ChemicalModelInput:
    control: object
    expression: object
    doses: object
    features: np.ndarray

    def fingerprints(self, rows):
        return self.features[rows]


class ChemicalInputs:
    """Reuse the frozen matched-control and expression readers for each arm."""

    def __init__(self, data, vectors, arm):
        self.data = data
        self.features = aligned_features(vectors, arm, data.parents)

    def __getattr__(self, name):
        return getattr(self.data, name)

    def model_input(self, phase):
        base = self.data.model_input(phase)
        return ChemicalModelInput(base.control, base.expression, base.doses, self.features)


@contextmanager
def open_inputs(training, representations):
    with axis_run.open_inputs(training, representations) as (data, vectors):
        accepted_parents = set(vectors.parents.drug_parent_inchikey)
        if len(data.parents) != 1383 or len(accepted_parents) != 459 or set(data.parents) != accepted_parents:
            raise ValueError("shared tenv releases must cover exactly 1,383 conditions and 459 parents")
        all_rows = np.arange(len(data.parents))
        for fold in range(3):
            partitions = [data.rows(fold, phase) for phase in axis_run.PARTITIONS]
            if not np.array_equal(np.sort(np.concatenate(partitions)), all_rows):
                raise ValueError(f"frozen cold-drug fold {fold} does not partition all pooled conditions")
        if not np.array_equal(aligned_features(vectors, REFERENCE, data.parents),
                              data.model_input("test").fingerprints(all_rows)):
            raise ValueError("ordered KPGT vectors differ from the alpha")
        yield data, vectors


def arguments(coordinate, smoke):
    arm, fold, seed = coordinate
    if type(fold) is not int or type(seed) is not int or coordinate not in coordinates(smoke):
        raise ValueError(f"unregistered Step 08 coordinate: {coordinate}")
    return {"arm": arm, "fold": fold, "seed": seed, "smoke": smoke, "axis": AXIS, "device": "cuda"}


def bind(coordinate, data, vectors, source, smoke):
    args = arguments(coordinate, smoke)
    arm = args["arm"]
    config = bind_measured({**args, "shape": MEASURED_SHARED, "source": str(source)}, data, vectors)
    width = vectors.width(arm)
    config.update(schema=SCHEMA, protocol=PROTOCOL, alpha_release_sha256=alpha_release_digest(),
                  representation={"arm": arm, "width": width, "projection_width": PROJECTION_WIDTH,
                                  "manifest_sha256": vectors.manifest_sha256,
                                  "vectors_sha256": vectors.manifest["files"][f"{arm}.npy"]})
    model = config["model"]
    model.update(molecule_feature=arm, features_dim=width, representation_width=width)
    if arm != REFERENCE:
        model["constructor_order"].append("replace the first chemical projection with isolated RNG")
        model["initialization"] = "; ".join(model["constructor_order"])
        model["adaptations"] = [
            f"frozen {arm} {width}-value input" if item == "frozen KPGT2304 input" else item
            for item in model["adaptations"]
        ]
    config["code_sha256"][Path(__file__).relative_to(axis_run.REPO).as_posix()] = axis_run.digest(Path(__file__))
    return axis_run.transigen.bind_experiment(args, config)


def register(output, coordinate, training, representations, source, smoke):
    verify_release()

    def definition():
        with open_inputs(training, representations) as (data, vectors):
            config = bind(coordinate, data, vectors, source, smoke)
        return {**config, "arguments": arguments(coordinate, smoke)}

    return registered.register(output, definition)


def initialization(model, arm):
    projection = model.feat_embeddings[0]
    if (projection.in_features, projection.out_features) != (WIDTHS[arm], PROJECTION_WIDTH):
        raise ValueError(f"alpha uses a noncanonical chemical projection for {arm}")
    downstream, projection_hash = sha256(), sha256()
    for name, tensor in sorted(model.state_dict().items()):
        digest = projection_hash if name.startswith("feat_embeddings.0.") else downstream
        digest.update(name.encode() + b"\0" + tensor.detach().cpu().numpy().tobytes())
    return {
        "downstream_sha256": downstream.hexdigest(),
        "projection_sha256": projection_hash.hexdigest(),
        "torch_cpu_rng_sha256": sha256(torch.get_rng_state().numpy().tobytes()).hexdigest(),
        "torch_cuda_rng_sha256": sha256(torch.cuda.get_rng_state().numpy().tobytes()).hexdigest(),
        "projection_shape": [projection.in_features, projection.out_features],
    }


def execute(output, phase, training, representations, source):
    verify_release()
    with registered.phase(output, phase, PHASE_FILES) as run:
        config = run.config
        args = config["arguments"]
        coordinate = args["arm"], args["fold"], args["seed"]
        if args != arguments(coordinate, args["smoke"]):
            raise ValueError("run arguments differ from its coordinate")
        os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
        torch.use_deterministic_algorithms(True)
        axis_run.kpgt.configure_backend("cuda")
        with open_inputs(training, representations) as (data, vectors):
            identity = {key: value for key, value in config.items() if key not in {"arguments", "nonce"}}
            current = bind(coordinate, data, vectors, source, args["smoke"])
            verify_replay_configuration(identity, current, replay=run.sealed is not None,
                                        repo=axis_run.REPO, revision=ORIGINAL_CODE)
            arm, _, seed = coordinate
            adapter = CPSTAdapter(source, config["model"]["transform"], chemical_width=WIDTHS[arm])
            model = adapter.build(seed, "cuda", n_genes=len(data.genes))
            model_data = ChemicalInputs(data, vectors, arm)
            if phase == "train":
                paired = initialization(model, arm)
                details = fit_with_diagnostics(run.scratch, adapter, model, model_data, config, run.registration_sha)
                details["initialization"] = paired
            else:
                values, rows = axis_run.predict_test(run.root, adapter, model, model_data, config,
                                                     run.registration_sha, run.training)
                details = axis_run.record_predictions(run.root, run.scratch, values, rows, model_data, config,
                                                      replay=phase == "replay")
        return run.publish(details, {**args, "schema": SCHEMA, "protocol": PROTOCOL, "objective": config["objective"]})


def read_run(path, smoke):
    registration = path / "registration.json"
    config = axis_run.read_json(registration)
    verify_sealed_code(config, axis_run.REPO, ORIGINAL_CODE)
    args = config["arguments"]
    key = args["arm"], args["fold"], args["seed"]
    registration_sha = axis_run.digest(registration)
    if (args != arguments(key, smoke) or config["schema"] != SCHEMA or config["protocol"] != PROTOCOL
            or config["seed"] != args["seed"] or config["alpha_release_sha256"] != alpha_release_digest()):
        raise ValueError(f"invalid sealed Step 08 registration: {path}")
    training, prediction, replay = (
        axis_run.verified_receipt(path, PHASES[phase], registration_sha, PHASE_FILES[phase]) for phase in PHASES
    )
    if (len({config["nonce"], training["process_token"], prediction["process_token"], replay["process_token"]}) != 4
            or any(record.get(field) != value for record in (training, prediction, replay)
                   for field, value in {**args, "schema": SCHEMA, "protocol": PROTOCOL,
                                       "objective": config["objective"]}.items())
            or any(record["cuda_peak_allocated_bytes"] <= 0 for record in (training, prediction, replay))
            or any(replay.get(field) != "exact" for field in
                   ("validation_replay", "validation_objective_replay", "test_replay"))
            or any(replay["outputs"][name] != digest for record in (training, prediction)
                   for name, digest in record["outputs"].items())):
        raise ValueError(f"invalid Step 08 independent phase evidence: {path}")
    trace = axis_evaluate._csv(path / "learning_curve.csv")
    if tuple(trace.columns) != TRACE_COLUMNS:
        raise ValueError(f"invalid Step 08 learning trace: {path}")
    check_trace(trace, training, config)
    if not axis_evaluate._csv(path / "kl-diagnostics.csv").equals(kl_diagnostics(trace, config["model"]["n_latent"])):
        raise ValueError(f"invalid Step 08 KL trace: {path}")
    return key, (path, config, replay, training, trace)


def add_inputs(parser):
    for name, default in (("training", DEFAULT_TRAINING), ("representations", DEFAULT_REPRESENTATIONS),
                          ("source", DEFAULT_SOURCE)):
        parser.add_argument(f"--{name}", type=Path, default=default)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    phases = parser.add_subparsers(dest="phase", required=True)
    registration = phases.add_parser("register")
    registration.add_argument("--arm", choices=ARMS, required=True)
    registration.add_argument("--fold", choices=range(3), type=int, required=True)
    registration.add_argument("--seed", choices=SEEDS, type=int, required=True)
    registration.add_argument("--smoke", action="store_true")
    for phase in PHASES:
        phases.add_parser(phase)
    for command in phases.choices.values():
        add_inputs(command)
        command.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    inputs = tuple(path.resolve() for path in (args.training, args.representations, args.source))
    if args.phase == "register":
        digest = register(args.output.resolve(), (args.arm, args.fold, args.seed), *inputs, args.smoke)
        result = {"registration_sha256": digest}
    else:
        result = execute(args.output.resolve(), args.phase, *inputs)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
