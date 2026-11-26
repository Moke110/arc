"""Verify the registration, objective trace and independent replay receipts of a fit."""

from composition_model import axis_run, objective_identity
from composition_roster import NO_X2_SHAPES, fit_coordinates
from run_composition import AXIS, PHASES, PHASE_FILES, PROTOCOL, SCHEMA, SEALED_FILES, architecture_config
from run_shape import kl_diagnostics
from evaluate_loss import axis_evaluate, check_trace
from evaluate_sweep import TRACE_COLUMNS


def read_run(path, smoke):
    registration = path / "registration.json"
    config = axis_run.read_json(registration)
    args = config["arguments"]
    shape, fold, seed = (args[field] for field in ("shape", "fold", "seed"))
    key = shape, fold, seed
    registration_sha = axis_run.digest(registration)
    if (type(fold) is not int or type(seed) is not int or key not in fit_coordinates(smoke)
            or (config["schema"], config["protocol"], config["shape"]) != (SCHEMA, PROTOCOL, shape)
            or (args["axis"], args["device"], args["smoke"]) != (AXIS, "cuda", smoke)
            or config["scenario"] != "cold-drug" or config["seed"] != seed
            or (config["validation_seed"], config["test_seed"]) != (axis_run.VALIDATION_SEED, axis_run.TEST_SEED)
            or config["patience"] != axis_run.PATIENCE
            or config["max_epochs"] != (2 if smoke else axis_run.MAX_EPOCHS)
            or config["objective"] != objective_identity(shape)
            or any(config["model"][name] != value
                   for name, value in architecture_config(axis_run.kpgt.configuration(config["model"]["transform"]),
                                                          shape).items())):
        raise ValueError(f"invalid Step 06.3 registration: {path}")
    training, prediction, replay = (axis_run.verified_receipt(path, PHASES[phase], registration_sha, PHASE_FILES[phase])
                                    for phase in PHASES)
    if ((replay["schema"], replay["protocol"], replay["shape"], replay["fold"]) != (SCHEMA, PROTOCOL, shape, fold)
            or (replay["axis"], replay["seed"], replay["smoke"]) != (AXIS, config["seed"], smoke)
            or replay["objective"] != config["objective"]
            or any(replay[name] != "exact" for name in ("test_replay", "validation_replay", "validation_objective_replay"))
            or training["parameters"] < 1 or training["cuda_peak_allocated_bytes"] < 1
            or set(replay["outputs"]) != set(SEALED_FILES)):
        raise ValueError(f"invalid Step 06.3 replay or training receipt: {path}")
    for phase in (training, prediction):
        if any(replay["outputs"][name] != value for name, value in phase["outputs"].items()):
            raise ValueError(f"phase outputs differ from the independent replay receipt: {path}")
    if len({config["nonce"], training["process_token"], prediction["process_token"], replay["process_token"]}) != 4:
        raise ValueError(f"run phases require independent processes: {path}")
    trace = axis_evaluate._csv(path / "learning_curve.csv")
    if tuple(trace.columns) != TRACE_COLUMNS:
        raise ValueError(f"missing Step 05.2 objective components: {path}")
    check_trace(trace, training, config)
    kl = axis_evaluate._csv(path / "kl-diagnostics.csv")
    if not kl.equals(kl_diagnostics(trace, config["model"]["n_latent"])):
        raise ValueError(f"KL diagnostics differ from the objective trace: {path}")
    if shape in NO_X2_SHAPES and any(trace[f"{partition}_treatment_mse"].ne(0).any()
                              for partition in ("train", "validation")):
        raise ValueError(f"no-X2 arm contains a treatment reconstruction term: {path}")
    return key, (path, config, replay, training, trace)
