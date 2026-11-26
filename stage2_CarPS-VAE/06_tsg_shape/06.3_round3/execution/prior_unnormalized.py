"""Verify unchanged arms from the sealed fits of the first Round 3 attempt."""

from itertools import product
from pathlib import Path

from composition_model import STEP, axis_run
from composition_roster import ROUND_TWO_SHAPES, SHAPES
from architecture_model import objective_identity
from evaluate_loss import axis_evaluate, check_trace
from evaluate_sweep import TRACE_COLUMNS
from run_architecture import architecture_config
from run_shape import kl_diagnostics
from run_sweep import PHASES
from verify_composition import AXIS, PHASE_FILES, SEALED_FILES

SOURCE_ROOT = STEP / "output/formal/prior-unnormalized"
SOURCE_REVISION = "7ce9e6259c6554072a7193b2830ce9fd03c54248"
SOURCE_SCHEMA = "stage2-tenv-transigen-architecture-composition-run-v1"
SOURCE_PROTOCOL = "matched-plate-control-transigen-architecture-composition-v1"
SOURCE_EXECUTION = "stage2/06_tsg_shape/06.3_round3/execution"
BLOCKED_RECEIPT = STEP / "log/356-blocked.json"
BLOCKED_RECEIPT_SHA256 = "5a60af60eaa99226929f56ed42d59c889740972f642ff9524473247f900c72bd"
SOURCE_CODE = {
    f"{SOURCE_EXECUTION}/run_composition.py": "56f00520b1767ead108c46ef963b41332dea91460f02ee2c89c4289adf843fdc",
    f"{SOURCE_EXECUTION}/composition_model.py": "d8b78226c9ed2e9b8e8b0b3d5ad945bb42990dcfc510e61b517638cd156a8983",
    f"{SOURCE_EXECUTION}/composition_roster.py": "5ce7bff7ef9b0d47e0999a35c45a4633704861bfe2af0351645415eed8a2b335",
}
RETAINED = tuple(product(SHAPES[:3], range(3), (17, 2026, 43)))


def select_reusable(completed):
    by_key = {tuple(fit["coordinate"]): fit for fit in completed}
    if len(by_key) != len(completed) or not set(RETAINED) <= by_key.keys():
        raise ValueError("historical ledger repeats a coordinate or lacks a retained fit")
    return {key: by_key[key] for key in RETAINED}


def verify_fit(path, key, fit):
    shape, fold, seed = key
    registration = path / "registration.json"
    registration_sha = axis_run.digest(registration)
    config = axis_run.read_json(registration)
    args = config["arguments"]
    prior_shape = ROUND_TWO_SHAPES[shape]
    old_model = architecture_config(axis_run.kpgt.configuration(config["model"]["transform"]), prior_shape)
    if (registration_sha != fit["registration_sha256"]
            or (args["shape"], args["fold"], args["seed"]) != key
            or (args["axis"], args["device"], args["smoke"]) != (AXIS, "cuda", False)
            or (config["schema"], config["protocol"], config["shape"]) != (SOURCE_SCHEMA, SOURCE_PROTOCOL, shape)
            or config["scenario"] != "cold-drug" or config["seed"] != seed
            or (config["validation_seed"], config["test_seed"]) != (axis_run.VALIDATION_SEED, axis_run.TEST_SEED)
            or (config["patience"], config["max_epochs"]) != (axis_run.PATIENCE, axis_run.MAX_EPOCHS)
            or config["objective"] != objective_identity(prior_shape)
            or any(config["model"][field] != value for field, value in old_model.items())
            or any(config["code_sha256"].get(name) != value for name, value in SOURCE_CODE.items())):
        raise ValueError(f"retained registration does not match the original Round 3 fit: {key}")
    manifest_sha = axis_run.digest(path / "manifest.json")
    if manifest_sha != fit["run_manifest_sha256"]:
        raise ValueError(f"retained replay receipt differs from the blocked ledger: {key}")
    train, prediction, replay = (axis_run.verified_receipt(path, PHASES[phase], registration_sha, PHASE_FILES[phase])
                                 for phase in PHASES)
    if (len({config["nonce"], train["process_token"], prediction["process_token"], replay["process_token"]}) != 4
            or (replay["shape"], replay["fold"], replay["seed"]) != key
            or (replay["schema"], replay["protocol"], replay["axis"], replay["smoke"])
            != (SOURCE_SCHEMA, SOURCE_PROTOCOL, AXIS, False)
            or replay["objective"] != config["objective"]
            or any(replay[name] != "exact" for name in ("test_replay", "validation_replay", "validation_objective_replay"))
            or set(replay["outputs"]) != set(SEALED_FILES)
            or any(replay["outputs"][name] != value
                   for phase in (train, prediction) for name, value in phase["outputs"].items())
            or train["peak_rss_kib"] != fit["train_peak_rss_kib"]
            or train["cuda_peak_allocated_bytes"] != fit["train_cuda_peak_allocated_bytes"]):
        raise ValueError(f"retained fit requires separately sealed exact replay: {key}")
    trace = axis_evaluate._csv(path / "learning_curve.csv")
    if tuple(trace.columns) != TRACE_COLUMNS:
        raise ValueError(f"retained fit objective trace schema differs: {key}")
    check_trace(trace, train, config)
    if not axis_evaluate._csv(path / "kl-diagnostics.csv").equals(kl_diagnostics(trace, config["model"]["n_latent"])):
        raise ValueError(f"retained fit KL diagnostics differ: {key}")
    if shape == SHAPES[2] and any(trace[f"{phase}_treatment_mse"].ne(0).any() for phase in ("train", "validation")):
        raise ValueError(f"retained no-X2 trace includes treated reconstruction: {key}")
    return (path, config, replay, train, trace), {
        **fit, "shape": shape, "fold": fold, "seed": seed, "source_shape": shape,
        "best_epoch": train["best_epoch"], "stopping_epoch": train["epochs"],
        "train_wall_seconds": train["wall_seconds"], "source_ticket": 356,
        "path": f"matrix/runs/{shape}/fold-{fold}/seed-{seed}",
    }


def load_unnormalized(root=SOURCE_ROOT):
    if axis_run.digest(BLOCKED_RECEIPT) != BLOCKED_RECEIPT_SHA256:
        raise ValueError("the first Round 3 attempt receipt changed")
    blocked = axis_run.read_json(BLOCKED_RECEIPT)
    if blocked["matrix_execution_revision"] != SOURCE_REVISION:
        raise ValueError("retained Round 3 execution revision changed")
    root = Path(root)
    partial = root / "partial-status.json"
    preflight = root / "matrix/preflight.json"
    partial_sha, preflight_sha = axis_run.digest(partial), axis_run.digest(preflight)
    if (partial_sha, preflight_sha) != (blocked["partial_ledger_sha256"], blocked["preflight_sha256"]):
        raise ValueError("the first Round 3 attempt preflight or partial ledger changed")
    status = axis_run.read_json(partial)
    if status["status"] != "BLOCKED" or len(status["completed"]) != blocked["new_fits_with_exact_replay"]:
        raise ValueError("retained fits must come from the recorded partial matrix")
    verified, coordinates = {}, []
    for key, fit in select_reusable(status["completed"]).items():
        shape, fold, seed = key
        relative = f"matrix/runs/{shape}/fold-{fold}/seed-{seed}"
        verified[key], evidence = verify_fit(root / relative, key, fit)
        coordinates.append(evidence)
    return verified, {
        "ticket": 356, "matrix_code_revision": SOURCE_REVISION, "blocked_receipt_sha256": BLOCKED_RECEIPT_SHA256,
        "preflight_sha256": preflight_sha, "partial_ledger_sha256": partial_sha,
        "registered_code_sha256": SOURCE_CODE, "coordinates": coordinates,
    }
