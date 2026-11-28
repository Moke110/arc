"""Score or independently verify all registered architecture comparisons."""

import argparse
from functools import partial
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from shape_model import STEP, axis_run
from shape_roster import CONTROL, PRIOR_NAMES, coordinates, dimensions
from prior_shape import PRIOR_MATRIX, PRIOR_PROTOCOL, PRIOR_RECEIPT, PRIOR_RUN_SCHEMA, load_prior
from run_shape import (AXIS, DEFAULT_TRAINING, OBJECTIVE, PHASES, PHASE_FILES, PROTOCOL,
                       SCHEMA, SEALED_FILES)
from evaluate_loss import axis_evaluate, check_trace, score_run
from evaluate_sweep import TRACE_COLUMNS
from tsg_axis.inputs import AxisInputs
from tsg_axis.metrics import METRICS
from training_data import TENV_DEFAULT, open_tenv_training
from util.eval.shape import pair_shape_scores
from util.training.transigen import kl_diagnostics

pair_scores = partial(pair_shape_scores, control=CONTROL)
REPORT_SCHEMA = "stage2-tenv-transigen-shape-width-extension-report-v1"


def read_run(path, smoke, *, prior=False):
    registration = path / "registration.json"
    config = axis_run.read_json(registration)
    args = config["arguments"]
    source_shape, fold = args["shape"], args["fold"]
    key = (PRIOR_NAMES.get(source_shape) if prior else source_shape), fold
    schema, protocol = (PRIOR_RUN_SCHEMA, PRIOR_PROTOCOL) if prior else (SCHEMA, PROTOCOL)
    registration_sha = axis_run.digest(registration)
    if (type(key[1]) is not int or key not in coordinates(smoke)
            or (config["schema"], config["protocol"], config["shape"]) != (schema, protocol, source_shape)
            or (args["axis"], args["device"], args["smoke"]) != (AXIS, "cuda", smoke)
            or config["scenario"] != "cold-drug" or config["seed"] != axis_run.SEED
            or (config["validation_seed"], config["test_seed"]) != (axis_run.VALIDATION_SEED, axis_run.TEST_SEED)
            or config["patience"] != axis_run.PATIENCE
            or config["max_epochs"] != (2 if smoke else axis_run.MAX_EPOCHS)
            or config["objective"] != OBJECTIVE.identity()
            or any(config["model"][name] != value for name, value in dimensions(key[0]).items())):
        raise ValueError(f"invalid shape registration: {path}")
    model = config["model"]
    replacement = f"replace both predictive heads Linear{model['predictive_head_dim']}→{model['n_latent']} with PyTorch defaults"
    constructor = ["native constructor and selective Xavier", replacement]
    if model["constructor_order"] != constructor or model["initialization"] != "; ".join(constructor):
        raise ValueError(f"invalid shape initialization: {path}")
    training, prediction, replay = (axis_run.verified_receipt(path, PHASES[phase], registration_sha, PHASE_FILES[phase])
                                    for phase in PHASES)
    if ((replay["schema"], replay["protocol"], replay["shape"], replay["fold"]) != (schema, protocol, source_shape, fold)
            or (replay["axis"], replay["seed"], replay["smoke"]) != (AXIS, config["seed"], smoke)
            or replay["objective"] != config["objective"]
            or any(replay[name] != "exact" for name in ("test_replay", "validation_replay", "validation_objective_replay"))
            or training["parameters"] < 1 or training["cuda_peak_allocated_bytes"] < 1
            or set(replay["outputs"]) != set(SEALED_FILES)):
        raise ValueError(f"invalid shape replay or training receipt: {path}")
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
    return key, (path, config, replay, training, kl)


def shared_configuration(config):
    shape_fields = {*dimensions(CONTROL), "constructor_order", "initialization"}
    release_fields = {"nonce", "arguments", "shape", "schema", "protocol", "model", "code_sha256"}
    execution = (STEP / "execution").relative_to(axis_run.REPO)
    shape_code = {(execution / name).as_posix() for name in ("run_shape.py", "shape_model.py", "shape_roster.py")}
    return {name: value for name, value in config.items() if name not in release_fields} | {
        "model": {name: value for name, value in config["model"].items() if name not in shape_fields},
        "code_sha256": {name: digest for name, digest in config["code_sha256"].items()
                        if name not in shape_code},
    }


def validate_coordinates(expected, paths):
    if len(paths) != len(set(paths)):
        raise ValueError("duplicate shape coordinate")
    missing, unexpected = set(expected) - set(paths), set(paths) - set(expected)
    if missing or unexpected:
        raise ValueError(f"incomplete shape comparison: missing {sorted(missing)}, unexpected {sorted(unexpected)}")


def verify_configurations(verified):
    common_run, fold_configs = None, {}
    for (shape, fold), (_, config, _, _, _) in verified.items():
        common = shared_configuration(config)
        shared = {name: value for name, value in common.items() if name not in {"partitions", "conditions"}}
        shared["model"] = {name: value for name, value in common["model"].items() if name != "transform"}
        if common_run is None:
            common_run = shared
        if common_run != shared or fold_configs.setdefault(fold, common) != common:
            raise ValueError(f"shape comparison has mismatched protocol or populations: {shape}, fold {fold}")


def verify_prior_scores(prior_report, fold_scores):
    source = pd.read_csv(prior_report / "fold-scores.csv", float_precision="round_trip")
    source["shape"] = source["shape"].map(PRIOR_NAMES)
    old_shapes = set(PRIOR_NAMES.values())
    current = fold_scores.loc[fold_scores["shape"].isin(old_shapes)].reset_index(drop=True)
    if not current.equals(source):
        raise ValueError("reproduced historical scores or eligibility differ from the sealed #347 report")


def verify_historical_runs(paths):
    historical = dict(read_run(path, False, prior=True) for path in paths.values())
    if historical.keys() != paths.keys():
        raise ValueError("historical run identifiers differ from the sealed matrix")
    return historical


def build_report(candidate, runs, training, bandwidth_training, smoke, prior_matrix=None, prior_receipt=PRIOR_RECEIPT):
    expected = coordinates(smoke)
    prior, _, prior_provenance = load_prior(prior_matrix, prior_receipt) if prior_matrix else ({}, {}, None)
    verified_current = [read_run(path, smoke) for path in runs]
    validate_coordinates(expected, [*prior, *(key for key, _ in verified_current)])
    historical = verify_historical_runs(prior)
    verified = historical | dict(verified_current)
    verify_configurations(verified)
    bandwidths = {}
    runs_table, conditions, fold_metrics, kl_tables = [], [], [], []
    with (open_tenv_training(training) as release, open_tenv_training(bandwidth_training) as baseline):
        hvg = axis_evaluate.open_axes(family="hvg").axis(5000)
        if (baseline.manifest["inputs"]["axis_sha256"] != hvg.sha256
                or baseline.manifest["outputs"]["genes.csv"] != hvg.sha256):
            raise ValueError("fold bandwidths require the released HVG-DEG 5k training axis")
        data, base = AxisInputs(release), AxisInputs(baseline)
        for shape, fold in expected:
            path, config, manifest, training_record, kl = verified[shape, fold]
            if fold not in bandwidths:
                bandwidths[fold] = axis_evaluate._bandwidth(base, fold)
            paired, metrics = score_run(path, config, manifest, data, base, fold, smoke, bandwidths[fold]["bandwidth"])
            key = dict(shape=shape, fold=fold, seed=config["seed"], run_sha256=manifest["registration_sha256"])
            runs_table.append({**key, "manifest_sha256": axis_run.digest(path / "manifest.json"),
                               "encoder": ",".join(map(str, config["model"]["encoder"])),
                               "decoder": ",".join(map(str, config["model"]["decoder"])),
                               "n_latent": config["model"]["n_latent"], "trainable_parameters": training_record["parameters"],
                               "best_epoch": training_record["best_epoch"], "stopping_epoch": training_record["epochs"],
                               "training_peak_rss_kib": training_record["peak_rss_kib"],
                               "training_cuda_peak_allocated_bytes": training_record["cuda_peak_allocated_bytes"]})
            for table, frame in ((conditions, paired), (fold_metrics, metrics), (kl_tables, kl)):
                table.append(pd.DataFrame({**key, **dict(frame)}))
    fold_scores, summary = pair_scores(pd.concat(fold_metrics, ignore_index=True))
    if prior_matrix:
        verify_prior_scores(prior_matrix.parent / "report", fold_scores)
    candidate.mkdir()
    tables = {
        "runs.csv": pd.DataFrame(runs_table),
        "condition-scores.csv": pd.concat(conditions, ignore_index=True),
        "fold-scores.csv": fold_scores,
        "delta-metrics.csv": summary,
        "bandwidths.csv": pd.DataFrame(bandwidths.values()),
        "kl-diagnostics.csv": pd.concat(kl_tables, ignore_index=True),
    }
    for name, table in tables.items():
        table.to_csv(candidate / name, index=False)
    axis_run.write_json(candidate / "manifest.json", {
        "status": "PASS", "schema": REPORT_SCHEMA, "axis": AXIS, "seed": axis_run.SEED,
        "smoke": smoke, "runs": len(verified), "fold_scores": len(fold_scores), "metrics": list(METRICS),
        "fold_reference": CONTROL, "selected_checkpoint_and_replay_verified": True,
        "aggregation": "equal eligible conditions per fold, equal eligible folds per shape, differences paired within fold",
        "scientific_replay_tolerances": {"rtol": 0., "atol": 0.},
        "prior_release": prior_provenance,
        "outputs": axis_run.file_hashes(candidate, tables),
        "code_sha256": {path.name: axis_run.digest(path) for path in (
            Path(__file__), Path(axis_evaluate.__file__), Path(axis_evaluate.STEP / "output/api/tsg_axis/metrics.py"))},
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("score", "verify"))
    parser.add_argument("--run", action="append", type=Path, required=True)
    parser.add_argument("--training", type=Path, default=DEFAULT_TRAINING)
    parser.add_argument("--bandwidth-training", type=Path, default=TENV_DEFAULT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--scratch", type=Path)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--prior-matrix", type=Path, default=PRIOR_MATRIX)
    parser.add_argument("--prior-receipt", type=Path, default=PRIOR_RECEIPT)
    args = parser.parse_args()
    output = args.output.resolve()
    scratch = args.scratch or (output.parent if args.command == "score" else axis_run.REPO / "tmp")
    with axis_run.publication_lock(output) as target, TemporaryDirectory(prefix=".shape-score-", dir=scratch) as directory:
        candidate = Path(directory) / "report"
        prior_matrix = None if args.smoke else args.prior_matrix.resolve()
        build_report(candidate, args.run, args.training, args.bandwidth_training, args.smoke,
                     prior_matrix, args.prior_receipt)
        if args.command == "verify":
            receipt = axis_run.read_json(candidate / "manifest.json")
            if axis_run.read_json(target / "manifest.json") != receipt:
                raise ValueError("reproduced shape score receipt differs")
            axis_run.check_outputs(target, receipt)
        else:
            if target.exists():
                raise ValueError("shape score candidate requires an unused output directory")
            candidate.rename(target)
    print(json.dumps({"status": "PASS", "report_sha256": axis_run.digest(output / "manifest.json")}, sort_keys=True))


if __name__ == "__main__":
    main()
