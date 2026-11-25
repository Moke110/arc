"""Score or independently verify all registered architecture comparisons."""

import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from shape_model import CONTROL, axis_run, coordinates, dimensions
from run_shape import (AXIS, DEFAULT_TRAINING, OBJECTIVE, PHASES, PHASE_FILES, PROTOCOL,
                       SCHEMA, SEALED_FILES, kl_diagnostics)
from evaluate_loss import axis_evaluate, check_trace, score_run
from evaluate_sweep import TRACE_COLUMNS
from tsg_axis.inputs import AxisInputs
from tsg_axis.metrics import METRICS, summarize_folds
from training_data import TENV_DEFAULT, open_tenv_training

REPORT_SCHEMA = "stage2-tenv-transigen-shape-round1-report-v1"


def read_run(path, smoke):
    registration = path / "registration.json"
    config = axis_run.read_json(registration)
    args = config["arguments"]
    key = args["shape"], args["fold"]
    registration_sha = axis_run.digest(registration)
    if (type(key[1]) is not int or key not in coordinates(smoke)
            or (config["schema"], config["protocol"], config["shape"]) != (SCHEMA, PROTOCOL, key[0])
            or (args["axis"], args["device"], args["smoke"]) != (AXIS, "cuda", smoke)
            or config["scenario"] != "cold-drug" or config["seed"] != axis_run.SEED
            or (config["validation_seed"], config["test_seed"]) != (axis_run.VALIDATION_SEED, axis_run.TEST_SEED)
            or config["patience"] != axis_run.PATIENCE
            or config["max_epochs"] != (2 if smoke else axis_run.MAX_EPOCHS)
            or config["objective"] != OBJECTIVE.identity()
            or any(config["model"][name] != value for name, value in dimensions(key[0]).items())):
        raise ValueError(f"invalid shape registration: {path}")
    training, prediction, replay = (axis_run.verified_receipt(path, PHASES[phase], registration_sha, PHASE_FILES[phase])
                                    for phase in PHASES)
    if ((replay["schema"], replay["protocol"], replay["shape"], replay["fold"]) != (SCHEMA, PROTOCOL, *key)
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


def pair_scores(fold_scores):
    baseline = fold_scores.loc[fold_scores["shape"] == CONTROL, ["fold", "metric", "value"]]
    pairs = fold_scores.merge(baseline.rename(columns={"value": "control_score"}),
                              on=["fold", "metric"], how="left", validate="many_to_one", indicator=True)
    if not pairs["_merge"].eq("both").all():
        raise ValueError("every score needs the fresh original architecture on the same fold and metric")
    pairs = pairs.drop(columns="_merge")
    pairs["score_difference"] = pairs.value - pairs.control_score
    summary = summarize_folds(pairs, ["shape", "metric"])
    difference = pairs.groupby(["shape", "metric"], sort=False).score_difference.agg(["mean", "count"])
    summary = summary.merge(difference.rename(columns={"mean": "score_difference", "count": "n_paired_folds"}),
                            on=["shape", "metric"], validate="one_to_one")
    return pairs, summary


def shared_configuration(config):
    shape_fields = {*dimensions(CONTROL), "constructor_order", "initialization"}
    return {name: value for name, value in config.items() if name not in {"nonce", "arguments", "shape", "model"}} | {
        "model": {name: value for name, value in config["model"].items() if name not in shape_fields}}


def build_report(candidate, runs, training, bandwidth_training, smoke):
    expected = coordinates(smoke)
    verified = {}
    for path in runs:
        key, detail = read_run(path, smoke)
        if key in verified:
            raise ValueError(f"duplicate shape run {key}")
        verified[key] = detail
    if verified.keys() != dict.fromkeys(expected).keys():
        raise ValueError(f"incomplete shape comparison: missing {sorted(set(expected) - verified.keys())}")
    common_run, fold_configs, bandwidths = None, {}, {}
    runs_table, conditions, fold_metrics, kl_tables = [], [], [], []
    with (open_tenv_training(training) as release, open_tenv_training(bandwidth_training) as baseline):
        hvg = axis_evaluate.open_axes(family="hvg").axis(5000)
        if (baseline.manifest["inputs"]["axis_sha256"] != hvg.sha256
                or baseline.manifest["outputs"]["genes.csv"] != hvg.sha256):
            raise ValueError("fold bandwidths require the released HVG-DEG 5k training axis")
        data, base = AxisInputs(release), AxisInputs(baseline)
        for shape, fold in expected:
            path, config, manifest, training_record, kl = verified[shape, fold]
            common = shared_configuration(config)
            shared_model = {name: value for name, value in common["model"].items() if name != "transform"}
            shared = {name: value for name, value in common.items() if name not in {"model", "partitions", "conditions"}}
            shared["model"] = shared_model
            if common_run is None:
                common_run = shared
            if common_run != shared or fold_configs.setdefault(fold, common) != common:
                raise ValueError(f"shape comparison has mismatched protocol or populations: {shape}, fold {fold}")
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
    args = parser.parse_args()
    output = args.output.resolve()
    scratch = args.scratch or (output.parent if args.command == "score" else axis_run.REPO / "tmp")
    with axis_run.publication_lock(output) as target, TemporaryDirectory(prefix=".shape-score-", dir=scratch) as directory:
        candidate = Path(directory) / "report"
        build_report(candidate, args.run, args.training, args.bandwidth_training, args.smoke)
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
