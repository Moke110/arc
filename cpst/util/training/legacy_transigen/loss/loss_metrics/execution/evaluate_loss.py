"""Score or verify the Step 05.1 loss comparison at the fixed-prediction boundary."""

import argparse
from fractions import Fraction
from itertools import product
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from objective import ALPHAS, COMPONENTS, FORMS, NATIVE_ALPHA, Objective
from run_loss import AXIS, DEFAULT_TRAINING, PROTOCOL, SCHEMA, SEALED_FILES, axis_run, load_step04
from tsg_axis.inputs import AxisInputs, score_fixed_prediction
from tsg_axis.metrics import (CONDITION, DISTRIBUTION, METRICS, condition_scores,
                              distribution_scores, summarize_conditions, summarize_folds)
from training_data import TENV_DEFAULT, open_tenv_training

from util.eval.fixed_prediction import read_paired_prediction, score_paired_prediction

axis_evaluate = load_step04("evaluate")
REPORT_SCHEMA = "stage2-tenv-transigen-loss-report-v1"
TRACE_COLUMNS = ("epoch", "train_loss", "validation_loss", "validation_mse",
                 *(f"{partition}_{component}" for partition in ("train", "validation") for component in COMPONENTS))


def coordinates(smoke=False):
    return tuple(product(FORMS, (NATIVE_ALPHA,) if smoke else ALPHAS,
                         (0,) if smoke else range(3)))


def check_trace(trace, training, config):
    epochs = training["epochs"]
    if (not 1 <= epochs <= config["max_epochs"]
            or trace.epoch.tolist() != list(range(1, epochs + 1))
            or not np.isfinite(trace.to_numpy()).all()):
        raise ValueError("objective trace differs from the observed stopping epoch")
    best = trace.loc[trace.validation_loss.idxmin()]
    if (best.epoch != training["best_epoch"]
            or best.validation_loss != training["validation_loss"]
            or best.validation_mse != training["validation_mse"]
            or (epochs < config["max_epochs"] and epochs - best.epoch != config["patience"])):
        raise ValueError("objective trace disagrees with the selected checkpoint")
    return int(best.epoch)


def read_run(path, smoke):
    registration, receipt = (path / name for name in ("registration.json", "manifest.json"))
    config, manifest = map(axis_run.read_json, (registration, receipt))
    args, objective = config["arguments"], config["objective"]
    key = args["loss"], args["alpha"], args["fold"]
    registration_sha = axis_run.digest(registration)
    if (type(key[2]) is not int or key not in coordinates(smoke)
            or (config["schema"], config["protocol"]) != (SCHEMA, PROTOCOL)
            or (manifest["schema"], manifest["protocol"]) != (SCHEMA, PROTOCOL)
            or args["axis"] != AXIS or manifest["axis"] != AXIS or manifest["fold"] != key[2]
            or config["scenario"] != "cold-drug" or config["seed"] != axis_run.SEED
            or manifest["seed"] != axis_run.SEED or args["smoke"] is not smoke or manifest["smoke"] is not smoke
            or config["patience"] != axis_run.PATIENCE
            or config["max_epochs"] != (2 if smoke else axis_run.MAX_EPOCHS)
            or (config["validation_seed"], config["test_seed"]) != (axis_run.VALIDATION_SEED, axis_run.TEST_SEED)
            or objective != Objective(key[0], Fraction(key[1])).identity() or manifest["objective"] != objective
            or any(manifest[name] != "exact" for name in ("test_replay", "validation_replay", "validation_objective_replay"))
            or manifest["registration_sha256"] != registration_sha or set(manifest["outputs"]) != set(SEALED_FILES)):
        raise ValueError(f"invalid sealed loss run: {path}")
    axis_run.check_outputs(path, manifest)
    training, prediction = map(axis_run.read_json, (path / "training.json", path / "predictions.json"))
    for phase, files in ((training, axis_run.TRAINING_FILES), (prediction, axis_run.PREDICTION_FILES)):
        if (phase["registration_sha256"] != registration_sha
                or phase["outputs"] != {name: manifest["outputs"][name] for name in files}):
            raise ValueError(f"phase receipt differs from its sealed run: {path}")
    if len({config["nonce"], training["process_token"], prediction["process_token"], manifest["process_token"]}) != 4:
        raise ValueError(f"run phases require independent processes: {path}")
    trace = axis_evaluate._csv(path / "learning_curve.csv")
    if tuple(trace.columns) != TRACE_COLUMNS:
        raise ValueError(f"loss component trace columns differ: {path}")
    check_trace(trace, training, config)
    return key, (path, config, manifest, training, trace)


def score_run(path, config, manifest, data, base, fold, smoke, bandwidth):
    conditions, _, paired = read_paired_prediction(path, config, manifest, data, base, AXIS, fold, smoke)
    return score_paired_prediction(conditions, paired, data.genes.gene_index.to_numpy(), bandwidth)


def summarize_loss_scores(fold_scores):
    summary = summarize_folds(fold_scores, ["loss", "alpha", "metric"])
    native = summary.loc[(summary.loss == "mse") & (summary.alpha == NATIVE_ALPHA)]
    summary["native_mse_score"] = summary.metric.map(native.set_index("metric").value)
    summary["score_difference"] = summary.value - summary.native_mse_score
    return summary


def build_report(candidate, paths, training_path, bandwidth_path, smoke):
    runs = {}
    for path in paths:
        key, run = read_run(path, smoke)
        if key in runs:
            raise ValueError(f"duplicate loss run: {key}")
        runs[key] = run
    expected = coordinates(smoke)
    if set(runs) != set(expected):
        raise ValueError(f"incomplete loss grid: missing {sorted(set(expected) - set(runs))}")
    candidate.mkdir()
    memberships, bandwidths, records, conditions, scores, epochs = {}, {}, [], [], [], []
    experiment = None
    with (open_tenv_training(training_path) as release,
          open_tenv_training(bandwidth_path) as baseline):
        hvg = axis_evaluate.open_axes(family="hvg").axis(5000)
        if (baseline.manifest["inputs"]["axis_sha256"] != hvg.sha256
                or baseline.manifest["outputs"]["genes.csv"] != hvg.sha256):
            raise ValueError("fold bandwidths require the released HVG-DEG 5k training axis")
        data, base = AxisInputs(release), AxisInputs(baseline)
        for loss, alpha, fold in expected:
            path, config, manifest, training, trace = runs[loss, alpha, fold]
            identity = {key: value for key, value in config.items() if key not in {"objective", "arguments", "nonce"}}
            common = {key: value for key, value in identity.items() if key not in {"model", "partitions", "conditions"}}
            common["model"] = {key: value for key, value in identity["model"].items() if key != "transform"}
            if experiment is None:
                experiment = common
            if common != experiment:
                raise ValueError("loss runs use different shared provenance, settings or runtime")
            if memberships.setdefault(fold, identity) != identity:
                raise ValueError(f"unmatched experiment configuration for fold {fold}")
            if fold not in bandwidths:
                bandwidths[fold] = axis_evaluate._bandwidth(base, fold)
            paired, metrics = score_run(path, config, manifest, data, base, fold, smoke, bandwidths[fold]["bandwidth"])
            key = {"loss": loss, "alpha": alpha, "fold": fold, "seed": config["seed"],
                   "run_sha256": manifest["registration_sha256"]}
            records.append({**key, "manifest_sha256": axis_run.digest(path / "manifest.json"),
                            "best_epoch": training["best_epoch"], "stopping_epoch": training["epochs"],
                            "validation_loss": training["validation_loss"],
                            "validation_mse": training["validation_mse"],
                            "checkpoint_sha256": manifest["outputs"]["checkpoint.pt"],
                            "trace_sha256": manifest["outputs"]["learning_curve.csv"]})
            conditions.append(pd.DataFrame({**key, **dict(paired)}))
            scores.append(pd.DataFrame({**key, **dict(metrics)}))
            epochs.append(pd.DataFrame({**key, "best_epoch": training["best_epoch"],
                                        "stopping_epoch": training["epochs"], **dict(trace)}))
    fold_scores = pd.concat(scores, ignore_index=True)
    tables = {
        "runs.csv": pd.DataFrame(records),
        "condition-scores.csv": pd.concat(conditions, ignore_index=True),
        "fold-scores.csv": fold_scores,
        "epoch-objectives.csv": pd.concat(epochs, ignore_index=True),
        "bandwidths.csv": pd.DataFrame(bandwidths.values()),
        "delta-metrics.csv": summarize_loss_scores(fold_scores),
    }
    for name, table in tables.items():
        table.to_csv(candidate / name, index=False)
    axis_run.write_json(candidate / "manifest.json", {
        "status": "PASS", "schema": REPORT_SCHEMA, "axis": AXIS, "seed": axis_run.SEED,
        "smoke": smoke, "runs": len(runs), "metrics": list(METRICS),
        "ordered_populations_verified": True, "selected_checkpoint_and_replay_verified": True,
        "aggregation": "equal eligible conditions per fold, equal eligible folds per loss/alpha",
        "scientific_replay_tolerances": {"rtol": 0., "atol": 0.},
        "outputs": axis_run.file_hashes(candidate, tables),
        "code_sha256": {"evaluate_loss.py": axis_run.digest(Path(__file__)),
                        "Step 04 evaluate.py": axis_run.digest(Path(axis_evaluate.__file__)),
                        "Step 04 metrics.py": axis_run.digest(Path(axis_evaluate.STEP / "output/api/tsg_axis/metrics.py"))},
    })


def report_scratch(command, output, scratch):
    if scratch is not None:
        return scratch
    return output.parent if command == "score" else axis_run.REPO / "tmp"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("score", "verify"))
    parser.add_argument("--run", action="append", type=Path, required=True)
    parser.add_argument("--training", type=Path, default=DEFAULT_TRAINING)
    parser.add_argument("--bandwidth-training", type=Path, default=TENV_DEFAULT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--scratch", type=Path, help="temporary report workspace")
    args = parser.parse_args()
    with axis_run.publication_lock(args.output) as output:
        if args.command == "score" and output.exists():
            raise ValueError("score requires an unused output directory")
        staging = report_scratch(args.command, output, args.scratch)
        staging.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(prefix=f".{output.name}-", dir=staging) as scratch:
            candidate = Path(scratch) / "candidate"
            build_report(candidate, args.run, args.training, args.bandwidth_training, args.smoke)
            if args.command == "score":
                candidate.rename(output)
            else:
                manifest = axis_run.read_json(candidate / "manifest.json")
                for name in ("manifest.json", *manifest["outputs"]):
                    if axis_run.digest(output / name) != axis_run.digest(candidate / name):
                        raise ValueError(f"score replay differs: {name}")
        print(json.dumps({"status": "PASS", "manifest_sha256": axis_run.digest(output / "manifest.json")}, sort_keys=True))


if __name__ == "__main__":
    main()
