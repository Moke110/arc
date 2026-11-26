"""Score or independently reproduce the Step 06.2 report from sealed fit receipts."""

import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from architecture_model import STEP, axis_run
from architecture_roster import CONTROL, HEAD_PARAMETERS, coordinates, objective_label
from evaluate_loss import axis_evaluate, score_run
from evaluate_shape import pair_scores, validate_coordinates
from matrix_architecture import SCHEMA as MATRIX_SCHEMA, fit_coordinates
from prior_shape import verify_matrix
from run_architecture import AXIS, DEFAULT_TRAINING, architecture_config
from training_data import TENV_DEFAULT, open_tenv_training
from tsg_axis.inputs import AxisInputs
from tsg_axis.metrics import METRICS
from verify_architecture import read_run

REPORT_SCHEMA = "stage2-tenv-transigen-prediction-architecture-report-v1"
STAGE = STEP.parents[1]


def shared_configuration(config):
    release_fields = {"nonce", "arguments", "shape", "schema", "protocol", "objective", "model", "code_sha256"}
    architecture_fields = set(architecture_config(config["model"], CONTROL))
    provenance_prefix = "stage2/06_tsg_shape/"
    return {name: value for name, value in config.items() if name not in release_fields} | {
        "model": {name: value for name, value in config["model"].items() if name not in architecture_fields},
        "code_sha256": {name: value for name, value in config["code_sha256"].items()
                        if not name.startswith(provenance_prefix)},
    }


def verify_configurations(verified):
    common, fold_configs = None, {}
    for (shape, fold), (_, config, _, _, _) in verified.items():
        identity = shared_configuration(config)
        shared = {name: value for name, value in identity.items() if name not in {"partitions", "conditions", "model"}}
        shared["model"] = {name: value for name, value in identity["model"].items() if name != "transform"}
        if common is None:
            common = shared
        if common != shared or fold_configs.setdefault(fold, identity) != identity:
            raise ValueError(f"architecture comparison has mismatched protocol or populations: {shape}, fold {fold}")


def build_report(candidate, runs, training, bandwidth_training, smoke, reused=None):
    expected = coordinates(smoke)
    reused = {} if reused is None else reused
    fresh = [read_run(path, smoke) for path in runs]
    validate_coordinates(expected, [*reused, *(key for key, _ in fresh)])
    verified = reused | dict(fresh)
    verify_configurations(verified)
    bandwidths = {}
    run_rows, conditions, fold_metrics, epochs = [], [], [], []
    with (open_tenv_training(training) as release, open_tenv_training(bandwidth_training) as baseline):
        hvg = axis_evaluate.open_axes(family="hvg").axis(5000)
        if (baseline.manifest["inputs"]["axis_sha256"] != hvg.sha256
                or baseline.manifest["outputs"]["genes.csv"] != hvg.sha256):
            raise ValueError("fold bandwidths require the released HVG-DEG 5k training axis")
        data, base = AxisInputs(release), AxisInputs(baseline)
        for shape, fold in expected:
            path, config, receipt, training_record, trace = verified[shape, fold]
            if fold not in bandwidths:
                bandwidths[fold] = axis_evaluate._bandwidth(base, fold)
            paired, metrics = score_run(path, config, receipt, data, base, fold, smoke, bandwidths[fold]["bandwidth"])
            key = dict(shape=shape, fold=fold, seed=config["seed"], run_sha256=receipt["registration_sha256"])
            run_rows.append({
                **key, "source": "Step 06.1" if (shape, fold) in reused else "fitted",
                "manifest_sha256": axis_run.digest(path / "manifest.json"),
                "checkpoint_sha256": receipt["outputs"]["checkpoint.pt"],
                "trace_sha256": receipt["outputs"]["learning_curve.csv"],
                "objective": objective_label(shape), "prediction_head_parameters": HEAD_PARAMETERS[shape],
                "trainable_parameters": training_record["parameters"],
                "best_epoch": training_record["best_epoch"], "stopping_epoch": training_record["epochs"],
                "training_peak_rss_kib": training_record["peak_rss_kib"],
                "training_cuda_peak_allocated_bytes": training_record["cuda_peak_allocated_bytes"],
            })
            epochs.append(pd.DataFrame({**key, "objective": objective_label(shape), **dict(trace)}))
            conditions.append(pd.DataFrame({**key, **dict(paired)}))
            fold_metrics.append(pd.DataFrame({**key, **dict(metrics)}))
    fold_scores, summary = pair_scores(pd.concat(fold_metrics, ignore_index=True))
    candidate.mkdir()
    tables = {
        "runs.csv": pd.DataFrame(run_rows),
        "condition-scores.csv": pd.concat(conditions, ignore_index=True),
        "fold-scores.csv": fold_scores,
        "delta-metrics.csv": summary,
        "bandwidths.csv": pd.DataFrame(bandwidths.values()),
        "epoch-objectives.csv": pd.concat(epochs, ignore_index=True),
    }
    for name, table in tables.items():
        table.to_csv(candidate / name, index=False)
    axis_run.write_json(candidate / "manifest.json", {
        "status": "PASS", "schema": REPORT_SCHEMA, "axis": AXIS, "seed": axis_run.SEED,
        "smoke": smoke, "runs": len(verified), "new_runs": len(fresh), "reused_runs": len(reused),
        "fold_scores": len(fold_scores), "metrics": list(METRICS), "fold_reference": CONTROL,
        "selected_checkpoint_and_replay_verified": True,
        "aggregation": "equal eligible conditions per fold, equal eligible folds per shape, differences paired within fold",
        "scientific_replay_tolerances": {"rtol": 0., "atol": 0.},
        "outputs": axis_run.file_hashes(candidate, tables),
        "code_sha256": {path.relative_to(STAGE).as_posix(): axis_run.digest(path) for path in (
            Path(__file__), Path(axis_evaluate.__file__), axis_evaluate.STEP / "output/api/tsg_axis/metrics.py",
            STEP.with_name("06.1_round1") / "execution/evaluate_shape.py",
        )},
    })


def score_or_verify(output, scratch, runs, training, bandwidth_training, smoke, verify, reused=None):
    with axis_run.publication_lock(output) as target, TemporaryDirectory(prefix=".architecture-score-", dir=scratch) as directory:
        candidate = Path(directory) / "report"
        build_report(candidate, runs, training, bandwidth_training, smoke, reused)
        if verify:
            receipt = axis_run.read_json(candidate / "manifest.json")
            if axis_run.read_json(target / "manifest.json") != receipt:
                raise ValueError("reproduced Step 06.2 score receipt differs")
            axis_run.check_outputs(target, receipt)
        else:
            if target.exists():
                raise ValueError("Step 06.2 score candidate requires an unused output directory")
            candidate.rename(target)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("score", "verify"))
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--training", type=Path, default=DEFAULT_TRAINING)
    parser.add_argument("--bandwidth-training", type=Path, default=TENV_DEFAULT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    args = parser.parse_args()
    paths, _ = verify_matrix(args.matrix, fit_coordinates(True), MATRIX_SCHEMA)
    score_or_verify(args.output, args.scratch, paths.values(), args.training, args.bandwidth_training,
                    True, args.command == "verify")
    print(json.dumps({"status": "PASS", "report_sha256": axis_run.digest(args.output / "manifest.json")}, sort_keys=True))


if __name__ == "__main__":
    main()
