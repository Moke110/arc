"""Score and independently reproduce the Step 05 report from sealed fits."""

import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from run_module import STEP, axis_run
from measured_roster import SHAPES
from module_roster import coordinates
from measured_scores import CONTRASTS, CONTRAST_METRICS, PRIMARY_METRICS, paired_report

STAGE = STEP.parent

from evaluate_loss import axis_evaluate, score_run
from matrix_measured import verify_matrix
from run_module import AXIS, DEFAULT_TRAINING
from training_data import TENV_DEFAULT, open_tenv_training
from tsg_axis.inputs import AxisInputs
from cpst_util.comparison import METRIC_NAMES
from verify_module_fit import read_run

REPORT_SCHEMA = "cpst-tenv-transigen-module-report-v1"


def verify_configurations(verified):
    common, fold_configs = None, {}
    variable_fields = {"arguments", "nonce", "schema", "protocol", "shape", "seed", "objective",
                       "conditions", "partitions", "model", "experiment"}
    model_fields = {"constructor_order", "initialization", "prediction_head_parameters",
                    "validation_objective", "transform", "measured_encoding", "delta_readout"}
    for (shape, fold, seed), (_, config, _, _, _) in verified.items():
        shared = {name: value for name, value in config.items() if name not in variable_fields}
        shared["model"] = {name: value for name, value in config["model"].items() if name not in model_fields}
        membership = {name: config[name] for name in ("partitions", "conditions")}
        membership["transform"] = config["model"]["transform"]
        if common is None:
            common = shared
        if (seed != config["seed"] or common != shared or fold_configs.setdefault(fold, membership) != membership):
            raise ValueError(f"mismatched protocol or population in {shape}, fold {fold}, seed {seed}")


def build_report(candidate, paths, training, bandwidth_training, smoke):
    expected = coordinates(smoke)
    verified = dict(read_run(path, smoke) for path in paths.values())
    if len(verified) != len(paths) or set(verified) != set(expected):
        raise ValueError("Expected one fresh fit for each registered coordinate")
    verify_configurations(verified)

    bandwidths = {}
    run_rows, conditions, fold_metrics, epochs = [], [], [], []
    with open_tenv_training(training) as release, open_tenv_training(bandwidth_training) as baseline:
        hvg = axis_evaluate.open_axes(family="hvg").axis(5000)
        if (baseline.manifest["inputs"]["axis_sha256"] != hvg.sha256
                or baseline.manifest["outputs"]["genes.csv"] != hvg.sha256):
            raise ValueError("fold bandwidths require the released HVG-DEG 5k training axis")
        data, base = AxisInputs(release), AxisInputs(baseline)
        protocol = axis_run.shared_protocol()
        data.bind_protocol(protocol)
        base.bind_protocol(protocol)
        for coordinate in expected:
            shape, fold, seed = coordinate
            path, config, receipt, training_record, trace = verified[coordinate]
            if fold not in bandwidths:
                bandwidths[fold] = axis_evaluate._bandwidth(base, fold)
            by_condition, metrics = score_run(
                path, config, receipt, data, base, fold, False, bandwidths[fold]["bandwidth"])
            key = dict(shape=shape, fold=fold, seed=seed, run_sha256=receipt["registration_sha256"])
            objective = "full registered module objective"
            run_rows.append({
                **key, "source": "Step 05", "protocol": config["protocol"],
                "manifest_sha256": axis_run.digest(path / "manifest.json"),
                "checkpoint_sha256": receipt["outputs"]["checkpoint.pt"],
                "trace_sha256": receipt["outputs"]["learning_curve.csv"],
                "prediction_sha256": receipt["outputs"]["prediction.npy"],
                "objective": objective, "prediction_head_parameters": config["model"].get("prediction_head_parameters", 106600),
                "trainable_parameters": training_record["parameters"],
                "best_epoch": training_record["best_epoch"], "stopping_epoch": training_record["epochs"],
                "training_peak_rss_kib": training_record["peak_rss_kib"],
                "training_cuda_peak_allocated_bytes": training_record["cuda_peak_allocated_bytes"],
            })
            epochs.append(pd.DataFrame({**key, "objective": objective, **dict(trace)}))
            by_condition = by_condition.drop(columns=["JAC_pos100", "JAC_neg100"])
            conditions.append(pd.DataFrame({**key, **dict(by_condition)}))
            metrics = metrics.loc[metrics.metric.isin(METRIC_NAMES)]
            fold_metrics.append(pd.DataFrame({**key, **dict(metrics)}))
    paired, summary, contrasts, contrast_by_fold, contrast_summary = paired_report(
        pd.concat(fold_metrics, ignore_index=True), expected, METRIC_NAMES)
    bandwidth_table = pd.DataFrame(bandwidths.values())
    candidate.mkdir()
    tables = {
        "runs.csv": pd.DataFrame(run_rows),
        "condition-scores.csv": pd.concat(conditions, ignore_index=True),
        "fold-scores.csv": paired,
        "delta-metrics.csv": summary,
        "primary-contrasts.csv": contrasts,
        "primary-contrast-by-fold.csv": contrast_by_fold,
        "primary-contrast-summary.csv": contrast_summary,
        "bandwidths.csv": bandwidth_table,
        "epoch-objectives.csv": pd.concat(epochs, ignore_index=True),
    }
    for name, table in tables.items():
        table.to_csv(candidate / name, index=False)
    axis_run.write_json(candidate / "manifest.json", {
        "status": "PASS", "schema": REPORT_SCHEMA, "axis": AXIS, "smoke": smoke,
        "training_seeds": list(dict.fromkeys(seed for _, _, seed in expected)),
        "held_out_folds": list(dict.fromkeys(fold for _, fold, _ in expected)),
        "runs": len(verified), "new_runs": len(verified), "reused_runs": 0,
        "fold_scores": len(paired), "metrics": list(METRIC_NAMES), "reference": SHAPES[0],
        "primary_metrics": list(PRIMARY_METRICS), "contrast_metrics": list(CONTRAST_METRICS),
        "contrasts": list(CONTRASTS), "primary_contrast_scores": len(contrasts),
        "primary_contrast_by_fold": len(contrast_by_fold), "primary_contrast_aggregates": len(contrast_summary),
        "selected_checkpoint_and_replay_verified": True,
        "aggregation": "equal eligible conditions per fold; equal defined fold and seed coordinates per arm",
        "scientific_replay_tolerances": {"rtol": 0., "atol": 0.},
        "outputs": axis_run.file_hashes(candidate, tables),
        "code_sha256": {path.relative_to(STAGE).as_posix(): axis_run.digest(path) for path in (
            Path(__file__), Path(__file__).with_name("measured_scores.py"),
            Path(axis_evaluate.__file__), axis_evaluate.STEP / "output/api/tsg_axis/metrics.py",
        )},
    })


def score_or_verify(output, scratch, paths, training, bandwidth_training, smoke, verify):
    with (axis_run.publication_lock(output) as target,
          TemporaryDirectory(prefix=".measured-score-", dir=scratch) as directory):
        candidate = Path(directory) / "report"
        build_report(candidate, paths, training, bandwidth_training, smoke)
        receipt = axis_run.read_json(candidate / "manifest.json")
        if verify:
            if axis_run.read_json(target / "manifest.json") != receipt:
                raise ValueError("reproduced Step 05 score receipt differs")
            axis_run.check_outputs(target, receipt)
        else:
            if target.exists():
                raise ValueError("Step 05 score candidate requires an unused output directory")
            candidate.rename(target)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("score", "verify"))
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--training", type=Path, default=DEFAULT_TRAINING)
    parser.add_argument("--bandwidth-training", type=Path, default=TENV_DEFAULT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    paths, _ = verify_matrix(args.matrix, args.smoke)
    score_or_verify(args.output, args.scratch, paths, args.training, args.bandwidth_training,
                    args.smoke, args.command == "verify")
    print(json.dumps({"status": "PASS", "report_sha256": axis_run.digest(args.output / "manifest.json")},
                     sort_keys=True))


if __name__ == "__main__":
    main()
