"""Score and independently reproduce the Step 06.3 report from sealed fits."""

import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from composition_model import STEP, axis_run
from composition_roster import CONTROL, HEAD_PARAMETERS, coordinates, objective_label
from composition_scores import CONTRASTS, paired_report

STAGE = STEP.parents[1]

from evaluate_loss import axis_evaluate, score_run
from matrix_composition import verify_matrix
from prior_composition import ROUND_ONE_MATRIX, ROUND_TWO_MATRIX, load_reused
from run_composition import AXIS, DEFAULT_TRAINING
from training_data import TENV_DEFAULT, open_tenv_training
from tsg_axis.inputs import AxisInputs
from tsg_axis.metrics import METRICS
from verify_composition import read_run

REPORT_SCHEMA = "stage2-tenv-transigen-architecture-composition-report-v1"
SOURCE_SCORE_FIELDS = [
    "shape", "fold", "seed", "run_sha256", "metric", "value", "n_conditions", "n_eligible", "n_undefined",
]


def verify_configurations(verified):
    common, fold_configs = None, {}
    variable_fields = {"arguments", "nonce", "schema", "protocol", "shape", "seed", "objective",
                       "conditions", "partitions", "model", "code_sha256"}
    model_fields = {"constructor_order", "initialization", "prediction_head_parameters",
                    "validation_objective", "transform"}
    for (shape, fold, seed), (_, config, _, _, _) in verified.items():
        shared = {name: value for name, value in config.items() if name not in variable_fields}
        shared["model"] = {name: value for name, value in config["model"].items() if name not in model_fields}
        shared["code_sha256"] = {name: value for name, value in config["code_sha256"].items()
                                 if not name.startswith("stage2/06_tsg_shape/")}
        membership = {name: config[name] for name in ("partitions", "conditions")}
        membership["transform"] = config["model"]["transform"]
        if common is None:
            common = shared
        if (seed != config["seed"] or common != shared or fold_configs.setdefault(fold, membership) != membership):
            raise ValueError(f"mismatched protocol or population in {shape}, fold {fold}, seed {seed}")


def verify_source_scores(scores, bandwidths, reused, source_scores, source_bandwidths):
    index = pd.MultiIndex.from_tuples(reused, names=["shape", "fold", "seed"])
    reconstructed = scores.loc[pd.MultiIndex.from_frame(scores[["shape", "fold", "seed"]]).isin(index),
                               SOURCE_SCORE_FIELDS]
    score_index = SOURCE_SCORE_FIELDS[:3] + ["metric"]
    found = reconstructed.set_index(score_index).sort_index()
    accepted = source_scores[SOURCE_SCORE_FIELDS].set_index(score_index).sort_index()
    if not found.equals(accepted):
        raise ValueError("reproduced seed-42 scores or eligibility differ from accepted Round 2")
    if not bandwidths.sort_values("fold").reset_index(drop=True).equals(
            source_bandwidths[bandwidths.columns].sort_values("fold").reset_index(drop=True)):
        raise ValueError("reproduced fold bandwidths differ from accepted Round 2")


def build_report(candidate, paths, training, bandwidth_training, smoke, prior):
    reused, provenance, source_scores, source_bandwidths = prior
    expected = tuple(coordinates(smoke))
    fresh = dict(read_run(path, smoke) for path in paths.values())
    if len(fresh) != len(paths) or set(fresh) != set(paths):
        raise ValueError("duplicated or mismatched fresh fit coordinate")
    verified = reused | fresh
    if set(reused) & set(fresh) or set(verified) != set(expected):
        raise ValueError("reused and fresh fits must supply exactly one fit per registered coordinate")
    verify_configurations(verified)

    bandwidths = {}
    run_rows, conditions, fold_metrics, epochs = [], [], [], []
    with open_tenv_training(training) as release, open_tenv_training(bandwidth_training) as baseline:
        hvg = axis_evaluate.open_axes(family="hvg").axis(5000)
        if (baseline.manifest["inputs"]["axis_sha256"] != hvg.sha256
                or baseline.manifest["outputs"]["genes.csv"] != hvg.sha256):
            raise ValueError("fold bandwidths require the released HVG-DEG 5k training axis")
        data, base = AxisInputs(release), AxisInputs(baseline)
        for coordinate in expected:
            shape, fold, seed = coordinate
            path, config, receipt, training_record, trace = verified[coordinate]
            if fold not in bandwidths:
                bandwidths[fold] = axis_evaluate._bandwidth(base, fold)
            by_condition, metrics = score_run(
                path, config, receipt, data, base, fold, smoke, bandwidths[fold]["bandwidth"])
            key = dict(shape=shape, fold=fold, seed=seed, run_sha256=receipt["registration_sha256"])
            objective = objective_label(shape)
            source = ("Step 06.1" if shape == CONTROL else "Step 06.2") if coordinate in reused else "Step 06.3"
            run_rows.append({
                **key, "source": source, "source_shape": config["shape"], "protocol": config["protocol"],
                "manifest_sha256": axis_run.digest(path / "manifest.json"),
                "checkpoint_sha256": receipt["outputs"]["checkpoint.pt"],
                "trace_sha256": receipt["outputs"]["learning_curve.csv"],
                "objective": objective, "prediction_head_parameters": HEAD_PARAMETERS[shape],
                "trainable_parameters": training_record["parameters"],
                "best_epoch": training_record["best_epoch"], "stopping_epoch": training_record["epochs"],
                "training_peak_rss_kib": training_record["peak_rss_kib"],
                "training_cuda_peak_allocated_bytes": training_record["cuda_peak_allocated_bytes"],
            })
            epochs.append(pd.DataFrame({**key, "objective": objective, **dict(trace)}))
            conditions.append(pd.DataFrame({**key, **dict(by_condition)}))
            fold_metrics.append(pd.DataFrame({**key, **dict(metrics)}))
    paired, summary, contrasts, contrast_summary = paired_report(
        pd.concat(fold_metrics, ignore_index=True), expected, METRICS)
    bandwidth_table = pd.DataFrame(bandwidths.values())
    if reused:
        verify_source_scores(paired, bandwidth_table, reused, source_scores, source_bandwidths)
    candidate.mkdir()
    tables = {
        "runs.csv": pd.DataFrame(run_rows),
        "condition-scores.csv": pd.concat(conditions, ignore_index=True),
        "fold-scores.csv": paired,
        "delta-metrics.csv": summary,
        "primary-contrasts.csv": contrasts,
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
        "runs": len(verified), "new_runs": len(fresh), "reused_runs": len(reused),
        "fold_scores": len(paired), "metrics": list(METRICS), "reference": CONTROL,
        "contrasts": list(CONTRASTS), "primary_contrast_scores": len(contrasts),
        "exact_reused_scores": len(reused) * len(METRICS), "prior_release": provenance,
        "selected_checkpoint_and_replay_verified": True,
        "aggregation": "equal eligible conditions per fold; equal defined fold and seed coordinates per arm",
        "scientific_replay_tolerances": {"rtol": 0., "atol": 0.},
        "outputs": axis_run.file_hashes(candidate, tables),
        "code_sha256": {path.relative_to(STAGE).as_posix(): axis_run.digest(path) for path in (
            Path(__file__), Path(__file__).with_name("composition_scores.py"),
            Path(axis_evaluate.__file__), axis_evaluate.STEP / "output/api/tsg_axis/metrics.py",
        )},
    })


def score_or_verify(output, scratch, paths, training, bandwidth_training, smoke, verify, prior):
    with (axis_run.publication_lock(output) as target,
          TemporaryDirectory(prefix=".composition-score-", dir=scratch) as directory):
        candidate = Path(directory) / "report"
        build_report(candidate, paths, training, bandwidth_training, smoke, prior)
        receipt = axis_run.read_json(candidate / "manifest.json")
        if verify:
            if axis_run.read_json(target / "manifest.json") != receipt:
                raise ValueError("reproduced Step 06.3 score receipt differs")
            axis_run.check_outputs(target, receipt)
        else:
            if target.exists():
                raise ValueError("Step 06.3 score candidate requires an unused output directory")
            candidate.rename(target)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("score", "verify"))
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--training", type=Path, default=DEFAULT_TRAINING)
    parser.add_argument("--bandwidth-training", type=Path, default=TENV_DEFAULT)
    parser.add_argument("--round-one-matrix", type=Path, default=ROUND_ONE_MATRIX)
    parser.add_argument("--round-two-matrix", type=Path, default=ROUND_TWO_MATRIX)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    paths, _ = verify_matrix(args.matrix, args.smoke)
    prior = load_reused(args.round_two_matrix, args.round_one_matrix)
    if args.smoke:
        prior = {}, prior[1], None, None
    score_or_verify(args.output, args.scratch, paths, args.training, args.bandwidth_training,
                    args.smoke, args.command == "verify", prior)
    print(json.dumps({"status": "PASS", "report_sha256": axis_run.digest(args.output / "manifest.json")},
                     sort_keys=True))


if __name__ == "__main__":
    main()
