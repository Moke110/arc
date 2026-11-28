"""Rescore sealed Step 08 predictions and summarize paired metrics against KPGT."""

import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np
import pandas as pd

from chem_matrix import verify_matrix
from chem_runs import REFERENCE, add_inputs, axis_run, bind, coordinates, open_inputs, read_run, verify_configuration
from evaluate_loss import axis_evaluate, score_run
from training_data import TENV_DEFAULT, open_tenv_training
from tsg_axis.inputs import AxisInputs
from tsg_axis.metrics import METRICS

PRIMARY_METRICS = ("Δ PCC", "Δ RMSE")
SCORE_KEYS = ["arm", "fold", "seed", "metric"]
SCHEMA = "cpst-alpha-chemical-report-v1"


def paired_scores(scores, roster):
    reference_coordinates = {(fold, seed) for arm, fold, seed in roster if arm == REFERENCE}
    if reference_coordinates != {(fold, seed) for _, fold, seed in roster}:
        raise ValueError("every scored fold and seed requires a KPGT reference")
    expected = pd.MultiIndex.from_tuples(
        [(*coordinate, metric) for coordinate, metric in product(roster, METRICS)], names=SCORE_KEYS)
    indexed = scores.set_index(SCORE_KEYS)
    if not indexed.index.is_unique or len(indexed) != len(expected) or not expected.difference(indexed.index).empty:
        raise ValueError("scores must supply each arm, fold, seed and metric exactly once")
    counts = indexed[["n_conditions", "n_eligible", "n_undefined"]]
    if (counts.isna().any().any() or (counts < 0).any().any()
            or not counts.eq(counts.round()).all().all()
            or not counts.n_conditions.eq(counts.n_eligible + counts.n_undefined).all()
            or not indexed.value.notna().eq(counts.n_eligible.gt(0)).all()
            or not np.isfinite(indexed.value.dropna()).all()):
        raise ValueError("score values or eligibility do not match the frozen roster")
    paired = indexed.reindex(expected).reset_index()
    groups = ["fold", "seed", "metric"]
    reference = paired.loc[paired.arm.eq(REFERENCE)].set_index(groups)
    paired["kpgt_value"] = reference.value.reindex(pd.MultiIndex.from_frame(paired[groups])).to_numpy()
    paired["score_difference"] = paired.value - paired.kpgt_value
    aggregates = {
        "value": ("value", "mean"), "sample_sd": ("value", "std"), "n_defined": ("value", "count"),
        "kpgt_value": ("kpgt_value", "mean"), "score_difference": ("score_difference", "mean"),
        "score_difference_sd": ("score_difference", "std"), "n_paired": ("score_difference", "count"),
        "n_coordinates": ("value", "size"), "n_eligible_conditions": ("n_eligible", "sum"),
        "n_undefined_conditions": ("n_undefined", "sum"),
    }
    by_fold = paired.groupby(["arm", "metric", "fold"], sort=False, as_index=False).agg(**aggregates)
    summary = paired.groupby(["arm", "metric"], sort=False, as_index=False).agg(**aggregates)
    return paired, by_fold, summary


def verify_initializations(fresh):
    """Check that all fresh arms retain the same nonprojection state and RNG per coordinate."""
    shared = {}
    for (arm, fold, seed), (_, config, _, training, _) in fresh.items():
        observed = training["initialization"]
        rep = config["representation"]
        if observed["projection_shape"] != [rep["width"], rep["projection_width"]]:
            raise ValueError(f"{arm} did not use its native chemical input width")
        paired = tuple(observed[key] for key in (
            "downstream_sha256", "torch_cpu_rng_sha256", "torch_cuda_rng_sha256"))
        if shared.setdefault((fold, seed), paired) != paired:
            raise ValueError(f"{arm} initialization or RNG was not paired at fold {fold}, seed {seed}")
    return shared


def verified_runs(matrix, smoke):
    paths = verify_matrix(matrix, smoke)
    return dict(read_run(path, smoke) for path in paths.values())


def report_tables(fresh, reused, training, representations, source, bandwidth_training, smoke):
    """Independently rescore each fit, including supplied verified reference checkpoints."""
    roster = coordinates(smoke)
    if (set(fresh) & set(reused) or set(fresh) | set(reused) != set(roster)
            or any(key[0] != REFERENCE for key in reused)):
        raise ValueError("fresh and reused fits do not provide exactly the registered coordinates")
    verify_initializations(fresh)
    bandwidths, runs, fold_metrics, conditions, traces = {}, [], [], [], []
    with open_inputs(training, representations) as (data, vectors), open_tenv_training(bandwidth_training) as release:
        hvg = axis_evaluate.open_axes(family="hvg").axis(5000)
        if release.manifest["inputs"]["axis_sha256"] != hvg.sha256 or release.manifest["outputs"]["genes.csv"] != hvg.sha256:
            raise ValueError("fold bandwidth requires the released HVG-DEG 5k training axis")
        base = AxisInputs(release)
        for key in roster:
            arm, fold, seed = key
            reused_fit = key in reused
            path, config, replay, trained, trace = (reused if reused_fit else fresh)[key]
            if not reused_fit:
                identity = {name: value for name, value in config.items() if name not in ("arguments", "nonce")}
                verify_configuration(identity, bind(key, data, vectors, source, smoke), replay=True)
            if fold not in bandwidths:
                bandwidths[fold] = axis_evaluate._bandwidth(base, fold)
            condition, metrics = score_run(path, config, replay, data, base, fold, smoke,
                                           bandwidths[fold]["bandwidth"])
            keys = dict(arm=arm, fold=fold, seed=seed, registration_sha256=replay["registration_sha256"])
            runs.append({**keys, "fit_origin": "Step 06.4" if reused_fit else "Step 08",
                         "native_width": vectors.width(arm), "manifest_sha256": axis_run.digest(path / "manifest.json"),
                         "checkpoint_sha256": replay["outputs"]["checkpoint.pt"],
                         "prediction_sha256": replay["outputs"]["prediction.npy"],
                         "trace_sha256": replay["outputs"]["learning_curve.csv"],
                         "best_epoch": trained["best_epoch"], "stopping_epoch": trained["epochs"],
                         "parameters": trained["parameters"]})
            fold_metrics.append(metrics.assign(**keys))
            conditions.append(condition.assign(**keys))
            traces.append(trace.assign(**keys))
    paired, by_fold, summary = paired_scores(pd.concat(fold_metrics, ignore_index=True), roster)
    return {"fold-scores.csv": paired, "paired-by-fold.csv": by_fold, "paired-summary.csv": summary,
            "condition-scores.csv": pd.concat(conditions, ignore_index=True), "runs.csv": pd.DataFrame(runs),
            "bandwidths.csv": pd.DataFrame(bandwidths.values()), "epoch-objectives.csv": pd.concat(traces, ignore_index=True)}


def write_report(output, tables, smoke):
    output.mkdir()
    for name, table in tables.items():
        table.to_csv(output / name, index=False)
    receipt = {"schema": SCHEMA, "smoke": smoke, "scores": len(tables["fold-scores.csv"]),
               "csv_sha256": axis_run.file_hashes(output, tables)}
    axis_run.write_json(output / "manifest.json", receipt)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_inputs(parser)
    parser.add_argument("--bandwidth-training", type=Path, default=TENV_DEFAULT)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    fresh = verified_runs(args.matrix, args.smoke)
    inputs = args.training, args.representations, args.source, args.bandwidth_training
    tables = report_tables(fresh, {}, *inputs, args.smoke)
    print(json.dumps(write_report(args.output, tables, args.smoke), sort_keys=True))


if __name__ == "__main__":
    main()
