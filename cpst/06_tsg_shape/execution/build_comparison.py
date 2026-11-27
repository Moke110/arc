"""Build the four canonical Step 06 figure sets from accepted fit and report evidence."""

import argparse
from itertools import product
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

SHAPE = Path(__file__).resolve().parents[1]
STAGE = SHAPE.parent
sys.path[:0] = [str(SHAPE / "figure/script"), *(str(STAGE / path) for path in (
    "util/output/api", "04_tsg_axis/output/api", "03_data_train/output/api",
))]

from architecture_comparison import RENDERER, STUDIES
from cpst_util.comparison import CURVE_FIELDS, METRIC_FIELDS, METRIC_NAMES, write_paired_metrics
from cpst_util.comparison_publication import render_receipt
from cpst_util.derivation.release import check_outputs, read_json, sha256, write_json
from training_data import open_tenv_training, tenv_axis_output
from tsg_axis.comparison_evidence import read_fixed_jaccard
from tsg_axis.inputs import AxisInputs
from tsg_axis.metrics import summarize_conditions

KEY = ["shape", "fold", "seed"]
TRACE = ["epoch", "train_loss", "validation_loss"]
OLD_METRICS = set(METRIC_NAMES) - {"wJAC-2000"}
SCORE_FIELDS = ("value", "n_conditions", "n_eligible", "n_undefined")


def _table(path):
    return pd.read_csv(path, float_precision="round_trip")


def _checked_receipt(path, digest):
    if sha256(path) != digest:
        raise ValueError(f"accepted receipt digest differs for {path}")
    return read_json(path)


def _matrix_paths(output, accepted):
    matrix = output / "matrix"
    receipt = _checked_receipt(matrix / "manifest.json", accepted["matrix_manifest_sha256"])
    return {(row["registration_sha256"], row["run_manifest_sha256"]): matrix / row["path"]
            for row in receipt["fits"]}


def _fit_paths(evidence, names):
    prior = read_json(SHAPE / "06.1_round1/log/347-formal.json")
    paths = _matrix_paths(evidence / "06.1_round1/output/formal/prior-347", prior)
    by_step = {}
    for name in STUDIES:
        output = evidence / name / "output/formal"
        receipt = read_json(SHAPE / name / "log" / f"{STUDIES[name].ticket}-formal.json")
        current = _matrix_paths(output, receipt)
        if name == "06.3_round3":
            first = output / "prior-unnormalized"
            blocked = read_json(SHAPE / name / "log/356-blocked.json")
            partial = _checked_receipt(first / "partial-status.json", blocked["partial_ledger_sha256"])
            for row in partial["completed"]:
                arm, fold, seed = row["coordinate"]
                if arm in ("0_original", "1_output-delta", "2_output-delta-no-x2"):
                    key = (row["registration_sha256"], row["run_manifest_sha256"])
                    paths[key] = first / f"matrix/runs/{arm}/fold-{fold}/seed-{seed}"
        paths.update(current)
        by_step[name] = paths.copy()
        if set(names) <= by_step.keys():
            break
    return by_step


def accepted_fit(paths, record):
    """Resolve a fit matching both accepted run-ledger hashes and its formal coordinate."""
    identity = record.shape, record.fold, record.seed
    key = record.run_sha256, record.manifest_sha256
    run = paths.get(key)
    if run is None:
        raise ValueError(f"missing accepted fit source for {identity}")
    registration, manifest_path = (run / name for name in ("registration.json", "manifest.json"))
    if (not registration.is_file() or not manifest_path.is_file()
            or (sha256(registration), sha256(manifest_path)) != key):
        raise ValueError(f"accepted registration or replay differs for {identity}: {run}")
    config, manifest = map(read_json, (registration, manifest_path))
    if (manifest["registration_sha256"] != record.run_sha256
            or (manifest["axis"], manifest["fold"], manifest["seed"], manifest["smoke"])
            != ("denoise-5k", record.fold, record.seed, False)
            or any(manifest[name] != "exact" for name in ("test_replay", "validation_replay",
                                                         "validation_objective_replay"))):
        raise ValueError(f"accepted run coordinate or replay differs for {identity}")
    return run, config, manifest


def _read_trace(run, record, manifest, epochs):
    path = run / "learning_curve.csv"
    trace_sha = manifest["outputs"]["learning_curve.csv"]
    identity = record.shape, record.fold, record.seed
    if sha256(path) != trace_sha or getattr(record, "trace_sha256", trace_sha) != trace_sha:
        raise ValueError(f"accepted trace digest differs for {identity}")
    trace = _table(path).loc[:, TRACE]
    if (trace.epoch.tolist() != list(range(1, record.stopping_epoch + 1))
            or (epochs is not None and not trace.equals(epochs.loc[:, TRACE].reset_index(drop=True)))):
        raise ValueError(f"accepted epoch history differs for {identity}")
    return trace


def _write_plot_data(report, fits, study, data, output):
    runs = _table(report / "runs.csv")
    expected = set(product(study.arms, range(3), study.seeds))
    if len(runs) != len(expected) or set(runs[KEY].itertuples(index=False, name=None)) != expected:
        raise ValueError(f"accepted fit coordinates differ for {report}")

    tables = [_table(report / name) for name in ("condition-scores.csv", "fold-scores.csv")]
    if (report / "epoch-objectives.csv").is_file():
        tables.append(_table(report / "epoch-objectives.csv"))
    ledger = set(runs[[*KEY, "run_sha256"]].itertuples(index=False, name=None))
    if any(set(table[[*KEY, "run_sha256"]].itertuples(index=False, name=None)) != ledger for table in tables):
        raise ValueError(f"accepted score or trace identities differ for {report}")
    groups = [table.groupby(KEY, sort=False) for table in tables]
    curves, metrics = [], []
    for record in runs.itertuples(index=False):
        identity = record.shape, record.fold, record.seed
        run, config, manifest = accepted_fit(fits, record)
        sealed, scores, *epochs = (group.get_group(identity) for group in groups)
        trace = _read_trace(run, record, manifest, next(iter(epochs), None))
        conditions, top200, weighted = read_fixed_jaccard(run, config, manifest, data)
        if (conditions != tuple(sealed.condition_id)
                or not np.array_equal(top200, sealed.JAC_top200, equal_nan=True)):
            raise ValueError(f"held-out JAC-top200 differs from accepted conditions for {identity}")
        scores = scores.set_index("metric", verify_integrity=True)
        if not OLD_METRICS <= set(scores.index):
            raise ValueError(f"accepted fold metrics differ for {identity}")
        jacc = scores.loc["JAC_top200"]
        observed = summarize_conditions(top200)
        if (not np.isclose(jacc.value, observed["value"], rtol=0, atol=1e-15, equal_nan=True)
                or any(observed[key] != jacc[key] for key in SCORE_FIELDS[1:])):
            raise ValueError(f"held-out JAC-top200 fold score differs for {identity}")
        weighted_mean = summarize_conditions(weighted)["value"]
        curves.extend((*identity, *row) for row in trace.itertuples(index=False, name=None))
        metrics.extend((*identity, metric, weighted_mean if metric == "wJAC-2000" else scores.at[metric, "value"])
                       for metric in METRIC_NAMES)

    output.mkdir()
    pd.DataFrame(curves, columns=CURVE_FIELDS).to_csv(output / "learning-curve.csv", index=False)
    pd.DataFrame(metrics, columns=METRIC_FIELDS).to_csv(output / "metrics.csv", index=False)
    write_paired_metrics(output, study.reference)
    return len(runs)


def build_study(source, output, fits, data):
    """Verify the accepted report and publish candidate plotted CSVs, PNGs and replay receipt."""
    study = STUDIES[source.name]
    step = SHAPE / source.name
    receipt_path = step / "log" / f"{study.ticket}-formal.json"
    accepted = read_json(receipt_path)
    report = source / "output/formal/report"
    manifest = _checked_receipt(report / "manifest.json", accepted["report_manifest_sha256"])
    if (accepted["status"] != "PASS" or manifest["status"] != "PASS"
            or manifest.get("reference", manifest.get("fold_reference")) != study.reference
            or check_outputs(report, manifest) != accepted["report_outputs_sha256"]):
        raise ValueError(f"accepted report or reference differs for {source.name}")

    output.mkdir()
    figure = output / "figure"
    figure.mkdir()
    plotted = figure / "data"
    count = _write_plot_data(report, fits, study, data, plotted)
    renderer = step / RENDERER
    dependencies = (SHAPE / "figure/script/architecture_comparison.py",)
    write_json(output / "log/figure-render.json", {
        **render_receipt(STAGE, step.name, output, renderer, dependencies=dependencies),
        "evidence": {
            "acceptance_receipt_sha256": sha256(receipt_path),
            "matrix_manifest_sha256": accepted["matrix_manifest_sha256"],
            "report_manifest_sha256": accepted["report_manifest_sha256"],
            "build_sha256": sha256(Path(__file__)),
            "matched_controls_sha256": data.release.manifest_sha256,
            "fixed_jaccard_code_sha256": sha256(STAGE / "04_tsg_axis/output/api/tsg_axis/comparison_evidence.py"),
            "shared_metrics_sha256": sha256(STAGE / "util/output/api/cpst_util/metrics.py"),
            "sealed_jac_top200_parity": True,
            "fits": count,
            "reference": study.reference,
        },
    })
    print(json.dumps({"status": "PASS", "step": step.name, "fits": count}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-root", type=Path, default=STAGE.parent,
                        help="materialized checkout containing the accepted Step 03 and Step 06 releases")
    parser.add_argument("--step", action="append", choices=STUDIES, help="rebuild only the named Step 06 comparison")
    parser.add_argument("--output", type=Path, required=True, help="new candidate directory under tmp/")
    args = parser.parse_args()
    candidate = args.output.resolve()
    evidence = args.evidence_root.resolve() / STAGE.name
    shape = evidence / SHAPE.name
    names = args.step or STUDIES
    paths = _fit_paths(shape, names)
    candidate.mkdir()
    with open_tenv_training(tenv_axis_output("denoise-5k", evidence / "03_data_train")) as release:
        data = AxisInputs(release)
        for name in names:
            build_study(shape / name, candidate / name, paths[name], data)


if __name__ == "__main__":
    main()
