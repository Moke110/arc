"""Build the three Step 05 figures from the independently verified report."""

import argparse
from pathlib import Path

import pandas as pd

from run_module import STAGE, STEP, axis_run
from measured_roster import CONTROL
from module_roster import coordinates
from cpst_util.comparison import CURVE_FIELDS, METRIC_FIELDS, METRIC_NAMES, write_paired_metrics
from cpst_util.comparison_publication import render_receipt


def build(report, output, smoke):
    manifest = report / "manifest.json"
    receipt = axis_run.read_json(manifest)
    roster = coordinates(smoke)
    if (receipt["status"] != "PASS" or receipt["smoke"] != smoke
            or receipt["runs"] != len(roster) or receipt["metrics"] != list(METRIC_NAMES)):
        raise ValueError("Report does not match the registered module study")
    axis_run.check_outputs(report, receipt)

    curves = pd.read_csv(report / "epoch-objectives.csv", float_precision="round_trip").rename(
        columns={"shape": "arm", "train_loss": "train", "validation_loss": "validation"})
    scores = pd.read_csv(report / "fold-scores.csv", float_precision="round_trip").rename(
        columns={"shape": "arm"})
    expected = {(arm, fold, seed, metric) for arm, fold, seed in roster for metric in METRIC_NAMES}
    keys = list(scores[["arm", "fold", "seed", "metric"]].itertuples(index=False, name=None))
    if set(keys) != expected or len(keys) != len(expected):
        raise ValueError("Plot scores must contain all registered fold-and-seed coordinates")
    data = output / "figure/data"
    data.mkdir(parents=True)
    curves[list(CURVE_FIELDS)].to_csv(data / "learning-curve.csv", index=False)
    scores[list(METRIC_FIELDS)].to_csv(data / "metrics.csv", index=False)
    write_paired_metrics(data, CONTROL)
    renderer = STEP / "figure/script/render.py"
    result = render_receipt(STAGE, STEP.name, output, renderer)
    result["evidence"] = {
        "report_manifest_sha256": axis_run.digest(manifest),
        "build_sha256": axis_run.digest(Path(__file__)),
        "fits": len(roster),
    }
    (output / "log").mkdir()
    axis_run.write_json(output / "log/figure-render.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    build(args.report.resolve(), args.output.resolve(), args.smoke)


if __name__ == "__main__":
    main()
