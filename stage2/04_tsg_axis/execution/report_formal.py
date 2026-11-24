"""Build and exactly replay the formal score report and its CSV-only figures."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
import time

from matrix import SCHEMA, STAGE, STEP, file_sha256, read, run_paths, verify_run

FIGURES = ("delta-metrics.png", "disease-f1.png")
FIGURE_CSVS = ("delta-metrics.csv", "disease-f1.csv")


def matrix_inputs(matrix, audit):
    """Verify the full matrix receipt before providing runs to the evaluator."""
    receipt = read(matrix / "manifest.json")
    paths = run_paths(matrix)
    records = {(item["axis"], item["fold"]): item for item in receipt["fits"]}
    if (receipt["schema"] != SCHEMA or set(records) != set(paths)
            or len(receipt["fits"]) != len(paths)
            or receipt["preflight_sha256"] != file_sha256(matrix / "preflight.json")
            or set(paths.values()) != {path for path in (matrix / "runs").glob("*/*") if path.is_dir()}):
        raise ValueError("formal matrix manifest must bind the exact 18-fit inventory")
    for (axis, fold), path in paths.items():
        registration = file_sha256(path / "registration.json")
        if records[axis, fold] != verify_run(path, axis, fold, registration):
            raise ValueError(f"formal matrix fit identity differs: {axis}, fold {fold}")
    logs = {path.relative_to(audit).as_posix(): file_sha256(path)
            for path in sorted((audit / "logs").rglob("*")) if path.is_file()}
    if not logs or logs != receipt["resource_logs"]:
        raise ValueError("formal matrix process logs differ from their receipt")
    return paths


def command(script, args, logs, name):
    attempt = logs / f"{name}-{time.time_ns()}"
    with attempt.with_suffix(".log").open("xb") as output:
        subprocess.run(["/usr/bin/time", "-v", "-o", str(attempt.with_suffix(".resources")),
                        sys.executable, str(script), *args],
                       stdout=output, stderr=subprocess.STDOUT, check=True)


def report(args):
    paths = matrix_inputs(args.matrix, args.matrix_audit or args.matrix)
    args.logs.mkdir(parents=True, exist_ok=True)
    evaluate = STEP / "execution/evaluate.py"
    render = STEP / "figure/script/render.py"
    options = ["--step03", str(args.step03), "--raw", str(args.raw), "--output", str(args.report)]
    for path in paths.values():
        options.extend(("--run", str(path)))
    if args.command == "build":
        if args.figures.exists():
            raise ValueError("report building requires an unused figures directory")
        command(evaluate, ["score", *options], args.logs, "score")
    command(evaluate, ["verify", *options], args.logs, "verify")
    args.figures.parent.mkdir(parents=True, exist_ok=True)
    staging = args.figures.parent if args.command == "build" else args.logs
    with TemporaryDirectory(prefix="csv-figures-", dir=staging) as scratch:
        scratch = Path(scratch)
        figures = scratch / "figures"
        csvs = scratch / "csv"
        csvs.mkdir()
        for name in FIGURE_CSVS:
            shutil.copyfile(args.report / "figure/data" / name, csvs / name)
        command(render, ["--data", str(csvs), "--output", str(figures)], args.logs, "csv-render")
        hashes = {name: file_sha256(figures / name) for name in FIGURES}
        if args.command == "build":
            replay = scratch / "replay"
            command(render, ["--data", str(csvs), "--output", str(replay)], args.logs, "csv-replay")
            if any(file_sha256(replay / name) != digest for name, digest in hashes.items()):
                raise ValueError("independent CSV-only PNG renders differ")
            csvs.rename(figures / "data")
            figures.rename(args.figures)
        else:
            if any(file_sha256(args.figures / name) != digest for name, digest in hashes.items()):
                raise ValueError("CSV-only figure replay differs from the accepted PNG")
            for name in FIGURE_CSVS:
                if file_sha256(args.figures / "data" / name) != file_sha256(csvs / name):
                    raise ValueError(f"published plot CSV differs from the accepted report: {name}")
    scoring_code = (STEP / "output/api/tsg_axis/metrics.py", STEP / "output/api/tsg_axis/disease_reference.py")
    return {"result": "PASS", "fits": len(paths),
            "matrix_manifest_sha256": file_sha256(args.matrix / "manifest.json"),
            "report_manifest_sha256": file_sha256(args.report / "manifest.json"),
            "csv_only_figure_sha256": hashes,
            "code_sha256": {path.relative_to(STEP).as_posix(): file_sha256(path)
                            for path in (evaluate, render, Path(__file__), *scoring_code)}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("build", "verify"))
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--matrix-audit", type=Path, help="directory holding the matrix process logs")
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--figures", type=Path, required=True)
    parser.add_argument("--logs", type=Path, required=True)
    parser.add_argument("--step03", type=Path, default=STAGE / "03_data_train")
    parser.add_argument("--raw", type=Path, default=STAGE / "01_data_raw/output")
    args = parser.parse_args()
    print(json.dumps(report(args), sort_keys=True))


if __name__ == "__main__":
    main()
