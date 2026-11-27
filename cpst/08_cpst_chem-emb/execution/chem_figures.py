"""Render twice from an isolated CSV copy and compare the complete figure set."""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

from chem_runs import STEP, axis_run
from cpst_util.comparison import FIGURES as COMPARISON_FIGURES

RENDERER = STEP / "execution/render_chem.py"
FIGURES = ("delta-metrics.png", "delta-metrics-paired-difference.png")
COMPARISON_RENDERER = STEP / "figure/script/render_comparison.py"
COMPARISON_TABLES = ("learning-curve.csv", "metrics.csv", "metrics-paired.csv")


def replay(data, output, scratch, verify=False, *, comparison=False):
    renderer, tables, figures = ((COMPARISON_RENDERER, COMPARISON_TABLES, COMPARISON_FIGURES)
                                  if comparison else (RENDERER, ("fold-scores.csv",), FIGURES))
    with TemporaryDirectory(prefix="step08-figures-", dir=scratch) as temporary:
        work = Path(temporary)
        csv_only = work / "data"
        csv_only.mkdir()
        for table in tables:
            shutil.copyfile(data / table, csv_only / table)
        rendered = []
        for name in ("first", "second"):
            destination = work / name
            subprocess.run([sys.executable, str(renderer), "--data", str(csv_only), "--output", str(destination)],
                           check=True, env={**os.environ, "MPLCONFIGDIR": str(work / f"{name}-mpl")})
            if {path.name for path in destination.iterdir()} != set(figures):
                raise ValueError("renderer emitted an unexpected Step 08 figure set")
            rendered.append(axis_run.file_hashes(destination, figures))
        if rendered[0] != rendered[1]:
            raise ValueError("independent CSV-only figure renders differ")
        if verify:
            if ({path.name for path in output.glob("*.png")} != set(figures)
                    or axis_run.file_hashes(output, figures) != rendered[0]):
                raise ValueError("published Step 08 figures differ from the CSV-only replay")
        else:
            if any((output / name).exists() for name in figures):
                raise ValueError("figure build requires unused PNG destinations")
            output.mkdir(exist_ok=True, parents=True)
            for name in figures:
                shutil.copyfile(work / "first" / name, output / name)
        return {"status": "PASS", "independent_csv_only_renders": 2,
                "csv_sha256": axis_run.file_hashes(csv_only, tables),
                "png_sha256": rendered[0]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("build", "verify"))
    for name in ("data", "output", "scratch"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--comparison", action="store_true", help="replay the canonical Step 08 figure set")
    args = parser.parse_args()
    print(json.dumps(replay(args.data, args.output, args.scratch, args.mode == "verify",
                            comparison=args.comparison), sort_keys=True))


if __name__ == "__main__":
    main()
