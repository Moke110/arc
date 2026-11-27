"""Regenerate the Step 10 panels twice with only the plotted CSVs available."""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

from dose_protocol import STEP
from dose_report import PLOT_TABLES
from dose_runs import axis_run
from stage2_util.comparison import FIGURES

RENDERER = STEP / "figure/script/render_doses.py"


def replay(data, output, scratch, *, verify=False):
    with TemporaryDirectory(prefix="step10-figures-", dir=scratch) as temporary:
        work = Path(temporary)
        csv_only = work / "data"
        csv_only.mkdir()
        for filename in PLOT_TABLES:
            shutil.copyfile(data / filename, csv_only / filename)
        renders = []
        for attempt in ("first", "second"):
            rendered = work / attempt
            subprocess.run([sys.executable, str(RENDERER), "--data", str(csv_only), "--output", str(rendered)],
                           check=True, cwd=work, env={**os.environ, "MPLCONFIGDIR": str(work / f"{attempt}-mpl")})
            if {path.name for path in rendered.iterdir()} != set(FIGURES):
                raise ValueError("Step 10 renderer emitted an unexpected set of figures")
            renders.append(axis_run.file_hashes(rendered, FIGURES))
        if renders[0] != renders[1]:
            raise ValueError("independent Step 10 CSV-only renders differ")
        if verify:
            if ({path.name for path in output.glob("*.png")} != set(FIGURES)
                    or axis_run.file_hashes(output, FIGURES) != renders[0]):
                raise ValueError("published Step 10 figures differ from their plotted CSVs")
        else:
            output.mkdir(parents=True, exist_ok=True)
            if any((output / filename).exists() for filename in FIGURES):
                raise ValueError("Step 10 figure build requires unused PNG destinations")
            for filename in FIGURES:
                shutil.copyfile(work / "first" / filename, output / filename)
        return {"status": "PASS", "independent_csv_only_renders": len(renders),
                "csv_sha256": axis_run.file_hashes(csv_only, PLOT_TABLES), "png_sha256": renders[0]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("build", "verify"))
    for name in ("data", "output", "scratch"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    result = replay(args.data, args.output, args.scratch, verify=args.mode == "verify")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
