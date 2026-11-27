"""Verify that the seven figures reproduce exactly from their isolated CSVs."""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

from axis_protocol import STEP
from axis_report import CSV_FIGURES
from axis_runs import axis_run

RENDERER = STEP / "figure/script/render_axes.py"
FIGURES = ("learning-curves", *CSV_FIGURES)


def replay(data, output, scratch, *, verify=False):
    """Render twice in fresh processes and compare outputs before publishing."""
    with TemporaryDirectory(prefix="step09-figures-", dir=scratch) as temporary:
        work = Path(temporary)
        csv_only = work / "data"
        csv_only.mkdir()
        for name in FIGURES:
            shutil.copyfile(data / f"{name}.csv", csv_only / f"{name}.csv")
        files = tuple(f"{name}.png" for name in FIGURES)
        rendered = []
        for name in ("first", "second"):
            place = work / name
            subprocess.run([sys.executable, str(RENDERER), "--data", str(csv_only), "--output", str(place)],
                           check=True, cwd=work, env={**os.environ, "MPLCONFIGDIR": str(work / f"{name}-mpl")})
            if {path.name for path in place.iterdir()} != set(files):
                raise ValueError("renderer emitted an unexpected set of figures")
            rendered.append(axis_run.file_hashes(place, files))
        if rendered[0] != rendered[1]:
            raise ValueError("independent CSV-only renders differ")
        if verify:
            if axis_run.file_hashes(output, files) != rendered[0]:
                raise ValueError("published figures differ from their CSV-only reproduction")
        else:
            output.mkdir(parents=True, exist_ok=True)
            if any((output / name).exists() for name in files):
                raise ValueError("figure build requires unused PNG destinations")
            for name in files:
                shutil.copyfile(work / "first" / name, output / name)
        return {"status": "PASS", "independent_csv_only_renders": 2,
                "csv_sha256": axis_run.file_hashes(csv_only, (f"{name}.csv" for name in FIGURES)),
                "png_sha256": rendered[0]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("build", "verify"))
    for name in ("data", "output", "scratch"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(replay(args.data, args.output, args.scratch, verify=args.mode == "verify"), sort_keys=True))


if __name__ == "__main__":
    main()
