"""Render Step 06.4 figures twice with only the published plotting CSVs."""

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory


STEP = Path(__file__).resolve().parents[1]
STAGE = STEP.parents[1]
RENDERER = STEP / "execution/render_measured.py"
CSV_NAMES = ("delta-metrics.csv", "fold-scores.csv", "epoch-objectives.csv", "primary-contrasts.csv",
             "primary-contrast-by-fold.csv", "primary-contrast-summary.csv")
PNG_NAMES = ("delta-metrics.png", "delta-metrics-paired-difference.png", "learning-curves.png",
             "primary-contrasts.png")


def hashes(root, names):
    return {str(name): sha256((root / name).read_bytes()).hexdigest() for name in names}


def replay(data, output, receipt_path, scratch, smoke, verify):
    with TemporaryDirectory(prefix="step064-figures-", dir=scratch) as temporary:
        work = Path(temporary)
        csv_only = work / "data"
        csv_only.mkdir()
        for name in CSV_NAMES:
            shutil.copyfile(data / name, csv_only / name)
        rendered = []
        for name in ("first", "second"):
            destination = work / name
            subprocess.run(
                [sys.executable, str(RENDERER), "--data", str(csv_only), "--output", str(destination),
                 *(["--smoke"] if smoke else [])],
                check=True,
                env={**os.environ, "MPLCONFIGDIR": str(work / f"{name}-matplotlib")},
            )
            if {path.name for path in destination.iterdir()} != set(PNG_NAMES):
                raise ValueError("renderer emitted an incomplete or unexpected figure set")
            rendered.append(hashes(destination, PNG_NAMES))
        first, second = rendered
        if first != second:
            raise ValueError("independent CSV-only renders differ")
        code = (RENDERER, Path(__file__), *(STEP / "execution" / name for name in
               ("measured_roster.py", "measured_scores.py")),
                *(STAGE / "util/output/api/cpst_util" / name for name in ("palette.py", "plotting.py")))
        receipt = {"status": "PASS", "smoke": smoke, "independent_csv_only_renders": 2,
                   "csv_sha256": hashes(csv_only, CSV_NAMES), "png_sha256": first,
                   "code_sha256": hashes(STAGE, [path.relative_to(STAGE) for path in code])}
        if verify:
            if first != hashes(output, PNG_NAMES) or receipt != json.loads(receipt_path.read_text(encoding="utf-8")):
                raise ValueError("published plots or figure receipt differ from CSV-only replay")
        else:
            if receipt_path.exists() or any((output / name).exists() for name in PNG_NAMES):
                raise ValueError("figure build requires unused figure and receipt destinations")
            output.mkdir(parents=True, exist_ok=True)
            for name in PNG_NAMES:
                shutil.copyfile(work / "first" / name, output / name)
            receipt_path.parent.mkdir(parents=True, exist_ok=True)
            with receipt_path.open("x", encoding="utf-8") as record:
                json.dump(receipt, record, indent=2, sort_keys=True)
                record.write("\n")
    return {"status": "PASS", "figures": len(PNG_NAMES), "independent_csv_only_renders": 2,
            "receipt_sha256": hashes(receipt_path.parent, (receipt_path.name,))[receipt_path.name]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("build", "verify"))
    for name in ("data", "output", "receipt", "scratch"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    print(json.dumps(replay(args.data, args.output, args.receipt, args.scratch, args.smoke,
                            args.command == "verify"), sort_keys=True))


if __name__ == "__main__":
    main()
