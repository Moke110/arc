"""Verify the 22 current Stage 2 PNGs with their render receipts and plot data."""

import argparse
from hashlib import sha256
from importlib.metadata import version
import json
import os
from pathlib import Path
from platform import python_version
import shutil
import subprocess
import sys

STAGE = Path(__file__).resolve().parents[3]
AXIS_REPORTS = {
    "02_gene_axis": (
        ("02.1_hvg-deg/figure", "02.1_hvg-deg/figure/script/render_coverage.py"),
        ("02.2_denoise-deg/figure", "02.2_denoise-deg/figure/script/render.py"),
    ),
    "04_tsg_axis": (("figure/canonical-axis", "figure/script/render.py"),),
}


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def check_hashes(root, hashes):
    for path, expected in hashes.items():
        if digest(root / path) != expected:
            raise ValueError(f"render receipt hash differs for {root / path}")


def verify_axis_report(step, renderers, work, env):
    receipt = json.loads((step / "log/342-figure-render.json").read_text())
    if (receipt["schema"], receipt["ticket"], receipt["independent_renders"]) != (
        "stage2-figure-render-v1", 342, 2,
    ):
        raise ValueError(f"unknown figure render receipt for {step}")
    environment = {name: version(name) for name in ("matplotlib", "numpy", "pandas")}
    if receipt["environment"] != {**environment, "python": python_version()}:
        raise ValueError(f"render environment differs from the receipt for {step}")

    csvs = receipt["csv_files_sha256"]
    pngs = receipt["png_files_sha256"]
    scripts = receipt["renderer_files_sha256"]
    if set(scripts) != {script for _, script in renderers}:
        raise ValueError(f"renderer inventory differs from the receipt for {step}")
    published = {path.relative_to(step).as_posix() for directory, _ in renderers
                 for path in (step / directory).glob("*.png")}
    if set(pngs) != published:
        raise ValueError(f"published PNG inventory differs from the receipt for {step}")
    for hashes in (csvs, pngs, scripts):
        check_hashes(step, hashes)

    for attempt in ("first", "second"):
        scratch = work / attempt
        for path in csvs:
            target = scratch / path
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(step / path, target)
        for directory, script in renderers:
            output = scratch / directory
            subprocess.run([sys.executable, str(step / script), "--data", str(output / "data"),
                            "--output", str(output)], cwd=work, env=env, check=True)
        actual = {path.relative_to(scratch).as_posix() for path in scratch.rglob("*.png")}
        if actual != set(pngs):
            raise ValueError(f"rendered PNG inventory differs from the receipt for {step}")
        check_hashes(scratch, pngs)
        check_hashes(scratch, csvs)
    check_hashes(step, csvs)
    return len(pngs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True,
                        help="unused scratch directory under the checkout's tmp/")
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir()
    for child in ("matplotlib", "scratch"):
        (work / child).mkdir()
    env = {**os.environ, "MPLCONFIGDIR": str(work / "matplotlib"), "TMPDIR": str(work / "scratch")}
    count = sum(verify_axis_report(STAGE / step, paths, work / step, env)
                for step, paths in AXIS_REPORTS.items())
    if count != 7:
        raise ValueError(f"expected seven Step 02 and Step 04 figures, verified {count}")
    loss_verifier = STAGE / "05_tsg_loss/figure/script/verify_figures.py"
    subprocess.run([sys.executable, str(loss_verifier), "--work", str(work / "loss")],
                   cwd=work, env=env, check=True)
    print("Verified all 22 current Stage 2 PNGs")


if __name__ == "__main__":
    main()
