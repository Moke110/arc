"""Verify every current Stage 2 PNG against its render receipt and plot data."""

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
REPORTS = {
    "02_gene_axis": ("342-figure-render.json", (
        ("02.1_hvg-deg/figure", "02.1_hvg-deg/figure/script/render_coverage.py"),
        ("02.2_denoise-deg/figure", "02.2_denoise-deg/figure/script/render.py"),
    )),
    "04_tsg_axis": ("342-figure-render.json", (("figure", "figure/script/render.py"),)),
    "06_tsg_shape/06.1_round1": ("figure-render.json", (("figure", "figure/script/render_shape.py"),)),
    "06_tsg_shape/06.2_round2": ("figure-render.json", (("figure", "figure/script/render_architecture.py"),)),
    "06_tsg_shape/06.3_round3": ("figure-render.json", (("figure", "figure/script/render_composition.py"),)),
    "06_tsg_shape/06.4_round4": ("figure-render.json", (("figure", "figure/script/render_measured.py"),)),
    "08_cpsv_chem-emb": ("figure-render.json", (("figure", "figure/script/render_chem.py"),)),
    "09_cpsv_axis": ("figure-render.json", (("figure", "figure/script/render_axes.py"),)),
}


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def check_hashes(root, hashes):
    for path, expected in hashes.items():
        if digest(root / path) != expected:
            raise ValueError(f"render receipt hash differs for {root / path}")


def verify_report(step, receipt_name, renderers, work, env):
    receipt = json.loads((step / "log" / receipt_name).read_text())
    if (receipt["schema"], receipt["status"], receipt["independent_renders"]) != (
        "stage2-figure-render-v1", "PASS", 2,
    ):
        raise ValueError(f"unknown figure render receipt for {step}")
    environment = {name: version(name) for name in receipt["environment"] if name != "python"}
    if receipt["environment"] != {**environment, "python": python_version()}:
        raise ValueError(f"render environment differs from the receipt for {step}")

    csvs = receipt["csv_files_sha256"]
    pngs = receipt["png_files_sha256"]
    scripts = receipt["renderer_files_sha256"]
    script_root = step if receipt_name == "342-figure-render.json" else STAGE
    paths = {(step / script).relative_to(script_root).as_posix() for _, script in renderers}
    if not paths <= set(scripts):
        raise ValueError(f"renderer inventory differs from the receipt for {step}")
    published = {path.relative_to(step).as_posix() for directory, _ in renderers
                 for path in (step / directory).glob("*.png")}
    if set(pngs) != published:
        raise ValueError(f"published PNG inventory differs from the receipt for {step}")
    for hashes in (csvs, pngs):
        check_hashes(step, hashes)
    check_hashes(script_root, scripts)

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


def check_layout(directory):
    extra = {path.name for path in directory.iterdir() if path.is_dir()} - {"data", "script"}
    if extra:
        raise ValueError(f"unexpected directories in {directory}: {sorted(extra)}")
    for child in ("data", "script"):
        if any(path.is_dir() and path.name != "__pycache__" for path in (directory / child).glob("*")):
            raise ValueError(f"nested figure set in {directory / child}")


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
    paths = [str(STAGE / "util/output/api")]
    if inherited := env.get("PYTHONPATH"):
        paths.append(inherited)
    env["PYTHONPATH"] = os.pathsep.join(paths)
    count = sum(verify_report(STAGE / step, receipt, renderers, work / step, env)
                for step, (receipt, renderers) in REPORTS.items())
    if count != 30:
        raise ValueError(f"expected 30 Step 02, Step 04, Step 06, Step 08 and Step 09 figures, verified {count}")
    loss_verifier = STAGE / "05_tsg_loss/figure/script/verify_figures.py"
    subprocess.run([sys.executable, str(loss_verifier), "--work", str(work / "loss")],
                   cwd=work, env=env, check=True)
    palette = STAGE / "util/figure/palette.png"
    for attempt in ("first", "second"):
        output = work / f"palette-{attempt}.png"
        subprocess.run([sys.executable, str(palette.parent / "script/render_palette.py"),
                        "--output", str(output)], cwd=work, env=env, check=True)
        if digest(output) != digest(palette):
            raise ValueError("published palette differs from its render")

    for directory in STAGE.glob("*/figure"):
        check_layout(directory)
    for directory in STAGE.glob("*/*/figure"):
        check_layout(directory)
    print("Verified all 46 current Stage 2 PNGs")


if __name__ == "__main__":
    main()
