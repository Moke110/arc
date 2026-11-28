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
AXIS_RENDERERS = (
    ("02.1_hvg-deg/figure", "02.1_hvg-deg/figure/script/render_coverage.py"),
    ("02.2_denoise-deg/figure", "02.2_denoise-deg/figure/script/render.py"),
)
COMPARISON_RENDERERS = {
    "04_tsg_axis": "render.py",
    "05_tsg_module": "render.py",
    "08_cpst_chem-emb": "render_comparison.py",
    "09_cpst_axis": "render_axes.py",
    "10_cpst_dose": "render_doses.py",
}
HISTORICAL_COMPARISON_SOURCE = "325be5dd45e2b670f1e01882e225aa5debbeaab7"
SHARED_COMPARISON = "util/output/api/cpst_util/comparison.py"
COMPARISON_NAMES = ("learning-curve", "metrics", "metrics-paired")


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def check_hashes(root, hashes):
    for path, expected in hashes.items():
        if digest(root / path) != expected:
            raise ValueError(f"render receipt hash differs for {root / path}")


def check_renderer_hashes(root, scripts):
    for path, expected in scripts.items():
        if digest(root / path) == expected:
            continue
        if root == STAGE and path == SHARED_COMPARISON:
            original = subprocess.check_output(
                ["git", "show", f"{HISTORICAL_COMPARISON_SOURCE}:cpst/{path}"], cwd=STAGE)
            if sha256(original).hexdigest() == expected:
                continue
        raise ValueError(f"render receipt hash differs for {root / path}")


def verify_report(step, receipt_name, renderers, work, env, *, comparison=False):
    receipt = json.loads((step / "log" / receipt_name).read_text())
    if (receipt["schema"], receipt["status"], receipt["independent_renders"]) != (
        "cpst-figure-render-v1", "PASS", 2,
    ):
        raise ValueError(f"unknown figure render receipt for {step}")
    environment = {name: version(name) for name in receipt["environment"] if name != "python"}
    if receipt["environment"] != {**environment, "python": python_version()}:
        raise ValueError(f"render environment differs from the receipt for {step}")

    csvs = receipt["csv_files_sha256"]
    pngs = receipt["png_files_sha256"]
    if comparison and (
        set(pngs) != {f"figure/{name}.png" for name in COMPARISON_NAMES}
        or set(csvs) != {f"figure/data/{name}.csv" for name in COMPARISON_NAMES}
    ):
        raise ValueError(f"canonical comparison files differ for {step}")
    if comparison:
        expected_renderers = {Path(script).name for _, script in renderers}
        actual_renderers = {path.name for path in (step / "figure/script").glob("*.py")}
        if actual_renderers != expected_renderers:
            raise ValueError(f"unexpected comparison renderer in {step}: {sorted(actual_renderers - expected_renderers)}")
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
    check_renderer_hashes(script_root, scripts)
    projected_tables = receipt.get("evidence", {}).get("projected_scientific_tables_sha256", {})
    check_hashes(step, projected_tables)

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
    axis_env = {name: value for name, value in env.items() if name != "PYTHONPATH"}
    paths = [str(STAGE / "util/output/api")]
    if inherited := env.get("PYTHONPATH"):
        paths.append(inherited)
    env["PYTHONPATH"] = os.pathsep.join(paths)
    axis_count = verify_report(STAGE / "02_gene_axis", "342-figure-render.json",
                               AXIS_RENDERERS, work / "02_gene_axis", axis_env)
    if axis_count != 5:
        raise ValueError(f"expected five Step 02 figures, verified {axis_count}")
    for step, script in COMPARISON_RENDERERS.items():
        renderers = (("figure", f"figure/script/{script}"),)
        verify_report(STAGE / step, "figure-render.json", renderers, work / step, env,
                      comparison=True)
    palette = STAGE / "util/figure/palette.png"
    for attempt in ("first", "second"):
        output = work / f"palette-{attempt}.png"
        subprocess.run([sys.executable, str(palette.parent / "script/render_palette.py"),
                        "--output", str(output)], cwd=work, env=env, check=True)
        if digest(output) != digest(palette):
            raise ValueError("published palette differs from its render")

    directories = [*STAGE.glob("*/figure"), *STAGE.glob("*/*/figure")]
    for directory in directories:
        check_layout(directory)
    published_count = sum(1 for directory in directories for _ in directory.glob("*.png"))
    if published_count != 21:
        raise ValueError(f"expected 21 published Stage 2 PNGs, found {published_count}")
    print("Verified all 21 current CPST PNGs, including five comparison sets")


if __name__ == "__main__":
    main()
