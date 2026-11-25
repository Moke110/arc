"""Verify current Step 05 figures by rendering twice from the sealed plot CSVs."""

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

LOSS = Path(__file__).resolve().parents[2]
RECEIPT = LOSS / "log/343-figure-render.json"
STUDIES = {
    "05.1_loss_metrics": "330-formal.json",
    "05.2_mse_l1": "333-formal.json",
    "05.3_mse-l1_mse-x": "336-formal.json",
    "05.4_mse-l1_pcc": "339-formal.json",
    "05.5_mse-l1_cos": "339-formal.json",
}
STANDARD_RENDERERS = ("render_metrics.py", "render_learning_curves.py")
SIMILARITY_RENDERERS = ("render.py",)
INPUTS = {
    "delta-metrics.csv": "figure/data/delta-metrics.csv",
    "epoch-objectives.csv": "figure/data/epoch-objectives.csv",
    "fold-scores.csv": "log/fold-scores.csv",
}
FIGURES = {"delta-metrics.png", "delta-metrics-improvement.png", "learning-curves.png"}


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def inventory(work):
    work.mkdir()
    for name in ("matplotlib", "scratch"):
        (work / name).mkdir()
    env = {**os.environ, "MPLCONFIGDIR": str(work / "matplotlib"), "TMPDIR": str(work / "scratch")}
    csv_hashes = {}
    figure_hashes = {}
    renderer_hashes = {}

    for study, formal_receipt in STUDIES.items():
        release = LOSS / study
        prior = json.loads((release / "log" / formal_receipt).read_text())
        data = work / study / "data"
        data.mkdir(parents=True)
        for name, formal_path in INPUTS.items():
            csv_path = release / "figure/data" / name
            key = csv_path.relative_to(LOSS).as_posix()
            csv_hashes[key] = digest(csv_path)
            if csv_hashes[key] != prior["published_files_sha256"][formal_path]:
                raise ValueError(f"plot CSV differs from the original formal receipt: {key}")
            shutil.copyfile(csv_path, data / name)

        renderers = SIMILARITY_RENDERERS if study.startswith(("05.4", "05.5")) else STANDARD_RENDERERS
        scripts = [release / "figure/script" / name for name in renderers]
        for script in scripts:
            renderer_hashes[script.relative_to(LOSS).as_posix()] = digest(script)
        for attempt in ("first", "second"):
            output = work / study / attempt
            for script in scripts:
                subprocess.run([sys.executable, str(script), "--data", str(data), "--output", str(output)],
                               cwd=work, env=env, check=True)
            if {path.name for path in output.iterdir()} != FIGURES:
                raise ValueError(f"unexpected figure inventory for {study}")
        for name in sorted(FIGURES):
            path = release / "figure" / name
            key = path.relative_to(LOSS).as_posix()
            figure_hashes[key] = digest(work / study / "first" / name)
            if figure_hashes[key] != digest(work / study / "second" / name):
                raise ValueError(f"independent CSV-only renders differ: {key}")
            if figure_hashes[key] != digest(path):
                raise ValueError(f"published PNG differs from the CSV-only render: {key}")
        for name in INPUTS:
            key = f"{study}/figure/data/{name}"
            if csv_hashes[key] != digest(data / name) or csv_hashes[key] != digest(LOSS / key):
                raise ValueError(f"plot CSV changed during rendering: {key}")

    for source in ("figure/script/loss_plot_style.py", "similarity/render_similarity.py",
                   "05.3_mse-l1_mse-x/execution/sweep_protocol.py"):
        renderer_hashes[source] = digest(LOSS / source)
    return {
        "schema": "stage2-figure-render-v2",
        "ticket": 343,
        "independent_renders": 2,
        "environment": {
            "python": python_version(),
            "matplotlib": version("matplotlib"),
            "numpy": version("numpy"),
        },
        "csv_files_sha256": csv_hashes,
        "png_files_sha256": figure_hashes,
        "renderer_files_sha256": renderer_hashes,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True,
                        help="unused scratch directory under the checkout's tmp/")
    args = parser.parse_args()
    actual = inventory(args.work.resolve())
    expected = json.loads(RECEIPT.read_text())
    mismatches = sorted(key for key in actual.keys() | expected.keys()
                        if actual.get(key) != expected.get(key))
    if mismatches:
        raise ValueError(f"Step 05 rendering receipt differs in {', '.join(mismatches)}")
    print(f"Verified {len(actual['png_files_sha256'])} Step 05 PNGs from "
          f"{len(actual['csv_files_sha256'])} sealed CSVs")


if __name__ == "__main__":
    main()
