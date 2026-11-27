"""Replay canonical comparison PNGs and bind their plotted data and renderers."""

from importlib.metadata import version
from pathlib import Path
from platform import python_version
import subprocess
import sys
from tempfile import TemporaryDirectory

from .comparison import FIGURES
from .derivation.release import sha256

PLOTTED_CSVS = ("learning-curve.csv", "metrics.csv", "metrics-paired.csv")


def render_receipt(stage, step, output, renderer, *, dependencies=()):
    """Render twice from the canonical CSVs, then return a figure replay receipt."""
    figure = output / "figure"
    data = figure / "data"
    csvs = {f"figure/data/{name}": sha256(data / name) for name in PLOTTED_CSVS}
    with TemporaryDirectory(prefix="figure-replay-", dir=output) as replay:
        for candidate in (figure, Path(replay)):
            subprocess.run([sys.executable, str(renderer), "--data", str(data),
                            "--output", str(candidate)], check=True)
        pngs = {f"figure/{name}": sha256(figure / name) for name in FIGURES}
        if pngs != {f"figure/{name}": sha256(Path(replay) / name) for name in FIGURES}:
            raise ValueError(f"independent figure renders differ for {step}")

    if csvs != {f"figure/data/{name}": sha256(data / name) for name in PLOTTED_CSVS}:
        raise ValueError(f"plot data changed during figure replay for {step}")
    renderers = (renderer, *dependencies, *(Path(__file__).with_name(f"{name}.py")
                                             for name in ("comparison", "palette", "plotting")))
    return {
        "schema": "cpst-figure-render-v1",
        "status": "PASS",
        "step": step,
        "independent_renders": 2,
        "environment": {"python": python_version(), **{name: version(name) for name in ("matplotlib", "numpy", "pandas")}},
        "csv_files_sha256": csvs,
        "png_files_sha256": pngs,
        "renderer_files_sha256": {path.relative_to(stage).as_posix(): sha256(path) for path in renderers},
    }
