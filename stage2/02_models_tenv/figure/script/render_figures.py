#!/usr/bin/env python3
"""Render the tenv comparison figures from their recorded plot values alone.

A publication script, not scientific code: it reads only the tracked
``figure/data/*.csv`` tables and the shared public plotting implementation under
``01_bases/output/api``.  It never opens this step's ``output/`` tree, matrices,
checkpoints or truth, and it recomputes no metric (ADR 0062, ADR 0063).

    python figure/script/render_figures.py \
        --figure-data stage2/02_models_tenv/figure/data \
        --out-dir stage2/02_models_tenv/figure
"""
import argparse
import sys
from pathlib import Path

STEP_ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = STEP_ROOT.parents[1]
API = REPO_ROOT / "stage2/01_bases/output/api"
if not API.is_dir():
    raise SystemExit(f"public plotting implementation not found: {API}")
sys.path.insert(0, str(API))

from bases import plots  # noqa: E402  (public implementation bound above)

FIGURES = {
    "metrics": ("metrics_values.csv", "metrics.png", plots.replay_metrics),
    "learning_curve": ("learning_curve_values.csv", "learning_curve.png",
                       plots.replay_learning_curve),
}


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--figure-data", default=str(STEP_ROOT / "figure" / "data"),
                        help="directory holding the recorded plot-value CSVs")
    parser.add_argument("--out-dir", default=str(STEP_ROOT / "figure"),
                        help="directory to write the PNGs")
    parser.add_argument("--only", choices=sorted(FIGURES), action="append",
                        help="render only the named figure (repeatable)")
    arguments = parser.parse_args(argv)
    data, out = Path(arguments.figure_data), Path(arguments.out_dir)
    rendered = {}
    for name, (csv_name, png_name, replay) in FIGURES.items():
        if arguments.only and name not in arguments.only:
            continue
        source = data / csv_name
        if not source.is_file():
            raise SystemExit(f"missing recorded plot values: {source}")
        values = replay(source, out / png_name)
        rendered[name] = str(out / png_name)
        print(f"rendered {out / png_name} from {source} "
              f"({len(values['panels'])} panel(s))")
    return 0


if __name__ == "__main__":
    sys.exit(main())
