#!/usr/bin/env python
"""Thin CLI: (re)generate ``model_card.md`` for a run directory.

    python execution/make_model_card.py --run-dir <run-dir> [--eval-dir <eval-dir>]

Reads only the run's own artifacts, so it can regenerate a card after an
evaluation without retraining.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ptlib import modelcard  # noqa: E402


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate a PRnet formal-transfer model card from a run directory.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--eval-dir", type=Path, default=None,
                        help="optional evaluation directory (aggregate.json)")
    parser.add_argument("--out", type=Path, default=None,
                        help="default <run-dir>/model_card.md")
    return parser


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ns = build_parser().parse_args(argv)
    path = modelcard.write_model_card(ns.run_dir, eval_dir=ns.eval_dir, out_path=ns.out)
    print(str(path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
