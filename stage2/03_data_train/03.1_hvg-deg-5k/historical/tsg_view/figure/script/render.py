#!/usr/bin/env python3
"""Canonical figure entry point for the `tsg_view` study.

The implementation lives in ``execution/render.py`` so the figure stage can be unit
tested with the rest of the execution code; this wrapper is the documented path::

    python figure/script/render.py --data figure/data --output figure --suite tenv

The renderer reads only ``figure/data/*.csv`` and draws stored values -- no matrices,
checkpoints or truth products.  A CSV-only replay copies ``figure/data/`` plus
``execution/render.py`` into a fresh tree and runs it under ``python -I`` with an
isolated ``MPLCONFIGDIR``; the PNGs are then byte-identical (see the renderer test).
"""
import sys
from pathlib import Path

STEP = Path(__file__).resolve().parents[2]
EXECUTION = STEP / "execution"
if str(EXECUTION) not in sys.path:
    sys.path.insert(0, str(EXECUTION))

import render  # noqa: E402

if __name__ == "__main__":
    sys.exit(render.main())
