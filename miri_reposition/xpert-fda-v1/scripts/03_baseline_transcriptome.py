#!/usr/bin/env python3
"""
Stage 03 — Baseline transcriptome.

Builds the consensus cardiomyocyte baseline transcriptome from DTOXS control
wells and caches it under data/baseline/. This is workflow-owned: the baseline
is generated and cached inside xpert-fda-v1 rather than in the model component.
The implementation lives in the moved predict/build_baseline.py, which this
stage invokes; the output is copied to data/baseline/cm_ctl_landmark.csv.

Usage:
    python scripts/03_baseline_transcriptome.py
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import WORKFLOW_ROOT, XPERT_ROOT  # noqa: E402

OUT_DIR = WORKFLOW_ROOT / "data" / "baseline"


def main():
    parser = argparse.ArgumentParser(description="Build the workflow-owned baseline transcriptome")
    parser.add_argument("--force", action="store_true",
                        help="Rebuild even if the cached baseline exists")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / "cm_ctl_landmark.csv"
    if out_path.exists() and not args.force:
        print(f"Cached baseline present: {out_path} — skipping")
        return

    # predict/build_baseline.py writes to predict/fda-cm/cm_ctl_landmark.csv
    build_script = WORKFLOW_ROOT / "predict" / "build_baseline.py"
    print(f"Running {build_script}")
    subprocess.run([sys.executable, str(build_script)], cwd=XPERT_ROOT, check=True)

    src = WORKFLOW_ROOT / "predict" / "fda-cm" / "cm_ctl_landmark.csv"
    if not src.exists():
        print("ERROR: baseline not produced by build_baseline.py")
        sys.exit(1)
    shutil.copy2(src, out_path)
    print(f"Stage 03 complete — baseline cached at {out_path}")


if __name__ == "__main__":
    main()
