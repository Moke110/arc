"""Verify the five current Step 05 comparison publications from their plotted CSVs."""

import argparse
import json
import os
from pathlib import Path
import sys

from loss_comparison import LOSS, RENDERER, STUDIES
from stage2_util.comparison import FIGURES

sys.path.insert(0, str(LOSS.parent / "util/figure/script"))

from verify_report_figures import digest, verify_report


def verify_figures(work):
    work.mkdir()
    for child in ("matplotlib", "scratch"):
        (work / child).mkdir()
    env = {**os.environ, "MPLCONFIGDIR": str(work / "matplotlib"), "TMPDIR": str(work / "scratch")}
    count = 0
    for name, study in STUDIES.items():
        step = LOSS / name
        receipt_name = "figure-render.json"
        receipt = json.loads((step / "log" / receipt_name).read_text())
        evidence = receipt["evidence"]
        if (set(receipt["png_files_sha256"]) != {f"figure/{figure}" for figure in FIGURES}
                or evidence["reference"] != study.reference
                or evidence["fits"] != len(study.arms) * 3
                or evidence["acceptance_receipt_sha256"] != digest(step / "log" / study.receipt)):
            raise ValueError(f"current figure receipt does not match the accepted study {name}")
        count += verify_report(step, receipt_name, (("figure", RENDERER),), work / name, env)
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True, help="unused directory under the checkout's tmp/")
    args = parser.parse_args()
    count = verify_figures(args.work.resolve())
    print(f"Verified {count} Step 05 PNGs from sealed plotted CSVs")


if __name__ == "__main__":
    main()
