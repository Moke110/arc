#!/usr/bin/env python3
"""
Stage 05 — Perturbed transcriptomes.

Runs XPert predict on the FDA drug list to produce predicted drug perturbation
(deg = perturbed - control) across the 978 L1000 landmark genes. Consumes the
workflow baseline (data/baseline/cm_ctl_landmark.csv), the workflow finetuned
model (finetune/output/dtoxs_dtoxs/), and the frozen drug embeddings
(data/embeddings/). XPert predict indexes the precomputed UniMol array by drug
index so finetune and predict consume identical drug features.

Output: data/deg/xpert_fda-cm_deg.csv

Usage:
    python scripts/05_perturbed_transcriptomes.py
"""

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import WORKFLOW_ROOT, XPERT_ROOT  # noqa: E402

OUTPUT = WORKFLOW_ROOT / "data" / "deg" / "xpert_fda-cm_deg.csv"


def main():
    parser = argparse.ArgumentParser(description="Predict perturbed transcriptomes for FDA drugs")
    parser.add_argument("--model-path",
                        default=str(WORKFLOW_ROOT / "finetune/output/dtoxs_dtoxs"),
                        help="Finetuned XPert model (dir of .pth or single .pth)")
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--dose-idx", type=int, default=0)
    parser.add_argument("--time-idx", type=int, default=2)
    args = parser.parse_args()

    prediction_root = WORKFLOW_ROOT / "data" / "deg"
    prediction_root.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable, "-m", "src.predict",
        "--prediction_root", str(prediction_root),
        "--ctl_csv", str(WORKFLOW_ROOT / "data/baseline/cm_ctl_landmark.csv"),
        "--drug_csv", str(WORKFLOW_ROOT / "data/drug_list.csv"),
        "--output_csv", "xpert_fda-cm_deg.csv",
        "--model_path", args.model_path,
        "--dose_idx", str(args.dose_idx),
        "--time_idx", str(args.time_idx),
        "--device", args.device,
        "--batch_size", str(args.batch_size),
        "--config", str(WORKFLOW_ROOT / "finetune/config/dtoxs"),
    ]
    print(f"Running predict:\n  {' '.join(cmd)}")
    subprocess.run(cmd, cwd=XPERT_ROOT, check=True)

    print(f"Stage 05 complete — {OUTPUT}")


if __name__ == "__main__":
    main()
