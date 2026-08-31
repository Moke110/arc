#!/usr/bin/env python3
"""
Step 07 — Perturbed transcriptome.

Runs XPert predict on the FDA drug list to produce predicted drug perturbation
(deg = perturbed - control) across the 978 L1000 landmark genes. Consumes the
baseline (02_cm_baseline/output/cm_ctl_landmark.csv), the drug list
(01_drug_embed/output/drug_list.csv), the finetuned model
(06_xpert_finetune/output/dtoxs_dtoxs/), and the frozen drug embeddings
(input/ + 01_drug_embed/output/ via the step-06 config). XPert predict
indexes the precomputed UniMol array by drug index so finetune and predict
consume identical drug features.

Output: 07_xpert_predict/output/xpert_fda-cm_deg.csv

Usage:
    python 07_xpert_predict/perturb_transcriptome.py
"""

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from utils import WORKFLOW_ROOT, XPERT_ROOT  # noqa: E402

STEP_DIR = WORKFLOW_ROOT / "07_xpert_predict"
OUTPUT = STEP_DIR / "output" / "xpert_fda-cm_deg.csv"


def main():
    parser = argparse.ArgumentParser(description="Predict perturbed transcriptomes for FDA drugs")
    parser.add_argument("--model-path",
                        default=str(WORKFLOW_ROOT / "06_xpert_finetune/output/dtoxs_dtoxs"),
                        help="Finetuned XPert model (dir of .pth or single .pth)")
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--dose-idx", type=int, default=0)
    parser.add_argument("--time-idx", type=int, default=2)
    parser.add_argument("--force", action="store_true",
                        help="Re-predict even if the output DEG table already exists")
    args = parser.parse_args()

    if not args.force and OUTPUT.exists():
        print(f"Cached DEG table present: {OUTPUT} — skipping "
              f"(pass --force to re-predict; re-prediction from retrained weights "
              f"is GPU-fp-nondeterministic across environments and breaks parity "
              f"with the archived reference)")
        return

    prediction_root = STEP_DIR / "output"
    prediction_root.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable, "-m", "src.predict",
        "--prediction_root", str(prediction_root),
        "--ctl_csv", str(WORKFLOW_ROOT / "02_cm_baseline/output/cm_ctl_landmark.csv"),
        "--drug_csv", str(WORKFLOW_ROOT / "01_drug_embed/output/drug_list.csv"),
        "--output_csv", "xpert_fda-cm_deg.csv",
        "--model_path", args.model_path,
        "--dose_idx", str(args.dose_idx),
        "--time_idx", str(args.time_idx),
        "--device", args.device,
        "--batch_size", str(args.batch_size),
        "--config", str(WORKFLOW_ROOT / "06_xpert_finetune/config/dtoxs"),
    ]
    print(f"Running predict:\n  {' '.join(cmd)}")
    subprocess.run(cmd, cwd=XPERT_ROOT, check=True)

    print(f"Step 07 complete — {OUTPUT}")


if __name__ == "__main__":
    main()
