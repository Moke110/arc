#!/usr/bin/env python3
"""
Stage 04 — Finetune XPert on DTOXS.

Runs the dtoxs finetune of the XPert model using the workflow-owned finetune
directory (configs, splits, and weights all live under finetune/). Consumes the
embedding inputs from data/embeddings/. The base model and training library are
imported from the model component (pert_model/xpert), which stays a pure
library — the workflow owns the finetuned weights.

Output: finetune/output/dtoxs_dtoxs/ (the all-train final prediction model).

Usage:
    python scripts/04_finetune_xpert.py \
        --pretrained-model pert_model/xpert/pretrain/output/pretrain_mdmt_full_50_epoch.pth
"""

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import WORKFLOW_ROOT, XPERT_ROOT  # noqa: E402

DEFAULT_PRETRAINED = (
    XPERT_ROOT / "pretrain/output/pretrain_mdmt_full_50_epoch.pth"
)


def main():
    parser = argparse.ArgumentParser(description="Finetune XPert on DTOXS (workflow-owned)")
    parser.add_argument("--pretrained-model",
                        default=str(DEFAULT_PRETRAINED),
                        help="Pretrained XPert checkpoint")
    parser.add_argument("--config", default="finetune/config/dtoxs",
                        help="Finetune config name (workflow-relative)")
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--epochs", type=int, default=4)
    args = parser.parse_args()

    config_path = WORKFLOW_ROOT / args.config
    if not Path(str(config_path) + ".yaml").exists():
        print(f"ERROR: config not found: {config_path}")
        sys.exit(1)
    if not Path(args.pretrained_model).exists():
        print(f"ERROR: pretrained model not found: {args.pretrained_model}")
        sys.exit(1)

    cmd = [
        sys.executable, "-m", "src.train",
        "--mode", "train",
        "--model", "XPert",
        "--config", str(config_path),
        "--drug_feat", "unimol",
        "--drug_embeds", "unimol,HG",
        "--split_method", "dtoxs",
        "--dataset", "dtoxs",
        "--use_gradscaler", "True",
        "--include_cell_idx", "True",
        "--pretrained_model", args.pretrained_model,
        "--num_epochs", str(args.epochs),
        "--device", args.device,
    ]
    print(f"Running finetune:\n  {' '.join(cmd)}")
    subprocess.run(cmd, cwd=XPERT_ROOT, check=True)

    out_dir = WORKFLOW_ROOT / "finetune/output/dtoxs_dtoxs"
    print(f"Stage 04 complete — finetuned weights in {out_dir}")


if __name__ == "__main__":
    main()
