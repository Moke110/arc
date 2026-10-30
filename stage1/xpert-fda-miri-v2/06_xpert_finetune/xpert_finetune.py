#!/usr/bin/env python3
"""
Step 06 — Finetune XPert on the union CM data.

Runs the cm finetune of the XPert model (cm_full all-train split) using the
step-owned config (06_xpert_finetune/config/cm.yaml). Consumes the embedding
inputs from input/ + 01_drug_embed/output/ and the union finetune dataset at
03_cm_finetune_data/output/cm_union.h5ad (built by step 03 from all five
pert_data datasets — the substantive v2 change vs v1's frozen DToxS
`stage1/xpert-fda-miri-v1/03_cm_finetune_data/output/dtoxs.h5ad`). The HG table is the author reference (step 01 copies it
byte-identical — no HG training), and the finetune uses the same 978-gene set
verified by _lib/verify_genes.py. The base model and training library are
imported from the model component (pert_model/XPert), which stays a pure
library — the workflow owns the pretrained base weight (04_xpert_pretrain/)
and the finetuned weights.

XPert writes training artifacts to <workflow>/finetune/output/cm_cm_full/
(the library derives this from the config root); this step moves them into
06_xpert_finetune/output/cm_cm_full/ when training finishes.

Output: 06_xpert_finetune/output/cm_cm_full/ (the all-train final prediction model).

Usage:
    python 06_xpert_finetune/xpert_finetune.py
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from utils import WORKFLOW_ROOT, XPERT_ROOT  # noqa: E402

STEP_DIR = WORKFLOW_ROOT / "06_xpert_finetune"
OUT_DIR = STEP_DIR / "output" / "cm_cm_full"
STAGING_DIR = WORKFLOW_ROOT / "finetune" / "output" / "cm_cm_full"
DEFAULT_PRETRAINED = (
    WORKFLOW_ROOT / "04_xpert_pretrain" / "pretrain_mdmt_full_50_epoch.pth"
)


def main():
    parser = argparse.ArgumentParser(description="Finetune XPert on the union CM data (workflow-owned)")
    parser.add_argument("--pretrained-model",
                        default=str(DEFAULT_PRETRAINED),
                        help="Pretrained XPert checkpoint")
    parser.add_argument("--config", default=str(STEP_DIR / "config/cm"),
                        help="Finetune config path (without .yaml)")
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--epochs", type=int, default=4)
    parser.add_argument("--force", action="store_true",
                        help="Retrain even if finetuned weights already exist")
    args = parser.parse_args()

    if not args.force and list(OUT_DIR.glob("*.pth")):
        print(f"Cached finetuned weights present: {OUT_DIR} — skipping "
              f"(pass --force to retrain; retraining is GPU-fp-nondeterministic "
              f"across environments and breaks parity with archived predictions)")
        return

    if not Path(str(args.config) + ".yaml").exists():
        print(f"ERROR: config not found: {args.config}")
        sys.exit(1)
    if not Path(args.pretrained_model).exists():
        print(f"ERROR: pretrained model not found: {args.pretrained_model}")
        sys.exit(1)

    shutil.rmtree(STAGING_DIR, ignore_errors=True)
    cmd = [
        sys.executable, "-m", "src.train",
        "--mode", "train",
        "--model", "XPert",
        "--config", str(args.config),
        "--drug_feat", "unimol",
        "--drug_embeds", "unimol,HG",
        "--split_method", "cm_full",
        "--dataset", "cm",
        "--use_gradscaler", "True",
        "--include_cell_idx", "True",
        "--pretrained_model", args.pretrained_model,
        "--num_epochs", str(args.epochs),
        "--batch_size", "32",
        "--init_epoch", "0",
        "--finetune_lr", "0.002",
        "--device", args.device,
    ]
    print(f"Running finetune:\n  {' '.join(cmd)}")
    subprocess.run(cmd, cwd=XPERT_ROOT, check=True)

    if not STAGING_DIR.exists() or not list(STAGING_DIR.glob("*.pth")):
        print(f"ERROR: training produced no weights in {STAGING_DIR}")
        sys.exit(1)
    OUT_DIR.parent.mkdir(parents=True, exist_ok=True)
    shutil.rmtree(OUT_DIR, ignore_errors=True)
    shutil.move(str(STAGING_DIR), str(OUT_DIR))
    shutil.rmtree(WORKFLOW_ROOT / "finetune", ignore_errors=True)

    print(f"Step 06 complete — finetuned weights in {OUT_DIR}")


if __name__ == "__main__":
    main()
