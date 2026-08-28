#!/usr/bin/env python3
"""
Step 01 — Drug embeddings (author reference, no HG training).

The finalized screen screens only L1000-contained FDA drugs (drug_idx < 8981),
so the HG biological embeddings are **the author reference table itself** — the
exact embedding the author's XPert was pretrained on. No HG training, no low-LR
finetune, no FDA-extension inference is ever run here.

  01_drug_embed/output/all_drugs_idx2HG.npy   (8981 x 256)
      == input/HG_drug_embeddings.npy  (byte-identical; author reference,
         workflow-owned in input/)

The UniMol 3D features stay the shared frozen array (input/all_drugs_unimol_arr.npy);
the model indexes it by drug_idx, and every screened drug has drug_idx < 8981.

The step is idempotent and parity-gated: if the output already exists AND is
byte-identical to the author reference it is reused; pass --force to rewrite.

Usage:
    python 01_drug_embed/drug_embeddings.py
    python 01_drug_embed/drug_embeddings.py --force
"""

import argparse
import shutil
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from utils import WORKFLOW_ROOT, AUTHOR_HG_REFERENCE, N_L1000_DRUGS  # noqa: E402

STEP_DIR = WORKFLOW_ROOT / "01_drug_embed"
INPUT_DIR = WORKFLOW_ROOT / "input"
OUTPUT_DIR = STEP_DIR / "output"
UNIMOL_ARR = "all_drugs_unimol_arr.npy"
SMI_MAPPING = "all_drugs_idx2smi.npy"
HG_EMBEDDINGS = "all_drugs_idx2HG.npy"


def _drug_count(arr) -> int:
    if getattr(arr, "ndim", 0) == 0 and arr.dtype == object:
        return len(arr.item())
    if hasattr(arr, "shape") and arr.ndim:
        return arr.shape[0]
    return len(arr)


def author_reference_intact() -> bool:
    return AUTHOR_HG_REFERENCE.exists()


def parity_with_author() -> bool:
    """Byte-identical parity gate against the author reference (rows 0-8980)."""
    out = OUTPUT_DIR / HG_EMBEDDINGS
    if not out.exists():
        return False
    a = np.load(AUTHOR_HG_REFERENCE)
    b = np.load(out, mmap_mode="r")
    if a.shape != b.shape:
        print(f"ERROR: shape {b.shape} != author reference {a.shape}")
        return False
    return np.array_equal(a, b)


def main():
    parser = argparse.ArgumentParser(description="Verify frozen inputs + author-reference HG")
    parser.add_argument("--force", action="store_true",
                        help="Rewrite the HG table even if an intact parity output exists")
    args = parser.parse_args()

    # 1. Static inputs (shared frozen artifacts, gitignored).
    for name in (UNIMOL_ARR, SMI_MAPPING):
        path = INPUT_DIR / name
        if not path.exists():
            print(f"ERROR: missing static input {path}")
            sys.exit(1)
        try:
            arr = np.load(path, mmap_mode="r", allow_pickle=True)
        except ValueError:
            arr = np.load(path, allow_pickle=True)
        print(f"Static input present: {path} ({_drug_count(arr)} drugs)")

    # 2. Author HG reference must exist (that is the entire screen HG table).
    if not author_reference_intact():
        print(f"ERROR: author HG reference not found: {AUTHOR_HG_REFERENCE}")
        sys.exit(1)
    ref = np.load(AUTHOR_HG_REFERENCE)
    if ref.shape != (N_L1000_DRUGS, 256):
        print(f"ERROR: author reference shape {ref.shape}, expected "
              f"({N_L1000_DRUGS}, 256)")
        sys.exit(1)

    # 3. Write / reuse the workflow copy with the parity gate.
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUTPUT_DIR / HG_EMBEDDINGS
    if not args.force and parity_with_author():
        print(f"Parity-intact author-reference HG reused: {out_path} "
              f"(byte-identical to the author reference)")
    else:
        shutil.copy2(AUTHOR_HG_REFERENCE, out_path)
        if not parity_with_author():
            print("ERROR: parity gate failed after copy")
            sys.exit(1)
        print(f"Wrote {out_path} (byte-identical to the author reference)")

    for path in (INPUT_DIR / UNIMOL_ARR, INPUT_DIR / SMI_MAPPING, out_path):
        arr = np.load(path, allow_pickle=True)
        print(f"  {path.name}: {_drug_count(arr)} drugs")
    print("Step 01 complete — no HG training performed (author-aligned screen).")


if __name__ == "__main__":
    main()
