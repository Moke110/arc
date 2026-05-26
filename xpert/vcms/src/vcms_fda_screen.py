"""
vcms_fda_screen.py — Ensemble virtual drug screening on cardiomyocytes.

Loads 5 fold models from finetune/output/dtoxs_warm-drug/, runs inference
for each FDA-approved drug against the unified control transcriptome X_cm,
averages deg predictions across folds, and saves results incrementally.

Usage:
  python vcms/vcms_fda_screen.py [--device cuda:0] [--sequential]
"""
import sys
import os
import time
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.model_xpert import XPertNet, get_unimol_drug_feat
from src.config import load_config

# ---------------------------------------------------------------------------
# Dummy logger — XPertNet uses logger.info() during __init__
# ---------------------------------------------------------------------------
class DummyLogger:
    def info(self, msg): pass
    def warning(self, msg): pass
    def error(self, msg): pass


# ---------------------------------------------------------------------------
# Minimal args namespace matching what XPertNet.__init__ expects
# ---------------------------------------------------------------------------
class ScreenArgs:
    dataset = "dtoxs"
    drug_feat = "unimol"
    pretrained_mode = "global"
    include_cell_idx = True
    wo_HG = False
    wo_atom = False
    wo_atom_HG = False
    wo_unimol = False
    wo_ppi = False
    use_gene_pos_emed = False
    mode = "infer"
    output_attention = False
    output_cls_embed = False


# ---------------------------------------------------------------------------
# Binning (mirrors src/dataset.py _digitize with side="one")
# ---------------------------------------------------------------------------
def bin_expression(x, min_val=1.854, max_val=15.993, n_bins=128):
    bins = np.quantile(np.array([min_val, max_val]), np.linspace(0, 1, n_bins - 1))
    return np.digitize(x, bins).astype(np.int64)


# ---------------------------------------------------------------------------
# Model loading
# ---------------------------------------------------------------------------
def load_model(fold_idx, config, device):
    model_path = (
        PROJECT_ROOT / "finetune" / "output" / "dtoxs_warm-drug"
        / f"dtoxs_warm-drug_fold-{fold_idx}_best.pth"
    )
    args = ScreenArgs()
    model = XPertNet(args, config, device, logger=DummyLogger())
    model.to(device)
    model.eval()

    ckpt = torch.load(str(model_path), map_location=device)
    model_dict = model.state_dict()
    pretrained_dict = {
        k: v
        for k, v in ckpt["model_state_dict"].items()
        if k in model_dict and v.shape == model_dict[k].shape
    }
    model_dict.update(pretrained_dict)
    model.load_state_dict(model_dict)

    print(f"  Loaded fold-{fold_idx} | best_val_loss={ckpt.get('best_val_loss', '?')}")
    return model


# ---------------------------------------------------------------------------
# Single-drug prediction across all models
# ---------------------------------------------------------------------------
def predict_one_drug(
    drug_idx, models, x_cm_tensor, x_cm_binned_tensor,
    drug_unimol_arr, device,
):
    """
    Run one drug through all loaded models, return average deg_output (978,).
    """
    batch_dim = 1

    trt = x_cm_tensor.unsqueeze(0)            # (1, 978)
    ctl = x_cm_tensor.unsqueeze(0)            # (1, 978)
    trt_b = x_cm_binned_tensor.unsqueeze(0)   # (1, 978)
    ctl_b = x_cm_binned_tensor.unsqueeze(0)   # (1, 978)

    drug_feat_np = drug_unimol_arr[drug_idx]  # (122, 514) float64
    drug_feat = torch.tensor(drug_feat_np, dtype=torch.float32).unsqueeze(0)

    pert_dose_idx = torch.tensor([0], dtype=torch.long)
    pert_time_idx = torch.tensor([2], dtype=torch.long)
    drug_idx_t = torch.tensor([drug_idx], dtype=torch.long)
    cell_idx = torch.tensor([0], dtype=torch.long)
    tissue_idx = torch.tensor([0], dtype=torch.long)

    fold_degs = []
    for model in models:
        data = (
            trt, ctl, trt_b, ctl_b, drug_feat,
            pert_dose_idx, pert_time_idx, drug_idx_t, cell_idx, tissue_idx,
        )
        with torch.no_grad():
            _, _, deg_output, _, _, _, _, _, _ = model(data, mode="infer")
        fold_degs.append(deg_output.squeeze(0).cpu().numpy())

    return np.mean(fold_degs, axis=0)  # (978,)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Cardiomyocyte FDA drug virtual screening")
    parser.add_argument("--device", default="cuda:0", help="Device (cuda:0, cpu)")
    parser.add_argument("--sequential", action="store_true",
                        help="Keep only 1 model on GPU at a time (low VRAM mode)")
    args_cli = parser.parse_args()

    device = torch.device(args_cli.device if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    # --- Load config ---
    config = load_config("finetune/config/dtox")
    print(f"Config loaded: {len(config)} top-level keys")

    # --- Load X_cm ---
    xcm_path = PROJECT_ROOT / "vcms" / "Xcm.csv"
    if not xcm_path.exists():
        print(f"ERROR: {xcm_path} not found. Run Xcm.py first.")
        sys.exit(1)
    xcm_df = pd.read_csv(xcm_path)
    gene_symbols = list(xcm_df.columns)
    x_cm = xcm_df.iloc[0].values.astype(np.float32)
    print(f"Loaded X_cm: {len(gene_symbols)} genes, range [{x_cm.min():.2f}, {x_cm.max():.2f}]")

    # --- Load FDA drug list ---
    fda_path = PROJECT_ROOT / "vcms" / "fda" / "fda_drugs.csv"
    if not fda_path.exists():
        print(f"ERROR: {fda_path} not found. Run fda_drugs.py first.")
        sys.exit(1)
    fda_df = pd.read_csv(fda_path)
    print(f"Loaded {len(fda_df)} FDA drugs from fda_drugs.csv")

    # --- Check existing output for resume ---
    output_path = PROJECT_ROOT / "vcms" / "fda" / "cm_fda_pert_deg.csv"
    completed_indices = set()
    if output_path.exists():
        existing = pd.read_csv(output_path, index_col=0)
        completed_indices = set(existing.index.tolist())
        print(f"Resume: {len(completed_indices)} drugs already completed")

    pending = fda_df[~fda_df["embedding_index"].isin(completed_indices)]
    if len(pending) == 0:
        print("All drugs already processed. Done.")
        return

    print(f"Drugs to process: {len(pending)}")

    # --- Pre-load UniMol embeddings ---
    print("Loading UniMol drug embeddings...")
    unimol_path = str(config["drug_embeddings"]["unimol"])
    drug_unimol_arr = np.load(unimol_path, allow_pickle=True)  # (8981, 122, 514), float64
    print(f"  UniMol shape: {drug_unimol_arr.shape}, dtype: {drug_unimol_arr.dtype}")

    # --- Bin X_cm ---
    x_cm_binned = bin_expression(x_cm)
    x_cm_tensor = torch.tensor(x_cm, dtype=torch.float32)
    x_cm_binned_tensor = torch.tensor(x_cm_binned, dtype=torch.long)

    # --- Load models ---
    print("Loading models...")
    if args_cli.sequential:
        models = []  # will load one at a time
    else:
        models = [load_model(i, config, device) for i in range(1, 6)]
    print(f"  {len(models)} models loaded")

    # --- Setup output file ---
    output_dir = output_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    if not output_path.exists():
        header = "embedding_index," + ",".join(gene_symbols)
        with open(output_path, "w") as f:
            f.write(header + "\n")

    # --- Screening loop ---
    total = len(pending)
    t_total_start = time.time()

    for i, (_, row) in enumerate(pending.iterrows()):
        t_start = time.time()
        drug_idx = int(row["embedding_index"])
        drug_name = str(row["drug_name"])

        if args_cli.sequential:
            # Load models one at a time, predict, then unload
            fold_degs = []
            for fold_i in range(1, 6):
                model = load_model(fold_i, config, device)
                deg = _predict_single_model(
                    model, drug_idx, x_cm_tensor, x_cm_binned_tensor,
                    drug_unimol_arr, device,
                )
                fold_degs.append(deg)
                del model
                torch.cuda.empty_cache()
            avg_deg = np.mean(fold_degs, axis=0)
        else:
            avg_deg = predict_one_drug(
                drug_idx, models, x_cm_tensor, x_cm_binned_tensor,
                drug_unimol_arr, device,
            )

        # Append to CSV
        row_str = f"{drug_idx}," + ",".join(f"{v:.6f}" for v in avg_deg)
        with open(output_path, "a") as f:
            f.write(row_str + "\n")
            f.flush()

        elapsed = time.time() - t_start
        eta = (time.time() - t_total_start) / (i + 1) * (total - i - 1)
        print(
            f"[{i + 1}/{total}] Drug {drug_idx} ({drug_name[:30]}): "
            f"deg mean={avg_deg.mean():+.4f} std={avg_deg.std():.4f} "
            f"| {elapsed:.1f}s | ETA {eta / 60:.0f}m{eta % 60:.0f}s"
        )

    t_total = time.time() - t_total_start
    print(f"\nDone. {total} drugs in {t_total / 60:.1f}m. Output: {output_path}")


def _predict_single_model(model, drug_idx, x_cm_tensor, x_cm_binned_tensor,
                          drug_unimol_arr, device):
    """Predict deg for one drug using a single model. Used in sequential mode."""
    trt = x_cm_tensor.unsqueeze(0)
    ctl = x_cm_tensor.unsqueeze(0)
    trt_b = x_cm_binned_tensor.unsqueeze(0)
    ctl_b = x_cm_binned_tensor.unsqueeze(0)

    drug_feat = torch.tensor(drug_unimol_arr[drug_idx], dtype=torch.float32).unsqueeze(0)
    pert_dose_idx = torch.tensor([0], dtype=torch.long)
    pert_time_idx = torch.tensor([2], dtype=torch.long)
    drug_idx_t = torch.tensor([drug_idx], dtype=torch.long)
    cell_idx = torch.tensor([0], dtype=torch.long)
    tissue_idx = torch.tensor([0], dtype=torch.long)

    data = (
        trt, ctl, trt_b, ctl_b, drug_feat,
        pert_dose_idx, pert_time_idx, drug_idx_t, cell_idx, tissue_idx,
    )
    with torch.no_grad():
        _, _, deg_output, _, _, _, _, _, _ = model(data, mode="infer")
    return deg_output.squeeze(0).cpu().numpy()


if __name__ == "__main__":
    main()
