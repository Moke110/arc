"""
Step 4: Ensemble model prediction for all FDA drugs on cardiomyocytes.
Uses all 7 fold checkpoints from dtoxs_warm, averages predictions.
Run in xpert conda env.
"""
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import sys
import os
from tqdm import tqdm

# Add xpert src to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from src.model_xpert import XPertNet, get_unimol_drug_feat
from src.model_utils import Embeddings, cell_Embeddings, unimol_Embeddings, Encoder, crossEncoder
from src.config import load_config, resolve_paths

FDA_DIR = "vcms/fda"
EMBED_DIR = "drug_embeddings"
CHECKPOINT_DIR = "finetune/output/dtoxs_warm"
N_FOLDS = 7

# DTOX data settings
PERT_DOSE_IDX = 0
PERT_TIME_IDX = 2
CELL_IDX = 0
TISSUE_IDX = 0

class Args:
    """Simulate argparse namespace for model initialization."""
    drug_feat = "unimol"
    dataset = "dtoxs"
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
    use_gradscaler = False

def load_config_dtoxs():
    """Load and resolve dtoxs finetune config."""
    config_path = "finetune/config/dtoxs"
    config = load_config(config_path)
    return config

def load_drug_feat(config):
    """Load all drug feature tables."""
    unimol = np.load(config["drug_embeddings"]["unimol"], allow_pickle=True)
    hg = np.load(config["drug_embeddings"]["HG"], allow_pickle=True)
    return unimol, hg

def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    # Load config
    config = load_config_dtoxs()
    args = Args()

    # Load drug features (merged: 8981 old + 2123 new = 11104 total)
    print("Loading drug features...")
    # Use merged files (Step 3 updated these in-place)
    unimol_arr = np.load(EMBED_DIR + "/all_drugs_unimol_arr.npy", allow_pickle=True)
    hg_arr = np.load(EMBED_DIR + "/all_drugs_idx2HG.npy", allow_pickle=True)
    # SMILES merged file (created as all_drugs_idx2smi.npy by step3)
    smi_merged_path = EMBED_DIR + "/all_drugs_idx2smi.npy"
    smi_dict = np.load(smi_merged_path, allow_pickle=True).item()
    print(f"  UniMol: {unimol_arr.shape}")
    print(f"  HG: {hg_arr.shape}")
    print(f"  SMILES dict: {len(smi_dict)}")

    # Load FDA SMILES and map to global indices
    fda_df = pd.read_csv(f"{FDA_DIR}/fda_drugs.csv")
    fda_smiles = fda_df["canonical_smiles"].dropna().str.strip().tolist()
    print(f"\nFDA drugs to predict: {len(fda_smiles)}")

    # Build reverse lookup: SMILES -> global drug index
    smi_to_global_idx = {v: k for k, v in smi_dict.items()}
    fda_global_indices = []
    fda_valid_smiles = []
    missing = 0
    for smi in fda_smiles:
        if smi in smi_to_global_idx:
            fda_global_indices.append(smi_to_global_idx[smi])
            fda_valid_smiles.append(smi)
        else:
            missing += 1
    print(f"FDA drugs mapped to global indices: {len(fda_global_indices)}")
    print(f"Missing: {missing}")

    # Get control expression (mean across dtoxs training data)
    import h5py
    from scipy import sparse
    with h5py.File(config["dataset"]["DTOX_data_root"], "r") as f:
        g = f["obsm/X_ctl"]
        X_ctl = sparse.csr_matrix(
            (g["data"][:], g["indices"][:], g["indptr"][:]),
            shape=(int(g.attrs.get("shape", [521, 978])[0]), int(g.attrs.get("shape", [521, 978])[1]))
        ).toarray()
        gene_symbols = f["var/gene_symbol"][:]
    ctl_mean = X_ctl.mean(axis=0).astype(np.float32)
    print(f"Control expr mean shape: {ctl_mean.shape}")
    n_genes = len(ctl_mean)

    # Binning parameters from config
    n_bins = config["dataset"]["n_bins"]
    max_value = config["dataset"]["max_value"]
    min_value = config["dataset"]["min_value"]
    bins = np.quantile(np.array([min_value, max_value]), np.linspace(0, 1, n_bins - 1))
    ctl_binned = np.digitize(ctl_mean, bins).astype(np.int64)

    # Load gene names for output columns
    gene_names = [g.decode("utf-8") if isinstance(g, bytes) else str(g) for g in gene_symbols]

    # Collect predictions from all folds
    all_y_preds = []  # list of (n_fda_drugs, n_genes)

    for fold in range(1, N_FOLDS + 1):
        print(f"\n{'='*60}")
        print(f"Processing Fold {fold}/{N_FOLDS}")
        print(f"{'='*60}")

        # Load checkpoint
        ckpt_path = f"{CHECKPOINT_DIR}/dtoxs_warm_fold-{fold}_best.pth"
        print(f"Loading checkpoint: {ckpt_path}")
        ckpt = torch.load(ckpt_path, map_location=device, weights_only=False)

        # Build model with dummy logger
        import logging
        dummy_logger = logging.getLogger(f"fold_{fold}")
        dummy_logger.setLevel(logging.WARNING)
        if not dummy_logger.handlers:
            dummy_logger.addHandler(logging.StreamHandler(sys.stdout))
        model = XPertNet(args, config, device, logger=dummy_logger)

        # Load state dict
        model.load_state_dict(ckpt["model_state_dict"])
        model.to(device)
        model.eval()

        # Load embeddings into model
        model.drug_HG_embed = torch.tensor(hg_arr, dtype=torch.float, device=device)

        # Predict for each FDA drug in batches
        fold_preds = []
        batch_size = 64
        n_drugs = len(fda_global_indices)

        with torch.no_grad():
            for start in tqdm(range(0, n_drugs, batch_size), desc=f"Fold {fold}"):
                end = min(start + batch_size, n_drugs)
                batch_indices = fda_global_indices[start:end]
                batch_n = len(batch_indices)

                # Get drug features for this batch
                batch_unimol = torch.tensor(unimol_arr[batch_indices], dtype=torch.float, device=device)
                batch_drug_idx = torch.tensor(batch_indices, dtype=torch.long, device=device)

                # Prepare inputs
                trt_raw_data = torch.tensor(np.tile(ctl_mean, (batch_n, 1)), dtype=torch.float, device=device)
                ctl_raw_data = torch.tensor(np.tile(ctl_mean, (batch_n, 1)), dtype=torch.float, device=device)
                trt_raw_data_binned = torch.tensor(np.tile(ctl_binned, (batch_n, 1)), dtype=torch.long, device=device)
                ctl_raw_data_binned = torch.tensor(np.tile(ctl_binned, (batch_n, 1)), dtype=torch.long, device=device)
                pert_dose_idx = torch.full((batch_n,), PERT_DOSE_IDX, dtype=torch.long, device=device)
                pert_time_idx = torch.full((batch_n,), PERT_TIME_IDX, dtype=torch.long, device=device)
                cell_idx = torch.full((batch_n,), CELL_IDX, dtype=torch.long, device=device)
                tissue_idx = torch.full((batch_n,), TISSUE_IDX, dtype=torch.long, device=device)

                # Forward pass (infer mode returns 9 values)
                data = (trt_raw_data, ctl_raw_data, trt_raw_data_binned, ctl_raw_data_binned,
                        batch_unimol, pert_dose_idx, pert_time_idx, batch_drug_idx, cell_idx, tissue_idx)

                trt_output, ctl_output, output_deg, _, _, attention_dict, _, _, _ = model(data)

                # y_pred = output_deg + ctl_raw_data (predicted perturbed expression)
                y_pred = output_deg + ctl_raw_data
                fold_preds.append(y_pred.cpu().numpy())

        fold_y_pred = np.concatenate(fold_preds, axis=0)
        print(f"Fold {fold} predictions shape: {fold_y_pred.shape}")
        all_y_preds.append(fold_y_pred)

        # Free GPU memory
        del model
        torch.cuda.empty_cache()

    # Ensemble: average predictions across folds
    print(f"\n{'='*60}")
    print("Computing ensemble (mean) predictions...")
    all_y_preds_stacked = np.stack(all_y_preds, axis=0)  # (n_folds, n_drugs, n_genes)
    ensemble_preds = all_y_preds_stacked.mean(axis=0)  # (n_drugs, n_genes)
    print(f"Ensemble predictions shape: {ensemble_preds.shape}")

    # Save as CSV
    print("Saving results...")
    df_pred = pd.DataFrame(ensemble_preds, columns=gene_names)
    df_pred.insert(0, "canonical_smiles", fda_valid_smiles)
    os.makedirs(FDA_DIR, exist_ok=True)
    outpath = f"{FDA_DIR}/fda_cm_perturb.csv"
    df_pred.to_csv(outpath, index=False)
    print(f"Saved to {outpath}")
    print(f"Shape: {df_pred.shape[0]} drugs x {df_pred.shape[1]-1} genes + SMILES column")

    # Also save ensemble std for uncertainty estimation
    ensemble_std = all_y_preds_stacked.std(axis=0)
    df_std = pd.DataFrame(ensemble_std, columns=gene_names)
    df_std.insert(0, "canonical_smiles", fda_valid_smiles)
    std_path = f"{FDA_DIR}/fda_cm_perturb_std.csv"
    df_std.to_csv(std_path, index=False)
    print(f"Saved std to {std_path}")

    print("\nDone!")

if __name__ == "__main__":
    main()
