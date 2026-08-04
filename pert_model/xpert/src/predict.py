"""
XPert drug perturbation prediction module.

Supports:
  - Single model or ensemble (directory of .pth files → averaged deg)
  - Fresh UniMol feature generation from SMILES via unimol_tools
  - HG embedding lookup from precomputed arrays
  - Batch inference with direct 10-tuple construction (no DataLoader needed)
  - prediction_root pattern: all paths relative to a task root directory
"""

import argparse
import logging
import os
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from tqdm import tqdm

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.config import load_config
from src.model_xpert import XPertNet

logger = logging.getLogger("predict")
logger.setLevel(logging.INFO)
if not logger.handlers:
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)


# ---------------------------------------------------------------------------
# Drug index resolution
# ---------------------------------------------------------------------------
def load_drug_smiles_mapping(smi_path: str) -> dict:
    """Load idx→SMILES dict, build reverse SMILES→idx mapping.

    Returns dict with keys: 'idx2smi' (dict), 'smi2idx' (dict)
    """
    idx2smi = np.load(smi_path, allow_pickle=True).item()
    smi2idx = {v: k for k, v in idx2smi.items()}
    logger.info(f"Loaded {len(idx2smi)} drug SMILES from {smi_path}")
    return {"idx2smi": idx2smi, "smi2idx": smi2idx}


# ---------------------------------------------------------------------------
# UniMol feature generation
# ---------------------------------------------------------------------------
def generate_unimol_features(smiles_list, max_atoms=122):
    """Generate UniMol (max_atoms, 514) features for a list of SMILES.

    Mimics src.utils.generate_unimol_feat but accepts a list of SMILES strings
    instead of a dict. Returns a dict {smiles: np.ndarray(122, 514)}.

    Each output row: [input_mask (122,1) | atom_symbols_id (122,1) | atom_feat (122,512)]
    """
    from unimol_tools import UniMolRepr

    atom2id = {
        'N': 2, 'Se': 3, 'I': 4, 'As': 5, 'Pt': 6, 'S': 7, 'Sn': 8,
        'Mg': 9, 'P': 10, 'Cl': 11, 'Au': 12, 'F': 13, 'Co': 14,
        'Br': 15, 'Si': 16, 'O': 17, 'H': 18, 'B': 19, 'C': 20,
        'Ca': 21, 'Hg': 22, 'Li': 23, "[UNK]": 24, "Na": 25,
        'hg_molecule': 0, 'molecule': 1,
    }

    unique_smiles = list(dict.fromkeys(smiles_list))  # dedup, preserve order
    clf = UniMolRepr(data_type='molecule', remove_hs=False)
    unimol_repr = clf.get_repr(unique_smiles, return_atomic_reprs=True)
    assert len(unimol_repr['cls_repr']) == len(unique_smiles), \
        f"Only {len(unimol_repr['cls_repr'])}/{len(unique_smiles)} molecules returned!"

    result = {}
    for idx in range(len(unique_smiles)):
        molecule_feat = unimol_repr['cls_repr'][idx]
        atom_feat = unimol_repr['atomic_reprs'][idx]
        atom_symbols = unimol_repr['atomic_symbol'][idx]

        # Prepend two virtual tokens (hg_molecule + molecule)
        atom_feat = np.vstack((np.array([molecule_feat, molecule_feat]), atom_feat))
        atom_symbols = ['hg_molecule', 'molecule'] + list(atom_symbols)
        atom_symbols_ids = [atom2id.get(a, 24) for a in atom_symbols]  # default [UNK]=24

        n_atoms = len(atom_feat)
        if n_atoms < max_atoms:
            pad_len = max_atoms - n_atoms
            atom_feat_pad = np.pad(atom_feat, ((0, pad_len), (0, 0)), constant_values=0)
            atom_symbols_ids_pad = atom_symbols_ids + [0] * pad_len
            input_mask = [1] * n_atoms + [0] * pad_len
        else:
            atom_feat_pad = atom_feat[:max_atoms]
            atom_symbols_ids_pad = atom_symbols_ids[:max_atoms]
            input_mask = [1] * max_atoms

        # Build (122, 514) = [mask | atom_id | atom_512dim_feat]
        feat = np.hstack((
            np.array(input_mask)[:, None],
            np.array(atom_symbols_ids_pad)[:, None],
            atom_feat_pad,
        )).astype(np.float32)

        result[unique_smiles[idx]] = feat

    return result


# ---------------------------------------------------------------------------
# Expression binning
# ---------------------------------------------------------------------------
def get_bin_edges(min_value, max_value, n_bins):
    """Compute bin edges for digitizing expression values.

    Mirrors MyDataset.__init__ binning: np.quantile([min, max], linspace(0,1,n_bins-1)).
    For a two-point distribution, this is equivalent to linspace(min, max, n_bins+1)[1:-1].
    """
    return np.linspace(min_value, max_value, n_bins + 1)[1:-1]


def bin_expression(expr, bin_edges):
    """Digitize a (N, 978) expression array using the given bin edges."""
    binned = np.digitize(expr, bin_edges).astype(np.int64)
    return binned


# ---------------------------------------------------------------------------
# Model loading
# ---------------------------------------------------------------------------
def load_model(checkpoint_path, args, config, device):
    """Load a single XPertNet model from a checkpoint."""
    ckpt = torch.load(checkpoint_path, map_location=device, weights_only=False)
    model = XPertNet(args, config, device, logger).to(device)
    model.load_state_dict(ckpt['model_state_dict'], strict=True)
    model.args.mode = "infer"
    model.eval()
    return model


def load_ensemble_models(model_path, args, config, device):
    """Load model(s) from a single .pth file or a directory of .pth files.

    Returns list of (model, name) tuples.
    """
    model_path = Path(model_path)
    if model_path.is_file():
        logger.info(f"Loading single model: {model_path}")
        model = load_model(str(model_path), args, config, device)
        return [(model, model_path.stem)]
    elif model_path.is_dir():
        pth_files = sorted(model_path.glob("*.pth"))
        if not pth_files:
            raise FileNotFoundError(f"No .pth files found in {model_path}")
        logger.info(f"Loading ensemble of {len(pth_files)} models from {model_path}")
        models = []
        for p in pth_files:
            model = load_model(str(p), args, config, device)
            models.append((model, p.stem))
        return models
    else:
        raise FileNotFoundError(f"Model path not found: {model_path}")


# ---------------------------------------------------------------------------
# Inference
# ---------------------------------------------------------------------------
def prepare_batch(ctl_expr, drug_indices, drug_feats, pert_dose_idx, pert_time_idx,
                  cell_idx, tissue_idx, bin_edges):
    """Build a 10-tuple batch for model.forward() from raw arrays.

    Args:
        ctl_expr: (N, 978) float32 — baseline expression repeated N times
        drug_indices: (N,) int64 — HG embedding lookup indices
        drug_feats: (N, 122, 514) float32 — UniMol features
        pert_dose_idx: (N,) int64 — dose index
        pert_time_idx: (N,) int64 — time index
        cell_idx: (N,) int64 — cell type index
        tissue_idx: (N,) int64 — tissue index
        bin_edges: (n_bins-1,) — expression bin edges

    Returns: 10-tuple of torch tensors
    """
    ctl_binned = bin_expression(ctl_expr, bin_edges)

    return (
        torch.tensor(ctl_expr, dtype=torch.float32),                    # trt_raw_data (unused, set to ctl)
        torch.tensor(ctl_expr, dtype=torch.float32),                    # ctl_raw_data
        torch.tensor(ctl_binned, dtype=torch.long),                     # trt_raw_data_binned
        torch.tensor(ctl_binned, dtype=torch.long),                     # ctl_raw_data_binned
        torch.tensor(drug_feats, dtype=torch.float32),                  # drug_feat
        torch.tensor(pert_dose_idx, dtype=torch.int64),                 # pert_dose_idx
        torch.tensor(pert_time_idx, dtype=torch.int64),                 # pert_time_idx
        torch.tensor(drug_indices, dtype=torch.long),                   # drug_idx
        torch.tensor(cell_idx, dtype=torch.long),                       # cell_idx
        torch.tensor(tissue_idx, dtype=torch.long),                     # tissue_idx
    )


def predict_deg(models, ctl_expr_1d, drug_items, pert_dose_idx, pert_time_idx,
                config, device, process_dir=None, batch_size=64):
    """Predict deg (perturbed − control) for a list of drugs.

    Args:
        models: list of (model, name) tuples from load_ensemble_models()
        ctl_expr_1d: (978,) float32 — single baseline expression vector
        drug_items: list of (drug_name, canonical_smiles, drug_idx) tuples
        pert_dose_idx: int — dose index for all drugs
        pert_time_idx: int — time index for all drugs
        config: loaded DTOXS config dict
        device: torch device
        process_dir: optional directory for caching UniMol features
        batch_size: batch size for inference

    Returns:
        deg_pred: (N_drugs, 978) float32 — ensemble-mean deg prediction
        valid_items: filtered list of drug_items that were successfully predicted
        gene_names: list of 978 gene column names
    """
    n_drugs = len(drug_items)
    drug_names = [d[0] for d in drug_items]
    drug_smiles = [d[1] for d in drug_items]
    drug_indices = np.array([d[2] for d in drug_items], dtype=np.int64)

    # Generate or load UniMol features
    unimol_cache = None
    if process_dir:
        os.makedirs(process_dir, exist_ok=True)
        unimol_cache = os.path.join(process_dir, "unimol_feats.npy")

    unimol_feat_dict = _get_unimol_features(drug_smiles, unimol_cache)

    # Build feature array, tracking which drugs succeeded
    valid_mask = np.zeros(n_drugs, dtype=bool)
    feat_list = []
    for i, smi in enumerate(drug_smiles):
        if smi in unimol_feat_dict:
            valid_mask[i] = True
            feat_list.append(unimol_feat_dict[smi])

    n_skipped = np.sum(~valid_mask)
    if n_skipped > 0:
        logger.warning(f"Skipping {n_skipped}/{n_drugs} drugs without UniMol features")
        skipped = [(drug_names[i], drug_smiles[i]) for i in np.where(~valid_mask)[0]]
        if process_dir:
            pd.DataFrame(skipped, columns=["drug_name", "smiles"]).to_csv(
                os.path.join(process_dir, "skipped_drugs.csv"), index=False)

    drug_feats = np.stack(feat_list, axis=0).astype(np.float32)
    valid_items = [drug_items[i] for i in np.where(valid_mask)[0]]
    valid_indices = drug_indices[valid_mask]
    valid_n = len(valid_items)

    if valid_n == 0:
        logger.error("No valid drugs to predict!")
        return None, [], []

    # Prepare inputs
    cell_idx = np.zeros(valid_n, dtype=np.int64)
    tissue_idx = np.zeros(valid_n, dtype=np.int64)
    dose_arr = np.full(valid_n, pert_dose_idx, dtype=np.int64)
    time_arr = np.full(valid_n, pert_time_idx, dtype=np.int64)

    bin_edges = get_bin_edges(
        config['dataset']['min_value'],
        config['dataset']['max_value'],
        config['dataset']['n_bins'],
    )

    # Run ensemble inference in batches
    ctl_1d = np.asarray(ctl_expr_1d, dtype=np.float32)
    all_degs = []  # list of (n_models, n_drugs, 978)

    for model, name in models:
        model.to(device)
        model.eval()
        deg_chunks = []

        for start in tqdm(range(0, valid_n, batch_size), desc=f"Predicting [{name}]", leave=False):
            end = min(start + batch_size, valid_n)
            batch_slice = slice(start, end)
            n_batch = end - start

            ctl_batch = np.tile(ctl_1d, (n_batch, 1))
            batch_data = prepare_batch(
                ctl_batch,
                valid_indices[batch_slice],
                drug_feats[batch_slice],
                dose_arr[batch_slice],
                time_arr[batch_slice],
                cell_idx[batch_slice],
                tissue_idx[batch_slice],
                bin_edges,
            )

            # Move to device and run
            batch_data = tuple(t.to(device) if isinstance(t, torch.Tensor) else t for t in batch_data)
            with torch.no_grad():
                _, _, deg_output, _, _, _, _, _, _ = model(batch_data)
            deg_chunks.append(deg_output.cpu().numpy())

        all_degs.append(np.concatenate(deg_chunks, axis=0))

    # Ensemble mean
    all_degs = np.stack(all_degs, axis=0)  # (n_models, n_drugs, 978)
    deg_mean = all_degs.mean(axis=0)

    # Gene column names from L1000 metadata
    l1000_meta = pd.read_csv(str(PROJECT_ROOT / "metadata" / "l1000_gene_info_978.csv"))
    gene_names = l1000_meta.sort_values("gene_idx")["gene_symbol"].tolist()

    return deg_mean, valid_items, gene_names


def _get_unimol_features(smiles_list, cache_path=None):
    """Get UniMol features, with optional disk cache."""
    if cache_path and os.path.exists(cache_path):
        logger.info(f"Loading cached UniMol features from {cache_path}")
        with open(cache_path, 'rb') as f:
            return pickle.load(f)

    logger.info(f"Generating UniMol features for {len(smiles_list)} SMILES...")
    feats = generate_unimol_features(smiles_list)

    if cache_path:
        with open(cache_path, 'wb') as f:
            pickle.dump(feats, f)
        logger.info(f"Cached UniMol features to {cache_path}")

    return feats


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="XPert drug perturbation prediction")
    p.add_argument("--prediction_root", required=True,
                   help="Root directory for this prediction task")
    p.add_argument("--ctl_csv", required=True,
                   help="Baseline transcriptome CSV path, relative to prediction_root")
    p.add_argument("--drug_csv", required=True,
                   help="Drug list CSV path (needs 'canonical_smiles' column), relative to prediction_root")
    p.add_argument("--output_csv", default="deg-perturb.csv",
                   help="Output deg CSV path, relative to prediction_root")
    p.add_argument("--dose_idx", type=int, default=0, help="Dose index (default 0)")
    p.add_argument("--time_idx", type=int, default=2, help="Time index: 0=3h, 1=6h, 2=24h (default 2)")
    p.add_argument("--model_path", required=True,
                   help="Single .pth file or directory of .pth files (absolute or relative to CWD)")
    p.add_argument("--device", default="cuda:0", help="Device (default cuda:0)")
    p.add_argument("--batch_size", type=int, default=64)
    p.add_argument("--config", default="finetune/config/dtoxs",
                   help="Config path relative to xpert/ (default finetune/config/dtoxs)")
    return p.parse_args()


def run_prediction(args):
    prediction_root = Path(args.prediction_root).resolve()
    process_dir = prediction_root / "process"
    ctl_path = prediction_root / args.ctl_csv
    drug_path = prediction_root / args.drug_csv
    output_path = prediction_root / args.output_csv

    logger.info(f"Prediction root: {prediction_root}")
    logger.info(f"Baseline: {ctl_path}")
    logger.info(f"Drug list: {drug_path}")
    logger.info(f"Output:   {output_path}")
    logger.info(f"Process:  {process_dir}")

    # Validate inputs
    for p, label in [(ctl_path, "Baseline CSV"), (drug_path, "Drug CSV")]:
        if not p.exists():
            raise FileNotFoundError(f"{label} not found: {p}")

    # Load config and drug mappings
    config = load_config(args.config)
    smi_path = config['drug_embeddings']['smi']
    if not os.path.isabs(smi_path):
        smi_path = str(PROJECT_ROOT / smi_path)
    mapping = load_drug_smiles_mapping(smi_path)

    # Load baseline expression
    ctl_df = pd.read_csv(ctl_path)
    if 'mean_log2_cpm1' in ctl_df.columns:
        ctl_expr = ctl_df['mean_log2_cpm1'].values.astype(np.float32)
    else:
        # Try to find a numeric expression column
        numeric_cols = ctl_df.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) >= 1:
            ctl_expr = ctl_df[numeric_cols[0]].values.astype(np.float32)
        else:
            raise ValueError(f"No numeric expression column found in {ctl_path}")

    if len(ctl_expr) != 978:
        raise ValueError(f"Baseline should have 978 genes, got {len(ctl_expr)}")

    logger.info(f"Loaded baseline: {len(ctl_expr)} genes, range [{ctl_expr.min():.3f}, {ctl_expr.max():.3f}]")

    # Load drug list
    drug_df = pd.read_csv(drug_path)
    if 'canonical_smiles' not in drug_df.columns:
        if 'smiles' in drug_df.columns:
            drug_df['canonical_smiles'] = drug_df['smiles']
        else:
            raise ValueError(f"Drug CSV must have 'canonical_smiles' or 'smiles' column")

    drug_name_col = 'pref_name' if 'pref_name' in drug_df.columns else (
        drug_df.columns[0] if drug_df.columns[0] != 'canonical_smiles' else drug_df.columns[1]
        if len(drug_df.columns) > 1 else None
    )

    # Resolve drug indices
    drug_items = []
    n_missing = 0
    for _, row in drug_df.iterrows():
        smi = row['canonical_smiles']
        name = row.get(drug_name_col, smi) if drug_name_col else smi
        idx = mapping['smi2idx'].get(smi, -1)
        if idx >= 0:
            drug_items.append((name, smi, idx))
        else:
            n_missing += 1

    logger.info(f"Drugs with valid embeddings: {len(drug_items)}/{len(drug_df)}")
    if n_missing:
        logger.warning(f"{n_missing} drugs not found in SMILES index — will be skipped")
        skipped = [(row.get(drug_name_col, row['canonical_smiles']), row['canonical_smiles'])
                   for _, row in drug_df.iterrows()
                   if mapping['smi2idx'].get(row['canonical_smiles'], -1) < 0]
        os.makedirs(process_dir, exist_ok=True)
        pd.DataFrame(skipped, columns=["drug_name", "smiles"]).to_csv(
            process_dir / "skipped_drugs.csv", index=False)

    # Build model args
    model_args = argparse.Namespace(
        dataset='dtoxs',
        drug_feat='unimol',
        drug_embeds='unimol,HG',
        mode='infer',
        pretrained_mode='global',
        include_cell_idx=True,
        output_attention=False,
        output_cls_embed=False,
        wo_HG=False,
        wo_atom=False,
        wo_atom_HG=False,
        wo_unimol=False,
        wo_ppi=False,
        use_gene_pos_emed=False,
        batch_size=args.batch_size,
    )

    device = torch.device(args.device if torch.cuda.is_available() else 'cpu')
    logger.info(f"Device: {device}")

    # Load models
    model_path = args.model_path
    if not os.path.isabs(model_path):
        model_path = str(Path(model_path).resolve())
    models = load_ensemble_models(model_path, model_args, config, device)
    logger.info(f"Loaded {len(models)} model(s)")

    # Run prediction
    deg_pred, valid_items, gene_names = predict_deg(
        models, ctl_expr, drug_items,
        pert_dose_idx=args.dose_idx,
        pert_time_idx=args.time_idx,
        config=config,
        device=device,
        process_dir=str(process_dir),
        batch_size=args.batch_size,
    )

    if deg_pred is None:
        logger.error("Prediction failed — no valid drugs")
        return

    # Save output
    logger.info(f"Saving predictions to {output_path}")
    output_dir = output_path.parent
    os.makedirs(output_dir, exist_ok=True)

    result = pd.DataFrame(deg_pred, columns=gene_names)
    result.insert(0, 'smiles', [item[1] for item in valid_items])
    result.insert(0, 'drug_name', [item[0] for item in valid_items])
    result.insert(2, 'drug_idx', [item[2] for item in valid_items])
    result.to_csv(output_path, index=False)

    n_genes = len(gene_names)
    logger.info(f"Done. Output: {len(result)} drugs × {n_genes} genes")
    logger.info(f"  Columns: drug_name, smiles, drug_idx, [{n_genes} gene symbols]")
    logger.info(f"  File: {output_path}")


if __name__ == "__main__":
    args = parse_args()
    run_prediction(args)
