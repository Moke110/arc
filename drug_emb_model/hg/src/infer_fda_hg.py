"""
Run HG inference on the extended heterogeneous graph to produce proper HG embeddings
for FDA-extension drugs (indices 8981-11103).

Standard HG process:
  1. Generate 512-dim input features for FDA drugs via Ridge(KPGT→512)
  2. Load existing graph (gene features, PPI, training drug features, DTI, DDS)
  3. Add FDA drug nodes with their features
  4. Add FDA DTI edges (drug→gene, gene→drug)
  5. Add FDA DDS edges (drug↔drug, between FDA-FDA and FDA-train)
  6. Full-graph forward pass → extract FDA drug embeddings

Usage:
  # After DTI data is collected:
  python -m src.infer_fda_hg --dti-fda data/DTI_edge_index_fda_ext.npy

  # Test with a subset:
  python -m src.infer_fda_hg --dti-fda data/DTI_edge_index_fda_ext.npy --test-embeddings
"""

import logging
import os
import sys
import time
from pathlib import Path
from collections import OrderedDict

import numpy as np
import torch
from torch_geometric.data import HeteroData
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

HG_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HG_ROOT))

from src.model_hg import HeteroGNN
from src.utils_hg import get_embeddings

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Step 1: Generate 512-dim input features for FDA drugs
# ---------------------------------------------------------------------------
class KPGT_to_HG512:
    """Learn KPGT(2304) → HG_input(512) mapping via Ridge regression.

    Fitted on the 8,981 training drugs (which have both KPGT and HG 512-dim features),
    then applied to FDA drugs to produce consistent 512-dim input features.
    """

    def __init__(self, alpha: float = 1000.0):
        self.alpha = alpha
        self.scaler = StandardScaler()
        self.ridge = Ridge(alpha=alpha, fit_intercept=True)

    def fit(self, kpgt_train: np.ndarray, hg_feat_512: np.ndarray):
        """Fit ridge regression on training drugs.

        Args:
            kpgt_train: (8981, 2304) KPGT features
            hg_feat_512: (8981, 512) original HG input features
        """
        kpgt_scaled = self.scaler.fit_transform(kpgt_train)
        self.ridge.fit(kpgt_scaled, hg_feat_512)

        # Log fit quality
        pred = self.ridge.predict(kpgt_scaled)
        ss_res = np.sum((hg_feat_512 - pred) ** 2)
        ss_tot = np.sum((hg_feat_512 - hg_feat_512.mean()) ** 2)
        r2 = 1 - ss_res / (ss_tot + 1e-10)
        logger.info(f"KPGT→512 Ridge fit: alpha={self.alpha}, R²={r2:.4f}")

    def transform(self, kpgt_fda: np.ndarray) -> np.ndarray:
        """Generate 512-dim features for FDA drugs.

        Args:
            kpgt_fda: (N_fda, 2304)

        Returns:
            (N_fda, 512) features
        """
        kpgt_scaled = self.scaler.transform(kpgt_fda)
        return self.ridge.predict(kpgt_scaled)


def generate_fda_node_features(
    kpgt_path: str,
    hg_train_feat_path: str,
    n_train: int = 8981,
    n_total: int = 11104,
) -> np.ndarray:
    """Generate 512-dim HG input features for FDA drugs.

    Args:
        kpgt_path: path to all_drugs_idx2KPGT.npy
        hg_train_feat_path: path to all_drug_node_feat_8981.npy
        n_train: number of training drugs (8981)
        n_total: total number of drugs (11104)

    Returns:
        (n_total - n_train, 512) FDA drug features
    """
    kpgt_all = np.load(kpgt_path, allow_pickle=True).item()
    hg_train_512 = np.load(hg_train_feat_path)

    kpgt_train = np.stack([kpgt_all[i] for i in range(n_train)])
    kpgt_fda = np.stack([kpgt_all[i] for i in range(n_train, n_total)])

    projector = KPGT_to_HG512()
    projector.fit(kpgt_train, hg_train_512)
    fda_feat = projector.transform(kpgt_fda)

    logger.info(f"Generated FDA features: {fda_feat.shape}, "
                f"range [{fda_feat.min():.3f}, {fda_feat.max():.3f}]")
    return fda_feat


# ---------------------------------------------------------------------------
# Step 2-5: Build extended graph and run inference
# ---------------------------------------------------------------------------
def build_extended_graph(
    data_dir: str,
    fda_node_feat: np.ndarray,
    dti_fda_path: str | None = None,
    dds_fda_path: str | None = None,
    include_edges: list[str] | None = None,
) -> HeteroData:
    """Build the extended HeteroData with FDA drugs added.

    Args:
        data_dir: path to hg/data/
        fda_node_feat: (N_fda, 512) — generated FDA drug features
        dti_fda_path: path to DTI_edge_index_fda_ext.npy (FDA DTI edges with LOCAL indices)
        dds_fda_path: path to DDS_edge_index_fda_ext.npy (FDA DDS edges with GLOBAL indices)
        include_edges: edge types to include (default: PPI, DTI, DDS)

    Returns:
        HeteroData with combined training + FDA graph
    """
    if include_edges is None:
        include_edges = ["PPI", "DTI", "DDS"]

    data_dir = Path(data_dir)
    n_train = 8981
    n_fda = fda_node_feat.shape[0]
    n_total = n_train + n_fda
    logger.info(f"Building extended graph: {n_train} train + {n_fda} FDA = {n_total} drugs")

    # Load existing node features
    gene_feat = torch.tensor(
        np.load(str(data_dir / "all_gene_node_feat_19392.npy")), dtype=torch.float)
    drug_feat_train = torch.tensor(
        np.load(str(data_dir / "all_drug_node_feat_8981.npy")), dtype=torch.float)
    fda_feat_tensor = torch.tensor(fda_node_feat, dtype=torch.float)

    # Combine drug features
    combined_drug_feat = torch.cat([drug_feat_train, fda_feat_tensor], dim=0)
    logger.info(f"  Drug features: {combined_drug_feat.shape}")
    logger.info(f"  Gene features: {gene_feat.shape}")

    # Build HeteroData
    data = HeteroData()
    data["gene"].x = gene_feat
    data["drug"].x = combined_drug_feat

    # --- Load existing edges ---
    edge_types = ["PPI", "GSEA", "PCC", "DTI", "DDS"]
    edge_files = [
        "PPI_edge_index_all_901260_pairs_weighted.npy",
        "GSEA_edge_indexes_all_pairs_426904_weighted.npy",
        "PCC_edge_Thr04_18674_pairs.npy",
        "DTI_edge_index_all_12890_pairs.npy",
        "DDS_edge_index_all_287834_pairs_weighted.npy",
    ]

    for edge_type, edge_file in zip(edge_types, edge_files):
        if edge_type not in include_edges:
            continue

        edge_path = data_dir / edge_file
        if not edge_path.exists():
            logger.warning(f"  Edge file not found: {edge_path}")
            continue

        edge_data = np.load(str(edge_path), allow_pickle=True).item()
        edge_index = torch.tensor(edge_data["edge_index"], dtype=torch.long)
        edge_attr = (
            torch.tensor(edge_data["edge_attr"], dtype=torch.float)
            if "weighted" in edge_file else None
        )

        if edge_type in ["PPI", "GSEA"]:
            data["gene", edge_type, "gene"].edge_index = edge_index
            data["gene", edge_type, "gene"].edge_attr = edge_attr
        elif edge_type == "PCC":
            data["gene", edge_type, "gene"].edge_index = edge_index
        elif edge_type == "DTI":
            data["drug", edge_type, "gene"].edge_index = edge_index
            data["gene", edge_type, "drug"].edge_index = torch.stack(
                [edge_index[1], edge_index[0]], dim=0)
        elif edge_type == "DDS":
            data["drug", edge_type, "drug"].edge_index = edge_index
            data["drug", edge_type, "drug"].edge_attr = edge_attr

    n_existing_dti = data["drug", "DTI", "gene"].edge_index.shape[1]
    logger.info(f"  Existing DTI edges: {n_existing_dti}")

    # --- Add FDA DTI edges ---
    if dti_fda_path and Path(dti_fda_path).exists():
        fda_dti_data = np.load(dti_fda_path, allow_pickle=True).item()
        fda_dti_edge_index = fda_dti_data["edge_index"]  # (2, N_edges) with LOCAL FDA indices

        if fda_dti_edge_index.shape[1] > 0:
            # Convert local FDA indices (0..N_fda-1) to global indices (n_train..n_total-1)
            global_edge_index = fda_dti_edge_index.copy()
            global_edge_index[0] = global_edge_index[0] + n_train  # offset FDA drug indices

            global_edge_index = torch.tensor(global_edge_index, dtype=torch.long)

            # Append to existing DTI edges
            existing_dti = data["drug", "DTI", "gene"].edge_index
            data["drug", "DTI", "gene"].edge_index = torch.cat(
                [existing_dti, global_edge_index], dim=1)

            # Add reverse edges too
            existing_dti_rev = data["gene", "DTI", "drug"].edge_index
            rev_edges = torch.stack([global_edge_index[1], global_edge_index[0]], dim=0)
            data["gene", "DTI", "drug"].edge_index = torch.cat(
                [existing_dti_rev, rev_edges], dim=1)

            logger.info(f"  Added {fda_dti_edge_index.shape[1]} FDA DTI edges")
            logger.info(f"  Total DTI edges: {data['drug', 'DTI', 'gene'].edge_index.shape[1]}")

    # --- Add FDA DDS edges ---
    if dds_fda_path and Path(dds_fda_path).exists():
        fda_dds_data = np.load(dds_fda_path, allow_pickle=True).item()
        fda_dds_edge_index = fda_dds_data["edge_index"]  # should be in GLOBAL indices

        if fda_dds_edge_index.shape[1] > 0:
            fda_dds_edge_index = torch.tensor(fda_dds_edge_index, dtype=torch.long)

            existing_dds = data["drug", "DDS", "drug"].edge_index
            data["drug", "DDS", "drug"].edge_index = torch.cat(
                [existing_dds, fda_dds_edge_index], dim=1)

            if "edge_attr" in fda_dds_data:
                fda_dds_attr = torch.tensor(fda_dds_data["edge_attr"], dtype=torch.float)
                existing_dds_attr = data["drug", "DDS", "drug"].edge_attr
                data["drug", "DDS", "drug"].edge_attr = torch.cat(
                    [existing_dds_attr, fda_dds_attr])

            logger.info(f"  Added {fda_dds_edge_index.shape[1]} FDA DDS edges")
            logger.info(f"  Total DDS edges: {data['drug', 'DDS', 'drug'].edge_index.shape[1]}")

    # --- Summary ---
    n_drugs = data["drug"].x.shape[0]
    n_genes = data["gene"].x.shape[0]
    n_dti = data["drug", "DTI", "gene"].edge_index.shape[1] if ("drug", "DTI", "gene") in data.edge_index_dict else 0
    n_dds = data["drug", "DDS", "drug"].edge_index.shape[1] if ("drug", "DDS", "drug") in data.edge_index_dict else 0
    n_ppi = data["gene", "PPI", "gene"].edge_index.shape[1] if ("gene", "PPI", "gene") in data.edge_index_dict else 0

    logger.info(f"Graph summary: {n_drugs} drugs, {n_genes} genes, "
                f"{n_dti} DTI edges, {n_dds} DDS edges, {n_ppi} PPI edges")

    return data


def run_inference(
    data: HeteroData,
    model_path: str,
    device: str = "cuda:0",
    hidden_dim: int = 256,
    out_dim: int = 256,
    num_layers: int = 3,
    include_edges: list[str] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Run HG forward pass and extract drug & gene embeddings.

    Args:
        data: HeteroData (extended graph)
        model_path: path to trained model checkpoint
        device: torch device
        hidden_dim, out_dim, num_layers: model architecture params

    Returns:
        (drug_embeddings, gene_embeddings)
        drug_embeddings: (N_total_drugs, out_dim)
        gene_embeddings: (N_genes, out_dim)
    """
    if include_edges is None:
        include_edges = ["PPI", "DTI", "DDS"]

    device = torch.device(device if torch.cuda.is_available() else "cpu")
    logger.info(f"Device: {device}")

    # Load model
    model = HeteroGNN(
        hidden_dim=hidden_dim,
        out_dim=out_dim,
        num_layers=num_layers,
        include_edges=include_edges,
    ).to(device)

    checkpoint = torch.load(model_path, map_location=device, weights_only=False)
    model.load_state_dict(checkpoint["model_state_dict"])
    logger.info(f"Loaded model from {model_path} (epoch {checkpoint.get('epoch', '?')})")

    model.eval()

    # Move data to device
    data = data.to(device)

    # Forward pass
    logger.info("Running HG forward pass on extended graph...")
    t0 = time.time()
    with torch.no_grad():
        drug_embeddings, gene_embeddings = get_embeddings(model, data)

    elapsed = time.time() - t0
    logger.info(f"Forward pass done in {elapsed:.1f}s")
    logger.info(f"  Drug embeddings: {drug_embeddings.shape}")
    logger.info(f"  Gene embeddings: {gene_embeddings.shape}")

    return drug_embeddings.cpu().numpy(), gene_embeddings.cpu().numpy()


# ---------------------------------------------------------------------------
# Analysis helpers
# ---------------------------------------------------------------------------
def analyze_embeddings(
    fda_embeddings_new: np.ndarray,
    fda_embeddings_old: np.ndarray,
    drug_indices: np.ndarray,
    gene_embeddings: np.ndarray,
    dti_data: dict | None = None,
):
    """Compare new (graph-based) vs old (ridge-only) FDA drug embeddings.

    Args:
        fda_embeddings_new: (N_fda, 256) from extended graph inference
        fda_embeddings_old: (N_fda, 256) from ridge-only approach
        drug_indices: (N_fda,) global HG drug indices
        gene_embeddings: (19392, 256) gene embeddings
        dti_data: optional DTI edge data for FDA drugs
    """
    n_fda = fda_embeddings_new.shape[0]

    logger.info("\n" + "=" * 60)
    logger.info("Embedding Analysis")
    logger.info("=" * 60)

    # 1. Per-drug cosine similarity (old vs new)
    cos_sims = []
    for i in range(n_fda):
        cos_sim = np.dot(fda_embeddings_new[i], fda_embeddings_old[i]) / (
            np.linalg.norm(fda_embeddings_new[i]) * np.linalg.norm(fda_embeddings_old[i]) + 1e-8)
        cos_sims.append(cos_sim)
    cos_sims = np.array(cos_sims)

    logger.info(f"Cosine similarity (new vs old FDA embeddings):")
    logger.info(f"  Mean: {cos_sims.mean():.4f}  Std: {cos_sims.std():.4f}")
    logger.info(f"  Min: {cos_sims.min():.4f}  Max: {cos_sims.max():.4f}")
    logger.info(f"  Fraction > 0.9: {(cos_sims > 0.9).mean():.2%}")
    logger.info(f"  Fraction < 0.5: {(cos_sims < 0.5).mean():.2%}")

    # 2. Distribution comparison
    logger.info(f"\nNew embeddings: mean={fda_embeddings_new.mean():.4f}, "
                f"std={fda_embeddings_new.std():.4f}")
    logger.info(f"Old embeddings: mean={fda_embeddings_old.mean():.4f}, "
                f"std={fda_embeddings_old.std():.4f}")

    # 3. Drug-gene cosine similarity (for drugs with DTI edges)
    if dti_data is not None:
        dti_edge_index = dti_data["edge_index"]  # (2, N_edges) local FDA drug→gene
        if dti_edge_index.shape[1] > 0:
            drug_local = dti_edge_index[0]  # local FDA index
            gene_hg = dti_edge_index[1]  # HG gene index

            dt_drug_emb = fda_embeddings_new[drug_local]
            dt_gene_emb = gene_embeddings[gene_hg]

            dti_cos = np.sum(dt_drug_emb * dt_gene_emb, axis=1) / (
                np.linalg.norm(dt_drug_emb, axis=1) *
                np.linalg.norm(dt_gene_emb, axis=1) + 1e-8)

            logger.info(f"\nDTI drug-gene cosine similarity:")
            logger.info(f"  N={len(dti_cos)}, mean={dti_cos.mean():.4f}, "
                        f"std={dti_cos.std():.4f}")

            # Compare with random drug-gene pairs
            n_random = min(10000, len(dti_cos) * 3)
            rand_drugs = np.random.randint(0, n_fda, n_random)
            rand_genes = np.random.randint(0, gene_embeddings.shape[0], n_random)
            rand_cos = np.sum(
                fda_embeddings_new[rand_drugs] * gene_embeddings[rand_genes], axis=1) / (
                np.linalg.norm(fda_embeddings_new[rand_drugs], axis=1) *
                np.linalg.norm(gene_embeddings[rand_genes], axis=1) + 1e-8)

            logger.info(f"Random drug-gene cosine similarity (n={n_random}):")
            logger.info(f"  mean={rand_cos.mean():.4f}, std={rand_cos.std():.4f}")
            logger.info(f"DTI vs Random ratio: {dti_cos.mean() / (rand_cos.mean() + 1e-8):.2f}x")

    return cos_sims


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    import argparse
    ap = argparse.ArgumentParser(description="HG inference for FDA drugs on extended graph")
    ap.add_argument("--data-dir", default=str(HG_ROOT / "data"))
    ap.add_argument("--kpgt-path", default=str(HG_ROOT.parent / "xpert/drug_embeddings/all_drugs_idx2KPGT.npy"))
    ap.add_argument("--model-path", default=str(HG_ROOT / "model/20260518_102121/global_fold_early_stop.pth"))
    ap.add_argument("--dti-fda", default=None, help="Path to FDA DTI edge .npy")
    ap.add_argument("--dds-fda", default=None, help="Path to FDA DDS edge .npy")
    ap.add_argument("--output", default=str(HG_ROOT / "model/20260518_102121/fda_hg_embeddings.npy"))
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--include-edges", nargs="+", default=["PPI", "DTI", "DDS"])
    ap.add_argument("--hidden-dim", type=int, default=256)
    ap.add_argument("--out-dim", type=int, default=256)
    ap.add_argument("--num-layers", type=int, default=3)
    ap.add_argument("--test-embeddings", action="store_true",
                    help="Analyze and compare new vs old embeddings")

    args = ap.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    # Step 1: Generate FDA drug node features
    logger.info("=" * 60)
    logger.info("Step 1: Generate 512-dim FDA drug features")
    logger.info("=" * 60)
    fda_feat_512 = generate_fda_node_features(
        args.kpgt_path,
        str(Path(args.data_dir) / "all_drug_node_feat_8981.npy"),
    )

    # Optionally save FDA features
    feat_path = Path(args.output).with_name("fda_node_feat_512.npy")
    np.save(str(feat_path), fda_feat_512)
    logger.info(f"Saved FDA node features to {feat_path}")

    # Step 2-4: Build extended graph
    logger.info("\n" + "=" * 60)
    logger.info("Step 2: Build extended heterogeneous graph")
    logger.info("=" * 60)
    data = build_extended_graph(
        args.data_dir,
        fda_feat_512,
        dti_fda_path=args.dti_fda,
        dds_fda_path=args.dds_fda,
        include_edges=args.include_edges,
    )

    # Step 5: Run inference
    logger.info("\n" + "=" * 60)
    logger.info("Step 3: Run HG inference")
    logger.info("=" * 60)
    drug_emb, gene_emb = run_inference(
        data,
        args.model_path,
        device=args.device,
        hidden_dim=args.hidden_dim,
        out_dim=args.out_dim,
        num_layers=args.num_layers,
        include_edges=args.include_edges,
    )

    # Extract FDA drug embeddings (last N_fda rows)
    n_fda = fda_feat_512.shape[0]
    fda_emb_new = drug_emb[-n_fda:]  # (N_fda, 256)

    # Save
    np.save(args.output, fda_emb_new)
    logger.info(f"\nSaved FDA HG embeddings to {args.output}")
    logger.info(f"  Shape: {fda_emb_new.shape}")

    # Optional: analyze
    if args.test_embeddings:
        old_hg_path = str(HG_ROOT / "embeddings/all_drugs_idx2HG.npy")
        old_hg = np.load(old_hg_path, allow_pickle=True)
        fda_emb_old = old_hg[-n_fda:]

        dti_data = None
        if args.dti_fda and Path(args.dti_fda).exists():
            dti_data = np.load(args.dti_fda, allow_pickle=True).item()

        analyze_embeddings(
            fda_emb_new, fda_emb_old,
            np.arange(8981, 8981 + n_fda),
            gene_emb,
            dti_data,
        )


if __name__ == "__main__":
    main()
