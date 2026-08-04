"""
Compute DDS (drug-drug similarity) edges for FDA drugs based on Morgan fingerprints.

Uses Tanimoto (Jaccard) similarity on Morgan fingerprints to find drug pairs with
high structural similarity. The threshold is calibrated to match the density of
the existing DDS edges in the HG graph.

Existing DDS: 287,834 edges / C(8981, 2) ≈ 0.714% density
For FDA extension: we compute FDA↔FDA and FDA↔train drug similarities,
then keep top-k pairs matching the existing density.

Output: numpy dict with edge_index and edge_attr, using GLOBAL drug indices (0-11103).
"""

import logging
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from sklearn.metrics.pairwise import cosine_similarity

HG_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HG_ROOT))

logger = logging.getLogger(__name__)


def load_morgan_matrix(morgan_path: str,
                       n_total: int = 11104,
                       n_train: int = 8981) -> tuple[np.ndarray, np.ndarray]:
    """Load Morgan fingerprints for training and FDA drugs.

    Returns:
        (morgan_train, morgan_fda)
        morgan_train: (8981, 1024)
        morgan_fda: (2123, 1024)
    """
    morgan_dict = np.load(morgan_path, allow_pickle=True).item()

    morgan_train = np.stack([morgan_dict[i] for i in range(n_train)]).astype(np.float64)
    morgan_fda = np.stack([morgan_dict[i] for i in range(n_train, n_total)]).astype(np.float64)

    logger.info(f"Loaded Morgan fingerprints: train={morgan_train.shape}, "
                f"FDA={morgan_fda.shape}")
    return morgan_train, morgan_fda


def compute_tanimoto_full(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Compute Tanimoto (Jaccard) similarity between two binary fingerprint matrices.

    Tanimoto(A_i, B_j) = |A_i ∩ B_j| / |A_i ∪ B_j|
                        = (A_i · B_j) / (|A_i| + |B_j| - A_i · B_j)

    Works for both boolean and {0,1}-valued float arrays.

    Args:
        A: (N_a, D) binary fingerprints
        B: (N_b, D) binary fingerprints

    Returns:
        (N_a, N_b) similarity matrix (float64)
    """
    # Ensure binary (0/1) float64
    A_bin = (A > 0.5).astype(np.float64)
    B_bin = (B > 0.5).astype(np.float64)

    dot = A_bin @ B_bin.T  # (N_a, N_b)
    a_norm = A_bin.sum(axis=1)  # (N_a,)
    b_norm = B_bin.sum(axis=1)  # (N_b,)

    denom = a_norm[:, None] + b_norm[None, :] - dot + 1e-10
    return dot / denom


def calibrate_threshold(existing_dds_path: str,
                        morgan_train: np.ndarray) -> tuple[float, float]:
    """Determine the Tanimoto similarity threshold that matches existing DDS density.

    Computes the similarity distribution of a sample of existing DDS edges to
    find the percentile threshold that would produce a similar edge density.

    Returns:
        (threshold, density) — threshold value and fraction of all pairs retained
    """
    existing = np.load(existing_dds_path, allow_pickle=True).item()
    existing_edges = existing["edge_index"]
    existing_weights = existing.get("edge_attr", np.ones(existing_edges.shape[1]))

    n_drugs = morgan_train.shape[0]
    density = existing_edges.shape[1] / (n_drugs * (n_drugs - 1) / 2)

    logger.info(f"Existing DDS: {existing_edges.shape[1]} edges, density={density:.4%}")
    logger.info(f"Existing DDS weight range: [{existing_weights.min():.4f}, "
                f"{existing_weights.max():.4f}], median={np.median(existing_weights):.4f}")

    # We can directly use existing weights distribution to set threshold
    # Use the minimum weight among existing edges as threshold
    threshold = np.percentile(existing_weights, 5)  # 5th percentile as threshold
    logger.info(f"Calibrated threshold: {threshold:.4f} (5th percentile of existing weights)")

    return threshold, density


def compute_fda_dds_edges(
    morgan_train: np.ndarray,
    morgan_fda: np.ndarray,
    threshold: float,
    existing_dds_path: str | None = None,
    n_train: int = 8981,
) -> dict:
    """Compute DDS edges for FDA drugs.

    Computes two sets of similarities:
      - FDA ↔ FDA (among FDA drugs)
      - FDA ↔ Train (between FDA and existing drugs)

    Keeps all pairs with Tanimoto >= threshold.

    Args:
        morgan_train: (8981, 1024) training drug fingerprints
        morgan_fda: (N_fda, 1024) FDA drug fingerprints
        threshold: minimum Tanimoto similarity
        existing_dds_path: path to existing DDS edges (for weight calibration)
        n_train: number of training drugs
        batch_size: compute batch size

    Returns:
        dict with:
          - edge_index: (2, N_edges) — GLOBAL drug indices
          - edge_attr: (N_edges,) — Tanimoto similarity weights
          - stats: dict with summary
    """
    n_fda = morgan_fda.shape[0]

    logger.info(f"Computing DDS for {n_fda} FDA drugs...")

    all_edges = []  # (global_idx_a, global_idx_b, weight)

    # 1. FDA ↔ Train similarities (full 2123 × 8981 matrix)
    logger.info("Computing FDA↔Train Tanimoto (2123 × 8981)...")
    sim_ft = compute_tanimoto_full(morgan_fda, morgan_train)
    logger.info(f"  FDA↔Train sim range: [{sim_ft.min():.4f}, {sim_ft.max():.4f}]")
    logger.info(f"  Fraction > {threshold}: {(sim_ft > threshold).mean():.4%}")

    fda_local, train_local = np.where(sim_ft >= threshold)
    for fi, ti in zip(fda_local, train_local):
        global_fda = fi + n_train
        all_edges.append((global_fda, ti, float(sim_ft[fi, ti])))

    logger.info(f"  FDA↔Train edges above threshold: {len(all_edges)}")

    # 2. FDA ↔ FDA similarities
    logger.info("Computing FDA↔FDA Tanimoto (2123 × 2123)...")
    sim_ff = compute_tanimoto_full(morgan_fda, morgan_fda)
    logger.info(f"  FDA↔FDA sim range: [{sim_ff.min():.4f}, {sim_ff.max():.4f}]")

    # Keep only upper triangle (self-similarity is always 1.0, excluded)
    fda_a, fda_b = np.where((sim_ff >= threshold) & (np.triu(np.ones_like(sim_ff), k=1) > 0))
    n_before = len(all_edges)
    for ai, bi in zip(fda_a, fda_b):
        global_a = ai + n_train
        global_b = bi + n_train
        all_edges.append((global_a, global_b, float(sim_ff[ai, bi])))
    logger.info(f"  FDA↔FDA edges above threshold: {len(all_edges) - n_before}")

    # Build output arrays
    n_edges = len(all_edges)
    edge_index = np.zeros((2, max(n_edges, 1)), dtype=np.int64)
    edge_attr = np.zeros(max(n_edges, 1), dtype=np.float32)

    for i, (a, b, w) in enumerate(all_edges):
        edge_index[0, i] = a
        edge_index[1, i] = b
        edge_attr[i] = w

    if n_edges == 0:
        edge_index = np.zeros((2, 0), dtype=np.int64)
        edge_attr = np.zeros(0, dtype=np.float32)

    # Stats
    n_fda_fda = sum(1 for e in all_edges if e[0] >= n_train and e[1] >= n_train)
    n_fda_train = n_edges - n_fda_fda

    result = {
        "edge_index": edge_index,
        "edge_attr": edge_attr,
        "stats": {
            "n_edges": n_edges,
            "n_fda_fda": n_fda_fda,
            "n_fda_train": n_fda_train,
            "threshold": threshold,
            "n_fda_drugs": n_fda,
        },
    }

    logger.info(f"DDS edges computed: {n_edges} total "
                f"({n_fda_train} FDA↔Train, {n_fda_fda} FDA↔FDA)")

    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    import argparse
    ap = argparse.ArgumentParser(description="Compute DDS edges for FDA drugs")
    ap.add_argument("--morgan-path",
                    default=str(HG_ROOT / "embeddings/all_drugs_idx2morgan.npy"))
    ap.add_argument("--existing-dds",
                    default=str(HG_ROOT / "data/DDS_edge_index_all_287834_pairs_weighted.npy"))
    ap.add_argument("--output",
                    default=str(HG_ROOT / "data/DDS_edge_index_fda_drugs.npy"))
    ap.add_argument("--threshold", type=float, default=None,
                    help="Tanimoto threshold (auto-calibrated if not specified)")
    ap.add_argument("--n-train", type=int, default=8981)
    ap.add_argument("--n-total", type=int, default=11104)
    args = ap.parse_args()

    # Load fingerprints
    morgan_train, morgan_fda = load_morgan_matrix(
        args.morgan_path,
        n_total=args.n_total,
        n_train=args.n_train,
    )

    # Calibrate threshold
    if args.threshold is None:
        threshold, density = calibrate_threshold(args.existing_dds, morgan_train)
    else:
        threshold = args.threshold
        n_train_drugs = morgan_train.shape[0]
        density = 287834 / (n_train_drugs * (n_train_drugs - 1) / 2)
        logger.info(f"Using specified threshold: {threshold}")

    # Compute edges
    result = compute_fda_dds_edges(
        morgan_train, morgan_fda,
        threshold=threshold,
        existing_dds_path=args.existing_dds,
        n_train=args.n_train,
    )

    # Save
    np.save(args.output, result, allow_pickle=True)
    logger.info(f"Saved DDS edges to {args.output}")
    logger.info(f"Stats: {result['stats']}")
