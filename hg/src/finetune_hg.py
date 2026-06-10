"""
Fine-tune HG on the extended graph (11,104 drugs) to produce properly-scaled
FDA drug embeddings while preserving training drug embeddings.

Key design:
  - Full-graph forward pass (no NeighborLoader) — graph fits in GPU memory
  - InfoNCE contrastive loss on sampled DTI/DDS/PPI edges
  - Very low LR (5e-5) to prevent training drug embedding drift
  - Cosine similarity monitoring against reference embeddings
  - Early stop if training drug cosim drops below threshold

Usage:
  cd hg/
  python -m src.finetune_hg \
      --device cuda:0 --lr 5e-5 --epochs 50

  # Test with fewer epochs:
  python -m src.finetune_hg --device cuda:0 --epochs 5 --check-interval 1
"""

import logging
import os
import sys
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F

HG_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HG_ROOT))

from src.model_hg import HeteroGNN
from src.utils_hg import ContrastiveLoss, get_embeddings
from src.infer_fda_hg import generate_fda_node_features, build_extended_graph

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Edge sampling for contrastive loss
# ---------------------------------------------------------------------------
def sample_positive_edges(
    edge_index: torch.Tensor,
    n_sample: int,
    device: torch.device,
) -> torch.Tensor:
    """Randomly sample positive edges from edge_index.

    Args:
        edge_index: (2, N_edges) tensor
        n_sample: number of edges to sample
        device: torch device

    Returns:
        (2, n_sample) sampled edge_index
    """
    n_total = edge_index.shape[1]
    if n_total <= n_sample:
        return edge_index
    idx = torch.randperm(n_total, device=device)[:n_sample]
    return edge_index[:, idx]


def generate_negatives_for_edges(
    pos_src: torch.Tensor,
    pos_dst: torch.Tensor,
    num_dst_nodes: int,
    num_neg_samples: int,
    device: torch.device,
) -> torch.Tensor:
    """Generate negative destination samples for each positive edge.

    Args:
        pos_src: (N_pos,) source indices
        pos_dst: (N_pos,) destination indices
        num_dst_nodes: total number of destination nodes
        num_neg_samples: number of negatives per positive edge
        device: torch device

    Returns:
        (N_pos, num_neg_samples) negative destination indices
    """
    n_pos = len(pos_src)
    neg_dst = torch.randint(0, num_dst_nodes, (n_pos, num_neg_samples), device=device)

    # Avoid sampling the positive dst as negative (unlikely but fix it)
    pos_dst_expanded = pos_dst.unsqueeze(1)  # (N_pos, 1)
    collision = (neg_dst == pos_dst_expanded)
    # Replace collisions with a different random index
    n_collisions = collision.sum().item()
    if n_collisions > 0:
        replacements = torch.randint(0, num_dst_nodes, (n_collisions,), device=device)
        neg_dst[collision] = replacements

    return neg_dst


def compute_edge_type_loss(
    x_dict: dict,
    edge_type: tuple,
    sampled_edges: torch.Tensor,
    num_neg_samples: int,
    contrastive_loss_fn: ContrastiveLoss,
) -> torch.Tensor:
    """Compute InfoNCE loss for a single edge type batch.

    Args:
        x_dict: node embeddings dict {'drug': (N,256), 'gene': (M,256)}
        edge_type: ('src_type', 'relation', 'dst_type')
        sampled_edges: (2, N_sample) positive edge_index
        num_neg_samples: negatives per positive
        contrastive_loss_fn: InfoNCE loss module

    Returns:
        scalar loss tensor
    """
    src_type, _, dst_type = edge_type

    pos_src = sampled_edges[0]
    pos_dst = sampled_edges[1]

    z_i = x_dict[src_type][pos_src]  # (N_pos, hidden_dim)
    z_j = x_dict[dst_type][pos_dst]  # (N_pos, hidden_dim)

    num_dst_nodes = x_dict[dst_type].size(0)
    neg_dst_idx = generate_negatives_for_edges(
        pos_src, pos_dst, num_dst_nodes, num_neg_samples, z_i.device
    )
    z_neg = x_dict[dst_type][neg_dst_idx]  # (N_pos, num_neg, hidden_dim)

    loss = contrastive_loss_fn(z_i, z_j, z_neg)
    return loss


# ---------------------------------------------------------------------------
# Cosine similarity monitoring
# ---------------------------------------------------------------------------
def compute_mean_cosine(u: np.ndarray, v: np.ndarray) -> float:
    """Compute mean per-row cosine similarity between two embeddings arrays.

    Args:
        u: (N, D) array
        v: (N, D) array

    Returns:
        scalar mean cosine similarity
    """
    u_norm = u / (np.linalg.norm(u, axis=1, keepdims=True) + 1e-10)
    v_norm = v / (np.linalg.norm(v, axis=1, keepdims=True) + 1e-10)
    cos_sims = np.sum(u_norm * v_norm, axis=1)
    return float(cos_sims.mean())


# ---------------------------------------------------------------------------
# Main fine-tuning function
# ---------------------------------------------------------------------------
def finetune_on_extended_graph(
    data_dir: str,
    model_path: str,
    dti_fda_path: str | None,
    dds_fda_path: str | None,
    kpgt_path: str,
    output_dir: str,
    device: str = "cuda:0",
    lr: float = 5e-5,
    epochs: int = 50,
    temperature: float = 0.5,
    num_neg_samples: int = 5,
    dds_sample_per_epoch: int = 20000,
    ppi_sample_per_epoch: int = 50000,
    cosim_threshold: float = 0.95,
    check_interval: int = 5,
    include_edges: list[str] | None = None,
    seed: int = 4242,
) -> dict:
    """Fine-tune HG model on the extended graph.

    Returns:
        dict with results (train_cosim, fda_std, train_std, dti_ratio, output_paths)
    """
    if include_edges is None:
        include_edges = ["PPI", "DTI", "DDS"]

    device_obj = torch.device(device if torch.cuda.is_available() else "cpu")
    os.makedirs(output_dir, exist_ok=True)

    # -----------------------------------------------------------------------
    # Step 1: Generate FDA drug node features
    # -----------------------------------------------------------------------
    logger.info("=" * 60)
    logger.info("Step 1: Generate 512-dim FDA drug features")
    logger.info("=" * 60)
    fda_feat_512 = generate_fda_node_features(
        kpgt_path,
        str(Path(data_dir) / "all_drug_node_feat_8981.npy"),
    )
    n_fda = fda_feat_512.shape[0]
    n_train = 8981
    n_total = n_train + n_fda
    logger.info(f"FDA drugs: {n_fda}, Train drugs: {n_train}, Total: {n_total}")

    # -----------------------------------------------------------------------
    # Step 2: Build extended graph
    # -----------------------------------------------------------------------
    logger.info("\n" + "=" * 60)
    logger.info("Step 2: Build extended heterogeneous graph")
    logger.info("=" * 60)
    data = build_extended_graph(
        data_dir,
        fda_feat_512,
        dti_fda_path=dti_fda_path,
        dds_fda_path=dds_fda_path,
        include_edges=include_edges,
    )
    data = data.to(device_obj)
    logger.info(f"Graph on device: {device_obj}")

    # -----------------------------------------------------------------------
    # Step 3: Load pretrained model
    # -----------------------------------------------------------------------
    logger.info("\n" + "=" * 60)
    logger.info("Step 3: Load pretrained model")
    logger.info("=" * 60)

    model = HeteroGNN(
        hidden_dim=256,
        out_dim=256,
        num_layers=3,
        include_edges=include_edges,
    ).to(device_obj)

    checkpoint = torch.load(model_path, map_location=device_obj, weights_only=False)
    model.load_state_dict(checkpoint["model_state_dict"])
    logger.info(f"Loaded checkpoint from {model_path} (epoch {checkpoint.get('epoch', '?')})")

    # -----------------------------------------------------------------------
    # Step 4: Save reference training embeddings
    # -----------------------------------------------------------------------
    logger.info("\n" + "=" * 60)
    logger.info("Step 4: Extract reference training embeddings")
    logger.info("=" * 60)

    model.eval()
    with torch.no_grad():
        ref_drug_emb, ref_gene_emb = get_embeddings(model, data)
    ref_train_emb = ref_drug_emb[:n_train].cpu().numpy()  # (8981, 256)
    ref_fda_emb = ref_drug_emb[n_train:].cpu().numpy()    # (2123, 256)

    logger.info(f"Reference train embeddings: shape={ref_train_emb.shape}, "
                f"mean={ref_train_emb.mean():.4f}, std={ref_train_emb.std():.4f}")
    logger.info(f"Reference FDA embeddings:   shape={ref_fda_emb.shape}, "
                f"mean={ref_fda_emb.mean():.4f}, std={ref_fda_emb.std():.4f}")
    logger.info(f"Scale ratio (train/fda std): {ref_train_emb.std() / (ref_fda_emb.std() + 1e-10):.2f}x")

    # Save reference embeddings
    np.save(os.path.join(output_dir, "ref_train_embeddings.npy"), ref_train_emb)
    np.save(os.path.join(output_dir, "ref_fda_embeddings.npy"), ref_fda_emb)

    # -----------------------------------------------------------------------
    # Step 5: Fine-tuning loop
    # -----------------------------------------------------------------------
    logger.info("\n" + "=" * 60)
    logger.info(f"Step 5: Fine-tuning (lr={lr}, epochs={epochs})")
    logger.info("=" * 60)

    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    contrastive_loss_fn = ContrastiveLoss(temperature)

    # Get edge index dict for sampling
    edge_index_dict = data.edge_index_dict

    history = {
        "epoch": [],
        "loss": [],
        "train_cosim": [],
        "fda_std": [],
        "train_std": [],
    }

    best_cosim = 1.0
    best_epoch = 0
    early_stop_reason = "completed"

    for epoch in range(1, epochs + 1):
        model.train()
        epoch_loss = 0.0
        n_edge_types = 0
        epoch_start = time.time()

        # ---------------------------------------------------------------
        # Full-graph forward pass
        # ---------------------------------------------------------------
        x_dict = model(data.x_dict, data.edge_index_dict)

        # ---------------------------------------------------------------
        # DTI loss (all edges — most important biological signal)
        # ---------------------------------------------------------------
        if ("drug", "DTI", "gene") in edge_index_dict:
            dti_edges = edge_index_dict[("drug", "DTI", "gene")]
            # Use all DTI edges (small enough: ~16,874)
            dti_loss = compute_edge_type_loss(
                x_dict, ("drug", "DTI", "gene"), dti_edges,
                num_neg_samples, contrastive_loss_fn,
            )
            if not torch.isnan(dti_loss):
                epoch_loss += dti_loss
                n_edge_types += 1

        # ---------------------------------------------------------------
        # DDS loss (sample 20K edges)
        # ---------------------------------------------------------------
        if ("drug", "DDS", "drug") in edge_index_dict:
            dds_edges = edge_index_dict[("drug", "DDS", "drug")]
            dds_sampled = sample_positive_edges(dds_edges, dds_sample_per_epoch, device_obj)
            dds_loss = compute_edge_type_loss(
                x_dict, ("drug", "DDS", "drug"), dds_sampled,
                num_neg_samples, contrastive_loss_fn,
            )
            if not torch.isnan(dds_loss):
                epoch_loss += dds_loss
                n_edge_types += 1

        # ---------------------------------------------------------------
        # PPI loss (sample 50K edges — maintains gene embedding structure)
        # ---------------------------------------------------------------
        if ("gene", "PPI", "gene") in edge_index_dict:
            ppi_edges = edge_index_dict[("gene", "PPI", "gene")]
            ppi_sampled = sample_positive_edges(ppi_edges, ppi_sample_per_epoch, device_obj)
            ppi_loss = compute_edge_type_loss(
                x_dict, ("gene", "PPI", "gene"), ppi_sampled,
                num_neg_samples, contrastive_loss_fn,
            )
            if not torch.isnan(ppi_loss):
                epoch_loss += ppi_loss
                n_edge_types += 1

        # ---------------------------------------------------------------
        # Backprop
        # ---------------------------------------------------------------
        if epoch_loss != 0.0 and not torch.isnan(epoch_loss):
            optimizer.zero_grad()
            epoch_loss.backward()
            optimizer.step()

        loss_val = epoch_loss.item() if epoch_loss != 0.0 else float("nan")
        epoch_time = time.time() - epoch_start

        # ---------------------------------------------------------------
        # Periodic cosine similarity check
        # ---------------------------------------------------------------
        if epoch % check_interval == 0 or epoch == epochs:
            model.eval()
            with torch.no_grad():
                curr_drug_emb, _ = get_embeddings(model, data)
            curr_train_emb = curr_drug_emb[:n_train].cpu().numpy()
            curr_fda_emb = curr_drug_emb[n_train:].cpu().numpy()

            train_cosim = compute_mean_cosine(curr_train_emb, ref_train_emb)
            fda_std = float(curr_fda_emb.std())
            train_std = float(curr_train_emb.std())

            history["epoch"].append(epoch)
            history["loss"].append(loss_val)
            history["train_cosim"].append(train_cosim)
            history["fda_std"].append(fda_std)
            history["train_std"].append(train_std)

            logger.info(
                f"Epoch {epoch:03d} | loss={loss_val:.4f} | "
                f"train_cosim={train_cosim:.4f} | "
                f"train_std={train_std:.1f} | fda_std={fda_std:.1f} | "
                f"time={epoch_time:.1f}s"
            )

            # Early stop if training embeddings drift too much
            if train_cosim < cosim_threshold:
                logger.warning(
                    f"Training drug cosim {train_cosim:.4f} < threshold {cosim_threshold} — "
                    f"early stop at epoch {epoch}"
                )
                early_stop_reason = f"cosim below threshold at epoch {epoch}"
                break

            # Track best
            if train_cosim > best_cosim:
                best_cosim = train_cosim
                best_epoch = epoch
                # Save best checkpoint
                torch.save({
                    "epoch": epoch,
                    "model_state_dict": model.state_dict(),
                    "optimizer_state_dict": optimizer.state_dict(),
                }, os.path.join(output_dir, "best_checkpoint.pth"))
        else:
            logger.info(f"Epoch {epoch:03d} | loss={loss_val:.4f} | time={epoch_time:.1f}s")

    # -----------------------------------------------------------------------
    # Step 6: Final evaluation + save
    # -----------------------------------------------------------------------
    logger.info("\n" + "=" * 60)
    logger.info("Step 6: Final evaluation and save")
    logger.info("=" * 60)

    # Load best checkpoint
    best_ckpt_path = os.path.join(output_dir, "best_checkpoint.pth")
    if os.path.exists(best_ckpt_path):
        best_ckpt = torch.load(best_ckpt_path, map_location=device_obj, weights_only=False)
        model.load_state_dict(best_ckpt["model_state_dict"])
        logger.info(f"Loaded best checkpoint from epoch {best_ckpt['epoch']}")

    model.eval()
    with torch.no_grad():
        final_drug_emb, final_gene_emb = get_embeddings(model, data)
    final_drug_emb = final_drug_emb.cpu().numpy()  # (11104, 256)
    final_gene_emb = final_gene_emb.cpu().numpy()  # (19392, 256)

    final_train_emb = final_drug_emb[:n_train]
    final_fda_emb = final_drug_emb[n_train:]

    # --- Compute final metrics ---
    train_cosim = compute_mean_cosine(final_train_emb, ref_train_emb)
    fda_cosim_vs_ref = compute_mean_cosine(final_fda_emb, ref_fda_emb)
    train_std = float(final_train_emb.std())
    fda_std = float(final_fda_emb.std())

    logger.info(f"Final Results:")
    logger.info(f"  Train drug cosim vs reference: {train_cosim:.4f} "
                f"({'OK' if train_cosim >= cosim_threshold else 'WARNING'})")
    logger.info(f"  FDA drug cosim vs reference:  {fda_cosim_vs_ref:.4f}")
    logger.info(f"  Train embedding std: {train_std:.1f}")
    logger.info(f"  FDA embedding std:   {fda_std:.1f}")
    logger.info(f"  FDA/Train std ratio:  {fda_std / (train_std + 1e-10):.3f}")
    logger.info(f"  Early stop reason: {early_stop_reason}")

    # --- DTI biological signal check ---
    dti_ratio = None
    if ("drug", "DTI", "gene") in edge_index_dict:
        dti_edges = edge_index_dict[("drug", "DTI", "gene")]
        dti_drug_idx = dti_edges[0].cpu().numpy()
        dti_gene_idx = dti_edges[1].cpu().numpy()

        dti_drug_emb = final_drug_emb[dti_drug_idx]
        dti_gene_emb = final_gene_emb[dti_gene_idx]

        dti_cos = np.sum(dti_drug_emb * dti_gene_emb, axis=1) / (
            np.linalg.norm(dti_drug_emb, axis=1) *
            np.linalg.norm(dti_gene_emb, axis=1) + 1e-10)

        n_random = min(50000, len(dti_cos) * 3)
        rand_drugs = np.random.RandomState(42).randint(0, n_total, n_random)
        rand_genes = np.random.RandomState(43).randint(0, final_gene_emb.shape[0], n_random)
        rand_cos = np.sum(
            final_drug_emb[rand_drugs] * final_gene_emb[rand_genes], axis=1) / (
            np.linalg.norm(final_drug_emb[rand_drugs], axis=1) *
            np.linalg.norm(final_gene_emb[rand_genes], axis=1) + 1e-10)

        dti_ratio = dti_cos.mean() / (rand_cos.mean() + 1e-10)
        logger.info(f"  DTI drug-gene cosim: {dti_cos.mean():.4f} "
                    f"vs random: {rand_cos.mean():.4f} "
                    f"(ratio: {dti_ratio:.2f}x)")

    # --- Save embeddings ---
    all_drug_path = os.path.join(output_dir, "all_drug_embeddings.npy")
    np.save(all_drug_path, final_drug_emb)
    logger.info(f"\nSaved all-drug embeddings to {all_drug_path}")

    all_gene_path = os.path.join(output_dir, "all_gene_embeddings.npy")
    np.save(all_gene_path, final_gene_emb)
    logger.info(f"Saved all-gene embeddings to {all_gene_path}")

    # Save final checkpoint
    torch.save({
        "epoch": epoch,
        "model_state_dict": model.state_dict(),
        "optimizer_state_dict": optimizer.state_dict(),
    }, os.path.join(output_dir, "final_checkpoint.pth"))

    # Save history
    import json
    with open(os.path.join(output_dir, "finetune_history.json"), "w") as f:
        json.dump(history, f, indent=2)

    return {
        "train_cosim": train_cosim,
        "fda_cosim_vs_ref": fda_cosim_vs_ref,
        "train_std": train_std,
        "fda_std": fda_std,
        "dti_ratio": dti_ratio,
        "early_stop_reason": early_stop_reason,
        "output_dir": output_dir,
        "all_drug_path": all_drug_path,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    import argparse
    ap = argparse.ArgumentParser(
        description="Fine-tune HG on extended graph for FDA drug embeddings"
    )
    ap.add_argument("--data-dir", default=str(HG_ROOT / "data"))
    ap.add_argument("--model-path",
                    default=str(HG_ROOT / "model/20260518_102121/global_fold_early_stop.pth"))
    ap.add_argument("--dti-fda", default=str(HG_ROOT / "data/DTI_edge_index_fda_ext.npy"))
    ap.add_argument("--dds-fda", default=str(HG_ROOT / "data/DDS_edge_index_fda_drugs.npy"))
    ap.add_argument("--kpgt-path",
                    default=str(HG_ROOT.parent / "xpert/drug_embeddings/all_drugs_idx2KPGT.npy"))
    ap.add_argument("--output-dir", default=str(HG_ROOT / "model/finetune"))
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--lr", type=float, default=5e-5)
    ap.add_argument("--epochs", type=int, default=50)
    ap.add_argument("--temperature", type=float, default=0.5)
    ap.add_argument("--num-neg-samples", type=int, default=5)
    ap.add_argument("--dds-sample-per-epoch", type=int, default=20000)
    ap.add_argument("--ppi-sample-per-epoch", type=int, default=50000)
    ap.add_argument("--cosim-threshold", type=float, default=0.95)
    ap.add_argument("--check-interval", type=int, default=5)
    ap.add_argument("--include-edges", nargs="+", default=["PPI", "DTI", "DDS"])
    ap.add_argument("--seed", type=int, default=4242)
    ap.add_argument("--update-xpert-embeddings", action="store_true",
                    help="Automatically update xpert/drug_embeddings/all_drugs_idx2HG.npy")

    args = ap.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    # Set seed
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    # Create output dir with timestamp
    timestamp = time.strftime("%Y%m%d_%H%M%S", time.localtime())
    output_dir = os.path.join(args.output_dir, timestamp)
    os.makedirs(output_dir, exist_ok=True)

    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Device: {args.device}")
    logger.info(f"LR: {args.lr}, Epochs: {args.epochs}, Cosim threshold: {args.cosim_threshold}")

    result = finetune_on_extended_graph(
        data_dir=args.data_dir,
        model_path=args.model_path,
        dti_fda_path=args.dti_fda,
        dds_fda_path=args.dds_fda,
        kpgt_path=args.kpgt_path,
        output_dir=output_dir,
        device=args.device,
        lr=args.lr,
        epochs=args.epochs,
        temperature=args.temperature,
        num_neg_samples=args.num_neg_samples,
        dds_sample_per_epoch=args.dds_sample_per_epoch,
        ppi_sample_per_epoch=args.ppi_sample_per_epoch,
        cosim_threshold=args.cosim_threshold,
        check_interval=args.check_interval,
        include_edges=args.include_edges,
        seed=args.seed,
    )

    # Optionally update XPert embeddings
    if args.update_xpert_embeddings and result["train_cosim"] >= args.cosim_threshold:
        xpert_hg_path = str(HG_ROOT.parent / "xpert/drug_embeddings/all_drugs_idx2HG.npy")
        backup_path = xpert_hg_path.replace(".npy", "_backup_graph.npy")
        import shutil
        shutil.copy2(xpert_hg_path, backup_path)
        logger.info(f"Backed up current HG embeddings to {backup_path}")

        # Load and update: keep train rows from fine-tuned, or all rows
        finetuned_all = np.load(result["all_drug_path"])
        np.save(xpert_hg_path, finetuned_all)
        logger.info(f"Updated {xpert_hg_path} with fine-tuned embeddings")

    # Summary
    logger.info("\n" + "=" * 60)
    logger.info("FINETUNING SUMMARY")
    logger.info("=" * 60)
    logger.info(f"  Train cosim:  {result['train_cosim']:.4f}  {'✓' if result['train_cosim'] >= args.cosim_threshold else '✗'}")
    logger.info(f"  Train std:    {result['train_std']:.1f}")
    logger.info(f"  FDA std:      {result['fda_std']:.1f}")
    if result["dti_ratio"]:
        logger.info(f"  DTI ratio:    {result['dti_ratio']:.2f}x  {'✓' if result['dti_ratio'] >= 2.0 else '✗'}")
    logger.info(f"  Stop reason:  {result['early_stop_reason']}")
    logger.info(f"  Output:       {result['output_dir']}")


if __name__ == "__main__":
    main()
