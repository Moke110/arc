"""Merge multiple DTI batch output files into a single consolidated file.

The DTI query script produces per-batch .npy files with LOCAL FDA drug indices
(relative to the batch start position). This script:

  1. Reads each batch file
  2. Offsets the drug indices to their global position within the 2123 FDA-extension drugs
  3. Concatenates all edges
  4. Saves the merged file

Usage:
  python -m src.merge_dti_batches \
      data/DTI_edge_index_fda_ext_b1.npy \
      data/DTI_edge_index_fda_ext_b2.npy \
      data/DTI_edge_index_fda_ext_b3.npy \
      --output data/DTI_edge_index_fda_ext.npy \
      --offsets 0 200 400
"""

import logging
import sys
from pathlib import Path

import numpy as np

HG_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HG_ROOT))

logger = logging.getLogger(__name__)


def merge_dti_batches(batch_files: list[str],
                      output_path: str,
                      offsets: list[int] | None = None) -> dict:
    """Merge multiple DTI batch .npy files.

    Args:
        batch_files: list of .npy file paths (each from query_dti_fda.py)
        output_path: where to save the merged file
        offsets: list of int offsets to apply to LOCAL drug indices in each batch.
                 If None, reads offset from batch metadata (if available).

    Returns:
        dict with edge_index, gene_hg_indices, stats
    """
    if offsets is None:
        offsets = [0] * len(batch_files)

    all_edge_index = []
    all_gene_indices = []
    all_sources = []
    total_drugs = 0
    total_edges = 0
    drugs_with_targets = 0
    total_drugs_in_batches = 0

    for i, (fpath, offset) in enumerate(zip(batch_files, offsets)):
        if not Path(fpath).exists():
            logger.warning(f"Skipping missing file: {fpath}")
            continue

        data = np.load(fpath, allow_pickle=True).item()
        ei = data["edge_index"]  # (2, N_edges) with LOCAL drug indices
        gene_idx = data.get("gene_hg_indices", ei[1])
        sources = data.get("edge_sources", [])

        # Apply offset to drug indices
        ei_offset = ei.copy()
        ei_offset[0] = ei_offset[0] + offset

        n_edges = ei.shape[1]
        if n_edges > 0:
            all_edge_index.append(ei_offset)
            all_gene_indices.append(gene_idx)
            all_sources.extend(sources)

        stats = data.get("stats", {})
        batch_drugs = stats.get("n_fda_drugs", 0)
        batch_edges = stats.get("n_dti_edges", n_edges)
        batch_with = stats.get("drugs_with_targets", 0)

        logger.info(f"  Batch {i} (offset={offset}): {batch_edges} edges, "
                    f"{batch_with}/{batch_drugs} drugs with targets")

        total_edges += batch_edges
        drugs_with_targets += batch_with
        total_drugs_in_batches += batch_drugs

    # Concatenate
    if all_edge_index:
        merged_edges = np.concatenate(all_edge_index, axis=1)
        merged_genes = np.concatenate(all_gene_indices, axis=0)
    else:
        merged_edges = np.zeros((2, 0), dtype=np.int64)
        merged_genes = np.zeros(0, dtype=np.int64)

    # Compute unique genes
    unique_genes = len(set(merged_genes)) if len(merged_genes) > 0 else 0

    result = {
        "edge_index": merged_edges,
        "gene_hg_indices": merged_genes,
        "edge_sources": np.array(all_sources, dtype=object),
        "stats": {
            "n_fda_drugs": total_drugs_in_batches,
            "n_dti_edges": total_edges,
            "drugs_with_targets": drugs_with_targets,
            "drugs_without_targets": total_drugs_in_batches - drugs_with_targets,
            "unique_target_genes": unique_genes,
            "avg_targets_per_drug": total_edges / total_drugs_in_batches if total_drugs_in_batches > 0 else 0,
        },
    }

    np.save(output_path, result, allow_pickle=True)
    logger.info(f"\nMerged DTI: {total_edges} edges for {total_drugs_in_batches} drugs")
    logger.info(f"  Drugs with targets: {drugs_with_targets}")
    logger.info(f"  Unique target genes: {unique_genes}")
    logger.info(f"  Saved to {output_path}")

    return result


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("batch_files", nargs="+", help="DTI batch .npy files")
    ap.add_argument("--output", required=True, help="Output merged .npy path")
    ap.add_argument("--offsets", type=int, nargs="+", default=None,
                    help="Drug index offsets for each batch")
    args = ap.parse_args()

    merge_dti_batches(args.batch_files, args.output, args.offsets)
