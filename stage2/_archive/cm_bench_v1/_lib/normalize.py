"""
Per-dataset raw normalization (log2 CPM+1).

Normalization is applied inside each dataset's own space *before* Δ; there is no
cross-dataset harmonization (spec: the per-plate baseline is part of the model
input and carries the dataset scale).
"""

import numpy as np


def log_cpm(counts, target_sum=1e6):
    """log2(CPM + 1) per sample; input (n_samples, n_genes) raw counts."""
    counts = np.asarray(counts, dtype=float)
    totals = counts.sum(axis=1, keepdims=True)
    totals[totals <= 0] = 1.0
    return np.log2(counts / totals * target_sum + 1.0)
