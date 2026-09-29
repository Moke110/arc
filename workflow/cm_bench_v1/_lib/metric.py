"""
Direction-aware AUPRC: the single Stage 2 v1 metric.

For a (dataset, drug, dose) group the true label vector over the common gene
universe is split into an **up set** (`mean Δ > 0`, significant) and a **down
set** (`mean Δ < 0`, significant). The model supplies a ranking (its predicted
Δ); for the up set the predicted values are the scores, for the down set the
negated values are. The two average-precision values are averaged.

Threshold-free on the predicted side (only a ranking is used); with a constant
prediction the average precision equals the true-set prevalence, which is
therefore the random baseline. AUPRC is only comparable within a fixed universe
and split.
"""

import numpy as np
from sklearn.metrics import average_precision_score


def direction_aware_auprc(pred, labels):
    """AUPRC recovering the up/down significant sets; NaN if no significant gene.

    pred: (n_genes,) predicted Δ (higher = more up).
    labels: (n_genes,) int8, 1 up, -1 down, 0 not significant.
    """
    pred = np.asarray(pred, dtype=float)
    labels = np.asarray(labels)
    vals = []
    up = labels == 1
    down = labels == -1
    if up.any():
        vals.append(average_precision_score(up.astype(int), pred))
    if down.any():
        vals.append(average_precision_score(down.astype(int), -pred))
    return float(np.mean(vals)) if vals else float("nan")


def prevalence(labels):
    """Random-baseline AUPRC for a label vector (mean of the two directions)."""
    labels = np.asarray(labels)
    n = labels.size
    rates = []
    for m in (labels == 1, labels == -1):
        rates.append(m.sum() / n)
    return float(np.mean(rates)) if n else float("nan")
