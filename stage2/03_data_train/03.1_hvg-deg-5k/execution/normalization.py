"""Normalize bounded count windows on their selected expression axis."""
import numpy as np


class ViewError(ValueError):
    """Counts or library sizes violate the normalization contract."""


def library_totals(counts):
    """Selected-axis library total of every profile (int64); a zero total fails closed."""
    counts = np.asarray(counts)
    if counts.ndim != 2:
        raise ViewError("library totals need a two-dimensional count block")
    if counts.size and counts.min() < 0:
        raise ViewError("corpus counts must be non-negative")
    totals = counts.sum(axis=1, dtype=np.int64)
    if (totals <= 0).any():
        raise ViewError("a corpus profile has a non-positive selected-axis library total")
    return totals


def log1p_cp10k(counts, totals):
    """``log1p(count * 10000 / library_total)`` on the supplied selected-axis totals."""
    counts = np.asarray(counts)
    totals = np.asarray(totals, dtype=np.float64)
    if counts.ndim != 2:
        raise ViewError("normalization needs a two-dimensional count block")
    if totals.shape != (counts.shape[0],):
        raise ViewError("one selected-axis library total per profile is required")
    if not np.isfinite(totals).all() or (totals <= 0).any():
        raise ViewError("library totals must be finite and positive")
    return np.log1p(counts.astype(np.float64) * (10000.0 / totals[:, None]))
