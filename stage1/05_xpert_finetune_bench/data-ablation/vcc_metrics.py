#!/usr/bin/env python3
"""VCC 2026 metric computation (PDS/JAC/FID/REACH), adapted for the XPert
data-ablation test sets.

Definitions follow Arc Institute's Virtual Cell Challenge 2026 scoring suite
(cell_eval2 v0.16.0, ``vcc2026`` profile); the formulas and constants are
documented in ``/home/chang/vcc/reference/vcc2026_metrics.html``. Two
adaptations are made for the ablation data:

  * Drug perturbations have no "target gene", so the panel-wide target-gene
    exclusion (``exclusion_scope=panel``) is a no-op and is dropped.
  * The ablation test set holds ~5 cells per drug, so a Wilcoxon ``p_adj <
    0.05`` differential-expression set is underpowered (a signed-rank test at
    n=5 cannot reach 0.05). The DE gene set for JAC/FID/REACH is therefore the
    top-k genes by |log2FC| (k = DE_TOP_K), matching the repo's existing
    Precision@k convention.

All four returned values are RAW and higher-is-better (no replicate-anchored
scaling).
"""

import numpy as np

DE_TOP_K = 50
REACH_PURITY_FLOOR = 0.9


def _midrank0(a, v):
    """0-based midrank of value ``v`` within 1-D array ``a``; ties share the
    average position. ``v`` is expected to be a member of ``a``."""
    lt = int(np.sum(a < v))
    eq = int(np.sum(a == v))
    return lt + 0.5 * (eq - 1)


def _mean_effect(deg, groups, drugs):
    """Per-drug mean effect vectors, shape (n_drugs, n_genes)."""
    out = np.zeros((len(drugs), deg.shape[1]), dtype=np.float64)
    for i, d in enumerate(drugs):
        m = groups == d
        out[i] = deg[m].mean(axis=0)
    return out


def pds(eff_pred, eff_real):
    """Perturbation discrimination: mean over drugs of
    ``1 - midrank(d_pp) / (n - 1)`` under cosine distance between predicted
    and real per-drug mean effects."""
    n = eff_real.shape[0]
    if n < 2:
        return np.nan
    pn = eff_pred / (np.linalg.norm(eff_pred, axis=1, keepdims=True) + 1e-12)
    rn = eff_real / (np.linalg.norm(eff_real, axis=1, keepdims=True) + 1e-12)
    cos = pn @ rn.T
    dist = 1.0 - cos
    scores = []
    for p in range(n):
        rank = _midrank0(dist[p], dist[p, p])
        scores.append(1.0 - rank / (n - 1))
    return float(np.mean(scores))


def _topk(deg_vec, k):
    return np.argsort(np.abs(deg_vec))[::-1][:k]


def jac(eff_pred, eff_real, k=DE_TOP_K):
    """DE significance overlap (Jaccard of top-k |log2FC| gene sets)."""
    n = eff_real.shape[0]
    vals = []
    for p in range(n):
        r = set(_topk(eff_real[p], k).tolist())
        q = set(_topk(eff_pred[p], k).tolist())
        union = len(r | q)
        vals.append(len(r & q) / union if union else 1.0)
    return float(np.mean(vals))


def fid(eff_pred, eff_real, k=DE_TOP_K):
    """DE direction fidelity: fraction of the top-k predicted genes whose
    direction (sign) matches the real log2FC."""
    n = eff_real.shape[0]
    vals = []
    for p in range(n):
        pred_top = _topk(eff_pred[p], k)
        r = eff_real[p][pred_top]
        f = eff_pred[p][pred_top]
        correct = int(np.sum((np.sign(f) == np.sign(r)) & (r != 0)))
        vals.append(correct / k)
    return float(np.mean(vals))


def reach(eff_pred, eff_real, k=DE_TOP_K, floor=REACH_PURITY_FLOOR):
    """DE direction reach: deepest prefix of the predicted |log2FC| ranking
    whose direction purity stays >= floor, divided by k."""
    n = eff_real.shape[0]
    vals = []
    for p in range(n):
        order = np.argsort(np.abs(eff_pred[p]))[::-1]
        r = eff_real[p][order]
        f = eff_pred[p][order]
        adjudicable = r != 0
        correct = np.sign(f) == np.sign(r)
        cnt = 0
        n_correct = 0
        k_star = 0
        for i in range(len(order)):
            if not adjudicable[i]:
                continue
            cnt += 1
            n_correct += int(correct[i])
            if n_correct / cnt >= floor:
                k_star = cnt
            if cnt >= k:
                break
        vals.append(k_star / k)
    return float(np.mean(vals))


def compute_vcc_metrics(deg_true, deg_pred, drugs):
    """Compute the four raw VCC metrics from per-cell DEG predictions.

    Args:
        deg_true: (n_cells, n_genes) true perturbation effect (deg = trt - ctl).
        deg_pred: (n_cells, n_genes) predicted perturbation effect.
        drugs: (n_cells,) integer drug index per cell.

    Returns:
        dict with keys PDS/JAC/FID/REACH (each averaged over drugs).
    """
    deg_true = np.asarray(deg_true, dtype=np.float64)
    deg_pred = np.asarray(deg_pred, dtype=np.float64)
    drugs = np.asarray(drugs)
    uniq = np.unique(drugs)
    eff_real = _mean_effect(deg_true, drugs, uniq)
    eff_pred = _mean_effect(deg_pred, drugs, uniq)
    return {
        "PDS": pds(eff_pred, eff_real),
        "JAC": jac(eff_pred, eff_real),
        "FID": fid(eff_pred, eff_real),
        "REACH": reach(eff_pred, eff_real),
    }
