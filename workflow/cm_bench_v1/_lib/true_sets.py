"""
Ground-truth true significant DEG sets per (dataset, drug, dose).

Per group (all treated wells of one dataset+drug+dose), a per-gene **moderated
one-sample t-test** of the per-well Δ vectors against 0 (limma/eBayes style,
per the spec) gives a p-value; Benjamini-Hochberg FDR < `DEG_FDR` marks
significance, and the sign of the mean Δ gives the direction. The up
(`mean Δ > 0`) and down (`mean Δ < 0`) significant sets are the labels the
direction-aware AUPRC recovers.

The `arc` conda env's bundled R/limma is broken (its R wrapper points at the
removed `xpert` env), so the moderated t is implemented directly in Python
following Smyth (2004) — the empirical-Bayes variance shrinkage is identical in
form; the p-values reduce to an ordinary one-sample t when the prior is diffuse.
"""

import numpy as np
from scipy import special
from scipy import stats


def benjamini_hochberg(p):
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)
    out = np.empty_like(q)
    out[order] = q
    return out


def _trigamma_inverse(x):
    if not np.isfinite(x) or x <= 0:
        return np.inf
    if x > 1e7:
        return 1.0 / np.sqrt(x)
    y = 0.5 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        dif = tri * (1.0 - special.polygamma(2, y) * (tri / x - 1.0) / tri)
        y = y + (x - tri) / dif
    return y


def fit_f_dist(s2, df):
    """limma `fitFDist` prior on per-gene variances (constant residual df).

    Returns (d0, s0sq). `df` is the common residual degrees of freedom.
    """
    s2 = np.asarray(s2, dtype=float)
    z = np.log(s2)
    e = z - special.digamma(df / 2.0) + np.log(df / 2.0)
    emean, evar = float(np.mean(e)), float(np.var(e, ddof=1))
    if evar > special.polygamma(1, df / 2.0):
        d0 = 2.0 * _trigamma_inverse(evar - special.polygamma(1, df / 2.0))
        s0sq = float(np.exp(emean + special.digamma(d0 / 2.0) - np.log(d0 / 2.0)))
    else:
        d0 = np.inf
        s0sq = float(np.exp(emean))
    return d0, s0sq


def moderated_ttest(delta):
    """Moderated one-sample t across replicate wells.

    `delta` is (n_wells, n_genes). Returns (t, p, df_total, mean).
    """
    x = np.asarray(delta, dtype=float)
    n = x.shape[0]
    df = n - 1
    mean = x.mean(axis=0)
    if n < 2:
        nan = np.full(x.shape[1], np.nan)
        return nan, nan, float(df), mean
    s2 = np.maximum(x.var(axis=0, ddof=1), 1e-12)
    d0, s0sq = fit_f_dist(s2, df)
    if not np.isfinite(d0) or d0 <= 0:
        s2_post, df_total = s2, float(df)
    else:
        s2_post = (d0 * s0sq + df * s2) / (d0 + df)
        df_total = df + d0
    t = mean / np.sqrt(s2_post / n)
    p = 2.0 * stats.t.sf(np.abs(t), df_total)
    return t, p, float(df_total), mean


def group_true_sets(delta, fdr=0.05, abs_mean_gate=0.0):
    """Up/down significant labels for one (dataset, drug, dose) group.

    Returns (labels int8 over genes: 1 up, -1 down, 0 n.s., n_sig_up, n_sig_down).
    """
    _, p, _, mean = moderated_ttest(delta)
    q = benjamini_hochberg(p)
    sig = np.isfinite(q) & (q < fdr) & (np.abs(mean) > abs_mean_gate)
    labels = np.zeros(mean.shape, dtype=np.int8)
    labels[sig & (mean > 0)] = 1
    labels[sig & (mean < 0)] = -1
    return labels, int((labels == 1).sum()), int((labels == -1).sum())
