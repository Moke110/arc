#!/usr/bin/env python3
"""Precision/recall curves + AUPRC for the v2 **official warm-start folds**.

Source data (task decision A): the 10 warm-start folds of the ADR 0018
*dtoxs-primary benchmark* in `05_xpert_finetune_bench/full/` — the same folds
the deck's `make_v2_warmstart_figure.py` reads from `full/output/benchmark.csv`.
That bench never dumped per-cell predictions, so step `dump` re-runs XPert in
`--mode test` on the *existing* fold weights (`full/output/cm_warm-start/`,
no retraining) with `--output_profile True` and collects the per-cell
`y_true / y_pred / ctl_true / deg_true / deg_pred` profiles.

Metric contract (task decision (a), matching XPert's own `precision_k`
convention in `pert_model/xpert/src/metrics.py`):

  * Per cell, the true set is the top-100 up genes (`deg_true` descending) and
    the bottom-100 down genes; every gene not in the top/bottom-100 is a
    negative. Both directions therefore have a fixed prevalence of
    100/978 = 0.1022, which is the random-ranking AUPRC baseline.
  * Up-direction score = `deg_pred`; down-direction score = `-deg_pred`
    (sign-aware). One average-precision value per direction.
  * Per-cell direction-aware AP = equal-weight mean of the two directions; the
    two components are reported separately too.
  * Aggregation: per cell -> per fold (mean over the fold's cells) -> across
    folds (mean +/- sd). AP is `sklearn.average_precision_score` (the step-wise
    average-precision convention, as in Stage 2's Direction-aware AUPRC), not
    the trapezoidal area under the PR curve; the curve itself is plotted.

A reproduction gate recomputes XPert's own `Precision@100 Positive/Negative_deg`
from the dumped profiles and compares against `full/output/benchmark.csv`.

Outputs under `warmstart_pr/output/`:
    profiles/cm_warm-start_fold-<i>_profile.npy
    auprc_per_cell.csv        per cell, both directions + direction-aware AP
    auprc_per_fold.csv        per fold summary
    pr_curves.csv             recall grid + across-fold mean/sd precision curves
    auprc_warmstart.png       PR curves + AUPRC bars (random baseline dashed)
    report.md                 numbers, recipe, reproduction check

Usage:
    python 05_xpert_finetune_bench/warmstart_pr/run_pr.py            # dump if needed, then metrics
    python 05_xpert_finetune_bench/warmstart_pr/run_pr.py dump
    python 05_xpert_finetune_bench/warmstart_pr/run_pr.py metrics
    python 05_xpert_finetune_bench/warmstart_pr/run_pr.py dump --force
"""

import argparse
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, precision_recall_curve

_HERE = Path(__file__).resolve().parent
STEP_DIR = _HERE.parent
sys.path.insert(0, str(STEP_DIR.parent / "_lib"))
sys.path.insert(0, str(STEP_DIR / "full"))
from utils import WORKFLOW_ROOT  # noqa: E402
from xpert_finetune_bench import DEFAULT_PRETRAINED, run_xpert  # noqa: E402

SPLIT = "warm-start"
N_FOLDS = 10
CONFIG = WORKFLOW_ROOT / "06_xpert_finetune" / "config" / "cm"
WEIGHTS_DIR = STEP_DIR / "full" / "output" / f"cm_{SPLIT}"
BENCH_CSV = STEP_DIR / "full" / "output" / "benchmark.csv"

OUT_DIR = _HERE / "output"
PROFILE_DIR = OUT_DIR / "profiles"
STAGING_ROOT = WORKFLOW_ROOT / "finetune" / "output"
STAGING = STAGING_ROOT / f"cm_{SPLIT}"

TOP_K = 100          # XPert precision_k num_pos = num_neg = 100
N_UNIVERSE = 978     # L1000 landmark space
RECALL_GRID = np.linspace(0.0, 1.0, 201)

BAR = {"up": "#2E6DA4", "down": "#C0504D"}


def _profile_path(fold):
    return PROFILE_DIR / f"cm_{SPLIT}_fold-{fold}_profile.npy"


# ---------------------------------------------------------------------------
# Step 1 — dump per-cell prediction profiles from the existing fold weights
# ---------------------------------------------------------------------------
def dump(device, force):
    pths = sorted(WEIGHTS_DIR.glob("*_best.pth"))
    if not pths:
        sys.exit(f"ERROR: no fold weights in {WEIGHTS_DIR} — train the full/ "
                 "warm-start bench first")
    if not force and all(_profile_path(f).exists() for f in range(1, N_FOLDS + 1)):
        print(f"All {N_FOLDS} profiles present in {PROFILE_DIR} — skipping dump")
        return
    PROFILE_DIR.mkdir(parents=True, exist_ok=True)
    shutil.rmtree(STAGING, ignore_errors=True)
    STAGING.mkdir(parents=True, exist_ok=True)
    # --pretrained_model selects the 'finetune' staging root (06 config ->
    # workflow root); the evaluated weights come from --saved_model_path.
    run_xpert("test", SPLIT, device,
              extra=["--saved_model_path", str(WEIGHTS_DIR),
                     "--pretrained_model", str(pths[0]),
                     "--output_profile", "True"],
              config=CONFIG)
    # `--mode test` runs validate(), which saves <fold_name>_predict_profile.npy
    # (fold_name = e.g. "warm-start_fold-1"), not the infer-mode name.
    produced = sorted(STAGING.glob("*_predict_profile.npy"))
    if not produced:
        sys.exit(f"ERROR: test mode produced no prediction profile in {STAGING}")
    for npy in produced:
        fold = npy.name.split("_fold-")[-1].split("_predict_profile")[0]
        dest = _profile_path(int(fold))
        if dest.exists():
            dest.unlink()
        shutil.move(str(npy), str(dest))
    got = sorted(_profile_path(f).name for f in range(1, N_FOLDS + 1)
                 if _profile_path(f).exists())
    print(f"Dumped {len(got)}/{N_FOLDS} profiles -> {PROFILE_DIR}")
    if len(got) < N_FOLDS:
        print(f"WARNING: expected {N_FOLDS} profiles")
    shutil.rmtree(STAGING, ignore_errors=True)


# ---------------------------------------------------------------------------
# Step 2 — direction-aware PR curves + AUPRC
# ---------------------------------------------------------------------------
def _true_sets(deg_true, k=TOP_K):
    """Top-k up and bottom-k down gene indices per cell (XPert convention)."""
    order = np.argsort(deg_true, axis=1)          # ascending
    down = order[:, :k]                            # bottom-k
    up = order[:, -k:]                             # top-k
    return up, down


def _mask(idx_rows, n_genes):
    m = np.zeros((idx_rows.shape[0], n_genes), dtype=bool)
    np.put_along_axis(m, idx_rows, True, axis=1)
    return m


def _pr_curve(labels, scores):
    """Precision interpolated on RECALL_GRID (ascending recall)."""
    precision, recall, _ = precision_recall_curve(labels, scores)
    r = recall[::-1]
    p = precision[::-1]
    grid = np.interp(RECALL_GRID, r, p)
    return grid, p, r


def score_fold(deg_true, deg_pred):
    """Per-cell direction-aware AP + averaged PR curves for one fold.

    Returns (per_cell DataFrame, up_curve, down_curve, p_at_100_up, p_at_100_down).
    """
    up_idx, down_idx = _true_sets(deg_true)
    up_mask = _mask(up_idx, deg_true.shape[1])
    down_mask = _mask(down_idx, deg_true.shape[1])
    rows, up_curves, down_curves = [], [], []
    for c in range(deg_true.shape[0]):
        ap_up = average_precision_score(up_mask[c].astype(int), deg_pred[c])
        ap_down = average_precision_score(down_mask[c].astype(int), -deg_pred[c])
        up_curves.append(_pr_curve(up_mask[c].astype(int), deg_pred[c])[0])
        down_curves.append(_pr_curve(down_mask[c].astype(int), -deg_pred[c])[0])
        rows.append({"cell": c, "ap_up": ap_up, "ap_down": ap_down,
                     "ap_direction_aware": 0.5 * (ap_up + ap_down)})
    per_cell = pd.DataFrame(rows)
    # reproduction gate: XPert's Precision@100 = |top100_pred ∩ top100_true| / 100
    pred_order = np.argsort(deg_pred, axis=1)
    pred_up = pred_order[:, -TOP_K:]
    pred_down = pred_order[:, :TOP_K]
    p100_up = np.mean([len(set(pred_up[i]) & set(up_idx[i])) / TOP_K
                       for i in range(deg_true.shape[0])])
    p100_down = np.mean([len(set(pred_down[i]) & set(down_idx[i])) / TOP_K
                         for i in range(deg_true.shape[0])])
    return (per_cell, np.vstack(up_curves).mean(axis=0),
            np.vstack(down_curves).mean(axis=0), p100_up, p100_down)


def metrics():
    per_cell_all, per_fold_rows = [], []
    up_by_fold, down_by_fold = [], []
    for fold in range(1, N_FOLDS + 1):
        path = _profile_path(fold)
        if not path.exists():
            sys.exit(f"ERROR: missing profile {path} — run the dump step first")
        d = np.load(path, allow_pickle=True).item()
        deg_true = np.asarray(d["deg_true"], dtype=float)
        deg_pred = np.asarray(d["deg_pred"], dtype=float)
        per_cell, up_c, down_c, p100_up, p100_down = score_fold(deg_true, deg_pred)
        per_cell.insert(0, "fold", fold)
        per_cell_all.append(per_cell)
        up_by_fold.append(up_c)
        down_by_fold.append(down_c)
        per_fold_rows.append({
            "fold": fold, "n_cells": len(per_cell),
            "ap_up": per_cell["ap_up"].mean(),
            "ap_down": per_cell["ap_down"].mean(),
            "ap_direction_aware": per_cell["ap_direction_aware"].mean(),
            "precision_at_100_up": p100_up,
            "precision_at_100_down": p100_down,
        })
        print(f"[fold {fold}] n_cells={len(per_cell):4d} "
              f"AP_up={per_cell['ap_up'].mean():.4f} "
              f"AP_down={per_cell['ap_down'].mean():.4f} "
              f"dir-aware={per_cell['ap_direction_aware'].mean():.4f}")

    per_cell_all = pd.concat(per_cell_all, ignore_index=True)
    per_fold = pd.DataFrame(per_fold_rows)

    # across-fold mean/sd PR curves (each fold = mean over its cells)
    up_mean = np.mean(up_by_fold, axis=0)
    up_sd = np.std(up_by_fold, axis=0, ddof=1)
    down_mean = np.mean(down_by_fold, axis=0)
    down_sd = np.std(down_by_fold, axis=0, ddof=1)
    curves = pd.DataFrame({
        "recall": RECALL_GRID,
        "up_mean": up_mean, "up_sd": up_sd,
        "down_mean": down_mean, "down_sd": down_sd,
    })

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    per_cell_all.to_csv(OUT_DIR / "auprc_per_cell.csv", index=False)
    per_fold.to_csv(OUT_DIR / "auprc_per_fold.csv", index=False)
    curves.to_csv(OUT_DIR / "pr_curves.csv", index=False)

    prevalence = TOP_K / N_UNIVERSE
    summary = {
        "ap_up": (per_fold["ap_up"].mean(), per_fold["ap_up"].std(ddof=1)),
        "ap_down": (per_fold["ap_down"].mean(), per_fold["ap_down"].std(ddof=1)),
        "ap_direction_aware": (per_fold["ap_direction_aware"].mean(),
                               per_fold["ap_direction_aware"].std(ddof=1)),
    }
    for k, (m, s) in summary.items():
        print(f"{k:20s} = {m:.4f} +/- {s:.4f}")
    print(f"random baseline      = {prevalence:.4f}")

    _plot(curves, per_fold, prevalence)
    _report(per_fold, summary, prevalence)
    return per_fold, summary, prevalence


# Font scale for the deck figure: axis fonts x2, every other font x1.75.
AXIS_FS = 20        # tick labels + axis titles (default ~10 -> x2)
OTHER_FS = 17.5     # general text (default 10 -> x1.75)
TITLE_FS = 22       # axes titles (~12.5 -> x1.75)
LEGEND_FS = 18      # legend (~10.5 -> x1.75)
NOTE_FS = 19        # bar annotations (~11 -> x1.75)


def _plot(curves, per_fold, prevalence):
    plt.rcParams.update({
        "font.size": OTHER_FS,
        "axes.labelsize": AXIS_FS,
        "xtick.labelsize": AXIS_FS,
        "ytick.labelsize": AXIS_FS,
        "axes.titlesize": TITLE_FS,
    })
    fig, axes = plt.subplots(1, 2, figsize=(16.0, 7.8),
                             constrained_layout=True)

    ax = axes[0]
    for tag, col in (("up", BAR["up"]), ("down", BAR["down"])):
        m = curves[f"{tag}_mean"].values
        s = curves[f"{tag}_sd"].values
        ax.fill_between(curves["recall"], m - s, m + s, color=col, alpha=0.18,
                        linewidth=0)
        ax.plot(curves["recall"], m, color=col, lw=2.2,
                label=f"{tag}  (AUPRC {per_fold[f'ap_{tag}'].mean():.3f})")
    ax.axhline(prevalence, ls="--", lw=1.6, color="#555555",
               label=f"random ({prevalence:.3f})")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_box_aspect(1)          # square plot area (x:y = 1:1)
    ax.set_xlabel("recall")
    ax.set_ylabel("precision")
    ax.set_title("PR curve", fontsize=TITLE_FS, fontweight="bold", pad=14)
    ax.legend(loc="upper right", frameon=False, fontsize=LEGEND_FS)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1]
    cats = [("Up", "ap_up", BAR["up"]), ("Down", "ap_down", BAR["down"]),
            ("Direction-aware", "ap_direction_aware", "#4F6228")]
    x = np.arange(len(cats))
    means = [per_fold[col].mean() for _, col, _ in cats]
    stds = [per_fold[col].std(ddof=1) for _, col, _ in cats]
    ax.bar(x, means, yerr=stds, capsize=6, width=0.6,
           color=[c for _, _, c in cats], alpha=0.9, edgecolor="white",
           zorder=2, error_kw=dict(ecolor="#333333", lw=1.5, zorder=4))
    rng = np.random.default_rng(2024)
    for xi, (_, col, _) in zip(x, cats):
        ax.scatter(xi + rng.uniform(-0.13, 0.13, size=len(per_fold)),
                   per_fold[col].values, color="black", s=24, zorder=3,
                   alpha=0.75)
    for xi, m, s in zip(x, means, stds):
        ax.text(xi, m + s + 0.012, f"{m:.3f}\n±{s:.3f}", ha="center",
                va="bottom", fontsize=NOTE_FS, linespacing=0.95)
    rand = ax.axhline(prevalence, ls="--", lw=1.6, color="#555555")
    ax.set_xticks(x)
    ax.set_xticklabels([c for c, _, _ in cats], fontsize=AXIS_FS)
    ax.set_ylabel("AUPRC")
    ax.set_ylim(0, max(0.45, max(means) + max(stds) + 0.12))
    ax.set_title("AUPRC", fontsize=TITLE_FS, fontweight="bold", pad=14)
    ax.legend([rand], [f"random {prevalence:.3f}"], loc="upper left",
              frameon=False, fontsize=LEGEND_FS)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    out = OUT_DIR / "auprc_warmstart.png"
    fig.savefig(out, dpi=200, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def _report(per_fold, summary, prevalence):
    bench = pd.read_csv(BENCH_CSV)
    ws = bench[bench["split"] == SPLIT]
    lines = [
        "# v2 warm-start folds — direction-aware PR curves & AUPRC",
        "",
        f"Source: the 10 official warm-start folds of the dtoxs-primary "
        f"benchmark (`05_xpert_finetune_bench/full/`), re-evaluated in "
        f"`--mode test --output_profile True` from the existing fold weights "
        f"(no retraining). Cells: "
        f"{int(per_fold['n_cells'].sum())} across 10 folds.",
        "",
        "## Metric",
        "",
        f"Per cell, the true set is the top-{TOP_K} up genes and bottom-{TOP_K} "
        f"down genes by `deg_true` (XPert `precision_k` convention), in the "
        f"{N_UNIVERSE}-gene landmark space; the up direction is scored by "
        f"`deg_pred`, the down direction by `-deg_pred`. AUPRC = "
        f"`sklearn.average_precision_score` (average precision, not the "
        f"trapezoidal PR area). Aggregation: per cell -> per fold -> across "
        f"folds (mean ± sd). Random baseline = {TOP_K}/{N_UNIVERSE} = "
        f"{prevalence:.4f}.",
        "",
        "## AUPRC (mean ± sd over 10 folds)",
        "",
        "| direction | AUPRC | sd |",
        "|---|---|---|",
    ]
    for k, lab in (("ap_up", "up"), ("ap_down", "down"),
                   ("ap_direction_aware", "direction-aware (equal-weight mean)")):
        m, s = summary[k]
        lines.append(f"| {lab} | {m:.4f} | {s:.4f} |")
    lines += ["", f"Random baseline: **{prevalence:.4f}**.", "",
              "## Per-fold", "",
              "| fold | cells | AP up | AP down | direction-aware |",
              "|---|---|---|---|---|"]
    for _, r in per_fold.iterrows():
        lines.append(f"| {int(r['fold'])} | {int(r['n_cells'])} | "
                     f"{r['ap_up']:.4f} | {r['ap_down']:.4f} | "
                     f"{r['ap_direction_aware']:.4f} |")

    lines += ["", "## Reproduction gate", "",
              "Recomputed XPert `Precision@100 Positive/Negative_deg` from the "
              "dumped profiles vs `full/output/benchmark.csv` (mean over folds):",
              "",
              "| metric | from profiles | benchmark.csv |",
              "|---|---|---|"]
    gp_up = per_fold["precision_at_100_up"].mean()
    gp_dn = per_fold["precision_at_100_down"].mean()
    bk_up = ws["Precision@100 Positive_deg"].astype(float).mean()
    bk_dn = ws["Precision@100 Negative_deg"].astype(float).mean()
    lines += [f"| Precision@100 Positive_deg | {gp_up:.4f} | {bk_up:.4f} |",
              f"| Precision@100 Negative_deg | {gp_dn:.4f} | {bk_dn:.4f} |", ""]
    d_up = abs(gp_up - bk_up)
    d_dn = abs(gp_dn - bk_dn)
    ok = d_up < 0.01 and d_dn < 0.01
    lines += [f"Match within 0.01: **{'yes' if ok else 'NO'}** "
              f"(Δ up {d_up:.4f}, Δ down {d_dn:.4f}); benchmark.csv rounds "
              f"Precision@k to 3 decimals.", ""]
    (OUT_DIR / "report.md").write_text("\n".join(lines) + "\n")
    print(f"Saved: {OUT_DIR / 'report.md'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("action", nargs="?", default="all",
                    choices=["all", "dump", "metrics"])
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--force", action="store_true",
                    help="re-dump profiles even if cached")
    args = ap.parse_args()
    if args.action in ("all", "dump"):
        dump(args.device, args.force)
    if args.action in ("all", "metrics"):
        metrics()


if __name__ == "__main__":
    main()
