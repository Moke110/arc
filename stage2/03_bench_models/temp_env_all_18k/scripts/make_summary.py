#!/usr/bin/env python
"""Render the Stage 2 benchmark summary: one metrics table + one bar panel per metric.

Reads ``results/leaderboard.csv`` and the model registry next to this script,
then writes ``summary_table.csv`` and ``summary.png`` into the same
``temp_env_all_18k/`` archive.  Ordering is identical in the table and every bar
panel: baselines first (alphabetical), then the real models (alphabetical).
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent
RESULTS = ARCHIVE / "results"
CONFIG = ARCHIVE / "configs" / "models.json"

BASE_COLOR = "#7b9fc4"   # baselines
MODEL_COLOR = "#e08a4a"  # real models
GRID_COLOR = "#d9d9d9"

# name -> (leaderboard column, panel title, higher_is_better, axis transform)
PANELS = [
    ("headline_auprc", "AUPRC@200 (headline)", True, None),
    ("mean_k100", "AUPRC@100", True, None),
    ("mean_k500", "AUPRC@500", True, None),
    ("pcc", "PCC", True, None),
    ("spearman", "Spearman", True, None),
    ("r2", "R²", True, "clip"),
    ("rmse", "RMSE", False, "log"),
    ("nrmse", "NRMSE", False, "log"),
]


def load_models() -> tuple[pd.DataFrame, list[str], list[str], dict[str, bool]]:
    frame = pd.read_csv(RESULTS / "leaderboard.csv")
    registry = json.loads(CONFIG.read_text())
    is_baseline = {name: bool(registry.get(name, {}).get("baseline", False)) for name in frame["model"]}
    baselines = sorted(name for name, flag in is_baseline.items() if flag)
    models = sorted(name for name, flag in is_baseline.items() if not flag)
    order = baselines + models
    missing = set(order) - set(frame["model"])
    if missing:
        raise SystemExit(f"leaderboard is missing models: {sorted(missing)}")
    frame = frame.set_index("model").loc[order].reset_index()
    return frame, baselines, models, is_baseline


def fmt(value: float, digits: int = 3) -> str:
    if value is None or not np.isfinite(float(value)):
        return "n/a"
    value = float(value)
    if value != 0 and (abs(value) >= 1e4 or abs(value) < 1e-3):
        return f"{value:.2e}"
    return f"{value:.{digits}f}"


def compact(value: float) -> str:
    """Short bar label: 4 significant digits, scientific only when unavoidable."""
    if not np.isfinite(value):
        return "n/a"
    if value == 0:
        return "0"
    magnitude = abs(value)
    if magnitude >= 1e4 or magnitude < 1e-3:
        mantissa, exponent = f"{value:.1e}".split("e")
        return f"{mantissa}e{int(exponent):+d}"
    if magnitude >= 100:
        return f"{value:.0f}"
    if magnitude >= 10:
        return f"{value:.1f}"
    return f"{value:.3f}".rstrip("0").rstrip(".")


def draw_table(ax, frame: pd.DataFrame, is_baseline: dict[str, bool]) -> None:
    ax.axis("off")
    columns = [
        ("Rank", 0.045),
        ("Model", 0.125),
        ("Type", 0.072),
        ("AUPRC@100", 0.083),
        ("AUPRC@200", 0.088),
        ("95% CI", 0.105),
        ("AUPRC@500", 0.088),
        ("PCC", 0.063),
        ("Spearman", 0.074),
        ("R²", 0.076),
        ("RMSE", 0.068),
        ("NRMSE", 0.068),
        ("Fidelity", 0.065),
    ]

    def cell(row, key: str) -> str:
        if key == "type":
            return "baseline" if is_baseline[str(row["model"])] else "model"
        if key == "ci":
            return f"[{float(row['ci_low']):.3f}, {float(row['ci_high']):.3f}]"
        value = row[key]
        if key == "rank":
            return str(int(value))
        if key == "model":
            return str(value)
        if key in {"headline_auprc", "mean_k100", "mean_k500"}:
            return f"{float(value):.4f}"
        if key == "ideal_fit":
            return str(value)
        return fmt(float(value))

    keys = ["rank", "model", "type", "mean_k100", "headline_auprc", "ci", "mean_k500", "pcc", "spearman", "r2", "rmse", "nrmse", "ideal_fit"]
    cell_text = [[cell(row, key) for key in keys] for _, row in frame.iterrows()]

    # best value per comparable metric (baselines and models together)
    best = {}
    for key, higher in [
        ("mean_k100", True), ("headline_auprc", True), ("mean_k500", True),
        ("pcc", True), ("spearman", True), ("r2", True),
        ("rmse", False), ("nrmse", False),
    ]:
        values = pd.to_numeric(frame[key], errors="coerce")
        if values.notna().any():
            best[key] = values.max() if higher else values.min()

    table = ax.table(
        cellText=cell_text,
        colLabels=[name for name, _ in columns],
        colWidths=[width for _, width in columns],
        cellLoc="center",
        loc="upper center",
        bbox=[0.0, 0.02, 1.0, 0.98],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(17.2)

    best_cells = {
        "mean_k100": 3, "headline_auprc": 4, "mean_k500": 6,
        "pcc": 7, "spearman": 8, "r2": 9, "rmse": 10, "nrmse": 11,
    }
    for (row, col), cell_obj in table.get_celld().items():
        cell_obj.set_edgecolor("#c8c8c8")
        cell_obj.set_linewidth(0.6)
        if row == 0:
            cell_obj.set_facecolor("#e9edf2")
            cell_obj.set_text_props(weight="bold")
            continue
        model = frame.iloc[row - 1]["model"]
        cell_obj.set_facecolor("#f2f6fa" if is_baseline[model] else "#fdf3ea")
        if col in best_cells.values():
            key = [k for k, c in best_cells.items() if c == col][0]
            value = pd.to_numeric(pd.Series([frame.iloc[row - 1][key]]), errors="coerce").iloc[0]
            if np.isfinite(value) and float(value) == float(best[key]):
                cell_obj.set_facecolor("#d6ecd4")
                cell_obj.set_text_props(weight="bold")

    ax.set_title(
        "Stage 2 benchmark — metric summary (ordering: baselines first, then models alphabetically)",
        fontsize=26,
        weight="bold",
        pad=42,
    )
    ax.text(
        0.0, 1.005,
        "Green = best value for that metric. AUPRC = direction-macro, group-macro; headline = k=200 (5 cold-drug folds).",
        transform=ax.transAxes, fontsize=18, color="#444444", va="bottom",
    )


def draw_panel(ax, frame: pd.DataFrame, is_baseline: dict[str, bool], spec) -> None:
    key, title, higher, transform = spec
    names = list(frame["model"])
    values = pd.to_numeric(frame[key], errors="coerce").to_numpy(dtype=float)
    colors = [BASE_COLOR if is_baseline[name] else MODEL_COLOR for name in names]
    x = np.arange(len(names))

    ax.set_facecolor("#fbfbfb")
    r2_floor = -15.0
    heights = np.where(np.isfinite(values), values, 0.0)
    if transform == "clip":
        heights = np.clip(heights, r2_floor, None)
    if transform == "log":
        bottom = np.nanmin(values[np.isfinite(values)]) / 4.0
        ax.set_yscale("log")
        ax.set_ylim(bottom=bottom)
    ax.bar(x, heights, color=colors, edgecolor="white", linewidth=0.7, width=0.72)
    ax.axhline(0.0, color="#555555", linewidth=0.8)
    ax.grid(axis="y", color=GRID_COLOR, linewidth=0.7, alpha=0.9)
    ax.set_axisbelow(True)

    if key == "headline_auprc":
        low = frame["ci_low"].to_numpy(dtype=float)
        high = frame["ci_high"].to_numpy(dtype=float)
        err = np.vstack([values - low, high - values])
        ax.errorbar(x, values, yerr=err, fmt="none", ecolor="#333333", elinewidth=1.0, capsize=2.5, capthick=1.0)

    high = frame["ci_high"].to_numpy(dtype=float) if key == "headline_auprc" else None
    for xi, value in zip(x, values):
        if not np.isfinite(value):
            ax.text(xi, 0.02, "n/a", ha="center", va="bottom", fontsize=14,
                    color="#666666", transform=ax.get_xaxis_transform())
            continue
        anchor = value if high is None else max(value, high[xi])
        if transform == "clip":
            anchor = max(anchor, r2_floor)
        va = "bottom" if value >= 0 else "top"
        offset = 2 if value >= 0 else -2
        rotation = 90
        if transform == "clip":
            if value < r2_floor:
                va, offset = "bottom", 5
            elif value >= 0:
                anchor, va, offset = 0.0, "bottom", 3
        ax.annotate(
            compact(value), (xi, anchor), xytext=(0, offset),
            textcoords="offset points", ha="center", va=va,
            fontsize=13.2, rotation=rotation,
        )

    finite = values[np.isfinite(values)]
    if transform == "log":
        ax.set_ylim(top=float(finite.max()) * 4.5)
    elif transform == "clip":
        ax.set_ylim(r2_floor - 3.0, 4.0)
    else:
        lo = min(0.0, float(finite.min()))
        hi = float(finite.max())
        if high is not None:
            hi = max(hi, float(np.nanmax(high)))
        span = max(hi - lo, 1e-9)
        ax.set_ylim(lo - (0.06 * span if lo < 0 else 0.0), hi + 0.32 * span)

    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=90, ha="center", va="top", fontsize=18)
    for label, name in zip(ax.get_xticklabels(), names):
        label.set_color("#2f5d8a" if is_baseline[name] else "#a95a1c")
    ax.set_title(f"{title}  ({'higher' if higher else 'lower'} is better)"
                 + (" · log scale" if transform == "log" else "")
                 + (f" · clipped at {r2_floor:.0f}" if transform == "clip" else ""),
                 fontsize=21, pad=12)
    ax.tick_params(axis="y", labelsize=16)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main() -> int:
    frame, baselines, models, is_baseline = load_models()

    out = frame.copy()
    out.insert(1, "type", ["baseline" if is_baseline[m] else "model" for m in out["model"]])
    out = out[[
        "rank", "model", "type", "mean_k100", "headline_auprc", "headline_sd",
        "ci_low", "ci_high", "mean_k500", "pcc", "spearman", "r2", "rmse",
        "nrmse", "ideal_fit", "status",
    ]]
    out.to_csv(ARCHIVE / "summary_table.csv", index=False, lineterminator="\n")

    fig = plt.figure(figsize=(32, 24))
    grid = GridSpec(2, 1, figure=fig, height_ratios=[0.95, 1.05], hspace=0.16)
    draw_table(fig.add_subplot(grid[0]), frame, is_baseline)

    sub = grid[1].subgridspec(2, 4, hspace=1.05, wspace=0.30)
    for index, spec in enumerate(PANELS):
        ax = fig.add_subplot(sub[index // 4, index % 4])
        draw_panel(ax, frame, is_baseline, spec)

    fig.suptitle(
        "Stage 2 cold-drug benchmark — 12 entrants × 5 folds (baselines vs. models)",
        fontsize=34, weight="bold", y=0.985,
    )
    fig.text(
        0.5, -0.055,
        "Baselines (blue): " + ", ".join(baselines)
        + "   |   Models (orange): " + ", ".join(models)
        + "\nRMSE is in raw log1p_cp10k Δ units; no_change has undefined PCC/Spearman (constant prediction). "
          "Squidiff's R² is clipped at −15 (true value −4.9e5); its R²/RMSE reflect the 20/100,000-step smoke budget.",
        ha="center", fontsize=19, color="#444444",
    )
    fig.savefig(ARCHIVE / "summary.png", dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {ARCHIVE / 'summary_table.csv'}")
    print(f"wrote {ARCHIVE / 'summary.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
