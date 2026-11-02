"""Figures for the first-round report — rendered from the frozen aggregation outputs.

Four small PNGs: per-fold headline AUPRC by init (with the `delta_mean` bar), the paired
per-drug gain distribution, the dose / Tahoe-seen strata, and per-fold seed stability. Each
figure is a pure function of the CSV tables written beside it, so it can be regenerated from
the tracked record alone.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")  # headless: never needs a display

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

__all__ = ["make_figures"]

_INIT_COLORS = {"random": "#6c757d", "tahoe": "#0b7285"}
_INIT_LABELS = {"random": "random init", "tahoe": "Tahoe init"}


def _primary_bar(summary: Mapping[str, Any]) -> float | None:
    value = ((summary.get("primary_bar") or {}).get("auprc_200"))
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def make_figures(
    *,
    fold_seed: pd.DataFrame,
    unit_results: pd.DataFrame,
    drug_gains: pd.DataFrame,
    strata: pd.DataFrame,
    summary: Mapping[str, Any],
    out_dir: Path | str,
) -> list[dict[str, Any]]:
    """Render the four report figures; returns ``[{"name","path","sha256"}]``."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[dict[str, Any]] = []
    bar = _primary_bar(summary)

    if not fold_seed.empty:
        folds = sorted(fold_seed["fold"].unique())
        fig, ax = plt.subplots(figsize=(6.4, 3.6))
        width = 0.36
        for offset, (init, color) in enumerate(
            (("random", _INIT_COLORS["random"]), ("tahoe", _INIT_COLORS["tahoe"]))
        ):
            subset = fold_seed[fold_seed["init"] == init]
            positions = np.arange(len(folds)) + (offset - 0.5) * width
            means = [
                float(subset.loc[subset["fold"] == fold, "auprc_200_mean"].iloc[0])
                if not subset.loc[subset["fold"] == fold].empty
                else np.nan
                for fold in folds
            ]
            sds = [
                float(subset.loc[subset["fold"] == fold, "auprc_200_sd"].iloc[0] or 0.0)
                if not subset.loc[subset["fold"] == fold].empty
                else 0.0
                for fold in folds
            ]
            ax.bar(
                positions, means, width,
                yerr=sds, capsize=2, color=color, label=_INIT_LABELS[init],
                error_kw={"elinewidth": 0.8},
            )
        if bar is not None:
            ax.axhline(bar, color="#c92a2a", linestyle="--", linewidth=1.0,
                       label=f"delta_mean baseline ({bar:.4f})")
        ax.set_xticks(np.arange(len(folds)), [f"fold {fold}" for fold in folds])
        ax.set_ylabel("direction AUPRC@200")
        ax.set_title("Test AUPRC@200 by fold (mean over 3 seeds, ±1 sd)")
        ax.legend(fontsize=7)
        ax.grid(axis="y", alpha=0.25)
        fig.tight_layout()
        written.append(_save(fig, out_dir / "auprc_by_fold.png"))

    if not drug_gains.empty and "gain_tahoe_minus_random_mean" in drug_gains.columns:
        values = (
            pd.to_numeric(drug_gains["gain_tahoe_minus_random_mean"], errors="coerce")
            .dropna().to_numpy(dtype=np.float64)
        )
        stats = ((summary.get("questions") or {}).get("q_a_tahoe_vs_random") or {}).get(
            "overall", {}
        )
        fig, ax = plt.subplots(figsize=(6.4, 3.6))
        ax.hist(values, bins=min(40, max(10, len(values) // 3)), color="#0b7285",
                alpha=0.85, edgecolor="white", linewidth=0.4)
        if stats.get("mean") is not None:
            ax.axvline(float(stats["mean"]), color="#c92a2a", linewidth=1.2,
                       label=f"mean {stats['mean']:+.4f}")
            if stats.get("ci_low") is not None:
                ax.axvline(float(stats["ci_low"]), color="#c92a2a", linestyle=":", linewidth=0.9)
                ax.axvline(float(stats["ci_high"]), color="#c92a2a", linestyle=":", linewidth=0.9)
            ax.set_title(
                f"Paired Tahoe−random gain by drug (n={stats.get('n')}, "
                f"95% CI [{stats.get('ci_low'):+.3f}, {stats.get('ci_high'):+.3f}], "
                f"verdict: {((summary.get('questions') or {}).get('q_a_tahoe_vs_random') or {}).get('verdict')})"
            )
            ax.legend(fontsize=7)
        ax.axvline(0.0, color="black", linewidth=0.8)
        ax.set_xlabel("drug-level AUPRC@200 gain (Tahoe − random)")
        ax.set_ylabel("drugs")
        ax.grid(axis="y", alpha=0.25)
        fig.tight_layout()
        written.append(_save(fig, out_dir / "paired_gain_by_drug.png"))

    if not strata.empty:
        known = strata[strata["label_known"]].copy()
        known["stratum_label"] = known["stratum_type"] + "=" + known["stratum_value"]
        fig, ax = plt.subplots(figsize=(7.2, 3.6))
        positions = np.arange(len(known))
        ax.bar(positions - 0.2, known["random_auprc_200"], width=0.4,
               color=_INIT_COLORS["random"], label="random init")
        ax.bar(positions + 0.2, known["tahoe_auprc_200"], width=0.4,
               color=_INIT_COLORS["tahoe"], label="Tahoe init")
        if bar is not None:
            ax.axhline(bar, color="#c92a2a", linestyle="--", linewidth=1.0,
                       label=f"delta_mean ({bar:.4f})")
        ax.set_xticks(positions, known["stratum_label"], rotation=20, ha="right", fontsize=7)
        ax.set_ylabel("direction AUPRC@200")
        ax.set_title("Dose and Tahoe-seen strata (test)")
        ax.legend(fontsize=7)
        ax.grid(axis="y", alpha=0.25)
        fig.tight_layout()
        written.append(_save(fig, out_dir / "strata.png"))

    if not unit_results.empty:
        folds = sorted(unit_results["fold"].unique())
        fig, ax = plt.subplots(figsize=(6.4, 3.6))
        width = 0.36
        for offset, init in enumerate(("random", "tahoe")):
            subset = unit_results[unit_results["init"] == init]
            for position, fold in enumerate(folds):
                group = subset[subset["fold"] == fold].sort_values("seed")
                ax.scatter(
                    np.full(len(group), position + (offset - 0.5) * width),
                    group["auprc_200"],
                    s=14, color=_INIT_COLORS[init],
                    label=_INIT_LABELS[init] if position == 0 else None,
                    zorder=3,
                )
        if bar is not None:
            ax.axhline(bar, color="#c92a2a", linestyle="--", linewidth=1.0,
                       label=f"delta_mean ({bar:.4f})")
        ax.set_xticks(np.arange(len(folds)), [f"fold {fold}" for fold in folds])
        ax.set_ylabel("direction AUPRC@200")
        ax.set_title("Seed stability: one point per (fold, seed)")
        ax.legend(fontsize=7)
        ax.grid(axis="y", alpha=0.25)
        fig.tight_layout()
        written.append(_save(fig, out_dir / "seed_stability.png"))
    return written


def _save(fig, path: Path) -> dict[str, Any]:
    fig.savefig(path, dpi=150)
    plt.close(fig)
    import common

    return {"name": path.name, "path": str(path), "sha256": common.sha256_file(path)}
