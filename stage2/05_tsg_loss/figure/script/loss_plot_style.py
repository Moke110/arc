"""Shared color assignments for the Step 05 report panels."""

from matplotlib.lines import Line2D

from stage2_util import palette

FOLD_HUES = (220, 20, 140)
LOWER_IS_BETTER = frozenset(("Δ RMSE", "Δ E-distance", "Δ MMD"))
LOWER_IS_BETTER_HUE = 340
HIGHER_IS_BETTER_HUE = 220


def metric_direction(axis, metric, label):
    """Apply the metric's background cue and label the preferred direction."""
    lower = metric in LOWER_IS_BETTER
    hue = LOWER_IS_BETTER_HUE if lower else HIGHER_IS_BETTER_HUE
    axis.set_facecolor(palette.color_pair(hue, 4).secondary)
    axis.set_title(f"{label} {'↓' if lower else '↑'}", fontweight="bold")
    return hue


def curve_legend_handles(folds, *, selected_epoch=False):
    """Build the shared fold and phase keys for a learning-curve legend."""
    handles = [
        Line2D([], [], color=palette.color_pair(FOLD_HUES[fold]).main, label=f"Fold {fold}")
        for fold in folds
    ]
    handles.extend((
        Line2D([], [], color="#333333", linestyle="--", label="Train"),
        Line2D([], [], color="#333333", label="Validation"),
    ))
    if selected_epoch:
        handles.append(Line2D([], [], color="#333333", linestyle="None",
                              marker="o", label="Selected epoch"))
    return handles
