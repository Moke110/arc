"""Shared color assignments for the Step 05 report panels."""

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
