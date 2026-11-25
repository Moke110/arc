"""Rendering primitives for figures whose plotted values are recorded in CSV."""

import csv
from pathlib import Path
import numpy as np

from . import palette


def read_rows(path):
    """Read the CSV's rows without loading any study execution assets."""
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def save_figure(figure, path, **kwargs):
    """Write a figure to a destination, creating its parent directory if needed."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, **kwargs)


def bar(axis, x, height, *, hue=220, depth=2, **kwargs):
    """Draw bars with secondary faces and main edges.

    A scalar hue/depth styles the entire series; sequences supply one
    coordinate per bar and may be mixed with a scalar.
    """
    if np.ndim(hue) or np.ndim(depth):
        shape = np.broadcast(x, height).shape
        hues = np.broadcast_to(hue, shape)
        depths = np.broadcast_to(depth, shape)
        colors = [palette.color_pair(h, d) for h, d in zip(hues, depths, strict=True)]
        facecolor = [pair.secondary for pair in colors]
        edgecolor = [pair.main for pair in colors]
    else:
        colors = palette.color_pair(hue, depth)
        facecolor = colors.secondary
        edgecolor = colors.main
    kwargs.setdefault("linewidth", .8)
    return axis.bar(x, height, facecolor=facecolor, edgecolor=edgecolor, **kwargs)


def grouped_bar(axis, x, groups, *, width=.8, **kwargs):
    """Center groups of (heights, hue, depth, label) on numeric x positions."""
    groups = list(groups)
    if not groups:
        return []
    positions = np.asarray(x)
    bar_width = width / len(groups)
    center = (len(groups) - 1) / 2
    return [
        bar(axis, positions + (index - center) * bar_width, heights,
            hue=hue, depth=depth, label=label, width=bar_width, **kwargs)
        for index, (heights, hue, depth, label) in enumerate(groups)
    ]


def line(axis, x, y, *, hue=220, depth=2, **kwargs):
    """Draw a series with the main color of its palette coordinate."""
    return axis.plot(x, y, color=palette.color_pair(hue, depth).main, **kwargs)


def learning_curve(axis, runs, *, train_hue=220, validation_hue=140):
    """Draw one train/validation pair per run for a single model's panel.

    Each run supplies (epochs, train_values, validation_values). Hue
    identifies the curve role, and both lines in a run share a stroke
    pattern. Return the (training, validation) artists in input order.
    """
    styles = ("-", "--", "-.", ":")
    pairs = []
    for index, (epochs, training, validation) in enumerate(runs):
        style = styles[index % len(styles)]
        train_line, = line(axis, epochs, training, hue=train_hue, linestyle=style,
                           linewidth=1, label="Training" if index == 0 else None)
        validation_line, = line(axis, epochs, validation, hue=validation_hue,
                                linestyle=style, linewidth=1.2,
                                label="Validation" if index == 0 else None)
        pairs.append((train_line, validation_line))
    return pairs


def errorbar(axis, x, y, *, hue=220, depth=2, **kwargs):
    """Draw errors and optional markers with the same main color."""
    color = palette.color_pair(hue, depth).main
    return axis.errorbar(x, y, color=color, ecolor=color, **kwargs)


def range_band(axis, x, lower, upper, *, hue=220, depth=2, **kwargs):
    """Fill the observed interval with the coordinate's secondary color."""
    kwargs.setdefault("linewidth", 0)
    return axis.fill_between(x, lower, upper,
                             color=palette.color_pair(hue, depth).secondary, **kwargs)


def scatter(axis, x, y, *, hue=220, depth=2, **kwargs):
    """Draw points with the main color of their palette coordinate."""
    return axis.scatter(x, y, color=palette.color_pair(hue, depth).main, **kwargs)


def replicate_dots(axis, x, values, *, hue=220, depth=2, spread=.18, **kwargs):
    """Draw each finite replicate at one group center with one hue-depth pair.

    The symmetric spread is measured in x-axis units. Sorting by value keeps
    the offsets independent of replicate identifiers and input order.
    """
    if not np.isfinite(spread) or spread < 0:
        raise ValueError("replicate dots need a finite nonnegative spread")
    y = np.sort(_finite_replicates(values))
    offsets = np.linspace(-spread, spread, size) if (size := len(y)) > 1 else np.zeros(size)
    return scatter(axis, x + offsets, y, hue=hue, depth=depth, **kwargs)


def _finite_replicates(values):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1:
        raise ValueError("replicates need a value vector")
    return values[np.isfinite(values)]


def replicate_errorbar(axis, x, values, *, hue=220, depth=2, spread=.18, **kwargs):
    """Draw finite replicate dots and sample SD about their mean.

    Keyword arguments style the dots. An empty group draws no marks and
    a singleton draws only its dot. Return the dots and error bar artists.
    """
    y = _finite_replicates(values)
    if not len(y):
        return None, None
    kwargs.setdefault("s", 15)
    kwargs.setdefault("zorder", 4)
    errors = None
    if len(y) > 1:
        errors = errorbar(axis, [x], [y.mean()], yerr=[y.std(ddof=1)], fmt="none",
                          hue=hue, depth=depth, capsize=3, linewidth=1,
                          zorder=kwargs["zorder"] - 1)
    dots = replicate_dots(axis, x, y, hue=hue, depth=depth, spread=spread, **kwargs)
    return dots, errors


def replicate_bar(axis, x, values, *, hue=220, depth=2, width=.8, spread=None,
                  label=None, dot_size=15, **kwargs):
    """Draw the replicate mean, individual dots and sample SD for one bar.

    Undefined replicates produce no marks. Return the bar, dots and error
    bar artists, with None for marks that were not drawn.
    """
    y = _finite_replicates(values)
    if not len(y):
        return None, None, None
    kwargs.setdefault("zorder", 2)
    dots, errors = replicate_errorbar(
        axis, x, y, hue=hue, depth=depth, spread=width * .23 if spread is None else spread,
        s=dot_size, zorder=kwargs["zorder"] + 2,
    )
    bars = bar(axis, [x], [y.mean()], hue=hue, depth=depth, width=width, label=label, **kwargs)
    return bars, dots, errors


def bubble(axis, x, y, sizes, *, hue=220, depth=2, **kwargs):
    """Draw size-encoded points with the main color for both faces and edges."""
    color = palette.color_pair(hue, depth).main
    kwargs.setdefault("linewidths", .8)
    return axis.scatter(x, y, s=sizes, facecolors=color, edgecolors=color, **kwargs)


def heatmap(axis, values, *, hue=220, vmin=None, vmax=None, colorbar_label=None, **kwargs):
    """Draw a scalar heatmap and its numeric color scale in the ordered depths of one hue."""
    from matplotlib.colors import ListedColormap

    colormap = ListedColormap(palette.sequential_colors(hue))
    kwargs.setdefault("interpolation", "nearest")
    image = axis.imshow(values, cmap=colormap, vmin=vmin, vmax=vmax, **kwargs)
    axis.get_figure().colorbar(image, ax=axis, label=colorbar_label)
    return image


def histogram(axis, values, *, hue=220, depth=2, **kwargs):
    """Draw a histogram with secondary faces and main edges."""
    colors = palette.color_pair(hue, depth)
    kwargs.setdefault("linewidth", .8)
    return axis.hist(values, facecolor=colors.secondary, edgecolor=colors.main, **kwargs)
