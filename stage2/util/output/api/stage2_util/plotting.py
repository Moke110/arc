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


def bubble(axis, x, y, sizes, *, hue=220, depth=2, **kwargs):
    """Draw size-encoded points with secondary faces and main edges."""
    colors = palette.color_pair(hue, depth)
    kwargs.setdefault("linewidths", .8)
    return axis.scatter(x, y, s=sizes, facecolors=colors.secondary,
                        edgecolors=colors.main, **kwargs)


def heatmap(axis, values, *, hue=220, vmin=None, vmax=None, colorbar_label=None, **kwargs):
    """Draw a scalar heatmap and its numeric color scale in five depths of one hue."""
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
