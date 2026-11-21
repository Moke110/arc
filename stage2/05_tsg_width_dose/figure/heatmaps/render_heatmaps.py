"""Render Step 05 heatmaps from its published metric CSVs."""

import argparse
from itertools import product
from pathlib import Path
import sys

import matplotlib

matplotlib.use('Agg')
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np

FIGURE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FIGURE / 'script'))
from plot_contract import FOLDS, METHODS, METHOD_LABELS, METRICS, OVERALL, PER_DOSE, WIDTHS
from plot_contract import read_rows
from render import checked_index, numeric, save


def metric_scale(values, metric):
    finite = values[np.isfinite(values)]
    low, high = (finite.min(), finite.max()) if finite.size else (0.0, 1.0)
    if low == high:
        margin = max(abs(low) * 0.001, 0.0001)
        low, high = low - margin, high + margin
    palette = ['#a51f31', '#ffffff', '#17569a']
    if metric == 'rmse':
        palette.reverse()
    cmap = colors.LinearSegmentedColormap.from_list(metric, palette, N=257)
    cmap.set_bad('#e6e8ec')
    return colors.Normalize(low, high), cmap


def panel(axis, values, labels, scale, *, show_methods=True, show_rows=True):
    norm, cmap = scale
    image = axis.imshow(values, norm=norm, cmap=cmap, interpolation='nearest', aspect='auto')
    axis.set_xticks(range(len(METHODS)), METHOD_LABELS if show_methods else [''] * len(METHODS),
                    rotation=23, ha='right', fontsize=8)
    axis.set_yticks(range(len(labels)), labels if show_rows else [''] * len(labels), fontsize=8)
    axis.tick_params(length=0)
    for spine in axis.spines.values():
        spine.set_visible(False)
    for (row, column), value in np.ndenumerate(values):
        red, green, blue, _ = cmap(norm(value))
        channels = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
                    for c in (red, green, blue)]
        luminance = sum(c * weight for c, weight in zip(channels, (0.2126, 0.7152, 0.0722)))
        axis.text(column, row, f'{value:.3f}' if np.isfinite(value) else 'NA',
                  ha='center', va='center', color='white' if luminance < 0.179 else 'black',
                  fontsize=9)
    return image


def colorbar(figure, image, axes):
    figure.colorbar(image, ax=axes, orientation='horizontal', shrink=0.78,
                    aspect=24).ax.tick_params(labelsize=8)


def render(data, output):
    keys = ('width', 'method', 'metric')
    overall = checked_index(read_rows(data, OVERALL), keys, product(WIDTHS, METHODS, METRICS), OVERALL)
    folds = tuple(map(str, FOLDS))
    rows = read_rows(data, PER_DOSE)
    held = checked_index(rows, (*keys, 'fold'), product(WIDTHS, METHODS, METRICS, folds), PER_DOSE)
    concentrations = {fold: {numeric(row['dose_um']) for row in rows if row['fold'] == fold}
                      for fold in folds}
    if any(len(values) != 1 for values in concentrations.values()):
        raise ValueError('Inconsistent held concentrations')
    doses = {fold: values.pop() for fold, values in concentrations.items()}
    if len(set(doses.values())) != len(folds) or any(
        not np.isfinite(dose) or dose <= 0 for dose in doses.values()
    ):
        raise ValueError('Expected distinct positive held concentrations')
    ordered = sorted(folds, key=doses.__getitem__)

    output.mkdir(parents=True, exist_ok=True)
    summary, summary_axes = plt.subplots(1, 4, figsize=(13, 3.9), layout='constrained')
    by_dose, dose_axes = plt.subplots(len(ordered), len(METRICS), figsize=(13, 8.5),
                                     layout='constrained')
    summary.suptitle('TranSiGen width and dose encoding · equal-fold means')
    by_dose.suptitle('TranSiGen width and dose encoding · means by held concentration')
    summary_axes[0].set_ylabel('Width')
    for column, (metric, label) in enumerate(METRICS.items()):
        axis = summary_axes[column]
        axis.set_title(label)
        values = np.array([[numeric(overall[width, method, metric]['mean'])
                            for method in METHODS] for width in WIDTHS])
        colorbar(summary, panel(axis, values, WIDTHS, metric_scale(values, metric),
                                show_rows=column == 0), axis)

        values = np.array([[[numeric(held[width, method, metric, fold]['mean'])
                             for method in METHODS] for width in WIDTHS] for fold in ordered])
        scale = metric_scale(values, metric)
        for row, fold in enumerate(ordered):
            axis = dose_axes[row, column]
            axis.set_title(label)
            if column == 0:
                axis.set_ylabel(f'{doses[fold]:g} µM\nWidth')
            image = panel(axis, values[row], WIDTHS, scale, show_rows=column == 0)
        colorbar(by_dose, image, dose_axes[:, column])

    save(summary, output, 'metrics_heatmap')
    save(by_dose, output, 'metrics_per_dose_heatmap')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=FIGURE / 'data')
    parser.add_argument('--output', type=Path, default=FIGURE)
    args = parser.parse_args()
    render(args.data, args.output)
