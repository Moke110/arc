"""Reproduce the four Step 05 figures using the published CSV tables alone."""

import argparse
from collections import defaultdict
from itertools import product
import math
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

from plot_contract import BAR, CURVES, FOLD_COLORS, FOLDS, METHOD_COLORS, METHOD_LABELS
from plot_contract import METHODS, METRICS, OVERALL, PER_DOSE, WIDTHS, read_rows


def numeric(value):
    if value == '':
        return float('nan')
    number = float(value)
    if not math.isfinite(number):
        raise ValueError('non-finite plotted value')
    return number


def checked_index(rows, keys, expected, name):
    index = {tuple(row[key] for key in keys): row for row in rows}
    if len(index) != len(rows) or set(index) != set(expected):
        raise ValueError(f'{name}.csv has an incomplete or duplicate plot grid')
    return index


def load(data):
    keys = ('width', 'method', 'metric')
    overall = checked_index(read_rows(data, OVERALL), keys, product(WIDTHS, METHODS, METRICS), OVERALL)
    bars = checked_index(read_rows(data, BAR), keys, product(WIDTHS, METHODS, METRICS), BAR)
    rows = read_rows(data, PER_DOSE)
    folds = tuple(map(str, FOLDS))
    per_dose = checked_index(rows, (*keys, 'fold'), product(WIDTHS, METHODS, METRICS, folds), PER_DOSE)
    doses = {fold: {numeric(row['dose_um']) for row in rows if row['fold'] == fold} for fold in folds}
    if any(len(values) != 1 for values in doses.values()):
        raise ValueError('inconsistent held dose across variants')
    concentrations = {fold: values.pop() for fold, values in doses.items()}
    if len(set(concentrations.values())) != len(folds) or any(
        not math.isfinite(dose) or dose <= 0 for dose in concentrations.values()
    ):
        raise ValueError('expected three distinct positive held doses')

    learning = defaultdict(list)
    for row in read_rows(data, CURVES):
        key = (row['width'], row['method'], row['fold'])
        curve = learning[key]
        epoch = int(row['epoch'])
        best = int(row['best_validation'])
        if epoch != len(curve) + 1 or best not in (0, 1):
            raise ValueError('invalid learning curve coordinates')
        values = [numeric(row['train_loss']), numeric(row['validation_mse'])]
        if not all(map(math.isfinite, values)):
            raise ValueError('missing learning curve losses')
        curve.append((epoch, *values, best))
    if set(learning) != set(product(WIDTHS, METHODS, folds)) or any(
        sum(row[3] for row in curve) != 1 for curve in learning.values()
    ):
        raise ValueError('incomplete learning curves or best validation epochs')
    return overall, bars, per_dose, concentrations, learning


def save(figure, output, name):
    figure.savefig(output / f'{name}.png', dpi=160, facecolor='white')
    plt.close(figure)


def methods_axis(axis, *, show_labels):
    axis.set_xticks(
        range(len(METHODS)), METHOD_LABELS if show_labels else [''] * len(METHODS),
        rotation=23, ha='right', fontsize=8,
    )
    axis.set_xlim(-0.55, len(METHODS) - 0.45)
    axis.set_axisbelow(True)
    axis.spines[['top', 'right']].set_visible(False)


def metric_color(rows, metric):
    values = [numeric(row['mean']) for row in rows if row['metric'] == metric and row['mean']]
    low, high = (min(values), max(values)) if values else (0.0, 1.0)
    if low == high:
        margin = max(abs(low) * 0.001, 0.0001)
        low, high = low - margin, high + margin
    return plt.cm.ScalarMappable(norm=plt.Normalize(low, high), cmap='viridis')


def balloon(axis, x, y, point, color):
    value = numeric(point['mean'])
    if math.isfinite(value):
        diameter = numeric(point['diameter_pt'])
        if not math.isfinite(diameter) or diameter <= 0:
            raise ValueError('a defined balloon must have a positive diameter')
        axis.scatter(x, y, s=diameter ** 2, color=color.to_rgba(value),
                     edgecolor='#364152', linewidth=0.35, zorder=3)
        axis.annotate(f'{value:.3f}', (x, y), xytext=(0, diameter / 2 + 2),
                      textcoords='offset points', ha='center', va='bottom', fontsize=7)
    else:
        axis.text(x, y, 'NA', color='#636b74', ha='center', va='center', fontsize=8)


def balloon_axes(axis, ylabels, *, show_methods, show_levels):
    methods_axis(axis, show_labels=show_methods)
    axis.set_yticks(range(len(ylabels)), ylabels if show_levels else [''] * len(ylabels))
    axis.set_ylim(-0.48, len(ylabels) - 0.35)
    axis.tick_params(axis='y', labelsize=8)
    axis.grid(alpha=0.16)


def overall_balloon(output, overall):
    figure, axes = plt.subplots(1, 4, figsize=(15, 4.45), layout='constrained')
    figure.suptitle('tenv cold-drug-dose · TranSiGen width and dose encoding', fontsize=15)
    for axis, (metric, label) in zip(axes, METRICS.items(), strict=True):
        color = metric_color(overall.values(), metric)
        axis.set_title(label, fontsize=11)
        for y, width in enumerate(WIDTHS):
            for x, method in enumerate(METHODS):
                balloon(axis, x, y, overall[width, method, metric], color)
        balloon_axes(axis, WIDTHS, show_methods=True, show_levels=True)
        figure.colorbar(color, ax=axis, orientation='horizontal', shrink=0.75, aspect=24).ax.tick_params(
            labelsize=7,
        )
    axes[0].set_ylabel('Width')
    figure.supxlabel(
        'Color and labels show equal-fold means. Diameter is 21 + 4 × z points, limited to 12–30 points. '
        'z uses all fold runs.',
        fontsize=9,
    )
    save(figure, output, OVERALL)


def metric_bars(output, bars):
    figure, axes = plt.subplots(4, 4, figsize=(13.8, 10.2), layout='constrained', sharey='col')
    figure.suptitle('tenv cold-drug-dose · outer-fold scores', fontsize=15)
    for row, (width, panels) in enumerate(zip(WIDTHS, axes, strict=True)):
        panels[0].set_ylabel(width)
        for axis, (metric, label) in zip(panels, METRICS.items(), strict=True):
            if row == 0:
                axis.set_title(label, fontsize=11)
            for x, (method, color) in enumerate(zip(METHODS, METHOD_COLORS, strict=True)):
                scores = bars[width, method, metric]
                average, sd = numeric(scores['mean']), numeric(scores['sd'])
                if math.isfinite(average):
                    axis.bar(x, average, width=0.67, color=color, alpha=0.30)
                if math.isfinite(sd):
                    axis.errorbar(x, average, yerr=sd, fmt='none', color='#303030',
                                  capsize=3, lw=1.1, zorder=5)
                dots = [numeric(scores[f'fold_{fold}']) for fold in FOLDS]
                axis.scatter(x + np.array([-.16, 0, .16]), dots, color=color, s=20,
                             edgecolor='white', linewidth=0.3, zorder=4)
            axis.yaxis.grid(alpha=0.20)
            axis.tick_params(axis='y', labelsize=8, labelleft=True)
            methods_axis(axis, show_labels=row == len(WIDTHS) - 1)
    figure.supxlabel('Bars show equal-fold means ± sample SD. Each dot is one outer test fold.', fontsize=9)
    save(figure, output, BAR)


def held_dose_balloons(output, per_dose, concentrations):
    figure, axes = plt.subplots(4, 4, figsize=(14.2, 12.0), layout='constrained')
    figure.suptitle('tenv cold-drug-dose · scores at each held concentration', fontsize=15)
    ordered = sorted(concentrations, key=concentrations.__getitem__)
    labels = [f'{concentrations[fold]:g}' for fold in ordered]
    for column, (metric, label) in enumerate(METRICS.items()):
        color = metric_color(per_dose.values(), metric)
        axes[0, column].set_title(label, fontsize=11)
        for row, width in enumerate(WIDTHS):
            axis = axes[row, column]
            if column == 0:
                axis.set_ylabel(f'{width}\nDose (µM)')
            for y, fold in enumerate(ordered):
                for x, method in enumerate(METHODS):
                    balloon(axis, x, y, per_dose[width, method, metric, fold], color)
            balloon_axes(axis, labels, show_methods=row == len(WIDTHS) - 1, show_levels=column == 0)
        figure.colorbar(color, ax=axes[:, column], orientation='horizontal',
                        shrink=0.75, aspect=24).ax.tick_params(labelsize=7)
    figure.supxlabel('Color, labels and diameter show the mean across test profiles at each held dose.', fontsize=9)
    save(figure, output, PER_DOSE)


def learning_curves(output, learning):
    figure, axes = plt.subplots(8, 4, figsize=(14.0, 15.5), layout='constrained')
    figure.suptitle('tenv training and validation histories by fold', fontsize=15)
    for row, panels in enumerate(axes):
        width = WIDTHS[row // 2]
        position = 1 + row % 2
        for axis, method, label in zip(panels, METHODS, METHOD_LABELS, strict=True):
            for fold, color in zip(FOLDS, FOLD_COLORS, strict=True):
                curve = np.asarray(learning[width, method, str(fold)])
                axis.plot(curve[:, 0], curve[:, position], lw=1, color=color, label=f'Fold {fold}')
                if position == 2:
                    best = curve[curve[:, 3] == 1][0]
                    axis.scatter(best[0], best[position], s=14, color=color, zorder=3)
            if row == 0:
                axis.set_title(label, fontsize=10)
            if row == len(axes) - 1:
                axis.set_xlabel('Epoch', fontsize=8)
            axis.xaxis.set_major_locator(MaxNLocator(integer=True))
            axis.grid(alpha=0.20)
            axis.tick_params(labelsize=7)
            axis.spines[['top', 'right']].set_visible(False)
        label = 'Training loss' if position == 1 else 'Validation MSE'
        panels[0].set_ylabel(f'{width}\n{label}', fontsize=8)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, loc='outside upper right', ncols=3, frameon=False, fontsize=9)
    save(figure, output, CURVES)


def render(data, output):
    overall, bars, doses, concentrations, learning = load(data)
    output.mkdir(parents=True, exist_ok=True)
    overall_balloon(output, overall)
    metric_bars(output, bars)
    held_dose_balloons(output, doses, concentrations)
    learning_curves(output, learning)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument('--data', default=root / 'data', type=Path)
    parser.add_argument('--output', default=root, type=Path)
    args = parser.parse_args()
    render(args.data, args.output)


if __name__ == '__main__':
    main()
