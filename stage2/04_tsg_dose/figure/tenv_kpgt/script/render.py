"""Reproduce the tenv KPGT figures using the published CSV tables."""

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

SCENARIOS = ('cold-drug', 'cold-drug-dose')
METHODS = (
    'no-dose',
    'scaling',
    'input-numeric',
    'latent-numeric',
    'input-embedding',
    'latent-embedding',
    'film',
    'KPGT_MLP',
    'Mean_cell',
)
LABELS = (
    'No dose',
    'Scaling',
    'Input numeric',
    'Latent numeric',
    'Input embedding',
    'Latent embedding',
    'FiLM',
    'KPGT_MLP',
    'Mean_cell',
)
COLORS = (
    '#386095',
    '#397e88',
    '#61943e',
    '#ad7745',
    '#5280be',
    '#8970b4',
    '#a65170',
    '#45535c',
    '#8b9297',
)
METRICS = {
    'pcc': 'Δ Pearson r',
    'spearman': 'Δ Spearman ρ',
    'rmse': 'Δ RMSE',
    'r2': 'Δ R²',
}
FIGURES = ('metrics.png', 'metrics_cold_drug_dose.png', 'learning_curves.png')


def read(data, name):
    with (data / f'{name}.csv').open(newline='') as stream:
        return list(csv.DictReader(stream))


def index(rows, fields, expected, table):
    result = {tuple(row[field] for field in fields): row for row in rows}
    if len(result) != len(rows) or set(result) != expected:
        raise ValueError(f'{table}.csv has an incomplete or duplicate plot grid')
    return result


def numeric(value):
    if value == '':
        return float('nan')
    number = float(value)
    if not math.isfinite(number):
        raise ValueError('non-finite plotted value')
    return number


def load(data):
    keys = {(scenario, method) for scenario in SCENARIOS for method in METHODS}
    folds = ('0', '1', '2')
    pooled = index(read(data, 'pooled'), ('scenario', 'method'), keys, 'pooled')
    dots = index(
        read(data, 'folds'),
        ('scenario', 'method', 'fold'),
        {(*key, fold) for key in keys for fold in folds},
        'folds',
    )
    coverage = index(
        read(data, 'coverage'),
        ('scenario', 'fold'),
        {(scenario, fold) for scenario in SCENARIOS for fold in folds},
        'coverage',
    )
    learning = {}
    for row in read(data, 'learning'):
        key = (row['scenario'], row['method'], row['fold'])
        curve = learning.setdefault(key, [])
        if int(row['epoch']) != len(curve) + 1:
            raise ValueError('learning.csv has invalid epoch ordering')
        curve.append(
            (
                int(row['epoch']),
                numeric(row['train_loss']),
                numeric(row['validation_mse']),
            )
        )
    if set(learning) != {key for key in dots if key[1] != 'Mean_cell'}:
        raise ValueError('learning.csv has an incomplete curve grid')
    return pooled, dots, coverage, learning


def save(figure, output, name):
    figure.savefig(output / name, dpi=160, facecolor='white')
    plt.close(figure)


def methods_axis(axis, *, labels):
    axis.set_xticks(
        range(len(METHODS)),
        LABELS if labels else [],
        rotation=32,
        ha='right',
        fontsize=8,
    )
    axis.set_xlim(-0.55, len(METHODS) - 0.45)
    axis.yaxis.grid(alpha=0.2)
    axis.set_axisbelow(True)
    axis.spines[['top', 'right']].set_visible(False)


def overall(output, pooled, dots):
    figure, axes = plt.subplots(2, 4, figsize=(17, 7.4), layout='constrained')
    figure.suptitle('tenv dose methods with frozen KPGT', fontsize=16)
    for use_labels, scenario, panels in zip(
        (False, True), SCENARIOS, axes, strict=True
    ):
        for axis, (metric, label) in zip(panels, METRICS.items(), strict=True):
            for pos, (method, color) in enumerate(zip(METHODS, COLORS, strict=True)):
                for fold, dx in zip(('0', '1', '2'), (-0.19, 0.0, 0.19), strict=True):
                    axis.scatter(
                        pos + dx,
                        numeric(dots[scenario, method, fold][metric]),
                        marker='o',
                        s=16,
                        color=color,
                        alpha=0.6,
                    )
                axis.scatter(
                    pos,
                    numeric(pooled[scenario, method][metric]),
                    marker='D',
                    s=27,
                    color=color,
                    edgecolor='white',
                    linewidth=0.5,
                )
            axis.set_title(f'{scenario}  ·  {label}', fontsize=10)
            methods_axis(axis, labels=use_labels)
    figure.supxlabel(
        'Diamonds show pooled sample-equal scores. Circles show individual folds.',
        fontsize=10,
    )
    save(figure, output, FIGURES[0])


def held_doses(output, dots, coverage):
    figure, axes = plt.subplots(
        4, 3, figsize=(15, 10), layout='constrained', sharey='row'
    )
    figure.suptitle('tenv cold-drug-dose test by held concentration', fontsize=16)
    scenario = 'cold-drug-dose'
    for level, ((metric, label), panels) in enumerate(
        zip(METRICS.items(), axes, strict=True)
    ):
        panels[0].set_ylabel(label)
        for fold, axis in enumerate(panels):
            meta = coverage[scenario, str(fold)]
            values = [
                numeric(dots[scenario, method, str(fold)][metric]) for method in METHODS
            ]
            axis.scatter(range(len(METHODS)), values, color=COLORS, s=32)
            methods_axis(axis, labels=level == len(METRICS) - 1)
            if level == 0:
                axis.set_title(
                    f'{float(meta["held_dose_um"]):g} µM  ·  n = {meta["n_test_samples"]}'
                )
    save(figure, output, FIGURES[1])


def curves(output, learning):
    figure, axes = plt.subplots(4, 8, figsize=(22, 9), layout='constrained')
    figure.suptitle('tenv training and validation histories by fold', fontsize=16)
    fold_colors = ('#23669d', '#c16041', '#43866a')
    for row, panels in enumerate(axes):
        scenario = SCENARIOS[row // 2]
        position = 1 + row % 2
        for method, label, axis in zip(METHODS[:-1], LABELS[:-1], panels, strict=True):
            for fold, color in enumerate(fold_colors):
                curve = np.asarray(learning[scenario, method, str(fold)])
                axis.plot(
                    curve[:, 0],
                    curve[:, position],
                    color=color,
                    linewidth=1,
                    label=f'Fold {fold}',
                )
                if position == 2:
                    best = curve[np.argmin(curve[:, position])]
                    axis.scatter(best[0], best[position], s=14, color=color, zorder=3)
            axis.grid(alpha=0.2)
            axis.spines[['top', 'right']].set_visible(False)
            if row == 0:
                axis.set_title(label, fontsize=9)
            if row == len(axes) - 1:
                axis.set_xlabel('Epoch', fontsize=8)
            axis.tick_params(labelsize=7)
        label = 'Training loss' if position == 1 else 'Validation MSE'
        panels[0].set_ylabel(f'{scenario}\n{label}', fontsize=9)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles, labels, loc='outside upper right', ncols=3, frameon=False, fontsize=9
    )
    save(figure, output, FIGURES[2])


def render(data, output):
    pooled, dots, coverage, learning = load(data)
    output.mkdir(exist_ok=True)
    overall(output, pooled, dots)
    held_doses(output, dots, coverage)
    curves(output, learning)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument('--data', type=Path, default=root / 'data')
    parser.add_argument('--output', type=Path, default=root)
    args = parser.parse_args()
    render(args.data, args.output)


if __name__ == '__main__':
    main()
