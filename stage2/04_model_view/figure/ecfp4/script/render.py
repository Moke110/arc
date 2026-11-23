"""Render the accepted 02 CSV values with canonical ECFP4 labels."""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / 'util/output/api'))
from stage2_util import palette

SEALED_MLP_LABEL = 'MLP_Morgan'
DISPLAY_MLP_LABEL = 'ECFP4_MLP'
HUES = dict(PRnet='blue', TranSiGen='orange', ECFP4_MLP='green',
            Mean='gray', Mean_cell='purple', Mean_drug='magenta')
SCENARIOS = ('cold-drug', 'warm-start')
METRICS = ('pcc', 'spearman', 'rmse', 'r2')


def read_values(path, labels):
    values = pd.read_csv(path)
    for column in labels:
        values[column] = values[column].replace(SEALED_MLP_LABEL, DISPLAY_MLP_LABEL)
    return values


def draw_metrics(values, output):
    figure, axes = plt.subplots(2, 4, figsize=(17, 7), constrained_layout=True)
    for scenario, row_axes in zip(SCENARIOS, axes, strict=True):
        for metric, axis in zip(METRICS, row_axes, strict=True):
            rows = values.loc[(values.scenario == scenario) & (values.metric == metric)].set_index('method')
            methods = [method for method in HUES if method in rows.index]
            axis.bar(range(len(methods)), [rows.loc[method, 'value'] for method in methods],
                     color=[palette.color(HUES[method], 6) for method in methods],
                     edgecolor=[palette.color(HUES[method], 2) for method in methods])
            labels = []
            for index, method in enumerate(methods):
                row = rows.loc[method]
                if bool(row.applicable):
                    status = (f'n={int(row.n_parents)} u={int(row.n_undefined)} '
                              f't={int(row.n_constant_truth)} p={int(row.n_constant_prediction)}')
                else:
                    status = 'N/A'
                    axis.text(index, 0, f'N/A\n{row.reason}', rotation=90, ha='center',
                              va='bottom', fontsize=6)
                labels.append(f'{method}\n{status}')
            axis.set_xticks(range(len(methods)), labels, rotation=50, ha='right', fontsize=7)
            axis.set_title(f'{scenario} · Δ {metric.upper()}')
            axis.axhline(float(rows.zero.iloc[0]), color=palette.color('gray', 4), lw=.6)
            axis.grid(axis='y', alpha=.15)
    figure.suptitle('Equal-weight parent metrics; n=defined parents, u=undefined; '
                    't=constant truth parents, p=constant prediction seed-parent vectors. IDs/reasons in CSV')
    figure.savefig(output / 'metrics.png')
    plt.close(figure)


def draw_paired_differences(pairs, output):
    figure, axes = plt.subplots(2, 4, figsize=(19, 10), constrained_layout=True)
    for scenario, row_axes in zip(SCENARIOS, axes, strict=True):
        for metric, axis in zip(METRICS, row_axes, strict=True):
            rows = pairs.loc[(pairs.scenario == scenario) & (pairs.metric == metric)]
            labels = []
            for index, row in enumerate(rows.itertuples()):
                color = palette.color(HUES[row.candidate], 2)
                axis.plot([row.low, row.high], [index, index], color=color, lw=2)
                axis.plot(row.estimate, index, 'o', color=color, ms=4)
                counts = f'{int(row.candidate_constant_prediction)}/{int(row.reference_constant_prediction)}'
                labels.append(f'{row.candidate} − {row.reference} [pred const {counts}]')
            axis.set_yticks(range(len(rows)), labels, fontsize=7)
            excluded = 0
            if not rows.empty:
                axis.axvline(float(rows.zero.iloc[0]), color=palette.color('gray', 4), lw=.7)
                excluded = int(rows.n_excluded_constant.iloc[0])
                unavailable = rows.unavailable_reference.iloc[0]
                if pd.notna(unavailable) and str(unavailable):
                    axis.set_xlabel(f'N/A: {unavailable}', fontsize=7)
            axis.set_title(f'{scenario} · Δ {metric.upper()}; excluded={excluded}', fontsize=9)
            axis.grid(axis='x', alpha=.15)
    figure.suptitle('Paired drug bootstrap: candidate minus reference, 95% CI (errors: negative is better)')
    figure.savefig(output / 'paired_differences.png')
    plt.close(figure)


def draw_learning_curves(curves, output):
    figure, axes = plt.subplots(2, 3, figsize=(17, 9), constrained_layout=True)
    for scenario, row_axes in zip(SCENARIOS, axes, strict=True):
        methods = ('PRnet', 'TranSiGen', DISPLAY_MLP_LABEL)
        for method, training_axis in zip(methods, row_axes, strict=True):
            validation_axis = training_axis.twinx()
            validation_axis.spines['right'].set_visible(True)
            rows = curves.loc[(curves.scenario == scenario) & (curves.method == method)]
            for index, (_, run) in enumerate(rows.groupby(['fold', 'seed'])):
                # Every recorded epoch is plotted; shades identify the nine fold/seed runs.
                validation_axis.plot(run.epoch, run.validation_mse,
                                     color=palette.color('green', index % 7), alpha=.65, lw=.8)
                training_axis.plot(run.epoch, run.train_loss, color=palette.color('blue', index % 7),
                                   alpha=.35, lw=.6, ls='--')
            training_axis.set_title(f'{scenario} · {method}')
            training_axis.set_xlabel('epoch')
            training_axis.set_ylabel('native train loss (dashed)', color=palette.color('blue', 2))
            validation_axis.set_ylabel('validation MSE (solid)', color=palette.color('green', 2))
            training_axis.tick_params(axis='y', labelcolor=palette.color('blue', 2))
            validation_axis.tick_params(axis='y', labelcolor=palette.color('green', 2))
            training_axis.grid(alpha=.15)
    figure.suptitle('Independent y-scales: native training loss (blue, left) and validation MSE '
                    '(green, right). Native losses are not comparable across models')
    figure.savefig(output / 'learning_curve.png')
    plt.close(figure)


def render(data, output):
    data, output = Path(data), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                         'savefig.dpi': 140})
    draw_metrics(read_values(data / 'metrics_values.csv', ('method',)), output)
    draw_paired_differences(read_values(data / 'paired_differences_values.csv',
                                        ('candidate', 'reference')), output)
    draw_learning_curves(read_values(data / 'learning_curve_values.csv', ('method',)), output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--figure-data', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    render(args.figure_data, args.out_dir)
