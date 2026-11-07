"""Render final display figures from the verified scorer tables and curves."""

import numpy as np

from artifacts import rows
from plots import METRICS, plt


MODEL_NAMES = {
    'no_change': 'No Change', 'delta_mean': 'Delta Mean', 'prnet': 'PRnet',
    'cpa': 'CPA', 'chemcpa': 'chemCPA', 'cmonge': 'CMonge', 'perturbdiff': 'PerturbDiff',
}


def model_name(model):
    return MODEL_NAMES.get(model, model)


def save(figure, destination, name):
    for extension in ('png', 'svg'):
        figure.savefig(destination / f'{name}.{extension}', dpi=150)
    plt.close(figure)


def render_metrics(output, destination, summary):
    runs = rows(output / 'per_run.csv')
    labels = [model_name(row['model']) for row in summary]
    figure, axes = plt.subplots(2, 4, figsize=(24, 12), layout='constrained')
    for metric, axis in zip(METRICS, axes.flat):
        for x, row in enumerate(summary):
            mean = float(row[f'{metric}_mean'])
            sd = float(row[f'{metric}_sample_sd'])
            values = np.array([float(run[metric]) for run in runs if run['model'] == row['model']])
            if np.isfinite(mean):
                axis.bar(x, mean, color='lightsteelblue', width=.7,
                         yerr=sd if np.isfinite(sd) else None, capsize=4)
                offsets = np.linspace(-.13, .13, len(values))
                finite = np.isfinite(values)
                axis.scatter(x + offsets[finite], values[finite], s=24, color='black', zorder=3)
            else:
                axis.text(x, .48, 'Undefined', rotation=90, ha='center', va='center',
                          transform=axis.get_xaxis_transform(), color='.35')
        axis.set_title(metric)
        axis.set_xticks(range(len(labels)), labels, rotation=42, ha='right', rotation_mode='anchor')
        axis.set_xlim(-.5, len(labels) - .5)
    figure.suptitle('Fold × seed scores  ·  finite-run mean ± sample SD\n'
                    'Dots show individual runs. Reference seeds are deterministic copies.')
    save(figure, destination, 'metrics')


def render_pr_panels(output, destination, summary):
    bands = rows(output / 'pr_bands.csv')
    figure, axes = plt.subplots(2, 4, figsize=(24, 12), layout='constrained')
    for row, axis in zip(summary, axes.flat):
        for direction, color in [('Positive', '#c62828'), ('Negative', '#1565c0')]:
            selected = [band for band in bands
                        if band['model'] == row['model'] and band['direction'] == direction.lower()]
            counts = [int(band['valid_runs']) for band in selected]
            if not counts or max(counts) == 0:
                axis.plot([], [], color=color, label=f'{direction} · undefined')
                continue
            recall, mean, low, high = ([float(band[key]) for band in selected]
                                       for key in ('recall', 'mean', 'min', 'max'))
            defined = str(min(counts)) if min(counts) == max(counts) else f'{min(counts)}–{max(counts)}'
            axis.plot(recall, mean, color=color, label=f'{direction} · {defined} runs')
            axis.fill_between(recall, low, high, color=color, alpha=.15)
        axis.set(xlabel='Recall', ylabel='Precision', title=model_name(row['model']),
                 xlim=(0, 1), ylim=(0, 1.01))
        axis.legend(loc='upper right', fontsize=16)
    for axis in list(axes.flat)[len(summary):]:
        axis.set_axis_off()
    figure.suptitle('Directional precision–recall  ·  run mean and min–max band\n'
                    'Groups are macro-averaged within each run at tied-score thresholds.')
    save(figure, destination, 'pr_panels')


def render_best_run(output, destination, summary):
    summaries = {row['model']: row for row in summary}
    figure, axis = plt.subplots(figsize=(16, 9), layout='constrained')
    with np.load(output / 'pr_sources.npz', allow_pickle=False) as source:
        identifiers = source['identifiers'].tolist()
        for selected in rows(output / 'best_runs.csv'):
            model = selected['model']
            row = summaries[model]
            name = model_name(model)
            completed = f"{row['completed_runs']}/{row['expected_runs']} completed"
            reference = ' · deterministic copies' if row['deterministic_reference'] == 'True' else ''
            if not selected['fold']:
                axis.plot([], [], label=f'{name} · undefined\n{completed}{reference}')
                continue
            key = f"{model}:{selected['fold']}:{selected['seed']}"
            curves = source['curves'][identifiers.index(key)]
            axis.plot(source['recall'], (curves[0] + curves[1]) / 2,
                      label=f"{name} · fold {selected['fold']}, seed {selected['seed']}\n{completed}{reference}")
    axis.set(xlabel='Recall', ylabel='Mean directional precision', xlim=(0, 1), ylim=(0, 1.01))
    figure.suptitle('Post-hoc best test run per model\n'
                    'Models may select different folds. The nine-run summaries are primary.')
    axis.legend(loc='center left', bbox_to_anchor=(1.01, .5), fontsize=15, labelspacing=.85)
    save(figure, destination, 'best_run')


def render(output, destination):
    summary = rows(output / 'summary.csv')
    with plt.rc_context({'font.size': 18, 'axes.titlesize': 21, 'axes.labelsize': 18,
                         'xtick.labelsize': 15, 'ytick.labelsize': 18, 'figure.titlesize': 22}):
        render_metrics(output, destination, summary)
        render_pr_panels(output, destination, summary)
        render_best_run(output, destination, summary)
