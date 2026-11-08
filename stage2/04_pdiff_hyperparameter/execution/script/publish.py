#!/usr/bin/env python3
"""Render the Step 04 figures from the complete compact formal records."""
import argparse
import math
from pathlib import Path
import sys
import tempfile

import numpy as np

from evaluate import PR_FILENAME, PR_MODELS, PR_SCHEMA
from train import REPO, STEP03
from training_policy import MONITOR, verify_best_checkpoint

sys.path.insert(0, str(REPO / STEP03 / 'execution/script'))
from artifacts import SEEDS, read_json, rows, sha256, write_json
from metrics import METRICS, RECALL_GRID

import matplotlib

matplotlib.use('Agg')
from matplotlib import colors, pyplot as plt, ticker
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

DATASETS = {'tahoe': 'Tahoe', 'tenv': 'tenv'}
MODELS = ('perturbdiff', 'no_change', 'delta_mean')
FOLDS = range(3)
SLOTS = {(fold, seed) for fold in FOLDS for seed in SEEDS}
CHECKPOINT_FIELDS = ('identity', 'recipe_sha256', 'checkpoint_sha256', 'best_step',
                     'best_validation_loss', 'model_state_sha256')
TRAINING, VALIDATION, BEST = '#2166ac', '#24804b', '#bd252d'
REFERENCE_COLORS = {'no_change': '#a5aab1', 'delta_mean': '#636c77'}
METRIC_TITLES = ('Positive AP@100', 'Negative AP@100', 'Positive recall at P80',
                 'Negative recall at P80', 'PCC', 'Spearman', 'RMSE', 'R²')
FIGURES = ('learning_curve.png', 'metrics.png', 'pr_curves.png')
PR_COLORS = {'perturbdiff': (TRAINING, '#8bb9e6'), 'delta_mean': ('#606771', '#b8bfc8')}
PR_STYLES = {'tahoe': '-', 'tenv': '--'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify_dataset(root, dataset, hashes):
    def tracked(name, expected=None):
        path = root / name
        digest = sha256(path)
        require(expected is None or digest == expected, f'{path} hash mismatch')
        hashes[f'{dataset}/{name}'] = digest
        return path

    evaluation = read_json(tracked('evaluation.json'))
    require(evaluation['schema'] == 'step04-evaluation-v1'
            and evaluation['datasets'] == [dataset] and evaluation['stage'] == 'formal'
            and evaluation['complete'] is True
            and evaluation['tenv_test_memberships_previously_scored_in_step03'] == (dataset == 'tenv'),
            f'{dataset} needs a complete formal evaluation')
    plan_sha = evaluation['plan_sha256']
    plan = read_json(tracked('plan.json', plan_sha))
    inventory = read_json(tracked('outcomes.json', evaluation['outcomes_sha256']))
    require(plan['schema'] == 'step04-evaluation-plan-v1'
            and plan['dataset'] == dataset and plan['stage'] == 'formal'
            and inventory['schema'] == 'step04-evaluation-outcomes-v1'
            and inventory['plan_sha256'] == plan_sha,
            f'{dataset} frozen plan/outcomes mismatch')
    folds = plan['folds']
    require(len(folds) == len(FOLDS),
            f'{dataset} needs three frozen folds')
    outcomes = {(entry['fold'], entry['seed']): entry for entry in inventory['outcomes']}
    require(outcomes.keys() == SLOTS and len(inventory['outcomes']) == len(SLOTS),
            f'{dataset} needs nine unique formal outcomes')

    runs = []
    manifests = {}
    for (fold, seed), outcome in sorted(outcomes.items()):
        slot = f'slots/fold-{fold}/seed-{seed}'
        require(outcome['dataset'] == dataset and outcome['status'] == 'completed',
                f'{dataset}/{slot} has an incomplete outcome')
        prediction = read_json(tracked(f'{slot}/prediction.json', outcome['prediction_manifest_sha256']))
        run = read_json(tracked(f'{slot}/run.json', prediction['run_sha256']))
        replay = read_json(tracked(f'{slot}/replay.json', prediction['replay_sha256']))
        identity = run['identity']
        checkpoint = {key: run[key] for key in CHECKPOINT_FIELDS}
        require((identity['dataset'], identity['fold'], identity['seed'], identity['stage'])
                == (dataset, fold, seed, 'formal')
                and run['status'] == 'completed' and run['monitor'] == MONITOR
                and run['recipe_sha256'] == folds[fold]['recipe_sha256_by_seed'][str(seed)]
                and prediction['schema'] == 'step04-heldout-prediction-v1'
                and prediction['plan_sha256'] == plan_sha
                and all(prediction[key] == value for key, value in checkpoint.items()),
                f'{dataset}/{slot} formal run identity mismatch')
        require(replay == {**run['reload'], **checkpoint, 'plan_sha256': plan_sha}
                and replay['validation_loss'] == run['best_validation_loss']
                and replay['prediction_max_abs_error'] == 0
                and replay['prediction_shape'][1] == 18164,
                f'{dataset}/{slot} exact EMA replay mismatch')
        verify_best_checkpoint(run, f'{dataset}/{slot}')
        history = run['history']
        require(history[-1]['step'] <= run['optimizer_steps']
                and all(math.isfinite(event['training_loss']) and event['training_loss'] > 0
                        and event['validation_loss'] > 0 for event in history),
                f'{dataset}/{slot} invalid history')
        runs.append(run)
        manifests[fold, seed] = prediction

    per_run = rows(tracked('per_run.csv', evaluation['files']['per_run.csv']))
    scores = read_scores(dataset, per_run,
                         rows(tracked('summary.csv', evaluation['files']['summary.csv'])), manifests)
    pr_path = tracked(PR_FILENAME, evaluation['files'].get(f'{dataset}/{PR_FILENAME}'))
    pr = read_pr_curves(read_json(pr_path), dataset, evaluation, per_run)
    checkpoints = []
    if dataset == 'tenv':
        index = read_json(tracked('checkpoints.json'))
        require(index['schema'] == 'step04-tenv-checkpoints-v1'
                and index['plan_sha256'] == plan_sha
                and index['outcomes_sha256'] == hashes[f'{dataset}/outcomes.json']
                and index['evaluation_sha256'] == hashes[f'{dataset}/evaluation.json']
                and index['selection_sha256'] == plan['selection_sha256'],
                'tenv checkpoint release identity mismatch')
        checkpoints = verify_checkpoint_index(index['checkpoints'], folds, manifests)
    return runs, scores, pr, checkpoints


def read_pr_curves(source, dataset, evaluation, per_run):
    require(source['schema'] == PR_SCHEMA and source['dataset'] == dataset
            and source['plan_sha256'] == evaluation['plan_sha256']
            and source['outcomes_sha256'] == evaluation['outcomes_sha256']
            and source['per_run_sha256'] == evaluation['files']['per_run.csv']
            and np.array_equal(source['recall'], RECALL_GRID),
            f'{dataset} PR source does not match the formal evaluation')
    expected = {(row['model'], int(row['fold']), row['seed']): row for row in per_run
                if row['model'] in PR_MODELS}
    slots = {(row['model'], row['fold'], '' if row['seed'] is None else str(row['seed'])): row
             for row in source['runs']}
    require(slots.keys() == expected.keys() and len(source['runs']) == len(expected),
            f'{dataset} PR source run inventory mismatch')
    curves = {model: [[], []] for model in PR_MODELS}
    for slot, row in slots.items():
        require(row['dataset'] == dataset and type(row['fold']) is int
                and len(row['precision']) == 2, f'{dataset}/{slot} invalid PR run')
        for direction, (metric, values) in enumerate(zip(METRICS[:2], row['precision'])):
            if int(expected[slot][f'{metric}_valid_groups']) == 0:
                require(values is None, f'{dataset}/{slot}/{metric} PR curve must be undefined')
                continue
            points = np.asarray(values)
            require(points.shape == RECALL_GRID.shape and points.dtype.kind in 'iuf'
                    and np.all((0 <= points) & (points <= 1)),
                    f'{dataset}/{slot}/{metric} invalid PR precision')
            curves[slot[0]][direction].append(points)
    return curves


def read_scores(dataset, per_run, summary, manifests):
    expected = {('perturbdiff', fold, str(seed)) for fold, seed in SLOTS}
    expected |= {(model, fold, '') for fold in FOLDS for model in MODELS[1:]}
    require(len(per_run) == len(expected)
            and {(row['model'], int(row['fold']), row['seed']) for row in per_run} == expected
            and all(row['dataset'] == dataset and row['status'] == 'completed' for row in per_run),
            f'{dataset} needs nine scored runs and six scored references')
    for fold in FOLDS:
        group_counts = {int(row['group_count']) for row in per_run if int(row['fold']) == fold}
        expected_counts = {manifests[fold, seed]['group_count'] for seed in SEEDS}
        require(len(group_counts) == 1 and min(group_counts) > 0 and group_counts == expected_counts,
                f'{dataset} fold {fold} held-out group counts do not match')

    require(len(summary) == len(MODELS) and {row['model'] for row in summary} == set(MODELS),
            f'{dataset} summary model inventory mismatch')
    scores = {}
    for row in summary:
        model = row['model']
        runs = [run for run in per_run if run['model'] == model]
        require(row['dataset'] == dataset and row['complete'] == 'True'
                and int(row['expected_runs']) == int(row['completed_runs']) == len(runs),
                f'{dataset}/{model} incomplete summary')
        scores[model] = {}
        for metric in METRICS:
            count = int(row[f'{metric}_valid_runs'])
            group_counts = [int(run[f'{metric}_valid_groups']) for run in runs]
            height = float(row[f'{metric}_mean']) if row[f'{metric}_mean'] else None
            require(all(0 <= n <= int(run['group_count']) and bool(n) == bool(run[metric])
                        for n, run in zip(group_counts, runs))
                    and count == sum(n > 0 for n in group_counts)
                    and ((height is None and count == 0)
                         or (height is not None and math.isfinite(height) and count > 0)),
                    f'{dataset}/{model}/{metric} invalid value or valid count')
            scores[model][metric] = {
                'height': height, 'status': 'undefined' if height is None else 'defined',
                'valid_runs': count, 'expected_runs': len(runs), 'valid_groups': sum(group_counts),
                'total_groups': sum(int(run['group_count']) for run in runs),
            }
    return scores


def verify_checkpoint_index(entries, folds, manifests):
    by_slot = {(row['identity']['fold'], row['identity']['seed']): row for row in entries}
    require(by_slot.keys() == SLOTS and len(entries) == len(SLOTS),
            'tenv release needs nine unique checkpoint identities')
    for (fold, seed), row in by_slot.items():
        manifest = manifests[fold, seed]
        fields = (*CHECKPOINT_FIELDS, 'run_sha256', 'replay_sha256')
        require(all(row[key] == manifest[key] for key in fields)
                and row['directory'] == f'fold-{fold}/seed-{seed}'
                and row['request_sha256'] == folds[fold]['request_sha256']
                and row['policy_sha256'] == folds[fold]['policy_sha256']
                and row['raw_recipe_sha256'] == folds[fold]['recipes'][str(seed)]['sha256'],
                f'tenv checkpoint index mismatch for fold {fold}, seed {seed}')
    return [by_slot[slot] for slot in sorted(by_slot)]


def verify_local_release(root, formal_root, index_path, checkpoints):
    require(sha256(root / 'checkpoints.json') == sha256(index_path),
            'local tenv checkpoint index differs from the recorded index')
    evidence = {}
    for row in checkpoints:
        slot = root / row['directory']
        for filename, field in [('run.json', 'run_sha256'), ('replay.json', 'replay_sha256'),
                                ('request.json', 'request_sha256'), ('recipe.json', 'raw_recipe_sha256')]:
            require(sha256(slot / filename) == row[field], f'local {slot / filename} hash mismatch')
        fold, seed = row['identity']['fold'], row['identity']['seed']
        fitted = formal_root / 'slots' / f'fold-{fold}-seed-{seed}' / 'fit/best.pt'
        checkpoint = slot / 'best.pt'
        require(checkpoint.is_file() and not checkpoint.is_symlink()
                and checkpoint.samefile(fitted) and checkpoint.stat().st_size > 0,
                f'{checkpoint} must hard-link to the verified formal fit')
        require(sha256(checkpoint) == row['checkpoint_sha256'], f'local {checkpoint} hash mismatch')
        evidence[row['directory']] = {
            'bytes': checkpoint.stat().st_size, 'same_file_as_formal_fit': True,
            'metadata_hashes_verified': True, 'checkpoint_hash_verified': True,
        }
    return evidence


def setup_axis(ax):
    ax.set_axisbelow(True)
    ax.grid(axis='y', color='#e4e8ed', linewidth=0.8)
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.spines['bottom'].set_color('#b0b7bf')
    ax.tick_params(length=0, pad=6)


def draw_curves(runs, destination):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharex=True, sharey=True)
    fig.subplots_adjust(left=.10, right=.975, top=.75, bottom=.18, wspace=.12)
    fig.suptitle('PerturbDiff formal learning curves', x=.10, ha='left', y=.96, fontweight='bold', fontsize=18)
    fig.text(.10, .875, 'Nine runs per dataset  •  Losses recorded at validation checks', fontsize=10)
    report = {}
    for ax, (dataset, label) in zip(axes, DATASETS.items()):
        setup_axis(ax)
        ax.set_title(label, loc='left', fontweight='bold', pad=10)
        ax.set_yscale('log')
        ax.set_yticks([.4, .6, 1, 2, 4, 8], labels=['0.4', '0.6', '1', '2', '4', '8'])
        ax.yaxis.set_minor_locator(ticker.NullLocator())
        ax.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
        ax.set_xlabel('Optimizer step')
        panel = []
        for run in runs[dataset]:
            history = run['history']
            fold, seed = run['identity']['fold'], run['identity']['seed']
            steps = [event['step'] for event in history]
            row = {'fold': fold, 'seed': seed, 'stop_step': run['optimizer_steps']}
            for name, color in [('training', TRAINING), ('validation', VALIDATION)]:
                line, = ax.plot(steps, [event[f'{name}_loss'] for event in history], color=color,
                                linestyle=('-', '--', ':')[fold], alpha=.66, linewidth=1.25)
                row[name] = {'color': colors.to_hex(line.get_color()),
                             'points': line.get_xydata().tolist()}
            marker, = ax.plot(run['best_step'], run['best_validation_loss'], color=BEST, marker='o',
                              markersize=4.2, markeredgecolor='white', markeredgewidth=.35, zorder=4)
            row['best'] = {'color': colors.to_hex(marker.get_color()),
                           'point': marker.get_xydata()[0].tolist()}
            panel.append(row)
        report[dataset] = panel
    axes[0].set_ylabel('Loss (log scale)')
    fig.legend(handles=[Line2D([], [], color=TRAINING, label='Training'),
                        Line2D([], [], color=VALIDATION, label='EMA validation'),
                        Line2D([], [], color=BEST, marker='o', linestyle='', label='Selected checkpoint')],
               loc='lower center', bbox_to_anchor=(.54, .026), ncols=3, frameon=False)
    fig.savefig(destination, dpi=160)
    plt.close(fig)
    return report


def draw_metrics(scores, destination):
    fig, axes = plt.subplots(2, 4, figsize=(15, 7.1))
    fig.subplots_adjust(left=.06, right=.984, top=.77, bottom=.135, wspace=.42, hspace=.48)
    fig.suptitle('PerturbDiff held-out metrics', x=.06, ha='left', y=.963, fontweight='bold', fontsize=18)
    fig.text(.06, .89, 'Mean across nine model runs or three reference folds per dataset', fontsize=10)
    fig.legend(handles=[
        Patch(facecolor=TRAINING, label='Tahoe PerturbDiff'),
        Patch(facecolor=VALIDATION, label='tenv PerturbDiff'),
        *[Patch(facecolor=color, label=model) for model, color in REFERENCE_COLORS.items()],
    ], loc='upper left', bbox_to_anchor=(.053, .864), frameon=False, ncols=4)
    fig.text(.06, .047, 'Dashes mark undefined scores (0 of 3 folds). Zero scores are labeled 0.000. '
             'Tenv reuses Step 03 held-out memberships.', fontsize=9)
    report = {}
    for ax, metric, title in zip(axes.flat, METRICS, METRIC_TITLES):
        setup_axis(ax)
        ax.set_title(title, loc='left', fontweight='bold', pad=10, fontsize=11)
        ax.axhline(0, linewidth=.8, color='#8d96a0')
        ax.set_xlim(-.68, 2.28)
        ax.set_xticks([0, 1.6], DATASETS.values())
        ax.yaxis.set_major_locator(ticker.MaxNLocator(5))
        ax.margins(y=.19)
        panel = []
        for dataset, center, model_color in zip(DATASETS, (0, 1.6), (TRAINING, VALIDATION)):
            for offset, model in zip((-.43, 0, .43), MODELS):
                value = scores[dataset][model][metric]
                height = value['height']
                color = REFERENCE_COLORS.get(model, model_color)
                x = center + offset
                if height is None:
                    mark = ax.annotate('—', (x, 0), xytext=(0, 3), textcoords='offset points',
                                       color=color, ha='center', va='bottom', fontsize=11)
                    color = colors.to_hex(mark.get_color())
                    x = mark.xy[0]
                else:
                    bar, = ax.bar(x, height, width=.36, color=color, zorder=2)
                    x = bar.get_x() + bar.get_width() / 2
                    height = bar.get_height()
                    color = colors.to_hex(bar.get_facecolor())
                    ax.annotate(f'{height:.3f}', (x, height), xytext=(0, 3 if height >= 0 else -3),
                                textcoords='offset points', ha='center',
                                va='bottom' if height >= 0 else 'top', fontsize=7)
                panel.append({'dataset': dataset, 'model': model, **value,
                              'x': x, 'height': height, 'color': color})
        report[metric] = panel
    fig.savefig(destination, dpi=160)
    plt.close(fig)
    return report


def draw_pr_curves(pr, destination, hashes):
    layers = [('tahoe', 'perturbdiff', 8), ('tenv', 'perturbdiff', 6),
              ('tahoe', 'delta_mean', 4), ('tenv', 'delta_mean', 4)]
    report = {'recall': RECALL_GRID.tolist()}
    dark = '#18212c'
    muted = '#45525f'
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.85), sharey=True, sharex=True)
    fig.subplots_adjust(left=.078, right=.974, top=.702, bottom=.19, wspace=.12)
    fig.suptitle('Held-out precision–recall curves', x=.078, y=.966, ha='left',
                 fontsize=20, fontweight='bold', color=dark)
    subtitle = 'PerturbDiff and delta_mean  •  Top 100 genes per direction'
    fig.text(.078, .861, subtitle, fontsize=10.5, color=muted)
    model_handles = [
        Line2D([], [], color=PR_COLORS['perturbdiff'][0], linewidth=2.5, label='PerturbDiff'),
        Line2D([], [], color=PR_COLORS['delta_mean'][0], linewidth=2.5, label='delta_mean'),
    ]
    dataset_handles = [
        Line2D([], [], color=dark, linewidth=1.8, linestyle=PR_STYLES[dataset], label=label)
        for dataset, label in DATASETS.items()
    ]
    for title, handles, anchor in [('Model', model_handles, .072),
                                   ('Dataset', dataset_handles, .36)]:
        legend = fig.legend(handles=handles, title=title, loc='upper left',
                            bbox_to_anchor=(anchor, .85), ncol=2, frameon=False,
                            fontsize=9, title_fontsize=9, handlelength=2.5,
                            columnspacing=1.5, alignment='left')
        legend.get_title().set_fontweight('bold')

    for direction, (ax, name) in enumerate(zip(axes, ('positive', 'negative'))):
        title = f'PR {name.capitalize()}'
        panel = []
        ax.set_title(title, loc='left', fontweight='bold', fontsize=12, color=dark, pad=9)
        ax.set(xlim=(0, 1), ylim=(-.015, 1.025), xlabel='Recall')
        ax.set_xticks(np.linspace(0, 1, 6))
        ax.set_yticks(np.linspace(0, 1, 6))
        ax.set_axisbelow(True)
        ax.grid(axis='y', color='#e4e8ed')
        ax.tick_params(colors=muted, length=3)
        ax.spines[['top', 'right']].set_visible(False)
        ax.spines[['bottom', 'left']].set_color('#b0b7bf')
        for dataset, model, layer in reversed(layers):
            curves = pr[dataset][model][direction]
            if not curves:
                continue
            color, fill = PR_COLORS[model]
            style = PR_STYLES[dataset]
            low, high = np.min(curves, axis=0), np.max(curves, axis=0)
            median = np.median(curves, axis=0)
            ax.fill_between(RECALL_GRID, low, high, step='pre', color=fill, alpha=.34,
                            linewidth=0, zorder=layer)
            ax.step(RECALL_GRID, median, where='pre', color=color, linestyle=style, linewidth=1.9,
                    zorder=layer + 1)
            panel.append({'dataset': dataset, 'model': model, 'color': color,
                          'linestyle': style, 'valid_runs': len(curves),
                          'min': low.tolist(), 'median': median.tolist(), 'max': high.tolist()})
        report[name] = list(reversed(panel))

    axes[0].set_ylabel('Precision')
    count_text = 'Per dataset: 9 PerturbDiff runs and 3 delta_mean folds.'
    fig.text(.078, .082, 'Dark lines show pointwise medians. Shaded bands span pointwise minima and maxima.',
             color=muted, fontsize=9)
    fig.text(.078, .043, count_text, color=muted, fontsize=9)
    description = (f'{subtitle}. {count_text} Blue PerturbDiff; gray delta_mean. Solid Tahoe; dashed tenv. '
                   'Pointwise median and min–max of group-macro-averaged run curves. ')
    evidence = [f'{dataset} evaluation SHA-256 {hashes[f"{dataset}/evaluation.json"]}' for dataset in DATASETS]
    fig.savefig(destination, dpi=180, facecolor='white', metadata={
        'Title': 'Step 04 held-out precision–recall curves',
        'Description': description + ' '.join(evidence),
    })
    plt.close(fig)
    return report


def publish(record, output, evidence, local):
    hashes = {}
    runs, scores, pr = {}, {}, {}
    for dataset in DATASETS:
        runs[dataset], scores[dataset], pr[dataset], checkpoints = verify_dataset(record / dataset, dataset, hashes)
    require(len({run['identity']['gene_sha256'] for panel in runs.values() for run in panel}) == 1,
            'formal runs require one shared gene axis')
    report = {'schema': 'step04-publication-v1', 'matplotlib': matplotlib.__version__,
              'renderer_sha256': sha256(__file__), 'inputs_sha256': hashes,
              'checkpoint_index': checkpoints}
    if local:
        report['local_checkpoint_release'] = verify_local_release(
            *local, record / 'tenv/checkpoints.json', checkpoints,
        )

    outputs = [output / name for name in FIGURES]
    destinations = (*outputs, evidence)
    require(len({path.resolve() for path in destinations}) == len(destinations)
            and not any(path.exists() or path.is_symlink() for path in destinations),
            'publication requires fresh figure and evidence paths')
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as temporary:
        staging = Path(temporary)
        figures = [staging / name for name in FIGURES]
        curve, metrics, pr_figure = figures
        with plt.rc_context({'font.family': 'DejaVu Sans', 'text.color': '#18212c',
                             'axes.labelcolor': '#18212c', 'axes.titlecolor': '#18212c'}):
            report['learning_curve'] = draw_curves(runs, curve)
            report['metrics'] = draw_metrics(scores, metrics)
        report['pr_curves'] = draw_pr_curves(pr, pr_figure, hashes)
        report['figures_sha256'] = {path.name: sha256(path) for path in figures}
        output.mkdir(parents=True, exist_ok=True)
        for path in figures:
            path.rename(output / path.name)
        write_json(evidence, report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('formal_record', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('evidence', type=Path)
    parser.add_argument('--local-checkpoints', type=Path, nargs=2, metavar=('RELEASE', 'FORMAL_ROOT'),
                        help='verify local metadata, checkpoint hashes and hard links')
    args = parser.parse_args()
    try:
        publish(args.formal_record, args.output, args.evidence, args.local_checkpoints)
    except (KeyError, OSError, TypeError, ValueError) as error:
        parser.exit(1, f'publication rejected: {error}\n')


if __name__ == '__main__':
    main()
