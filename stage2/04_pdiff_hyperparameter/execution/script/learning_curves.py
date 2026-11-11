"""Replay original Step 04 loss histories without loading retrieval scores."""
import argparse
import math
from pathlib import Path
import sys
import tempfile
REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / '.git').exists())
STEP = REPO / 'stage2/04_pdiff_hyperparameter'
sys.path.insert(0, str(REPO / 'stage2/03_modes_tenv/execution/script'))
from artifacts import SEEDS, read_json, sha256, write_json
from training_policy import MONITOR, verify_best_checkpoint
from checkpoint_release import CHECKPOINT_FIELDS, SLOTS, require, verify_checkpoint_index
import matplotlib
matplotlib.use('Agg')
from matplotlib import colors, pyplot as plt, ticker
from matplotlib.lines import Line2D
DATASETS = {'tahoe': 'Tahoe', 'tenv': 'tenv'}
TRAINING, VALIDATION, BEST = '#2166ac', '#24804b', '#bd252d'


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


def read_histories(root, dataset, hashes):
    def tracked(name, expected=None):
        path = root / name
        digest = sha256(path)
        require(expected is None or digest == expected, f'{path} hash mismatch')
        hashes[f'{dataset}/{name}'] = digest
        return path
    plan = read_json(tracked('plan.json'))
    plan_sha = hashes[f'{dataset}/plan.json']
    inventory = read_json(tracked('outcomes.json'))
    require(plan['schema'] == 'step04-evaluation-plan-v1' and plan['dataset'] == dataset
            and inventory['plan_sha256'] == plan_sha, f'{dataset} frozen plan/outcomes mismatch')
    outcomes = {(row['fold'], row['seed']): row for row in inventory['outcomes']}
    require(len(outcomes) == len(inventory['outcomes']) == len(SLOTS) and outcomes.keys() == SLOTS,
            f'{dataset} needs nine unique formal outcomes')
    runs, manifests = [], {}
    for (fold, seed), outcome in sorted(outcomes.items()):
        slot = f'slots/fold-{fold}/seed-{seed}'
        require(outcome['status'] == 'completed', f'{dataset}/{slot} incomplete outcome')
        prediction = read_json(tracked(f'{slot}/prediction.json', outcome['prediction_manifest_sha256']))
        run = read_json(tracked(f'{slot}/run.json', prediction['run_sha256']))
        replay = read_json(tracked(f'{slot}/replay.json', prediction['replay_sha256']))
        require((run['identity']['dataset'], run['identity']['fold'], run['identity']['seed']) == (dataset, fold, seed)
                and run['status'] == 'completed' and run['monitor'] == MONITOR
                and run['recipe_sha256'] == plan['folds'][fold]['recipe_sha256_by_seed'][str(seed)]
                and prediction['plan_sha256'] == plan_sha
                and all(prediction[key] == run[key] for key in CHECKPOINT_FIELDS),
                f'{dataset}/{slot} selected checkpoint identity mismatch')
        verify_best_checkpoint(run, f'{dataset}/{slot}')
        require(replay['prediction_max_abs_error'] == 0
                and replay['validation_loss'] == run['best_validation_loss']
                and replay['checkpoint_sha256'] == run['checkpoint_sha256'],
                f'{dataset}/{slot} EMA replay differs')
        require(all(math.isfinite(e['training_loss']) and e['training_loss'] > 0 for e in run['history']),
                f'{dataset}/{slot} invalid training history')
        runs.append(run); manifests[fold, seed] = prediction
    if dataset == 'tenv':
        index = read_json(tracked('checkpoints.json'))
        require(index['plan_sha256'] == plan_sha
                and index['outcomes_sha256'] == hashes[f'{dataset}/outcomes.json'],
                'tenv checkpoint release differs from selected histories')
        verify_checkpoint_index(index['checkpoints'], plan['folds'], manifests)
    return runs


def publish(record, output, evidence):
    hashes = {}
    runs = {dataset: read_histories(record / dataset, dataset, hashes) for dataset in DATASETS}
    require(len({run['identity']['gene_sha256'] for panel in runs.values() for run in panel}) == 1,
            'formal histories require one shared gene axis')
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as directory:
        temporary = Path(directory) / 'learning_curve.png'
        with plt.rc_context({'font.family': 'DejaVu Sans', 'text.color': '#18212c',
                             'axes.labelcolor': '#18212c', 'axes.titlecolor': '#18212c'}):
            plotted = draw_curves(runs, temporary)
        payload = {'schema': 'step04-learning-curves-v1', 'sources': hashes,
                   'renderer_sha256': sha256(Path(__file__)), 'plotted': plotted,
                   'figure_sha256': sha256(temporary)}
        temporary.replace(output)
        write_json(evidence, payload)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record', type=Path, default=STEP / 'execution/formal_record')
    parser.add_argument('--output', type=Path, default=STEP / 'output/learning_curve.png')
    parser.add_argument('--evidence', type=Path, default=STEP / 'execution/formal_record/learning_curves.json')
    args = parser.parse_args()
    publish(args.record, args.output, args.evidence)


if __name__ == '__main__':
    main()
