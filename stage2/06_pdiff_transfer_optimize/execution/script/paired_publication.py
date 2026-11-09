#!/usr/bin/env python3
"""Publish verified paired Transfer scores or replay the figure from their numeric record."""
import argparse
import math
from pathlib import Path
from statistics import mean
import tempfile

import paired_reference as reference
from paired_evidence import verify_matrix
from artifacts import SEEDS, read_json, sha256, write_json
from formal import slot_lock
from formal_matrix import bind
from training_inputs import json_hash

SCHEMA = 'step06-paired-publication-v1'
STEP = Path(__file__).resolve().parents[2]
FIGURE = STEP / 'output/03-paired-transfer-optimizations.png'
RECORD = STEP / 'execution/paired_publication.json'
PANELS = (
    ('Dose representation', (reference.ARM, 'random', 'log')),
    ('Learning rate and trainable parameters', (reference.ARM, 'all-half', 'main-half', 'adaLN-only')),
)
ARMS = tuple(dict.fromkeys(arm for _, arms in PANELS for arm in arms))
METRICS = (('mean_ap100', 'AUPRC (mean AP@100)'), ('pcc', 'PCC'))
COLORS = dict(zip(ARMS, ('#303843', '#006fac', '#00836c', '#ac4c0a', '#8457ad', '#a83353'), strict=True))


def plot_points(runs):
    if set(runs) != set(ARMS):
        raise ValueError('publication requires all five arms and original Transfer')
    plotted = {}
    for arm in ARMS:
        records = runs[arm]
        if tuple(run['seed'] for run in records) != SEEDS:
            raise ValueError(f'{arm} requires the three paired seeds in order')
        for run in records:
            if tuple(float(row['dose_um']) for row in run['per_dose']) != reference.DOSES:
                raise ValueError(f'{arm} seed {run["seed"]} requires all three doses in order')
        points = []
        for dose, *rows in zip(reference.DOSES, *(run['per_dose'] for run in records), strict=True):
            point = {'dose_um': dose}
            for metric, _ in METRICS:
                values = [row[metric] for row in rows]
                if any(value is not None and not math.isfinite(value) for value in values):
                    raise ValueError(f'{arm} dose {dose} has a non-finite {metric}')
                undefined = [seed for seed, value in zip(SEEDS, values, strict=True) if value is None]
                point[metric] = None if undefined else mean(values)
                point[f'{metric}_undefined_seeds'] = undefined
            points.append(point)
        plotted[arm] = points
    return plotted


def build_figure(points):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.ticker import NullLocator

    fig = Figure(figsize=(10, 7.4), dpi=180, facecolor='white')
    FigureCanvasAgg(fig)
    axes = fig.subplots(2, 2)
    fig.subplots_adjust(left=.09, right=.98, top=.78, bottom=.09, hspace=.74, wspace=.23)
    fig.suptitle('Paired PerturbDiff Transfer optimizations', y=.975, fontsize=18, weight='bold')
    fig.text(.5, .927, 'Cold-cm-drug fold 0 validation  |  Equal mean of seeds 17, 42 and 2026',
             ha='center', fontsize=10, color='#424b56')
    for (heading, arms), panels, header_y in zip(PANELS, axes, (.882, .425), strict=True):
        fig.text(.5, header_y, heading, ha='center', fontsize=12, weight='bold')
        for (metric, label), ax in zip(METRICS, panels, strict=True):
            for arm in arms:
                ax.plot(reference.DOSES, [point[metric] for point in points[arm]],
                        label='Original Transfer' if arm == reference.ARM else arm,
                        color=COLORS[arm], marker='o', markersize=4.5, linewidth=1.8)
            ax.set_xscale('log')
            ax.set_xticks(reference.DOSES, ['0.2', '1', '10'])
            ax.xaxis.set_minor_locator(NullLocator())
            ax.set_xlabel('Tenv dose (µM)')
            ax.set_ylabel(label)
            ax.grid(axis='y', alpha=.22)
            ax.spines[['top', 'right']].set_visible(False)
        fig.legend(handles=panels[0].lines, loc='upper center', bbox_to_anchor=(.5, header_y - .012),
                   ncol=len(arms), frameon=False, fontsize=9)
    return fig


def relative(binding, runtime):
    return {**binding, 'path': str(Path(binding['path']).relative_to(runtime))}


def fit_record(fit, runtime):
    return {
        'selected_step': fit['metadata']['best_step'],
        'validation_diffusion_loss': fit['metadata']['best_validation_loss'],
        'optimizer_steps': read_json(fit['history']['path'])['optimizer_steps'],
        **{key: relative(fit[key], runtime) for key in ('checkpoint', 'history', 'fit_seal')},
    }


def run_record(run, source, seed, runtime):
    return {
        'seed': seed,
        'source': fit_record(source, runtime),
        'target': fit_record(run['target'], runtime),
        **{key: relative(run[key], runtime) for key in ('plan', 'scores', 'prediction')},
        'score_files': {name: relative(binding, runtime) for name, binding in run['score_files'].items()},
        'per_dose': run['per_dose'],
    }


def resources(root, runtime):
    expected = {f'{arm}-{seed}-{phase}' for arm in ARMS[1:] for seed in SEEDS
                for phase in ('verify', 'fit', 'predict', 'score')}
    expected.update(f'reference-{seed}-{phase}' for seed in SEEDS for phase in ('verify', 'predict', 'score'))
    expected.update(f'log-source-{seed}-fit' for seed in SEEDS)
    expected.update(f'{arm}-matrix-verify' for arm in ARMS[1:])
    files = {}
    ledger = {}
    measurements = ('memory_max_gib', 'memory_swap_max_gib', 'peak_host_rss_gib',
                    'peak_workstation_gpu_used_gib', 'peak_gpu_utilization_percent')
    for status in sorted((root / 'schedule').glob('*/scheduler/status.json')):
        jobs = read_json(status)['jobs']
        tracked = expected & jobs.keys()
        if not tracked:
            continue
        slot = status.parent.parent
        files[slot.name] = [relative(bind(path), runtime) for path in
                            (slot / 'queue.json', status, status.with_name('telemetry.jsonl'))]
        for job in sorted(tracked):
            entry = jobs[job]
            if job in ledger or entry['status'] != 'completed':
                raise ValueError(f'{job} needs one completed scheduler record')
            ledger[job] = {'slot': slot.name, 'attempts': [
                {'outcome': item['outcome'], 'wall_seconds': item['finished_at'] - item['started_at'],
                 **{key: item.get(key) for key in measurements}}
                for item in entry['attempts']
            ]}
    if missing := expected - ledger.keys():
        raise ValueError(f'missing scheduler evidence: {", ".join(sorted(missing))}')
    return {'jobs': ledger, 'files': files}


def collect(root, runtime):
    matrices = {}
    originals = {}
    runs = {arm: [] for arm in ARMS}
    for arm in ARMS[1:]:
        matrix = verify_matrix(root, arm)
        path = root / arm / 'matrix.json'
        binding = bind(path)
        if matrix != read_json(path):
            raise ValueError(f'{arm} matrix differs from the verified cells')
        matrices[arm] = relative(binding, runtime)
        for cell in matrix['cells']:
            seed = cell['seed']
            original = (matrix['protocol_sha256'], cell['reference'])
            if originals.setdefault(seed, original) != original:
                raise ValueError(f'{arm} differs from the shared seed {seed} Transfer reference')
            runs[arm].append({**run_record(cell, cell['source'], seed, runtime),
                              'paired': relative(cell['paired'], runtime)})
    protocols = {protocol for protocol, _ in originals.values()}
    if len(protocols) != 1:
        raise ValueError('all arms and seeds must share one frozen validation protocol')
    for seed in SEEDS:
        _, original = originals[seed]
        plan = read_json(original['plan']['path'])
        verified = reference.verify_baseline(runtime, seed)
        if json_hash(plan) != json_hash(verified):
            changed = [key for key in sorted(plan.keys() | verified.keys())
                       if json_hash(plan.get(key)) != json_hash(verified.get(key))]
            raise ValueError(f'seed {seed} reference differs from published Step 05 in {", ".join(changed)}')
        runs[reference.ARM].append(run_record(original, plan['source'], seed, runtime))
    return {'schema': SCHEMA, 'regime': reference.REGIME, 'fold': 0,
            'protocol_sha256': protocols.pop(), 'matrices': matrices, 'runs': runs,
            'plotted': plot_points(runs), 'resources': resources(root, runtime)}


def publish(root, runtime, record, figure):
    data = collect(root, runtime)
    figure.parent.mkdir(parents=True, exist_ok=True)
    record.parent.mkdir(parents=True, exist_ok=True)
    with slot_lock(figure.parent / '.paired-publication-lock'), tempfile.TemporaryDirectory(dir=figure.parent) as tmp:
        png = Path(tmp) / 'figure.png'
        build_figure(data['plotted']).savefig(png)
        data['figure_sha256'] = sha256(png)
        if figure.exists():
            bind(figure, data['figure_sha256'])
        if record.exists() and read_json(record) != data:
            raise ValueError('published record differs from verified inputs')
        if not figure.exists():
            png.replace(figure)
        if not record.exists():
            write_json(record, data)


def replay(record, figure):
    data = read_json(record)
    if data['schema'] != SCHEMA or data['plotted'] != plot_points(data['runs']):
        raise ValueError('numeric record does not match its plotted points')
    figure.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=figure.parent) as tmp:
        png = Path(tmp) / 'figure.png'
        build_figure(data['plotted']).savefig(png)
        bind(png, data['figure_sha256'])
        if figure.exists():
            bind(figure, data['figure_sha256'])
        else:
            png.replace(figure)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    phases = parser.add_subparsers(dest='phase', required=True)
    publisher = phases.add_parser('publish')
    publisher.add_argument('root', type=Path)
    publisher.add_argument('--runtime', type=Path, default=reference.base.REPO)
    publisher.add_argument('--record', type=Path, default=RECORD)
    replayer = phases.add_parser('replay')
    replayer.add_argument('record', type=Path)
    for phase in (publisher, replayer):
        phase.add_argument('--figure', type=Path, default=FIGURE)
    args = parser.parse_args(argv)
    if args.phase == 'publish':
        publish(args.root.resolve(), args.runtime.resolve(), args.record.resolve(), args.figure.resolve())
    else:
        replay(args.record.resolve(), args.figure.resolve())
    print(f'{args.phase}: {args.figure} SHA-256 {sha256(args.figure)}', flush=True)


if __name__ == '__main__':
    main()
