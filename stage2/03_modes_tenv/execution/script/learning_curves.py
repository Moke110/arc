#!/usr/bin/env python
"""Render persisted native training histories without importing training/evaluation code."""
import argparse
import hashlib
import json
import math
from pathlib import Path

SEEDS = (17, 42, 2026)
TERMINAL = {'completed', 'failed', 'timeout', 'excluded'}
# x, training key, validation key, native objective label, sampling limitation
MAPPINGS = {
    'prnet': ('epoch', 'training_loss', 'validation_loss',
              'Gaussian NLL (native mean; dimensionless)', 'Epoch means; fixed validation noise stream.'),
    'cpa': ('epoch', 'training_loss', 'validation_loss',
            'Reconstruction NLL (gene sum / cell; dimensionless)',
            'Reconstruction component only; native adversarial training also occurs.'),
    'chemcpa': ('epoch', 'autoencoder_objective', 'validation_autoencoder_objective',
                'Gaussian NLL − weighted drug BCE (dimensionless)',
                'Combined native autoencoder objective; alternating adversarial updates.'),
    'cmonge': ('iteration', 'training_loss', 'validation_loss',
               'Sinkhorn divergence + Monge gap (native log-expression units)',
               'Event samples: one train condition / one seeded held-out condition.'),
    'perturbdiff': ('step', 'training_loss', 'validation_loss',
                   'Weighted diffusion MSE + energy MMD (dimensionless)',
                   'Last training batch vs bounded validation batches; EMA validation.'),
}
MODELS = (*MAPPINGS, 'no_change', 'delta_mean')
TITLES = dict(zip(MODELS, ('PRnet', 'CPA', 'chemCPA', 'CMonge', 'PerturbDiff', 'no_change', 'delta_mean')))


def read(path):
    return json.loads(path.read_text())


def points(history, xkey, ykey):
    """Preserve actual native coordinates; reject malformed/nonfinite history."""
    result = [(float(row[xkey]), float(row[ykey])) for row in history]
    if any(not math.isfinite(v) for pair in result for v in pair):
        raise ValueError('nonfinite native history')
    if any(b[0] <= a[0] for a, b in zip(result, result[1:])):
        raise ValueError('native history coordinates must increase')
    return result


def collect(output):
    records = []
    for model in MODELS:
        for fold in range(3):
            for seed in SEEDS:
                directory = output / 'runs' / model / f'fold_{fold}' / f'seed_{seed}'
                run_path = directory / 'run.json'
                run = read(run_path) if run_path.exists() else {}
                record = dict(model=model, fold=fold, seed=seed, status=run.get('status', 'not_dispatched'),
                              train=[], validation=[], history_status='not applicable' if model not in MAPPINGS else 'missing')
                checkpoint = directory / 'checkpoint.json'
                if model in MAPPINGS and checkpoint.exists():
                    raw = checkpoint.read_bytes()
                    data = json.loads(raw)
                    identity = data.get('identity', {})
                    if any(identity.get(k) != v for k, v in [('model', model), ('fold', fold), ('seed', seed)]):
                        raise ValueError(f'checkpoint identity mismatch: {checkpoint}')
                    x, train, validation, _, _ = MAPPINGS[model]
                    record.update(train=points(data.get('history', []), x, train),
                                  validation=points(data.get('history', []), x, validation),
                                  source=str(checkpoint), sha256=hashlib.sha256(raw).hexdigest())
                    record['history_status'] = 'available' if record['train'] else 'empty'
                records.append(record)
    return records


def render(output, publish, final=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams['svg.fonttype'] = 'none'
    records = collect(output)
    terminal = sum(r['status'] in TERMINAL for r in records)
    if final and terminal != 63:
        raise ValueError(f'final learning curves require 63 terminal outcomes, found {terminal}')
    fig, axes = plt.subplots(4, 2, figsize=(16, 18), layout='constrained')
    for ax, model in zip(axes.flat, MODELS):
        subset = [r for r in records if r['model'] == model]
        ax.set_title(TITLES[model], fontweight='bold')
        if model not in MAPPINGS:
            ax.text(.5, .55, 'Not applicable', ha='center', transform=ax.transAxes, fontsize=16)
            ax.text(.5, .40, 'Deterministic reference; no fitted training objective', ha='center', transform=ax.transAxes)
            ax.set_axis_off()
            continue
        x, _, _, ylabel, note = MAPPINGS[model]
        counts = []
        for key, color in [('train', '#1565c0'), ('validation', '#228b22')]:
            counts.append(sum(bool(r[key]) for r in subset))
            for index, record in enumerate(subset):
                if record[key]:
                    xs, ys = zip(*record[key])
                    line, = ax.plot(xs, ys, color=color, alpha=.45 + .25 * record['fold'], linewidth=1,
                                    linestyle=('-', '--', ':')[index % 3], marker='.' if len(xs) == 1 else None)
                    line.set_gid(f'{model}-fold-{record["fold"]}-seed-{record["seed"]}-{key}')
        ax.set_xlabel(x.capitalize() + ' (actual native coordinate)')
        ax.set_ylabel(ylabel, fontsize=9)
        ax.grid(alpha=.18)
        ax.text(.02, .98, f'Train {counts[0]}/9 · validation {counts[1]}/9 histories',
                transform=ax.transAxes, va='top', fontsize=10,
                bbox=dict(facecolor='white', alpha=.85, edgecolor='none'))
        ax.text(0, -.23, note, transform=ax.transAxes, fontsize=9, wrap=True)
    axes.flat[-1].set_axis_off()
    axes.flat[-1].text(0, .85,
        '3 folds × 3 seeds = 9 expected histories per fitted model\n\n'
        'Missing / unfinished histories are omitted, never reconstructed.\n'
        'Active adapters may persist history only after training + inference.\n'
        'Line styles: solid / dashed / dotted = seeds 17 / 42 / 2026.\n'
        'Opacity: light / medium / dark = folds 0 / 1 / 2.\n'
        'Exact fold × seed identities also appear in SVG IDs and JSON.\n\n'
        'Native objectives and reductions differ across panels.\n'
        'No smoothing, normalization, extrapolation or test-truth access.',
        transform=axes.flat[-1].transAxes, va='top', fontsize=12, linespacing=1.6)
    fig.suptitle(f'Native learning curves — {"final" if final else "interim snapshot"} ({terminal}/63 terminal outcomes)', fontsize=19)
    fig.legend(handles=[Line2D([0], [0], color='#1565c0', label='Training'),
                        Line2D([0], [0], color='#228b22', label='Validation')], loc='outside lower center', ncol=2)
    publish.mkdir(parents=True, exist_ok=True)
    for ext in ('png', 'svg'):
        fig.savefig(publish / f'learning_curves.{ext}', dpi=160)
    plt.close(fig)
    provenance = dict(schema='native-learning-curves-v1', final=final, terminal_outcomes=terminal,
                      renderer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), records=records)
    (publish / 'learning_curves.json').write_text(json.dumps(provenance, indent=2) + '\n')
    return provenance


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('stage2/03_modes_tenv/output/formal'))
    parser.add_argument('--publish', type=Path, default=Path('stage2/03_modes_tenv/output'))
    parser.add_argument('--final', action='store_true', help='require all 63 terminal outcomes')
    args = parser.parse_args()
    render(args.output, args.publish, args.final)
