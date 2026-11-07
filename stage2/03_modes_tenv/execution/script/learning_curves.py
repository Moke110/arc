#!/usr/bin/env python
"""Render persisted native training histories without importing training/evaluation code."""
import argparse
import hashlib
import json
import math
from pathlib import Path

V2_PROTOCOL = 'perturbdiff-v2-all-validation-ema'
V1_PROTOCOL = 'perturbdiff-v1-legacy'
POPULATIONS = {V2_PROTOCOL: 'all 138 groups / 46 drugs',
               V1_PROTOCOL: 'first 40 groups / 14 of 46 drugs'}

V2_MODELS = {'perturbdiff': V2_PROTOCOL, 'cmonge': 'cmonge-v2-full-validation-fixed-draws'}
LEGACY_POPULATIONS = {'perturbdiff': POPULATIONS[V1_PROTOCOL], 'cmonge': 'one seeded held-out condition per event'}

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


def collect(output, *, legacy=False, authenticated_v2=(), only=None):
    records = []
    for model in (V2_MODELS if legacy else MODELS):
        for fold in range(3):
            for seed in SEEDS:
                slot = (model, fold, seed)
                if only is not None and slot != only:
                    continue
                live = output / 'runs' / model / f'fold_{fold}' / f'seed_{seed}'
                archived = output / 'legacy' / f'{model}-v1/runs' / f'fold_{fold}' / f'seed_{seed}'
                migrated = (output / 'legacy' / f'{model}-v1/runs').exists()
                directory = archived if legacy and migrated else live
                run_path = directory / 'run.json'
                run = read(run_path) if run_path.exists() else {}
                record = dict(model=model, fold=fold, seed=seed, status=run.get('status', 'not_dispatched'),
                              train=[], validation=[], history_status='not applicable' if model not in MAPPINGS else 'missing')
                checkpoint = directory / 'checkpoint.json'
                # CMonge publishes durable events before final checkpoint/inference.
                use_checkpoint = checkpoint.exists()
                if slot in authenticated_v2 and run.get('status') != 'completed':
                    use_checkpoint = False
                source = checkpoint if use_checkpoint else directory / 'history.json'
                data = read(source) if model in MAPPINGS and source.exists() else None
                if data is not None:
                    identity = data.get('identity', {})
                    if any(identity.get(k) != v for k, v in [('model', model), ('fold', fold), ('seed', seed)]):
                        raise ValueError(f'checkpoint identity mismatch: {source}')
                if model in V2_MODELS:
                    v2 = V2_MODELS[model]
                    protocol = data.get('protocol', f'{model}-v1-legacy') if data else run.get('protocol')
                    if data is None and slot in authenticated_v2:
                        protocol = v2
                    if protocol and protocol not in (v2, f'{model}-v1-legacy'):
                        raise ValueError(f'unknown {model} protocol: {source}')
                    is_v2 = protocol == v2
                    record.update(protocol=f'{model}-v1-legacy' if legacy else v2,
                                  validation_population=LEGACY_POPULATIONS[model] if legacy else 'all 138 groups / 46 drugs')
                    if legacy and (is_v2 or (data is None and not run_path.exists())):
                        continue
                    if not legacy and not is_v2:
                        if data is not None or run_path.exists():
                            record['status'] = 'legacy_v1_only'
                        records.append(record)
                        continue
                    if is_v2 and data is not None:
                        history = data.get('history', [])
                        if model == 'perturbdiff' and data.get('validation_group_count') != 138:
                            raise ValueError(f'v2 validation population mismatch: {source}')
                        if model == 'cmonge':
                            if any(row.get('validation_conditions_evaluated') != 138 for row in history):
                                raise ValueError(f'v2 validation population mismatch: {source}')
                            if history:
                                record['validation_condition_count'] = 138
                        record['validation_events'] = len(history)
                        for key in ('validation_group_count', 'validation_group_ids', 'validation_group_ids_sha256', 'best_iteration', 'best_step', 'stop_reason'):
                            if key in data:
                                record[key] = data[key]
                if data is not None:
                    x, train, validation, _, _ = MAPPINGS[model]
                    record.update(train=points(data.get('history', []), x, train),
                                  validation=points(data.get('history', []), x, validation),
                                  source=str(source), sha256=hashlib.sha256(source.read_bytes()).hexdigest())
                    record['history_status'] = 'available' if record['train'] else 'empty'
                records.append(record)
    return records


def collect_legacy(output):
    return collect(output, legacy=True)


def bind_attempt(record, hashes, attempt_id, status):
    if record['status'] != status:
        raise ValueError(f'inconsistent attempt status: {attempt_id}')
    if 'source' in record and record['sha256'] != hashes.get(Path(record['source']).name):
        raise ValueError(f'history changed since authentication: {attempt_id}')
    record.update(attempt_id=attempt_id, run_sha256=hashes['run.json'],
                  validation_events=len(record['validation']))

def render(output, publish, final=False, *, supplement=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams['svg.fonttype'] = 'none'
    if supplement is None and (output / 'selection.json').exists():
        raise ValueError('use the authenticated supplement overlay publication command')
    authenticated_v2 = ()
    original_records = []
    if supplement:
        from supplement_overlay import TARGET, slot_key
        selection, proof, _ = supplement
        authenticated_v2 = {slot_key(row) for row in selection['selected'] if row['model'] in V2_MODELS}
        original_output = Path(selection['original_output'])
        original_records = collect(original_output, authenticated_v2={TARGET}, only=TARGET)
        bind_attempt(original_records[0], proof['original_timeout']['files'],
                     'original-cumulative-budget', 'timeout')
    else:
        original_output = output
    records = collect(output, authenticated_v2=authenticated_v2)
    if supplement:
        selected = next(row for row in records if slot_key(row) == TARGET)
        selected_attempt = proof['selected_attempt']
        selected_files = next(row['files'] for row in selection['selected'] if slot_key(row) == TARGET)
        bind_attempt(selected, selected_files, selected_attempt['attempt_id'], selected_attempt['status'])
    legacy_records = collect_legacy(original_output)
    terminal = sum(r['status'] in TERMINAL for r in records)
    if final and terminal != 63:
        raise ValueError(f'final learning curves require 63 terminal outcomes, found {terminal}')
    panels = (*MODELS, 'perturbdiff_legacy', 'cmonge_legacy')
    fig, axes = plt.subplots(6 if supplement else 5, 2, figsize=(16, 27 if supplement else 23), layout='constrained')
    if supplement:
        panels += ('cmonge_original',)
    for ax, model in zip(axes.flat, panels):
        is_legacy = model.endswith('_legacy')
        is_original = model.endswith('_original')
        base_model = model.removesuffix('_legacy').removesuffix('_original')
        source_records = original_records if is_original else legacy_records if is_legacy else records
        subset = [r for r in source_records if r['model'] == base_model]
        model = base_model
        protocol = f'{model}-v1-legacy' if is_legacy else V2_MODELS.get(model, model)
        title = TITLES[model]
        if model in V2_MODELS:
            if is_original:
                title += ' — original v2 timeout, fold 0 seed 17 (excluded)'
            elif is_legacy:
                title += ' — legacy v1 (excluded)'
            else:
                title += ' — formal v2'
                if supplement and model == 'cmonge':
                    title += f" · supplemental slot {selected['status']}"

        ax.set_title(title, fontweight='bold')
        if model not in MAPPINGS:
            ax.text(.5, .55, 'Not applicable', ha='center', transform=ax.transAxes, fontsize=16)
            ax.text(.5, .40, 'Deterministic reference; no fitted training objective', ha='center', transform=ax.transAxes)
            ax.set_axis_off()
            continue
        x, _, _, ylabel, note = MAPPINGS[model]
        if model in V2_MODELS:
            note = ('Validation: ' + (LEGACY_POPULATIONS[model] if is_legacy else 'all 138 groups / 46 drugs') + '. ' +
                    ('Historical operational context; excluded from formal slots.' if is_legacy else
                     'Seed restart; deterministic group-weighted validation.'))
        counts = []
        for key, color in [('train', '#1565c0'), ('validation', '#228b22')]:
            counts.append(sum(bool(r[key]) for r in subset))
            for record in subset:
                if record[key]:
                    xs, ys = zip(*record[key])
                    line, = ax.plot(xs, ys, color=color, alpha=.45 + .25 * record['fold'], linewidth=1,
                                    linestyle=('-', '--', ':')[SEEDS.index(record['seed'])], marker='.' if len(xs) == 1 else None)
                    attempt = record.get('attempt_id', protocol)
                    line.set_gid(f'{attempt}-fold-{record["fold"]}-seed-{record["seed"]}-{key}')
        ax.set_xlabel(x.capitalize() + ' (actual native coordinate)')
        ax.set_ylabel(ylabel, fontsize=9)
        ax.grid(alpha=.18)
        count_label = (f'Legacy histories: train {counts[0]} · validation {counts[1]}' if is_legacy else
                       f'Train {counts[0]}/9 · validation {counts[1]}/9 histories; missing {9-counts[0]}/{9-counts[1]}')
        if is_original:
            count_label = f'Original attempt histories: train {counts[0]}/1 · validation {counts[1]}/1'
            note = 'Separate original v2 attempt; validation averages all 138 conditions.'
        if supplement and model == 'cmonge' and not (is_original or is_legacy):
            note = (f"Fold 0 seed 17 selected supplement ({selected['status']}): "
                    'fresh 14,400 s training + 1,800 s inference from seed initialization.\n'
                    'Original cumulative-budget timeout shown separately.')
        status_counts = {status: sum(r['status'] == status for r in subset) for status in sorted({r['status'] for r in subset})}
        status_label = ' · '.join(f'{status}: {count}' for status, count in status_counts.items())
        ax.text(.02, .98, count_label + '\n' + status_label,
                transform=ax.transAxes, va='top', fontsize=10,
                bbox=dict(facecolor='white', alpha=.85, edgecolor='none'))
        ax.text(0, -.23, note, transform=ax.transAxes, fontsize=9, wrap=True)
    axes.flat[-1].set_axis_off()
    axes.flat[-1].text(0, .85,
        '3 folds × 3 seeds = 9 expected formal histories per fitted model\n\n'
        'Both models’ v1 is separate and never counts toward nine v2 slots.\n'
        'Different validation populations; no v1/v2 history splicing.\n\n'
        'Missing / unfinished histories are omitted, never reconstructed.\n'
        'CMonge v2 durable events may precede final inference.\n'
        'Line styles: solid / dashed / dotted = seeds 17 / 42 / 2026.\n'
        'Opacity: light / medium / dark = folds 0 / 1 / 2.\n'
        'Exact fold × seed identities also appear in SVG IDs and JSON.\n\n'
        'Native objectives and reductions differ across panels.\n'
        'No smoothing, normalization, extrapolation or test-truth access.',
        transform=axes.flat[-1].transAxes, va='top', fontsize=12, linespacing=1.6)
    title = f'Native learning curves — {"final" if final else "interim snapshot"} ({terminal}/63 terminal outcomes)'
    if supplement:
        from deliver import disclosure_lines
        title += '\n' + '\n'.join(disclosure_lines(proof))
    fig.suptitle(title, fontsize=15 if supplement else 19)
    fig.legend(handles=[Line2D([0], [0], color='#1565c0', label='Training'),
                        Line2D([0], [0], color='#228b22', label='Validation')], loc='outside lower center', ncol=2)
    publish.mkdir(parents=True, exist_ok=True)
    for ext in ('png', 'svg'):
        fig.savefig(publish / f'learning_curves.{ext}', dpi=160)
    plt.close(fig)
    provenance = dict(schema='native-learning-curves-v2', final=final, terminal_outcomes=terminal,
                      renderer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), records=records,
                      legacy_records=legacy_records,
                      perturbdiff_v2=dict(expected_slots=9,
                          available_histories=sum(bool(r['train']) for r in records if r['model']=='perturbdiff'),
                          missing_histories=sum(not r['train'] for r in records if r['model']=='perturbdiff')))
    provenance['cmonge_v2'] = dict(expected_slots=9, available_histories=sum(bool(r['train']) for r in records if r['model']=='cmonge'), missing_histories=sum(not r['train'] for r in records if r['model']=='cmonge'))
    if supplement:
        provenance.update(original_v2_records=original_records, selection_sha256=proof['selection_sha256'],
                          evaluation_provenance_sha256=supplement[2], budget_exception=proof)
    (publish / 'learning_curves.json').write_text(json.dumps(provenance, indent=2) + '\n')
    return provenance


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('stage2/03_modes_tenv/output/formal'))
    parser.add_argument('--publish', type=Path, default=Path('stage2/03_modes_tenv/output'))
    parser.add_argument('--final', action='store_true', help='require all 63 terminal outcomes')
    args = parser.parse_args()
    render(args.output, args.publish, args.final)
