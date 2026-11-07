#!/usr/bin/env python
"""Render persisted native training histories without importing training/evaluation code."""
import argparse
import hashlib
import json
import math
import os
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


def points(history, xkey, ykey):
    """Preserve actual native coordinates; reject malformed/nonfinite history."""
    result = [(float(row[xkey]), float(row[ykey])) for row in history]
    if any(not math.isfinite(v) for pair in result for v in pair):
        raise ValueError('nonfinite native history')
    if any(b[0] <= a[0] for a, b in zip(result, result[1:])):
        raise ValueError('native history coordinates must increase')
    return result


def curve_file(path, hashes):
    if hashes is None:
        return (path.read_bytes(), None) if path.exists() else (None, None)
    from supplement_overlay import checked_path

    if path.name not in hashes:
        if os.path.lexists(path):
            raise ValueError(f'unreviewed curve source: {path}')
        return None, None
    raw = checked_path(path).read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != hashes[path.name]:
        raise ValueError(f'history or status changed since authentication: {path}')
    return raw, digest


def curve_data(raw, slot, source):
    data = json.loads(raw)
    if not isinstance(data, dict):
        raise ValueError(f'expected native history object: {source}')
    identity = data.get('identity')
    if not isinstance(identity, dict) or any(identity.get(k) != v for k, v in zip(('model', 'fold', 'seed'), slot)):
        raise ValueError(f'checkpoint identity mismatch: {source}')
    return data


def curve_points(data, model, protocol, source, *, authenticated=False):
    history = data['history'] if authenticated else data.get('history', [])
    if not isinstance(history, list) or any(not isinstance(row, dict) for row in history):
        raise ValueError(f'invalid native history events: {source}')
    result = {}
    if model in V2_MODELS and protocol == V2_MODELS[model]:
        if model == 'perturbdiff' and data.get('validation_group_count') != 138:
            raise ValueError(f'v2 validation population mismatch: {source}')
        if model == 'cmonge':
            if any(row.get('validation_conditions_evaluated') != 138 for row in history):
                raise ValueError(f'v2 validation population mismatch: {source}')
            if history:
                result['validation_condition_count'] = 138
        result['validation_events'] = len(history)
        for key in ('validation_group_count', 'validation_group_ids', 'validation_group_ids_sha256',
                    'best_iteration', 'best_step', 'stop_reason'):
            if key in data:
                result[key] = data[key]
    x, train, validation, _, _ = MAPPINGS[model]
    if authenticated and any(type(row[key]) not in (int, float)
                             for row in history for key in (x, train, validation)):
        raise ValueError(f'expected numeric native history events: {source}')
    result.update(train=points(history, x, train), validation=points(history, x, validation))
    result['history_status'] = 'available' if result['train'] else 'empty'
    return result


def collect(output, *, legacy=False, authenticated_v2=(), only=None, inventories=None, tolerate_invalid=()):
    records = []
    for model in (V2_MODELS if legacy else MODELS):
        for fold in range(3):
            for seed in SEEDS:
                slot = (model, fold, seed)
                if only is not None and slot != only:
                    continue
                live = output / 'runs' / model / f'fold_{fold}' / f'seed_{seed}'
                archive = output / 'legacy' / f'{model}-v1/runs'
                archived = legacy and (inventories is not None or archive.exists())
                directory = archive / f'fold_{fold}' / f'seed_{seed}' if archived else live
                hashes = inventories[slot] if inventories is not None else None
                run_raw, run_sha = curve_file(directory / 'run.json', hashes)
                run = json.loads(run_raw) if run_raw is not None else {}
                record = dict(model=model, fold=fold, seed=seed, status=run.get('status', 'not_dispatched'),
                              train=[], validation=[], history_status='not applicable' if model not in MAPPINGS else 'missing')
                if legacy and run_sha is not None:
                    record['run_sha256'] = run_sha
                if model not in MAPPINGS:
                    records.append(record)
                    continue
                checkpoint = directory / 'checkpoint.json'
                # CMonge publishes durable events before final checkpoint/inference.
                use_checkpoint = (checkpoint.exists() if hashes is None else 'checkpoint.json' in hashes)
                if slot in authenticated_v2 and record['status'] != 'completed':
                    use_checkpoint = False
                source = checkpoint if use_checkpoint else directory / 'history.json'
                raw, digest = curve_file(source, hashes)
                if hashes is not None:
                    # Also reject unreviewed competing inputs and changed archived JSON.
                    alternate = directory / ('history.json' if use_checkpoint else 'checkpoint.json')
                    curve_file(alternate, hashes)
                if model in V2_MODELS:
                    v2 = V2_MODELS[model]
                    v1 = f'{model}-v1-legacy'
                    record.update(protocol=v1 if legacy else v2,
                                  validation_population=LEGACY_POPULATIONS[model] if legacy else 'all 138 groups / 46 drugs')
                try:
                    data = curve_data(raw, slot, source) if raw is not None else None
                    protocol = None
                    if model in V2_MODELS:
                        protocol = data.get('protocol', v1) if data else run.get('protocol')
                        if data is None and slot in authenticated_v2:
                            protocol = v2
                        if protocol and protocol not in (v2, v1):
                            raise ValueError(f'unknown {model} protocol: {source}')
                        if slot in authenticated_v2 and protocol != v2:
                            raise ValueError(f'expected authenticated {model} v2 history: {source}')
                        if legacy and hashes is not None and protocol == v2:
                            raise ValueError(f'unexpected v2 history in the legacy archive: {source}')
                        if legacy and (protocol == v2 or (data is None and run_raw is None)):
                            continue
                        if not legacy and protocol != v2:
                            if data is not None or run_raw is not None:
                                record['status'] = 'legacy_v1_only'
                            records.append(record)
                            continue
                    if data is not None:
                        record.update(curve_points(data, model, protocol, source, authenticated=digest is not None),
                                      source=str(source),
                                      sha256=digest or hashlib.sha256(raw).hexdigest())
                except (ValueError, KeyError, TypeError, OverflowError):
                    if (slot not in tolerate_invalid or digest is None
                            or record['status'] not in ('failed', 'timeout')):
                        raise
                    record.update(history_status='invalid', source=str(source), sha256=digest, validation_events=0)
                records.append(record)
    return records


def collect_legacy(output, *, journal_sha256=None):
    if journal_sha256 is None:
        return collect(output, legacy=True)
    from supplement_overlay import checked_path, pinned_json

    journal = pinned_json(output / 'joint-v2-apply.json', journal_sha256)
    if journal['schema'] != 'joint-v2-apply-v1' or journal['status'] != 'committed':
        raise ValueError('expected the committed legacy migration journal')
    inventories = {}
    for model in V2_MODELS:
        archive = checked_path(output / 'legacy' / f'{model}-v1/runs', directory=True)
        matches = [move for move in journal['moves'] if move['destination'] == str(archive)]
        if len(matches) != 1 or matches[0]['source'] != str(output / 'runs' / model):
            raise ValueError(f'legacy archive is missing its reviewed move: {archive}')
        inventory = matches[0]['inventory']
        for fold in range(3):
            for seed in SEEDS:
                relative = Path(f'fold_{fold}') / f'seed_{seed}'
                for directory in (relative.parent, relative):
                    if str(directory) in inventory:
                        checked_path(archive / directory, directory=True)
                inventories[model, fold, seed] = {
                    name: inventory[str(relative / name)] for name in ('run.json', 'history.json', 'checkpoint.json')
                    if str(relative / name) in inventory
                }
    return collect(output, legacy=True, inventories=inventories)


def bind_attempt(record, hashes, attempt_id, status):
    if record['status'] != status:
        raise ValueError(f'inconsistent attempt status: {attempt_id}')
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
    original_records = []
    if supplement:
        from supplement_overlay import TARGET, slot_key
        selection, proof, _ = supplement
        original_output = Path(selection['original_output'])
        original_files = proof['original_timeout']['files']
        original_records = collect(original_output, authenticated_v2={TARGET}, only=TARGET,
                                   inventories={TARGET: original_files}, tolerate_invalid={TARGET})
        bind_attempt(original_records[0], original_files, 'original-cumulative-budget', 'timeout')
        inventories = {slot_key(row): row['files'] for row in selection['selected']}
        authenticated_v2 = {slot for slot in inventories if slot[0] in V2_MODELS}
        selected_attempt = proof['selected_attempt']
        if selected_attempt['account'].get('stop_reason') == 'invalid_native_artifacts':
            raise ValueError('invalid native artifacts cannot settle the supplemental slot')
        records = collect(output, authenticated_v2=authenticated_v2, inventories=inventories,
                          tolerate_invalid={TARGET})
        selected = next(row for row in records if slot_key(row) == TARGET)
        bind_attempt(selected, inventories[TARGET], selected_attempt['attempt_id'], selected_attempt['status'])
        legacy_records = collect_legacy(original_output, journal_sha256=selection['migration_journal_sha256'])
    else:
        records = collect(output)
        legacy_records = collect_legacy(output)
    terminal = sum(r['status'] in TERMINAL for r in records)
    if final and terminal != 63:
        raise ValueError(f'final learning curves require 63 terminal outcomes, found {terminal}')
    panels = (*MODELS, 'perturbdiff_legacy', 'cmonge_legacy')
    fig, axes = plt.subplots(6 if supplement else 5, 2, figsize=(16, 27 if supplement else 23), layout='constrained')
    if supplement:
        panels += ('cmonge_original',)
        axes.flat[len(panels)].set_axis_off()
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
        invalid = sum(r['history_status'] == 'invalid' for r in subset)
        if invalid:
            count_label += f' · invalid histories {invalid}'
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
