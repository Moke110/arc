#!/usr/bin/env python3
"""Publish and replay the Step 06 diagnosis from frozen score and run evidence."""
import argparse
from pathlib import Path
import tempfile

import baseline
import confirmation
import dose_reset
import learning_rate as lr
import replication
from artifacts import SEEDS, read_json, rows, sha256, write_csv, write_json
from formal import slot_lock
from formal_matrix import bind, check_binding
from screen import THRESHOLD, screen_decision

METRICS = lr.METRICS
ROOTS = {'baseline': 'baseline', 'dose_reset': 'dose_reset', 'learning_rate': 'learning-rate'}
FIGURES = ('01-paired-validation-effects.png', '02-tahoe-exposure.png')


def diagnosis(screen, arm):
    """Apply the attribution gate to the declared primary mode only."""
    gate = screen_decision({int(seed): gain for seed, gain in screen['gains'].items()})
    mode = arm['primary_mode'] if arm else None
    pooled = arm['comparisons'][mode]['pooled'] if mode else None
    status = gate['status']
    if status == 'promoted':
        metric = pooled['mean_ap100'] if pooled else None
        if metric is None or metric['paired_seeds'] != list(SEEDS):
            status = 'incomplete_confirmation'
        elif metric['gain'] <= 0:
            status = 'nonpositive_confirmation_effect'
        elif metric['cm_gap'] is None or not metric['matched_cm_seeds']:
            status = 'matching_cm_unavailable'
        elif metric['cm_gap'] <= 0:
            status = 'nonpositive_cm_gap'
        else:
            status = 'major_contribution' if metric['recovery'] >= 0.5 else 'below_half_gap'
    return {'screen': screen, 'confirmation_mode': mode, 'confirmation': pooled,
            'attribution_status': status, 'major_contribution': status == 'major_contribution'}


def validation_scores(output, screens):
    baseline_seeds = {seed for seed in SEEDS
                      if (output / 'baseline' / f'fold-0-seed-{seed}/scores/scores.json').is_file()}
    if 17 not in baseline_seeds:
        raise ValueError('publication requires the shared seed 17 validation baseline')
    for arm, root in ROOTS.items():
        score_seeds = (baseline_seeds if arm == 'baseline'
                       else {int(seed) for seed in screens[arm]['scores']})
        for seed in sorted(score_seeds, key=SEEDS.index):
            slot = output / root / f'fold-0-seed-{seed}'
            score = bind(slot / 'scores/scores.json')
            if arm != 'baseline' and score != screens[arm]['scores'][str(seed)]:
                raise ValueError(f'changed frozen score for {arm}, seed {seed}')
            report = read_json(score['path'])
            plan = bind(slot / 'plan.json', report['plan_sha256'])
            if (report['fold'], report['seed']) != (0, seed):
                raise ValueError(f'unexpected validation score in {slot}')
            for name, digest in report['files'].items():
                bind(slot / 'scores' / name, digest)
            yield arm, slot, report, plan


def replay_validation(runtime, output, arm, slot, score, binding):
    digest = binding['sha256']
    plan = read_json(binding['path'])
    if arm == 'baseline':
        baseline.score_baseline(runtime, slot, plan, digest)
    elif arm == 'dose_reset':
        base = read_json(check_binding(plan['baseline'])['path'])
        dose_reset.score_reset(runtime, output, plan, base, digest)
    else:
        lr.score(runtime, slot, lr.load_plan(slot, score['seed'], digest), digest)


def summary_by_arm(rows_in):
    arms = {}
    for row in rows_in:
        arms.setdefault(row['arm'], []).append(row)
    return arms


def pair_rows(fold, seed, arm, mode, result, original=None, treatment=None, cm=None,
              steps=(None, None, None), cm_reason=None):
    """Keep each paired direction together with its values and support."""
    if original is None:
        yield {'fold': fold, 'seed': seed, 'arm': arm, 'mode': mode, 'dose_um': 'pooled',
               'status': result['status'], 'reason': result.get('reason'), 'cm_reason': cm_reason,
               'metric': 'mean_ap100'}
        return
    for original_row, treatment_row, cm_row in zip(
            original, treatment or [None] * len(original), cm or [None] * len(original), strict=True):
        dose, count = original_row['dose_um'], original_row['group_count']
        if any((row['dose_um'], row['group_count']) != (dose, count)
               for row in (treatment_row, cm_row) if row is not None):
            raise ValueError('paired dose or group count differs')
        paired = lr.contrast(cm_row, original_row, treatment_row or dict.fromkeys(METRICS))
        for metric in METRICS:
            support = f'{metric}_valid_groups'
            row = {'fold': fold, 'seed': seed, 'arm': arm, 'mode': mode, 'dose_um': dose,
                   'group_count': count, 'status': result['status'], 'reason': result.get('reason'),
                   'cm_reason': cm_reason, 'metric': metric, 'transfer_step': steps[0],
                   'intervention_step': steps[1], 'cm_step': steps[2]}
            for key, values in zip(('transfer', 'intervention', 'cm'),
                                   (original_row, treatment_row, cm_row), strict=True):
                row[key] = values.get(metric) if values else None
                row[f'{key}_valid_groups'] = values.get(support) if values else None
            yield {**row, **{key: paired[metric][key]
                             for key in ('gain', 'cm_gap', 'residual_gap', 'recovery')}}


def validation_pairs(arm, score):
    seed = score['seed']
    if arm != 'learning_rate':
        summaries = summary_by_arm(score['per_dose'])
        completed = score.get('overall_gain') is not None or arm == 'baseline'
        result = {'status': 'completed' if completed else 'incomplete',
                  'reason': None if completed else 'pooled directional mean AP@100 undefined'}
        steps = tuple(score['checkpoints'].get(key, {}).get('metadata', {}).get('best_step')
                      for key in ('Transfer', arm, 'CM'))
        yield from pair_rows(0, seed, arm, 'loss-selected', result, summaries['Transfer'],
                             summaries.get(arm), summaries['CM'], steps)
        return
    summaries = score['per_dose']
    checkpoints = score['baseline_checkpoints']
    for mode, result in score['comparisons'].items():
        matched = result.get('cm_status') == 'available'
        cm_reason = None if matched or mode == 'native' else 'matching CM EMA unavailable at frozen optimizer step'
        if not result['paired']:
            yield from pair_rows(0, seed, arm, mode, result, cm_reason=cm_reason)
            continue
        steps = (checkpoints['Transfer']['metadata']['best_step'],
                 score['intervention_checkpoints'][mode]['metadata']['best_step'],
                 checkpoints['CM']['metadata']['best_step'] if matched else None)
        if mode == 'common' and (steps[0] != steps[1] or matched and steps[0] != steps[2]):
            raise ValueError('equal-step validation comparison contains different optimizer steps')
        yield from pair_rows(0, seed, arm, mode, result, summaries['Transfer'], summaries[mode],
                             summaries['CM'] if matched else None, steps, cm_reason)


def confirmation_pairs(result, score_by_seed):
    for arm, arm_result in result['arms'].items():
        for mode, mode_result in arm_result['comparisons'].items():
            for key, comparison in mode_result['per_seed'].items():
                seed = int(key)
                cm_reason = comparison['cm_reason']
                if not comparison['contrasts']:
                    yield from pair_rows(1, seed, arm, mode, comparison, cm_reason=cm_reason)
                    continue
                score = score_by_seed[key]
                summaries = score['summaries']
                steps = (score['fits']['Transfer']['steps'][mode], score['fits'][arm]['steps'][mode],
                         score['cm_native_step'] if cm_reason is None else None)
                if mode == 'common' and (steps[0] != steps[1] or steps[2] not in (None, steps[0])):
                    raise ValueError('equal-step confirmation comparison contains different optimizer steps')
                yield from pair_rows(1, seed, arm, mode, comparison, summaries[f'Transfer/{mode}'],
                                     summaries[f'{arm}/{mode}'],
                                     summaries['CM/native'] if cm_reason is None else None, steps, cm_reason)


def audit_evidence(output):
    source = output / 'source-audit'
    lineage = read_json(source / 'lineage.json')
    if lineage['schema'] != 'step06-source-audit-v1' or lineage['descriptive_only'] is not True:
        raise ValueError('unexpected source-domain audit lineage')
    expected = {'exposure.csv', 'delta_panel.csv', 'exact_smiles_panel.csv'}
    if set(lineage['outputs']) != expected:
        raise ValueError('unexpected source-domain audit tables')
    for name, digest in lineage['outputs'].items():
        bind(source / name, digest)
        if (output / 'source-audit-replay').is_dir():
            bind(output / 'source-audit-replay' / name, digest)
    replayed = (output / 'source-audit-replay').is_dir()
    if replayed:
        bind(output / 'source-audit-replay/lineage.json', sha256(source / 'lineage.json'))
    panel = rows(source / 'delta_panel.csv')
    exact = rows(source / 'exact_smiles_panel.csv')
    return {'independent_byte_replay': replayed, 'shared_parents': len(lineage['shared_training_parents']),
            'exact_smiles_parents': len(lineage['exact_parent_smiles_training_parents']),
            'tahoe_cell_lines': len({row['tahoe_cell_line_id'] for row in panel if row['tahoe_cell_line_id']}),
            'delta_panel_rows': len(panel), 'exact_smiles_panel_rows': len(exact),
            'brackets_um': {'0.2': [0.05, 0.5], '1.0': [0.5, 5.0], '10.0': []},
            'interpretation': 'Descriptive panel. Assay, normalization, dose, exposure time, cell background '
                              'and coverage vary together. Tahoe per-sample exposure time is unavailable.'}


def resources(output, replay_record):
    root = output / 'resources'
    ids = (replay_record['jobs'] if replay_record else
           [path.parent.parent.name for path in sorted(root.glob('*/scheduler/status.json'))
            if not path.parent.parent.name.startswith('ticket96')
            or path.parent.parent.name == 'ticket96-repo-verify'])
    ledger = {}
    for job in ids:
        slot = root / job
        scheduler = slot / 'scheduler'
        status = read_json(scheduler / 'status.json')
        entry = status['jobs'][job]
        if entry['status'] not in ('completed', 'failed', 'skipped', 'timeout', 'excluded'):
            raise ValueError(f'scheduler run {job} requires a terminal disposition')
        attempts = [{**item, 'wall_seconds': (item['finished_at'] - item['started_at']
                                             if 'finished_at' in item and 'started_at' in item else None)}
                    for item in entry['attempts']]
        files = [slot / 'queue.json', scheduler / 'status.json', scheduler / 'telemetry.jsonl',
                 *(Path(item['log']) for item in attempts if item.get('log'))]
        disposition = slot / 'disposition.json'
        if disposition.is_file():
            files.append(disposition)
        ledger[job] = {'status': entry['status'], 'identity': entry['identity'], 'attempts': attempts,
                       'disposition': read_json(disposition) if disposition.is_file() else None,
                       'files': [bind(path) for path in files]}
    probe = root / 'between-phases.tsv'
    return {'jobs': ledger, 'between_phases': bind(probe) if probe.is_file() else None}


def evidence_inventory(output, documents):
    """Verify and list referenced files, including raw predictions and EMA checkpoints."""
    inventory = {}

    def add(path, digest=None):
        path = Path(path)
        key = str(path)
        if key not in inventory:
            inventory[key] = bind(path, digest)['sha256']
        elif digest is not None and digest != inventory[key]:
            raise ValueError(f'conflicting digests for {path}')

    def visit(value):
        if isinstance(value, dict):
            if {'path', 'sha256'} <= value.keys():
                add(value['path'], value['sha256'])
            for item in value.values():
                visit(item)
        elif isinstance(value, list):
            for item in value:
                visit(item)

    for document in documents:
        add(document)
        if document.suffix == '.json':
            visit(read_json(document))
    for path in sorted(output.rglob('prediction.json')):
        record = read_json(path)
        if record.get('schema') != baseline.PREDICTION_SCHEMA:
            raise ValueError(f'unrecognized prediction seal at {path}')
        add(path)
        visit(record)
        add(path.with_name('cases.json'), record['cases_sha256'])
        add(path.with_name('prediction.npy'), record['prediction_sha256'])
    for path in sorted((output / 'checkpoints').rglob('*')):
        if path.is_file() and path.suffix in ('.pt', '.json', '.csv'):
            if path.suffix == '.pt' and (not path.relative_to(output / 'checkpoints').parts[0][:2].isdigit()
                                         or not path.name[:2].isdigit()):
                raise ValueError(f'unnumbered local checkpoint {path}')
            add(path)
    for directory in sorted(output.glob('source-audit*')):
        for path in sorted(directory.iterdir()):
            if path.suffix in ('.csv', '.json'):
                add(path)
    return dict(sorted(inventory.items()))


def render(pairs, exposure, destination):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    figure = Figure(figsize=(9, 4.8), layout='constrained')
    FigureCanvasAgg(figure)
    ax = figure.subplots()
    visible = [row for row in pairs if row['metric'] == 'mean_ap100' and row['gain'] is not None]
    if visible:
        labels = {'dose_reset': 'Dose reset', 'learning_rate': 'Learning rate'}
        cohorts = {}
        for row in visible:
            key = (row['fold'], row['seed'], row['arm'], row['mode'])
            cohorts.setdefault(key, []).append(row)
        styles = {'loss-selected': 'o-', 'native': 's--', 'common': '^-'}
        dose_order = list(dict.fromkeys(row['dose_um'] for row in visible))
        for (fold, seed, arm, mode), values in cohorts.items():
            ax.plot([dose_order.index(row['dose_um']) for row in values], [row['gain'] for row in values],
                    styles[mode], label=f'{labels[arm]} · fold {fold} / seed {seed} · {mode}')
        ax.set_xticks(range(len(dose_order)), ['Pooled' if dose == 'pooled' else f'{float(dose):g} µM'
                                              for dose in dose_order])
        if any(row['fold'] == 0 and row['seed'] == 17 and row['dose_um'] == 'pooled' for row in visible):
            ax.plot(0, float(THRESHOLD), 'x', color='black', label='Fold 0 pooled promotion gate')
        ax.legend(frameon=False, fontsize=8)
    else:
        ax.text(.5, .5, 'No scored intervention pairs', transform=ax.transAxes, ha='center')
    ax.axhline(0, color='#64748b', linewidth=.8)
    ax.set(ylabel='Intervention − original Transfer (mean AP@100)', title='Paired PerturbDiff effects by dose')
    figure.savefig(destination / FIGURES[0], dpi=144)

    figure = Figure(figsize=(8.5, 4.1), layout='constrained')
    FigureCanvasAgg(figure)
    ax = figure.subplots()
    roles = ('train', 'validation', 'test')
    cohorts = {(row['source_release'], row['tahoe_role']) for row in exposure}
    counts = {(row['source_release'], row['tahoe_role'], row['tenv_role']): int(row['shared_parent_count'])
              for row in exposure}
    for index, (release, role) in enumerate(sorted(cohorts, reverse=True)):
        offset = (index - (len(cohorts) - 1) / 2) * .19
        ax.bar([i + offset for i in range(3)], [counts[release, role, key] for key in roles], width=.18,
               label=f'{release.replace("_", " ").capitalize()} Tahoe {role}')
    ax.set_xticks(range(3), [key.capitalize() for key in roles])
    ax.set(xlabel='Tenv fold 0 role', ylabel='Shared parent-connectivity blocks',
           title='Tahoe fitting exposure before and after global exclusion')
    ax.legend(frameon=False, fontsize=8)
    figure.savefig(destination / FIGURES[1], dpi=144)


PAIR_FIELDS = (
    'fold', 'seed', 'arm', 'mode', 'dose_um', 'metric', 'status', 'reason', 'group_count',
    'transfer_step', 'intervention_step', 'cm_step', 'cm_reason',
    'transfer', 'intervention', 'cm', 'transfer_valid_groups', 'intervention_valid_groups',
    'cm_valid_groups', 'gain', 'cm_gap', 'residual_gap', 'recovery',
)
RECORD_FILES = ('decision.json', 'evidence.json', 'resources.json', 'per_run.csv', 'per_dose.csv',
                'source_exposure.csv')


def collect(runtime, output, replay):
    path = output / replication.PLAN
    digest = sha256(path)
    plan = replication.load_plan(runtime, output, digest)
    frozen = {arm: confirmation.verified_screen(output, output / ROOTS[arm], arm)
              for arm in plan['screens']}
    if frozen != plan['screens']:
        raise ValueError('replication plan differs from the verified screening decisions')
    result_path = output / replication.RESULT
    result = read_json(result_path)
    score_by_seed = {key: read_json(check_binding(item)['path'])
                     for key, item in result.get('scores', {}).items()}
    if replay:
        for key in score_by_seed:
            replication.score(output, plan, digest, int(key))
    unavailable = {key: (item['status'], item['reason'])
                   for key, item in result.get('seed_dispositions', {}).items() if key not in score_by_seed}
    replication.finalize(output, plan, digest, unavailable)

    pairs = []
    documents = [path, result_path, output / confirmation.PLAN_PATH]
    seed_runs = []
    seen = set()
    for arm, slot, score, binding in validation_scores(output, frozen):
        seen.add((arm, score['seed']))
        if arm == 'baseline':
            seed_runs.append({'fold': 0, 'seed': score['seed'], 'arm': arm,
                              'status': 'completed', 'reason': None})
        documents.extend([Path(binding['path']), slot / 'scores/scores.json'])
        if replay:
            replay_validation(runtime, output, arm, slot, score, binding)
        pairs.extend(validation_pairs(arm, score))
    for arm, screen in frozen.items():
        for seed in SEEDS:
            key = str(seed)
            if (arm, seed) in seen:
                status = 'incomplete' if key in screen['incomplete_seeds'] else 'completed'
                reason = screen['incomplete_seeds'].get(key)
            elif key in screen['incomplete_seeds']:
                status, reason = 'incomplete', screen['incomplete_seeds'][key]
            else:
                status, reason = 'skipped', screen['reason']
                if (output / ROOTS[arm] / f'fold-0-seed-{seed}/plan.json').is_file():
                    raise ValueError(f'unexpected execution plan for skipped {arm} seed {seed}')
            seed_runs.append({'fold': 0, 'seed': seed, 'arm': arm, 'status': status, 'reason': reason})
    for seed in SEEDS:
        if result['status'] == 'skipped':
            seed_runs.append({'fold': 1, 'seed': seed, 'arm': 'replication',
                              'status': 'skipped', 'reason': result['reason']})
        else:
            seed_runs.append({'fold': 1, 'seed': seed, 'arm': 'replication',
                              **result['seed_dispositions'][str(seed)]})
    pairs.extend(confirmation_pairs(result, score_by_seed))
    summary = {'schema': 'step06-publication-v1', 'promotion_threshold': float(THRESHOLD),
               'major_contribution_minimum_recovery': 0.5,
               'arms': {arm: diagnosis(screen, result['arms'].get(arm)) for arm, screen in frozen.items()},
               'replication': {'status': result['status'], 'reason': result['reason'],
                               'historical_scoring': result['historical_scoring']},
               'runs': seed_runs, 'source_audit': audit_evidence(output),
               'pooled_validation': [row for row in pairs if row['fold'] == 0
                                     and row['dose_um'] == 'pooled' and row['metric'] == 'mean_ap100'],
               'interpretation': 'Residual CM gaps are paired observations. Each intervention was tested '
                                 'separately. The source-domain audit does not enter attribution.'}
    documents.extend(Path(item['path']) for item in result.get('scores', {}).values())
    for screen in frozen.values():
        documents.append(Path(screen['decision']['path']))
    documents.extend(sorted(output.glob('baseline/*/scores/per_*.csv')))
    documents.extend(sorted(output.glob('dose_reset/*/scores/per_*.csv')))
    documents.extend(sorted(output.glob('learning-rate/*/scores/per_*.csv')))
    return summary, pairs, evidence_inventory(output, documents)


def write_record(runtime, output, record, replay):
    summary, pairs, inventory = collect(runtime, output, replay)
    exposure = rows(output / 'source-audit/exposure.csv')
    ledger = resources(output, read_json(record / 'resources.json') if replay else None)
    with tempfile.TemporaryDirectory(dir=output, prefix='.publication-') as temporary:
        staging = Path(temporary)
        render(pairs, exposure, staging)
        write_csv(staging / 'per_run.csv', [row for row in pairs if row['dose_um'] == 'pooled'],
                  fields=PAIR_FIELDS)
        write_csv(staging / 'per_dose.csv', [row for row in pairs if row['dose_um'] != 'pooled'],
                  fields=PAIR_FIELDS)
        (staging / 'source_exposure.csv').write_bytes((output / 'source-audit/exposure.csv').read_bytes())
        write_json(staging / 'decision.json', summary)
        write_json(staging / 'resources.json', ledger)
        generated = (*FIGURES, *(name for name in RECORD_FILES if name != 'evidence.json'))
        write_json(staging / 'evidence.json', {
            'inputs': inventory,
            'figures': {f'output/{name}': sha256(staging / name) for name in FIGURES},
            'record': {name: sha256(staging / name) for name in generated if name not in FIGURES},
        })
        destinations = {name: record / name for name in RECORD_FILES}
        destinations.update({name: output / name for name in FIGURES})
        with slot_lock(output / '.publication-lock'):
            for name, path in destinations.items():
                if replay or path.exists():
                    bind(path, sha256(staging / name))
            if not replay:
                record.mkdir(parents=True, exist_ok=True)
                for name, path in destinations.items():
                    if not path.exists():
                        (staging / name).replace(path)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('publish', 'replay'))
    parser.add_argument('output', type=Path)
    parser.add_argument('--runtime', type=Path, default=baseline.REPO)
    parser.add_argument('--record', required=True, type=Path)
    args = parser.parse_args()
    summary = write_record(args.runtime.resolve(), args.output.resolve(), args.record.resolve(),
                           args.phase == 'replay')
    print(f"{args.phase}: fold 1 {summary['replication']['status']}; "
          + ', '.join(f"{arm} {item['attribution_status']}" for arm, item in summary['arms'].items()))


if __name__ == '__main__':
    main()
