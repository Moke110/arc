"""Recompute descriptive benchmark summaries from persisted predictions, without fitting."""
import hashlib
import json
from pathlib import Path

import numpy as np
from metrics import CLASS_ORDER, METRIC_CONVENTION, evaluate_predictions


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def _key(identity):
    return tuple(identity[k] for k in ('condition', 'phenotype', 'setting', 'fold', 'seed'))


def _counts(rows):
    return dict(measurements=len(rows), parents=len({r['parent_id'] for r in rows}),
                combinations=len({(r['parent_id'], r['dose_um']) for r in rows}),
                source_labels=len({r['source_label'] for r in rows if 'source_label' in r}),
                classes={str(c): sum(r['direction'] == c for r in rows) for c in CLASS_ORDER})


def _aggregate(items, axes):
    result = dict(n_runs=len(items), axes=axes, standard_deviation='sample; ddof=1; null for n<2')
    for metric in ('macro_f1', 'log_loss', 'precision', 'recall', 'f1', 'confusion_matrix'):
        values = np.asarray([x[metric] for x in items], dtype=float)
        result[metric] = dict(mean=values.mean(axis=0).tolist(),
                              std=values.std(axis=0, ddof=1).tolist() if len(items) > 1 else None)
    return result


def summarize(prepared, runs, require_complete=False):
    """Summarize run directories or ledger entries with identity/status/run_path.

    Missing entries are pending. Failed/pending ledger entries carry no predictions.
    Completed entries must reference an immutable run directory. Smoke and formal
    aggregates are separate; only complete formal coverage completes a benchmark.
    """
    if not isinstance(prepared, dict):
        prepared = json.loads(Path(prepared).read_text())
    rows, splits = prepared['rows'], prepared['splits']
    for field, content in [('input_fingerprint', rows), ('split_fingerprint', splits)]:
        if prepared[field] != _digest(content):
            raise ValueError('Prepared fingerprint mismatch: ' + field)
    if len({r['measurement_id'] for r in rows}) != len(rows):
        raise ValueError('Duplicate prepared measurement identity')
    expected = {_key(i): i for i in prepared['expected_runs']}
    if len(expected) != len(prepared['expected_runs']):
        raise ValueError('Duplicate expected run identity')
    seen, outcomes, validated = set(), [], []
    for entry in runs:
        if isinstance(entry, dict):
            identity, status = entry['identity'], entry['status']
            path = entry.get('run_path')
            if status not in ('pending', 'failed', 'completed'):
                raise ValueError('Unknown run status')
        else:
            path, status, identity = entry, 'completed', None
        records = None
        if status == 'completed':
            if path is None:
                raise ValueError('Completed ledger entry needs run_path')
            path = Path(path)
            path = path.parent if path.name == 'run.json' else path
            manifest = json.loads((path / 'run.json').read_text())
            if identity is not None and identity != manifest['identity']:
                raise ValueError('Ledger/run identity mismatch')
            identity = manifest['identity']
            if manifest['status'] != 'completed' or manifest['mode'] not in ('smoke', 'formal'):
                raise ValueError('Invalid completed run mode/status')
            if (manifest['mode'] == 'formal') != (manifest.get('smoke_epochs') is None):
                raise ValueError('Inconsistent smoke/formal mode')
            for field in ('input_fingerprint', 'split_fingerprint', 'features'):
                if manifest[field] != prepared[field]:
                    raise ValueError('Incompatible run ' + field)
            prediction_path = path / 'predictions.json'
            if hashlib.sha256(prediction_path.read_bytes()).hexdigest() != manifest['artifact_hashes']['predictions.json']:
                raise ValueError('Prediction artifact hash mismatch')
            records = json.loads(prediction_path.read_text())
        key = _key(identity)
        if key not in expected or expected[key] != identity:
            raise ValueError('Unexpected run identity')
        if key in seen:
            raise ValueError('Duplicate run identity')
        seen.add(key)
        outcome = dict(identity=identity, status=status)
        if records is not None:
            split = next(s for s in splits if s['setting'] == identity['setting'] and s['fold'] == identity['fold'])
            if len(split['roles']) != len(rows):
                raise ValueError('Partition length mismatch')
            test = [r for r, role in zip(rows, split['roles']) if role == 'test'
                    and (r['condition'], r['phenotype']) == (identity['condition'], identity['phenotype'])]
            ids = [r['measurement_id'] for r in records]
            if len(set(ids)) != len(ids):
                raise ValueError('Duplicate prediction identity')
            actual = [{k: v for k, v in r.items() if k not in ('run_identity', 'probabilities', 'prediction')} for r in records]
            if sorted(actual, key=lambda r: r['measurement_id']) != sorted(test, key=lambda r: r['measurement_id']):
                raise ValueError('Prediction truth/provenance/identities differ from prepared test partition')
            if any(r['run_identity'] != identity for r in records):
                raise ValueError('Prediction run identity mismatch')
            metrics = evaluate_predictions([r['direction'] for r in records], [r['probabilities'] for r in records], [r['prediction'] for r in records])
            outcome.update(mode=manifest['mode'], metrics=metrics)
            validated.append((identity, manifest['mode'], records, metrics, split))
        outcomes.append(outcome)
    outcomes.extend(dict(identity=i, status='pending') for k, i in expected.items() if k not in seen)
    outcomes.sort(key=lambda x: _key(x['identity']))
    formal = sum(x.get('mode') == 'formal' for x in outcomes)
    requested_scope_complete = bool(expected) and all(x['status'] == 'completed' for x in outcomes)
    formal_scope_complete = formal == len(expected) and bool(expected)
    # ADR 0045 defines the full benchmark independently of a prepared subset.
    canonical = {(condition, 'cell_viability', setting, fold, seed)
                 for condition in ('baseline', 'H-R')
                 for setting, folds in (('cold_drug', 3), ('cold_dose', 4), ('warm_start', 3))
                 for fold in range(folds) for seed in (17, 42, 2026)}
    complete = formal_scope_complete and set(expected) == canonical
    if require_complete and not complete:
        raise ValueError('Incomplete formal run coverage')
    targets = []
    for condition, phenotype in sorted({(r['condition'], r['phenotype']) for r in rows}):
        selected = [r for r in rows if (r['condition'], r['phenotype']) == (condition, phenotype)]
        roles = [dict(setting=s['setting'], fold=s['fold'], roles={role: _counts([r for r, assigned in zip(rows, s['roles']) if assigned == role and (r['condition'], r['phenotype']) == (condition, phenotype)]) for role in ('train', 'validation', 'test')}) for s in splits]
        targets.append(dict(condition=condition, phenotype=phenotype, coverage=_counts(selected), splits=roles))
    aggregates, pooled = [], []
    groups = sorted({(i['condition'], i['phenotype'], i['setting'], mode) for i, mode, *_ in validated})
    for condition, phenotype, setting, mode in groups:
        group = [x for x in validated if (x[0]['condition'], x[0]['phenotype'], x[0]['setting'], x[1]) == (condition, phenotype, setting, mode)]
        base = dict(condition=condition, phenotype=phenotype, setting=setting, mode=mode)
        if setting == 'cold_dose':
            for fold in sorted({x[0]['fold'] for x in group}):
                dose_group = [x for x in group if x[0]['fold'] == fold]
                split = dose_group[0][4]
                aggregates.append(dict(base, scope='held_dose', held_dose=split.get('held_dose_um', split.get('held_dose')), generalization_type=split.get('dose_generalization', split.get('generalization_type')), **_aggregate([x[3] for x in dose_group], ['seed'])))
        aggregates.append(dict(base, scope='overall_fold_mean', **_aggregate([x[3] for x in group], ['fold', 'seed'])))
        if setting == 'cold_drug':
            for seed in sorted({x[0]['seed'] for x in group}):
                parts = [x for x in group if x[0]['seed'] == seed]
                if {x[0]['fold'] for x in parts} != {0, 1, 2}:
                    continue
                records = [r for x in parts for r in x[2]]
                ids = [r['measurement_id'] for r in records]
                target_ids = {r['measurement_id'] for r in rows if (r['condition'], r['phenotype']) == (condition, phenotype)}
                if len(ids) != len(set(ids)) or set(ids) != target_ids:
                    raise ValueError('Cold-drug pooled OOF must cover each target measurement exactly once')
                metrics = evaluate_predictions([r['direction'] for r in records], [r['probabilities'] for r in records], [r['prediction'] for r in records])
                pooled.append(dict(base, seed=seed, metrics=metrics))
    pooled_aggregates = []
    for condition, phenotype, setting, mode in groups:
        items = [x['metrics'] for x in pooled if (x['condition'], x['phenotype'], x['setting'], x['mode']) == (condition, phenotype, setting, mode)]
        if items:
            pooled_aggregates.append(dict(condition=condition, phenotype=phenotype, setting=setting, mode=mode, **_aggregate(items, ['seed; each metric pools three outer folds'])))
    return dict(schema_version=1, complete=complete, requested_scope_complete=requested_scope_complete,
                formal_scope_complete=formal_scope_complete, expected_runs=len(expected), completed_formal=formal,
                status_counts={s: sum(x['status'] == s for x in outcomes) for s in ('completed', 'pending', 'failed')},
                metric_convention=METRIC_CONVENTION, targets=targets, runs=outcomes, aggregates=aggregates,
                cold_drug_pooled_oof=pooled, cold_drug_pooled_oof_aggregates=pooled_aggregates,
                source_coverage=prepared.get('source_coverage', prepared.get('provenance_metadata')))


def render_report(summary):
    """Render compact scientific tables; JSON artifacts retain exhaustive detail."""
    lines = ['# Morgan direction benchmark', '',
             'Status: ' + ('complete formal benchmark' if summary['complete'] else 'INCOMPLETE formal benchmark'), '',
             'Requested scope: ' + ('complete' if summary.get('requested_scope_complete') else 'incomplete')
             + '; formal scope: ' + ('complete' if summary.get('formal_scope_complete') else 'incomplete') + '.', '',
             'Machine-readable detail: [summary.json](summary.json) and [run ledger](benchmark.json).', '']

    def table(title, headers, rows):
        lines.extend(['## ' + title, '', '| ' + ' | '.join(headers) + ' |',
                      '| ' + ' | '.join('---' for _ in headers) + ' |'])
        for row in rows:
            lines.append('| ' + ' | '.join(str(x).replace('|', r'\|').replace('\n', ' ') for x in row) + ' |')
        if not rows:
            lines.append('| ' + ' | '.join(['No results available'] + ['—'] * (len(headers) - 1)) + ' |')
        lines.append('')

    def number(value):
        return '—' if value is None else f'{value:.4f}'

    def mean_sd(item, metric, index=None):
        value = item[metric]
        mean, sd = value['mean'], value['std']
        if index is not None:
            mean, sd = mean[index], None if sd is None else sd[index]
        return number(mean) + ' ± ' + number(sd)

    def target(item):
        return item['condition'] + '/' + item['phenotype']

    counts = summary['status_counts']
    smoke = sum(r.get('mode') == 'smoke' for r in summary['runs'])
    table('Run status', ['Expected', 'Completed formal', 'Completed smoke', 'Pending', 'Failed'],
          [[summary['expected_runs'], summary['completed_formal'], smoke, counts['pending'], counts['failed']]])
    lines.extend(['Smoke and formal results are reported separately. Smoke fits do not complete the formal benchmark.', ''])
    source = summary.get('source_coverage') or {}
    coverage = source.get('coverage', source)
    source_rows = [[condition, c.get('source_values', '—'), c.get('admitted_values', '—'),
                    c.get('pending_values', '—'), c.get('resolved_labels', '—')]
                   for condition, c in sorted(coverage.items()) if isinstance(c, dict)]
    table('Source coverage', ['Condition', 'Source measurements', 'Admitted', 'Pending identity', 'Resolved labels'], source_rows)
    lines.extend(['Pending identities are excluded from admitted measurements; — means coverage was not supplied.', ''])
    table('Admitted coverage and class distribution',
          ['Target', 'Measurements', 'Parents', 'Parent–dose combinations', 'Source labels', '-1', '0', '+1'],
          [[target(t), *[t['coverage'][k] for k in ('measurements', 'parents', 'combinations', 'source_labels')],
            *[t['coverage']['classes'][str(c)] for c in CLASS_ORDER]] for t in summary['targets']])
    table('Split role counts', ['Target', 'Setting', 'Fold', 'Role', 'Measurements', 'Parents', 'Combinations'],
          [[target(t), s['setting'], s['fold'], role, c['measurements'], c['parents'], c['combinations']]
           for t in summary['targets'] for s in t['splits'] for role, c in s['roles'].items()])
    held = [a for a in summary['aggregates'] if a['scope'] == 'held_dose']
    table('Cold-dose held levels', ['Target', 'Mode', 'Held dose (µM)', 'Generalization', 'n seeds', 'Macro-F1 mean ± SD', 'Log loss mean ± SD'],
          [[target(a), a['mode'], a['held_dose'], a['generalization_type'], a['n_runs'], mean_sd(a, 'macro_f1'), mean_sd(a, 'log_loss')] for a in held])
    overall = [a for a in summary['aggregates'] if a['scope'] == 'overall_fold_mean']
    table('Unweighted fold × seed means', ['Target', 'Setting', 'Mode', 'n runs', 'Macro-F1 mean ± SD', 'Log loss mean ± SD'],
          [[target(a), a['setting'], a['mode'], a['n_runs'], mean_sd(a, 'macro_f1'), mean_sd(a, 'log_loss')] for a in overall])
    table('Cold-drug pooled OOF per seed', ['Target', 'Mode', 'Seed', 'Macro-F1', 'Log loss'],
          [[target(a), a['mode'], a['seed'], number(a['metrics']['macro_f1']), number(a['metrics']['log_loss'])] for a in summary['cold_drug_pooled_oof']])
    table('Cold-drug pooled OOF across seeds', ['Target', 'Mode', 'n seeds', 'Macro-F1 mean ± SD', 'Log loss mean ± SD'],
          [[target(a), a['mode'], a['n_runs'], mean_sd(a, 'macro_f1'), mean_sd(a, 'log_loss')] for a in summary['cold_drug_pooled_oof_aggregates']])
    table('Per-class precision and recall (fold × seed)', ['Target', 'Setting', 'Mode', 'Class', 'n runs', 'Precision mean ± SD', 'Recall mean ± SD'],
          [[target(a), a['setting'], a['mode'], c, a['n_runs'], mean_sd(a, 'precision', j), mean_sd(a, 'recall', j)]
           for a in overall for j, c in enumerate(CLASS_ORDER)])
    table('Every expected run', ['Target', 'Setting', 'Fold', 'Seed', 'Status', 'Mode'],
          [[target(r['identity']), r['identity']['setting'], r['identity']['fold'], r['identity']['seed'], r['status'], r.get('mode', '—')] for r in summary['runs']])
    lines.extend(['## Interpretation and limitations', '',
        'Dose normalization uses training-only log10(µM) mean and standard deviation. Directions describe magnitude relative to each source arm’s assumed vehicle=100: -1 below 80, 0 from 80 through 120 inclusive, +1 above 120. Class 0 is within-band, not demonstrated equivalence.', '',
        'Source measurements are not established independent biological replicates. Training seeds describe optimization variation, not biological replication. Direction does not establish therapeutic benefit or clinical efficacy. Softmax does not establish probability calibration.', '',
        'All metrics retain classes [-1, 0, +1], with zero division set to zero. SD is sample SD (ddof=1); — denotes unavailable SD for fewer than two contributors. Fold means are unweighted. Pooled OOF metrics concatenate all three outer test folds within a seed and are distinct from fold means. Cold-dose interior levels are interpolation; boundary levels are extrapolation.', '',
        'Confusion matrices, per-run class metrics, numerical log-loss conventions and split details are retained in [summary.json](summary.json). Saved predictions retain measurement provenance; checkpoint selection uses validation data only.', ''])
    return '\n'.join(lines)
