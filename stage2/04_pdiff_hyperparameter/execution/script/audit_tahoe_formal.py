"""Audit Tahoe held-out coverage, defined counts, and exact scoring replay."""
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
import sys

METRICS = ('AUPRC_positive@100', 'AUPRC_negative@100', 'Recall_positive@P80',
           'Recall_negative@P80', 'PCC', 'Spearman', 'RMSE', 'R2')
MODELS = {'perturbdiff': ('17', '42', '2026'), 'no_change': ('',), 'delta_mean': ('',)}
FAILURES = {'failed', 'incomplete', 'interrupted', 'timeout', 'cancelled'}


def read(path):
    return json.loads(path.read_text())


def rows(path):
    with path.open(newline='') as stream:
        yield from csv.DictReader(stream)


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def check(condition, explanation):
    if not condition:
        raise ValueError(explanation)


def verify(path, expected):
    check(digest(path) == expected, f'digest mismatch: {path}')


def audit(root, plan_sha):
    plan_file, outcomes_file = root / 'plan.json', root / 'outcomes.json'
    verify(plan_file, plan_sha)
    plan, outcomes = read(plan_file), read(outcomes_file)
    check(plan['dataset'] == 'tahoe' and len(plan['folds']) == 3, 'expected frozen Tahoe folds')
    check(outcomes['plan_sha256'] == plan_sha, 'outcomes do not match the frozen plan')
    primary, replay = root / 'evaluate', root / 'replay'
    record_file = primary / 'evaluation.json'
    record_sha = digest(record_file)
    verify(replay / record_file.name, record_sha)
    record = read(record_file)
    outcomes_sha = digest(outcomes_file)
    check(record['schema'] == 'step04-evaluation-v1' and record['stage'] == 'formal'
          and record['datasets'] == ['tahoe'] and record['plan_sha256'] == plan_sha
          and record['outcomes_sha256'] == outcomes_sha
          and record['tenv_test_memberships_previously_scored_in_step03'] is False,
          'unexpected scoring scope or identities')
    expected_files = set(record['files']) | {record_file.name}
    check({'per_group.csv', 'per_run.csv', 'summary.csv'} <= expected_files, 'missing scored table')
    for directory in (primary, replay):
        check({str(path.relative_to(directory)) for path in directory.rglob('*') if path.is_file()}
              == expected_files, f'scoring file coverage mismatch: {directory}')
        for name, checksum in record['files'].items():
            verify(directory / name, checksum)

    groups = {}
    for fold, entry in enumerate(plan['folds']):
        request_file = root / entry['request']
        verify(request_file, entry['request_sha256'])
        request = read(request_file)
        cases = Path(request['predictor_files']['test_cases'])
        verify(cases, request['file_sha256']['test_cases'])
        groups[str(fold)] = Counter(row['group_id'] for row in rows(cases))
    expected = {(fold, model, seed) for fold in groups
                for model, seeds in MODELS.items() for seed in seeds}
    model_keys = {key for key in expected if key[1] == 'perturbdiff'}
    terminal = {}
    for row in outcomes['outcomes']:
        key = (str(row['fold']), 'perturbdiff', str(row['seed']))
        check(row['dataset'] == 'tahoe' and key in model_keys and key not in terminal,
              f'invalid terminal slot {key}')
        if row['status'] == 'completed':
            verify(root / row['prediction'], row['prediction_manifest_sha256'])
        else:
            check(row['status'] in FAILURES and row['reason'] and 'prediction' not in row,
                  f'failed slot requires a terminal cause: {key}')
        terminal[key] = row
    check(set(terminal) == model_keys, 'formal terminal inventory must cover all nine slots')

    def slot(row):
        check(row['dataset'] == 'tahoe', 'unexpected dataset in score row')
        return row['fold'], row['model'], row['seed']

    runs, remaining, defined, undefined = {}, {}, defaultdict(Counter), defaultdict(Counter)
    for row in rows(primary / 'per_run.csv'):
        key = slot(row)
        check(key in expected and key not in runs, f'invalid run slot {key}')
        terminal_row = terminal.get(key, {'status': 'completed', 'reason': ''})
        check((row['status'], row['reason']) == (terminal_row['status'], terminal_row.get('reason', '')),
              f'run outcome differs from terminal inventory: {key}')
        runs[key] = row
        remaining[key] = groups[key[0]].copy() if row['status'] == 'completed' else {}
        check(int(row['group_count']) == len(remaining[key]), f'run group count differs: {key}')
    check(set(runs) == expected, 'run inventory must cover nine model and six reference slots')
    for row in rows(primary / 'per_group.csv'):
        key = slot(row)
        check(row['group_id'] in remaining[key] and remaining[key].pop(row['group_id']) == int(row['n_wells']),
              f'invalid, mismatched, or duplicate scored group: {key} {row["group_id"]}')
        for metric in METRICS:
            value, reason = row[metric], row[f'{metric}_reason']
            check(bool(value) != bool(reason), f'score or undefined reason required: {key} {metric}')
            if value:
                check(math.isfinite(float(value)), f'nonfinite score: {key} {metric}')
                defined[key][metric] += 1
            else:
                undefined[key[1], metric][reason] += 1
    check(not any(remaining.values()), 'a completed run lacks frozen test groups')
    for key, row in runs.items():
        for metric in METRICS:
            count = defined[key][metric]
            check(int(row[f'{metric}_valid_groups']) == count and bool(row[metric]) == bool(count)
                  and (not count or math.isfinite(float(row[metric]))),
                  f'run metric coverage differs: {key} {metric}')

    summary = {}
    for row in rows(primary / 'summary.csv'):
        model = row['model']
        check(row['dataset'] == 'tahoe' and model in MODELS and model not in summary,
              f'invalid summary model: {model}')
        model_runs = [run for key, run in runs.items() if key[1] == model]
        completed = sum(run['status'] == 'completed' for run in model_runs)
        check(int(row['expected_runs']) == len(model_runs) and int(row['completed_runs']) == completed,
              f'summary run counts differ: {model}')
        for metric in METRICS:
            count = sum(bool(run[metric]) for run in model_runs)
            check(int(row[f'{metric}_valid_runs']) == count
                  and bool(row[f'{metric}_mean']) == bool(count)
                  and bool(row[f'{metric}_sample_sd']) == (count > 1),
                  f'summary metric coverage differs: {model} {metric}')
        summary[model] = completed
    check(summary.keys() == MODELS.keys(), 'summary model coverage differs')
    failed = [row for row in outcomes['outcomes'] if row['status'] != 'completed']
    check(record['complete'] is (not failed), 'completion marker differs from terminal outcomes')
    result = {'plan_sha256': plan_sha, 'outcomes_sha256': outcomes_sha, 'evaluation_sha256': record_sha,
              'replay_identical': True, 'complete': not failed, 'failed_slots': failed,
              'model_slots_scored': summary['perturbdiff'], 'reference_slots_scored': 6,
              'heldout_groups_by_fold': {fold: len(members) for fold, members in groups.items()},
              'heldout_wells_by_fold': {fold: sum(members.values()) for fold, members in groups.items()},
              'metrics': list(METRICS), 'scored_files': record['files'],
              'defined_groups_by_model': {
                  model: {metric: sum(defined[key][metric] for key in runs if key[1] == model)
                          for metric in METRICS}
                  for model in MODELS},
              'undefined_groups_by_model': {
                  model: {metric: dict(undefined[model, metric]) for metric in METRICS}
                  for model in MODELS}}
    pending = root / '.audit.pending'
    pending.write_text(json.dumps(result, indent=2) + '\n')
    pending.replace(root / 'audit.json')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    audit(Path(sys.argv[1]), sys.argv[2])
