"""Accept replayed formal origins and publish candidate scientific tables."""
from collections import defaultdict
import math
from pathlib import Path
from statistics import mean, stdev
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from auxiliary import audit_coverage
from formal import SCENARIOS, read_json, run_path, sha256, verify, write_json
from tsg_round1.inputs import Inputs, SOURCE_POPULATIONS
from tsg_round1.matrix import METHODS, registration_path, selection_path
from tsg_round1.protocol import METRICS, summary

IDENTITIES = ['profile_id', 'drug_parent', 'dose_um']
SCORE_COLUMNS = {'constant_truth', 'constant_prediction', 'nonfinite', *METRICS}
PHASES = {('validation', 'train'), ('validation', 'replay'),
          ('test', 'evaluate'), ('test', 'replay')}


def summarize_folds(scores_by_fold):
    if set(scores_by_fold) != {0, 1, 2}:
        raise ValueError('three distinct outer folds required')
    pooled = summary(pd.concat(scores_by_fold.values(), ignore_index=True))
    dots = [dict(fold=fold, **summary(scores_by_fold[fold])) for fold in sorted(scores_by_fold)]
    spread = []
    for metric in METRICS:
        values = [row[metric] for row in dots if row[metric] is not None]
        spread.append(dict(metric=metric, n_defined_folds=len(values),
                           mean=mean(values) if values else None,
                           sd=stdev(values) if len(values) > 1 else None,
                           estimand='equal-fold mean and sample SD (ddof=1)'))
    return pooled, dots, spread


def validate_scores(path, inputs, origin, partition):
    scores = pd.read_csv(path / f'{partition}_per_sample.csv',
                         dtype={'profile_id': str, 'drug_parent': str, 'cell_line': str},
                         float_precision='round_trip')
    if set(scores) != set(IDENTITIES) | {'cell_line'} | SCORE_COLUMNS:
        raise ValueError('score file column coverage mismatch')
    rows = inputs.rows(origin['scenario'], origin['fold'], partition)
    obs = inputs.obs.iloc[rows]
    expected = pd.DataFrame(dict(profile_id=obs.index.astype(str),
                                 drug_parent=obs.drug_parent_inchikey.to_numpy(),
                                 dose_um=obs.dose_um.to_numpy(),
                                 cell_line=obs.cell_line.astype(str).to_numpy()))
    if not scores[[*IDENTITIES, 'cell_line']].equals(expected):
        raise ValueError(f'{partition} scores differ from exact frozen partition rows')
    actual = summary(scores)
    reported = read_json(path / f'{partition}_metrics.json')
    if set(actual) != set(reported) or any(
        not math.isclose(value, reported[name], rel_tol=1e-12, abs_tol=1e-15)
        if value is not None and reported[name] is not None else value != reported[name]
        for name, value in actual.items()
    ):
        raise ValueError(f'{partition} score summary does not match sealed per-sample rows')
    predictions = np.load(path / f'{partition}_prediction.npy', mmap_mode='r')
    if predictions.shape != (len(scores), 5000) or predictions.dtype != np.float32:
        raise ValueError(f'{partition} prediction dimensions/precision changed')
    for start in range(0, len(scores), 512):
        if not np.isfinite(predictions[start:start + 512]).all():
            raise ValueError(f'{partition} prediction contains nonfinite outputs')
    return scores


def cold_drug_coverage(inputs, scores_by_fold):
    pooled = pd.concat(scores_by_fold.values(), ignore_index=True)
    if pooled.profile_id.duplicated().any() or set(pooled.profile_id) != set(inputs.obs.index):
        raise ValueError('cold-drug test must cover every treated sample once')
    held_parents = set()
    records = []
    for fold in range(3):
        rows = inputs.splits.loc[inputs.splits.scenario.eq('cold-drug') & inputs.splits.fold.eq(fold)]
        parts = {label: rows.loc[rows.partition.eq(label)] for label in ('train', 'validation', 'test')}
        parents = {label: set(frame.drug_parent) for label, frame in parts.items()}
        if (set(rows.partition) != set(parts) or held_parents & parents['test']
                or parents['train'] & parents['test'] or parents['validation'] & parents['test']
                or parents['train'] & parents['validation']
                or len(parents['validation']) != math.ceil(.1 * len(parents['train'] | parents['validation']))):
            raise ValueError('cold-drug parent partition or validation coverage changed')
        held_parents.update(parents['test'])
        records.append(dict(fold=fold, n_population_samples=len(inputs.obs),
                            n_held_parents=len(parents['test']),
                            n_test_parents=len(parents['test']),
                            n_test_samples=len(parts['test']),
                            n_train_samples=len(parts['train']),
                            n_validation_samples=len(parts['validation']),
                            test_population_fraction=len(parts['test']) / len(inputs.obs),
                            test_scope='all dose samples for held parents'))
    if held_parents != set(inputs.obs.drug_parent_inchikey):
        raise ValueError('cold-drug test parent inventory changed')
    return records


def resources(root, origin):
    """Verify the four timed phases plus archived attempts for one origin."""
    path = run_path(root, origin)
    base = root / 'resources' / path.relative_to(root / 'runs')
    registration_digest = sha256(registration_path(root, origin))
    successful = []
    records = []
    locked_at = selection_path(root, origin).stat().st_mtime
    for attempt in sorted(base.iterdir()):
        phase = read_json(attempt / 'phase.json')
        result = read_json(attempt / 'result.json')
        phase_pair = (phase['phase'], phase['action'])
        if (phase_pair not in PHASES or phase['registration_sha256'] != registration_digest
                or phase['origin'] != origin or phase['available_ram_kib'] < 12 * 1024 * 1024
                or result['time_sha256'] != sha256(attempt / 'time.txt')):
            raise ValueError('formal phase resource evidence changed')
        if phase['phase'] == 'test' and phase['started_utc_unix_seconds'] < locked_at:
            raise ValueError('formal outer test started before selection lock')
        if result['returncode'] == 0:
            successful.append(phase_pair)
        timing = (attempt / 'time.txt').read_text().splitlines()
        rss = next(int(line.rsplit(': ', 1)[1]) for line in timing if 'Maximum resident set size' in line)
        records.append(dict(origin, phase=phase['phase'], action=phase['action'],
                            attempt=attempt.name, available_ram_kib=phase['available_ram_kib'],
                            started_utc_unix_seconds=phase['started_utc_unix_seconds'],
                            seconds=result['elapsed_seconds'], peak_rss_kib=rss,
                            returncode=result['returncode']))
    if len(successful) != len(PHASES) or set(successful) != PHASES:
        raise ValueError('missing or duplicate successful timed formal phases')
    return records


def record_details(root, inputs, record, tables):
    origin = {key: record[key] for key in ('source', 'scenario', 'fold', 'method', 'seed')}
    path = run_path(root, origin)
    validate_scores(path, inputs, origin, 'validation')
    scores = validate_scores(path, inputs, origin, 'test')
    learning = pd.read_csv(path / 'learning_curve.csv', float_precision='round_trip')
    terminal = record['terminal']
    best = learning['validation_mse'].idxmin()
    if (learning.epoch.tolist() != list(range(1, terminal['epochs'] + 1))
            or not np.isfinite(learning.select_dtypes('number').to_numpy()).all()
            or int(learning.loc[best, 'epoch']) != terminal['best_epoch']
            or not math.isclose(float(learning.loc[best, 'validation_mse']), terminal['validation_mse'],
                                rel_tol=1e-12, abs_tol=1e-15)):
        raise ValueError('formal learning history differs from restored validation best epoch')
    tables['learning'].append(learning.assign(**origin))
    tables['test_per_sample'].append(scores.assign(**origin))
    tables['fit_resources'].extend(resources(root, origin))
    return scores


def source_tables(root, inputs, records, tables):
    if len(records) != len(METHODS) * 3:
        raise ValueError('source acceptance requires all methods in all three folds')
    for fold in range(3):
        origin = next(record for record in records if record['fold'] == fold)
        register_folder = registration_path(root, origin).parent
        if {path.name for path in register_folder.iterdir()} != {f'{method}.json' for method in METHODS}:
            raise ValueError('formal registration inventory mismatch')
    by_method = defaultdict(dict)
    selected = {}
    for record in records:
        scores = record_details(root, inputs, record, tables)
        method, fold = record['method'], record['fold']
        by_method[method][fold] = scores
        if record['selected_method'] == method:
            selected[fold] = scores.assign(selected_method=method)
    source, scenario = records[0]['source'], records[0]['scenario']
    common = dict(source=source, scenario=scenario, seed=records[0]['seed'])
    for method, scores_by_fold in by_method.items():
        label = dict(common, method=method)
        coverage = (cold_drug_coverage(inputs, scores_by_fold) if scenario == 'cold-drug'
                    else audit_coverage(inputs.splits, scores_by_fold))
        tables['coverage'].extend(dict(label, **row) for row in coverage)
        pooled, dots, spread = summarize_folds(scores_by_fold)
        tables['pooled'].append(dict(label, **pooled))
        tables['folds'].extend(dict(label, **row) for row in dots)
        tables['fold_spread'].extend(dict(label, **row) for row in spread)
        if scenario == 'cold-drug':
            pooled_scores = pd.concat(scores_by_fold.values(), ignore_index=True)
            for dose, rows in pooled_scores.groupby('dose_um'):
                tables['per_dose'].append(dict(label, dose_um=dose, **summary(rows)))
    if scenario == 'cold-drug':
        if set(selected) != {0, 1, 2}:
            raise ValueError('missing locked cold-drug selected procedure fold')
        selected_scores = pd.concat(
            [rows.assign(fold=fold) for fold, rows in sorted(selected.items())], ignore_index=True)
        tables['selected_per_sample'].append(selected_scores.assign(**common))
        tables['selected_pooled'].append(dict(common, method='selected-procedure', **summary(selected_scores)))
    elif selected or any(record['selected_method'] is not None for record in records):
        raise ValueError('auxiliary cold-drug-dose fold must not choose a method')


def publish_scenario(root, scenario, tables, records):
    destination = root / 'reports' / scenario
    destination.parent.mkdir(exist_ok=True)
    with TemporaryDirectory(prefix=f'.{scenario}-', dir=destination.parent) as temporary:
        candidate = Path(temporary)
        for name, rows in tables.items():
            frame = pd.concat(rows, ignore_index=True) if isinstance(rows[0], pd.DataFrame) else pd.DataFrame(rows)
            frame.to_csv(candidate / f'{name}.csv', index=False)
        receipt = dict(schema='tsg-dose-round1-axis-local-scenario-acceptance-v1',
                       scenario=scenario, formal_fits=len(records),
                       plan_sha256=sha256(root / 'formal_plan.json'),
                       origins=[{key: value for key, value in record.items() if key != 'terminal'}
                                for record in records],
                       outputs={file.name: sha256(file) for file in sorted(candidate.iterdir())})
        write_json(candidate / 'acceptance.json', receipt)
        if destination.exists():
            current = {file.name: sha256(file) for file in destination.iterdir()}
            expected = {file.name: sha256(file) for file in candidate.iterdir()}
            if current != expected:
                raise ValueError(f'existing {scenario} report differs from accepted evidence')
        else:
            candidate.rename(destination)
    return receipt


def accept(args):
    """Accept one 42-origin scenario or both scenarios and seal 84 origins."""
    records = verify(args, args.scenario)
    root = args.candidate
    scenarios = SCENARIOS if args.scenario is None else (args.scenario,)
    reports = {scenario: defaultdict(list) for scenario in scenarios}
    for source in SOURCE_POPULATIONS:
        inputs = Inputs(args.inputs / source)
        for scenario in scenarios:
            subset = [record for record in records
                      if record['source'] == source and record['scenario'] == scenario]
            source_tables(root, inputs, subset, reports[scenario])
    receipts = {
        scenario: publish_scenario(root, scenario, reports[scenario],
                                   [record for record in records if record['scenario'] == scenario])
        for scenario in scenarios
    }
    if args.scenario is not None:
        return receipts[args.scenario]
    combined = dict(schema='tsg-dose-round1-axis-local-formal-acceptance-v1',
                    formal_fits=len(records), plan_sha256=sha256(root / 'formal_plan.json'),
                    scenario_receipt_sha256={
                        scenario: sha256(root / 'reports' / scenario / 'acceptance.json')
                        for scenario in scenarios
                    })
    path = root / 'formal_acceptance.json'
    if path.exists() and read_json(path) != combined:
        raise ValueError('existing combined report receipt differs from accepted evidence')
    write_json(path, combined)
    return combined
