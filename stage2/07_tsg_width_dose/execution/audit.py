"""Verify the Step 05 held profiles and compare original-width fits to Step 04."""

import argparse
from collections import defaultdict
import csv
from decimal import Decimal
from pathlib import Path

from run import ACCEPTED, METHODS, SCENARIO, STEP04, digest, load_json, run_path, save_json, verify


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_csv(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))


def population(rows):
    result = {row['profile_id']: (row['drug_parent'], Decimal(row['dose_um'])) for row in rows}
    require(len(result) == len(rows), 'Duplicate profile IDs')
    return result


def frozen_test_population(inputs):
    manifest = inputs / 'tenv/manifest.json'
    splits = inputs / 'tenv/splits.csv'
    require(digest(manifest) == load_json(inputs / 'manifest.json')['files']['tenv/manifest.json'],
            'tenv input manifest digest mismatch')
    require(digest(splits) == load_json(manifest)['outputs']['splits.csv'], 'Split digest mismatch')
    by_fold = defaultdict(list)
    for row in read_csv(splits):
        if row['scenario'] == SCENARIO and row['partition'] == 'test':
            by_fold[int(row['fold'])].append(row)
    expected = {fold: population(rows) for fold, rows in by_fold.items()}
    doses = {fold: {dose for _, dose in rows.values()} for fold, rows in expected.items()}
    require(doses == {0: {Decimal('0.2')}, 1: {Decimal('1')}, 2: {Decimal('10')}},
            'Frozen test concentrations differ from the registered folds')
    return expected


def compare(run, prior, source_manifest, expected, origin, registered, terminal, metrics):
    require(digest(prior / 'manifest.json') == source_manifest, f'Step 04 run manifest mismatch: {prior}')
    before = load_json(prior / 'manifest.json')['outputs']
    for name in ('terminal.json', 'test_metrics.json', 'test_per_sample.csv',
                 'learning_curve.csv', 'validation_prediction.npy'):
        require(digest(prior / name) == before[name], f'Step 04 artifact digest mismatch: {prior / name}')
    require(population(read_csv(prior / 'test_per_sample.csv')) == expected,
            f'Step 04 held profiles differ from the frozen split: {prior}')
    prior_terminal = load_json(prior / 'terminal.json')
    prior_metrics = load_json(prior / 'test_metrics.json')
    require(prior_metrics.keys() == metrics.keys(), f'Step 04 metric keys differ: {prior}')
    differences = {key: {'step04': prior_metrics[key], 'step05': value}
                   for key, value in metrics.items() if value != prior_metrics[key]}
    for key in ('epochs', 'best_epoch', 'validation_mse', 'parameter_count'):
        if terminal[key] != prior_terminal[key]:
            differences[key] = {'step04': prior_terminal[key], 'step05': terminal[key]}
    return dict(
        fold=origin['fold'], method=origin['method'], step04_manifest_sha256=source_manifest,
        step05_manifest_sha256=digest(run / 'manifest.json'), held_rows=len(expected),
        identical_per_sample_csv=before['test_per_sample.csv'] == registered['test_per_sample.csv'],
        identical_validation_prediction=before['validation_prediction.npy'] == registered['validation_prediction.npy'],
        identical_learning_curve=before['learning_curve.csv'] == registered['learning_curve.csv'],
        differences=differences,
    )


def audit(candidate, freeze):
    origins = verify(candidate)
    input_digest = digest(freeze / 'inputs/manifest.json')
    source_receipt = STEP04 / 'log/275-kpgt-formal/acceptance.json'
    release = load_json(source_receipt)
    require(release['status'] == 'passed', 'Step 04 release did not pass')
    require(release['input_manifest_sha256'] == input_digest == digest(ACCEPTED),
            'Step 04 input digest does not match the accepted freeze')
    expected = frozen_test_population(freeze / 'inputs')
    source = {(row['origin']['method'], row['origin']['fold']): row for row in release['origins']
              if row['origin']['scenario'] == SCENARIO}
    comparisons = []
    for origin in sorted(origins, key=lambda value: run_path(candidate, value)):
        run = run_path(candidate, origin)
        terminal = load_json(run / 'terminal.json')
        require(terminal['status'] == 'complete' and terminal['validation_replay'] == 'exact'
                and terminal['evaluation_pid'] != terminal['training_pid'],
                f'Missing independent exact validation replay: {run}')
        samples = read_csv(run / 'test_per_sample.csv')
        metrics = load_json(run / 'test_metrics.json')
        fold = origin['fold']
        require(len(samples) == terminal['test_samples'] == metrics['n_samples']
                and population(samples) == expected[fold], f'Held profiles differ from the frozen split: {run}')
        if origin['width'] == 'original':
            method = METHODS[origin['method']]
            prior = freeze / 'retry/formal/runs/tenv' / SCENARIO / f'fold-{fold}' / method
            manifest = source[method, fold]['run_manifest_sha256']
            registered = load_json(run / 'manifest.json')['files']
            comparisons.append(compare(run, prior, manifest, expected[fold], origin,
                                       registered, terminal, metrics))
    return dict(verified_coordinates=len(origins), held_rows={fold: len(rows) for fold, rows in expected.items()},
                step04_receipt_sha256=digest(source_receipt), original_comparison=comparisons)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--step04-freeze', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    save_json(args.output, audit(args.candidate, args.step04_freeze))


if __name__ == '__main__':
    main()
