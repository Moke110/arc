"""Verify the Stage 03 baseline matrix and prepare its CSV-only report successor."""
import argparse
from bisect import bisect_left, bisect_right
import csv
from datetime import datetime
import re
import runpy
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import numpy as np
import pandas as pd

from baselines import (ECFP4_MLP, METHODS, PARTITIONS, SCHEMA, SEED, STEP,
                       configuration, open_ecfp4, run_path)
from formal import SCENARIOS, read_json, sha256, write_json
from formal_report import IDENTITIES, summarize_folds
from publish_axis_local_report import figure_hashes, snapshot, verify_science
from tsg_round1.inputs import SOURCE_POPULATIONS
from tsg_round1.protocol import summary

SOURCE_RENDERER = STEP / 'figure/script/render_axis_local.py'
RENDERER = STEP / 'figure/baselines/script/render.py'
RELEASE_SCHEMA = 'tsg-dose-round1-baselines-report-v1'
TABLES = dict(folds='round1_runs.csv', fold_spread='round1_summary.csv',
              per_dose='round1_dose.csv', learning='round1_curves.csv',
              pooled='round1_pooled.csv')
ID_FIELDS = [*IDENTITIES, 'cell_line']
ID_TYPES = {'profile_id': str, 'drug_parent': str, 'cell_line': str}


def validate_fit(path, source, scenario, fold, method, matrix_plan):
    expected_partitions = PARTITIONS if method == ECFP4_MLP else ('test',)
    manifest = path / 'manifest.json'
    result = read_json(manifest)
    config = result['config']
    root = path.parents[4]
    inputs = Path(matrix_plan['inputs'])
    source_directory = inputs / source
    fingerprints = root / 'ecfp4' / source
    open_ecfp4(fingerprints, source_manifest=source_directory / 'manifest.json')
    args = SimpleNamespace(inputs=inputs, scenario=scenario, fold=fold, method=method, device='cuda')
    inputs_identity = SimpleNamespace(source=source, directory=source_directory)
    expected = configuration(args, inputs_identity, fingerprints)
    if (config != expected or config['input_manifest_sha256'] != matrix_plan['input_source_sha256'][source]
            or config['source_pair_sha256'] != matrix_plan['input_pair_sha256']):
        raise ValueError(f'baseline configuration differs from accepted comparison: {path}')
    required = {f'{part}_{kind}.{suffix}' for part in expected_partitions
                for kind, suffix in [('prediction', 'npy'), ('per_sample', 'csv'), ('metrics', 'json')]}
    if method == ECFP4_MLP:
        required |= {'checkpoint.pt', 'learning_curve.csv'}
    if (set(result['outputs']) != required or {f.name for f in path.iterdir()} != required | {'manifest.json'}
            or any(sha256(path / name) != digest for name, digest in result['outputs'].items())):
        raise ValueError(f'baseline output digest or coverage mismatch: {path}')
    replay_path = path.parent / f'{method}-replay.json'
    replay = read_json(replay_path)
    if replay != dict(schema=SCHEMA, status='exact fresh-process replay',
                      partitions=list(expected_partitions), manifest_sha256=sha256(manifest)):
        raise ValueError(f'baseline fresh-process replay mismatch: {path}')
    resource_dir = root / 'resources' / path.relative_to(root / 'runs')
    timing = {phase: resource(resource_dir / f'{phase}.time') for phase in ('fit', 'replay')}
    return result, dict(manifest_sha256=sha256(manifest), replay_sha256=sha256(replay_path),
                        resources=timing)


def resource(path):
    lines = path.read_text().splitlines()
    fields = dict(line.strip().split(': ', 1) for line in lines if ': ' in line)
    if int(fields['Exit status']) != 0:
        raise ValueError(f'failed timed baseline command: {path}')
    wall_time = fields['Elapsed (wall clock) time (h:mm:ss or m:ss)']
    end = path.stat().st_mtime
    seconds = sum(float(part) * 60 ** position
                  for position, part in enumerate(reversed(wall_time.split(':'))))
    return dict(sha256=sha256(path), peak_rss_kib=int(fields['Maximum resident set size (kbytes)']),
                wall_time=wall_time, start_utc_unix_seconds=end - seconds, end_utc_unix_seconds=end)


def monitor_evidence(log, timings):
    text = log.read_text()
    samples = []
    for block in re.finditer(r'^(\S+T\S+)\n[^\n]+\nMem:\s+([^\n]+)\nSwap:[^\n]+\n'
                             r'(\d+) MiB, (\d+) MiB, (\d+) %$', text, re.MULTILINE):
        stamp, memory, _, gpu_memory, gpu_busy = block.groups()
        samples.append((datetime.fromisoformat(stamp).timestamp(), int(memory.split()[-1]),
                        int(gpu_memory), int(gpu_busy)))
    timestamps = [sample[0] for sample in samples]
    if timestamps != sorted(set(timestamps)):
        raise ValueError('resource monitor must contain ordered, unique samples')
    used = set()
    for timing in timings:
        first = bisect_right(timestamps, timing['start_utc_unix_seconds']) - 1
        last = bisect_left(timestamps, timing['end_utc_unix_seconds'])
        if first < 0 or last >= len(samples):
            raise ValueError('resource monitor must cover the complete driver and unit windows')
        if any(second - first > 20
               for first, second in zip(timestamps[first:last], timestamps[first + 1:last + 1])):
            raise ValueError('incomplete 10-second resource evidence')
        used.update(range(first, last + 1))
    window = [samples[index] for index in sorted(used)]
    return dict(sha256=sha256(log), samples=len(window), minimum_available_ram_mib=min(s[1] for s in window),
                peak_gpu_used_mib=max(s[2] for s in window), peak_gpu_utilization_percent=max(s[3] for s in window),
                first_timestamp=int(window[0][0]), last_timestamp=int(window[-1][0]))


def accepted_science(matrix):
    evidence = read_json(STEP / 'log/round1-axis-local-installation.json')
    acceptance = matrix / 'report-candidate-receipt.json'
    if sha256(acceptance) != evidence['acceptance_sha256']:
        raise ValueError('predecessor does not match the installed 84-fit receipt')
    receipt = read_json(acceptance)
    scientific = verify_science(matrix)
    if scientific != receipt['scientific_artifacts']:
        raise ValueError('predecessor scientific artifacts changed')
    return scientific

def legacy_files():
    evidence = read_json(STEP / 'log/round1-axis-local-installation.json')
    inventories = ((STEP / 'figure', evidence['active_report_files']),
                   (STEP / evidence['previous_report_archive'], evidence['previous_report_files']))
    for directory, files in inventories:
        observed = snapshot(directory)
        if directory == STEP / 'figure':
            observed = {name: digest for name, digest in observed.items()
                        if not name.startswith('baselines/')}
        if observed != files:
            raise ValueError(f'sealed scientific report changed: {directory}')
    return {field: sha256(STEP / 'log' / field) for field in
            ('round1-axis-local-installation.json', '242-round1-report-acceptance.json')}


def calculate(matrix, baselines):
    tables = {key: [] for key in TABLES}
    provenance = []
    matrix_plan = read_json(matrix / 'formal_plan.json')
    expected_runs = {run_path(baselines, source, scenario, fold, method)
                     for source in SOURCE_POPULATIONS for scenario in SCENARIOS
                     for fold in range(3) for method in METHODS}
    if set((baselines / 'runs').glob('*/*/*/*/manifest.json')) != {
            path / 'manifest.json' for path in expected_runs}:
        raise ValueError('baseline run inventory differs from the 24 formal coordinates')
    for scenario in SCENARIOS:
        historical = pd.read_csv(matrix / 'reports' / scenario / 'test_per_sample.csv',
                                 usecols=[*ID_FIELDS, 'source', 'fold', 'method'],
                                 dtype=ID_TYPES, float_precision='round_trip')
        for source in SOURCE_POPULATIONS:
            for method in METHODS:
                label = dict(source=source, scenario=scenario, method=method, seed=SEED)
                scores_by_fold = {}
                for fold in range(3):
                    path = run_path(baselines, source, scenario, fold, method)
                    result, evidence = validate_fit(path, source, scenario, fold, method, matrix_plan)
                    provenance.append(dict(label, fold=fold, fit=result['fit'], **evidence))
                    scores = pd.read_csv(path / 'test_per_sample.csv', dtype=ID_TYPES,
                                         float_precision='round_trip')
                    reference = historical.loc[historical.source.eq(source) & historical.fold.eq(fold)
                                               & historical.method.eq('no-dose'), ID_FIELDS]
                    if (not scores[ID_FIELDS].equals(reference.reset_index(drop=True))
                            or summary(scores) != read_json(path / 'test_metrics.json')
                            or scores.nonfinite.any()):
                        raise ValueError(f'unmatched or invalid scored baseline rows: {path}')
                    scores_by_fold[fold] = scores
                    if method == ECFP4_MLP:
                        learning = pd.read_csv(path / 'learning_curve.csv', float_precision='round_trip')
                        if (learning.epoch.tolist() != list(range(1, result['fit']['epochs'] + 1))
                                or not np.isfinite(learning.to_numpy()).all()
                                or int(learning.loc[learning.validation_mse.idxmin(), 'epoch']) != result['fit']['best_epoch']
                                or float(learning.validation_mse.min()) != result['fit']['validation_mse']):
                            raise ValueError(f'invalid best validation epoch: {path}')
                        best = float('inf')
                        best_epoch = 0
                        for row in learning.to_dict('records'):
                            if row['validation_mse'] < best:
                                best = row['validation_mse']
                                best_epoch = row['epoch']
                            tables['learning'].append(dict(row, **label, fold=fold, best_mse=best,
                                                           stale=row['epoch'] - best_epoch, lr=1e-4))
                pooled, folds, spreads = summarize_folds(scores_by_fold)
                tables['pooled'].append(dict(label, **pooled))
                tables['folds'].extend(dict(label, **row) for row in folds)
                tables['fold_spread'].extend(dict(label, **row) for row in spreads)
                if scenario == 'cold-drug':
                    for dose, scores in pd.concat(scores_by_fold.values()).groupby('dose_um'):
                        tables['per_dose'].append(dict(label, dose_um=float(dose), **summary(scores)))
    return tables, provenance


def export_tables(matrix, destination, additions):
    destination.mkdir()
    for scientific, filename in TABLES.items():
        scenarios = ('cold-drug',) if scientific == 'per_dose' else SCENARIOS
        paths = [matrix / 'reports' / scenario / f'{scientific}.csv' for scenario in scenarios]
        with paths[0].open(newline='') as stream:
            fields = next(csv.reader(stream))
        with (destination / filename).open('x', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for path in paths:
                with path.open(newline='') as source:
                    writer.writerows(csv.DictReader(source))
            writer.writerows(additions[scientific])


def reproduce(candidate):
    with TemporaryDirectory(prefix='arc-baseline-csv-only-') as temporary:
        isolated = Path(temporary)
        shutil.copy2(candidate / 'script/render.py', isolated / 'render.py')
        shutil.copy2(SOURCE_RENDERER, isolated / 'reference.py')
        shutil.copytree(candidate / 'data', isolated / 'data')
        subprocess.run([sys.executable, '-I', str(isolated / 'render.py'),
                        '--reference', str(isolated / 'reference.py'), '--data', str(isolated / 'data'),
                        '--output', str(isolated / 'output')], check=True)
        return figure_hashes(isolated / 'output')


def prepare(matrix, baselines, destination, timing, monitor):
    scientific = accepted_science(matrix)
    history = legacy_files()
    driver = resource(timing)
    tables, origins = calculate(matrix, baselines)
    timings = [driver, *(timing for origin in origins for timing in origin['resources'].values())]
    monitored = monitor_evidence(monitor, timings)
    if destination.exists():
        raise ValueError('report candidate must be fresh')
    destination.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix='.baseline-report-', dir=destination.parent) as temporary:
        candidate = Path(temporary) / 'report'
        (candidate / 'script').mkdir(parents=True)
        shutil.copy2(RENDERER, candidate / 'script/render.py')
        export_tables(matrix, candidate / 'data', tables)
        render = runpy.run_path(str(RENDERER))['render']
        render(SOURCE_RENDERER, candidate / 'data', candidate, prepare_axes=True)
        figures = figure_hashes(candidate)
        if reproduce(candidate) != figures or accepted_science(matrix) != scientific or legacy_files() != history:
            raise ValueError('scientific evidence changed or CSV-only replay differs')
        receipt = dict(schema=RELEASE_SCHEMA, status='candidate-verified', csv_only_figures=True,
                       visual_inspection=False, predecessor_fits=84, ecfp4_mlp_fits=12,
                       mean_cell_evaluations=12, scientific_artifacts=scientific,
                       historical_receipts=history, files=snapshot(candidate), figures=figures,
                       origins=origins, driver=driver, resource_monitor=monitored,
                       baseline_runtime_sha256=sha256(Path(__file__).with_name('baselines.py')))
        candidate.rename(destination)
    return receipt


def accept(candidate, receipt):
    if (receipt.get('schema') != RELEASE_SCHEMA or receipt.get('status') != 'candidate-verified'
            or snapshot(candidate) != receipt['files'] or reproduce(candidate) != receipt['figures']
            or legacy_files() != receipt['historical_receipts']):
        raise ValueError('candidate or legacy report changed before visual acceptance')
    return dict(receipt, status='passed', visual_inspection=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='action', required=True)
    for name in ('prepare', 'accept'):
        command = commands.add_parser(name)
        command.add_argument('--candidate', type=Path, required=True)
        command.add_argument('--receipt', type=Path, required=True)
    for name in ('matrix', 'baselines', 'timing', 'monitor'):
        commands.choices['prepare'].add_argument(f'--{name}', type=Path, required=True)
    commands.choices['accept'].add_argument('--visual-inspection', action='store_true', required=True)
    args = parser.parse_args()
    if args.candidate.resolve() in args.receipt.resolve().parents:
        parser.error('write the receipt outside the report candidate')
    if args.action == 'prepare':
        result = prepare(args.matrix, args.baselines, args.candidate, args.timing, args.monitor)
    else:
        result = accept(args.candidate, read_json(args.receipt))
    write_json(args.receipt, result)
