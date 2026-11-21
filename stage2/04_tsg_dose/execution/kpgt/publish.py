"""Verify scientific evidence and prepare or accept the tenv KPGT report."""

import argparse
from collections import defaultdict
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import numpy as np
import pandas as pd

from formal import (
    STEP,
    SCENARIOS,
    candidate_lock,
    planned_runs,
    read_json,
    sha256,
    verify,
    write_json,
)
from baseline import (
    MLP,
    PHASES as BASELINE_PHASES,
    configuration,
    paths as baseline_paths,
    verify_fit,
    verify_replay,
)
from render import FIGURES
from tsg_kpgt.inputs import Inputs
from tsg_round1.protocol import summary

from auxiliary import audit_coverage
from formal_report import (
    PHASES as TRANSIGEN_PHASES,
    cold_drug_coverage,
    source_tables,
    summarize_folds,
    validate_scores,
)
from install_axis_local_report import (
    AXIS_LOCAL_ARCHIVE,
    installed_tree_matches,
    safe_tree,
    successor_files,
    verify_baseline_presentation,
)
from publish_baselines import verify_mean_cell
from release_resources import samples, verify_window

SCHEMA = 'tsg-dose-tenv-kpgt-report-v1'
RENDERER = Path(__file__).with_name('render.py')


def predecessor_inventory(step):
    verify_baseline_presentation(step, step / 'figure')
    installation = step / 'log/round1-axis-local-installation.json'
    evidence = read_json(installation)
    sealed = step / AXIS_LOCAL_ARCHIVE
    files = evidence['active_report_files'] | successor_files(
        step, sealed, installation
    )
    if not installed_tree_matches(sealed, files) or not installed_tree_matches(
        step / evidence['previous_report_archive'], evidence['previous_report_files']
    ):
        raise ValueError('sealed predecessor artifacts changed')
    names = (
        'round1-axis-local-installation.json',
        'round1-baselines-acceptance.json',
        '242-round1-report-acceptance.json',
        'baselines-gray-acceptance.json',
    )
    return {name: sha256(step / 'log' / name) for name in names}


def baseline_resources(directory, origin):
    records = []
    for attempt in sorted(directory.iterdir()):
        phase = read_json(attempt / 'phase.json')
        result = read_json(attempt / 'result.json')
        timing = attempt / 'time.txt'
        if (
            phase['origin'] != origin
            or phase['phase'] not in BASELINE_PHASES
            or phase['available_ram_kib'] < 12 * 1024 * 1024
            or result['time_sha256'] != sha256(timing)
        ):
            raise ValueError('baseline fit/replay timing changed')
        rss = re.search(
            r'Maximum resident set size \(kbytes\):\s*(\d+)', timing.read_text()
        )
        if rss is None:
            raise ValueError('baseline fit/replay has no peak RSS')
        records.append(
            dict(
                origin,
                phase=phase['phase'],
                action=phase['phase'],
                attempt=attempt.name,
                available_ram_kib=phase['available_ram_kib'],
                started_utc_unix_seconds=phase['started_utc_unix_seconds'],
                seconds=result['elapsed_seconds'],
                peak_rss_kib=int(rss[1]),
                returncode=result['returncode'],
            )
        )
    return records


def baseline_scores(args, inputs, scenario, tables, coverage):
    options = SimpleNamespace(inputs=args.inputs, device=args.device)
    planned = [
        origin
        for origin in planned_runs()['baselines']
        if origin['scenario'] == scenario
    ]
    by_method = defaultdict(dict)
    evidence = []
    for origin in planned:
        method, fold = origin['method'], origin['fold']
        path, resources = baseline_paths(args.baselines, origin)
        fit = verify_fit(path, configuration(options, origin))
        replay = verify_replay(path, origin, fit)
        scores = validate_scores(path, inputs, origin, 'test')
        if scores.nonfinite.any():
            raise ValueError('baseline test metrics include nonfinite sample flags')
        result = dict(
            origin=origin,
            manifest_sha256=replay['manifest_sha256'],
            replay_sha256=sha256(path.parent / f'{method}-replay.json'),
            fit=fit['fit'],
        )
        if method == MLP:
            validate_scores(path, inputs, origin, 'validation')
            curve = pd.read_csv(
                path / 'learning_curve.csv', float_precision='round_trip'
            )
            key = curve.validation_mse.idxmin()
            if (
                not np.isfinite(curve.to_numpy()).all()
                or curve.epoch.tolist() != list(range(1, fit['fit']['epochs'] + 1))
                or fit['fit']['weights_updated'] is not True
                or int(curve.loc[key, 'epoch']) != fit['fit']['best_epoch']
                or float(curve.loc[key, 'validation_mse'])
                != fit['fit']['validation_mse']
            ):
                raise ValueError(
                    'baseline learning curve does not identify the best restored epoch'
                )
            tables['learning'].append(curve.assign(**origin))
        else:
            result['independent_replay'] = verify_mean_cell(
                inputs, scenario, fold, path
            )
        evidence.append(result)
        tables['fit_resources'].extend(baseline_resources(resources, origin))
        by_method[method][fold] = scores
    for method, by_fold in by_method.items():
        actual = (
            cold_drug_coverage(inputs, by_fold)
            if scenario == 'cold-drug'
            else audit_coverage(inputs.splits, by_fold)
        )
        if actual != coverage:
            raise ValueError(
                'baseline test population differs from the TranSiGen test population'
            )
        label = dict(source='tenv', scenario=scenario, method=method, seed=42)
        pooled, folds, spread = summarize_folds(by_fold)
        tables['pooled'].append(dict(label, **pooled))
        tables['folds'].extend(dict(label, **row) for row in folds)
        tables['fold_spread'].extend(dict(label, **row) for row in spread)
        if scenario == 'cold-drug':
            for dose, scores in pd.concat(by_fold.values()).groupby('dose_um'):
                tables['per_dose'].append(dict(label, dose_um=dose, **summary(scores)))
    return evidence


def scientific_tables(args):
    inputs = Inputs(args.inputs)
    records = verify(
        SimpleNamespace(
            inputs=args.inputs,
            source=args.source,
            candidate=args.formal,
            device=args.device,
        ),
        None,
    )
    actual = {path.parent for path in args.baselines.glob('runs/*/*/*/manifest.json')}
    expected = {
        baseline_paths(args.baselines, origin)[0]
        for origin in planned_runs()['baselines']
    }
    if actual != expected:
        raise ValueError('baseline fit inventory differs from the 12 tenv coordinates')
    combined = defaultdict(list)
    baseline_evidence = []
    for scenario in SCENARIOS:
        subset = [
            dict(record, **record['origin'])
            for record in records
            if record['origin']['scenario'] == scenario
        ]
        tables = defaultdict(list)
        source_tables(args.formal, inputs, subset, tables)
        coverage = [
            {
                key: value
                for key, value in row.items()
                if key not in ('source', 'scenario', 'seed', 'method')
            }
            for row in tables['coverage']
            if row['method'] == 'no-dose'
        ]
        baseline_evidence.extend(
            baseline_scores(args, inputs, scenario, tables, coverage)
        )
        tables['coverage'] = [dict(row, scenario=scenario) for row in coverage]
        # Raw per-sample score bytes remain bound by each fit manifest.
        del tables['test_per_sample']
        if scenario == 'cold-drug':
            selected = pd.concat(tables.pop('selected_per_sample'), ignore_index=True)
            tables['selection'] = [
                dict(
                    scenario=scenario,
                    fold=int(fold),
                    method=rows.selected_method.iat[0],
                    **summary(rows),
                )
                for fold, rows in selected.groupby('fold')
            ]
        for name, rows in tables.items():
            combined[name].extend(rows)
    origins = [
        {key: item for key, item in record.items() if key != 'terminal'}
        for record in records
    ]
    return combined, origins, baseline_evidence


def check_resource_window(records, family, cap, timing, monitor, driver_log):
    """Bind scope limits and consecutive samples to every fit phase."""
    limits = dict(
        line.split('=', 1) for line in cap.read_text().splitlines() if '=' in line
    )
    invocation = limits.get('InvocationID', '')
    launch = f'Running as unit: {limits.get("Id")}; invocation ID: {invocation}'
    log = driver_log.read_text()
    if (
        not re.fullmatch(r'[0-9a-f]{32}', invocation)
        or not limits.get('ActiveEnterTimestamp')
        or limits.get('MemoryMax') != str(10 * 1024**3)
        or limits.get('MemorySwapMax') != str(8 * 1024**3)
        or launch not in log
    ):
        raise ValueError(f'{family} launch does not bind the required 10/8 GiB scope')
    driver_time = timing.read_text()
    rss = re.search(r'Maximum resident set size \(kbytes\):\s*([1-9]\d*)', driver_time)
    if rss is None or not re.search(
        r'^\s*Exit status:\s*0\s*$', driver_time, re.MULTILINE
    ):
        raise ValueError(f'{family} driver has no successful wall-time/RSS receipt')
    success = [row for row in records if row['returncode'] == 0]
    phase_pairs = (
        TRANSIGEN_PHASES
        if family == 'transigen'
        else {(phase, phase) for phase in BASELINE_PHASES}
    )
    coordinates = ('source', 'scenario', 'fold', 'method', 'seed')
    expected = {
        (*(origin[key] for key in coordinates), *phase)
        for origin in planned_runs()[family]
        for phase in phase_pairs
    }
    observed = {
        (*(row[key] for key in coordinates), row['phase'], row['action'])
        for row in success
    }
    if observed != expected or len(success) != len(expected):
        raise ValueError(f'incomplete or duplicate {family} resource receipts')
    bounds = []
    for row in success:
        begin, elapsed = row['started_utc_unix_seconds'], row['seconds']
        if not math.isfinite(begin) or not math.isfinite(elapsed) or elapsed <= 0:
            raise ValueError(f'invalid {family} phase time bounds')
        bounds.append((begin, begin + elapsed))
    bounds.sort()
    if any(end > begin for (_, end), (begin, _) in zip(bounds, bounds[1:])):
        raise ValueError(f'{family} fits exceeded serial concurrency')
    begin, end = bounds[0][0], bounds[-1][1]
    points = verify_window(samples(monitor, gpu=True), begin, end, family)
    return dict(
        scope=limits['Id'],
        invocation_id=invocation,
        memory_max=int(limits['MemoryMax']),
        memory_swap_max=int(limits['MemorySwapMax']),
        first_phase=begin,
        last_phase=end,
        samples=len(points),
        min_available_ram_kib=min(point['MemAvailable'] for point in points),
        max_gpu_used_mib=max(point['gpu'][0] for point in points),
        max_driver_rss_kib=int(rss[1]),
        successful_phases=len(success),
        cap_sha256=sha256(cap),
        driver_log_sha256=sha256(driver_log),
        timing_sha256=sha256(timing),
        telemetry_sha256=sha256(monitor),
    )


def report_resources(args, records):
    resources = {}
    for family, stem in (('transigen', 'formal'), ('baselines', 'baseline')):
        work = args.resources
        relevant = [
            row
            for row in records
            if (row['method'] in ('KPGT_MLP', 'Mean_cell')) == (family == 'baselines')
        ]
        resources[family] = check_resource_window(
            relevant,
            family,
            work / f'{stem}-cap-actual.txt',
            work / f'{stem}-driver.time',
            args.monitor,
            work / f'{stem}-driver.log',
        )
    if resources['transigen']['last_phase'] > resources['baselines']['first_phase']:
        raise ValueError('TranSiGen and baseline execution overlapped')
    return resources


def render_subprocess(script, data, output):
    subprocess.run(
        [
            sys.executable,
            '-I',
            str(script),
            '--data',
            str(data),
            '--output',
            str(output),
        ],
        check=True,
    )
    return {name: sha256(output / name) for name in FIGURES}


def reproduce(candidate):
    with TemporaryDirectory(prefix='arc-tenv-kpgt-csv-') as temporary:
        isolated = Path(temporary)
        shutil.copyfile(candidate / 'script/render.py', isolated / 'render.py')
        shutil.copytree(candidate / 'data', isolated / 'data')
        return render_subprocess(
            isolated / 'render.py', isolated / 'data', isolated / 'figures'
        )


def prepare(args):
    if args.candidate.exists() or args.receipt.exists():
        raise FileExistsError('publication requires a fresh candidate and receipt')
    history = predecessor_inventory(STEP)
    with candidate_lock(args.formal), candidate_lock(args.baselines):
        tables, origins, baselines = scientific_tables(args)
        resources = report_resources(args, tables['fit_resources'])
    args.candidate.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(
        prefix='.tenv-kpgt-', dir=args.candidate.parent
    ) as temporary:
        candidate = Path(temporary)
        (candidate / 'data').mkdir()
        (candidate / 'script').mkdir()
        for name, rows in sorted(tables.items()):
            frame = (
                pd.concat(rows, ignore_index=True)
                if isinstance(rows[0], pd.DataFrame)
                else pd.DataFrame(rows)
            )
            frame.to_csv(candidate / 'data' / f'{name}.csv', index=False)
        shutil.copyfile(RENDERER, candidate / 'script/render.py')
        figures = render_subprocess(
            candidate / 'script/render.py', candidate / 'data', candidate
        )
        if reproduce(candidate) != figures or predecessor_inventory(STEP) != history:
            raise ValueError('CSV-only render replay or predecessor inventory changed')
        evidence = dict(
            schema=SCHEMA,
            status='candidate-verified',
            source='tenv',
            representation='KPGT',
            formal_fits=len(origins),
            kpgt_mlp_fits=sum(row['origin']['method'] == MLP for row in baselines),
            mean_cell_evaluations=sum(
                row['origin']['method'] != MLP for row in baselines
            ),
            input_manifest_sha256=sha256(args.inputs / 'manifest.json'),
            formal_plan_sha256=sha256(args.formal / 'formal_plan.json'),
            origins=origins,
            baselines=baselines,
            resources=resources,
            historical_receipts=history,
            csv_only_replay=True,
            figures=figures,
            files=safe_tree(candidate),
        )
        candidate.rename(args.candidate)
    write_json(args.receipt, evidence)


def accept(candidate, receipt_path, acceptance_path):
    receipt = read_json(receipt_path)
    if (
        receipt['schema'] != SCHEMA
        or receipt['status'] != 'candidate-verified'
        or safe_tree(candidate) != receipt['files']
        or reproduce(candidate) != receipt['figures']
        or predecessor_inventory(STEP) != receipt['historical_receipts']
    ):
        raise ValueError('scientific or CSV-only publication evidence changed')
    result = dict(
        receipt,
        status='passed',
        visual_inspection=True,
        candidate_receipt_sha256=sha256(receipt_path),
    )
    if acceptance_path.exists():
        raise FileExistsError('acceptance path already exists')
    write_json(acceptance_path, result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_subparsers(dest='action', required=True)
    for action in ('prepare', 'accept'):
        command = actions.add_parser(action)
        for name in ('candidate', 'receipt'):
            command.add_argument(f'--{name}', type=Path, required=True)
    command = actions.choices['prepare']
    for name in ('inputs', 'source', 'formal', 'baselines', 'resources', 'monitor'):
        command.add_argument(f'--{name}', type=Path, required=True)
    command.add_argument('--device', default='cuda')
    command = actions.choices['accept']
    command.add_argument('--acceptance', type=Path, required=True)
    command.add_argument('--visual-inspection', action='store_true', required=True)
    args = parser.parse_args()
    if args.action == 'prepare':
        for name in ('inputs', 'source', 'formal', 'baselines'):
            setattr(args, name, getattr(args, name).resolve())
        args.candidate = args.candidate.resolve()
        prepare(args)
    else:
        accept(args.candidate, args.receipt, args.acceptance)


if __name__ == '__main__':
    main()
