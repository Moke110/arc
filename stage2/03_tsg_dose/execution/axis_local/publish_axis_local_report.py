"""Prepare and accept the dual-source figure candidate after formal acceptance."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

STEP = Path(__file__).resolve().parents[2]
RENDERER = Path(__file__).with_name('render_axis_local.py').resolve()
sys.path.insert(0, str(STEP / 'execution'))
sys.path.insert(0, str(RENDERER.parent))
from publish_report import digest, snapshot
from render_axis_local import METHODS, PNG_NAMES, SCENARIOS, SOURCES
from release_resources import OPTIONAL_RESOURCE_FILES, RESOURCE_FILES, verify_resource_evidence

COMBINED_SCHEMA = 'tsg-dose-round1-axis-local-formal-acceptance-v1'
SCENARIO_SCHEMA = 'tsg-dose-round1-axis-local-scenario-acceptance-v1'
CANDIDATE_SCHEMA = 'tsg-dose-round1-axis-local-figure-acceptance-v1'
TABLES = frozenset(('coverage.csv', 'fit_resources.csv', 'folds.csv', 'fold_spread.csv',
                    'learning.csv', 'pooled.csv', 'test_per_sample.csv'))
COLD_DRUG_TABLES = frozenset(('per_dose.csv', 'selected_per_sample.csv', 'selected_pooled.csv'))


def verify_science(matrix):
    """Verify sealed scalar plot inputs and source-local origin identities."""
    matrix = Path(matrix)
    combined_path = matrix / 'formal_acceptance.json'
    if combined_path.is_symlink() or not combined_path.is_file():
        raise ValueError('complete 84-origin scientific acceptance required')
    combined = json.loads(combined_path.read_text())
    plan = matrix / 'formal_plan.json'
    if (combined.get('schema') != COMBINED_SCHEMA or combined.get('formal_fits') != 84
            or set(combined.get('scenario_receipt_sha256', {})) != set(SCENARIOS)
            or plan.is_symlink() or not plan.is_file()
            or combined.get('plan_sha256') != digest(plan)):
        raise ValueError('complete 84-origin scientific acceptance required')
    artifacts = {
        'formal_acceptance.json': digest(combined_path),
        'formal_plan.json': combined['plan_sha256'],
    }

    def require(path, sha256):
        relatives = (path, *[parent for parent in path.parents if matrix in parent.parents])
        if any(item.is_symlink() for item in relatives) or not path.is_file() or digest(path) != sha256:
            raise ValueError(f'accepted scientific artifact changed: {path.relative_to(matrix)}')
        artifacts[str(path.relative_to(matrix))] = sha256

    for scenario in SCENARIOS:
        directory = matrix / 'reports' / scenario
        receipt_path = directory / 'acceptance.json'
        require(receipt_path, combined['scenario_receipt_sha256'][scenario])
        receipt = json.loads(receipt_path.read_text())
        origins = receipt.get('origins', [])
        expected = {(source, scenario, fold, method, 42)
                    for source in SOURCES for fold in range(3) for method in METHODS}
        fields = ('source', 'scenario', 'fold', 'method', 'seed')
        if (receipt.get('schema') != SCENARIO_SCHEMA or receipt.get('scenario') != scenario
                or receipt.get('formal_fits') != 42 or receipt.get('plan_sha256') != combined['plan_sha256']
                or len(origins) != 42 or {tuple(origin[key] for key in fields) for origin in origins} != expected):
            raise ValueError(f'incomplete 42-origin scientific acceptance for {scenario}')
        for origin in origins:
            source, _, fold, method, _ = (origin[key] for key in fields)
            require(matrix / 'runs' / source / scenario / f'fold-{fold}' / method / 'manifest.json',
                    origin['run_manifest_sha256'])
            require(matrix / 'selection' / source / scenario / f'fold-{fold}.json',
                    origin['selection_sha256'])
        tables = TABLES | (COLD_DRUG_TABLES if scenario == 'cold-drug' else set())
        if set(receipt.get('outputs', {})) != tables or {path.name for path in directory.iterdir()} != tables | {'acceptance.json'}:
            raise ValueError(f'scientific output inventory mismatch for {scenario}')
        for table in sorted(tables):
            require(directory / table, receipt['outputs'][table])
    resource_hashes = {name: digest(matrix / name) for name in (*RESOURCE_FILES, *OPTIONAL_RESOURCE_FILES) if (matrix / name).is_file()}
    verify_resource_evidence(matrix, SCENARIOS)
    for name, expected in resource_hashes.items():
        require(matrix / name, expected)
    return artifacts


def run_renderer(script, *args, python=sys.executable, mplconfig):
    subprocess.run([str(python), '-I', str(script), *map(str, args)],
                   check=True, env=dict(os.environ, MPLCONFIGDIR=str(mplconfig)))


def figure_hashes(directory):
    files = {path.name: digest(path) for path in Path(directory).glob('*.png')}
    if set(files) != set(PNG_NAMES):
        raise ValueError('unexpected figure inventory')
    return files


def reproduce(script, data, python):
    """Render with only a script and plotting CSVs in a fresh working tree."""
    with TemporaryDirectory(prefix='arc-axis-local-csv-') as temporary:
        isolated = Path(temporary)
        shutil.copy2(script, isolated / 'render.py')
        shutil.copytree(data, isolated / 'data')
        run_renderer(isolated / 'render.py', '--data', isolated / 'data',
                     '--output', isolated / 'output', python=python, mplconfig=isolated / 'mpl')
        return figure_hashes(isolated / 'output')


def prepare_candidate(matrix, destination, python=sys.executable):
    matrix, destination = Path(matrix).resolve(), Path(destination).resolve()
    scientific_artifacts = verify_science(matrix)
    if destination.exists():
        raise ValueError('report candidate destination must be fresh')
    destination.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix='.axis-local-figure-', dir=destination.parent) as temporary:
        work = Path(temporary)
        candidate = work / 'figure'
        script = candidate / 'script/render_axis_local.py'
        script.parent.mkdir(parents=True)
        shutil.copy2(RENDERER, script)
        data = candidate / 'data'
        run_renderer(script, '--reports', matrix / 'reports', '--data', data,
                     python=python, mplconfig=work / 'mpl')
        run_renderer(script, '--data', data, '--output', candidate,
                     python=python, mplconfig=work / 'mpl')
        figures = figure_hashes(candidate)
        if reproduce(script, data, python) != figures:
            raise ValueError('CSV-only reproduction differs from candidate figures')
        if verify_science(matrix) != scientific_artifacts:
            raise ValueError('scientific evidence changed during plotting')
        receipt = dict(schema=CANDIDATE_SCHEMA, status='candidate-verified', formal_fits=84,
                       csv_only_figures=True, visual_inspection=False, matrix=str(matrix),
                       scientific_artifacts=scientific_artifacts,
                       resources=verify_resource_evidence(matrix, SCENARIOS), files=snapshot(candidate),
                       figures=figures)
        candidate.rename(destination)
    return receipt


def accept_candidate(candidate, receipt, *, visually_inspected):
    if (not visually_inspected or receipt.get('schema') != CANDIDATE_SCHEMA
            or receipt.get('status') != 'candidate-verified' or receipt.get('formal_fits') != 84
            or receipt.get('csv_only_figures') is not True):
        raise ValueError('visual inspection of a verified dual-source candidate required')
    if snapshot(candidate) != receipt['files'] or figure_hashes(candidate) != receipt['figures']:
        raise ValueError('verified candidate changed before visual acceptance')
    if (verify_science(receipt['matrix']) != receipt['scientific_artifacts']
            or verify_resource_evidence(Path(receipt['matrix']), SCENARIOS) != receipt.get('resources')):
        raise ValueError('scientific evidence changed before visual acceptance')
    return dict(receipt, status='passed', visual_inspection=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    for name in ('prepare', 'accept'):
        sub = commands.add_parser(name)
        sub.add_argument('--candidate', type=Path, required=True)
        sub.add_argument('--receipt', type=Path, required=True)
    commands.choices['prepare'].add_argument('--matrix', type=Path, required=True)
    commands.choices['accept'].add_argument('--visual-inspection', action='store_true', required=True)
    args = parser.parse_args()
    if args.candidate.resolve() in args.receipt.resolve().parents:
        parser.error('store the acceptance receipt outside the candidate')
    if args.command == 'prepare':
        result = prepare_candidate(args.matrix, args.candidate)
    else:
        result = accept_candidate(args.candidate, json.loads(args.receipt.read_text()),
                                  visually_inspected=args.visual_inspection)
    args.receipt.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
