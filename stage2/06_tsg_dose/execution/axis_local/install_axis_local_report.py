"""Install a visually accepted report while retaining the sealed 126-fit release."""
import argparse
from contextlib import contextmanager
import fcntl
import json
import os
from pathlib import Path
import re
import subprocess
from tempfile import NamedTemporaryFile

from publish_axis_local_report import (
    CANDIDATE_SCHEMA, PNG_NAMES, SCENARIOS, digest, figure_hashes, snapshot, verify_resource_evidence,
    verify_science,
)

STEP = Path(__file__).resolve().parents[2]
AXIS_LOCAL_ARCHIVE = Path('execution/historical/report_axis_local')
GRAY_ARCHIVE = Path('execution/historical/report_baselines_gray')
PRESENTATION_RECEIPT = Path('log/baselines-gray-acceptance.json')
KPGT_ACCEPTANCE = Path('log/275-kpgt-formal/acceptance.json')
KPGT_PRESENTATION = Path('log/275-kpgt-formal/bar-presentation.json')
KPGT_PRESENTATION_FILES = {
    'data/held_parents.csv', 'metrics.png', 'metrics_cold_drug_dose.png', 'script/render.py',
}
CHECKS = ('check_paths.sh', 'check_tracked.sh', 'verify_all.sh')
HEALTH_SCHEMA = 'tsg-dose-round1-axis-local-health-v1'
SUCCESSOR_SCHEMA = 'tsg-dose-round1-baselines-report-v1'
SUCCESSOR_FILES = {*PNG_NAMES, 'script/render.py',
                   *(f'data/round1_{table}.csv' for table in
                     ('axis_segments', 'curves', 'dose', 'pooled', 'runs', 'summary'))}


def read(path):
    return json.loads(path.read_text())


def write_atomic(path, contents):
    with NamedTemporaryFile('w', dir=path.parent, prefix='.report-receipt-', delete=False) as stream:
        temporary = Path(stream.name)
        try:
            json.dump(contents, stream, indent=2, sort_keys=True)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
            temporary.replace(path)
        finally:
            temporary.unlink(missing_ok=True)


def safe_tree(directory):
    if directory.is_symlink() or not directory.is_dir():
        raise ValueError(f'report directory is missing or a symlink: {directory}')
    if any(path.is_symlink() or not (path.is_file() or path.is_dir()) for path in directory.rglob('*')):
        raise ValueError(f'report tree contains an unsafe entry: {directory}')
    return snapshot(directory)


def health(repo, destination):
    """Run the established offline project checks in the materialized checkout."""
    destination.mkdir(parents=True)
    outputs = {}
    for name in CHECKS:
        script = repo / 'script' / name
        output = destination / f'{name}.log'
        with output.open('w') as stream:
            status = subprocess.run(['bash', str(script)], cwd=repo, stdout=stream,
                                    stderr=subprocess.STDOUT, check=False)
        outputs[name] = dict(log_sha256=digest(output), script_sha256=digest(script),
                             returncode=status.returncode)
        status.check_returncode()
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    write_atomic(destination / 'health.json', dict(schema=HEALTH_SCHEMA, repo=str(repo),
                                                    commit=head, checks=outputs))


def verify_health(repo, receipt_path):
    receipt = read(receipt_path)
    if (receipt.get('schema') != HEALTH_SCHEMA
            or set(receipt.get('checks', {})) != set(CHECKS)
            or not re.fullmatch(r'[0-9a-f]{40}', receipt.get('commit', ''))):
        raise ValueError('three completed offline health checks required')
    for name, record in receipt['checks'].items():
        script = repo / 'script' / name
        log = receipt_path.parent / f'{name}.log'
        if (log.is_symlink() or record['returncode'] != 0
                or digest(log) != record['log_sha256'] or digest(script) != record['script_sha256']):
            raise ValueError(f'offline health evidence changed: {name}')
    return digest(receipt_path)


def locations(step):
    return (step / 'figure', step / 'execution/historical/report_126',
            step / 'log/round1-axis-local-installation.json')


@contextmanager
def deploy_lock(step):
    directory = step / 'execution/checkpoints'
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / '.axis-local-report-install.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        yield


def legacy_snapshot(step):
    previous = read(step / 'log/242-round1-report-acceptance.json')
    expected = previous['files']
    if previous.get('status') != 'passed' or len(expected) != 19:
        raise ValueError('missing complete sealed 126-fit report')
    return expected


def verify_paths(step, candidate):
    active, archive, installed = locations(step)
    reserved = (active, archive, installed)
    if any(path.is_symlink() for path in (candidate, *reserved)):
        raise ValueError('symlinks are forbidden at report deployment paths')
    if any(candidate == path or candidate in path.parents or path in candidate.parents for path in reserved):
        raise ValueError('candidate overlaps deployment paths')
    if candidate.parent.stat().st_dev != active.parent.stat().st_dev or archive.parent.stat().st_dev != active.parent.stat().st_dev:
        raise ValueError('candidate, active and historical report must share one filesystem')
    return reserved


def recover_swap(active, archive, candidate, previous, replacement):
    """Restore the old report after an interrupted rename, without deleting bytes."""
    if not archive.exists():
        if safe_tree(active) != previous:
            raise ValueError('active report no longer matches the historical seal')
        return
    if safe_tree(archive) != previous:
        raise ValueError('historical report changed during installation')
    if active.exists():
        if candidate.exists():
            raise ValueError('ambiguous release recovery: both active and candidate exist')
        if safe_tree(active) != replacement:
            raise ValueError('cannot recover unrecognized active bytes')
        active.rename(candidate)
    elif safe_tree(candidate) != replacement:
        raise ValueError('cannot recover unrecognized candidate bytes')
    archive.rename(active)
    if safe_tree(active) != previous or safe_tree(candidate) != replacement:
        raise ValueError('report recovery did not restore both seals')


def install(step, candidate, receipt_path, health_path):
    receipt = read(receipt_path)
    if (receipt.get('schema') != CANDIDATE_SCHEMA or receipt.get('status') != 'passed'
            or receipt.get('formal_fits') != 84 or receipt.get('visual_inspection') is not True
            or receipt.get('csv_only_figures') is not True):
        raise ValueError('verified CSV-only and visually inspected report required')
    active, archive, installed = verify_paths(step, candidate)
    with deploy_lock(step):
        previous = legacy_snapshot(step)
        if archive.exists() or installed.exists():
            raise ValueError('existing release or interrupted installation; verify or recover first')
        if (safe_tree(candidate) != receipt['files'] or figure_hashes(candidate) != receipt['figures']
                or safe_tree(active) != previous):
            raise ValueError('historical or accepted candidate report bytes changed')
        if (verify_science(receipt['matrix']) != receipt['scientific_artifacts']
                or verify_resource_evidence(Path(receipt['matrix']), SCENARIOS) != receipt['resources']):
            raise ValueError('accepted science or resources changed before installation')
        health_digest = verify_health(step.parents[1], health_path)
        evidence = dict(schema='tsg-dose-round1-axis-local-installation-v1',
                        acceptance_sha256=digest(receipt_path), health_sha256=health_digest,
                        previous_report_receipt='log/242-round1-report-acceptance.json',
                        previous_report_archive=str(archive.relative_to(step)),
                        previous_report_files=previous, active_report_files=receipt['files'])
        active.rename(archive)
        try:
            candidate.rename(active)
            if safe_tree(archive) != previous or safe_tree(active) != receipt['files']:
                raise ValueError('report tree seal failed after placement')
            write_atomic(installed, evidence)
        except BaseException:
            if not installed.exists():
                recover_swap(active, archive, candidate, previous, receipt['files'])
            raise
    return evidence


def successor_files(step, active, predecessor_receipt):
    successor = active / 'baselines'
    acceptance_path = step / 'log/round1-baselines-acceptance.json'
    if not any(path.exists() or path.is_symlink() for path in (successor, acceptance_path)):
        return {}
    if acceptance_path.is_symlink() or not acceptance_path.is_file():
        raise ValueError('missing or unsafe baseline successor receipt')
    acceptance = read(acceptance_path)
    historical = {
        path.name: digest(path) for path in
        (predecessor_receipt, step / 'log/242-round1-report-acceptance.json')
    }
    if (acceptance.get('schema') != SUCCESSOR_SCHEMA or acceptance.get('status') != 'passed'
            or acceptance.get('visual_inspection') is not True
            or acceptance.get('csv_only_figures') is not True
            or acceptance.get('predecessor_fits') != 84 or acceptance.get('ecfp4_mlp_fits') != 12
            or acceptance.get('mean_cell_evaluations') != 12
            or acceptance.get('historical_receipts') != historical
            or set(acceptance.get('files', {})) != SUCCESSOR_FILES):
        raise ValueError('invalid baseline successor receipt')
    files = safe_tree(successor)
    if files != acceptance['files'] or figure_hashes(successor) != acceptance.get('figures'):
        raise ValueError('accepted baseline successor report bytes changed')
    return {f'baselines/{name}': value for name, value in files.items()}


def installed_tree_matches(active, files):
    directories = {parent for name in files for parent in Path(name).parents if parent != Path('.')}
    return (safe_tree(active) == files
            and {path.relative_to(active) for path in active.rglob('*') if path.is_dir()} == directories)


def verify_baseline_presentation(step, active):
    receipt_path = step / PRESENTATION_RECEIPT
    if receipt_path.is_symlink():
        raise ValueError('unsafe baseline presentation receipt')
    receipt = read(receipt_path)
    successor_receipt = step / 'log/round1-baselines-acceptance.json'
    successor_data = {name.replace('data/round1_', 'data/'): value
                      for name, value in read(successor_receipt)['files'].items()
                      if name.startswith('data/')}
    expected = {'script/render.py', *successor_data,
                *(name.removeprefix('round1_') for name in PNG_NAMES)}
    files = receipt.get('files', {})
    published = step / KPGT_ACCEPTANCE
    gray = step / GRAY_ARCHIVE
    active_files = files.copy()
    if published.is_symlink():
        raise ValueError('unsafe tenv KPGT presentation receipt')
    if gray.exists() or gray.is_symlink():
        if not installed_tree_matches(gray, files) or not published.is_file():
            raise ValueError('archived baseline presentation bytes changed or KPGT receipt missing')
        active_files = {}
    if published.exists():
        successor = read(published)
        if (successor.get('schema') != 'tsg-dose-tenv-kpgt-report-v1'
                or successor.get('status') != 'passed'
                or successor.get('visual_inspection') is not True or not successor.get('files')):
            raise ValueError('invalid tenv KPGT presentation receipt')
        if gray.exists():
            active_files = successor['files'].copy()
            presentation = step / KPGT_PRESENTATION
            if presentation.is_symlink():
                raise ValueError('unsafe tenv KPGT figure presentation receipt')
            if presentation.exists():
                figures = read(presentation)
                if (figures.get('schema') != 'tsg-dose-tenv-kpgt-bar-presentation-v1'
                        or figures.get('acceptance_sha256') != digest(published)
                        or set(figures.get('files', {})) != KPGT_PRESENTATION_FILES):
                    raise ValueError('invalid tenv KPGT figure presentation receipt')
                active_files.update(figures['files'])
        else:
            active_files.update({f'tenv_kpgt/{name}': value for name, value in successor['files'].items()})
    if (receipt.get('schema') != 'tsg-dose-baselines-gray-v2'
            or receipt.get('csv_only_replay_byte_identical') is not True
            or receipt.get('visual_inspection') is not True
            or receipt.get('predecessor_sha256') != digest(successor_receipt)
            or receipt.get('predecessor_report_archive') != str(AXIS_LOCAL_ARCHIVE)
            or set(files) != expected
            or any(files[name] != value for name, value in successor_data.items())
            or not installed_tree_matches(active, active_files)):
        raise ValueError('accepted baseline presentation bytes changed')


def recover(step, candidate, receipt_path):
    receipt = read(receipt_path)
    active, archive, installed = verify_paths(step, candidate)
    with deploy_lock(step):
        if installed.exists():
            evidence = read(installed)
            sealed = active
            if (step / PRESENTATION_RECEIPT).exists() or (step / AXIS_LOCAL_ARCHIVE).exists():
                verify_baseline_presentation(step, active)
                sealed = step / AXIS_LOCAL_ARCHIVE
            active_files = evidence['active_report_files'] | successor_files(step, sealed, installed)
            if (evidence['acceptance_sha256'] != digest(receipt_path)
                    or not installed_tree_matches(sealed, active_files)
                    or not installed_tree_matches(archive, evidence['previous_report_files'])):
                raise ValueError('installed report evidence changed')
            return
        recover_swap(active, archive, candidate, legacy_snapshot(step), receipt['files'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='action', required=True)
    check = commands.add_parser('health')
    check.add_argument('--repo', type=Path, required=True)
    check.add_argument('--destination', type=Path, required=True)
    for action in ('install', 'recover'):
        command = commands.add_parser(action)
        command.add_argument('--step', type=Path, default=STEP)
        command.add_argument('--candidate', type=Path, required=True)
        command.add_argument('--receipt', type=Path, required=True)
        if action == 'install':
            command.add_argument('--health', type=Path, required=True)
    args = parser.parse_args()
    if args.action == 'health':
        health(args.repo.resolve(), args.destination.resolve())
    elif args.action == 'install':
        print(json.dumps(install(args.step.resolve(), args.candidate.absolute(), args.receipt, args.health)))
    else:
        recover(args.step.resolve(), args.candidate.absolute(), args.receipt)


if __name__ == '__main__':
    main()
