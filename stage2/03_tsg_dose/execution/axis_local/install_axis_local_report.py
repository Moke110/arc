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
    CANDIDATE_SCHEMA, SCENARIOS, digest, figure_hashes, snapshot, verify_resource_evidence,
    verify_science,
)

STEP = Path(__file__).resolve().parents[2]
CHECKS = ('check_paths.sh', 'check_tracked.sh', 'verify_all.sh')
HEALTH_SCHEMA = 'tsg-dose-round1-axis-local-health-v1'


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


def recover(step, candidate, receipt_path):
    receipt = read(receipt_path)
    active, archive, installed = verify_paths(step, candidate)
    with deploy_lock(step):
        if installed.exists():
            evidence = read(installed)
            if (evidence['acceptance_sha256'] != digest(receipt_path)
                    or safe_tree(active) != evidence['active_report_files']
                    or safe_tree(archive) != evidence['previous_report_files']):
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
