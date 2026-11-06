"""Review-gated v2 migration transaction; never starts training.

Dry-run re-stages all identities and inventories. Apply additionally requires the
reviewed inventory SHA, actual host PID namespace and stopped systemd service.
An interrupted transaction resumes forward from its durable journal; it never
credits time, overwrites legacy bytes, or dispatches workers. A scheduler must
refuse a transaction whose journal status is not committed.
"""
import argparse
import ctypes
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

from perturbdiff_v2_boundary import digest, sync_directory
from perturbdiff_v2_migration import read, stage_locked

NAMES = {'queue.json', 'commands.json', 'expected_identities.json', 'budget-ledger.json', 'plan.json'}


def publish(path, value):
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('w') as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)
    sync_directory(path.parent)


def inventory(root):
    if root.is_symlink():
        raise ValueError('symlink refused: ' + str(root))
    if root.is_file():
        return {'.': digest(root)}
    if not root.is_dir():
        raise ValueError('missing inventory root: ' + str(root))
    result = {}
    for path in sorted(root.rglob('*')):
        if path.is_symlink() or not (path.is_dir() or path.is_file()):
            raise ValueError('nonregular entry refused: ' + str(path))
        result[str(path.relative_to(root))] = digest(path) if path.is_file() else None
    return result


def mkdir_durable(path):
    missing = []
    cursor = path
    while not cursor.exists():
        missing.append(cursor)
        cursor = cursor.parent
    for directory in reversed(missing):
        directory.mkdir()
        sync_directory(directory)
        sync_directory(directory.parent)


def sync_tree(root):
    for path in root.rglob('*'):
        if path.is_file():
            with path.open('rb') as stream:
                os.fsync(stream.fileno())
    for path in sorted((p for p in root.rglob('*') if p.is_dir()), reverse=True):
        sync_directory(path)
    sync_directory(root)
    sync_directory(root.parent)


def exchange(left, right):
    # Single atomic namespace operation publishes queue/commands/identities and
    # ledger together. Fail closed if the filesystem/kernel lacks RENAME_EXCHANGE.
    libc = ctypes.CDLL(None, use_errno=True)
    call = libc.renameat2
    call.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
    call.restype = ctypes.c_int
    if call(-100, os.fsencode(left), -100, os.fsencode(right), 2):
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error))
    sync_directory(left.parent)
    if right.parent != left.parent:
        sync_directory(right.parent)


def host_stopped(stopped, namespace, service):
    if (not namespace or os.readlink('/proc/self/ns/pid') != namespace or
            Path('/proc/1/comm').read_text().strip() != 'systemd'):
        raise ValueError('actual host PID namespace required')
    request = read(stopped / 'request.json')
    if Path('/proc/sys/kernel/random/boot_id').read_text().strip() != request['target']['boot_id']:
        raise ValueError('stopped boundary boot changed')
    report = subprocess.check_output(['systemctl', '--user', 'show', service,
        '--property=Restart,MainPID,ActiveState'], text=True)
    values = dict(line.split('=', 1) for line in report.splitlines() if '=' in line)
    if values != {'Restart': 'no', 'MainPID': '0', 'ActiveState': 'inactive'}:
        raise ValueError('scheduler service must be inactive with Restart=no')
    for pid, ticks in request['tracked'].items():
        path = Path('/proc') / pid / 'stat'
        if path.exists() and int(path.read_text().rsplit(')', 1)[1].split()[19]) == ticks:
            raise ValueError('old tracked process still exists')
    # The host namespace proof makes absence meaningful, unlike sandbox /proc.
    for path in Path('/proc').glob('[0-9]*/cmdline'):
        try:
            arguments = path.read_bytes().split(b'\0')
        except FileNotFoundError:
            continue
        if any(Path(os.fsdecode(a)).name in ('schedule.py', 'perturbdiff_v2_runtime.py') for a in arguments if a):
            raise ValueError('scheduler process still exists')
    return {'pid_namespace': namespace, 'service': service, 'properties': values}


def reviewed_stage(stage, sha):
    if digest(stage / 'staging-sha256.json') != sha:
        raise ValueError('reviewed staging inventory SHA mismatch')
    hashes = read(stage / 'staging-sha256.json')
    if set(hashes) != NAMES or any(digest(stage / n) != h for n, h in hashes.items()):
        raise ValueError('staged document hash mismatch')
    if any(p.is_symlink() for p in stage.iterdir()):
        raise ValueError('staging symlink refused')
    return hashes


def verify_identities(formal):
    from artifacts import run_identity
    queue = read(formal / 'queue.json')
    expected = read(formal / 'expected_identities.json')
    for row in expected:
        job, = [j for j in queue['jobs'] if (j['model'], j['fold'], j['seed']) ==
                (row['model'], row['fold'], row['seed'])]
        if run_identity(queue['input_root'], job['model'], job['fold'], job['seed'], job.get('command')) != row['identity']:
            raise ValueError('source identity changed')


def fresh_check(formal, stopped, staged, hashes, output):
    from artifacts import run_identity
    with tempfile.TemporaryDirectory(prefix='tenv-v2-recheck-') as temporary:
        destination = Path(temporary) / 'stage'
        actual = stage_locked(formal, stopped, destination, read(formal / 'queue.json'), output, run_identity)
        if actual != hashes:
            raise ValueError('fresh migration differs from reviewed proposal')


def finish(journal_path, journal, host_check, fault=lambda phase: None):
    """Idempotent forward recovery; every rename is checked against durable bytes."""
    formal = Path(journal['formal'])
    replacement = Path(journal['replacement'])
    for operation in journal['moves']:
        source, destination = Path(operation['source']), Path(operation['destination'])
        if any(p.is_symlink() for p in (source, destination, *source.parents, *destination.parents)):
            raise ValueError('archive path symlink refused')
        if source.exists() and not destination.exists():
            if inventory(source) != operation['inventory']:
                raise ValueError('archive source changed')
            mkdir_durable(destination.parent)
            os.rename(source, destination)
            sync_directory(source.parent)
            sync_directory(destination.parent)
        elif not source.exists() and destination.exists():
            if inventory(destination) != operation['inventory']:
                raise ValueError('archive destination changed')
        else:
            raise ValueError('ambiguous archive state; no overwrite allowed')
        fault('archive_move')
    # Detect exchange completion by full directory inventories, including all
    # unchanged formal records. Crash between exchange and journal is recoverable.
    current, alternate = inventory(formal), inventory(replacement)
    if current == journal['old_formal'] and alternate == journal['new_formal']:
        verify_identities(replacement)
        host_check()
        exchange(formal, replacement)
        fault('exchange')
    elif current != journal['new_formal'] or alternate != journal['old_formal']:
        raise ValueError('formal generation changed or ambiguous')
    verify_identities(formal)
    journal['status'] = 'committed'
    publish(journal_path, journal)
    return journal


def transaction(formal, stopped, staged, sha, apply=False, namespace=None,
                service='arc-tenv-formal-scheduler.service'):
    formal, stopped, staged = formal.resolve(), stopped.resolve(), staged.resolve()
    hashes = reviewed_stage(staged, sha)
    output = Path(read(staged / 'queue.json')['output_root']).resolve()
    journal_path = output / 'perturbdiff-v2-apply.json'
    with (output / 'scheduling/.scheduler.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        host = lambda: host_stopped(stopped, namespace, service)
        if journal_path.exists():
            journal = read(journal_path)
            if (journal['reviewed_sha256'] != sha or journal['formal'] != str(formal) or
                    journal['stopped'] != str(stopped)):
                raise ValueError('transaction review or boundary mismatch')
            if not apply:
                return {'status': journal['status'], 'recovery_required': journal['status'] != 'committed'}
            host()
            if journal['status'] == 'committed':
                if inventory(formal) != journal['new_formal']:
                    raise ValueError('committed formal generation changed')
                verify_identities(formal)
                return journal
            return finish(journal_path, journal, host)
        fresh_check(formal, stopped, staged, hashes, output)
        if not apply:
            return {'status': 'verified_dry_run', 'reviewed_sha256': sha, 'activation_authorized': False}
        host_evidence = host()
        replacement = formal.with_name('.' + formal.name + '-v2-generation')
        if replacement.exists():
            raise ValueError('replacement generation exists without journal; inspect manually')
        old_formal = inventory(formal)
        shutil.copytree(formal, replacement)
        for name in NAMES - {'plan.json'}:
            shutil.copyfile(staged / name, replacement / name)
        shutil.copyfile(staged / 'plan.json', replacement / 'v2-migration-plan.json')
        sync_tree(replacement)
        moves = read(staged / 'plan.json')['archive_operations']
        for operation in moves:
            source, destination = Path(operation['source']), Path(operation['destination'])
            if destination.exists() or destination.is_symlink():
                raise ValueError('archive destination already exists')
            operation['inventory'] = inventory(source)
        journal = {'schema': 'perturbdiff-v2-apply-v1', 'status': 'prepared',
            'formal': str(formal), 'stopped': str(stopped), 'replacement': str(replacement),
            'reviewed_sha256': sha, 'staged_hashes': hashes, 'host_evidence': host_evidence,
            'old_formal': old_formal, 'new_formal': inventory(replacement), 'moves': moves,
            'restart_authorized': False, 'resume_supported': False}
        publish(journal_path, journal)
        return finish(journal_path, journal, host)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--formal', type=Path, required=True)
    parser.add_argument('--stopped', type=Path, required=True)
    parser.add_argument('--staged', type=Path, required=True)
    parser.add_argument('--reviewed-stage-sha256', required=True)
    parser.add_argument('--apply', action='store_true')
    parser.add_argument('--host-pid-namespace')
    parser.add_argument('--service', default='arc-tenv-formal-scheduler.service')
    args = parser.parse_args()
    print(json.dumps(transaction(args.formal, args.stopped, args.staged,
        args.reviewed_stage_sha256, args.apply, args.host_pid_namespace, args.service), sort_keys=True))
