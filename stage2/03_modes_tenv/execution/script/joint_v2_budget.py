"""Offline-prepared v2 cumulative account and detached worker watchdog.

Not integrated into the live queue. `initialize` requires a reviewed, final
migration ledger digest. Historical and worker wall time conservatively include
setup/inference, because an independently authenticated training phase is absent.
An interrupted watchdog leaves the entire remaining budget reserved (fail closed).
"""
import argparse
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

CAP = 14400


def debit(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise ValueError('missing or invalid cumulative debit')
    return float(value)


def write(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('x') as stream:
        json.dump(value, stream, sort_keys=True, indent=2)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)
    fd = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def read(path):
    return json.loads(Path(path).read_text())


def initialize(source, digest, key, account):
    raw = Path(source).read_bytes()
    if hashlib.sha256(raw).hexdigest() != digest:
        raise ValueError('migration ledger digest mismatch')
    migration = json.loads(raw)
    if migration.get('schema') != 'joint-v2-migration-budget-v1' or migration.get('boundary_verified') is not True:
        raise ValueError('reviewed final migration ledger required')
    row = migration['ledger'][key]
    if row.get('reconciliation_required', True) or row.get('cap_seconds') != CAP:
        raise ValueError('unreconciled debit or cap reset')
    consumed = debit(row.get('consumed_seconds'))
    if not row.get('evidence'):
        raise ValueError('durable debit evidence required, including zero debit')
    account = Path(account)
    if str(account.resolve()) != row.get('account_path'):
        raise ValueError('account path must match protected migration ledger')
    if not isinstance(row.get('identity'), dict):
        raise ValueError('reviewed v2 identity required')
    account.mkdir(exist_ok=False)
    # Immutable local copy plus caller-pinned digest is retained for audit.
    with (account / 'migration.json').open('xb') as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())
    write(account / 'account.json', dict(schema='joint-v2-account-v1', key=key,
        cap_seconds=CAP, migration_sha256=digest, consumed_seconds=consumed,
        attempts=[], accounting='conservative full worker wall; includes setup and inference'))


def load_account(account):
    state = read(account / 'account.json')
    if state['schema'] != 'joint-v2-account-v1' or state['cap_seconds'] != CAP:
        raise ValueError('invalid account or cap reset')
    raw = (account / 'migration.json').read_bytes()
    if hashlib.sha256(raw).hexdigest() != state['migration_sha256']:
        raise ValueError('protected migration ledger changed')
    migration_row = json.loads(raw)['ledger'][state['key']]
    if str(account.resolve()) != migration_row['account_path']:
        raise ValueError('account moved from protected migration path')
    initial = debit(migration_row['consumed_seconds'])
    consumed = debit(state['consumed_seconds'])
    expected = initial
    for attempt in state['attempts']:
        if attempt['prior_training_seconds'] != expected:
            raise ValueError('attempt debit chain reset')
        if attempt['status'] == 'exited':
            expected = min(CAP, expected + debit(attempt['wall_seconds']))
        else:
            expected = CAP
    if consumed != expected:
        raise ValueError('historical or attempt debit reset')
    return state


def process_identity(pid):
    fields = Path(f'/proc/{pid}/stat').read_text().rsplit(')', 1)[1].split()
    return dict(pid=pid, start_ticks=int(fields[19]),
        boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip())


def matches(identity):
    try:
        return process_identity(identity['pid']) == identity
    except (FileNotFoundError, ProcessLookupError):
        return False


def supervise(account, request_path, command):
    """Run from a detached session/service; hold one account lock across the attempt.

    Full remaining credit is reserved durably BEFORE spawn. SIGKILL, boot loss or
    disk failure therefore cannot create reusable credit. No resume is supported.
    """
    account = Path(account)
    with (account / 'lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        state = load_account(account)
        consumed = debit(state['consumed_seconds'])
        remaining = CAP - consumed
        if remaining <= 0:
            raise ValueError('cumulative budget exhausted or outstanding reservation')
        request = read(request_path)
        if any(request.get(k) is not None for k in ('resume', 'resume_checkpoint', 'restart_checkpoint')):
            raise ValueError('native restart state unsupported')
        if 'prior_training_seconds' in request:
            raise ValueError('scheduler must not supply its own debit')
        migration_row = read(account / 'migration.json')['ledger'][state['key']]
        if request['identity'] != migration_row['identity']:
            raise ValueError('request identity differs from protected migration')
        if request['identity'].get('model') not in ('perturbdiff','cmonge'):
            raise ValueError('unsupported v2 model')
        expected = f"{request['identity']['model']}_fold{request['identity']['fold']}_seed{request['identity']['seed']}"
        if expected != state['key']:
            raise ValueError('request fold/seed mismatch')
        inference_limit = float(request.get('inference_timeout_seconds', 1800))
        if not math.isfinite(inference_limit) or not 0 < inference_limit <= 1800:
            raise ValueError('inference timeout must be finite and at most 1800 seconds')
        request['prior_training_seconds'] = consumed
        attempt_index = len(state['attempts'])
        injected = account / f'request-{attempt_index}.json'
        write(injected, request)
        if command.count('{request}') != 1:
            raise ValueError('exactly one literal {request} command argument required')
        units = [a.split('=',1)[1] for a in command if a.startswith('--unit=')]
        scope = units[0] if command[:3] == ['systemd-run','--user','--scope'] and len(units) == 1 else None
        if scope is not None and (not scope.startswith('arc-tenv-v2-') or not scope.endswith('.scope') or '/' in scope):
            raise ValueError('unowned worker scope')
        def kill_scope():
            if scope:
                result = subprocess.run(['systemctl','--user','kill','--kill-whom=all','--signal=SIGKILL',scope],
                    capture_output=True,timeout=5)
                attempt['scope_kill_returncode'] = result.returncode
                deadline = time.monotonic() + 5
                while True:
                    report = subprocess.run(['systemctl','--user','show',scope,'--property=ActiveState,LoadState'],
                        capture_output=True,text=True,timeout=5)
                    values = dict(line.split('=',1) for line in report.stdout.splitlines() if '=' in line)
                    if values.get('ActiveState') in ('inactive','failed') or values.get('LoadState') == 'not-found':
                        break
                    if time.monotonic() >= deadline:
                        raise ValueError('owned scope exit unverified; reservation remains spent')
                    time.sleep(.05)
        argv = [str(injected.resolve()) if arg == '{request}' else arg for arg in command]
        attempt = dict(prior_training_seconds=consumed, reserved_seconds=remaining,
                       status='reserved', request=str(injected), command=argv)
        state['attempts'].append(attempt)
        state['consumed_seconds'] = CAP
        write(account / 'account.json', state)
        start = time.monotonic()
        worker = None
        try:
            worker = subprocess.Popen(argv, start_new_session=True)
            identity = process_identity(worker.pid)
            attempt.update(status='running', process=identity)
            write(account / 'account.json', state)
            # Popen owns an unreaped child. A live child cannot have its PID
            # recycled; poll returning None does not reap it. Kill the group
            # before wait, preserving this ownership guarantee.
            inference_started = None
            inference_expired = False
            while worker.poll() is None and time.monotonic() - start < remaining:
                if request.get('output_directory'):
                    try:
                        progress = read(Path(request['output_directory']) / 'progress.json')
                    except (OSError, ValueError):
                        progress = {}
                    if progress.get('phase') == 'inference':
                        inference_started = inference_started or time.monotonic()
                    if inference_started is not None and time.monotonic() - inference_started >= inference_limit:
                        inference_expired = True
                        break
                time.sleep(min(.1, max(.001, remaining - (time.monotonic() - start))))
            alive = worker.poll() is None
            timed_out = alive or inference_expired or time.monotonic() - start >= remaining
            if alive:
                if not matches(identity):
                    raise ValueError('worker process identity changed')
                os.killpg(worker.pid, signal.SIGKILL)
                kill_scope()
            returncode = worker.wait()
            elapsed = time.monotonic() - start
            # A timeout remains terminal even if an adapter emitted artifacts.
            attempt.update(status='timeout' if timed_out else 'exited', returncode=returncode,
                           wall_seconds=elapsed, stop_reason='inference_timeout' if inference_expired else ('timeout' if timed_out else 'process_exit'))
            state['consumed_seconds'] = CAP if timed_out else min(CAP, consumed + elapsed)
            write(account / 'account.json', state)
            return 124 if timed_out else returncode
        except BaseException:
            # Reservation remains spent on any ambiguous outcome; kill only owned child.
            if worker is not None and worker.poll() is None:
                os.killpg(worker.pid, signal.SIGKILL)
                kill_scope()
                worker.wait()
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    init = sub.add_parser('initialize')
    for name in ('source', 'sha256', 'key', 'account'):
        init.add_argument('--' + name, required=True)
    check = sub.add_parser('inspect')
    check.add_argument('--account', required=True)
    launch = sub.add_parser('launch')
    launch.add_argument('--account', required=True)
    launch.add_argument('--request', required=True)
    launch.add_argument('command', nargs=argparse.REMAINDER)
    run = sub.add_parser('supervise')
    run.add_argument('--account', required=True)
    run.add_argument('--request', required=True)
    run.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.action == 'initialize':
        initialize(args.source, args.sha256, args.key, args.account)
    elif args.action == 'inspect':
        print(json.dumps(load_account(Path(args.account)), sort_keys=True))
    elif args.action == 'launch':
        account = Path(args.account).resolve()
        load_account(account)
        command = args.command[1:] if args.command[:1] == ['--'] else args.command
        with (account / 'watchdog.log').open('ab', buffering=0) as log:
            child = subprocess.Popen([sys.executable, str(Path(__file__).resolve()),
                'supervise', '--account', str(account), '--request', str(Path(args.request).resolve()),
                '--', *command], stdin=subprocess.DEVNULL, stdout=log, stderr=log,
                start_new_session=True, close_fds=True)
        print(json.dumps({'watchdog_pid': child.pid, 'dispatch_confirmed': False}))
    else:
        # Detach from scheduler session; caller must launch with stdin/out/err files
        # and close_fds=True (or a separate systemd service). Never inherit pipes.
        if os.getsid(0) != os.getpid():
            raise ValueError('watchdog must be launched in its own session')
        command = args.command[1:] if args.command[:1] == ['--'] else args.command
        sys.exit(supervise(args.account, args.request, command))


if __name__ == '__main__':
    main()
