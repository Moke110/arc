"""Durable capture and explicitly invoked SIGINT boundary stop; no migration or dispatch.

Use in the actual host PID namespace. A running process's entire lifetime is
charged conservatively, including setup/inference. Completed attempts likewise
charge wall time. This is an upper bound, not exact native training time.
"""
import argparse
import ctypes
import hashlib
import json
import math
import os
import signal
from pathlib import Path
import time


def sync_directory(path):
    fd = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def seconds(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise ValueError('invalid debit')
    return value


def process_identity(pid):
    # /proc stat command may itself contain spaces or parentheses.
    stat = Path(f'/proc/{int(pid)}/stat').read_text().rsplit(')', 1)[1].split()
    start = int(stat[19])
    uptime = float(Path('/proc/uptime').read_text().split()[0])
    return dict(pid=int(pid), start_ticks=start,
                boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
                lifetime_seconds=max(0, uptime - start / os.sysconf('SC_CLK_TCK')))


def capture(queue_path, destination):
    queue_bytes = queue_path.read_bytes()
    queue = json.loads(queue_bytes)
    output = Path(queue['output_root']).resolve()
    destination = destination.resolve()
    if destination == output or output in destination.parents:
        raise ValueError('capture destination must be outside live output')
    jobs = queue['jobs']
    selected = [j for j in jobs if j['model'] == 'perturbdiff']
    expected = {(f, s) for f in range(3) for s in (17, 42, 2026)}
    if len(selected) != 9 or {(j['fold'], j['seed']) for j in selected} != expected:
        raise ValueError('exactly nine unique PerturbDiff slots required')
    ledger = {}
    sources = {}
    state_copies = {}
    for job in selected:
        key = f"perturbdiff_fold{job['fold']}_seed{job['seed']}"
        path = output / 'scheduling' / (key + '.json')
        if not path.exists():
            ledger[key] = dict(consumed_seconds=None, remaining_seconds=None,
                               reconciliation_required=True, reason='never-dispatched status unverified')
            continue
        raw = path.read_bytes()
        state = json.loads(raw)
        sources[str(path.relative_to(output))] = hashlib.sha256(raw).hexdigest()
        state_copies[path.name] = raw
        sampled = seconds(state['training_seconds'])
        attempts = state.get('attempts', [])
        wall = sum(seconds(a['wall_seconds']) for a in attempts)
        process = None
        if state['status'] == 'running':
            process = process_identity(state['pid'])  # fails closed when PID is hidden/dead
            wall += process['lifetime_seconds']
        consumed = max(sampled, wall)
        ledger[key] = dict(consumed_seconds=consumed, remaining_seconds=max(0, 14400-consumed),
            sampled_training_seconds=sampled, attempt_wall_upper_bound_seconds=wall,
            legacy_identity=state['identity'], process=process,
            reconciliation_required=True, accounting='conservative wall upper bound; includes inference')
    # Hash binaries and histories without deserializing. A live capture remains provisional.
    inventory = {}
    legacy = output / 'runs' / 'perturbdiff'
    for path in sorted(legacy.rglob('*')):
        if path.is_symlink():
            raise ValueError('legacy symlink refused: ' + str(path))
        if path.is_file():
            before = path.stat()
            sha = digest(path)
            after = path.stat()
            if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                raise ValueError('legacy artifact changed during capture: ' + str(path))
            inventory[str(path.relative_to(output))] = dict(sha256=sha, bytes=after.st_size)
    result = dict(schema='perturbdiff-boundary-capture-v1', activation_authorized=False,
        captured_unix_seconds=time.time(), queue_sha256=hashlib.sha256(queue_bytes).hexdigest(),
        ledger=ledger, source_state_hashes=sources, legacy_inventory=inventory,
        other_jobs=[j for j in jobs if j['model'] != 'perturbdiff'],
        limitations=['Live snapshot is provisional; repeat after verified process exit.',
                    'Absent state is not proof of zero debit.',
                    'No SIGTERM handler exists in legacy scheduler; SIGTERM does not run finally.',
                    'No migration transaction or budget watchdog is provided by this command.'])
    destination.mkdir(parents=True, exist_ok=False)
    with (destination / 'queue.json').open('xb') as stream:
        stream.write(queue_bytes)
        stream.flush()
        os.fsync(stream.fileno())
    archive = destination / 'states'
    archive.mkdir()
    for name, raw in state_copies.items():
        with (archive / name).open('xb') as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
    sync_directory(archive)
    with (destination / 'boundary.json').open('x') as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    sync_directory(destination)
    sync_directory(destination.parent)


def process_table():
    records = {}
    for directory in Path('/proc').iterdir():
        if not directory.name.isdecimal():
            continue
        try:
            fields = (directory / 'stat').read_text().rsplit(')', 1)[1].split()
            records[int(directory.name)] = dict(state=fields[0], parent=int(fields[1]),
                group=int(fields[2]), start_ticks=int(fields[19]))
        except FileNotFoundError:
            continue
    return records


def verify_request(request):
    current = process_identity(request['pid'])
    for key in ('pid', 'start_ticks', 'boot_id'):
        if current[key] != request[key]:
            raise ValueError('scheduler identity mismatch: ' + key)
    proc = Path('/proc') / str(request['pid'])
    if digest(proc / 'cmdline') != request['cmdline_sha256']:
        raise ValueError('scheduler command identity mismatch')
    status = dict(line.split(':', 1) for line in (proc / 'status').read_text().splitlines() if ':' in line)
    # A shell background process may inherit ignored SIGINT. Never claim graceful
    # interruption unless SIGINT is caught and unblocked in the scheduler.
    bit = 1 << (signal.SIGINT - 1)
    if int(status['SigIgn'], 16) & bit or int(status['SigBlk'], 16) & bit or not int(status['SigCgt'], 16) & bit:
        raise ValueError('scheduler SIGINT is not caught and unblocked')
    return current


def stop(queue_path, destination, request, timeout):
    if not math.isfinite(timeout) or timeout <= 0 or timeout > 60:
        raise ValueError('timeout must be in (0, 60]')
    if request['pid'] == os.getpid():
        raise ValueError('refusing self interruption')
    verify_request(request)
    # pidfd binds signal delivery to this process, preventing PID-reuse races.
    libc = ctypes.CDLL(None, use_errno=True)
    open_pidfd = libc.pidfd_open
    open_pidfd.argtypes = [ctypes.c_int, ctypes.c_uint]
    open_pidfd.restype = ctypes.c_int
    pidfd = open_pidfd(request['pid'], 0)
    if pidfd < 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error))
    try:
        verify_request(request)
        destination.mkdir(parents=True, exist_ok=False)
        debit_started = time.monotonic()
        capture(queue_path, destination / 'before')
        verify_request(request)
        table = process_table()
        owned = {request['pid']}
        while True:
            expanded = owned | {pid for pid, row in table.items() if row['parent'] in owned}
            if expanded == owned:
                break
            owned = expanded
        tracked = {pid: table[pid]['start_ticks'] for pid in owned if pid in table}
        groups = {table[pid]['group'] for pid in owned if pid != request['pid']}
        # Persist authorization target and initial ownership before signal delivery.
        with (destination / 'request.json').open('x') as stream:
            json.dump(dict(target=request, tracked=tracked, groups=sorted(groups)), stream, indent=2)
            stream.flush()
            os.fsync(stream.fileno())
        sync_directory(destination)
        sync_directory(destination.parent)
        started = time.monotonic()
        send_signal = libc.pidfd_send_signal
        send_signal.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_void_p, ctypes.c_uint]
        send_signal.restype = ctypes.c_int
        if send_signal(pidfd, signal.SIGINT, None, 0) != 0:
            error = ctypes.get_errno()
            raise OSError(error, os.strerror(error))
        while True:
            now = process_table()
            survivors = [pid for pid, row in now.items() if row['state'] != 'Z' and
                (tracked.get(pid) == row['start_ticks'] or row['group'] in groups)]
            if not survivors:
                break
            if time.monotonic() - started >= timeout:
                raise RuntimeError('stop unverified; live processes remain: ' + repr(survivors))
            time.sleep(.05)
        elapsed = time.monotonic() - debit_started
        capture(queue_path, destination / 'after')
        before = json.loads((destination / 'before/boundary.json').read_text())
        after = json.loads((destination / 'after/boundary.json').read_text())
        if before['queue_sha256'] != after['queue_sha256']:
            raise ValueError('queue changed during stop')
        ledger = {}
        for key, prior in before['ledger'].items():
            final = after['ledger'][key]
            values = [v for v in (prior['consumed_seconds'], final['consumed_seconds']) if v is not None]
            if prior.get('process'):
                values.append(prior['consumed_seconds'] + elapsed)
            debit = max(values) if values else None
            ledger[key] = dict(consumed_seconds=debit,
                remaining_seconds=None if debit is None else max(0, 14400 - debit),
                reconciliation_required=True)
        with (destination / 'stop.json').open('x') as stream:
            json.dump(dict(process_exit_verified=True, activation_authorized=False,
                signal='SIGINT', stop_elapsed_seconds=elapsed, ledger=ledger,
                accounting='max pre/post sampled or attempt-wall debit; live pre debit plus stop elapsed; includes setup/inference; unresolved absent states',
                limitations=['No lock is held after exit; coordinator must prevent service restart and acquire scheduler lock before migration.',
                             'No budget watchdog or migration is installed.']), stream, indent=2)
            stream.flush()
            os.fsync(stream.fileno())
        sync_directory(destination)
    finally:
        os.close(pidfd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['capture', 'stop'])
    parser.add_argument('--queue', required=True, type=Path)
    parser.add_argument('--destination', required=True, type=Path)
    parser.add_argument('--request', type=Path)
    parser.add_argument('--timeout', type=float, default=45)
    args = parser.parse_args()
    if args.action == 'stop':
        if args.request is None:
            parser.error('stop requires --request with pinned process identity')
        stop(args.queue, args.destination, json.loads(args.request.read_text()), args.timeout)
    else:
        capture(args.queue, args.destination)
