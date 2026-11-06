"""Read-only durable boundary evidence. Never signals, migrates, or authorizes dispatch.

Use in the actual host PID namespace. A running process's entire lifetime is
charged conservatively, including setup/inference. Completed attempts likewise
charge wall time. This is an upper bound, not exact native training time.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import time


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
    (destination / 'queue.json').write_bytes(queue_bytes)
    archive = destination / 'states'
    archive.mkdir()
    for name, raw in state_copies.items():
        with (archive / name).open('xb') as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
    with (destination / 'boundary.json').open('x') as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    fd = os.open(destination, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['capture'])
    parser.add_argument('--queue', required=True, type=Path)
    parser.add_argument('--destination', required=True, type=Path)
    args = parser.parse_args()
    capture(args.queue, args.destination)
