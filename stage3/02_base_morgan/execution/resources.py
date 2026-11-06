"""Conservative host telemetry and monitoring of a runner-owned subprocess.

The caller must apply the repository's systemd memory cap before invoking this
module. Process reports deliberately contain no command arguments or environment.
"""
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import time

GIB = 1024 ** 3
RAM_FLOOR = 12 * GIB


class ResourceError(RuntimeError):
    """Resource policy or visibility prevents a safe launch."""


def _workload(args):
    """Inspect entry points, never inline code or arbitrary option values."""
    if not args:
        return None, False, False
    safe_path = re.compile(r'[A-Za-z0-9_./+-]+\.py')
    executable = Path(args[0]).name.lower()
    entry = None
    remaining = []
    if re.fullmatch(r'python(?:\d+(?:\.\d+)*)?', executable):
        skip_value = False
        for index, token in enumerate(args[1:], 1):
            if skip_value:
                skip_value = False
                continue
            if token == '-c':
                break  # Inline source may contain secrets or misleading paths.
            if token == '-m':
                if index + 1 < len(args) and re.fullmatch(r'[A-Za-z_][A-Za-z0-9_.]*', args[index + 1]):
                    entry = args[index + 1]
                    remaining = args[index + 2:]
                break
            if token in ('-W', '-X'):
                # These consume an option value; do not mistake it for a script.
                skip_value = True
                continue
            if token.startswith('-'):
                continue
            if safe_path.fullmatch(token):
                entry, remaining = token, args[index + 1:]
            break
    elif safe_path.fullmatch(args[0]):
        entry, remaining = args[0], args[1:]
    elif executable in ('torchrun', 'accelerate', 'deepspeed'):
        return None, True, True
    if entry is None:
        return None, False, False
    training = bool(re.search(r'(^|[/_.-])(train(?:er|ing)?|fit|worker|pretrain)([._/-]|$)', entry.lower()))
    # Public workflows use an immediate subcommand. Do not inspect option values.
    training |= bool(remaining and remaining[0] in ('train', 'fit', 'pretrain', 'benchmark'))
    repository_job = 'stage2/' in entry or 'stage3/' in entry
    return entry, repository_job or training, training


def _processes(proc_root):
    records = []
    readable = True
    for directory in proc_root.iterdir():
        if not directory.name.isdigit():
            continue
        try:
            args = directory.joinpath('cmdline').read_bytes().decode(errors='replace').split('\0')
            status = directory.joinpath('status').read_text()
            rss = re.search(r'^VmRSS:\s+(\d+)', status, re.M)
            parent = re.search(r'^PPid:\s+(\d+)', status, re.M)
            code, heavy, training = _workload(args)
            records.append(dict(pid=int(directory.name), ppid=int(parent[1]) if parent else 0,
                                code_path=code, rss_bytes=int(rss[1]) * 1024 if rss else 0,
                                heavy=heavy, long_training=training))
        except FileNotFoundError:
            continue  # Process exited during the scan.
        except (PermissionError, OSError):
            readable = False
    return records, readable


def collect_snapshot(exclude_pids=(), proc_root='/proc'):
    """Read bounded, argument-redacted telemetry; unknown visibility fails closed."""
    root = Path(proc_root)
    errors = []
    available = None
    try:
        match = re.search(r'^MemAvailable:\s+(\d+)\s+kB', (root/'meminfo').read_text(), re.M)
        available = int(match[1]) * 1024 if match else None
    except OSError:
        errors.append('memory_unobservable')
    try:
        processes, readable = _processes(root)
        init = (root/'1/comm').read_text().strip()
        host_visible = readable and init in ('systemd', 'init')
    except OSError:
        processes, host_visible = [], False
    excluded = set(exclude_pids) | {os.getpid()}
    # Exclude descendants of explicitly owned processes as well.
    while True:
        children = {p['pid'] for p in processes if p['ppid'] in excluded}
        if children <= excluded:
            break
        excluded |= children
    peers = [p for p in processes if p['pid'] not in excluded and p['heavy']]
    gpu = []
    gpu_visible = False
    try:
        result = subprocess.run(['nvidia-smi', '--query-gpu=index,memory.used,memory.free,utilization.gpu',
                                 '--format=csv,noheader,nounits'], capture_output=True, text=True,
                                timeout=3, check=True)
        for line in result.stdout.splitlines():
            index, used, free, utilization = [int(v.strip()) for v in line.split(',')]
            gpu.append(dict(index=index, used_mib=used, free_mib=free, utilization_percent=utilization))
        gpu_visible = bool(gpu)
    except (OSError, subprocess.SubprocessError, ValueError):
        errors.append('gpu_unobservable')
    return dict(timestamp=time.time(), memory_available_bytes=available,
                host_processes_observable=host_visible, gpu_observable=gpu_visible,
                gpu=gpu, peers=peers, processes=processes, errors=errors)


def preflight(snapshot, formal=True):
    """Enforce the RAM floor and repository-wide training serialization."""
    if snapshot.get('memory_available_bytes') is None or snapshot['memory_available_bytes'] < RAM_FLOOR:
        raise ResourceError('available_ram_below_12_gib')
    if not snapshot.get('host_processes_observable'):
        raise ResourceError('host_process_inventory_unobservable')
    if not snapshot.get('gpu_observable'):
        raise ResourceError('gpu_telemetry_unobservable')
    peers = snapshot.get('peers', [])
    if any(p['long_training'] for p in peers):
        raise ResourceError('external_long_training')
    if formal and peers:
        raise ResourceError('formal_training_requires_serial_execution')
    if len(peers) >= 2:
        raise ResourceError('two_heavy_jobs_already_running')


def _tree_rss(snapshot, pid):
    owned = {pid}
    while True:
        children = {p['pid'] for p in snapshot.get('processes', []) if p['ppid'] in owned}
        if children <= owned:
            return sum(p['rss_bytes'] for p in snapshot.get('processes', []) if p['pid'] in owned)
        owned |= children


def _terminate_owned(process):
    """Signal only our new session, tolerating a concurrent natural exit."""
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()


def run_guarded(command, log_dir, formal=True, interval=5, snapshot_fn=None):
    """Monitor and terminate only the process group created by this invocation.

    Commands are passed directly to Popen (no shell); their arguments are never
    written to telemetry. A hard memory cap remains the caller's responsibility.
    """
    if interval <= 0:
        raise ValueError('monitor interval must be positive')
    snapshot_fn = snapshot_fn or collect_snapshot
    directory = Path(log_dir)
    directory.mkdir(parents=True, exist_ok=True)
    peak = 0
    reason = None
    with (directory/'resources.jsonl').open('w') as log:
        initial = snapshot_fn(exclude_pids=())
        log.write(json.dumps(initial, sort_keys=True) + '\n')
        log.flush()
        preflight(initial, formal=formal)
        process = subprocess.Popen(command, start_new_session=True)
        try:
            while True:
                snapshot = snapshot_fn(exclude_pids=(process.pid,))
                peak = max(peak, _tree_rss(snapshot, process.pid))
                log.write(json.dumps(snapshot, sort_keys=True) + '\n')
                log.flush()
                try:
                    preflight(snapshot, formal=formal)
                except ResourceError as error:
                    reason = str(error)
                    _terminate_owned(process)
                    break
                if process.poll() is not None:
                    break
                time.sleep(interval)
        finally:
            if process.poll() is None:
                _terminate_owned(process)
    summary = dict(exit_code=process.returncode, stop_reason=reason,
                   peak_tree_rss_bytes=peak, formal=formal)
    (directory/'resource_summary.json').write_text(json.dumps(summary, sort_keys=True, indent=2) + '\n')
    return summary
