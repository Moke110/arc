#!/usr/bin/env python3
"""Dispatch Step 04 workers and record workstation admission and recovery."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import fcntl
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time
import uuid

REPO = next(path for path in Path(__file__).resolve().parents if (path / '.git').exists())
STEP03 = REPO / 'stage2/03_modes_tenv/execution'
# The Step 03 module has the same filename. Import it by file identity.
spec = importlib.util.spec_from_file_location('step03_scheduler', STEP03 / 'script/schedule.py')
step03 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(step03)
GIB = 1024 ** 3
SAMPLE_SECONDS = 2
PRESSURE_SAMPLES = 3


def read(path):
    return json.loads(path.read_text())


def number(value, minimum=0):
    return type(value) in (int, float) and math.isfinite(value) and value >= minimum


def load_queue(path):
    value = read(path)
    record = (path.parent / value['record_dir']).resolve()
    jobs = value['jobs']
    ids = set()
    for job in jobs:
        key = job['id']
        if not re.fullmatch(r'[A-Za-z0-9_-]{1,100}', key) or key in ids:
            raise ValueError('invalid or duplicate job id: ' + key)
        ids.add(key)
        command = job['command']
        if not isinstance(command, list) or not command or not all(
                isinstance(arg, str) and '\0' not in arg for arg in command):
            raise ValueError('job requires a command argument list: ' + key)
        if job['kind'] not in ('cpu', 'gpu', 'training'):
            raise ValueError('unknown worker kind: ' + key)
        memory = job.setdefault('memory_gib', 10)
        if memory not in (10, 18) or (memory == 18 and job['kind'] != 'training'):
            raise ValueError('18 GiB scope is reserved for exclusive formal training: ' + key)
        job['cwd'] = str((path.parent / job.get('cwd', str(REPO))).resolve())
        if not number(job.setdefault('disk_gib', 1)):
            raise ValueError('invalid required disk space: ' + key)
        if not number(job.setdefault('timeout_seconds', 14400), 1):
            raise ValueError('invalid worker timeout: ' + key)
        profile = job.get('profile')
        if profile is not None and (
                not number(profile.get('host_peak_gib'), 0.001)
                or profile['host_peak_gib'] > memory
                or not number(profile.get('gpu_peak_gib'))
                or (job['kind'] == 'cpu' and profile['gpu_peak_gib'] != 0)):
            raise ValueError('invalid measured profile: ' + key)
    if not jobs:
        raise ValueError('queue requires at least one job')
    return record, jobs


def live_sample(record, pids):
    snapshot = step03.process_snapshot(pids)
    if not snapshot['observable']:
        return snapshot
    gpu_util = subprocess.check_output(
        ['nvidia-smi', '--query-gpu=utilization.gpu', '--format=csv,noheader,nounits'],
        text=True, timeout=1)
    snapshot['gpu_utilization_percent'] = int(gpu_util.strip())
    mem = dict(line.split(':', 1) for line in Path('/proc/meminfo').read_text().splitlines())
    for label, key in (('host_total_gib', 'MemTotal'), ('swap_total_gib', 'SwapTotal'),
                       ('swap_free_gib', 'SwapFree')):
        snapshot[label] = int(mem[key].split()[0]) / 1024 ** 2
    disk = next(path for path in (record, *record.parents) if path.exists())
    snapshot['disk_free_gib'] = shutil.disk_usage(disk).free / GIB
    return snapshot


def snapshot(record, fixture=None, pids=()):
    started = time.time()
    started_monotonic = time.monotonic()
    try:
        if fixture:
            sample = read(fixture)
        else:
            # Bound the entire inherited Step 03 probe, including slow GPU calls.
            result = subprocess.check_output(
                [sys.executable, __file__, '_sample', str(record), *map(str, pids)],
                text=True, timeout=2.5)
            sample = json.loads(result)
        age = sample.pop('age_seconds', None)
        if age is None:
            age = started - sample.get('timestamp', started) if fixture else 0
        sample.setdefault('timestamp', started - age)
        sample['sample_started_monotonic'] = started_monotonic - age
    except (OSError, ValueError, TypeError, subprocess.SubprocessError) as exc:
        sample = dict(observable=False, error=str(exc), timestamp=started,
                      sample_started_monotonic=started_monotonic)
    sample['observed_at'] = time.time()
    sample['observed_monotonic'] = time.monotonic()
    return sample


def unsafe(sample):
    required = ('available_gib', 'gpu_free_gib', 'gpu_used_gib', 'gpu_utilization_percent',
                'disk_free_gib', 'external_heavy', 'external_peak_gib', 'timestamp',
                'sample_started_monotonic', 'observed_monotonic')
    if not sample.get('observable') or any(not number(sample.get(key)) for key in required):
        return 'missing telemetry'
    if not 0 <= sample['observed_monotonic'] - sample['sample_started_monotonic'] <= 5:
        return 'stale telemetry'
    if sample['available_gib'] < 12:
        return 'less than 12 GiB available RAM'
    if sample['gpu_free_gib'] < 0.5:
        return 'GPU memory below 0.5 GiB free'
    if sample['disk_free_gib'] < 1:
        return 'less than 1 GiB free disk'
    return None


def decision(sample, job, active=()):
    reason = unsafe(sample)
    if reason is None and sample['external_heavy']:
        reason = 'unprofiled external heavy job'
    if reason is None and active and (job['kind'] == 'training' or
                                     any(entry['kind'] == 'training' for entry in active)):
        reason = 'long training requires an exclusive slot'
    if reason is None and job['kind'] != 'cpu' and sample['gpu_utilization_percent'] >= 90:
        reason = 'GPU utilization at or above 90 percent'
    if reason is None and sample['disk_free_gib'] < 1 + job['disk_gib'] + sum(
            entry['disk_gib'] for entry in active):
        reason = 'insufficient free disk'
    if reason:
        return dict(admit=False, reason=reason)
    if not job.get('profile'):
        if active:
            return dict(admit=False, reason='unprofiled worker requires an exclusive slot')
        # Reserve the whole hard cap until an observed profile is available.
        gpu_reserve = 0 if job['kind'] == 'cpu' else 4
        return step03.admit(sample, dict(host_peak_gib=job['memory_gib'],
                                         gpu_peak_gib=gpu_reserve))
    return step03.admit(sample, job['profile'], active)


@contextmanager
def workstation_lock():
    common = subprocess.check_output(['git', '-C', str(REPO), 'rev-parse',
                                      '--path-format=absolute', '--git-common-dir'],
                                     text=True).strip()
    lock_path = Path(common).parent / 'stage2/03_modes_tenv/execution/.scheduler.lock'
    with lock_path.open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError('another Stage 2 scheduler holds the workstation lock') from exc
        yield


def load_status(record, jobs):
    path = record / 'status.json'
    identity = {job['id']: job for job in jobs}
    if not path.exists():
        return dict(jobs={key: dict(identity=job, status='pending', attempts=[])
                          for key, job in identity.items()}, last_decision={})
    status = read(path)
    if identity != {key: row['identity'] for key, row in status['jobs'].items()}:
        raise ValueError('queue identity changed; use a new record directory')
    status['jobs'] = {key: status['jobs'][key] for key in identity}
    return status


def publish(record, status):
    step03.save(record / 'status.json', status)


def handoff(record, status):
    return dict(status, paused=(record / 'pause').exists(),
                active=[key for key, row in status['jobs'].items() if row['status'] == 'running'],
                queue=[key for key, row in status['jobs'].items() if row['status'] == 'pending'])


def scope_state(unit):
    result = subprocess.run(['systemctl', '--user', 'show', unit, '-p', 'ActiveState', '--value'],
                            capture_output=True, text=True, timeout=2)
    state = result.stdout.strip()
    if state not in ('inactive', 'failed', 'active', 'activating', 'deactivating'):
        raise ValueError('cannot verify previous owned scope: ' + unit)
    return state


def process_identity(pid):
    try:
        proc = Path(f'/proc/{pid}/stat').read_text()
        # Linux stat field 22, after the last ')' that terminates comm.
        return proc.rsplit(')', 1)[1].split()[19]
    except FileNotFoundError:
        return None


def recover(record, status):
    for row in status['jobs'].values():
        if row['status'] != 'running':
            continue
        attempt = row['attempts'][-1]
        pid = attempt.get('pid')
        if (pid and process_identity(pid) == attempt['process_start']
                or scope_state(attempt['scope']) not in ('inactive', 'failed')):
            raise ValueError('previous owned worker still alive: ' + attempt['scope'])
        row['status'] = attempt['outcome'] = 'interrupted'
        attempt.update(reason='scheduler interrupted before worker outcome', finished_at=time.time())
    publish(record, status)


def prepare(record, row):
    job = row['identity']
    attempt = dict(scope=f"arc-step04-{uuid.uuid4().hex}.scope",
                   started_at=time.time(), memory_max_gib=job['memory_gib'],
                   memory_swap_max_gib=8, peak_host_rss_gib=0,
                   peak_workstation_gpu_used_gib=0, peak_gpu_utilization_percent=0,
                   log=str(record / f"{job['id']}.attempt{len(row['attempts']) + 1}.log"))
    row['attempts'].append(attempt)
    row['status'] = 'running'
    return attempt


def launch(job, attempt):
    command = ['systemd-run', '--user', '--scope', '--unit=' + attempt['scope'],
               '-p', f"MemoryMax={attempt['memory_max_gib']}G",
               '-p', f"MemorySwapMax={attempt['memory_swap_max_gib']}G", '--',
               '/usr/bin/time', '-v', '--', *job['command']]
    environment = dict(os.environ, OMP_NUM_THREADS='2', MKL_NUM_THREADS='2',
                       OPENBLAS_NUM_THREADS='2', NUMEXPR_NUM_THREADS='2')
    with Path(attempt['log']).open('w') as log:
        process = subprocess.Popen(command, cwd=job['cwd'], env=environment,
                                   stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
    attempt.update(pid=process.pid, process_start=process_identity(process.pid))
    return dict(process=process, started=time.monotonic())


def request_stop(entry, attempt, reason):
    if 'stop_at' in entry:
        return
    entry['stop_at'] = time.monotonic()
    attempt['reason'] = reason
    try:
        os.killpg(entry['process'].pid, signal.SIGTERM)
    except ProcessLookupError:
        pass


def finalize(row, returncode):
    attempt = row['attempts'][-1]
    reason = attempt.get('reason')
    outcome = ('interrupted' if reason and reason != 'worker timeout' else
               'timeout' if reason == 'worker timeout' else
               'completed' if returncode == 0 else 'failed')
    log = Path(attempt['log'])
    if log.exists():
        with log.open('rb') as stream:
            stream.seek(max(0, log.stat().st_size - 8192))
            tail = stream.read().decode(errors='replace')
        matches = re.findall(r'Maximum resident set size \(kbytes\):\s*(\d+)', tail)
        if matches:
            attempt['peak_host_rss_gib'] = max(attempt['peak_host_rss_gib'],
                                               int(matches[-1]) / 1024 ** 2)
    attempt.update(outcome=outcome, returncode=returncode, finished_at=time.time(),
                   reason=reason or f'process exit {returncode}')
    row['status'] = outcome


def run(record, jobs, fixture):
    with workstation_lock():
        status = load_status(record, jobs)
        rows = status['jobs']
        recover(record, status)
        active = {}
        pressure_streak = 0
        blocked_since = None
        interrupted = False

        def interrupt(_signum, _frame):
            nonlocal interrupted
            interrupted = True

        signal.signal(signal.SIGTERM, interrupt)
        signal.signal(signal.SIGINT, interrupt)
        status['scheduler_pid'] = os.getpid()
        try:
            while True:
                tick = time.monotonic()
                sample = snapshot(record, fixture,
                                  [entry['process'].pid for entry in active.values()])
                sample['active'] = list(active)
                with (record / 'telemetry.jsonl').open('a') as stream:
                    stream.write(json.dumps(sample) + '\n')
                pressure = unsafe(sample)
                if pressure is None:
                    external = sample['external_heavy']
                    owned_rss = sum(sample.get('owned_rss_gib', {}).values())
                    if external + len(active) > 2:
                        pressure = 'repo-wide heavy worker limit exceeded'
                    elif external and any(rows[key]['identity']['kind'] == 'training' for key in active):
                        pressure = 'external heavy worker conflicts with exclusive training'
                    elif owned_rss + sample['external_peak_gib'] > 28:
                        pressure = 'observed heavy worker RSS exceeds 28 GiB'
                pressure_streak = pressure_streak + 1 if pressure else 0
                status['pressure'] = dict(reason=pressure, consecutive_samples=pressure_streak)
                if active and pressure_streak >= PRESSURE_SAMPLES:
                    victim = next(reversed(active))
                    request_stop(active[victim], rows[victim]['attempts'][-1], pressure)
                    pressure_streak = 0
                for key, entry in list(active.items()):
                    row = rows[key]
                    attempt = row['attempts'][-1]
                    process = entry['process']
                    attempt['peak_host_rss_gib'] = max(
                        attempt['peak_host_rss_gib'],
                        sample.get('owned_rss_gib', {}).get(str(process.pid), 0))
                    for peak, metric in (
                            ('peak_workstation_gpu_used_gib', 'gpu_used_gib'),
                            ('peak_gpu_utilization_percent', 'gpu_utilization_percent')):
                        if number(sample.get(metric)):
                            attempt[peak] = max(attempt[peak], sample[metric])
                    code = process.poll()
                    if code is None and interrupted:
                        request_stop(entry, attempt, 'scheduler interrupted')
                    if code is None and tick - entry['started'] >= row['identity']['timeout_seconds']:
                        request_stop(entry, attempt, 'worker timeout')
                    if code is None and tick - entry.get('stop_at', tick) >= 15:
                        try:
                            os.killpg(process.pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                    if code is not None:
                        finalize(row, code)
                        del active[key]

                paused = (record / 'pause').exists()
                pending = [key for key, row in rows.items() if row['status'] == 'pending']
                for key in pending:
                    job = rows[key]['identity']
                    paused = (record / 'pause').exists()
                    if paused or interrupted:
                        verdict = dict(admit=False, reason='paused' if paused else 'scheduler interrupted')
                    elif pressure:
                        verdict = dict(admit=False, reason=pressure)
                    else:
                        verdict = decision(sample, job, [rows[name]['identity'] for name in active])
                    status['last_decision'] = dict(verdict, job=key, timestamp=time.time())
                    if not verdict['admit']:
                        continue
                    attempt = prepare(record, rows[key])
                    publish(record, status)
                    try:
                        active[key] = launch(job, attempt)
                    except OSError as exc:
                        attempt['launch_error'] = str(exc)
                        finalize(rows[key], 126)
                    # Use a new live observation before another dispatch.
                    break
                publish(record, status)
                pending = any(row['status'] == 'pending' for row in rows.values())
                if not active:
                    if paused or interrupted or not pending:
                        return 130 if interrupted else 0
                    blocked_since = tick if blocked_since is None else blocked_since
                    if tick - blocked_since >= 300:
                        return 2
                else:
                    blocked_since = None
                time.sleep(max(0, SAMPLE_SECONDS - (time.monotonic() - tick)))
        finally:
            for key, entry in active.items():
                row = rows[key]
                request_stop(entry, row['attempts'][-1], 'scheduler interrupted')
                process = entry['process']
                try:
                    code = process.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    code = process.wait(timeout=5)
                finalize(row, code)
            status['scheduler_pid'] = None
            publish(record, status)


def main():
    if sys.argv[1:2] == ['_sample']:
        try:
            print(json.dumps(live_sample(Path(sys.argv[2]), tuple(map(int, sys.argv[3:])))))
        except (OSError, ValueError, subprocess.SubprocessError) as exc:
            print(json.dumps(dict(observable=False, error=str(exc))))
        return 0
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    for name in ('preflight', 'run', 'status', 'pause', 'resume', 'retry'):
        command = sub.add_parser(name)
        command.add_argument('queue', type=Path)
        if name in ('preflight', 'run'):
            command.add_argument('--fixture', type=Path, help=argparse.SUPPRESS)
        if name == 'retry':
            command.add_argument('job_id')
    args = parser.parse_args()
    try:
        record, jobs = load_queue(args.queue.resolve())
        if args.action == 'preflight':
            sample = snapshot(record, args.fixture)
            decisions = {job['id']: decision(sample, job) for job in jobs}
            print(json.dumps(dict(sample=sample, decisions=decisions, jobs=jobs), indent=2))
            return 0 if any(item['admit'] for item in decisions.values()) else 2
        if args.action == 'run':
            return run(record, jobs, args.fixture)
        status = load_status(record, jobs)
        if args.action == 'pause':
            record.mkdir(parents=True, exist_ok=True)
            (record / 'pause').touch()
        elif args.action == 'resume':
            (record / 'pause').unlink(missing_ok=True)
        elif args.action == 'retry':
            with workstation_lock():
                status = load_status(record, jobs)
                recover(record, status)
                row = status['jobs'][args.job_id]
                if row['status'] not in ('interrupted', 'failed', 'timeout'):
                    raise ValueError('only an unsuccessful attempt can be retried')
                row['status'] = 'pending'
                publish(record, status)
        print(json.dumps(handoff(record, status), indent=2))
        return 0
    except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as exc:
        parser.exit(2, f'{parser.prog}: {exc}\n')


if __name__ == '__main__':
    sys.exit(main())
