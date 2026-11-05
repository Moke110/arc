#!/usr/bin/env python
"""Own experiment processes and resource admission; telemetry is injectable."""
from __future__ import annotations

import argparse
import fcntl
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

TERMINAL = {'completed', 'failed', 'timeout', 'excluded'}


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    temporary.replace(path)


def transition(state, reason, seconds):
    """Account attempts against a single cumulative four-hour allowance."""
    state = dict(state)
    state.setdefault('training_seconds', 0)
    state.setdefault('environment_retries', 0)
    state.setdefault('oom_retries', 0)
    state.setdefault('batch_scale', 1.0)
    if state.get('status') in TERMINAL:
        return state
    if seconds < 0:
        raise ValueError('elapsed seconds must be nonnegative')
    state['training_seconds'] += seconds
    state['stop_reason'] = reason
    state['status'] = 'failed'
    if state['training_seconds'] >= 14400 or reason in ('timeout', 'inference_timeout'):
        state['status'] = 'timeout'
    elif reason == 'completed':
        state['status'] = 'completed'
    elif reason == 'environment' and state['environment_retries'] < 1:
        state['environment_retries'] += 1
        state['status'] = 'pending'
    elif reason == 'oom' and state['oom_retries'] < 1:
        state['oom_retries'] += 1
        state['batch_scale'] *= 0.5
        state['status'] = 'pending'
    return state


def admit(snapshot, profile=None, active=()):
    """Require known resources before overlapping work, including unrelated jobs."""
    def decision(ok, reason):
        return {'admit': ok, 'reason': reason}
    if not snapshot.get('observable', False):
        return decision(False, 'host/GPU process telemetry unavailable')
    if snapshot['available_gib'] < 12:
        return decision(False, 'less than 12 GiB available RAM')
    external = snapshot['external_heavy']
    if len(active) + external >= 2:
        return decision(False, 'two heavy tasks already active')
    if not profile:
        return decision(not active and not external and snapshot['gpu_free_gib'] >= 1,
                        'unknown profile starts alone with GPU headroom')
    if any(not task.get('profile') for task in active):
        return decision(False, 'active task has unknown profile')
    peak = profile['host_peak_gib'] * 1.2
    combined = snapshot['external_peak_gib'] + peak
    combined += sum(task['profile']['host_peak_gib'] * 1.2 for task in active)
    if combined > 28 or snapshot['available_gib'] - peak < 12:
        return decision(False, 'host peak plus 20 percent headroom exceeds budget')
    if profile.get('gpu_peak_gib', 0) * 1.2 + 1 > snapshot['gpu_free_gib']:
        return decision(False, 'insufficient GPU headroom')
    return decision(True, 'resource budgets satisfied')


def pressure(history, tasks):
    if len(history) < 2:
        return {'stop': None}
    previous, current = history[-2:]
    worsening = (current['available_gib'] < 12 and
                 current['available_gib'] < previous['available_gib'])
    worsening |= (current['gpu_free_gib'] < 0.5 and
                  current['gpu_free_gib'] < previous['gpu_free_gib'])
    owned = [task for task in tasks if task.get('owned')]
    target = max(owned, key=lambda task: task['started']) if worsening and owned else None
    return {'stop': target['id'] if target else None}


def process_snapshot(owned_pids=(), idle_pids=()):
    """Read the real workstation; a restricted PID namespace is not evidence."""
    try:
        init = Path('/proc/1/comm').read_text().strip()
        if init not in ('systemd', 'init', 'wsl-init'):
            raise RuntimeError('restricted PID namespace: ' + init)
        memory = dict(line.split(':', 1) for line in Path('/proc/meminfo').read_text().splitlines())
        available = int(memory['MemAvailable'].split()[0]) / 1024 ** 2
        gpu = subprocess.check_output(
            ['nvidia-smi', '--query-gpu=memory.total,memory.used',
             '--format=csv,noheader,nounits'], text=True, timeout=5)
        total, used = map(float, gpu.splitlines()[0].split(','))
        compute = subprocess.check_output(
            ['nvidia-smi', '--query-compute-apps=pid', '--format=csv,noheader,nounits'],
            text=True, timeout=5)
        gpu_pids = {int(line.strip()) for line in compute.splitlines() if line.strip().isdigit()}
        processes = {}
        parents = {}
        for directory in Path('/proc').iterdir():
            if not directory.name.isdigit():
                continue
            try:
                status = dict(line.split(':', 1) for line in
                              (directory / 'status').read_text().splitlines())
                pid = int(directory.name)
                parents[pid] = int(status['PPid'])
                rss = int(status.get('VmRSS', '0 kB').split()[0]) / 1024 ** 2
                command = (directory / 'cmdline').read_bytes().replace(b'\0', b' ').decode(errors='replace')
                # Compute interpreters, GPU clients, and capped scopes are tasks;
                # workstation browsers are included in MemAvailable, not task count.
                executable = command.split(' ', 1)[0].rsplit('/', 1)[-1]
                heavy = pid in gpu_pids or (rss >= 0.5 and
                        (executable.startswith(('python', 'R', 'julia'))))
                processes[pid] = (rss, heavy)
            except (OSError, ValueError, KeyError):
                continue  # processes may exit between directory and status reads
        def belongs(pid, roots):
            seen = set()
            while pid and pid not in seen:
                if pid in roots:
                    return True
                seen.add(pid)
                pid = parents.get(pid, 0)
            return False
        excluded = set(owned_pids) | set(idle_pids)
        external = [pid for pid, (_, heavy) in processes.items()
                    if heavy and not belongs(pid, excluded)]
        # Child workers of a training interpreter form one heavy task.
        roots = [pid for pid in external if not belongs(parents.get(pid, 0), set(external))]
        external_peak = sum(rss for pid, (rss, _) in processes.items() if belongs(pid, set(roots)))
        owned_rss = {str(root):sum(rss for pid, (rss, _) in processes.items()
                                  if belongs(pid, {root})) for root in owned_pids}
        return dict(observable=True, available_gib=available, owned_rss_gib=owned_rss,
                    gpu_used_gib=used / 1024,
                    gpu_free_gib=(total - used) / 1024,
                    external_heavy=len(roots), external_peak_gib=external_peak,
                    external_pids=roots, idle_pids=list(idle_pids))
    except (OSError, ValueError, KeyError, subprocess.SubprocessError, RuntimeError) as exc:
        return dict(observable=False, error=str(exc), available_gib=0,
                    gpu_free_gib=0, external_heavy=0, external_peak_gib=0)


def stop_owned(process):
    """Give the adapter's TERM handler time to save, then reap its own group."""
    if process.poll() is not None:
        return
    os.killpg(process.pid, signal.SIGTERM)
    try:
        process.wait(timeout=15)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait(timeout=10)


def scope_failure_evidence(log_text):
    """Scope OOM may surface as SIGTERM after systemd kills the whole group."""
    match = re.search(r'Running as unit:\s*([\w.-]+\.scope)', log_text)
    if not match:
        return ''
    unit = match.group(1)
    try:
        result = subprocess.run(['journalctl', '--user', '-u', unit, '--no-pager', '-n', '30'],
                                capture_output=True, text=True, timeout=5)
        return result.stdout + result.stderr
    except (OSError, subprocess.SubprocessError) as exc:
        return 'scope journal unavailable: ' + str(exc)


def worker(job_path):
    from benchmark import run_one
    job = json.loads(Path(job_path).read_text())
    run_one(job['input_root'], job['output_root'], job['model'], job['fold'],
            job['seed'], command=job.get('command'))


def schedule(queue_path):
    # Admission and launch form one workstation-wide critical section. Separate
    # qualification queues must not both observe an empty host and race imports.
    lock_path = Path(__file__).resolve().parent.parent / '.scheduler.lock'
    with lock_path.open('a') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        return schedule_queue(queue_path)


def schedule_queue(queue_path):
    """Run a frozen queue, recording telemetry and durable attempt accounting."""
    from artifacts import run_directory, validate_run, run_identity
    queue_path = Path(queue_path).resolve()
    queue = json.loads(queue_path.read_text())
    output = Path(queue['output_root']).resolve()
    input_root = str(Path(queue['input_root']).resolve())
    jobs = queue['jobs']
    keys = [(job['model'], job['fold'], job['seed']) for job in jobs]
    if len(keys) != len(set(keys)):
        raise ValueError('duplicate run in scheduler queue')
    active = {}
    states = {}
    pending = []
    record = output / 'scheduling'
    record.mkdir(parents=True, exist_ok=True)
    def mirror_terminal(job, state):
        # Recovery/shutdown must publish the same durable outcome to the
        # prediction barrier as the scheduler, including pre-existing failures.
        if state.get('status') in TERMINAL and state['status'] != 'completed':
            directory = run_directory(output, job['model'], job['fold'], job['seed'])
            save(directory / 'run.json', dict(identity=state['identity'], status=state['status'],
                 stop_reason=state.get('stop_reason', 'recovered_terminal_state'),
                 training_seconds=state.get('training_seconds', 0)))
    for job in jobs:
        key = f"{job['model']}_fold{job['fold']}_seed{job['seed']}"
        job = dict(job, input_root=input_root, output_root=str(output), id=key)
        identity = run_identity(input_root, job['model'], job['fold'], job['seed'], job.get('command'))
        state_path = record / (key + '.json')
        state = json.loads(state_path.read_text()) if state_path.exists() else {}
        if state and state.get('identity') != identity:
            raise ValueError('stale scheduler identity: ' + key)
        state['identity'] = identity
        states[key] = state
        if state.get('status') == 'completed':
            validate_run(input_root, output, job['model'], job['fold'], job['seed'], job.get('command'))
        elif state.get('status') in TERMINAL:
            mirror_terminal(job, state)
            continue
        elif state.get('status') == 'running':
            # Reconnecting to a live orphan risks dispatching a duplicate fit.
            # A recovered dead attempt has already accrued durable sampled time.
            if Path(f"/proc/{state['pid']}").exists():
                raise ValueError('previous owned process still alive; refuse duplicate: ' + key)
            state = transition(state, 'interrupted', 0)
            states[key] = state
            save(state_path, state)
            mirror_terminal(job, state)
        else:
            pending.append(job)
    history = []
    blocked_since = None
    try:
        while pending or active:
            snapshot = process_snapshot([entry['process'].pid for entry in active.values()],
                                        queue.get('verified_idle_pids', []))
            snapshot['timestamp'] = time.time()
            with (record / 'telemetry.jsonl').open('a') as stream:
                stream.write(json.dumps(snapshot) + '\n')
            history.append(snapshot)
            history = history[-2:]
            victim = pressure(history, [dict(id=key, owned=True, started=entry['started'])
                                       for key, entry in active.items()])['stop']
            if victim:
                active[victim]['forced_reason'] = 'resource_pressure'
                stop_owned(active[victim]['process'])
            for key, entry in list(active.items()):
                state = states[key]
                state['peak_host_rss_gib'] = max(state.get('peak_host_rss_gib', 0),
                    snapshot.get('owned_rss_gib', {}).get(str(entry['process'].pid), 0))
                state['peak_workstation_gpu_used_gib'] = max(state.get('peak_workstation_gpu_used_gib', 0),
                    snapshot.get('gpu_used_gib', 0))
                now = time.monotonic()
                delta = now - entry['last_sample']
                entry['last_sample'] = now
                directory = run_directory(output, entry['job']['model'], entry['job']['fold'], entry['job']['seed'])
                progress_path = directory / 'progress.json'
                try:
                    progress = json.loads(progress_path.read_text())
                except (OSError, ValueError):
                    progress = {}
                if progress.get('phase') != 'inference':
                    state['training_seconds'] = state.get('training_seconds', 0) + delta
                else:
                    state['inference_seconds'] = state.get('inference_seconds', 0) + delta
                if state['training_seconds'] >= 14400:
                    entry['forced_reason'] = 'timeout'
                    stop_owned(entry['process'])
                elif state.get('inference_seconds', 0) >= entry['job'].get('inference_timeout_seconds', 1800):
                    entry['forced_reason'] = 'inference_timeout'
                    stop_owned(entry['process'])
                status = entry['process'].poll()
                if status is None:
                    save(record / (key + '.json'), state)
                    continue
                entry['log'].close()
                log_text = entry['log_path'].read_text(errors='replace')
                if status:
                    scope_evidence = scope_failure_evidence(log_text)
                    entry['log_path'].with_suffix('.journal.log').write_text(scope_evidence)
                    log_text += scope_evidence
                rss_matches = re.findall(r'Maximum resident set size \(kbytes\):\s*(\d+)', log_text)
                if rss_matches:
                    state['peak_host_rss_gib'] = max(state.get('peak_host_rss_gib', 0),
                                                   max(map(int, rss_matches)) / 1024 ** 2)
                if (directory / 'adapter.log').exists():
                    log_text += (directory / 'adapter.log').read_text(errors='replace')
                reason = entry.get('forced_reason')
                if reason is None:
                    if status == 0:
                        try:
                            validate_run(input_root, output, entry['job']['model'], entry['job']['fold'],
                                         entry['job']['seed'], entry['job'].get('command'))
                            reason = 'completed'
                        except ValueError:
                            reason = 'invalid_artifacts'
                    elif 'TENV_RETRY_ENVIRONMENT' in log_text:
                        reason = 'environment'
                    elif any(marker in log_text for marker in
                             ('CUDA out of memory', 'Memory cgroup out of memory',
                              'oom-kill', 'killed by the OOM killer')):
                        reason = 'oom'
                    else:
                        reason = 'process_exit_' + str(status)
                if reason == 'completed':
                    checkpoint = json.loads((directory / 'checkpoint.json').read_text())
                    if 'training_seconds' in checkpoint and 'inference_seconds' in checkpoint:
                        state['training_seconds'] = entry['training_base'] + checkpoint['training_seconds']
                        state['inference_seconds'] = entry['inference_base'] + checkpoint['inference_seconds']
                state = transition(state, reason, 0)
                state.setdefault('attempts', []).append(dict(returncode=status, reason=reason,
                    wall_seconds=now - entry['started'], log=str(entry['log_path'])))
                states[key] = state
                save(record / (key + '.json'), state)
                if state['status'] == 'pending':
                    pending.append(entry['job'])
                elif state['status'] != 'completed':
                    save(directory / 'run.json', dict(identity=state['identity'], status=state['status'],
                         stop_reason=reason, training_seconds=state['training_seconds']))
                del active[key]
            dispatched = False
            for job in list(pending):
                policy = admit(snapshot, job.get('profile'),
                               [entry['job'] for entry in active.values()])
                if not policy['admit']:
                    continue
                key = job['id']
                state = states[key]
                job_path = record / (key + '.job.json')
                save(job_path, job)
                attempt = len(state.get('attempts', []))
                log_path = record / f'{key}.attempt{attempt}.log'
                log = log_path.open('w')
                command = ['systemd-run', '--user', '--scope', '-p', 'MemoryMax=10G',
                           '-p', 'MemorySwapMax=8G', '--', '/usr/bin/time', '-v',
                           sys.executable, str(Path(__file__).resolve()), 'worker', str(job_path)]
                environment = dict(os.environ, OMP_NUM_THREADS='2', MKL_NUM_THREADS='2',
                                   OPENBLAS_NUM_THREADS='2', NUMEXPR_NUM_THREADS='2',
                                   TENV_BATCH_SCALE=str(state.get('batch_scale', 1.0)))
                process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                           start_new_session=True, env=environment)
                now = time.monotonic()
                state.update(status='running', pid=process.pid)
                state.setdefault('training_seconds', 0)
                save(record / (key + '.json'), state)
                active[key] = dict(process=process, job=job, started=now, last_sample=now,
                                   log=log, log_path=log_path,
                                   training_base=state['training_seconds'],
                                   inference_base=state.get('inference_seconds', 0))
                pending.remove(job)
                dispatched = True
                # Refresh real telemetry before any second dispatch.
                break
            if pending and not active and not dispatched:
                blocked_since = blocked_since or time.monotonic()
                if time.monotonic() - blocked_since >= queue.get('admission_timeout_seconds', 300):
                    save(record / 'blocked.json', dict(reason='resource admission', snapshot=snapshot,
                         pending=[job['id'] for job in pending]))
                    return 2
            else:
                blocked_since = None
            if pending or active:
                time.sleep(5)
    finally:
        for key, entry in active.items():
            stop_owned(entry['process'])
            entry['log'].close()
            states[key] = transition(states[key], 'scheduler_interrupted', 0)
            save(record / (key + '.json'), states[key])
            mirror_terminal(entry['job'], states[key])
    save(record / 'summary.json', {key: state for key, state in states.items()})
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    check = sub.add_parser('admit')
    check.add_argument('request', type=Path)
    advance = sub.add_parser('transition')
    advance.add_argument('state', type=Path)
    advance.add_argument('--reason', required=True)
    advance.add_argument('--seconds', type=float, required=True)
    press = sub.add_parser('pressure')
    press.add_argument('request', type=Path)
    sub.add_parser('telemetry')
    for name in ('run', 'worker'):
        action = sub.add_parser(name)
        action.add_argument('queue', type=Path)
    args = parser.parse_args()
    if args.action == 'admit':
        request = json.loads(args.request.read_text())
        print(json.dumps(admit(**request)))
    elif args.action == 'transition':
        state = json.loads(args.state.read_text()) if args.state.exists() else {}
        state = transition(state, args.reason, args.seconds)
        save(args.state, state)
        print(json.dumps(state))
    elif args.action == 'pressure':
        print(json.dumps(pressure(**json.loads(args.request.read_text()))))
    elif args.action == 'telemetry':
        print(json.dumps(process_snapshot()))
    elif args.action == 'run':
        sys.exit(schedule(args.queue))
    elif args.action == 'worker':
        worker(args.queue)


if __name__ == '__main__':
    main()
