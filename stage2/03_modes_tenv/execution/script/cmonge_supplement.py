"""Review-gated, one-shot CMonge v2 supplement for fold 0, seed 17."""
import argparse
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import stat
import subprocess
import sys
import threading
import time
import uuid

from artifacts import (TERMINAL, checkpoint_artifacts, command_evidence, read_json,
                       rows, run_directory, sha256, validate_prediction)
from benchmark import predictor_files
from cmonge_v2_protocol import EVAL_EVERY, PROTOCOL, Selection
from joint_v2_apply import mkdir_durable, sync_directory
from joint_v2_budget import matches, process_identity, write
from joint_v2_runtime import key, outcome, preflight
from schedule import admit, pressure, process_snapshot

SCRIPT = Path(__file__).resolve()
EXECUTION = SCRIPT.parent.parent
FORMAL = EXECUTION / 'formal'
ROOT = EXECUTION.parents[2]
LOCK = EXECUTION / '.scheduler.lock'
RECORD = EXECUTION / 'supplemental_record'
OUTPUT = EXECUTION.parent / 'output' / 'supplement-cmonge-v2-fold0-seed17'
ATTEMPT = 'cmonge-v2-fold0-seed17-supplement-1'
SLOT = ('cmonge', 0, 17)
KEY = 'cmonge_fold0_seed17'
MODELS = ('prnet', 'cpa', 'chemcpa', 'cmonge', 'perturbdiff', 'no_change', 'delta_mean')
TRAIN, INFERENCE = 14400, 1800
STAGE_SHA = '73fb21036a2d04a12c5d8d9ce712ed544b1c38eb70fbfff884dc6753c8b60028'
ORIGINAL_UNIT = 'arc-tenv-formal-scheduler.service'
ORIGINAL_INVOCATION = '6efea6fb957340f4a68754a09dff1940'
ORIGINAL_MAIN_PID = 1619811
SUPPLEMENT_UNIT = 'arc-tenv-cmonge-supplement.service'
SCOPE_PREFIX = 'arc-tenv-v2-cmonge-supplement-'
MEMORY_MAX, SWAP_MAX = 10 * 1024**3, 8 * 1024**3
CGROUP_ROOT = Path('/sys/fs/cgroup')
ENVIRONMENT = {
    'HOME': '/home/chang', 'USER': 'chang', 'LOGNAME': 'chang', 'LANG': 'C.UTF-8',
    'PATH': '/home/chang/miniconda3/envs/arc/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/lib/wsl/lib',
    **dict.fromkeys(('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'), '2'),
}
CONFIG = dict(protocol=PROTOCOL, iterations=10000, batch_size=512, eval_every=50, qualification=False)


def safe(path, *, directory=False):
    """Require a canonical path and a single-link regular file or real directory."""
    path = Path(path)
    if path.resolve() != path:
        raise ValueError(f'noncanonical path: {path}')
    info = path.lstat()
    if directory:
        if not stat.S_ISDIR(info.st_mode):
            raise ValueError(f'expected real directory: {path}')
    elif not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
        raise ValueError(f'expected single-link regular file: {path}')
    return path


def files_at(directory):
    return {p.name: sha256(safe(p)) for p in sorted(safe(directory, directory=True).iterdir())}


def read_pinned(path, digest):
    path = safe(Path(path).absolute())
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != digest:
        raise ValueError(f'reviewed digest mismatch: {path}')
    return json.loads(raw)


def frozen_job(formal, reviewed):
    if Path(formal).absolute() != FORMAL or reviewed != STAGE_SHA:
        raise ValueError('the committed joint v2 production generation is required')
    safe(OUTPUT.parent / 'formal' / 'joint-v2-apply.json')
    queue, expected = preflight(safe(FORMAL, directory=True), reviewed)
    slots = {(j['model'], j['fold'], j['seed']) for j in queue['jobs']}
    if slots != {(m, f, s) for m in MODELS for f in range(3) for s in (17, 42, 2026)}:
        raise ValueError('unexpected formal matrix')
    job = next(j for j in queue['jobs'] if key(j) == KEY)
    command = job['command']
    config = EXECUTION / 'models/cmonge/v2.config.json'
    if (job['inference_timeout_seconds'] != INFERENCE or
            command[1:] != [str(SCRIPT.with_name('cmonge_v2_gpu_adapter.py')), '--config', str(config)]
            or read_json(safe(config)) != CONFIG):
        raise ValueError('unexpected frozen CMonge recipe')
    return queue, expected, job


def original_records(queue, expected):
    output = safe(Path(queue['output_root']), directory=True)
    schedule = output / 'scheduling'
    summary_path = safe(schedule / 'summary.json')
    summary = read_json(summary_path)
    ledger_path = safe(FORMAL / 'budget-ledger.json')
    ledger = read_json(ledger_path)['ledger']
    if set(summary) != set(expected):
        raise ValueError('original scheduler summary is incomplete')
    for job in queue['jobs']:
        slot = key(job)
        directory = run_directory(output, job['model'], job['fold'], job['seed'])
        run = read_json(safe(directory / 'run.json'))
        record = read_json(safe(schedule / (slot + '.json')))
        status = run.get('status')
        if (status not in TERMINAL or run.get('identity') != expected[slot] or
                record.get('status') != status or record.get('identity') != expected[slot] or
                summary[slot] != record):
            raise ValueError(f'original slot has inconsistent terminal records: {slot}')
        if slot in ledger:
            row = ledger[slot]
            account = safe(Path(row['account_path']), directory=True)
            if account != output / 'joint_v2_budget' / slot or row['identity'] != expected[slot]:
                raise ValueError(f'original account path or identity mismatch: {slot}')
            for name in ('account.json', 'migration.json'):
                safe(account / name)
            state_status, state = outcome(account, ledger_path, job)
            if (state_status != status or record.get('attempts') != state['attempts']
                    or record.get('training_seconds') != state['consumed_seconds']):
                raise ValueError(f'original scheduler disagrees with budget account: {slot}')
            if slot == KEY:
                target = original_timeout(output, directory, job, run, record, row, account, state)
    return dict(target, summary_sha256=sha256(summary_path), ledger_sha256=sha256(ledger_path))


def original_timeout(output, directory, job, run, record, row, account, state):
    attempt, = state['attempts']
    prior = row['consumed_seconds']
    remaining = TRAIN - prior
    request_path = safe(account / 'request-0.json')
    request = read_json(request_path)
    original_request = read_json(safe(directory / 'request.json'))
    argv = attempt['command']
    if (run != record or run['status'] != 'timeout' or run['stop_reason'] != 'timeout' or
            row['cap_seconds'] != TRAIN or row['remaining_seconds'] != remaining or
            prior != 13975.509756658023 or remaining != 424.4902433419775 or
            attempt['prior_training_seconds'] != prior or attempt['reserved_seconds'] != remaining or
            attempt['wall_seconds'] < remaining or state['consumed_seconds'] != TRAIN or
            attempt['request'] != str(request_path) or
            argv[-len(job['command'])-2:] != [*job['command'], '--request', str(request_path)] or
            request != dict(original_request, prior_training_seconds=prior) or
            original_request['identity'] != row['identity'] or
            original_request['output_directory'] != str(directory)):
        raise ValueError('original target timeout does not match the reviewed exhausted attempt')
    return dict(status='timeout', stop_reason=run['stop_reason'], cap_seconds=TRAIN,
                prior_training_seconds=prior, remaining_seconds=remaining,
                worker_wall_seconds=attempt['wall_seconds'],
                consumed_seconds=state['consumed_seconds'], account=str(account),
                account_files=files_at(account), run_directory=str(directory),
                run_files=files_at(directory),
                scheduler_sha256=sha256(safe(output / 'scheduling' / (KEY + '.json'))))


def systemctl(*args, check=True):
    return subprocess.run(['/usr/bin/systemctl', '--user', *args], check=check,
                          capture_output=True, text=True, timeout=5)


def unit_properties(unit, *properties):
    report = systemctl('show', unit, '--property=' + ','.join(properties), check=False)
    values = dict(line.split('=', 1) for line in report.stdout.splitlines() if '=' in line)
    if not values or (report.returncode and values.get('LoadState') != 'not-found'):
        raise ValueError(f'cannot verify systemd unit {unit}')
    return values


def original_exit():
    if Path('/proc/1/comm').read_text().strip() not in ('systemd', 'init', 'wsl-init'):
        raise ValueError('actual host PID namespace required')
    service = unit_properties(ORIGINAL_UNIT, 'ActiveState', 'LoadState', 'MainPID',
                              'Result', 'ExecMainCode', 'ExecMainStatus', 'InvocationID')
    if service.get('ActiveState') != 'inactive' or service.get('MainPID', '0') != '0':
        raise ValueError('original scheduler has not cleanly stopped')
    if service.get('LoadState') != 'not-found':
        if (service.get('InvocationID'), service.get('Result'), service.get('ExecMainCode'),
                service.get('ExecMainStatus')) != (ORIGINAL_INVOCATION, 'success', '1', '0'):
            raise ValueError('original service invocation or exit status mismatch')
    # A collected transient unit has no ExecMainStatus. The original service's
    # GNU time main PID logged the scheduler's exit status to this invocation.
    result = subprocess.check_output([
        '/usr/bin/journalctl', '--user', '--no-pager', '-o', 'json',
        '_SYSTEMD_USER_UNIT=' + ORIGINAL_UNIT, '_SYSTEMD_INVOCATION_ID=' + ORIGINAL_INVOCATION,
        '_PID=' + str(ORIGINAL_MAIN_PID), '_UID=' + str(os.getuid()),
    ], text=True, timeout=10)
    exits = [entry for line in result.splitlines() if
             (entry := json.loads(line)).get('MESSAGE', '').strip().startswith('Exit status:')]
    if len(exits) != 1 or exits[0]['MESSAGE'].strip() != 'Exit status: 0':
        raise ValueError('original invocation has no unique authenticated zero exit')
    for line in systemctl('list-units', '--all', '--type=scope', '--plain', '--no-legend',
                          'arc-tenv-*').stdout.splitlines():
        if line and line.split()[2] not in ('inactive', 'failed'):
            raise ValueError(f'an owned TEN-V worker scope remains: {line}')
    # Catch detached watchdogs and workers even if their scope was collected.
    names = {b'joint_v2_runtime.py', b'joint_v2_budget.py', b'perturbdiff_v2_budget.py',
             b'cmonge_v2_gpu_adapter.py', b'cmonge_v2_adapter.py', b'perturbdiff_v2_adapter.py'}
    for item in Path('/proc').iterdir():
        if not item.name.isdecimal():
            continue
        try:
            args = (item / 'cmdline').read_bytes().split(b'\0')
            if any(arg.rsplit(b'/', 1)[-1] in names for arg in args):
                raise ValueError(f'owned scheduler or worker still visible: PID {item.name}')
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
    return dict(unit=ORIGINAL_UNIT, invocation=ORIGINAL_INVOCATION, main_pid=ORIGINAL_MAIN_PID,
                boot_id=exits[0]['_BOOT_ID'], exit_cursor=exits[0]['__CURSOR'], exit_status=0)


def policy_for(formal, reviewed, host_exit):
    queue, expected, job = frozen_job(formal, reviewed)
    output = Path(queue['output_root'])
    if output != OUTPUT.parent / 'formal':
        raise ValueError('unexpected original output root')
    original = original_records(queue, expected)
    journal = safe(output / 'joint-v2-apply.json')
    return dict(schema='tenv-cmonge-supplement-policy-v1', attempt_id=ATTEMPT,
                formal=str(FORMAL), reviewed_stage_sha256=reviewed,
                journal_sha256=sha256(journal),
                plan_sha256=sha256(safe(FORMAL / 'v2-migration-plan.json')),
                original=original, original_exit=host_exit, original_output=str(output),
                input_root=queue['input_root'], output_root=str(OUTPUT), account=str(OUTPUT / 'account'),
                terminal=str(RECORD / 'terminal.json'), job=job, identity=expected[KEY],
                caps=dict(training_seconds=TRAIN, inference_seconds=INFERENCE),
                worker_environment=ENVIRONMENT, runtime_sha256=sha256(SCRIPT))


def checked_policy(formal, reviewed, path, digest, *, live):
    policy = read_pinned(path, digest)
    exit_evidence = original_exit() if live else policy['original_exit']
    if policy != policy_for(formal, reviewed, exit_evidence):
        raise ValueError('policy differs from the fixed reviewed slot or preserved original evidence')
    return policy


def request_for(policy, *, verify_inputs):
    root = safe(Path(policy['input_root']), directory=True)
    manifest_path = safe(root / 'manifest.json')
    if sha256(manifest_path) != policy['identity']['input_sha256']:
        raise ValueError('predictor manifest hash mismatch')
    manifest = read_json(manifest_path)
    files = predictor_files(root, SLOT[1])
    for name in files.values():
        path = safe(Path(name))
        expected = manifest['files'][str(path.relative_to(root))]
        if verify_inputs and sha256(path) != expected:
            raise ValueError(f'predictor input hash mismatch: {path}')
    if len(rows(files['genes'])) != 18164 or manifest['gene_count'] != 18164:
        raise ValueError('the frozen gene axis is required')
    return dict(identity=policy['identity'], predictor_files=files, target_kind='delta',
                n_genes=18164, output_directory=str(run_directory(OUTPUT, *SLOT)),
                command_evidence=command_evidence(policy['job']['command']),
                inference_timeout_seconds=INFERENCE)


def lock_fd():
    fd = os.open(safe(LOCK), os.O_RDONLY | os.O_NOFOLLOW)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        check_lock(fd)
    except BaseException:
        os.close(fd)
        raise
    return fd


def check_lock(fd):
    with safe(LOCK).open('rb') as probe:
        if not os.path.samestat(os.fstat(fd), os.fstat(probe.fileno())):
            raise ValueError('scheduler lock inode changed or was not inherited')
        try:
            fcntl.flock(probe, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            try:
                fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise ValueError('inherited descriptor does not own the scheduler lock') from exc
        else:
            raise ValueError('scheduler lock is unheld')


def process_argv(pid):
    raw = Path(f'/proc/{pid}/cmdline').read_bytes()
    return list(map(os.fsdecode, raw.rstrip(b'\0').split(b'\0')))


def process_cgroup(pid):
    return Path(f'/proc/{pid}/cgroup').read_text().strip()


def run_command(path, digest):
    return [sys.executable, str(SCRIPT), 'run', '--formal', str(FORMAL),
            '--reviewed-stage-sha256', STAGE_SHA, '--policy', str(path),
            '--reviewed-policy-sha256', digest]


def validate_service(state):
    if (state.get('ActiveState') != 'active' or state.get('Restart') != 'no' or
            state.get('Type') != 'exec' or state.get('ExitType') != 'main' or
            not re.fullmatch(r'[0-9a-f]{32}', state.get('InvocationID', '')) or
            state.get('MemoryMax') != str(MEMORY_MAX) or state.get('MemorySwapMax') != str(SWAP_MAX) or
            state.get('WorkingDirectory') != str(ROOT) or
            not state.get('ControlGroup', '').endswith('/app.slice/' + SUPPLEMENT_UNIT)):
        raise ValueError('launch requires the reviewed dedicated systemd service and resource limits')


def check_service(supervisor, path, digest):
    state = unit_properties(SUPPLEMENT_UNIT, 'MainPID', 'InvocationID', 'ActiveState',
                            'Type', 'ExitType', 'Restart', 'MemoryMax', 'MemorySwapMax',
                            'ControlGroup', 'WorkingDirectory')
    validate_service(state)
    pid = supervisor['pid']
    main_pid = int(state['MainPID'])
    proc = Path(f'/proc/{pid}')
    parent = int(next(line for line in (proc / 'status').read_text().splitlines()
                      if line.startswith('PPid:')).split()[1])
    env = (proc / 'environ').read_bytes().split(b'\0')
    command = run_command(path, digest)
    if (not matches(supervisor) or parent != main_pid or
            process_cgroup(pid) != '0::' + state['ControlGroup'] or
            process_cgroup(main_pid) != '0::' + state['ControlGroup'] or
            process_argv(pid) != command or process_argv(main_pid) != ['/usr/bin/time', '-v', *command] or
            f'INVOCATION_ID={state["InvocationID"]}'.encode() not in env):
        raise ValueError('reserving process is not the reviewed supervisor invocation')
    return state


def admission(job):
    free = subprocess.check_output(['/usr/bin/free', '-g'], text=True, timeout=5)
    available = int(next(line for line in free.splitlines() if line.startswith('Mem:')).split()[-1])
    snapshot = process_snapshot()
    decision = admit(snapshot, job['profile'], (), max_active_jobs=1)
    if available < 12 or not decision['admit']:
        raise ValueError(f'host admission refused: {available} GiB available, {decision["reason"]}')
    return snapshot


def scope_command(scope, fd, executed):
    return ['/usr/bin/systemd-run', '--user', '--scope', '--slice=app.slice', '--unit=' + scope,
            '-p', f'MemoryMax={MEMORY_MAX}', '-p', f'MemorySwapMax={SWAP_MAX}',
            '-p', f'RuntimeMaxSec={TRAIN + INFERENCE}', '-p', 'KillSignal=SIGKILL',
            '-p', 'KillMode=control-group', '-p', 'TimeoutStopSec=1',
            '-p', 'BindsTo=' + SUPPLEMENT_UNIT, '-p', 'After=' + SUPPLEMENT_UNIT,
            '--', sys.executable, str(SCRIPT), 'worker', str(fd), scope,
            *executed]


def manager_usec(unit, interface, prop):
    encoded = re.sub(r'[^a-zA-Z0-9]', lambda match: f'_{ord(match[0]):02x}', unit)
    result = subprocess.check_output(
        ['/usr/bin/busctl', '--user', 'get-property', 'org.freedesktop.systemd1',
         '/org/freedesktop/systemd1/unit/' + encoded, 'org.freedesktop.systemd1.' + interface, prop],
        text=True, timeout=5)
    kind, value = result.split()
    if kind != 't':
        raise ValueError(f'expected integer microseconds for {unit} {prop}')
    return int(value)


def cutoff_record(scope, phase, deadline):
    unit = scope.removesuffix('.scope') + '-' + phase
    return dict(timer=unit + '.timer', service=unit + '.service', deadline_usec=round(deadline * 1e6))


def verify_cutoff(scope, record):
    timer, service, due = (record[k] for k in ('timer', 'service', 'deadline_usec'))
    props = unit_properties(timer, 'ActiveState', 'SubState', 'Unit', 'AccuracyUSec',
                            'WakeSystem', 'RandomizedDelayUSec', 'BindsTo', 'After')
    trigger = unit_properties(service, 'Conflicts', 'After', 'ExecStart', 'MemoryMax', 'MemorySwapMax')
    expected = dict(ActiveState='active', SubState='waiting', Unit=service,
                    AccuracyUSec='1us', WakeSystem='no', RandomizedDelayUSec='0')
    if (any(props.get(k) != v for k, v in expected.items()) or
            SUPPLEMENT_UNIT not in props.get('BindsTo', '').split() or
            SUPPLEMENT_UNIT not in props.get('After', '').split() or
            scope not in trigger.get('Conflicts', '').split() or scope not in trigger.get('After', '').split() or
            not trigger.get('ExecStart', '').startswith('{ path=/usr/bin/true ; argv[]=/usr/bin/true ;') or
            trigger.get('MemoryMax') != str(MEMORY_MAX) or trigger.get('MemorySwapMax') != str(SWAP_MAX) or
            manager_usec(timer, 'Timer', 'NextElapseUSecMonotonic') != due or
            time.monotonic_ns() // 1000 >= due):
        raise ValueError(f'manager cutoff could not be verified: {timer}')


def arm_cutoff(scope, record):
    # The user manager interprets OnBootSec relative to its own startup.
    offset = int(systemctl('show', '--value', '-p', 'UserspaceTimestampMonotonic').stdout)
    due = record['deadline_usec'] - offset
    subprocess.run(
        ['/usr/bin/systemd-run', '--user', '--slice=app.slice', '--unit=' + record['service'],
         f'--on-boot={due}us', '--timer-property=AccuracySec=1us',
         '--timer-property=RemainAfterElapse=no', f'--timer-property=BindsTo={SUPPLEMENT_UNIT}',
         f'--timer-property=After={SUPPLEMENT_UNIT}', '-p', f'Conflicts={scope}', '-p', f'After={scope}',
         '-p', f'MemoryMax={MEMORY_MAX}', '-p', f'MemorySwapMax={SWAP_MAX}', '/usr/bin/true'],
        check=True, capture_output=True, timeout=5)
    verify_cutoff(scope, record)


def stop_cutoff(unit):
    systemctl('stop', unit, check=False)
    if unit_properties(unit, 'ActiveState').get('ActiveState') != 'inactive':
        raise ValueError(f'manager cutoff still active: {unit}')


def scope_proof(scope, service):
    props = unit_properties(scope, 'ControlGroup', 'ActiveState', 'InvocationID', 'BindsTo',
                            'After', 'KillMode', 'KillSignal', 'MemoryMax', 'MemorySwapMax')
    expected = str(Path(service['ControlGroup']).parent / scope)
    cgroup = CGROUP_ROOT / expected.lstrip('/')
    if (props.get('ControlGroup') != expected or process_cgroup('self') != '0::' + expected or
            props.get('ActiveState') != 'active' or
            not re.fullmatch(r'[0-9a-f]{32}', props.get('InvocationID', '')) or
            SUPPLEMENT_UNIT not in props.get('BindsTo', '').split() or
            SUPPLEMENT_UNIT not in props.get('After', '').split() or
            props.get('KillMode') != 'control-group' or props.get('KillSignal') != str(signal.SIGKILL) or
            props.get('MemoryMax') != str(MEMORY_MAX) or props.get('MemorySwapMax') != str(SWAP_MAX) or
            (cgroup / 'memory.max').read_text().strip() != str(MEMORY_MAX) or
            (cgroup / 'memory.swap.max').read_text().strip() != str(SWAP_MAX) or Path.cwd() != ROOT):
        raise ValueError('worker scope or working-directory limits differ from the reviewed recipe')
    runtime = manager_usec(scope, 'Scope', 'RuntimeMaxUSec')
    stop = manager_usec(scope, 'Scope', 'TimeoutStopUSec')
    if runtime != (TRAIN + INFERENCE) * 10**6 or stop != 10**6:
        raise ValueError('worker scope does not enforce hard stops')
    return scope_record(expected, props['InvocationID'])


def scope_record(cgroup, invocation):
    return dict(control_group=cgroup, invocation=invocation, runtime_usec=(TRAIN + INFERENCE) * 10**6,
                timeout_stop_usec=10**6, kill_signal=signal.SIGKILL, kill_mode='control-group',
                binds_to=SUPPLEMENT_UNIT, after=SUPPLEMENT_UNIT, memory_max=MEMORY_MAX, swap_max=SWAP_MAX)


def lock_records(pid, fd):
    return [line.split(': ', 2)[-1] for line in
            Path(f'/proc/{pid}/fdinfo/{fd}').read_text().splitlines() if line.startswith('lock:')]


def worker_exec(fd, scope, command):
    """Authenticate the live reservation before entering the fixed native command."""
    if fd < 3 or not re.fullmatch(SCOPE_PREFIX + r'[0-9a-f]{32}\.scope', scope):
        raise ValueError('invalid scope or lock descriptor')
    check_lock(fd)
    account = OUTPUT / 'account'
    claim = safe(RECORD / 'launch-claim.json')
    claim_digest = sha256(claim)
    reservation = read_pinned(claim, claim_digest)
    if (reservation != read_json(safe(account / 'account.json')) or reservation['status'] != 'reserved'
            or reservation['scope'] != scope or reservation['lock_fd'] != fd):
        raise ValueError('scope launch is not the durable reserved worker')
    service = check_service(reservation['supervisor'], reservation['policy_path'], reservation['policy_sha256'])
    if service != reservation['service']:
        raise ValueError('reserving service invocation changed')
    own_locks = lock_records('self', fd)
    if len(own_locks) != 1 or own_locks != lock_records(reservation['supervisor']['pid'], fd):
        raise ValueError('lock does not match the reserving supervisor')
    scope_state = scope_proof(scope, service)
    policy = checked_policy(FORMAL, STAGE_SHA, reservation['policy_path'],
                            reservation['policy_sha256'], live=False)
    request = request_for(policy, verify_inputs=False)
    validate_reservation(policy, reservation['policy_sha256'], reservation, request)
    if command != reservation['executed_argv']:
        raise ValueError('worker command differs from the reviewed execution envelope')
    training = reservation['training_cutoff']
    verify_cutoff(scope, training)
    if sha256(claim) != claim_digest or os.path.lexists(account / 'worker-start.json'):
        raise ValueError('worker launch claim changed or was already used')
    guard_pid = os.fork()
    if guard_pid == 0:
        try:
            guard_phase(reservation, Path(request['output_directory']), account)
        except Exception as exc:
            try:
                write(account / 'phase-guard-error.json', dict(error=str(exc)))
            finally:
                systemctl('kill', '--kill-whom=all', '--signal=SIGKILL', scope, check=False)
        finally:
            os._exit(0)
    write(account / 'worker-start.json',
          dict(process=process_identity(os.getpid()), guard=process_identity(guard_pid),
               scope=scope, scope_state=scope_state,
               supervisor=reservation['supervisor'], invocation=service['InvocationID'],
               claim_sha256=claim_digest, policy_sha256=reservation['policy_sha256'],
               training_cutoff=training, injected_request_sha256=reservation['injected_request_sha256']))
    os.set_inheritable(fd, True)
    os.execve(command[0], command, ENVIRONMENT)


def cleanup(worker, scope, identity):
    # Killing the local launch session prevents a late scope registration.
    # The systemd stop then covers any setsid descendant of the frozen command.
    if worker.poll() is None:
        if not matches(identity):
            raise ValueError('worker process identity changed')
        os.killpg(worker.pid, signal.SIGKILL)
    worker.wait(timeout=10)
    systemctl('kill', '--kill-whom=all', '--signal=SIGKILL', scope, check=False)
    stopped = systemctl('stop', scope, check=False)
    state = unit_properties(scope, 'ActiveState', 'LoadState', 'ControlGroup')
    cgroup = state.get('ControlGroup')
    if (state.get('ActiveState') not in ('inactive', 'failed') or
            (stopped.returncode and state.get('LoadState') != 'not-found') or
            (cgroup and (Path('/sys/fs/cgroup') / cgroup.lstrip('/') / 'cgroup.events').exists()
             and 'populated 0' not in (Path('/sys/fs/cgroup') / cgroup.lstrip('/') / 'cgroup.events').read_text().splitlines())):
        raise ValueError('owned scope exit unverified; one-shot reservation remains spent')
    return dict(scope=scope, active_state=state['ActiveState'], load_state=state.get('LoadState'))


class PhaseDeadline:
    def __init__(self, start):
        self.start = start
        self.transition = None
        self.training_seconds = None

    @property
    def deadline(self):
        if self.transition is None:
            return self.start + TRAIN
        return min(self.start + TRAIN + INFERENCE, self.transition + INFERENCE)

    def observe(self, now, progress):
        # Refuse late transitions before processing any phase information.
        if now >= self.deadline:
            return 'training_timeout' if self.transition is None else 'inference_timeout'
        if progress is None:
            if self.training_seconds is not None:
                raise ValueError('adapter progress disappeared')
            return None
        phase = progress['phase']
        seconds = progress['training_seconds']
        if (type(seconds) not in (int, float) or not math.isfinite(seconds) or
                not 0 <= seconds <= min(TRAIN, now - self.start) or
                (self.training_seconds is not None and seconds < self.training_seconds) or
                phase not in ('training', 'inference')):
            raise ValueError('invalid or regressed native progress')
        if self.transition is not None:
            if phase != 'inference' or seconds != self.training_seconds:
                raise ValueError('native inference transition changed')
        elif phase == 'inference':
            self.transition = now
        self.training_seconds = seconds
        return None

    def timing(self, end):
        split = end if self.transition is None else self.transition
        return dict(worker_wall_seconds=end - self.start,
                    pre_inference_watchdog_seconds=split - self.start,
                    observed_inference_seconds=None if self.transition is None else end - split,
                    adapter_training_seconds=self.training_seconds)


def read_progress(path):
    try:
        return read_json(safe(path))
    except FileNotFoundError:
        return None


def guard_phase(state, directory, account):
    """Transfer the manager cutoff independently of the reserving supervisor."""
    scope = state['scope']
    clock = PhaseDeadline(state['started_usec'] / 1e6)
    while clock.transition is None:
        progress = read_progress(directory / 'progress.json')
        now = time.monotonic()
        if clock.observe(now, progress):
            return
        if clock.transition is not None:
            cutoff = cutoff_record(scope, 'inference', clock.deadline)
            arm_cutoff(scope, cutoff)
            if time.monotonic() < clock.start + TRAIN:
                stop_cutoff(state['training_cutoff']['timer'])
            if time.monotonic() >= clock.start + TRAIN:
                systemctl('kill', '--kill-whom=all', '--signal=SIGKILL', scope, check=False)
                return
            transition = dict(clock.timing(now), manager_cutoff=cutoff)
            write(account / 'inference-transition.json', transition)
            return
        time.sleep(min(.05, max(0, clock.deadline - now)))


def watch(worker, clock, progress_path, stop, transition_path, scope=None):
    pending = None
    while True:
        code = worker.poll()
        progress = read_progress(progress_path)
        now = time.monotonic()
        old = clock.transition
        if scope is not None and old is None:
            transition = read_progress(transition_path)
            if transition is not None:
                observed = clock.start + transition['pre_inference_watchdog_seconds']
                if (progress is None or progress['phase'] != 'inference' or observed > now or
                        clock.observe(observed, progress) or
                        transition != dict(clock.timing(observed),
                                           manager_cutoff=cutoff_record(scope, 'inference', clock.deadline))):
                    raise ValueError('phase guardian recorded an invalid inference transition')
            elif progress is not None and progress['phase'] == 'inference' and now < clock.deadline:
                pending = now if pending is None else pending
                if now - pending > 12:
                    raise ValueError('phase guardian did not authenticate the inference transition')
                time.sleep(min(.05, clock.deadline - now))
                continue
        reason = clock.observe(now, progress)
        if scope is None and old is None and clock.transition is not None:
            write(transition_path, clock.timing(now))
        end = time.monotonic()
        if reason or end >= clock.deadline:
            return 'timeout', reason or ('training_timeout' if clock.transition is None else 'inference_timeout'), end, code
        if code == 124 and clock.transition is None:
            return 'timeout', 'native_training_timeout', end, code
        if code is not None:
            if code != 0:
                return 'failed', 'worker_nonzero_exit', end, code
            if clock.transition is None:
                return 'failed', 'worker_exit_without_inference', end, code
            return 'completed', 'worker_zero_exit', end, code
        if stop.is_set():
            return 'failed', 'host_resource_pressure', end, code
        time.sleep(min(.05, max(0, clock.deadline - end)))


def owned_rss(owned):
    peak = 0
    live = False
    for pid in owned:
        try:
            for line in Path(f'/proc/{pid}/status').read_text().splitlines():
                if line.startswith('VmRSS:'):
                    peak += int(line.split()[1])
                    live = True
                    break
        except FileNotFoundError:
            continue
    return peak / 1024**2 if live else None


def sample_resources(done, stop, path, pid, scope, scope_cgroup):
    history = []
    try:
        with path.open('x', buffering=1) as stream:
            while True:
                # A scope can reparent its command away from the systemd-run client.
                process_file = scope_cgroup / 'cgroup.procs'
                try:
                    owned = {pid, *map(int, process_file.read_text().split())}
                except FileNotFoundError:
                    owned = {pid}
                snapshot = process_snapshot(owned_pids=sorted(owned))
                stream.write(json.dumps(dict(snapshot, scope=scope, timestamp=time.time(),
                                             sampled_owned_rss_gib=owned_rss(owned))) + '\n')
                history = [*history[-1:], snapshot]
                if (not snapshot['observable'] or snapshot['external_heavy'] or
                        pressure(history, [dict(owned=True, started=0, id=scope)])['stop']):
                    stop.set()
                if done.wait(5):
                    break
            os.fsync(stream.fileno())
    except (OSError, ValueError):
        stop.set()


def validate_completion(directory, request, state):
    checkpoint = read_json(safe(directory / 'checkpoint.json'))
    history = read_json(safe(directory / 'history.json'))
    if set(checkpoint['blobs']) != {'best.msgpack', 'latest.msgpack', 'history.json'}:
        raise ValueError('unexpected native checkpoint blobs')
    for name in checkpoint['blobs']:
        safe(directory / name)
    checkpoint_artifacts(directory, checkpoint, require_blobs=True)
    timing = state['timing']
    progress = read_json(safe(directory / 'progress.json'))
    inference = checkpoint['inference_seconds']
    fields = dict(identity=request['identity'], protocol=PROTOCOL, prior_training_seconds=0,
                  gene_count=18164, backend='gpu', validation_groups=138, batch_size=512,
                  maximum_iterations=10000, qualification_only=False,
                  snapshot_semantics='prediction-only-not-resumable',
                  early_stopping=dict(minimum_iterations=500, patience=10, relative_min_delta=.001,
                                      earliest_stop_iteration=1000))
    if (any(checkpoint.get(k) != v for k, v in fields.items()) or
            history != {k: checkpoint[k] for k in ('protocol', 'identity', 'history', 'snapshot_semantics')} or
            progress != dict(phase='inference', training_seconds=checkpoint['training_seconds']) or
            timing['adapter_training_seconds'] != checkpoint['training_seconds'] or
            type(inference) not in (int, float) or not math.isfinite(inference) or
            not 0 <= inference <= INFERENCE or
            inference > timing['worker_wall_seconds'] - checkpoint['training_seconds']):
        raise ValueError('native completion evidence or timing mismatch')
    events = checkpoint['history']
    iteration = checkpoint['iterations']
    if (type(iteration) is not int or not 1000 <= iteration <= 10000 or
            [event['iteration'] for event in events] != list(range(EVAL_EVERY, iteration + 1, EVAL_EVERY)) or
            events[-1]['iteration'] != iteration):
        raise ValueError('native validation cadence mismatch')
    selection = Selection()
    for event in events:
        if event['validation_conditions_evaluated'] != 138:
            raise ValueError('native validation excluded conditions')
        _, stop = selection.observe(event['iteration'], event['validation_loss'])
        if stop and event is not events[-1]:
            raise ValueError('native validation continued past early stopping')
    if (checkpoint['best_iteration'] != selection.best_iteration or
            checkpoint['best_validation_loss'] != selection.best_loss or
            checkpoint['stop_reason'] != ('validation_patience' if stop else 'maximum_iterations') or
            (not stop and iteration != 10000)):
        raise ValueError('native checkpoint selection evidence mismatch')
    artifacts = {name: sha256(safe(directory / name)) for name in
                 ('prediction.npy', 'case_ids.json', 'checkpoint.json', 'request.json')}
    prediction = validate_prediction(directory, rows(request['predictor_files']['test_cases']), 18164)
    if prediction.dtype.name != 'float32':
        raise ValueError('unexpected prediction dtype')
    return artifacts, inference


def untouched():
    for path in (OUTPUT, RECORD / 'launch-claim.json', RECORD / 'terminal.json'):
        if os.path.lexists(path):
            raise ValueError(f'supplement was already claimed: {path}')
    safe(OUTPUT.parent, directory=True)
    if os.path.lexists(RECORD):
        safe(RECORD, directory=True)


def reserve(policy, path, digest, request, fd, service, snapshot):
    untouched()
    OUTPUT.mkdir(mode=0o700)
    sync_directory(OUTPUT.parent)
    mkdir_durable(RECORD)
    account = OUTPUT / 'account'
    account.mkdir(mode=0o700)
    directory = Path(request['output_directory'])
    mkdir_durable(directory)
    write(directory / 'request.json', request)
    injected = account / 'request-0.json'
    write(injected, dict(request, prior_training_seconds=0))
    write(account / 'policy.json', policy)
    scope = SCOPE_PREFIX + uuid.uuid4().hex + '.scope'
    executed = executed_command(policy, injected)
    argv = scope_command(scope, fd, executed)
    start = time.monotonic_ns() // 1000
    state = dict(schema='tenv-cmonge-supplement-account-v1', status='reserved',
                 attempt_id=ATTEMPT, policy_path=str(path), policy_sha256=digest, identity=policy['identity'],
                 caps=policy['caps'], prior_training_seconds=0, credit_available_seconds=0,
                 reservation_count=1, supervisor=process_identity(os.getpid()), scope=scope,
                 lock_fd=fd, started_usec=start,
                 training_cutoff=cutoff_record(scope, 'training', start / 1e6 + TRAIN),
                 scope_command=argv, executed_argv=executed,
                 request_sha256=sha256(directory / 'request.json'),
                 injected_request_sha256=sha256(injected), service=service, admission=snapshot)
    # The protected claim and distinct account are durable before worker creation.
    write(RECORD / 'launch-claim.json', state)
    write(account / 'account.json', state)
    return state


def executed_command(policy, injected):
    return ['/usr/bin/time', '-v', *policy['job']['command'], '--request', str(injected)]


def validate_reservation(policy, digest, claim, request):
    account = OUTPUT / 'account'
    original = safe(Path(request['output_directory']) / 'request.json')
    injected = safe(account / 'request-0.json')
    scope = claim['scope']
    command = executed_command(policy, injected)
    validate_service(claim['service'])
    if (claim['status'] != 'reserved' or claim['schema'] != 'tenv-cmonge-supplement-account-v1' or
            claim['attempt_id'] != ATTEMPT or claim['policy_sha256'] != digest or
            claim['identity'] != policy['identity'] or claim['caps'] != policy['caps'] or
            claim['reservation_count'] != 1 or claim['credit_available_seconds'] != 0 or
            claim['prior_training_seconds'] != 0 or not claim['admission']['observable'] or
            not re.fullmatch(SCOPE_PREFIX + r'[0-9a-f]{32}\.scope', scope) or
            type(claim['lock_fd']) is not int or claim['lock_fd'] < 3 or
            type(claim['started_usec']) is not int or claim['started_usec'] <= 0 or
            claim['training_cutoff'] != cutoff_record(scope, 'training', claim['started_usec'] / 1e6 + TRAIN) or
            claim['executed_argv'] != command or claim['scope_command'] != scope_command(scope, claim['lock_fd'], command) or
            read_json(safe(account / 'policy.json')) != policy or
            read_json(original) != request or read_json(injected) != dict(request, prior_training_seconds=0) or
            sha256(original) != claim['request_sha256'] or sha256(injected) != claim['injected_request_sha256']):
        raise ValueError('invalid one-shot reservation or execution envelope')


def supervise(state, directory, account, fd):
    scope = state['scope']
    worker = None
    identity = None
    done, stop = threading.Event(), threading.Event()
    sampler = None
    scope_cgroup = (CGROUP_ROOT / state['service']['ControlGroup'].lstrip('/')).parent / scope
    clock = PhaseDeadline(state['started_usec'] / 1e6)
    with (account / 'watchdog.log').open('xb', buffering=0) as log:
        try:
            arm_cutoff(scope, state['training_cutoff'])
            worker = subprocess.Popen(state['scope_command'], pass_fds=(fd,), start_new_session=True,
                                      stdin=subprocess.DEVNULL, stdout=log, stderr=log,
                                      env=dict(ENVIRONMENT, XDG_RUNTIME_DIR=f'/run/user/{os.getuid()}'))
            try:
                identity = process_identity(worker.pid)
            except FileNotFoundError:
                if worker.poll() is None:
                    raise
            sampler = threading.Thread(target=sample_resources,
                args=(done, stop, account / 'telemetry.jsonl', worker.pid, scope, scope_cgroup), daemon=True)
            sampler.start()
            status, reason, end, code = watch(worker, clock, directory / 'progress.json',
                                              stop, account / 'inference-transition.json', scope)
        finally:
            done.set()
            if worker is not None:
                cleanup_proof = cleanup(worker, scope, identity)
            for phase in ('training', 'inference'):
                cutoff = cutoff_record(scope, phase, 0)
                stop_cutoff(cutoff['timer'])
                stop_cutoff(cutoff['service'])
            if sampler is not None:
                sampler.join(timeout=12)
                if sampler.is_alive():
                    raise ValueError('resource telemetry did not settle')
        timing = clock.timing(end)
        log.flush()
        os.fsync(log.fileno())
    # Failure of the isolation handshake is a failed launch, never a settled
    # scientific outcome eligible for substitution.
    started = read_json(safe(account / 'worker-start.json'))
    if matches(started['process']):
        raise ValueError('owned scoped worker start or exit not authenticated')
    resources = resource_evidence(account, scope, started['scope_state']['invocation'])
    return dict(status=status, stop_reason=reason, timing=timing, returncode=worker.returncode,
                observed_returncode=code, process=identity, worker=started,
                cleanup=cleanup_proof, **resources)


def publish_terminal(policy, digest, state, request):
    account = Path(policy['account'])
    directory = Path(request['output_directory'])
    artifacts, inference = None, None
    if state['status'] == 'completed':
        artifacts, inference = validate_completion(directory, request, state)
    state['timing']['adapter_inference_seconds'] = inference
    run = dict(status=state['status'], identity=policy['identity'], attempt_id=ATTEMPT,
               stop_reason=state['stop_reason'], training_seconds=state['timing']['adapter_training_seconds'],
               inference_seconds=inference, policy_sha256=digest)
    if artifacts is not None:
        run['artifacts'] = artifacts
    validate_execution(policy, digest, state, request, run)
    write(directory / 'run.json', run)
    write(account / 'account.json', state)
    terminal = dict(schema='tenv-cmonge-supplement-terminal-v1', status=run['status'],
                    attempt_id=ATTEMPT, identity=policy['identity'], run_directory=str(directory),
                    policy_sha256=digest, original=policy['original'], account=state,
                    account_files=files_at(account), run_sha256=sha256(directory / 'run.json'),
                    claim_sha256=sha256(safe(RECORD / 'launch-claim.json')), files=files_at(directory))
    write(RECORD / 'terminal.json', terminal)
    return terminal


def seconds(value):
    return type(value) in (int, float) and math.isfinite(value) and value >= 0


def resource_evidence(account, scope, invocation):
    log = safe(account / 'watchdog.log').read_text(errors='replace')
    if f'Running as unit: {scope}; invocation ID: {invocation}' not in log.splitlines():
        raise ValueError('watchdog log lacks the owned scope invocation')
    peaks = [int(value) / 1024**2 for value in
             re.findall(r'Maximum resident set size \(kbytes\):\s*(\d+)', log)]
    telemetry = [json.loads(line) for line in safe(account / 'telemetry.jsonl').read_text().splitlines()]
    if not telemetry:
        raise ValueError('missing worker resource telemetry')
    rss, gpu = [], []
    for sample in telemetry:
        owned = sample['sampled_owned_rss_gib']
        if (sample['scope'] != scope or not seconds(sample['timestamp']) or
                type(sample['observable']) is not bool or
                (owned is not None and not seconds(owned))):
            raise ValueError('invalid resource sample for owned scope')
        if owned is not None:
            rss.append(owned)
        if sample['observable']:
            if (not all(seconds(sample[k]) for k in ('gpu_used_gib', 'gpu_free_gib', 'available_gib')) or
                    type(sample['external_heavy']) is not int or sample['external_heavy'] < 0):
                raise ValueError('invalid observable host or GPU resource sample')
            gpu.append(sample['gpu_used_gib'])
    pressured = any(not sample['observable'] or sample['external_heavy'] for sample in telemetry)
    pressured |= any(pressure(telemetry[i-1:i+1], [dict(owned=True, started=0, id=scope)])['stop']
                     is not None for i in range(1, len(telemetry)))
    return dict(peak_host_rss_gib=max(rss + peaks, default=None),
                sampled_peak_host_rss_gib=max(rss, default=None),
                time_peak_rss_gib=max(peaks, default=None),
                peak_workstation_gpu_used_gib=max(gpu, default=None), observed_resource_pressure=pressured)


def validate_terminal_timing(state, run):
    timing = state['timing']
    wall = timing['worker_wall_seconds']
    pre = timing['pre_inference_watchdog_seconds']
    observed = timing['observed_inference_seconds']
    native_training = timing['adapter_training_seconds']
    native_inference = timing['adapter_inference_seconds']
    if (not seconds(wall) or not seconds(pre) or pre > wall or
            (observed is None and pre != wall) or
            (observed is not None and (not seconds(observed) or observed != wall - pre
                                       and not math.isclose(observed, wall - pre, abs_tol=1e-9) or pre >= TRAIN)) or
            (native_training is not None and
             (not seconds(native_training) or native_training > min(TRAIN, pre))) or
            run['training_seconds'] != native_training or run['inference_seconds'] != native_inference):
        raise ValueError('inconsistent supplemental watchdog or native timing')
    phase_elapsed = wall >= TRAIN if observed is None else observed >= INFERENCE or wall >= TRAIN + INFERENCE
    code = state['observed_returncode']
    rc = state['returncode']
    if type(rc) is not int or (code is not None and (type(code) is not int or code != rc)):
        raise ValueError('supplemental exit status is inconsistent')
    reason = state['stop_reason']
    status = state['status']
    if status == 'completed':
        if (observed is None or observed >= INFERENCE or wall >= TRAIN + INFERENCE or
                code != 0 or reason != 'worker_zero_exit' or
                not seconds(native_inference) or native_inference > INFERENCE):
            raise ValueError('completed attempt has inconsistent phase timing or exit status')
        return
    elif native_inference is not None:
        raise ValueError('unsuccessful attempt reported completed inference')
    if status == 'timeout':
        valid = ((reason == ('training_timeout' if observed is None else 'inference_timeout') and phase_elapsed)
                 or (reason == 'native_training_timeout' and not phase_elapsed and observed is None and code == 124))
    elif status == 'failed':
        valid = not phase_elapsed and (
            (reason == 'worker_nonzero_exit' and code not in (None, 0) and (code != 124 or observed is not None)) or
            (reason == 'host_resource_pressure' and code is None and rc != 0 and state['observed_resource_pressure']))
    else:
        valid = False
    if not valid:
        raise ValueError('ineligible supplemental status, reason, return code or deadline')


def validate_execution(policy, digest, state, request, run):
    claim = safe(RECORD / 'launch-claim.json')
    reservation = read_json(claim)
    validate_reservation(policy, digest, reservation, request)
    if any(state[k] != v for k, v in reservation.items() if k != 'status'):
        raise ValueError('terminal account diverges from the reserved execution')
    scope = state['scope']
    account = OUTPUT / 'account'
    started = read_json(safe(account / 'worker-start.json'))
    scope_state = started['scope_state']
    invocation = scope_state['invocation']
    process = started['process']
    guard = started['guard']
    expected_cgroup = str(Path(reservation['service']['ControlGroup']).parent / scope)
    if (os.path.lexists(account / 'phase-guard-error.json') or started != state['worker'] or
            started != dict(process=process, guard=guard, scope=scope, scope_state=scope_state,
                            supervisor=state['supervisor'], invocation=state['service']['InvocationID'],
                            claim_sha256=sha256(claim), policy_sha256=digest,
                            training_cutoff=reservation['training_cutoff'],
                            injected_request_sha256=reservation['injected_request_sha256']) or
            process != state['process'] or process['pid'] == state['supervisor']['pid'] or
            type(process['pid']) is not int or process['pid'] <= 1 or
            process['boot_id'] != state['supervisor']['boot_id'] or
            process['start_ticks'] < state['supervisor']['start_ticks'] or
            guard['pid'] in (process['pid'], state['supervisor']['pid']) or
            type(guard['pid']) is not int or guard['pid'] <= 1 or
            guard['boot_id'] != process['boot_id'] or guard['start_ticks'] < process['start_ticks'] or
            not re.fullmatch(r'[0-9a-f]{32}', invocation) or
            scope_state != scope_record(expected_cgroup, invocation) or
            state['cleanup']['scope'] != scope or state['cleanup']['active_state'] not in ('inactive', 'failed')):
        raise ValueError('terminal worker identity, scope or cleanup is inconsistent')
    resources = resource_evidence(account, scope, invocation)
    if (any(state[k] != v for k, v in resources.items()) or
            (state['status'] == 'completed' and
             (resources['time_peak_rss_gib'] is None or resources['peak_workstation_gpu_used_gib'] is None))):
        raise ValueError('terminal peak resources differ from the captured evidence')
    validate_terminal_timing(state, run)
    timing = state['timing']
    observed = timing['observed_inference_seconds']
    transition = account / 'inference-transition.json'
    progress = read_progress(Path(request['output_directory']) / 'progress.json')
    native = timing['adapter_training_seconds']
    if observed is None:
        if os.path.lexists(transition):
            raise ValueError('training outcome unexpectedly includes a committed inference transition')
        if ((progress is None and native is not None) or
                (progress is not None and
                 (not seconds(progress['training_seconds']) or
                  (native is not None and progress['training_seconds'] < native) or
                  (progress['phase'] != 'training' and
                   not (progress['phase'] == 'inference' and run['stop_reason'] == 'training_timeout'))))):
            raise ValueError('terminal training evidence disagrees with native progress')
    else:
        pre = timing['pre_inference_watchdog_seconds']
        start = state['started_usec'] / 1e6
        if (native is None or progress != dict(phase='inference', training_seconds=native) or
                read_json(safe(transition)) != dict(worker_wall_seconds=pre,
                    pre_inference_watchdog_seconds=pre, observed_inference_seconds=0,
                    adapter_training_seconds=timing['adapter_training_seconds'],
                    manager_cutoff=cutoff_record(scope, 'inference', start + pre + INFERENCE))):
            raise ValueError('terminal inference evidence differs from the committed transition')


def verify_terminal(formal, reviewed_stage_sha256, policy_path, reviewed_policy_sha256,
                    reviewed_terminal_sha256):
    """Authenticate the reviewed terminal outcome and both isolated histories."""
    policy = checked_policy(formal, reviewed_stage_sha256, policy_path, reviewed_policy_sha256, live=False)
    terminal = read_pinned(Path(policy['terminal']), reviewed_terminal_sha256)
    directory = run_directory(OUTPUT, *SLOT)
    account = OUTPUT / 'account'
    state = read_json(safe(account / 'account.json'))
    claim = safe(RECORD / 'launch-claim.json')
    request = request_for(policy, verify_inputs=False)
    run = read_json(safe(directory / 'run.json'))
    status = state['status']
    if (terminal != dict(schema='tenv-cmonge-supplement-terminal-v1', status=status,
                        attempt_id=ATTEMPT, identity=policy['identity'], run_directory=str(directory),
                        policy_sha256=reviewed_policy_sha256, original=policy['original'], account=state,
                        account_files=files_at(account), run_sha256=sha256(directory / 'run.json'),
                        claim_sha256=sha256(claim), files=files_at(directory)) or
            run['status'] != status or run['identity'] != policy['identity'] or
            run['attempt_id'] != ATTEMPT or run['policy_sha256'] != reviewed_policy_sha256 or
            run['stop_reason'] != state['stop_reason']):
        raise ValueError('terminal evidence does not match the one-shot supplemental account')
    validate_execution(policy, reviewed_policy_sha256, state, request, run)
    if status == 'completed':
        artifacts, inference = validate_completion(directory, request, state)
        if (artifacts != run['artifacts'] or state['returncode'] != 0 or
                run['inference_seconds'] != inference):
            raise ValueError('completed supplemental artifacts disagree')
    elif 'artifacts' in run:
        raise ValueError('unsuccessful supplement cannot supply scored artifacts')
    return policy, terminal


def run_supplement(formal, reviewed, path, digest):
    service = check_service(process_identity(os.getpid()), path, digest)
    fd = lock_fd()
    try:
        untouched()
        policy = checked_policy(formal, reviewed, path, digest, live=True)
        request = request_for(policy, verify_inputs=True)
        snapshot = admission(policy['job'])
        state = reserve(policy, path, digest, request, fd, service, snapshot)
        directory = Path(request['output_directory'])
        account = Path(policy['account'])
        with (account / 'lock').open('xb') as account_lock:
            fcntl.flock(account_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            state.update(supervise(state, directory, account, fd))
        check_lock(fd)
        checked_policy(formal, reviewed, path, digest, live=True)
        terminal = publish_terminal(policy, digest, state, request)
        print(json.dumps(dict(status=terminal['status'], attempt_id=ATTEMPT, terminal=policy['terminal'],
                              terminal_sha256=sha256(Path(policy['terminal'])))))
    finally:
        os.close(fd)


def interrupted(signum, _frame):
    raise InterruptedError(f'supplement supervisor received signal {signum}')


def main():
    # Internal scope command. It can run only with the inherited canonical lock
    # and live supervisor's matching reservation, inside the exact owned scope.
    if sys.argv[1:2] == ['worker']:
        worker_exec(int(sys.argv[2]), sys.argv[3], sys.argv[4:])
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('draft-policy', 'preflight', 'run', 'verify'))
    parser.add_argument('--formal', required=True)
    parser.add_argument('--reviewed-stage-sha256', required=True)
    parser.add_argument('--policy')
    parser.add_argument('--reviewed-policy-sha256')
    parser.add_argument('--reviewed-terminal-sha256')
    args = parser.parse_args()
    formal, stage = args.formal, args.reviewed_stage_sha256
    if args.action == 'draft-policy':
        fd = lock_fd()
        try:
            untouched()
            print(json.dumps(policy_for(formal, stage, original_exit()), indent=2, sort_keys=True))
        finally:
            os.close(fd)
        return
    if not args.policy or not args.reviewed_policy_sha256:
        parser.error('reviewed policy and digest required')
    path, digest = args.policy, args.reviewed_policy_sha256
    if args.action == 'preflight':
        fd = lock_fd()
        try:
            untouched()
            policy = checked_policy(formal, stage, path, digest, live=True)
            request_for(policy, verify_inputs=False)
            print(json.dumps(dict(status='preflight_passed', dispatch=False, attempt_id=ATTEMPT,
                                  admission=admission(policy['job']))))
        finally:
            os.close(fd)
    elif args.action == 'verify':
        if not args.reviewed_terminal_sha256:
            parser.error('reviewed terminal digest required')
        _, terminal = verify_terminal(formal, stage, path, digest, args.reviewed_terminal_sha256)
        print(json.dumps(dict(status=terminal['status'], attempt_id=ATTEMPT, dispatch=False)))
    else:
        signal.signal(signal.SIGTERM, interrupted)
        signal.signal(signal.SIGINT, interrupted)
        run_supplement(formal, stage, path, digest)


if __name__ == '__main__':
    main()
