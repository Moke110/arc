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
            return
        raise ValueError('scheduler lock is unheld')


def check_service():
    state = unit_properties(SUPPLEMENT_UNIT, 'MainPID', 'InvocationID', 'ActiveState',
                            'Restart', 'MemoryMax', 'MemorySwapMax', 'ControlGroup', 'WorkingDirectory')
    if (state.get('ActiveState') != 'active' or state.get('Restart') != 'no' or
            state.get('InvocationID') != os.environ.get('INVOCATION_ID') or
            state.get('MainPID') not in (str(os.getpid()), str(os.getppid())) or
            state.get('MemoryMax') != str(10 * 1024**3) or
            state.get('MemorySwapMax') != str(8 * 1024**3) or
            state.get('WorkingDirectory') != str(ROOT) or
            not state.get('ControlGroup', '').endswith('/app.slice/' + SUPPLEMENT_UNIT) or
            '0::' + state['ControlGroup'] not in Path('/proc/self/cgroup').read_text().splitlines()):
        raise ValueError('launch requires the reviewed dedicated systemd service and resource limits')
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
            '-p', 'MemoryMax=10G', '-p', 'MemorySwapMax=8G', '-p', 'RuntimeMaxSec=16200',
            '-p', 'BindsTo=' + SUPPLEMENT_UNIT, '-p', 'After=' + SUPPLEMENT_UNIT,
            '--', sys.executable, str(SCRIPT), 'worker', str(fd), scope,
            *executed]


def worker_exec(fd, scope, command):
    """The scope must retain the same scheduler lock even if its supervisor dies."""
    check_lock(fd)
    reservation = read_json(safe(OUTPUT / 'account' / 'account.json'))
    if (reservation['status'] != 'reserved' or reservation['scope'] != scope or
            command != reservation['executed_argv'] or not matches(reservation['supervisor']) or
            not re.fullmatch(SCOPE_PREFIX + r'[0-9a-f]{32}\.scope', scope) or
            not Path('/proc/self/cgroup').read_text().strip().endswith('/' + scope)):
        raise ValueError('scope launch is not the reserved worker')
    cgroup = Path('/sys/fs/cgroup') / Path('/proc/self/cgroup').read_text().strip()[3:].lstrip('/')
    if ((cgroup / 'memory.max').read_text().strip() != str(10 * 1024**3) or
            (cgroup / 'memory.swap.max').read_text().strip() != str(8 * 1024**3) or Path.cwd() != ROOT):
        raise ValueError('worker scope or working-directory limits differ from the reviewed recipe')
    request = safe(OUTPUT / 'account' / 'request-0.json')
    if sha256(request) != reservation['injected_request_sha256']:
        raise ValueError('injected request changed before exec')
    write(OUTPUT / 'account' / 'worker-start.json',
          dict(process=process_identity(os.getpid()), scope=scope,
               injected_request_sha256=sha256(request)))
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


def watch(worker, clock, progress_path, stop, transition_path):
    while True:
        code = worker.poll()
        try:
            progress = read_json(safe(progress_path))
        except FileNotFoundError:
            progress = None
        now = time.monotonic()
        old = clock.transition
        reason = clock.observe(now, progress)
        if old is None and clock.transition is not None:
            write(transition_path, clock.timing(now))
        end = time.monotonic()
        if reason or end >= clock.deadline:
            return 'timeout', reason or ('training_timeout' if clock.transition is None else 'inference_timeout'), end
        if code == 124:
            return 'timeout', 'native_training_timeout', end
        if code is not None:
            if code != 0:
                return 'failed', 'worker_nonzero_exit', end
            if clock.transition is None:
                return 'failed', 'worker_exit_without_inference', end
            return 'completed', 'worker_zero_exit', end
        if stop.is_set():
            return 'failed', 'host_resource_pressure', end
        time.sleep(min(.05, max(0, clock.deadline - end)))


def sample_resources(done, stop, path, pid, scope, scope_cgroup):
    history = []
    try:
        with path.open('x', buffering=1) as stream:
            while not done.is_set():
                # A scope can reparent its command away from the systemd-run client.
                process_file = scope_cgroup / 'cgroup.procs'
                try:
                    owned = {pid, *map(int, process_file.read_text().split())}
                except FileNotFoundError:
                    owned = {pid}
                snapshot = process_snapshot(owned_pids=sorted(owned))
                stream.write(json.dumps(dict(snapshot, timestamp=time.time())) + '\n')
                history = [*history[-1:], snapshot]
                if (not snapshot['observable'] or snapshot['external_heavy'] or
                        pressure(history, [dict(owned=True, started=0, id=scope)])['stop']):
                    stop.set()
                done.wait(5)
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


def reserve(policy, digest, request, fd, service, snapshot):
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
    executed = ['/usr/bin/time', '-v', *policy['job']['command'], '--request', str(injected)]
    argv = scope_command(scope, fd, executed)
    state = dict(schema='tenv-cmonge-supplement-account-v1', status='reserved',
                 attempt_id=ATTEMPT, policy_sha256=digest, identity=policy['identity'],
                 caps=policy['caps'], prior_training_seconds=0, credit_available_seconds=0,
                 reservation_count=1, supervisor=process_identity(os.getpid()), scope=scope,
                 scope_command=argv, executed_argv=executed,
                 request_sha256=sha256(directory / 'request.json'),
                 injected_request_sha256=sha256(injected), service=service, admission=snapshot)
    # The protected claim and distinct account are durable before worker creation.
    write(RECORD / 'launch-claim.json', state)
    write(account / 'account.json', state)
    return state


def supervise(state, directory, account, fd):
    scope = state['scope']
    worker = None
    identity = None
    done, stop = threading.Event(), threading.Event()
    sampler = None
    scope_cgroup = (Path('/sys/fs/cgroup') / state['service']['ControlGroup'].lstrip('/')).parent / scope
    with (account / 'watchdog.log').open('xb', buffering=0) as log:
        try:
            clock = PhaseDeadline(time.monotonic())
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
            status, reason, end = watch(worker, clock, directory / 'progress.json',
                                        stop, account / 'inference-transition.json')
        finally:
            done.set()
            if worker is not None:
                cleanup_proof = cleanup(worker, scope, identity)
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
    if (started['scope'] != scope or started['injected_request_sha256'] != state['injected_request_sha256']
            or matches(started['process'])):
        raise ValueError('owned scoped worker start or exit not authenticated')
    peaks = re.findall(r'Maximum resident set size \(kbytes\):\s*(\d+)',
                      (account / 'watchdog.log').read_text(errors='replace'))
    telemetry = [json.loads(line) for line in safe(account / 'telemetry.jsonl').read_text().splitlines()]
    if not telemetry:
        raise ValueError('missing worker resource telemetry')
    sampled_peak = max((sum(t.get('owned_rss_gib', {}).values()) for t in telemetry), default=0)
    return dict(status=status, stop_reason=reason, timing=timing, returncode=worker.returncode,
                process=identity, worker=started, cleanup=cleanup_proof,
                peak_host_rss_gib=max(sampled_peak, max(map(int, peaks), default=0) / 1024**2),
                time_peak_rss_gib=max(map(int, peaks)) / 1024**2 if peaks else None,
                peak_workstation_gpu_used_gib=max(t.get('gpu_used_gib', 0) for t in telemetry))


def publish_terminal(policy, digest, state, request):
    account = Path(policy['account'])
    directory = Path(request['output_directory'])
    if (read_json(safe(directory / 'request.json')) != request or
            read_json(safe(account / 'request-0.json')) != dict(request, prior_training_seconds=0)):
        raise ValueError('request changed during the attempt')
    artifacts, inference = None, None
    if state['status'] == 'completed':
        try:
            artifacts, inference = validate_completion(directory, request, state)
        except (OSError, ValueError, KeyError, TypeError, IndexError) as exc:
            state.update(status='failed', stop_reason='invalid_native_artifacts', artifact_error=str(exc))
    state['timing']['adapter_inference_seconds'] = inference
    run = dict(status=state['status'], identity=policy['identity'], attempt_id=ATTEMPT,
               stop_reason=state['stop_reason'], training_seconds=state['timing']['adapter_training_seconds'],
               inference_seconds=inference, policy_sha256=digest)
    if artifacts is not None:
        run['artifacts'] = artifacts
    validate_terminal_timing(state, run)
    write(directory / 'run.json', run)
    write(account / 'account.json', state)
    terminal = dict(schema='tenv-cmonge-supplement-terminal-v1', status=run['status'],
                    attempt_id=ATTEMPT, identity=policy['identity'], run_directory=str(directory),
                    policy_sha256=digest, original=policy['original'], account=state,
                    account_files=files_at(account), run_sha256=sha256(directory / 'run.json'),
                    claim_sha256=sha256(safe(RECORD / 'launch-claim.json')), files=files_at(directory))
    write(RECORD / 'terminal.json', terminal)
    return terminal


def validate_terminal_timing(state, run):
    def seconds(value):
        return type(value) in (int, float) and math.isfinite(value) and value >= 0

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
    if state['status'] == 'completed':
        if (observed is None or observed >= INFERENCE or wall >= TRAIN + INFERENCE or
                not seconds(native_inference) or native_inference > INFERENCE):
            raise ValueError('completed attempt exceeded a watchdog phase limit')
    elif native_inference is not None:
        raise ValueError('unsuccessful attempt reported completed inference')


def verify_terminal(formal, reviewed_stage_sha256, policy_path, reviewed_policy_sha256,
                    reviewed_terminal_sha256):
    """Authenticate the reviewed terminal outcome and both isolated histories."""
    policy = checked_policy(formal, reviewed_stage_sha256, policy_path, reviewed_policy_sha256, live=False)
    terminal = read_pinned(Path(policy['terminal']), reviewed_terminal_sha256)
    directory = run_directory(OUTPUT, *SLOT)
    account = OUTPUT / 'account'
    state = read_json(safe(account / 'account.json'))
    claim = safe(RECORD / 'launch-claim.json')
    reservation = read_json(claim)
    request = request_for(policy, verify_inputs=False)
    run = read_json(safe(directory / 'run.json'))
    status = state['status']
    if (terminal != dict(schema='tenv-cmonge-supplement-terminal-v1', status=status,
                        attempt_id=ATTEMPT, identity=policy['identity'], run_directory=str(directory),
                        policy_sha256=reviewed_policy_sha256, original=policy['original'], account=state,
                        account_files=files_at(account), run_sha256=sha256(directory / 'run.json'),
                        claim_sha256=sha256(claim), files=files_at(directory)) or
            status not in {'completed', 'failed', 'timeout'} or
            state['reservation_count'] != 1 or state['credit_available_seconds'] != 0 or
            state['prior_training_seconds'] != 0 or state['caps'] != policy['caps'] or
            state['policy_sha256'] != reviewed_policy_sha256 or state['identity'] != policy['identity'] or
            state['attempt_id'] != ATTEMPT or reservation['status'] != 'reserved' or
            any(state[k] != v for k, v in reservation.items() if k != 'status') or
            state['cleanup']['scope'] != state['scope'] or state['cleanup']['active_state'] not in ('inactive', 'failed') or
            read_json(safe(account / 'policy.json')) != policy or
            read_json(safe(directory / 'request.json')) != request or
            read_json(safe(account / 'request-0.json')) != dict(request, prior_training_seconds=0) or
            sha256(directory / 'request.json') != state['request_sha256'] or
            sha256(account / 'request-0.json') != state['injected_request_sha256'] or
            state['executed_argv'] != ['/usr/bin/time', '-v', *policy['job']['command'],
                                       '--request', str(account / 'request-0.json')] or
            run['status'] != status or run['identity'] != policy['identity'] or
            run['attempt_id'] != ATTEMPT or run['policy_sha256'] != reviewed_policy_sha256 or
            run['stop_reason'] != state['stop_reason']):
        raise ValueError('terminal evidence does not match the one-shot supplemental account')
    validate_terminal_timing(state, run)
    if status == 'completed':
        artifacts, inference = validate_completion(directory, request, state)
        if (artifacts != run['artifacts'] or state['returncode'] != 0 or
                run['inference_seconds'] != inference):
            raise ValueError('completed supplemental artifacts disagree')
    elif 'artifacts' in run:
        raise ValueError('unsuccessful supplement cannot supply scored artifacts')
    return policy, terminal


def run_supplement(formal, reviewed, path, digest):
    service = check_service()
    fd = lock_fd()
    try:
        untouched()
        policy = checked_policy(formal, reviewed, path, digest, live=True)
        request = request_for(policy, verify_inputs=True)
        snapshot = admission(policy['job'])
        state = reserve(policy, digest, request, fd, service, snapshot)
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
