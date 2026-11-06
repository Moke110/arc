"""Isolated serial v2 continuation. Preparation never implies restart authorization."""
import argparse
import fcntl
import json
import os
import re
from pathlib import Path
import subprocess
import sys
import time

from perturbdiff_v2_apply import inventory, publish
from perturbdiff_v2_budget import initialize, load_account
from artifacts import (read_json, run_identity, run_directory, validate_run,
                       sha256, rows, checkpoint_artifacts, validate_prediction, TERMINAL)

SCRIPT = Path(__file__).resolve().parent


def key(job):
    return f"{job['model']}_fold{job['fold']}_seed{job['seed']}"


def preflight(formal, reviewed):
    formal = Path(formal).resolve()
    queue = read_json(formal / 'queue.json')
    output = Path(queue['output_root']).resolve()
    journal = read_json(output / 'perturbdiff-v2-apply.json')
    if (journal.get('status') != 'committed' or journal.get('reviewed_sha256') != reviewed
            or journal.get('formal') != str(formal)):
        raise ValueError('committed reviewed migration required')
    if inventory(formal) != journal['new_formal']:
        raise ValueError('committed formal generation changed')
    jobs = queue['jobs']
    if queue.get('max_active_jobs') != 1 or len(jobs) != 63 or len({key(j) for j in jobs}) != 63:
        raise ValueError('exact serial formal matrix required')
    expected = {key(e): e['identity'] for e in read_json(formal / 'expected_identities.json')}
    commands = read_json(formal / 'commands.json')
    for job in jobs:
        if job.get('command') != commands.get(job['model']):
            raise ValueError('queue command mismatch')
        identity = run_identity(queue['input_root'], job['model'], job['fold'], job['seed'], job.get('command'))
        if expected.get(key(job)) != identity:
            raise ValueError('frozen source/config/input/environment identity changed: ' + key(job))
    if sum(j['model'] == 'perturbdiff' for j in jobs) != 9:
        raise ValueError('nine v2 outcomes required')
    if not any('perturbdiff_v2_adapter.py' in a for a in commands['perturbdiff']):
        raise ValueError('v2 adapter required')
    return queue, expected


def outcome(account, ledger_path=None, job=None):
    """Only authenticated terminal account outcomes permit artifact acceptance."""
    state = load_account(Path(account))
    if ledger_path is not None:
        if state['migration_sha256'] != sha256(ledger_path):
            raise ValueError('account ledger binding mismatch')
        if state['key'] != key(job):
            raise ValueError('account job binding mismatch')
    if not state['attempts']:
        return None, state
    attempt = state['attempts'][-1]
    if attempt['status'] == 'timeout':
        return 'timeout', state
    if attempt['status'] == 'exited':
        return ('completed' if attempt['returncode'] == 0 else 'failed'), state
    return None, state


def finish(queue, job, identity, account, ledger_path=None):
    from artifacts import write_json
    directory = run_directory(queue['output_root'], job['model'], job['fold'], job['seed'])
    status, account_state = outcome(account, ledger_path, job)
    if status is None:
        raise ValueError('watchdog is not terminal')
    reason = account_state['attempts'][-1].get('stop_reason', 'watchdog_' + status)
    if status == 'completed':
        try:
            request = read_json(directory / 'request.json')
            validate_prediction(directory, rows(request['predictor_files']['test_cases']), request['n_genes'])
            checkpoint = read_json(directory / 'checkpoint.json')
            checkpoint_artifacts(directory, checkpoint, require_blobs=True)
            if checkpoint.get('identity') != identity:
                raise ValueError('checkpoint identity mismatch')
            run = dict(status='completed', identity=identity,
                       training_seconds=account_state['consumed_seconds'],
                       inference_seconds=checkpoint.get('inference_seconds'),
                       accounting=account_state['accounting'],
                       artifacts={name: sha256(directory / name) for name in
                                  ('prediction.npy', 'case_ids.json', 'checkpoint.json', 'request.json')})
            write_json(directory / 'run.json', run)
            validate_run(queue['input_root'], queue['output_root'], job['model'], job['fold'], job['seed'], job['command'])
        except (ValueError, OSError, KeyError):
            status, reason = 'failed', 'invalid_artifacts'
    state = dict(identity=identity, status=status, stop_reason=reason,
                 training_seconds=account_state['consumed_seconds'],
                 accounting=account_state['accounting'], attempts=account_state['attempts'])
    log = Path(account) / 'watchdog.log'
    if log.exists():
        peaks = re.findall(r'Maximum resident set size \(kbytes\):\s*(\d+)', log.read_text(errors='replace'))
        if peaks:
            state['peak_host_rss_gib'] = max(map(int, peaks)) / 1024 ** 2
    if status != 'completed':
        write_json(directory / 'run.json', state)
    publish(Path(queue['output_root']) / 'scheduling' / (key(job) + '.json'), state)
    return state


def run_v2(formal, reviewed, queue, job, identity):
    from benchmark import predictor_files
    from artifacts import command_evidence, verify_input_files
    from schedule import process_snapshot, admit
    output = Path(queue['output_root']); record = output / 'scheduling'
    directory = run_directory(output, job['model'], job['fold'], job['seed'])
    state_path = record / (key(job) + '.json')
    ledger_path = Path(formal) / 'budget-ledger.json'
    ledger = read_json(ledger_path); account = Path(ledger['ledger'][key(job)]['account_path'])
    if not account.exists():
        initialize(str(ledger_path), sha256(ledger_path), key(job), str(account))
    _, account_state = outcome(account, ledger_path, job)
    if (directory / 'run.json').exists():
        previous = read_json(directory / 'run.json')
        if previous.get('identity') != identity:
            raise ValueError('stale v2 artifact identity')
        if previous.get('status') == 'completed':
            validate_run(queue['input_root'], output, job['model'], job['fold'], job['seed'], job['command'])
            return
        if previous.get('status') in TERMINAL:
            return
    if not account_state['attempts'] and account_state['consumed_seconds'] >= 14400:
        state = dict(status='timeout', identity=identity, training_seconds=14400, stop_reason='historical_budget_exhausted')
        directory.mkdir(parents=True, exist_ok=True)
        publish(directory / 'run.json', state)
        publish(state_path, state)
        return
    if account_state['attempts']:
        status, _ = outcome(account, ledger_path, job)
        if status is None:
            raise ValueError('outstanding watchdog reservation; refuse duplicate and require review')
        finish(queue, job, identity, account, ledger_path)
        return
    blocked = time.monotonic()
    while True:
        snapshot = process_snapshot()
        with (record / 'telemetry.jsonl').open('a') as stream:
            stream.write(json.dumps(dict(snapshot, timestamp=time.time())) + '\n')
        if admit(snapshot, job.get('profile'), (), max_active_jobs=1)['admit']:
            break
        if time.monotonic() - blocked > queue.get('admission_timeout_seconds', 300):
            raise ValueError('resource admission timed out')
        time.sleep(5)
    preflight(formal, reviewed)
    outcome(account, ledger_path, job)
    verify_input_files(queue['input_root'])
    files = predictor_files(queue['input_root'], job['fold'])
    directory.mkdir(parents=True, exist_ok=True)
    request = dict(identity=identity, predictor_files=files, target_kind='delta',
                   n_genes=len(rows(files['genes'])), output_directory=str(directory.resolve()),
                   command_evidence=command_evidence(job['command']),
                   inference_timeout_seconds=job.get('inference_timeout_seconds', 1800))
    publish(directory / 'request.json', request)
    publish(directory / 'run.json', dict(status='running', identity=identity))
    # The short-lived launch CLI exits before polling. Its watchdog session survives
    # scheduler death; an outstanding reservation can never dispatch another worker.
    command = ['systemd-run', '--user', '--scope', '-p', 'MemoryMax=10G', '-p',
               'MemorySwapMax=8G', '--', '/usr/bin/time', '-v', *job['command'], '--request', '{request}']
    env = dict(os.environ, OMP_NUM_THREADS='2', MKL_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2', NUMEXPR_NUM_THREADS='2')
    subprocess.run([sys.executable, str(SCRIPT / 'perturbdiff_v2_budget.py'), 'launch',
                    '--account', str(account), '--request', str(directory / 'request.json'), '--', *command],
                   check=True, env=env)
    launch_deadline = time.monotonic() + 30
    while True:
        status, account_state = outcome(account, ledger_path, job)
        if status is not None:
            finish(queue, job, identity, account, ledger_path)
            return
        if not account_state['attempts'] and time.monotonic() > launch_deadline:
            raise ValueError('watchdog launch not confirmed')
        snapshot = process_snapshot()
        with (record / 'telemetry.jsonl').open('a') as stream:
            stream.write(json.dumps(dict(snapshot, timestamp=time.time())) + '\n')
        state = dict(status='running', identity=identity, training_seconds=account_state['consumed_seconds'],
                     accounting='remaining cap reserved; terminal debit pending', account=str(account))
        publish(state_path, state)
        time.sleep(5)


def run(formal, reviewed):
    from schedule import schedule_queue
    formal = Path(formal).resolve()
    with (SCRIPT.parent / '.scheduler.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        queue, expected = preflight(formal, reviewed)
        record = Path(queue['output_root']) / 'scheduling'; record.mkdir(exist_ok=True)
        # Validate all preserved completed outcomes before any new dispatch.
        for job in queue['jobs']:
            if job['model'] == 'perturbdiff':
                continue
            directory = run_directory(queue['output_root'], job['model'], job['fold'], job['seed'])
            if (directory / 'run.json').exists() and read_json(directory / 'run.json').get('status') == 'completed':
                validate_run(queue['input_root'], queue['output_root'], job['model'], job['fold'], job['seed'], job.get('command'))
            elif job['model'] != 'cmonge':
                raise ValueError('preserved completed outcome missing: ' + key(job))
        for job in queue['jobs']:
            if job['model'] == 'perturbdiff':
                run_v2(formal, reviewed, queue, job, expected[key(job)])
        preflight(formal, reviewed)
        continuation = dict(queue, jobs=[j for j in queue['jobs'] if j['model'] != 'perturbdiff'])
        continuation_path = record / 'v2-cmonge-continuation.json'
        publish(continuation_path, continuation)
        result = schedule_queue(continuation_path)
        # Include v2 terminal states in the standard scheduler summary.
        summary = {key(j): read_json(record / (key(j) + '.json')) for j in queue['jobs']}
        publish(record / 'summary.json', summary)
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--formal', required=True)
    parser.add_argument('--reviewed-stage-sha256', required=True)
    parser.add_argument('--run', action='store_true', help='explicit approved production restart')
    args = parser.parse_args()
    if args.run:
        return run(args.formal, args.reviewed_stage_sha256)
    queue, _ = preflight(args.formal, args.reviewed_stage_sha256)
    print(json.dumps(dict(status='preflight_passed', jobs=len(queue['jobs']), dispatch=False)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
