"""Bind one reviewed terminal CMonge supplement to the unchanged evaluator."""
import argparse
from contextlib import contextmanager
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import stat
import tempfile

from artifacts import SEEDS, TERMINAL, read_json, run_directory, sha256, verify_input_files, write_json
import benchmark
from joint_v2_runtime import key, outcome, preflight

MODELS = ('no_change', 'delta_mean', 'prnet', 'cpa', 'chemcpa', 'cmonge', 'perturbdiff')
SLOTS = tuple((model, fold, seed) for model in MODELS for fold in range(3) for seed in SEEDS)
TARGET = ('cmonge', 0, 17)
ATTEMPT = 'cmonge-v2-fold0-seed17-supplement-1'
FITTED = set(MODELS) - {'no_change', 'delta_mean'}
FITTED_V2 = {'cmonge', 'perturbdiff'}
EVALUATION_FILES = ('per_group.csv', 'per_run.csv', 'summary.csv', 'best_runs.csv',
                    'pr_bands.csv', 'pr_sources.npz', 'predictions.json')


def checked_path(path, *, directory=False):
    path = Path(os.path.abspath(path))
    if path.resolve() != path:
        raise ValueError(f'linked path: {path}')
    info = path.lstat()
    if directory:
        if not stat.S_ISDIR(info.st_mode):
            raise ValueError(f'expected directory: {path}')
    elif not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
        raise ValueError(f'expected privately owned regular file: {path}')
    return path


def digest(path):
    return sha256(checked_path(path))


def pinned_json(path, reviewed):
    data = checked_path(path).read_bytes()
    if hashlib.sha256(data).hexdigest() != reviewed:
        raise ValueError(f'reviewed digest mismatch: {path}')
    return json.loads(data)


def files(directory):
    """Hash the complete flat run directory, including native histories."""
    return {path.name: digest(path) for path in sorted(checked_path(directory, directory=True).iterdir())}


def run_files(directory, run, command):
    inventory = files(directory)
    if run['status'] == 'completed':
        required = {'request.json', 'checkpoint.json', 'case_ids.json', 'prediction.npy'}
        if set(run['artifacts']) != required or any(run['artifacts'][name] != inventory.get(name) for name in required):
            raise ValueError(f'completed artifact mismatch: {directory}')
        checkpoint = read_json(directory / 'checkpoint.json')
        blobs = checkpoint.get('blobs', {})
        if (checkpoint.get('identity') != run['identity'] or (command and not blobs)
                or any(inventory.get(name) != value for name, value in blobs.items())):
            raise ValueError(f'checkpoint identity/blob mismatch: {directory}')
    return inventory


def original_records(formal, queue, expected):
    original = checked_path(queue['output_root'], directory=True)
    commands = read_json(checked_path(formal / 'commands.json'))
    jobs = {(job['model'], job['fold'], job['seed']): job for job in queue['jobs']}
    if len(queue['jobs']) != len(SLOTS) or set(jobs) != set(SLOTS) or set(commands) != FITTED:
        raise ValueError('expected the fixed seven-entry matrix')
    if set(expected) != {key(job) for job in jobs.values()}:
        raise ValueError('unexpected expected-identity matrix')
    ledger_path = checked_path(formal / 'budget-ledger.json')
    ledger = read_json(ledger_path)['ledger']
    summary_path = checked_path(original / 'scheduling' / 'summary.json')
    summary = read_json(summary_path)
    if set(summary) != set(expected):
        raise ValueError('original scheduler summary is incomplete')
    records = []
    for slot in SLOTS:
        job = jobs[slot]
        job_key = key(job)
        directory = run_directory(original, *slot)
        run = read_json(checked_path(directory / 'run.json'))
        scheduler_path = checked_path(original / 'scheduling' / (job_key + '.json'))
        state = read_json(scheduler_path)
        status = run.get('status')
        identity = expected[job_key]
        if (status not in TERMINAL or state != summary[job_key] or state.get('status') != status
                or run.get('identity') != identity or state.get('identity') != identity
                or job.get('command') != commands.get(slot[0])
                or (status != 'completed' and state.get('stop_reason') != run.get('stop_reason'))):
            raise ValueError(f'original run/scheduler disagreement: {job_key}')
        record = dict(model=slot[0], fold=slot[1], seed=slot[2], identity=identity, status=status,
                      command=job.get('command'), stop_reason=state.get('stop_reason'),
                      directory=str(directory), files=run_files(directory, run, job.get('command')),
                      scheduler_sha256=digest(scheduler_path))
        if slot[0] in FITTED_V2:
            row = ledger[job_key]
            account = checked_path(row['account_path'], directory=True)
            if account != original / 'joint_v2_budget' / job_key or row['identity'] != identity:
                raise ValueError('original account path mismatch')
            account_files = {name: digest(account / name) for name in ('account.json', 'migration.json')}
            account_status, account_state = outcome(account, ledger_path, job)
            # A failed artifact validation can downgrade a successful process exit.
            if (account_status != status
                    and not (account_status == 'completed' and status == 'failed'
                             and state.get('stop_reason') == 'invalid_artifacts')):
                raise ValueError(f'original account outcome mismatch: {job_key}')
            if (account_state['attempts'] != state.get('attempts')
                    or account_state['consumed_seconds'] != state.get('training_seconds')):
                raise ValueError(f'original debit mismatch: {job_key}')
            record['account_files'] = account_files
            if slot == TARGET:
                attempt = account_state['attempts'][-1]
                prior = row['consumed_seconds']
                remaining = row['cap_seconds'] - prior
                if (status != 'timeout' or row['cap_seconds'] != 14400 or run['training_seconds'] != 14400
                        or job.get('inference_timeout_seconds') != 1800
                        or attempt['prior_training_seconds'] != prior
                        or not math.isclose(attempt['reserved_seconds'], remaining, abs_tol=1e-8)):
                    raise ValueError('original target cumulative cap proof mismatch')
                record['budget'] = dict(cap_seconds=14400, prior_debit_seconds=prior,
                                        remaining_at_launch_seconds=remaining, cumulative_debit_seconds=14400,
                                        attempt_wall_seconds=attempt['wall_seconds'],
                                        peak_host_rss_gib=state.get('peak_host_rss_gib'))
        records.append(record)
    return records, digest(summary_path)


def selection(formal, reviewed_stage, policy_path, reviewed_policy, reviewed_terminal):
    from cmonge_supplement import verify_terminal

    formal = checked_path(formal, directory=True)
    queue, expected = preflight(formal, reviewed_stage)
    policy, terminal = verify_terminal(formal, reviewed_stage, policy_path, reviewed_policy, reviewed_terminal)
    original, summary_sha = original_records(formal, queue, expected)
    previous = next(row for row in original if slot_key(row) == TARGET)
    selected_root = checked_path(policy['output_root'], directory=True)
    directory = run_directory(selected_root, *TARGET)
    run = read_json(checked_path(directory / 'run.json'))
    selected_files = run_files(directory, run, previous['command'])
    original_proof = terminal['original']
    if (original_proof['run_files'] != previous['files']
            or original_proof['scheduler_sha256'] != previous['scheduler_sha256']
            or original_proof['summary_sha256'] != summary_sha
            or any(original_proof['account_files'].get(name) != value
                   for name, value in previous['account_files'].items())):
        raise ValueError('original timeout differs from reviewed terminal provenance')
    if (selected_root != Path(queue['output_root']).parent / 'supplement-cmonge-v2-fold0-seed17'
            or Path(terminal['run_directory']) != directory
            or terminal['attempt_id'] != ATTEMPT or terminal['identity'] != previous['identity']
            or policy['original'] != terminal['original'] or terminal['policy_sha256'] != reviewed_policy
            or terminal['status'] not in {'completed', 'timeout', 'failed'}
            or terminal['status'] != run['status'] or run['identity'] != previous['identity']
            or terminal['files'] != selected_files or terminal['run_sha256'] != selected_files['run.json']):
        raise ValueError('supplement terminal selection mismatch')
    selected = {field: previous[field] for field in ('model', 'fold', 'seed', 'identity', 'command')}
    selected.update(status=run['status'], stop_reason=run.get('stop_reason'), directory=str(directory),
                    files=selected_files, attempt_id=ATTEMPT)
    return dict(schema='tenv-cmonge-supplement-overlay-v1', reviewed_stage_sha256=reviewed_stage,
                reviewed_policy_sha256=reviewed_policy, reviewed_terminal_sha256=reviewed_terminal,
                original_output=queue['output_root'], input_root=queue['input_root'],
                input_sha256=digest(Path(queue['input_root']) / 'manifest.json'),
                scheduler_summary_sha256=summary_sha, original=original, terminal=terminal,
                selected=[selected if row is previous else row for row in original])


def slot_key(row):
    return row['model'], row['fold'], row['seed']


def overlay_root(selection_record):
    return Path(selection_record['original_output']).parent / 'final-with-cmonge-supplement'


def materialize(record):
    output = overlay_root(record)
    checked_path(output.parent, directory=True)
    if os.path.lexists(output):
        raise ValueError('evaluation overlay already exists')
    with tempfile.TemporaryDirectory(prefix='.cmonge-evaluation-', dir=output.parent) as temporary:
        staged = Path(temporary) / 'overlay'
        staged.mkdir(mode=0o700)
        for row in record['selected']:
            source = Path(row['directory'])
            target = run_directory(staged, *slot_key(row))
            target.mkdir(parents=True)
            for name in row['files']:
                shutil.copyfile(checked_path(source / name), target / name)
            if files(target) != row['files'] or files(source) != row['files']:
                raise ValueError(f'run changed during copy: {source}')
        original_timeout = next(row for row in record['original'] if slot_key(row) == TARGET)
        if files(original_timeout['directory']) != original_timeout['files']:
            raise ValueError('original timeout changed during overlay preparation')
        write_json(staged / 'selection.json', record)
        staged.rename(output)
    return output


def verify_copy(output, record):
    actual = set()
    for model_dir in checked_path(output / 'runs', directory=True).iterdir():
        for fold_dir in checked_path(model_dir, directory=True).iterdir():
            for run_dir in checked_path(fold_dir, directory=True).iterdir():
                actual.add(run_dir)
    expected = {run_directory(output, *slot_key(row)): row['files'] for row in record['selected']}
    if actual != set(expected) or any(files(path) != hashes for path, hashes in expected.items()):
        raise ValueError('private overlay slot/file mismatch')


def verify_seal(output, record):
    entries = [{field: row[field] for field in ('model', 'fold', 'seed', 'status', 'command')}
               | {'run_sha256': row['files']['run.json']} for row in record['selected']]
    expected = dict(input_sha256=record['input_sha256'], models=list(MODELS), runs=entries,
                    seeds=list(SEEDS), folds=[0, 1, 2], sealed=True)
    if read_json(checked_path(output / 'predictions.json')) != expected:
        raise ValueError('frozen seal differs from the reviewed attempt selection')


def evaluation_proof(output, record, selection_sha):
    verify_seal(output, record)
    evaluation = read_json(checked_path(output / 'evaluation.json'))
    sources = {name: digest(output / name) for name in EVALUATION_FILES}
    if evaluation['artifacts'] != sources or evaluation['workers'] != 2:
        raise ValueError('scored artifact checksum or worker count mismatch')
    return dict(schema='tenv-cmonge-supplement-evaluation-v1', selection_sha256=selection_sha,
                original_timeout=next(row for row in record['original'] if slot_key(row) == TARGET),
                selected_attempt=record['terminal'], evaluation_sha256=digest(output / 'evaluation.json'),
                seal_sha256=sources['predictions.json'], artifacts=sources)


@contextmanager
def scheduler_idle(formal):
    lock_path = checked_path(Path(formal).parent / '.scheduler.lock')
    with lock_path.open('r') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield


def evaluate(output, record, selection_sha):
    from deliver import FIGURES

    if any(os.path.lexists(output / name) for name in (*EVALUATION_FILES, *FIGURES, 'evaluation.json',
                                                       'evaluation-provenance.json')):
        raise ValueError('evaluation files already exist; preserve them for review')
    commands = {row['model']: row['command'] for row in record['selected'] if row['command'] is not None}
    benchmark.freeze_predictions(record['input_root'], output, list(MODELS), commands)
    verify_seal(output, record)
    benchmark.verify_barrier(record['input_root'], output)
    benchmark.evaluate(record['input_root'], output, workers=2)
    verify_copy(output, record)
    write_json(output / 'evaluation-provenance.json', evaluation_proof(output, record, selection_sha))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'verify', 'evaluate', 'publish'))
    parser.add_argument('--formal', type=Path, required=True)
    parser.add_argument('--reviewed-stage-sha256', required=True)
    parser.add_argument('--policy', type=Path, required=True)
    parser.add_argument('--reviewed-policy-sha256', required=True)
    parser.add_argument('--reviewed-terminal-sha256', required=True)
    parser.add_argument('--reviewed-selection-sha256')
    parser.add_argument('--reviewed-evaluation-sha256')
    args = parser.parse_args()
    if args.action != 'prepare' and not args.reviewed_selection_sha256:
        parser.error('a separately reviewed selection digest is required')
    if args.action == 'publish' and not args.reviewed_evaluation_sha256:
        parser.error('a separately reviewed evaluation digest is required')

    with scheduler_idle(args.formal):
        record = selection(args.formal, args.reviewed_stage_sha256, args.policy,
                           args.reviewed_policy_sha256, args.reviewed_terminal_sha256)
        output = overlay_root(record)
        selection_path = output / 'selection.json'
        proof_path = output / 'evaluation-provenance.json'
        if args.action == 'prepare':
            verify_input_files(record['input_root'])
            materialize(record)
            print(json.dumps(dict(output=str(output), selection_sha256=digest(selection_path))))
            return
        if pinned_json(selection_path, args.reviewed_selection_sha256) != record:
            raise ValueError('selection evidence changed since review')
        verify_copy(output, record)
        if args.action == 'evaluate':
            evaluate(output, record, args.reviewed_selection_sha256)
            print(json.dumps(dict(evaluation_provenance_sha256=digest(proof_path))))
        elif args.action == 'publish':
            from deliver import publish
            proof = pinned_json(proof_path, args.reviewed_evaluation_sha256)
            if proof != evaluation_proof(output, record, args.reviewed_selection_sha256):
                raise ValueError('evaluation evidence changed since review')
            publish(output, output.parent, supplement=(record, proof, args.reviewed_evaluation_sha256))
        else:
            verify_input_files(record['input_root'])
            if os.path.lexists(output / 'predictions.json'):
                verify_seal(output, record)
                benchmark.verify_barrier(record['input_root'], output)
            print(json.dumps(dict(status='verified', slots=len(SLOTS))))


if __name__ == '__main__':
    main()
