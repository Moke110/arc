"""Stage a stopped-boundary migration for independent review; never apply or launch.

Run with the formal scheduler's interpreter so environment identities stay fixed.
The resulting budget ledger is a proposal until its exact hash is reviewed.
"""
import argparse
import copy
import fcntl
import json
import math
from pathlib import Path

from perturbdiff_v2_boundary import digest


def read(path):
    return json.loads(path.read_text())


def key(job):
    return f"{job['model']}_fold{job['fold']}_seed{job['seed']}"


def debit(value):
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or value < 0:
        raise ValueError('invalid historical debit')
    return value


def reconcile(queue, stopped, output):
    """Reconcile only durable stopped data and explicit absent-dispatch evidence."""
    stop = read(stopped / 'stop.json')
    proof = read(stopped / 'post-stop-proof.json')
    after = read(stopped / 'after/boundary.json')
    if (not stop.get('process_exit_verified') or not proof.get('process_exit_verified') or
            proof.get('stop_sha256') != digest(stopped / 'stop.json') or
            not proof.get('nonblocking_scheduler_lock_acquired') or
            proof.get('old_scheduler_exists') is not False or proof.get('old_worker_exists') is not False or
            'Restart=no' not in proof.get('service', '') or 'ActiveState=inactive' not in proof.get('service', '')):
        raise ValueError('verified stopped boundary required')
    selected = [j for j in queue['jobs'] if j['model'] == 'perturbdiff']
    expected = {(f, s) for f in range(3) for s in (17, 42, 2026)}
    if len(selected) != 9 or {(j['fold'], j['seed']) for j in selected} != expected:
        raise ValueError('exactly nine unique PerturbDiff slots required')
    scheduler_documents = {}
    for path in sorted((output / 'scheduling').rglob('*.json')):
        if path.is_symlink():
            raise ValueError('scheduler symlink refused')
        scheduler_documents[str(path.relative_to(output))] = (digest(path), read(path))

    def mentions(value, job):
        if isinstance(value, dict):
            if (value.get('model'), value.get('fold'), value.get('seed')) == ('perturbdiff', job['fold'], job['seed']):
                return value != {'status': 'not_dispatched', 'model': 'perturbdiff',
                                 'fold': job['fold'], 'seed': job['seed']}
            return any(mentions(v, job) for v in value.values())
        if isinstance(value, list):
            return any(mentions(v, job) for v in value)
        if not isinstance(value, str):
            return False
        run_suffix = f"runs/perturbdiff/fold_{job['fold']}/seed_{job['seed']}"
        return key(job) in value or run_suffix in value

    result = {}
    for job in selected:
        name = key(job)
        state = output / 'scheduling' / (name + '.json')
        run = output / 'runs/perturbdiff' / f"fold_{job['fold']}" / f"seed_{job['seed']}"
        row = stop['ledger'][name]
        evidence = {'stop_sha256': digest(stopped / 'stop.json'),
                    'post_stop_proof_sha256': digest(stopped / 'post-stop-proof.json')}
        if row['consumed_seconds'] is None:
            # A missing state alone is insufficient. Refuse any attempt/log/job or
            # run artifact, including an empty run directory or symlink.
            attempts = list((output / 'scheduling').glob(name + '*'))
            if (state.exists() or attempts or run.exists() or run.is_symlink() or
                    'scheduling/' + name + '.json' in after['source_state_hashes'] or
                    any(mentions(doc, job) for _, doc in scheduler_documents.values())):
                raise ValueError('zero debit contradicted by dispatch evidence: ' + name)
            consumed = 0
            evidence['zero_debit_proof'] = {'queue_member': job, 'state_absent': str(state),
                'all_scheduler_prefix_entries_absent': name, 'run_directory_absent': str(run),
                'stopped_state_inventory_absent': True,
                'all_scheduler_json_inspected_sha256': {p: h for p, (h, _) in scheduler_documents.items()}}
        else:
            consumed = debit(row['consumed_seconds'])
            relative = str(state.relative_to(output))
            if digest(state) != after['source_state_hashes'][relative]:
                raise ValueError('state changed after stop: ' + name)
            current = read(state)
            if current['status'] == 'running' or consumed < debit(current['training_seconds']):
                raise ValueError('live or undercharged state: ' + name)
            evidence['state_sha256'] = digest(state)
            evidence['accounting'] = 'conservative full wall upper bound from verified stop; includes setup/inference'
        result[name] = {'consumed_seconds': consumed, 'cap_seconds': 14400,
                        'remaining_seconds': max(0, 14400-consumed),
                        'reconciliation_required': False, 'evidence': evidence}
    return result, after


def stage(formal, stopped, destination):
    # Frozen shared identity implementation stays untouched.
    from artifacts import run_identity
    formal = formal.resolve()
    queue_path = formal / 'queue.json'
    queue = read(queue_path)
    output = Path(queue['output_root']).resolve()
    with (output / 'scheduling/.scheduler.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        return stage_locked(formal, stopped.resolve(), destination.resolve(), queue, output, run_identity)


def stage_locked(formal, stopped, destination, queue, output, identity_function):
    if destination.exists():
        raise ValueError('staging destination already exists')
    if digest(formal / 'queue.json') != read(stopped / 'after/boundary.json')['queue_sha256']:
        raise ValueError('queue changed since stopped boundary')
    ledger, after = reconcile(queue, stopped, output)
    frozen = read(formal / 'expected_identities.json')
    commands = read(formal / 'commands.json')
    script = Path(__file__).resolve().parent
    command = [commands['perturbdiff'][0], str(script / 'perturbdiff_v2_adapter.py'),
               '--config', str(script.parent / 'models/perturbdiff/v2.config.json')]
    proposed = copy.deepcopy(queue)
    new_frozen = copy.deepcopy(frozen)
    other_jobs = [j for j in queue['jobs'] if j['model'] != 'perturbdiff']
    if len(other_jobs) != 54 or len(frozen) != 63:
        raise ValueError('63-slot matrix with other54 required')
    for job in proposed['jobs']:
        old = next(e['identity'] for e in frozen if key(e) == key(job))
        current = identity_function(queue['input_root'], job['model'], job['fold'], job['seed'], job.get('command'))
        if current != old:
            raise ValueError('frozen live identity mismatch: ' + key(job))
        if job['model'] != 'perturbdiff':
            continue
        job['command'] = command
        new = identity_function(queue['input_root'], job['model'], job['fold'], job['seed'], command)
        if new == old:
            raise ValueError('v2 identity unchanged')
        next(e for e in new_frozen if key(e) == key(job))['identity'] = new
        ledger[key(job)]['identity'] = new
        ledger[key(job)]['account_path'] = str(output / 'perturbdiff_v2_budget' / key(job))
    # Rehash native legacy bytes; no deserializing or rewriting.
    actual = {}
    for path in sorted((output / 'runs/perturbdiff').rglob('*')):
        if path.is_symlink():
            raise ValueError('legacy symlink refused')
        if path.is_file():
            actual[str(path.relative_to(output))] = {'sha256': digest(path), 'bytes': path.stat().st_size}
    if actual != after['legacy_inventory']:
        raise ValueError('legacy artifact inventory changed')
    archived = output / 'legacy/perturbdiff-v1'
    if archived.exists():
        raise ValueError('legacy archive destination exists')
    plan = {'schema': 'perturbdiff-v2-migration-plan-v1', 'activation_authorized': False,
        'independent_review_required': True, 'resume_supported': False,
        'initialization': 'restart_from_seed_at_remaining_original_budget',
        'source_hashes': {n: digest(formal / n) for n in ('queue.json', 'commands.json', 'expected_identities.json')},
        'legacy_inventory': actual,
        'old_new9': [{'run': key(e), 'old_identity': e['identity'], 'new_identity': ledger[key(e)]['identity']}
                     for e in frozen if e['model'] == 'perturbdiff'],
        'other54_exact_equality': {
            'jobs': other_jobs == [j for j in proposed['jobs'] if j['model'] != 'perturbdiff'],
            'identities': [e for e in frozen if e['model'] != 'perturbdiff'] ==
                          [e for e in new_frozen if e['model'] != 'perturbdiff'],
            'commands': all(commands.get(j['model']) == j.get('command') for j in other_jobs if j.get('command'))},
        'other54_jobs': other_jobs,
        'other54_identities': [e for e in frozen if e['model'] != 'perturbdiff'],
        'archive_operations': [{'source': str(output / 'runs/perturbdiff'), 'destination': str(archived / 'runs')},
            *[{'source': str(p), 'destination': str(archived / 'scheduling' / p.name)}
              for p in sorted((output / 'scheduling').glob('perturbdiff_*'))]],
        'limitations': ['No apply transaction implemented; no live file changed.',
            'Independent approval must pin this entire staging inventory before apply.',
            'Runtime scheduler integration and native GPU acceptance remain separate gates.',
            'Legacy v1 predictions/checkpoints never become v2 outcomes.']}
    commands['perturbdiff'] = command
    documents = {'queue.json': proposed, 'commands.json': commands, 'expected_identities.json': new_frozen,
        'budget-ledger.json': {'schema': 'perturbdiff-v2-migration-budget-v1', 'boundary_verified': True,
            'activation_authorized': False, 'ledger': ledger}, 'plan.json': plan}
    destination.mkdir(parents=True)
    for name, doc in documents.items():
        with (destination / name).open('x') as stream:
            json.dump(doc, stream, indent=2, sort_keys=True)
            stream.write('\n')
    hashes = {name: digest(destination / name) for name in documents}
    (destination / 'staging-sha256.json').write_text(json.dumps(hashes, indent=2, sort_keys=True) + '\n')
    return hashes


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--formal', type=Path, required=True)
    parser.add_argument('--stopped', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(stage(args.formal, args.stopped, args.destination), sort_keys=True))
