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
    """Retain reviewed PD boundary and independently sum all CMonge attempt walls."""
    from perturbdiff_v2_migration import reconcile as pd_reconcile
    result, after = pd_reconcile(queue, stopped, output)
    jobs = [j for j in queue['jobs'] if j['model'] == 'cmonge']
    if len(jobs) != 9 or {(j['fold'],j['seed']) for j in jobs} != {(f,s) for f in range(3) for s in (17,42,2026)}:
        raise ValueError('exactly nine unique CMonge slots required')
    docs = {}
    for path in sorted((output/'scheduling').rglob('*.json')):
        if path.is_symlink(): raise ValueError('scheduler symlink refused')
        docs[str(path.relative_to(output))] = (digest(path),read(path))
    def mentions(value, job):
        if isinstance(value,dict):
            if (value.get('model'),value.get('fold'),value.get('seed')) == ('cmonge',job['fold'],job['seed']):
                return value != dict(status='not_dispatched',model='cmonge',fold=job['fold'],seed=job['seed'])
            return any(mentions(v,job) for v in value.values())
        if isinstance(value,list): return any(mentions(v,job) for v in value)
        return isinstance(value,str) and (key(job) in value or f"runs/cmonge/fold_{job['fold']}/seed_{job['seed']}" in value)
    for job in jobs:
        name=key(job); state=output/'scheduling'/(name+'.json')
        run=output/'runs/cmonge'/f"fold_{job['fold']}"/f"seed_{job['seed']}"
        evidence={'stop_sha256':digest(stopped/'stop.json'),
                  'post_stop_proof_sha256':digest(stopped/'post-stop-proof.json'),
                  'all_scheduler_json_inspected_sha256':{p:h for p,(h,_) in docs.items()}}
        if state.exists():
            doc=read(state)
            if doc.get('status') not in ('completed','failed','timeout') or not doc.get('attempts'):
                raise ValueError('terminal CMonge attempt evidence required')
            walls=[debit(a.get('wall_seconds')) for a in doc['attempts']]
            consumed=max(sum(walls),debit(doc.get('training_seconds')))
            evidence.update(state_sha256=digest(state),attempt_wall_seconds=walls,
                accounting='max of all durable attempt walls and reported cumulative training; includes setup/inference')
        else:
            if (state.is_symlink() or list((output/'scheduling').glob(name+'*')) or run.exists() or run.is_symlink()
                    or any(mentions(doc,job) for _,doc in docs.values())):
                raise ValueError('zero debit contradicted by dispatch evidence: '+name)
            consumed=0
            evidence['zero_debit_proof']={'state_and_prefix_and_run_absent':True,'queue_member':job}
        result[name]={'consumed_seconds':consumed,'cap_seconds':14400,
            'remaining_seconds':max(0,14400-consumed),'reconciliation_required':False,'evidence':evidence}
    return result,after


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
    v2_commands = {m: [commands[m][0], str(script / adapter), '--config', str(script.parent / 'models' / m / 'v2.config.json')] for m, adapter in [('perturbdiff','perturbdiff_v2_adapter.py'), ('cmonge','cmonge_v2_gpu_adapter.py')]}
    proposed = copy.deepcopy(queue)
    new_frozen = copy.deepcopy(frozen)
    other_jobs = [j for j in queue['jobs'] if j['model'] not in ('perturbdiff','cmonge')]
    if len(other_jobs) != 45 or len(frozen) != 63:
        raise ValueError('63-slot matrix with other45 required')
    for job in proposed['jobs']:
        old = next(e['identity'] for e in frozen if key(e) == key(job))
        current = identity_function(queue['input_root'], job['model'], job['fold'], job['seed'], job.get('command'))
        if current != old:
            raise ValueError('frozen live identity mismatch: ' + key(job))
        if job['model'] not in ('perturbdiff','cmonge'):
            continue
        command = v2_commands[job['model']]
        job['command'] = command
        new = identity_function(queue['input_root'], job['model'], job['fold'], job['seed'], command)
        if new == old:
            raise ValueError('v2 identity unchanged')
        next(e for e in new_frozen if key(e) == key(job))['identity'] = new
        ledger[key(job)]['identity'] = new
        ledger[key(job)]['account_path'] = str(output / 'joint_v2_budget' / key(job))
    # Rehash native legacy bytes; no deserializing or rewriting.
    actual = {}
    for path in sorted(p for model in ('perturbdiff','cmonge') for p in (output / 'runs' / model).rglob('*')):
        if path.is_symlink():
            raise ValueError('legacy symlink refused')
        if path.is_file():
            actual[str(path.relative_to(output))] = {'sha256': digest(path), 'bytes': path.stat().st_size}
    if {p:v for p,v in actual.items() if p.startswith('runs/perturbdiff/')} != after['legacy_inventory']:
        raise ValueError('legacy artifact inventory changed')
    archived = output / 'legacy/perturbdiff-v1'
    if any((output / 'legacy' / (m+'-v1')).exists() for m in ('perturbdiff','cmonge')):
        raise ValueError('legacy archive destination exists')
    plan = {'schema': 'joint-v2-migration-plan-v1', 'activation_authorized': False,
        'independent_review_required': True, 'resume_supported': False,
        'initialization': 'restart_from_seed_at_remaining_original_budget',
        'source_hashes': {n: digest(formal / n) for n in ('queue.json', 'commands.json', 'expected_identities.json')},
        'legacy_inventory': actual,
        'old_new18': [{'run': key(e), 'old_identity': e['identity'], 'new_identity': ledger[key(e)]['identity']}
                     for e in frozen if e['model'] in ('perturbdiff','cmonge')],
        'other45_exact_equality': {
            'jobs': other_jobs == [j for j in proposed['jobs'] if j['model'] not in ('perturbdiff','cmonge')],
            'identities': [e for e in frozen if e['model'] not in ('perturbdiff','cmonge')] ==
                          [e for e in new_frozen if e['model'] not in ('perturbdiff','cmonge')],
            'commands': all(commands.get(j['model']) == j.get('command') for j in other_jobs if j.get('command'))},
        'other45_jobs': other_jobs,
        'other45_identities': [e for e in frozen if e['model'] not in ('perturbdiff','cmonge')],
        'archive_operations': [operation for model in ('perturbdiff','cmonge') for operation in
            [{'source':str(output/'runs'/model),'destination':str(output/'legacy'/(model+'-v1')/'runs')},
             *[{'source':str(p),'destination':str(output/'legacy'/(model+'-v1')/'scheduling'/p.name)}
               for p in sorted((output/'scheduling').glob(model+'_*'))]]],
        'limitations': ['Staging only; joint_v2_apply.py is a separate reviewed transaction; no live file changed.',
            'Independent approval must pin this entire staging inventory before apply.',
            'Runtime scheduler integration and native GPU acceptance remain separate gates.',
            'Legacy v1 predictions/checkpoints never become v2 outcomes.']}
    commands.update(v2_commands)
    documents = {'queue.json': proposed, 'commands.json': commands, 'expected_identities.json': new_frozen,
        'budget-ledger.json': {'schema': 'joint-v2-migration-budget-v1', 'boundary_verified': True,
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
