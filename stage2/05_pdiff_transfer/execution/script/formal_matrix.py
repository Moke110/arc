#!/usr/bin/env python3
"""Freeze and run the Step 05 learned-arm matrix without opening test truth."""
import argparse
from contextlib import contextmanager
import fcntl
import json
from pathlib import Path
import sys
import time

import transfer
from artifacts import SEEDS, read_json, sha256, write_json

REGIMES = transfer.REGIMES
ARMS = ('CM', 'Tahoe', 'Transfer', 'Mix')
PHASES = ('initial-fits', 'tahoe-predictions', 'transfer-fits', 'target-predictions')
SCHEMA = 'step05-formal-matrix-v1'
INVENTORY_SCHEMA = 'step05-prediction-inventory-v1'


def fixed_matrix():
    slots, workers = [], {}

    def job(key, phase, operation, owner, regime, seed, requires=()):
        workers[key] = {'phase': phase, 'operation': operation, 'owner': owner,
                        'regime': regime, 'seed': seed, 'requires': list(requires)}
        return key

    for regime in REGIMES:
        for arm in ARMS:
            for seed in SEEDS:
                key = f'{regime}-{arm.lower()}-{seed}'
                slot = {'id': key, 'regime': regime, 'fold': 0, 'arm': arm, 'seed': seed}
                slots.append(slot)
                if arm == 'CM' and regime != 'warm-start':
                    slot['reuse'] = str(seed)
                    continue
                if arm == 'Tahoe':
                    source = []
                    if regime == 'cold-drug':
                        slot['fit'] = job(key + '-fit', PHASES[0], 'source-fit', 'transfer', regime, seed)
                        source.append(slot['fit'])
                    slot['prediction'] = job(key + '-predict', PHASES[1], 'tahoe-predict', 'transfer',
                                             regime, seed, source)
                    continue
                owner = {'CM': 'warm', 'Transfer': 'transfer', 'Mix': 'mix'}[arm]
                requires = [f'{regime}-tahoe-{seed}-predict'] if arm == 'Transfer' else []
                phase = PHASES[2] if requires else PHASES[0]
                slot['fit'] = job(key + '-fit', phase, 'fit', owner, regime, seed, requires)
                slot['prediction'] = job(key + '-predict', PHASES[3], 'predict', owner, regime, seed, [slot['fit']])
    return slots, workers


def bind(path, expected=None):
    path = Path(path)
    if path.is_file():
        digest = sha256(path)
    else:
        sys.path.insert(0, str(transfer.REPO))
        from stage2.evaluation.history import metadata_digest
        digest = metadata_digest(path)
    if expected is not None and digest != expected:
        raise ValueError(f'frozen artifact hash mismatch: {path}')
    return {'path': str(path), 'sha256': digest}


def check_binding(binding):
    return bind(binding['path'], binding['sha256'])


def case_identity(request):
    from artifacts import rows
    from training_inputs import json_hash

    identity = {}
    for role in ('train', 'validation', 'test'):
        name = f'{role}_cases'
        cases = rows(request['predictor_files'][name])
        order = [[case['case_id'], case['group_id']] for case in cases]
        identity[role] = {**bind(request['predictor_files'][name], request['file_sha256'][name]),
                          'case_count': len(order), 'group_count': len({pair[1] for pair in order}),
                          'ordered_case_group_sha256': json_hash(order)}
    return identity


def target_identity(request, split, inference):
    return {'identity': request['identity'], 'split': split, 'cases': case_identity(request),
            'predictor_files': {key: {'path': path, 'sha256': request['file_sha256'][key]}
                                for key, path in request['predictor_files'].items()},
            'inference': inference, 'n_genes': request['n_genes']}


def reuse_cm(runtime, directory):
    import cold_tenv

    formal = runtime / transfer.STEP04 / 'execution/output/tenv-formal-ticket72'
    release = runtime / transfer.STEP04 / 'output/checkpoints/tenv'
    historical = cold_tenv.RECORD
    qualified = read_json(transfer.RECORD / 'cold_record/replay.json')['local_artifact_sha256']
    comparison = read_json(cold_tenv.pinned(directory, 'comparison.json',
                                           qualified[f'{directory.name}/comparison.json']))
    plan_sha = sha256(historical / 'plan.json')
    plan = read_json(cold_tenv.pinned(formal, 'plan.json', plan_sha))
    index = read_json(cold_tenv.pinned(release, 'checkpoints.json', sha256(historical / 'checkpoints.json')))
    inventory = read_json(cold_tenv.pinned(formal, 'outcomes.json', sha256(historical / 'outcomes.json')))
    entry = plan['folds'][0]
    fold = cold_tenv.fold_inputs(formal, entry, 'formal')
    rows_by_seed = {row['seed']: row for row in inventory['outcomes'] if row['fold'] == 0}
    index_by_seed = {row['identity']['seed']: row for row in index['checkpoints'] if row['identity']['fold'] == 0}
    reuse = {}
    for seed in SEEDS:
        values, row = cold_tenv.verify_reused_cm(formal, release, index_by_seed[seed], entry, fold,
                                                 rows_by_seed[seed], plan_sha)
        del values
        if row != comparison['cm'][str(seed)]:
            raise ValueError(f'historical CM reuse differs from its qualification for seed {seed}')
        reuse[str(seed)] = row
    return reuse


def qualification_profiles():
    paths = [transfer.RECORD / path for path in (
        'warm_record/replay.json', 'cold_tahoe_record/qualification.json',
        'tahoe_zero_shot_record/qualification.json', 'transfer_record/qualification.json',
        'mix_record/qualification.json',
    )]
    for dataset in ('tenv', 'tahoe'):
        record = transfer.REPO / transfer.STEP04 / 'execution/formal_record' / dataset
        paths.extend(record / f'slots/fold-0/seed-{seed}/run.json' for seed in SEEDS)
    return {'evidence': [bind(path) for path in paths], 'profiles': {
        'source-fit': {'host_peak_gib': 7, 'gpu_peak_gib': 7},
        'fit': {'host_peak_gib': 5, 'gpu_peak_gib': 7},
        'tahoe-predict': {'host_peak_gib': 4, 'gpu_peak_gib': 4},
        'predict': {'host_peak_gib': 4, 'gpu_peak_gib': 4},
    }, 'measured_peak_basis_gib': {'historical_formal_tahoe_rss': 5.901,
                                   'historical_formal_training_workstation_vram': 6.216,
                                   'bounded_mix_fit_rss': 3.148,
                                   'bounded_tahoe_prediction_workstation_vram_two_jobs': 3.301}}


def prepare(runtime, root):
    import subprocess
    import mix
    from formal import executable, launch_identity
    from training_inputs import json_hash

    if (root / 'plan.json').exists():
        raise ValueError('formal plan already exists; verify the frozen plan')
    if Path(__file__).resolve().parent != runtime / transfer.STEP05 / 'execution/script':
        raise ValueError('prepare must run from the qualification checkout')
    formal = runtime / transfer.STEP04 / 'execution/output'
    output, _, _ = transfer.records(runtime)
    worker_roots = {name: root / name for name in ('warm', 'mix', 'transfer')}
    if not worker_roots['warm'].exists():
        transfer.warm.prepare(runtime, formal / 'tenv-formal-ticket72', worker_roots['warm'], False)
    if not worker_roots['mix'].exists():
        mix.prepare(runtime, worker_roots['mix'], 'formal')
    if not worker_roots['transfer'].exists():
        transfer.prepare(runtime, worker_roots['transfer'], 'formal')
    plans = {name: bind(path / 'plan.json') for name, path in worker_roots.items()}
    warm_plan, warm_request = transfer.warm.verified_plan(worker_roots['warm'], plans['warm']['sha256'], runtime)
    mix_plan, mix_context = mix.verify(runtime, worker_roots['mix'], plans['mix']['sha256'])
    transfer_plan, context, source, source_sha = transfer.verify(
        runtime, worker_roots['transfer'], plans['transfer']['sha256'])
    if any(plan['stage'] != 'formal' for plan in (warm_plan, mix_plan, transfer_plan)):
        raise ValueError('formal matrix requires all three formal worker plans')
    source_path = output / 'tahoe-zero-shot-ticket78/source.json'
    source_check = root / 'historical_tahoe_verified.json'
    if not source_check.exists():
        transfer.tahoe.verify_source(runtime, formal / 'tahoe-formal-ticket73', source_check)
    bind(source_check, source_sha)
    historical_source = bind(source_path, source_sha)
    cold_inputs = transfer.cold_inputs(runtime)
    _, filtered, filtered_request, filtered_evidence = cold_inputs
    filtered_release = read_json(filtered / 'release.json')
    cm_reuse = reuse_cm(runtime, output / 'cold-tenv-ticket76')
    targets, recipient_plans = {}, {}
    for regime in REGIMES:
        request = context[regime][0]
        recipient, recipient_sha = transfer.recipient_slot(output, regime)
        name = 'plan.json' if regime == 'warm-start' else 'split.json'
        recipient_plans[regime] = bind(recipient / name, recipient_sha)
        cm_request = warm_request if regime == 'warm-start' else request
        cm_split = (bind(worker_roots['warm'] / 'split.json', warm_plan['split_sha256'])
                    if regime == 'warm-start' else recipient_plans[regime])
        target = target_identity(cm_request, cm_split, transfer_plan['targets'][regime]['inference'])
        expected_cases = (1402, 460) if regime == 'warm-start' else (1401, 462)
        test = target['cases']['test']
        compared = case_identity(request)
        if ((test['case_count'], test['group_count']) != expected_cases
                or cm_request['identity'] != request['identity']
                or cm_request['file_sha256'] != request['file_sha256']
                or mix_context[regime]['file_sha256'] != request['file_sha256']
                or any(row['ordered_case_group_sha256'] != compared[role]['ordered_case_group_sha256']
                       for role, row in target['cases'].items())
                or warm_plan['inference'] != mix_plan['inference']
                or mix_plan['inference'] != target['inference']):
            raise ValueError(f'{regime} arms do not share the frozen tenv split, files and cases')
        target['recipient_plan'] = recipient_plans[regime]
        target['request_sha256_by_arm'] = {
            'CM': json_hash(cm_request), 'Tahoe': json_hash(request), 'Transfer': json_hash(request),
            'Mix': json_hash(mix_plan['requests'][regime]),
        }
        target['request_files'] = {
            'CM': bind(worker_roots['warm'] / 'request.json') if regime == 'warm-start' else
            bind(recipient / 'predictor_request.json'), 'Tahoe': bind(
                recipient / ('request.json' if regime == 'warm-start' else 'predictor_request.json')),
        }
        pooled = mix_plan['requests'][regime]
        target['mix_source'] = {key: {'path': path, 'sha256': pooled['predictor_sha256']['tahoe'][key]}
                                for key, path in pooled['tahoe'].items()}
        targets[regime] = target
    # Mix's request differs because its Tahoe source is filtered.
    for field in ('identity', 'split', 'cases', 'predictor_files', 'inference', 'n_genes'):
        if targets['cold-drug'][field] != targets['cold-cm-drug'][field]:
            raise ValueError('the two cold regimes must use exactly the same target split')

    scripts = Path(__file__).resolve().parent
    python = executable(runtime)
    launch = launch_identity(runtime)
    extra_files = [*scripts.glob('*.py'), runtime / transfer.STEP04 / 'execution/script/evaluate.py',
                   runtime / 'stage2/03_modes_tenv/execution/script/schedule.py', Path(sys.executable)]
    launch.update({str(path): sha256(path) for path in sorted(extra_files)})
    revision = subprocess.check_output(['git', '-C', str(runtime), 'rev-parse', 'HEAD'], text=True).strip()
    slots, workers = fixed_matrix()
    admission = qualification_profiles()
    for key, worker in workers.items():
        op, owner, regime, seed = (worker[field] for field in ('operation', 'owner', 'regime', 'seed'))
        fit = op in ('fit', 'source-fit')
        worker['queue_job'] = {
            'id': key, 'kind': 'training' if fit else 'gpu', 'cwd': str(runtime), 'memory_gib': 10,
            'disk_gib': 5 if fit else 2, 'timeout_seconds': 15000 if fit else 1800,
            'profile': admission['profiles'][op],
            'command': [str(python), str(scripts / {'warm': 'warm_tenv', 'transfer': 'transfer',
                                                    'mix': 'mix'}[owner]) + '.py',
                        '--runtime', str(runtime), op, str(worker_roots[owner]),
                        *([regime] if owner != 'warm' and op != 'source-fit' else []),
                        str(seed), '--plan-sha256', plans[owner]['sha256']],
        }
        if op == 'source-fit':
            directory = worker_roots[owner] / 'source' / f'seed-{seed}' / 'fit'
        elif op == 'tahoe-predict':
            directory = transfer.tahoe_directory(runtime, worker_roots[owner], regime, seed, 'formal')
        else:
            slot_path = worker_roots[owner] / 'slots'
            if owner != 'warm':
                slot_path /= regime
            directory = slot_path / f'seed-{seed}' / ('fit' if fit else 'prediction')
        worker['directory'] = str(directory)

    for slot in slots:
        regime, arm, seed = (slot[key] for key in ('regime', 'arm', 'seed'))
        raw = warm_plan['formal_recipes'][str(seed)]
        slot['target_recipe'] = {**bind(worker_roots['warm'] / raw['path'], raw['sha256']),
                                 'resolved_sha256': raw['config_sha256']}
        if arm in ('Tahoe', 'Transfer'):
            if regime == 'cold-drug':
                source_fit = workers[f'{regime}-tahoe-{seed}-fit']['directory']
                item = filtered_release['recipes'][str(seed)]
                slot['source'] = {'fit_directory': source_fit, 'seed': seed,
                                  'release': bind(filtered / 'release.json',
                                                  filtered_evidence[f'{filtered.name}/release.json']),
                                  'request': bind(filtered / 'request.json',
                                                  filtered_evidence[f'{filtered.name}/request.json']),
                                  'recipe': {**bind(filtered / item['path'], item['sha256']),
                                             'resolved_sha256': item['resolved_sha256']}}
            else:
                entry = source['checkpoints'][str(seed)]
                producer = entry['producer']
                source_root = formal / 'tahoe-formal-ticket73'
                tahoe_plan = read_json(source_root / 'plan.json')
                source_fold = tahoe_plan['folds'][0]
                slot['source'] = {'source_manifest': historical_source, **entry,
                                  'request': bind(source_root / source_fold['request'], producer['request_sha256']),
                                  'recipe': bind(source_root / source_fold['recipes'][str(seed)]['path'],
                                                 producer['raw_recipe_sha256']),
                                  'history': bind(Path(entry['path']).parent / 'run.json', producer['run_sha256'])}
        if arm == 'Tahoe':
            del slot['target_recipe']

    plan = {'schema': SCHEMA, 'runtime': str(runtime), 'root': str(root), 'producer_revision': revision,
            'launch_files': launch, 'worker_plans': plans, 'targets': targets, 'slots': slots,
            'workers': workers, 'reuse_cm': cm_reuse, 'admission': admission, 'phases': list(PHASES),
            'filtered_source_files': {key: {'path': path, 'sha256': filtered_request['file_sha256'][key]}
                                      for key, path in filtered_request['predictor_files'].items()},
            'attempt_limit': 2, 'minimum_initial_free_disk_gib': 30}
    write_json(root / 'plan.json', plan)
    print(f'Frozen Step 05 matrix SHA-256 {sha256(root / "plan.json")}', flush=True)


def read_plan(root, digest):
    from formal import pinned

    plan = read_json(pinned(root, 'plan.json', digest))
    if (plan['schema'], plan['root'], plan['phases']) != (SCHEMA, str(root), list(PHASES)):
        raise ValueError('unexpected formal matrix identity')
    return plan


def verify_bindings(items):
    checked = {}
    for item in items:
        path, digest = item['path'], item['sha256']
        if path not in checked:
            checked[path] = bind(path)['sha256']
        if checked[path] != digest:
            raise ValueError(f'frozen artifact hash mismatch: {path}')


def cm_bindings(reuse):
    checkpoint = Path(reuse['checkpoint_directory'])
    prediction = Path(reuse['prediction_manifest'])
    bindings = [{'path': str(prediction), 'sha256': reuse['prediction_manifest_sha256']}]
    for directory, names in ((checkpoint, {
        'best.pt': 'checkpoint_sha256', 'run.json': 'run_sha256', 'replay.json': 'replay_sha256',
        'request.json': 'request_sha256', 'recipe.json': 'raw_recipe_sha256',
    }), (prediction.parent, {'prediction.npy': 'prediction_sha256', 'cases.json': 'cases_sha256'})):
        bindings.extend({'path': str(directory / name), 'sha256': reuse[field]} for name, field in names.items())
    return bindings


def verify_plan(plan):
    bindings = [{'path': path, 'sha256': digest} for path, digest in plan['launch_files'].items()]
    bindings.extend(plan['worker_plans'].values())
    bindings.extend(plan['admission']['evidence'])
    bindings.extend(plan['filtered_source_files'].values())
    for target in plan['targets'].values():
        bindings.extend((target['split'], target['recipient_plan'], *target['request_files'].values(),
                         *target['predictor_files'].values(), *target['mix_source'].values()))
    for slot in plan['slots']:
        if 'target_recipe' in slot:
            bindings.append(slot['target_recipe'])
        if 'source' in slot:
            source = slot['source']
            bindings.extend((source['request'], source['recipe']))
            if 'source_manifest' in source:
                bindings.extend((source['source_manifest'], source['history'],
                                 {'path': source['path'], 'sha256': source['checkpoint_sha256']}))
            else:
                bindings.append(source['release'])
    for row in plan['reuse_cm'].values():
        bindings.extend(cm_bindings(row))
    verify_bindings(bindings)


def verify_worker(plan, key, prior):
    from artifacts import rows
    from evaluate import load_array, np, prediction
    from formal import completed_fit, verify_fit_seal
    from training_inputs import json_hash

    worker = plan['workers'][key]
    op, owner, regime, seed = (worker[field] for field in ('operation', 'owner', 'regime', 'seed'))
    directory = Path(worker['directory'])
    plan_sha = plan['worker_plans'][owner]['sha256']
    arm = {'warm': 'CM', 'mix': 'Mix', 'transfer': 'Tahoe' if op in ('source-fit', 'tahoe-predict')
           else 'Transfer'}[owner]
    slot = next(row for row in plan['slots'] if (row['regime'], row['arm'], row['seed']) == (regime, arm, seed))
    target = plan['targets'][regime]

    if op in ('fit', 'source-fit'):
        verify_fit_seal(directory, plan_sha)
        run = read_json(directory / 'run.json')
        if op == 'source-fit':
            identity = read_json(slot['source']['request']['path'])['identity']
            recipe = slot['source']['recipe']
        else:
            identity = target['identity']
            if arm == 'Mix':
                identity = read_json(plan['worker_plans']['mix']['path'])['requests'][regime]['identity']
            elif arm == 'Transfer':
                identity = {**identity, 'arm': arm, 'regime': regime}
            recipe = slot['target_recipe']
        expected = {**identity, 'seed': seed, 'stage': 'formal'}
        completed_fit(directory, expected, run['recipe'])
        if run['recipe_sha256'] != recipe['resolved_sha256']:
            raise ValueError(f'{key} did not use its frozen recipe')
        result = {'identity': expected, 'recipe_sha256': run['recipe_sha256'], 'best_step': run['best_step'],
                  'reload': run['reload'], 'fit': bind(directory / 'fit.json'),
                  'checkpoint': bind(directory / 'best.pt', run['checkpoint_sha256']),
                  'history': bind(directory / 'run.json')}
        if arm == 'Transfer':
            bridge = prior[worker['requires'][0]]['artifacts']
            initialization = run['initialization']
            if (initialization['bridge_sha256'] != bridge['bridge']['sha256']
                    or initialization['mapped_state_sha256'] != bridge['mapped_state_sha256']
                    or initialization['source_checkpoint_sha256'] != bridge['source_checkpoint_sha256']
                    or initialization['recipient_request_sha256'] != target['request_sha256_by_arm']['Transfer']):
                raise ValueError(f'{key} did not initialize from its matching Tahoe prediction')
            result['initialization'] = initialization
        return result

    path = directory / 'prediction.json'
    manifest = read_json(path)
    cases = rows(target['cases']['test']['path'])
    case_order = [[case['case_id'], case['group_id']] for case in cases]
    case_binding = bind(directory / 'cases.json', manifest['cases_sha256'])
    if (read_json(case_binding['path']) != case_order
            or json_hash(case_order) != target['cases']['test']['ordered_case_group_sha256']
            or manifest['case_count'] != target['cases']['test']['case_count']
            or manifest['group_count'] != target['cases']['test']['group_count']):
        raise ValueError(f'{key} prediction does not cover the frozen tenv cases in order')
    result = {'prediction': bind(path), 'cases': case_binding}
    if op == 'tahoe-predict':
        if regime == 'cold-drug':
            fit = prior[worker['requires'][0]]['artifacts']
            expected_sha = fit['checkpoint']['sha256']
            expected_history = fit['history']['sha256']
            expected_recipe = fit['recipe_sha256']
        else:
            source = slot['source']
            expected_sha = source['checkpoint_sha256']
            expected_history = source['producer']['run_sha256']
            expected_recipe = source['recipe_sha256']
            if manifest['source_sha256'] != source['source_manifest']['sha256']:
                raise ValueError(f'{key} differs from the pinned historical Tahoe source')
        bridge = bind(directory / 'bridge.json', manifest['bridge_sha256'])
        mapped = read_json(bridge['path'])
        source = mapped['source']
        if (manifest['schema'] != transfer.tahoe.PREDICTION_SCHEMA
                or (manifest['regime'], manifest['arm'], manifest['seed']) != (regime, 'Tahoe', seed)
                or manifest['recipient_sha256'] != target['recipient_plan']['sha256']
                or mapped['recipient_identity'] != target['identity']
                or mapped['recipient_request_sha256'] != target['request_sha256_by_arm']['Tahoe']
                or manifest['recipient_request_sha256'] != mapped['recipient_request_sha256']
                or mapped['mapped_state_sha256'] != manifest['mapped_state_sha256']
                or source['identity']['seed'] != seed
                or source['checkpoint_sha256'] != expected_sha
                or source['recipe_sha256'] != expected_recipe
                or source['producer']['run_sha256'] != expected_history):
            raise ValueError(f'{key} bridge differs from the frozen Tahoe and recipient identities')
        if regime == 'cold-drug' and any(manifest['source_parent_overlap'][role][heldout]
                                         for role in ('train', 'validation') for heldout in ('validation', 'test')):
            raise ValueError(f'{key} exposed a held-out tenv parent to the filtered source')
        result.update(bridge=bridge, mapped_state_sha256=mapped['mapped_state_sha256'],
                      source_checkpoint_sha256=expected_sha)
        values = load_array(directory / 'prediction.npy', (len(cases), target['n_genes']))
    else:
        fit = prior[worker['requires'][0]]['artifacts']
        if (manifest['run_sha256'] != fit['history']['sha256']
                or manifest['checkpoint_sha256'] != fit['checkpoint']['sha256']):
            raise ValueError(f'{key} prediction did not use its verified sealed fit')
        values = prediction(directory, {'prediction': path.name, 'seed': seed,
                                       'prediction_manifest_sha256': result['prediction']['sha256']},
                            {'identity': {k: v for k, v in fit['identity'].items() if k not in ('seed', 'stage')},
                             'recipes': {seed: fit['recipe_sha256']}, 'cases': cases,
                             'n_genes': target['n_genes']}, 'formal', plan_sha)
        result['replay'] = bind(directory / 'replay.json', manifest['replay_sha256'])
    if any(not np.isfinite(values[start:start + 32]).all() for start in range(0, len(values), 32)):
        raise ValueError(f'{key} has nonfinite prediction values')
    result['values'] = bind(directory / 'prediction.npy', manifest['prediction_sha256'])
    return result


@contextmanager
def locked(path, nonblocking=False):
    with path.open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX | (fcntl.LOCK_NB if nonblocking else 0))
        yield


def append_event(path, event):
    with path.open('a') as stream:
        stream.write(json.dumps({'timestamp': time.time(), **event}, allow_nan=False) + '\n')


def phase_paths(root, phase):
    directory = root / 'phases' / phase
    return directory / 'queue.json', directory / 'scheduler', directory / 'outcomes.json'


def control(root, action, reason=None):
    with locked(root / 'control.lock'):
        state_path = root / 'control.json'
        state = read_json(state_path) if state_path.exists() else {'paused': False, 'phase': None}
        phase = state['phase']
        report = {}
        if phase:
            _, record, _ = phase_paths(root, phase)
            path = record / 'status.json'
            jobs = read_json(path)['jobs'] if path.exists() else {}
            report['active'] = [key for key, row in jobs.items() if row['status'] == 'running']
            report['pending'] = [key for key, row in jobs.items() if row['status'] == 'pending']
        sealed = (root / 'prediction_inventory.json').exists()
        if action in ('pause', 'resume'):
            if sealed:
                raise ValueError('prediction inventory is sealed')
            if not reason or not reason.strip():
                raise ValueError('pause and resume require a recorded operator reason')
            state['paused'] = action == 'pause'
            if phase:
                if state['paused']:
                    (record / 'pause').touch()
                else:
                    (record / 'pause').unlink(missing_ok=True)
            write_json(state_path, state)
            append_event(root / 'control.jsonl', {'action': action, 'reason': reason, **state, **report,
                                                'active_disposition': 'drain under pressure limits'
                                                if report.get('active') and state['paused'] else None})
        return {**state, **report, 'sealed': sealed}


class CappedScheduler:
    """Reuse Step 04's execution loop and add formal admission evidence and pause gating."""

    def __init__(self, runtime):
        import importlib.util

        path = runtime / transfer.STEP04 / 'execution/script/schedule.py'
        spec = importlib.util.spec_from_file_location('step05_capped_step04', path)
        self.scheduler = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.scheduler)

    def sample(self, root):
        return self.scheduler.snapshot(root)

    def read(self, queue):
        record, jobs = self.scheduler.load_queue(queue)
        return self.scheduler.load_status(record, jobs)['jobs']

    def retry(self, queue, keys):
        step04 = self.scheduler
        record, jobs = step04.load_queue(queue)
        with step04.workstation_lock():
            status = step04.load_status(record, jobs)
            step04.recover(record, status)
            for key in keys:
                row = status['jobs'][key]
                if row['status'] not in ('failed', 'timeout', 'interrupted'):
                    raise ValueError(f'cannot retry {key} from {row["status"]}')
                row['status'] = 'pending'
                append_event(record / 'retry.jsonl', {'job': key, 'previous_attempt': row['attempts'][-1],
                                                       'command': row['identity']['command']})
            step04.publish(record, status)

    def run(self, root, queue):
        from contextlib import ExitStack
        import signal
        from unittest.mock import patch

        step04 = self.scheduler
        record, jobs = step04.load_queue(queue)
        decide, launch, publish = step04.decision, step04.launch, step04.publish
        last = None
        with ExitStack() as dispatch_lock:
            def decide_guarded(sample, job, active, **_):
                dispatch_lock.enter_context(locked(root / 'control.lock'))
                if (record / 'pause').exists():
                    verdict = {'admit': False, 'reason': 'paused'}
                else:
                    verdict = decide(sample, job, active, reserve_active_gpu=True)
                if not verdict['admit']:
                    dispatch_lock.close()
                return verdict

            def launch_guarded(job, attempt):
                try:
                    return launch(job, attempt)
                finally:
                    dispatch_lock.close()

            def publish_logged(path, status):
                nonlocal last
                publish(path, status)
                decision = status.get('last_decision') or {}
                timestamp = decision.get('timestamp')
                if timestamp is not None and timestamp != last:
                    append_event(record / 'admission.jsonl', {'decision': decision})
                    last = timestamp

            with (patch.object(step04, 'decision', decide_guarded),
                  patch.object(step04, 'launch', launch_guarded),
                  patch.object(step04, 'publish', publish_logged)):
                signals = (signal.SIGTERM, signal.SIGINT)
                handlers = [signal.getsignal(sig) for sig in signals]
                try:
                    return step04.run(record, jobs, {}, None)
                finally:
                    for sig, handler in zip(signals, handlers):
                        signal.signal(sig, handler)


def run_phase(root, phase, phase_workers, plan, done, backend, validate):
    queue, record, terminal = phase_paths(root, phase)
    ready, outcomes = [], {}
    for key, worker in phase_workers.items():
        blocked = [dep for dep in worker['requires'] if done[dep]['status'] != 'completed']
        if blocked:
            causes = sorted({cause for dep in blocked for cause in done[dep].get('root_failures', [dep])})
            outcomes[key] = {'status': 'blocked', 'failed_prerequisites': blocked,
                             'root_failures': causes, 'attempts': []}
        else:
            ready.append(worker['queue_job'])
    payload = {'record_dir': str(record), 'jobs': ready}
    with locked(root / 'control.lock'):
        control_path = root / 'control.json'
        state = read_json(control_path) if control_path.exists() else {'paused': False, 'phase': None}
        if state['paused']:
            return None
        if terminal.exists():
            sealed = read_json(terminal)
            if ready:
                bind(queue, sealed['queue']['sha256'])
                if read_json(queue) != payload:
                    raise ValueError(f'sealed phase queue differs from the frozen plan: {phase}')
                check_binding(sealed['status'])
            if sealed['phase'] != phase or sealed['outcomes'].keys() != phase_workers.keys():
                raise ValueError(f'sealed phase outcomes differ from the frozen plan: {phase}')
            return sealed
        if ready:
            record.mkdir(parents=True, exist_ok=True)
            if queue.exists():
                if read_json(queue) != payload:
                    raise ValueError(f'phase queue changed after it was frozen: {phase}')
            else:
                write_json(queue, payload)
        if state['phase'] != phase:
            write_json(control_path, {**state, 'phase': phase})
            append_event(root / 'control.jsonl', {'action': 'phase-ready', 'phase': phase, 'jobs': len(ready)})
    while ready:
        if control(root, 'status')['paused']:
            return None
        rows = backend.read(queue)
        retriable = [key for key, row in rows.items()
                     if row['status'] in ('failed', 'timeout', 'interrupted')
                     and len(row['attempts']) < plan['attempt_limit']]
        if retriable:
            backend.retry(queue, retriable)
        pending = any(row['status'] in ('pending', 'running') for row in rows.values())
        if retriable or pending:
            code = backend.run(root, queue)
            if code:
                return None
            continue
        for key, row in rows.items():
            attempts = [{**attempt, 'log_sha256': sha256(Path(attempt['log']))}
                        if Path(attempt['log']).is_file() else attempt for attempt in row['attempts']]
            outcome = {'status': row['status'], 'attempts': attempts}
            if row['status'] == 'completed':
                try:
                    outcome['artifacts'] = validate(plan, key, done)
                except (ArithmeticError, KeyError, ValueError, TypeError, IndexError, OSError, EOFError) as error:
                    outcome.update(status='invalid', reason=str(error))
            outcomes[key] = outcome
        break
    seal = {'phase': phase, 'outcomes': outcomes}
    if ready:
        seal.update(queue=bind(queue), status=bind(record / 'status.json'), evidence={
            name: bind(record / name) for name in ('telemetry.jsonl', 'admission.jsonl', 'retry.jsonl')
            if (record / name).is_file()})
    write_json(terminal, seal)
    return seal


def inventory_slots(plan, done):
    slots = []
    for slot in plan['slots']:
        row = {field: slot[field] for field in ('id', 'regime', 'fold', 'arm', 'seed')}
        if 'reuse' in slot:
            row.update(status='reused', reuse=plan['reuse_cm'][slot['reuse']])
        else:
            prediction = done[slot['prediction']]
            row.update(status=prediction['status'], prediction=prediction)
            if 'fit' in slot:
                row['fit'] = done[slot['fit']]
        slots.append(row)
    return slots


def run_matrix(root, digest, plan, backend, validate=verify_worker):
    with locked(root / 'run.lock', nonblocking=True):
        inventory = root / 'prediction_inventory.json'
        if inventory.exists():
            verify_inventory(root, digest, plan)
            return True
        done, seals = {}, {}
        for phase in PHASES:
            phase_workers = {key: row for key, row in plan['workers'].items() if row['phase'] == phase}
            result = run_phase(root, phase, phase_workers, plan, done, backend, validate)
            if result is None:
                return False
            done.update(result['outcomes'])
            seals[phase] = bind(phase_paths(root, phase)[2])
        slots = inventory_slots(plan, done)
        with locked(root / 'control.lock'):
            state = read_json(root / 'control.json')
            if state['paused']:
                return False
            write_json(root / 'control.json', {**state, 'phase': None})
            append_event(root / 'control.jsonl', {'action': 'sealed', 'slots': len(slots)})
            write_json(inventory, {'schema': INVENTORY_SCHEMA, 'plan_sha256': digest, 'slots': slots,
                                   'phases': seals, 'control_evidence': bind(root / 'control.jsonl')})
        print(f'Sealed Step 05 prediction inventory SHA-256 {sha256(inventory)}', flush=True)
        return True


def verify_inventory(root, digest, plan):
    inventory = read_json(root / 'prediction_inventory.json')
    if ((inventory['schema'], inventory['plan_sha256']) != (INVENTORY_SCHEMA, digest)
            or set(inventory['phases']) != set(PHASES)):
        raise ValueError('inventory does not match the frozen matrix plan')
    bindings = [inventory['control_evidence']]
    done = {}
    for name, phase in inventory['phases'].items():
        check_binding(phase)
        seal = read_json(phase['path'])
        expected = {key for key, worker in plan['workers'].items() if worker['phase'] == name}
        if seal['phase'] != name or seal['outcomes'].keys() != expected:
            raise ValueError(f'phase seal does not account for the planned workers: {name}')
        done.update(seal['outcomes'])
        bindings.extend(seal[key] for key in ('queue', 'status') if key in seal)
        bindings.extend(seal.get('evidence', {}).values())
        for outcome in seal['outcomes'].values():
            if outcome['status'] == 'completed':
                bindings.extend(item for item in outcome['artifacts'].values()
                                if isinstance(item, dict) and item.keys() == {'path', 'sha256'})
            for attempt in outcome['attempts']:
                if 'log_sha256' in attempt:
                    bindings.append({'path': attempt['log'], 'sha256': attempt['log_sha256']})
    if inventory['slots'] != inventory_slots(plan, done):
        raise ValueError('inventory slots differ from the frozen plan and terminal worker evidence')
    for row in plan.get('reuse_cm', {}).values():
        bindings.extend(cm_bindings(row))
    verify_bindings(bindings)
    return inventory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime', type=Path, default=transfer.REPO)
    commands = parser.add_subparsers(dest='command', required=True)
    for name in ('prepare', 'verify', 'preflight', 'run', 'status', 'pause', 'resume', 'verify-inventory'):
        command = commands.add_parser(name)
        command.add_argument('root', type=Path)
        if name != 'prepare':
            command.add_argument('--plan-sha256', required=True)
        if name in ('pause', 'resume'):
            command.add_argument('--reason', required=True)
    args = parser.parse_args()
    root, runtime = args.root.resolve(), args.runtime.resolve()
    try:
        if args.command == 'prepare':
            root.parent.mkdir(parents=True, exist_ok=True)
            with locked(root.parent / f'.{root.name}-prepare.lock'):
                prepare(runtime, root)
            return 0
        plan = read_plan(root, args.plan_sha256)
        if runtime != Path(plan['runtime']):
            raise ValueError('runtime differs from the frozen qualification checkout')
        if args.command in ('status', 'pause', 'resume'):
            print(json.dumps(control(root, args.command, getattr(args, 'reason', None)), indent=2))
        elif args.command == 'verify-inventory':
            verify_inventory(root, args.plan_sha256, plan)
        else:
            verify_plan(plan)
            if args.command in ('preflight', 'run'):
                backend = CappedScheduler(runtime)
                if args.command == 'preflight' or not (root / 'control.json').exists():
                    sample = backend.sample(root)
                    print(json.dumps({'preflight': sample}), flush=True)
                    if (not sample.get('observable')
                            or sample['available_gib'] < 12
                            or sample['disk_free_gib'] < plan['minimum_initial_free_disk_gib']):
                        raise ValueError('formal preflight requires observable telemetry, 12 GiB RAM and 30 GiB disk')
                if args.command == 'run' and not run_matrix(root, args.plan_sha256, plan, backend):
                    print(json.dumps(control(root, 'status'), indent=2))
                    return 0 if control(root, 'status')['paused'] else 2
        return 0
    except (ArithmeticError, KeyError, ValueError, TypeError, IndexError, OSError, EOFError) as error:
        parser.exit(1, f'formal matrix {args.command} rejected: {error}\n')


if __name__ == '__main__':
    sys.exit(main())
