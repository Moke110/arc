#!/usr/bin/env python3
"""Plan and analyze validation-only Step 04 PerturbDiff pilots."""
import argparse
from itertools import product
from pathlib import Path
import sys

from train import REPO, STEP03

sys.path.insert(0, str(REPO / STEP03 / 'execution/script'))
from artifacts import SEEDS, read_json, sha256, write_json
from pilot_evidence import match_trajectory, read_profile, strict_best
from training_inputs import json_hash
from training_policy import MONITOR, Selection

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE.parent / 'pilot_record' / 'protocol.json'


def training_command(runtime, trainer, *arguments):
    return [str(runtime / 'stage2/_envs/perturbdiff/bin/python'),
            str(trainer), '--runtime-root', str(runtime), *map(str, arguments)]


def trainer_identity(runtime, trainer):
    trainer_files = [trainer, *(trainer.with_name(name) for name in (
        'native_training.py', 'training_inputs.py', 'training_policy.py', 'tahoe_inputs.py',
    ))]
    runtime_files = [runtime / STEP03 / 'execution' / relative for relative in (
        'models/perturbdiff/v2.config.json', 'models/perturbdiff/env.actual.lock',
        'script/perturbdiff_v2_adapter.py', 'script/perturbdiff_v2_protocol.py',
        'script/predictor_inputs.py', 'script/artifacts.py',
    )]
    vendor = runtime / 'pert_model/PerturbDiff/manifest.sha256'
    return {str(path): sha256(path) for path in [*trainer_files, *runtime_files, vendor] if path.is_file()}


def queue_job(runtime, job_id, kind, seconds, disk_gib, command):
    return {'id': job_id, 'kind': kind, 'cwd': str(runtime), 'timeout_seconds': seconds,
            'disk_gib': disk_gib, 'command': command}


def profile_recipe(protocol, fold, interval, steps, fit_seconds):
    return {'stage': 'pilot', 'fold': fold, 'seed': protocol['pilot_seed'], 'device': 'cuda',
            'max_steps': steps, 'validation_interval': interval, 'minimum_steps': steps + 1,
            'patience_steps': steps + 1, 'min_delta': 0, 'max_training_seconds': fit_seconds}


def fit_job(runtime, trainer, phase, fold, interval, request, steps, fit_seconds, protocol):
    job_id = f'fold-{fold}-v{interval}'
    recipe = phase / 'recipes' / f'{job_id}.json'
    run = phase / 'runs' / job_id
    write_json(recipe, profile_recipe(protocol, fold, interval, steps, fit_seconds))
    command = training_command(runtime, trainer, 'fit', request, recipe, run)
    return queue_job(runtime, job_id, 'training', fit_seconds + protocol['scheduler_grace_seconds'], 2, command)


def validate_protocol(protocol):
    if (protocol['schema'] != 'step04-pilot-protocol-v1' or protocol['folds'] != list(range(3))
            or protocol['pilot_seed'] not in SEEDS or protocol['validation_intervals'] != [250, 500]):
        raise ValueError('invalid pilot protocol inventory')
    short, long = protocol['profile_steps'], protocol['calibration_steps']
    if (type(short) is not int or type(long) is not int or short < 1000 or short % 500
            or long <= short or long % 500):
        raise ValueError('pilot horizons must be matched multiples of 500 steps')
    for name in ('minimum_steps', 'patience_steps', 'min_delta_fraction_of_first_loss'):
        values = protocol[name]
        if not values or values != sorted(set(values)) or values[0] <= 0:
            raise ValueError(f'{name} must be increasing and positive')


def plan_profile(runtime, trainer, output, dataset, release, requests, protocol_path):
    protocol = read_json(protocol_path)
    validate_protocol(protocol)
    folds = protocol['folds']
    phase = output / 'profile'
    jobs = []
    if release is not None:
        requests = {fold: output / 'requests' / f'fold-{fold}.json' for fold in folds}
        for fold, request in requests.items():
            command = training_command(runtime, trainer, f'prepare-{dataset}', release, request, '--fold', fold)
            jobs.append(queue_job(runtime, f'prepare-fold-{fold}', 'cpu', 300, 1, command))
    else:
        identities = [read_json(path)['identity'] for path in requests]
        datasets = {identity['dataset'] for identity in identities}
        if len(datasets) != 1 or sorted(identity['fold'] for identity in identities) != folds:
            raise ValueError('supply one request per fold of one dataset')
        dataset = datasets.pop()
        requests = {identity['fold']: path for identity, path in zip(identities, requests)}
    output.mkdir(parents=True)
    for fold in folds:
        for interval in protocol['validation_intervals']:
            jobs.append(fit_job(runtime, trainer, phase, fold, interval, requests[fold], protocol['profile_steps'],
                                protocol['profile_fit_seconds'], protocol))
    write_json(output / 'plan.json', {'dataset': dataset, 'runtime_root': str(runtime), 'trainer': str(trainer),
                                     'requests': {str(fold): str(path) for fold, path in requests.items()},
                                     'trainer_files': trainer_identity(runtime, trainer),
                                     'protocol': protocol})
    write_json(phase / 'queue.json', {'record_dir': 'scheduler', 'jobs': jobs})


def checked_profiles(output, phase_name, plan, intervals):
    phase = output / phase_name
    protocol = plan['protocol']
    steps = protocol[f'{phase_name}_steps']
    budget = protocol[f'{phase_name}_fit_seconds']
    runtime = Path(plan['runtime_root'])
    trainer = Path(plan.get('trainer', HERE / 'train.py'))
    queue = read_json(phase / 'queue.json')['jobs']
    fits = [job for job in queue if job['kind'] == 'training']
    jobs = {job['id']: job for job in fits}
    status = read_json(phase / 'scheduler/status.json')['jobs']
    expected = {f'fold-{fold}-v{interval}' for fold, grid in intervals.items() for interval in grid}
    if jobs.keys() != expected or len(fits) != len(expected):
        raise ValueError(f'{phase_name}: queue has an unexpected fit inventory')
    traces = {}
    for fold, grid in intervals.items():
        request = Path(plan['requests'][str(fold)])
        identity = read_json(request)['identity']
        if (identity['dataset'], identity['fold']) != (plan['dataset'], fold):
            raise ValueError(f'{phase_name}: request dataset or fold mismatch')
        traces[fold] = []
        for interval in grid:
            job_id = f'fold-{fold}-v{interval}'
            job = jobs[job_id]
            state = status.get(job_id)
            if not state or state['status'] != 'completed':
                raise ValueError(f'{phase_name}/{job_id}: scheduler has no completed attempt')
            command = training_command(runtime, trainer, 'fit', request, phase / 'recipes' / f'{job_id}.json',
                                       phase / 'runs' / job_id)
            if (any(state['identity'].get(key) != value for key, value in job.items())
                    or job['command'] != command
                    or job['timeout_seconds'] != budget + protocol['scheduler_grace_seconds']):
                raise ValueError(f'{phase_name}/{job_id}: scheduler command does not match frozen queue')
            recipe = profile_recipe(protocol, fold, interval, steps, budget)
            if read_json(Path(job['command'][-2])) != recipe:
                raise ValueError(f'{phase_name}/{job_id}: recipe differs from frozen protocol')
            traces[fold].append(read_profile(phase, job_id, request, recipe, state))
    return traces


def cadence_choice(fold, dense, sparse, protocol):
    match_trajectory(dense, sparse, protocol)
    gain = sparse['best_validation_loss'] - dense['best_validation_loss']
    tolerance = abs(sparse['best_validation_loss']) * protocol['loss_tolerance_fraction']
    overhead = dense['overhead_fraction']
    if overhead <= protocol['overhead_target_fraction_of_training'] and gain >= -tolerance:
        chosen, reason = dense, 'dense_cadence_within_cost_and_loss_targets'
    elif overhead <= protocol['overhead_limit_fraction_of_training'] and gain > tolerance:
        chosen, reason = dense, 'dense_cadence_improves_validation_within_cost_limit'
    else:
        chosen, reason = sparse, 'sparse_cadence_preferred_by_cost_and_loss_rule'
    return {'fold': fold, 'validation_interval': chosen['recipe']['validation_interval'],
            'reason': reason, 'dense_best_loss_gain': gain, 'loss_tolerance': tolerance,
            'profiles': [dense, sparse]}


def cadence_report(output, plan):
    protocol = plan['protocol']
    traces = checked_profiles(output, 'profile', plan,
                              {fold: protocol['validation_intervals'] for fold in protocol['folds']})
    return {'dataset': plan['dataset'], 'protocol_sha256': json_hash(protocol),
            'plan_sha256': sha256(output / 'plan.json'),
            'queue_sha256': sha256(output / 'profile/queue.json'),
            'folds': [cadence_choice(fold, *traces[fold], protocol) for fold in protocol['folds']]}


def calibration_plan(output, plan, cadence):
    runtime = Path(plan['runtime_root'])
    trainer = Path(plan.get('trainer', HERE / 'train.py'))
    fingerprints = trainer_identity(runtime, trainer)
    if 'trainer_files' in plan and fingerprints != plan['trainer_files']:
        raise ValueError('pilot trainer changed since profile planning')
    return {'protocol_sha256': cadence['protocol_sha256'], 'profile_plan_sha256': cadence['plan_sha256'],
            'cadence_sha256': sha256(output / 'cadence.json'),
            'calibration_queue_sha256': sha256(output / 'calibration/queue.json'),
            'trainer_files': fingerprints,
            'trainer_hash_capture': ('profile_planning' if 'trainer_files' in plan else 'calibration_planning'),
            'profile_run_sha256': {str(choice['fold']): [run['run_sha256'] for run in choice['profiles']]
                                   for choice in cadence['folds']}}


def plan_calibration(output, plan):
    report = cadence_report(output, plan)
    protocol = plan['protocol']
    runtime = Path(plan['runtime_root'])
    trainer = Path(plan.get('trainer', HERE / 'train.py'))
    if 'trainer_files' in plan and trainer_identity(runtime, trainer) != plan['trainer_files']:
        raise ValueError('pilot trainer changed since profile planning')
    phase = output / 'calibration'
    phase.mkdir()
    jobs = []
    for choice in report['folds']:
        fold, interval = choice['fold'], choice['validation_interval']
        jobs.append(fit_job(runtime, trainer, phase, fold, interval, Path(plan['requests'][str(fold)]),
                            protocol['calibration_steps'], protocol['calibration_fit_seconds'], protocol))
    write_json(output / 'cadence.json', report)
    write_json(phase / 'queue.json', {'record_dir': 'scheduler', 'jobs': jobs})
    write_json(phase / 'plan.json', calibration_plan(output, plan, report))


def replay_candidate(history, minimum, patience, delta, tolerance):
    policy = Selection(minimum, patience, delta)
    for event in history:
        observation = policy.observe(event['step'], event['validation_loss'])
        if observation['early_stop']:
            break
    before_horizon = observation['early_stop'] and observation['step'] < history[-1]['step']
    loss_gap = policy.best_loss - strict_best(history)['validation_loss']
    return {'minimum_steps': minimum, 'patience_steps': patience, 'min_delta': delta,
            'stop_step': observation['step'], 'best_step': policy.best_step,
            'best_validation_loss': policy.best_loss, 'loss_gap': loss_gap,
            'right_censored': not before_horizon, 'eligible': before_horizon and loss_gap <= tolerance}


def candidate_rank(candidate):
    return tuple(candidate[key] for key in (
        'stop_step', 'best_validation_loss', 'minimum_steps', 'patience_steps', 'min_delta',
    ))


def calibrate_fold(protocol, trace):
    history = trace['history']
    loss_scale = abs(history[0]['validation_loss'])
    tolerance = abs(trace['best_validation_loss']) * protocol['loss_tolerance_fraction']
    thresholds = [loss_scale * fraction for fraction in protocol['min_delta_fraction_of_first_loss']]
    candidates = [replay_candidate(history, minimum, patience, delta, tolerance)
                  for minimum, patience, delta in product(protocol['minimum_steps'],
                                                         protocol['patience_steps'], thresholds)]
    eligible = [candidate for candidate in candidates if candidate['eligible']]
    if eligible:
        status = 'eligible_stop_observed'
        selected = min(eligible, key=candidate_rank)
    else:
        status = 'right_censored'
        selected = replay_candidate(history, protocol['censored_fallback_minimum_steps'],
                                    protocol['censored_fallback_patience_steps'], thresholds[0], tolerance)
    return {'calibration_status': status, 'selected_candidate': selected, 'candidates': candidates,
            'loss_tolerance': tolerance, 'initial_validation_loss_magnitude': loss_scale}


def formal_recipes(protocol, fold, interval, selected, fixed):
    policy = {key: selected[key] for key in ('minimum_steps', 'patience_steps', 'min_delta')}
    return {seed: {'stage': 'formal', 'fold': fold, 'seed': seed, 'device': fixed['device'],
                   'max_steps': protocol['formal_max_steps'], 'validation_interval': interval,
                   'max_training_seconds': protocol['formal_fit_seconds'], **policy} for seed in SEEDS}


def replay_selection(output, cadence, protocol, calibrations):
    choices, recipe_files = [], {}
    calibration_plan_sha = sha256(output / 'calibration/plan.json')
    for choice in cadence['folds']:
        fold, interval = choice['fold'], choice['validation_interval']
        calibration = calibrations[fold]
        profile = next(profile for profile in choice['profiles']
                       if profile['recipe']['validation_interval'] == interval)
        match_trajectory(profile, calibration, protocol)
        result = calibrate_fold(protocol, calibration)
        formal = {}
        for seed, recipe in formal_recipes(protocol, fold, interval, result['selected_candidate'],
                                           profile['recipe']).items():
            relative = f'fold-{fold}/seed-{seed}.json'
            recipe_files[relative] = recipe
            formal[str(seed)] = {'path': relative, 'recipe_sha256': json_hash({**profile['recipe'], **recipe})}
        frozen = {'identity': {key: value for key, value in profile['identity'].items()
                               if key not in ('seed', 'stage')}, 'validation_interval': interval,
                  'cadence_reason': choice['reason'], 'formal_recipes': formal, **result}
        provenance = {'protocol_sha256': cadence['protocol_sha256'],
                      'calibration_plan_sha256': calibration_plan_sha,
                      'calibration_run_sha256': calibration['run_sha256'],
                      'profile_run_sha256': [item['run_sha256'] for item in choice['profiles']]}
        choices.append({**frozen, 'provenance': provenance,
                        'policy_sha256': json_hash({**frozen, 'provenance': provenance}),
                        'calibration': calibration})
    return {'schema': 'step04-pilot-selection-v1', 'monitor': MONITOR,
            'dataset': cadence['dataset'], 'protocol': protocol, 'plan_sha256': cadence['plan_sha256'],
            'cadence_sha256': sha256(output / 'cadence.json'),
            'calibration_plan_sha256': calibration_plan_sha, 'folds': choices}, recipe_files


def select_recipes(output, plan):
    protocol = plan['protocol']
    cadence = cadence_report(output, plan)
    if cadence != read_json(output / 'cadence.json'):
        raise ValueError('cadence evidence changed after scheduling calibration')
    if calibration_plan(output, plan, cadence) != read_json(output / 'calibration/plan.json'):
        raise ValueError('calibration plan or training sources changed after scheduling')
    grids = {choice['fold']: [choice['validation_interval']] for choice in cadence['folds']}
    traces = checked_profiles(output, 'calibration', plan, grids)
    report, recipe_files = replay_selection(output, cadence, protocol,
                                            {fold: items[0] for fold, items in traces.items()})
    selection = output / 'selection'
    selection.mkdir()
    for relative, recipe in recipe_files.items():
        write_json(selection / relative, recipe)
    write_json(selection / 'selection.json', report)


def verify_selection(output):
    selection = output / 'selection'
    report = read_json(selection / 'selection.json')
    cadence = read_json(output / 'cadence.json')
    plan = read_json(output / 'calibration/plan.json')
    protocol = report['protocol']
    validate_protocol(protocol)
    if (cadence['protocol_sha256'] != json_hash(protocol)
            or cadence['protocol_sha256'] != plan['protocol_sha256']
            or cadence['plan_sha256'] != plan['profile_plan_sha256']
            or report['cadence_sha256'] != plan['cadence_sha256']
            or [choice['fold'] for choice in cadence['folds']] != protocol['folds']):
        raise ValueError('compact pilot evidence does not match the frozen plan')
    for choice in cadence['folds']:
        if (cadence_choice(choice['fold'], *choice['profiles'], protocol) != choice
                or [trace['run_sha256'] for trace in choice['profiles']]
                != plan['profile_run_sha256'][str(choice['fold'])]):
            raise ValueError('cadence decision does not match profile evidence')
    calibrations = {choice['identity']['fold']: choice['calibration'] for choice in report['folds']}
    expected, recipe_files = replay_selection(output, cadence, protocol, calibrations)
    if report != expected or any(read_json(selection / path) != recipe for path, recipe in recipe_files.items()):
        raise ValueError('formal recipes do not match replayed pilot evidence')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime-root', type=Path, default=REPO)
    commands = parser.add_subparsers(dest='command', required=True)
    profile = commands.add_parser('profile')
    profile.add_argument('output', type=Path)
    profile.add_argument('--protocol', type=Path, default=PROTOCOL)
    profile.add_argument('--trainer', type=Path, default=HERE / 'train.py')
    source = profile.add_mutually_exclusive_group(required=True)
    source.add_argument('--tenv-release', type=Path)
    source.add_argument('--tahoe-release', type=Path)
    source.add_argument('--requests', type=Path, nargs=3)
    commands.add_parser('calibrate').add_argument('output', type=Path)
    commands.add_parser('select').add_argument('output', type=Path)
    commands.add_parser('verify').add_argument('output', type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    if args.command == 'verify':
        verify_selection(output)
        return
    if args.command != 'profile':
        plan = read_json(output / 'plan.json')
        validate_protocol(plan['protocol'])
        if args.command == 'calibrate':
            plan_calibration(output, plan)
        else:
            select_recipes(output, plan)
        return
    dataset = 'tahoe' if args.tahoe_release is not None else 'tenv'
    release = args.tahoe_release or args.tenv_release
    requests = [path.resolve() for path in args.requests] if args.requests else None
    plan_profile(args.runtime_root.resolve(), args.trainer.resolve(), output, dataset,
                 release.resolve() if release else None, requests, args.protocol)


if __name__ == '__main__':
    main()
