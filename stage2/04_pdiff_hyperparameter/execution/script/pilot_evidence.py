"""Validate completed pilot artifacts and compare paired validation traces."""
import csv
import math

from artifacts import read_json, sha256
from training_inputs import json_hash
from training_policy import MONITOR


def nonnegative(value):
    return type(value) in (int, float) and math.isfinite(value) and value >= 0


def strict_best(history):
    return min(history, key=lambda event: event['validation_loss'])


def read_profile(phase, job_id, request_path, recipe, job):
    run_path = phase / 'runs' / job_id / 'run.json'
    request = read_json(request_path)
    run = read_json(run_path)
    attempt = job['attempts'][-1]
    steps, interval = recipe['max_steps'], recipe['validation_interval']
    if (job['status'] != 'completed' or attempt['outcome'] != 'completed' or attempt['returncode'] != 0
            or run['status'] != 'completed' or run['stop_reason'] != 'max_steps'
            or run['optimizer_steps'] != steps or run['monitor'] != MONITOR):
        raise ValueError(f'{phase.name}/{job_id}: completed scheduled pilot required')
    if (run['identity'] != {**request['identity'], 'stage': 'pilot', 'seed': recipe['seed']}
            or any(run['recipe'].get(key) != value for key, value in recipe.items())
            or run['recipe_sha256'] != json_hash(run['recipe'])):
        raise ValueError(f'{phase.name}/{job_id}: request or recipe identity mismatch')

    cases = request['predictor_files']['validation_cases']
    if 'file_sha256' in request and sha256(cases) != request['file_sha256']['validation_cases']:
        raise ValueError(f'{phase.name}/{job_id}: validation cases differ from frozen request')
    with open(cases, newline='') as source:
        groups = list(dict.fromkeys(row['group_id'] for row in csv.DictReader(source)))
    batching = run['validation_batching']
    width = run['recipe']['batch_groups']
    sizes = [min(width, len(groups) - offset) for offset in range(0, len(groups), width)]
    if (not groups or batching['group_ids'] != groups or batching['group_ids_sha256'] != json_hash(groups)
            or batching['batch_groups'] != width or batching['population_size'] != run['recipe']['population_size']
            or batching['random_seed'] != recipe['seed'] + 100000
            or batching['remainder'] != 'keep and weight as one batch'):
        raise ValueError(f'{phase.name}/{job_id}: validation batching mismatch')
    history = run['history']
    if [event['step'] for event in history] != list(range(interval, steps + 1, interval)):
        raise ValueError(f'{phase.name}/{job_id}: missing or unexpected validation steps')

    best_loss, train_so_far = math.inf, 0
    for event in history:
        loss = event['validation_loss']
        batch_losses = event['batch_losses']
        if (not math.isfinite(loss) or event['batch_sizes'] != sizes or len(batch_losses) != len(sizes)
                or not all(map(math.isfinite, batch_losses))
                or not math.isclose(loss, math.fsum(batch_losses) / len(sizes), rel_tol=1e-6, abs_tol=1e-7)
                or event['checkpoint_refreshed'] != (loss < best_loss)
                or not nonnegative(event['training_seconds']) or event['training_seconds'] <= train_so_far
                or not nonnegative(event['validation_seconds']) or not nonnegative(event['checkpoint_seconds'])
                or (not event['checkpoint_refreshed'] and event['checkpoint_seconds'] != 0)):
            raise ValueError(f'{phase.name}/{job_id}: invalid validation event at {event["step"]}')
        best_loss = min(loss, best_loss)
        train_so_far = event['training_seconds']
    best = strict_best(history)
    if (run['best_step'] != best['step'] or run['best_validation_loss'] != best_loss
            or run['reload']['validation_loss'] != best_loss
            or run['reload']['prediction_max_abs_error'] != 0
            or not nonnegative(run['training_seconds']) or run['training_seconds'] != train_so_far
            or not nonnegative(run['wall_seconds']) or run['wall_seconds'] > recipe['max_training_seconds']
            or not nonnegative(run['reload_seconds'])):
        raise ValueError(f'{phase.name}/{job_id}: invalid fit, best checkpoint or reload evidence')
    return compact(run, request_path, run_path, job, attempt, sizes)


def compact(run, request_path, run_path, job, attempt, sizes):
    validation_seconds = math.fsum(event['validation_seconds'] for event in run['history'])
    checkpoint_seconds = math.fsum(event['checkpoint_seconds'] for event in run['history'])
    training_seconds = run['training_seconds']
    fields = ('step', 'validation_loss', 'training_loss', 'training_seconds',
              'validation_seconds', 'checkpoint_seconds', 'checkpoint_refreshed')
    attempt_fields = ('started_at', 'finished_at', 'memory_max_gib', 'memory_swap_max_gib',
                      'peak_host_rss_gib', 'peak_workstation_gpu_used_gib', 'outcome', 'returncode')
    evidence = {key: run[key] for key in (
        'identity', 'recipe', 'recipe_sha256', 'monitor', 'best_step', 'best_validation_loss',
        'checkpoint_sha256', 'model_state_sha256', 'reload', 'resources', 'wall_seconds',
        'reload_seconds', 'environment', 'upstream_commit',
    )}
    evidence.update(
        run_sha256=sha256(run_path), request_sha256=sha256(request_path),
        validation_batching={key: value for key, value in run['validation_batching'].items()
                             if key != 'group_ids'},
        validation_batch_sizes=sizes,
        training_seconds=training_seconds, validation_seconds=validation_seconds,
        checkpoint_seconds=checkpoint_seconds,
        overhead_fraction=(validation_seconds + checkpoint_seconds) / training_seconds,
        unallocated_wall_seconds=run['wall_seconds'] - training_seconds - validation_seconds
        - checkpoint_seconds - run['reload_seconds'],
        history=[{key: event[key] for key in fields} for event in run['history']],
        scheduler={'job': job['identity'], 'attempt_count': len(job['attempts']),
                   'latest_attempt': {key: attempt[key] for key in attempt_fields}},
    )
    return evidence


def match_trajectory(first, second, protocol):
    varied = {'validation_interval', 'max_steps', 'minimum_steps', 'patience_steps', 'max_training_seconds'}
    for key in ('identity', 'request_sha256', 'validation_batching', 'validation_batch_sizes',
                'upstream_commit', 'environment'):
        if first[key] != second[key]:
            raise ValueError(f'pilot trajectories differ in {key}')
    if any(first['recipe'].get(key) != second['recipe'].get(key)
           for key in first['recipe'].keys() | second['recipe'].keys() if key not in varied):
        raise ValueError('pilot trajectories differ in fixed training settings')
    values = {event['step']: event['validation_loss'] for event in first['history']}
    for event in second['history']:
        step = event['step']
        if step in values and not math.isclose(
                values[step], event['validation_loss'],
                rel_tol=protocol['matched_grid_loss_rtol'], abs_tol=protocol['matched_grid_loss_atol']):
            raise ValueError(f'fold {first["identity"]["fold"]}: pilot validation trajectories diverge at {step}')
