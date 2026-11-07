#!/usr/bin/env python3
"""Collect terminal scheduler outcomes for one frozen Step 04 dataset."""
import argparse
import errno
import os
from pathlib import Path
import shutil
import tempfile

from schedule import load_queue, load_status
from formal import CHECKPOINT_FIELDS, SEEDS, pinned, read_json, read_plan, sha256, write_json

FAILURES = {'failed', 'timeout', 'interrupted', 'cancelled'}


def queue_status(root, phase, plan_sha):
    queue = root / f'{phase}-queue.json'
    record, jobs = load_queue(queue)
    expected = {f'fold-{fold}-seed-{seed}': [phase, str(root / 'plan.json'), plan_sha, str(fold), str(seed)]
                for fold in range(3) for seed in SEEDS}
    if (record != root / f'{phase}-scheduler'
            or {job['id']: job['command'][-5:] for job in jobs} != expected):
        raise ValueError(f'{phase} queue does not match the frozen plan')
    status_file = record / 'status.json'
    status_hash = sha256(status_file) if status_file.is_file() else None
    status = load_status(record, jobs)
    if status_hash is not None:
        pinned(record, 'status.json', status_hash)
    evidence = {'queue': queue.name, 'queue_sha256': sha256(queue),
                'status': str(status_file.relative_to(root)) if status_hash else None,
                'status_sha256': status_hash}
    return status['jobs'], evidence


def failure(root, row, phase, slot):
    status = row['status']
    if status not in FAILURES:
        raise ValueError(f'{slot.name} {phase} is {status}. Wait for the scheduler to finish.')
    attempt = row['attempts'][-1]
    reason = attempt['reason']
    log = Path(attempt['log'])
    if log.parent != root / f'{phase}-scheduler' or log.resolve() != log:
        raise ValueError('scheduler log is outside the formal evidence root')
    if phase == 'fit' and log.is_file():
        with log.open(errors='replace') as stream:
            for line in stream:
                if line.startswith('Fit attempt evidence: '):
                    run = Path(line.removeprefix('Fit attempt evidence: ').strip()) / 'run.json'
                    if run.is_relative_to(slot) and run.resolve() == run and run.is_file():
                        metadata = read_json(run)
                        if metadata['status'] in ('failed', 'incomplete'):
                            if status == 'failed':
                                status = metadata['status']
                            detail = metadata.get('error') or metadata.get('stop_reason')
                            if detail:
                                reason = f'{reason}, {detail}'
                    break
    return {'status': status, 'phase': phase, 'reason': reason,
            'scheduler_attempts': len(row['attempts']), 'worker_log': str(log.relative_to(root))}


def collect(plan_path, plan_sha, output):
    plan = read_plan(plan_path, plan_sha)
    root = plan_path.parent
    if output.parent != root or output.exists() or output.is_symlink():
        raise ValueError('outcomes require a fresh file alongside the frozen plan')
    fits, fit_evidence = queue_status(root, 'fit', plan_sha)
    predictions, prediction_evidence = queue_status(root, 'predict', plan_sha)
    outcomes = []
    for fold, entry in enumerate(plan['folds']):
        request = read_json(pinned(root, entry['request'], entry['request_sha256']))
        for seed in SEEDS:
            job = f'fold-{fold}-seed-{seed}'
            slot = root / 'slots' / job
            fit, prediction = fits[job], predictions[job]
            outcome = {'dataset': plan['dataset'], 'fold': fold, 'seed': seed}
            if fit['status'] != 'completed':
                outcome.update(failure(root, fit, 'fit', slot))
            elif prediction['status'] != 'completed':
                outcome.update(failure(root, prediction, 'predict', slot))
            else:
                path = slot / 'prediction/prediction.json'
                manifest = read_json(path)
                run = read_json(pinned(path.parent, 'run.json', manifest['run_sha256']))
                if (manifest['schema'] != 'step04-heldout-prediction-v1'
                        or manifest['plan_sha256'] != plan_sha
                        or manifest['identity'] != {**request['identity'], 'seed': seed, 'stage': 'formal'}
                        or manifest['recipe_sha256'] != entry['recipe_sha256_by_seed'][str(seed)]
                        or any(manifest[key] != run[key] for key in CHECKPOINT_FIELDS)):
                    raise ValueError(f'{job} prediction identity mismatch')
                outcome.update(status='completed', prediction=str(path.relative_to(root)),
                               prediction_manifest_sha256=sha256(path))
            outcomes.append(outcome)
    write_json(output, {'schema': 'step04-evaluation-outcomes-v1', 'plan_sha256': plan_sha,
                        'scheduler': {'fit': fit_evidence, 'predict': prediction_evidence},
                        'outcomes': outcomes})


def publish_checkpoints(plan_path, plan_sha, outcomes_path, evaluation, destination):
    plan = read_plan(plan_path, plan_sha)
    root = plan_path.parent
    if outcomes_path.parent != root:
        raise ValueError('checkpoint outcomes must be alongside the frozen plan')
    outcomes_sha = sha256(outcomes_path)
    outcomes = read_json(outcomes_path)
    evaluation_path = evaluation / 'evaluation.json'
    audit = read_json(evaluation_path)
    expected = {(fold, seed) for fold in range(3) for seed in SEEDS}
    slots = {(row['fold'], row['seed']): row for row in outcomes['outcomes']}
    if (plan['dataset'] != 'tenv' or outcomes['schema'] != 'step04-evaluation-outcomes-v1'
            or outcomes['plan_sha256'] != plan_sha or audit['schema'] != 'step04-evaluation-v1'
            or audit['stage'] != 'formal' or audit['datasets'] != ['tenv']
            or not audit['tenv_test_memberships_previously_scored_in_step03']
            or audit['complete'] is not True or audit['plan_sha256'] != plan_sha
            or audit['outcomes_sha256'] != outcomes_sha or set(slots) != expected
            or len(outcomes['outcomes']) != len(expected)
            or any(row['dataset'] != 'tenv' or row['status'] != 'completed' for row in slots.values())):
        raise ValueError('checkpoint release requires a complete evaluated tenv matrix')
    for path, digest in audit['files'].items():
        pinned(evaluation, path, digest)
    if destination.exists() or destination.is_symlink():
        raise ValueError('checkpoint release requires a fresh output directory')
    destination.parent.mkdir(parents=True, exist_ok=True)
    index = {'schema': 'step04-tenv-checkpoints-v1', 'plan_sha256': plan_sha,
             'selection_sha256': plan['selection_sha256'], 'outcomes_sha256': outcomes_sha,
             'evaluation_sha256': sha256(evaluation_path), 'checkpoints': []}
    with tempfile.TemporaryDirectory(dir=destination.parent) as temporary:
        release = Path(temporary) / 'tenv'
        release.mkdir()
        for (fold, seed), outcome in sorted(slots.items()):
            entry = plan['folds'][fold]
            request_path = pinned(root, entry['request'], entry['request_sha256'])
            request = read_json(request_path)
            manifest_path = pinned(root, outcome['prediction'], outcome['prediction_manifest_sha256'])
            manifest = read_json(manifest_path)
            artifact = manifest_path.parent
            run_path = pinned(artifact, 'run.json', manifest['run_sha256'])
            run = read_json(run_path)
            replay_path = pinned(artifact, 'replay.json', manifest['replay_sha256'])
            replay = read_json(replay_path)
            checkpoint = {key: run[key] for key in CHECKPOINT_FIELDS}
            if (manifest['schema'] != 'step04-heldout-prediction-v1' or run['status'] != 'completed'
                    or manifest['plan_sha256'] != plan_sha
                    or any(manifest[key] != value for key, value in checkpoint.items())
                    or checkpoint['identity'] != {**request['identity'], 'seed': seed, 'stage': 'formal'}
                    or checkpoint['recipe_sha256'] != entry['recipe_sha256_by_seed'][str(seed)]
                    or replay != {**run['reload'], **checkpoint, 'plan_sha256': plan_sha}
                    or replay['validation_loss'] != checkpoint['best_validation_loss']
                    or replay['prediction_max_abs_error'] != 0):
                raise ValueError(f'checkpoint or exact EMA replay mismatch for fold {fold} seed {seed}')
            path = release / f'fold-{fold}' / f'seed-{seed}'
            path.mkdir(parents=True)
            source = pinned(artifact, 'best.pt', checkpoint['checkpoint_sha256'])
            try:
                os.link(source, path / 'best.pt')
            except OSError as error:
                if error.errno == errno.EXDEV:
                    raise ValueError('checkpoint release must share the formal artifacts filesystem') from error
                raise
            recipe_entry = entry['recipes'][str(seed)]
            recipe = pinned(root, recipe_entry['path'], recipe_entry['sha256'])
            for item, name in [(request_path, 'request.json'), (recipe, 'recipe.json'),
                               (run_path, 'run.json'), (replay_path, 'replay.json')]:
                shutil.copyfile(item, path / name)
            index['checkpoints'].append({
                **checkpoint, 'directory': str(path.relative_to(release)),
                'policy_sha256': entry['policy_sha256'], 'request_sha256': entry['request_sha256'],
                'raw_recipe_sha256': recipe_entry['sha256'], 'run_sha256': manifest['run_sha256'],
                'replay_sha256': manifest['replay_sha256'], 'upstream_commit': run['upstream_commit'],
                'environment': run['environment'],
            })
        write_json(release / 'checkpoints.json', index)
        release.rename(destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('plan', type=Path)
    parser.add_argument('plan_sha256')
    parser.add_argument('output', type=Path)
    parser.add_argument('--publish-checkpoints', type=Path, nargs=2, metavar=('EVALUATION', 'DESTINATION'),
                        help='publish all nine tenv checkpoints from existing outcomes and evaluation')
    args = parser.parse_args()
    try:
        plan, output = args.plan.resolve(), args.output.absolute()
        if args.publish_checkpoints:
            evaluation, destination = args.publish_checkpoints
            publish_checkpoints(plan, args.plan_sha256, output, evaluation.resolve(), destination.absolute())
        else:
            collect(plan, args.plan_sha256, output)
    except (KeyError, OSError, TypeError, ValueError) as error:
        parser.exit(1, f'formal outcomes rejected: {error}\n')


if __name__ == '__main__':
    main()
