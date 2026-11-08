#!/usr/bin/env python3
"""Freeze a Step 04 dataset plan, then fit and predict through independent jobs."""
import argparse
from contextlib import contextmanager
import fcntl
import os
from pathlib import Path
import shutil
import signal
import sys
import tempfile
import time

from train import REPO, STEP03

for name in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(name, '2')
sys.dont_write_bytecode = True
sys.path.insert(0, str(REPO / STEP03 / 'execution/script'))

from artifacts import SEEDS, read_json, sha256, write_json
from pilot import queue_job, trainer_identity, training_command, verify_selection
from training_inputs import json_hash, verify_request
from training_policy import MONITOR, verify_best_checkpoint

SCRIPT = Path(__file__).resolve()
DATASETS = ('tahoe', 'tenv')
PLAN_SCHEMA = 'step04-evaluation-plan-v1'
CHECKPOINT_FIELDS = ('identity', 'recipe_sha256', 'checkpoint_sha256', 'best_step',
                     'best_validation_loss', 'model_state_sha256')


def pinned(root, relative, expected):
    path = root / relative
    if (Path(relative).is_absolute() or '..' in Path(relative).parts or path.resolve() != path
            or sha256(path) != expected):
        raise ValueError(f'frozen artifact path/hash mismatch: {relative}')
    return path


def executable(runtime):
    return runtime / 'stage2/_envs/perturbdiff/bin/python'


def launch_identity(runtime):
    files = trainer_identity(runtime, SCRIPT.with_name('train.py'))
    for path in (SCRIPT, executable(runtime), SCRIPT.with_name('pilot.py'),
                 SCRIPT.with_name('pilot_evidence.py'), SCRIPT.with_name('schedule.py'),
                 runtime / STEP03 / 'execution/input_record/verification.json'):
        files[str(path)] = sha256(path)
    return files


def use_runtime(runtime):
    if Path(sys.executable).absolute() != executable(runtime):
        raise ValueError('use the pinned isolated PerturbDiff interpreter')
    sys.path.insert(0, str(runtime / STEP03 / 'execution/script'))


def freeze(runtime, pilot_root, requests, output):
    verify_selection(pilot_root)
    selection = pilot_root / 'selection/selection.json'
    report = read_json(selection)
    dataset = report['dataset']
    if dataset not in DATASETS:
        raise ValueError('invalid formal dataset')
    use_runtime(runtime)
    from native_training import load_recipe

    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() or output.is_symlink():
        raise ValueError('formal plan requires a fresh output directory')
    with tempfile.TemporaryDirectory(dir=output.parent) as temporary:
        staging = Path(temporary)
        plan = {'schema': PLAN_SCHEMA, 'stage': 'formal', 'dataset': dataset,
                'runtime_root': str(runtime), 'launch_files': launch_identity(runtime),
                'selection_sha256': sha256(selection), 'folds': [],
                'inference': {'order': 'first appearance by population size, then group',
                              'seed_offset': 200000, 'max_groups': 8, 'max_profiles': 24}}
        shutil.copyfile(selection, staging / 'selection.json')
        fits, predictions = [], []
        sources, axes = set(), set()
        selection_folds = {choice['identity']['fold']: choice for choice in report['folds']}
        if set(selection_folds) != set(range(3)) or len(selection_folds) != len(report['folds']):
            raise ValueError('expected three frozen pilot folds')
        for fold, request_path in enumerate(requests):
            request = read_json(request_path)
            verify_request(request, runtime)
            identity = request['identity']
            choice = selection_folds[fold]
            if identity != choice['identity'] or (identity['dataset'], identity['fold']) != (dataset, fold):
                raise ValueError('request identity differs from frozen pilot fold')
            sources.add(identity['source_input_sha256'])
            axes.add(identity['gene_sha256'])
            request_name = f'requests/fold-{fold}.json'
            (staging / 'requests').mkdir(exist_ok=True)
            shutil.copyfile(request_path, staging / request_name)
            entry = {'request': request_name, 'request_sha256': sha256(staging / request_name),
                     'policy_sha256': choice['policy_sha256'], 'recipes': {}, 'recipe_sha256_by_seed': {}}
            for seed in SEEDS:
                frozen = choice['formal_recipes'][str(seed)]
                relative = f'fold-{fold}/seed-{seed}.json'
                if frozen['path'] != relative:
                    raise ValueError('unexpected pilot recipe path')
                recipe_name = f'recipes/{relative}'
                recipe_path = staging / recipe_name
                recipe_path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(selection.parent / relative, recipe_path)
                config = load_recipe(runtime, recipe_path, request)
                normalized = json_hash(config)
                if normalized != frozen['recipe_sha256'] or (config['stage'], config['seed']) != ('formal', seed):
                    raise ValueError('normalized recipe differs from frozen pilot selection')
                entry['recipes'][str(seed)] = {'path': recipe_name, 'sha256': sha256(recipe_path)}
                entry['recipe_sha256_by_seed'][str(seed)] = normalized
                fits.append((fold, seed, config['max_training_seconds'] + report['protocol']['scheduler_grace_seconds']))
                predictions.append((fold, seed, 7200 if dataset == 'tenv' else 14400))
            plan['folds'].append(entry)
        if len(sources) != 1 or len(axes) != 1:
            raise ValueError('frozen folds must share their source and gene axis')
        write_json(staging / 'plan.json', plan)
        plan_sha = sha256(staging / 'plan.json')
        for phase, slots in [('fit', fits), ('predict', predictions)]:
            jobs = [queue_job(
                runtime, f'fold-{fold}-seed-{seed}', 'training' if phase == 'fit' else 'gpu',
                timeout, 2, training_command(runtime, SCRIPT, phase, output / 'plan.json', plan_sha, fold, seed))
                for fold, seed, timeout in slots]
            write_json(staging / f'{phase}-queue.json', {'record_dir': f'{phase}-scheduler', 'jobs': jobs})
        # The sealed plan is published together with the exact worker arguments.
        ready = staging / 'ready'
        ready.mkdir()
        for name in ('plan.json', 'selection.json', 'requests', 'recipes', 'fit-queue.json', 'predict-queue.json'):
            (staging / name).rename(ready / name)
        ready.rename(output)
    print(plan_sha, flush=True)


def read_plan(path, expected):
    if sha256(path) != expected:
        raise ValueError('frozen plan hash mismatch')
    plan = read_json(path)
    if (plan['schema'] != PLAN_SCHEMA or plan['stage'] != 'formal'
            or plan['dataset'] not in DATASETS or len(plan['folds']) != 3):
        raise ValueError('expected frozen formal dataset plan')
    return plan


def slot_inputs(path, plan, fold, seed):
    runtime = Path(plan['runtime_root'])
    use_runtime(runtime)
    if plan['launch_files'] != launch_identity(runtime):
        raise ValueError('formal runtime or producer identity changed')
    root = path.parent
    pinned(root, 'selection.json', plan['selection_sha256'])
    entry = plan['folds'][fold]
    request = read_json(pinned(root, entry['request'], entry['request_sha256']))
    verify_request(request, runtime)
    identity = {**request['identity'], 'seed': seed, 'stage': 'formal'}
    if (identity['dataset'], identity['fold']) != (plan['dataset'], fold):
        raise ValueError('formal slot request identity mismatch')
    recipe_entry = entry['recipes'][str(seed)]
    recipe = pinned(root, recipe_entry['path'], recipe_entry['sha256'])
    from native_training import load_recipe
    config = load_recipe(runtime, recipe, request)
    recipe_hash = entry['recipe_sha256_by_seed'][str(seed)]
    if (config['stage'], config['seed'], json_hash(config)) != ('formal', seed, recipe_hash):
        raise ValueError('formal slot recipe identity mismatch')
    return runtime, request, identity, recipe, config


@contextmanager
def slot_lock(slot):
    slot.mkdir(parents=True, exist_ok=True)
    with (slot / 'worker.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield


def completed_fit(run, identity, config):
    metadata = read_json(run / 'run.json')
    if (metadata['status'] != 'completed' or metadata['monitor'] != MONITOR
            or metadata['identity'] != identity or metadata['recipe'] != config
            or metadata['recipe_sha256'] != json_hash(config)):
        raise ValueError('completed fit differs from formal identity or recipe')
    best = verify_best_checkpoint(metadata, 'fit')
    if (metadata['reload']['validation_loss'] != best['validation_loss']
            or metadata['reload']['prediction_max_abs_error'] != 0):
        raise ValueError('completed fit lacks an exact EMA replay')
    return metadata


def fit_slot(slot, plan_sha, runtime, request, identity, recipe, config):
    from native_training import fit

    destination = slot / 'fit'
    if not destination.exists():
        attempt = Path(tempfile.mkdtemp(prefix='fit-attempt-', dir=slot))
        run = attempt / 'run'
        print(f'Fit attempt evidence: {run}', flush=True)
        fit(runtime, request, recipe, run)
        completed_fit(run, identity, config)
        write_json(run / 'fit.json', {'plan_sha256': plan_sha, 'run_sha256': sha256(run / 'run.json')})
        run.rename(destination)
        attempt.rmdir()
    verify_fit_seal(destination, plan_sha)
    metadata = completed_fit(destination, identity, config)
    pinned(destination, 'best.pt', metadata['checkpoint_sha256'])


def verify_fit_seal(directory, plan_sha):
    seal = read_json(directory / 'fit.json')
    if seal['plan_sha256'] != plan_sha:
        raise ValueError('fit did not run against this frozen plan')
    pinned(directory, 'run.json', seal['run_sha256'])


def predict_slot(slot, plan, plan_sha, request, identity, config):
    from native_training import (JointPerturbDiff, PredictorInputs, checkpoint_reload,
                                 groups_for, np, resources, seed_randomness, torch,
                                 validation_random_stream)

    fit_dir, destination = slot / 'fit', slot / 'prediction'
    verify_fit_seal(fit_dir, plan_sha)
    metadata = completed_fit(fit_dir, identity, config)
    if destination.exists():
        manifest = read_json(destination / 'prediction.json')
        if (manifest['plan_sha256'] != plan_sha
                or any(manifest[key] != metadata[key] for key in CHECKPOINT_FIELDS)):
            raise ValueError('existing prediction differs from frozen slot')
        for name, key in (('run.json', 'run_sha256'), ('best.pt', 'checkpoint_sha256'),
                          ('prediction.npy', 'prediction_sha256'), ('cases.json', 'cases_sha256'),
                          ('replay.json', 'replay_sha256')):
            pinned(destination, name, manifest[key])
        return
    staging = Path(tempfile.mkdtemp(prefix='predict-attempt-', dir=slot))
    print(f'Prediction attempt evidence: {staging}', flush=True)

    def stop(_signal, _frame):
        raise InterruptedError('held-out prediction interrupted')

    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    started = time.monotonic()
    torch.set_num_threads(config['threads'])
    seed_randomness(config['seed'])
    device = torch.device(config['device'])
    torch.cuda.reset_peak_memory_stats(device)
    inputs = PredictorInputs(request)
    model = JointPerturbDiff(inputs, config, staging).to(device)
    model.native.to(device)
    validation = groups_for(inputs.cases['validation'])
    replay = checkpoint_reload(fit_dir, metadata, model, inputs, validation, config, device, lambda: None)
    if replay != metadata['reload']:
        raise ValueError('held-out predictor did not reproduce the validation-selected EMA')
    write_json(staging / 'replay.json', {**replay, **{key: metadata[key] for key in CHECKPOINT_FIELDS},
                                        'plan_sha256': plan_sha})
    cases = inputs.cases['test']
    groups = groups_for(cases)
    buckets = {}
    for group in groups:
        buckets.setdefault(len(group), []).append(group)
    inference = plan['inference']
    seen = np.zeros(len(cases), dtype=bool)
    predictions = np.lib.format.open_memmap(staging / 'prediction.npy', mode='w+', dtype='float32',
                                           shape=(len(cases), inputs.n_genes))
    inference_seed = config['seed'] + inference['seed_offset']
    generator = torch.Generator(device=device).manual_seed(inference_seed)
    with torch.no_grad(), validation_random_stream(inference_seed, device):
        for size, population_groups in buckets.items():
            batch_size = min(inference['max_groups'], max(1, inference['max_profiles'] // size))
            for offset in range(0, len(population_groups), batch_size):
                chosen = population_groups[offset:offset + batch_size]
                batch, indices = model.batch(inputs, 'test', chosen, device)
                if 'pert_emb' in batch or not np.array_equal(indices, chosen):
                    raise ValueError('held-out predictor changed frozen case membership')
                noise = torch.randn(batch['cont_emb'].shape, device=device, generator=generator)
                delta = (model.sample(batch, config, noise) - batch['cont_emb']) * config['expression_scale']
                values = delta.cpu().numpy()
                indices = indices.reshape(-1)
                if (values.shape != (len(chosen), size, inputs.n_genes) or values.dtype != np.float32
                        or not np.isfinite(values).all() or seen[indices].any()):
                    raise ValueError('invalid or repeated held-out prediction rows')
                predictions[indices] = values.reshape(-1, inputs.n_genes)
                seen[indices] = True
    if not seen.all():
        raise ValueError('held-out predictions did not cover every frozen case')
    predictions.flush()
    del predictions
    write_json(staging / 'cases.json', [[case['case_id'], case['group_id']] for case in cases])
    for name in ('run.json', 'best.pt'):
        os.link(fit_dir / name, staging / name)
    write_json(staging / 'prediction.json', {
        'schema': 'step04-heldout-prediction-v1', 'plan_sha256': plan_sha,
        **{key: metadata[key] for key in CHECKPOINT_FIELDS},
        'run_sha256': sha256(staging / 'run.json'), 'cases_sha256': sha256(staging / 'cases.json'),
        'prediction_sha256': sha256(staging / 'prediction.npy'), 'replay_sha256': sha256(staging / 'replay.json'),
        'case_count': len(cases), 'group_count': len(groups), 'maximum_population': max(buckets),
        'wall_seconds': time.monotonic() - started, 'resources': resources(device),
    })
    staging.rename(destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime-root', type=Path, default=REPO)
    commands = parser.add_subparsers(dest='command', required=True)
    freeze_cmd = commands.add_parser('freeze')
    freeze_cmd.add_argument('pilot_record', type=Path)
    freeze_cmd.add_argument('output', type=Path)
    freeze_cmd.add_argument('requests', type=Path, nargs=3, metavar='REQUEST')
    for phase in ('fit', 'predict'):
        worker = commands.add_parser(phase)
        worker.add_argument('plan', type=Path)
        worker.add_argument('plan_sha256')
        worker.add_argument('fold', type=int, choices=range(3))
        worker.add_argument('seed', type=int, choices=SEEDS)
    args = parser.parse_args()
    try:
        if args.command == 'freeze':
            freeze(args.runtime_root.resolve(), args.pilot_record.resolve(), args.requests,
                   args.output.resolve())
        else:
            path = args.plan.resolve()
            plan = read_plan(path, args.plan_sha256)
            runtime, request, identity, recipe, config = slot_inputs(path, plan, args.fold, args.seed)
            slot = path.parent / 'slots' / f'fold-{args.fold}-seed-{args.seed}'
            with slot_lock(slot):
                if args.command == 'fit':
                    fit_slot(slot, args.plan_sha256, runtime, request, identity, recipe, config)
                else:
                    predict_slot(slot, plan, args.plan_sha256, request, identity, config)
    except (KeyError, OSError, TypeError, ValueError) as error:
        parser.exit(1, f'formal {args.command} rejected: {error}\n')


if __name__ == '__main__':
    main()
