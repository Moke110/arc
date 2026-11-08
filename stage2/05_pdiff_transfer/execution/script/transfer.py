#!/usr/bin/env python3
"""Freeze Step 05 Transfer slots and use a shared Tahoe bridge for both arms."""
import argparse
from contextlib import nullcontext
from pathlib import Path
import tempfile
from unittest.mock import patch

import cold_tahoe as cold
import tahoe_zero_shot as tahoe
import warm_tenv as warm
from artifacts import SEEDS, read_json, rows, sha256, write_json
from evaluate import grouped_cases, load_schema, prediction
from formal import (CHECKPOINT_FIELDS, completed_fit, fit_slot, pinned, predict_slot, slot_lock,
                    use_runtime, verify_fit_seal)
from tahoe_bridge import DOSE_ROWS, LIMITATION, load_initialization, request_digest, state_digest
from training_inputs import STEP03, json_hash

REPO, STEP04, STEP05 = cold.REPO, cold.STEP04, cold.STEP05
REGIMES = (*tahoe.REGIMES, 'cold-drug')
SCHEMA = 'step05-transfer-plan-v1'
RECORD = REPO / STEP05 / 'execution'


def records(runtime):
    output = runtime / STEP05 / 'execution/output'
    historical = read_json(RECORD / 'tahoe_zero_shot_record/qualification.json')
    cold_record = read_json(RECORD / 'cold_tahoe_record/qualification.json')['local_artifact_sha256']
    return output, historical, cold_record


def recipient_slot(output, regime):
    warm_start = regime == 'warm-start'
    name, plan, record = ('warm-tenv-ticket77', 'plan.json', 'warm') if warm_start else (
        'cold-tenv-ticket76', 'split.json', 'cold')
    hashes = read_json(RECORD / f'{record}_record/replay.json')['local_artifact_sha256']
    return output / name, hashes[f'{name}/{plan}']


def recipes(runtime, stage):
    base = read_json(tahoe.historical_record('tenv') / 'plan.json')
    config = STEP03 / 'execution/models/perturbdiff/v2.config.json'
    digest, = (value for path, value in base['launch_files'].items() if path.endswith(str(config)))
    fixed = read_json(pinned(runtime, str(config), digest))
    fixed_keys = ('batch_groups', 'population_size', 'chemical_dimension', 'learning_rate', 'weight_decay',
                  'warmup_steps', 'sampling_start_time', 'threads', 'expression_scale',
                  'schedule_training_steps', 'ema_decay', 'ema_update_steps')
    defaults = {key: fixed[key] for key in fixed_keys}
    root = runtime / STEP04 / 'execution/output/tenv-formal-ticket72'
    fold = base['folds'][0]
    formal = {}
    for seed in map(str, SEEDS):
        entry = fold['recipes'][seed]
        raw = read_json(pinned(root, entry['path'], entry['sha256']))
        config = {**defaults, **raw}
        digest = fold['recipe_sha256_by_seed'][seed]
        if json_hash(config) != digest:
            raise ValueError('target recipe differs from the matching Step 04 CM recipe')
        formal[seed] = {'recipe': raw, 'config': config, 'recipe_sha256': digest}
    if stage == 'formal':
        return formal, formal
    smoke = {**formal['17']['recipe'], **warm.SMOKE}
    config = {**defaults, **smoke}
    return formal, {'17': {'recipe': smoke, 'config': config, 'recipe_sha256': json_hash(config)}}


def manifest(runtime, stage, selected=REGIMES):
    output, historical, cold_record = records(runtime)
    historical_root = output / 'tahoe-zero-shot-ticket78'
    source_sha = historical['source_artifact']['tahoe-zero-shot-ticket78/source.json']
    source = load_schema(pinned(historical_root, 'source.json', source_sha), tahoe.SOURCE_SCHEMA)
    cold_root = output / 'cold-tahoe-ticket79'
    release_sha = cold_record['cold-tahoe-ticket79/release.json']
    release = read_json(pinned(cold_root, 'release.json', release_sha))
    formal, active = recipes(runtime, stage)
    targets, slots, context = {}, [], {}
    for regime in selected:
        target_root, target_sha = recipient_slot(output, regime)
        request, inference, fold, warm_source = tahoe.recipient(runtime, target_root, regime, target_sha)
        tahoe.verify_pair(source, request)
        targets[regime] = {'split_sha256': target_sha, 'request_sha256': request_digest(request),
                           'identity': request['identity'], 'inference': inference}
        for seed in SEEDS:
            if regime == 'cold-drug':
                lineage = {'source_release_sha256': release_sha,
                           'source_request_sha256': cold_record['cold-tahoe-ticket79/request.json'],
                           'source_recipe_sha256': release['recipes'][str(seed)]['resolved_sha256'],
                           'source_fit': f'source/seed-{seed}/fit'}
            else:
                item = source['checkpoints'][str(seed)]
                lineage = {'source_sha256': source_sha, 'source_request_sha256': source['request_file_sha256'],
                           'source_checkpoint_sha256': item['checkpoint_sha256'],
                           'source_recipe_sha256': item['recipe_sha256']}
            slots.append({'regime': regime, 'fold': 0, 'arm': 'Transfer', 'seed': seed, 'stage': 'formal',
                          'recipient_request_sha256': targets[regime]['request_sha256'],
                          'recipe_sha256': formal[str(seed)]['recipe_sha256'], **lineage})
        context[regime] = request, fold, warm_source
    plan = {'schema': SCHEMA, 'stage': stage, 'formal_slots': slots, 'targets': targets,
            'formal_recipes': formal, 'active_recipes': active}
    return plan, context, source, source_sha


def verify(runtime, root, digest, regime=None):
    plan = load_schema(pinned(root, 'plan.json', digest), SCHEMA)
    if plan['stage'] not in ('formal', 'smoke') or set(plan['targets']) != set(REGIMES):
        raise ValueError('unexpected Transfer stage or regime inventory')
    expected, context, source, source_sha = manifest(runtime, plan['stage'], (regime,) if regime else REGIMES)
    subset = plan if regime is None else {
        **plan, 'targets': {regime: plan['targets'][regime]},
        'formal_slots': [slot for slot in plan['formal_slots'] if slot['regime'] == regime],
    }
    if subset != expected:
        raise ValueError('Transfer plan differs from frozen requests, recipes or source slots')
    return plan, context, source, source_sha


def prepare(runtime, root, stage):
    if root.exists() or root.is_symlink():
        raise ValueError('prepare requires a fresh output directory')
    plan, _, _, _ = manifest(runtime, stage)
    root.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=root.parent) as temporary:
        ready = Path(temporary) / 'ready'
        write_json(ready / 'plan.json', plan)
        ready.rename(root)
    print(f'Frozen Transfer plan SHA-256 {sha256(root / "plan.json")}', flush=True)


def cold_inputs(runtime):
    output, _, evidence = records(runtime)
    release = output / 'cold-tahoe-ticket79'
    formal = runtime / STEP04 / 'execution/output/tahoe-formal-ticket73'
    split = output / 'cold-tenv-ticket76/split.json'
    request = cold.verify_release(runtime, formal, split, release)
    pinned(release, 'request.json', evidence['cold-tahoe-ticket79/request.json'])
    return output, release, request, evidence


def cold_source(runtime, root, plan, digest, seed, inputs, fit=False):
    output, release, request, evidence = inputs
    stage = plan['stage']
    identity = {**request['identity'], 'stage': stage, 'seed': seed}
    if stage == 'formal':
        slot = root / 'source' / f'seed-{seed}'
        directory = slot / 'fit'
        entry = next(slot for slot in plan['formal_slots'] if slot['regime'] == 'cold-drug' and slot['seed'] == seed)
        recipe = release / f'recipes/seed-{seed}.json'
        config_sha = entry['source_recipe_sha256']
        if fit:
            use_runtime(runtime)
            from native_training import load_recipe

            config = load_recipe(runtime, recipe, request)
            if json_hash(config) != config_sha:
                raise ValueError('cold Tahoe training recipe differs from the frozen source recipe')
            with slot_lock(slot):
                fit_slot(slot, digest, runtime, request, identity, recipe, config)
        verify_fit_seal(directory, digest)
    else:
        if fit or seed != 17:
            raise ValueError('qualification uses the previously sealed seed-17 cold Tahoe source')
        directory = output / 'cold-tahoe-smoke-ticket79'
        pinned(directory, 'run.json', evidence[f'{directory.name}/run.json'])
        replay_name = 'cold-tahoe-smoke-replay-ticket79.json'
        replay = read_json(pinned(output, replay_name, evidence[replay_name]))
        config_sha = read_json(directory / 'run.json')['recipe_sha256']
        recipe = cold.SMOKE_RECIPE
    run = read_json(directory / 'run.json')
    completed_fit(directory, identity, run['recipe'])
    validation_groups = list(dict.fromkeys(row['group_id'] for row in rows(request['predictor_files']['validation_cases'])))
    if (run['recipe_sha256'] != config_sha
            or any(run['recipe'][key] != value for key, value in read_json(recipe).items())
            or run['validation_batching']['group_ids_sha256'] != json_hash(validation_groups)):
        raise ValueError('cold Tahoe source did not use its frozen recipe and validation groups')
    if stage == 'smoke' and replay != {'identity': identity, 'best_step': run['best_step'],
                                      'checkpoint_sha256': run['checkpoint_sha256'], **run['reload']}:
        raise ValueError('cold Tahoe qualification replay differs from its selected source state')
    checkpoint = pinned(directory, 'best.pt', run['checkpoint_sha256'])
    cases = tahoe.role_cases(request)
    parents = tahoe.parent_blocks(request, cases)
    doses, indices = tahoe.dose_layout(tahoe.key_layout(cases))
    if doses != sorted({float(row['dose_um']) for row in cases['train']}):
        raise ValueError('cold Tahoe training does not cover the source dose vocabulary')
    selected = {key: run[key] for key in CHECKPOINT_FIELDS}
    source = {
        'request_identity': request['identity'], 'request_file_sha256': sha256(release / 'request.json'),
        'release_sha256': sha256(release / 'release.json'),
        'roles': {role: {'parent_connectivity': parents[role], 'case_count': len(cases[role]),
                         'group_count': len(grouped_cases(cases[role])[0]),
                         'cases_sha256': request['file_sha256'][f'{role}_cases']} for role in cold.ROLES},
        'checkpoints': {str(seed): {**selected, 'path': str(checkpoint), 'doses_um': doses, 'dose_indices': indices,
                                   'producer': {'run_sha256': sha256(directory / 'run.json'),
                                                'raw_recipe_sha256': sha256(recipe),
                                                'request_sha256': sha256(release / 'request.json'),
                                                'replay_sha256': json_hash(run['reload'])}}},
    }
    return source, json_hash(source)


def tahoe_directory(runtime, root, regime, seed):
    output, historical, _ = records(runtime)
    relative = f'tahoe-zero-shot-ticket78/{regime}/seed-{seed}/prediction.json'
    digest = historical['qualification'].get(regime, {}).get('artifacts', {}).get(relative)
    if digest is not None:
        return pinned(output, relative, digest).parent
    return root / 'tahoe' / regime / f'seed-{seed}'


def verify_tahoe(runtime, root, regime, seed, request, target, source, source_sha, *, load_values=True):
    directory = tahoe_directory(runtime, root, regime, seed)
    manifest, values = tahoe.verify_prediction(directory, source, source_sha, request, regime,
                                               target['split_sha256'], seed, load_values=load_values)
    return directory, manifest, values


def transfer_lineage(source, source_sha, target, item, bridge, seed):
    checkpoint = source['checkpoints'][str(seed)]
    return {'source_sha256': source_sha, 'source_request_sha256': source['request_file_sha256'],
            'source_checkpoint_sha256': checkpoint['checkpoint_sha256'],
            'source_recipe_sha256': checkpoint['recipe_sha256'],
            'source_history_sha256': checkpoint['producer']['run_sha256'],
            'recipient_request_sha256': target['request_sha256'], 'recipient_split_sha256': target['split_sha256'],
            'recipe_sha256': item['recipe_sha256'], 'bridge_sha256': bridge['bridge_sha256'],
            'mapped_state_sha256': bridge['mapped_state_sha256'], 'dose_limitation': LIMITATION}


def verify_transfer_fit(slot, digest, identity, item, lineage):
    directory = slot / 'fit'
    verify_fit_seal(directory, digest)
    run = completed_fit(directory, identity, item['config'])
    if (run['initialization'] != lineage or run['training_initial_state'] != {
            'optimizer_state_entries': 0, 'scheduler_start_step': 0, 'ema_start_step': 0,
            'mapped_dose_rows_trainable': True, 'mapped_dose_gradient_observed': True}):
        raise ValueError('Transfer fit did not consume the frozen Tahoe bridge with fresh training state')
    pinned(directory, 'best.pt', run['checkpoint_sha256'])
    return run


def transfer_worker(runtime, root, plan, digest, regime, seed, command, request, directory, lineage):
    use_runtime(runtime)
    import native_training as native
    from warm_inputs import WarmInputs

    item = plan['active_recipes'][str(seed)]
    worker_request = {**request, 'identity': transfer_identity(request, regime)}
    identity = {**worker_request['identity'], 'seed': seed, 'stage': plan['stage']}
    slot = root / 'slots' / regime / f'seed-{seed}'
    inputs = patch.object(native, 'PredictorInputs', WarmInputs) if regime == 'warm-start' else nullcontext()
    with slot_lock(slot), inputs:
        if command == 'fit':
            fit_transfer(slot, digest, runtime, worker_request, request, identity, item, directory, lineage, native)
        else:
            verify_transfer_fit(slot, digest, identity, item, lineage)
            predict_slot(slot, plan['targets'][regime], digest, worker_request, identity, item['config'])


def transfer_identity(request, regime):
    return {**request['identity'], 'arm': 'Transfer', 'regime': regime}


def fit_transfer(slot, digest, runtime, request, recipient, identity, item, directory, lineage, native):
    native_fit = native.fit
    initial = {}

    class InitializedModel(native.JointPerturbDiff):
        def __init__(self, inputs, config, directory_out):
            super().__init__(inputs, config, directory_out)
            bridge = load_initialization(self, directory, lineage['bridge_sha256'], recipient)
            if bridge['mapped_state_sha256'] != lineage['mapped_state_sha256']:
                raise ValueError('Transfer initialization differs from Tahoe prediction')
            parameter = self.get_parameter(DOSE_ROWS)
            if not parameter.requires_grad:
                raise ValueError('mapped dose embeddings must be trainable')
            initial['mapped_dose_rows_trainable'] = True

            def observe(gradient):
                if gradient[1:].count_nonzero().item():
                    initial['mapped_dose_gradient_observed'] = True
                    hook.remove()

            hook = parameter.register_hook(observe)

    native_optimizer = native.get_optimizer
    native_schedule = native.get_cosine_schedule_with_warmup

    def optimizer(*args, **kwargs):
        result = native_optimizer(*args, **kwargs)
        initial['optimizer_state_entries'] = len(result.state)
        return result

    def schedule(*args, **kwargs):
        result = native_schedule(*args, **kwargs)
        initial['scheduler_start_step'] = result.last_epoch
        return result

    class FreshEMA(native.EMAWeightAveraging):
        def setup(self, trainer, model, stage):
            super().setup(trainer, model, stage)
            if (self._average_model.n_averaged.item()
                    or state_digest(self._average_model.module.state_dict()) != state_digest(model.state_dict())):
                raise ValueError('fresh EMA must start at the mapped model state')
            initial['ema_start_step'] = trainer.global_step

    def fit(runtime, request, recipe, run_dir):
        try:
            with patch.multiple(native, JointPerturbDiff=InitializedModel, get_optimizer=optimizer,
                                get_cosine_schedule_with_warmup=schedule, EMAWeightAveraging=FreshEMA):
                native_fit(runtime, request, recipe, run_dir)
            if not initial.get('mapped_dose_gradient_observed'):
                raise ValueError('no training gradient reached the mapped positive dose rows')
        finally:
            run_file = run_dir / 'run.json'
            if run_file.exists():
                run = read_json(run_file)
                write_json(run_file, {**run, 'initialization': lineage, 'training_initial_state': initial})

    with tempfile.TemporaryDirectory(dir=slot) as temporary:
        recipe = Path(temporary) / 'recipe.json'
        write_json(recipe, item['recipe'])
        if native.load_recipe(runtime, recipe, request) != item['config']:
            raise ValueError('Transfer training config differs from the frozen recipe')
        with patch.object(native, 'fit', fit):
            fit_slot(slot, digest, runtime, request, identity, recipe, item['config'])
    verify_transfer_fit(slot, digest, identity, item, lineage)


def score(runtime, root, plan, digest, regime, output, context, sources, arm):
    request, fold, warm_source = context[regime]
    target = plan['targets'][regime]
    cases = rows(request['predictor_files']['test_cases'])
    groups, _ = grouped_cases(cases)
    evidence, values = {}, {}
    for seed in map(int, plan['active_recipes']):
        source, source_sha = sources(seed)
        directory, bridge, zero_shot = verify_tahoe(runtime, root, regime, seed, request, target,
                                                    source, source_sha, load_values=arm == 'Tahoe')
        item = plan['active_recipes'][str(seed)]
        lineage = transfer_lineage(source, source_sha, target, item, bridge, seed)
        overlap = bridge['source_parent_overlap']
        entry = {'source_sha256': source_sha, 'source_parent_overlap': overlap,
                 'bridge_sha256': lineage['bridge_sha256'],
                 'mapped_state_sha256': lineage['mapped_state_sha256']}
        if arm == 'Tahoe':
            manifest = directory / 'prediction.json'
            values[seed] = zero_shot
        else:
            slot = root / 'slots' / regime / f'seed-{seed}'
            transfer = transfer_identity(request, regime)
            identity = {**transfer, 'seed': seed, 'stage': plan['stage']}
            run = verify_transfer_fit(slot, digest, identity, item, lineage)
            manifest = slot / 'prediction/prediction.json'
            seal = read_json(manifest)
            if (seal['run_sha256'] != sha256(slot / 'fit/run.json')
                    or seal['case_count'] != len(cases) or seal['group_count'] != len(groups)):
                raise ValueError('Transfer prediction differs from the sealed fit or frozen cases')
            record = {'prediction': str(manifest.relative_to(root)), 'seed': seed,
                      'prediction_manifest_sha256': sha256(manifest)}
            values[seed] = prediction(root, record, {
                'identity': transfer, 'cases': cases, 'n_genes': request['n_genes'],
                'recipes': {seed: item['recipe_sha256']},
            }, plan['stage'], digest)
            entry['target_checkpoint_sha256'] = run['checkpoint_sha256']
        evidence[str(seed)] = {**entry, 'prediction_manifest_sha256': sha256(manifest)}
    tahoe.write_scores(runtime, request, fold, warm_source, output, arm, values, {
        'schema': 'step05-transfer-evaluation-v1', 'arm': arm, 'regime': regime,
        'stage': plan['stage'], 'plan_sha256': digest, 'recipient_identity': request['identity'],
        'recipient_request_sha256': target['request_sha256'], 'recipient_sha256': target['split_sha256'],
        'complete': plan['stage'] == 'formal', 'dose_limitation': LIMITATION, 'slots': evidence,
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime', type=Path, default=REPO)
    commands = parser.add_subparsers(dest='command', required=True)
    build = commands.add_parser('prepare')
    build.add_argument('root', type=Path)
    build.add_argument('--smoke', action='store_true')
    for name in ('verify', 'source-fit', 'tahoe-predict', 'fit', 'predict', 'tahoe-score', 'score'):
        command = commands.add_parser(name)
        command.add_argument('root', type=Path)
        command.add_argument('--plan-sha256', required=True)
        if name not in ('verify', 'source-fit'):
            command.add_argument('regime', choices=REGIMES)
        if name in ('source-fit', 'tahoe-predict', 'fit', 'predict'):
            command.add_argument('seed', type=int, choices=SEEDS)
        if name in ('score', 'tahoe-score'):
            command.add_argument('output', type=Path)
    args = parser.parse_args()
    runtime, root = args.runtime.resolve(), args.root.absolute()
    try:
        if args.command == 'prepare':
            prepare(runtime, root, 'smoke' if args.smoke else 'formal')
            return
        regime = 'cold-drug' if args.command == 'source-fit' else getattr(args, 'regime', None)
        digest = args.plan_sha256
        plan, context, historical, historical_sha = verify(runtime, root, digest, regime)
        if args.command == 'verify':
            return
        if hasattr(args, 'seed') and str(args.seed) not in plan['active_recipes']:
            raise ValueError('seed is outside the active Transfer plan')
        request = context[regime][0]
        target = plan['targets'][regime]
        cold_run_inputs = cold_inputs(runtime) if regime == 'cold-drug' else None

        def sources(seed):
            if cold_run_inputs is None:
                return historical, historical_sha
            source, source_sha = cold_source(runtime, root, plan, digest, seed, cold_run_inputs)
            overlap = tahoe.overlap(source, request)
            if any(overlap[role][heldout] for role in cold.FITTED for heldout in ('validation', 'test')):
                raise ValueError('globally cold Tahoe fit exposed held-out tenv parents')
            return source, source_sha

        if args.command == 'source-fit':
            cold_source(runtime, root, plan, digest, args.seed, cold_run_inputs, fit=True)
        elif args.command in ('score', 'tahoe-score'):
            score(runtime, root, plan, digest, regime, args.output.absolute(), context, sources,
                  'Tahoe' if args.command == 'tahoe-score' else 'Transfer')
        else:
            source, source_sha = sources(args.seed)
            if args.command == 'tahoe-predict':
                tahoe.predict_recipient(runtime, source, source_sha, request, target['inference'], regime,
                                        target['split_sha256'], args.seed,
                                        tahoe_directory(runtime, root, regime, args.seed))
            else:
                directory, bridge, _ = verify_tahoe(runtime, root, regime, args.seed, request, target,
                                                     source, source_sha, load_values=False)
                item = plan['active_recipes'][str(args.seed)]
                lineage = transfer_lineage(source, source_sha, target, item, bridge, args.seed)
                transfer_worker(runtime, root, plan, digest, regime, args.seed, args.command, request, directory, lineage)
    except (ArithmeticError, KeyError, ValueError, TypeError, IndexError, OSError, EOFError) as error:
        parser.exit(1, f'Transfer {args.command} rejected: {error}\n')


if __name__ == '__main__':
    main()
