#!/usr/bin/env python3
"""Replay original Transfer and test its halved target-stage learning rate."""
import argparse
from pathlib import Path
import tempfile
from unittest.mock import patch

import baseline as base
from artifacts import SEEDS, read_json, sha256, write_csv, write_json
from formal import slot_lock, use_runtime
from formal_matrix import bind, check_binding
from native_lr import NATIVE_CHECKPOINT, common_checkpoint, fit_diagnostic
from screen import paired_gain, screen_decision
from training_inputs import json_hash
from training_policy import MONITOR, verify_best_checkpoint

SCHEMA = 'step06-learning-rate-v1'
ARM = 'Transfer-lr-half'
COMMON_STEPS = {17: 6500, 42: 6000, 2026: 6500}
METRICS = ('positive_ap100', 'negative_ap100', 'mean_ap100')


def arm_config(original):
    if original['learning_rate'] != 0.0002:
        raise ValueError('expected the original Transfer base learning rate of 2e-4')
    return {**original, 'learning_rate': 0.0001}


def publish(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with slot_lock(path.parent / f'.{path.name}-lock'), tempfile.TemporaryDirectory(dir=path.parent) as temp:
        pending = Path(temp) / path.name
        write_json(pending, data)
        if path.exists():
            bind(path, sha256(pending))
        else:
            pending.rename(path)
    return bind(path)


def seed_root(root, seed):
    return root / f'fold-0-seed-{seed}'


def checkpoint_root(runtime, seed):
    number = SEEDS.index(seed) + 4
    return runtime / 'stage2/06_pdiff_transfer_optimize/output/checkpoints' / f'{number:02}-learning-rate-fold-0-seed-{seed}'


def freeze(root, seed, baseline_path, baseline_sha):
    primary = None
    if seed != 17:
        primary_path = seed_root(root, 17) / 'scores/scores.json'
        primary = bind(primary_path)
        result = read_json(primary_path)
        if screen_decision({17: primary_gain(result)})['status'] != 'pending_replication':
            raise ValueError('additional fits require a seed 17 matched-step screening pass')
    baseline_plan = bind(baseline_path, baseline_sha)
    plan = read_json(baseline_path)
    if (plan['schema'] != base.SCHEMA or plan['seed'] != seed
            or plan['request']['role'] != 'validation' or plan['request']['identity']['fold'] != 0):
        raise ValueError('expected the verified Step 06 validation baseline for this seed')
    step = COMMON_STEPS[seed]
    for selected in plan['checkpoints'].values():
        check_binding(selected['checkpoint'])
        check_binding(selected['history'])
    matching_transfer = plan['checkpoints']['Transfer']['metadata']['best_step'] == step
    bridge = plan['source']['bridge']
    bridge_record = read_json(check_binding(bridge)['path'])
    if bridge_record['source']['checkpoint_sha256'] != plan['source']['checkpoint']['sha256']:
        raise ValueError('baseline bridge differs from the verified Tahoe source EMA')
    release = plan['release']
    matrix = read_json(check_binding(release['step05_matrix'])['path'])
    recipient = check_binding(matrix['targets'][base.REGIME]['request_files']['Tahoe'])
    request = read_json(recipient['path'])
    if (base.inference_request(request) != plan['request']
            or base.transfer.request_digest(request) != bridge_record['recipient_request_sha256']):
        raise ValueError('baseline validation and bridge recipient differ')
    producer = read_json(check_binding(release['step05_transfer'])['path'])
    original = producer['active_recipes'][str(seed)]
    if original['config'] != plan['config'] or original['recipe_sha256'] != json_hash(plan['config']):
        raise ValueError('baseline and producer recipes differ')
    config = arm_config(original['config'])
    if step > config['max_steps'] or (step % config['validation_interval'] and step != config['max_steps']):
        raise ValueError('published original Transfer step falls outside the native validation cadence')
    return {'schema': SCHEMA, 'seed': seed, 'fold': 0, 'baseline': baseline_plan, 'primary_screen': primary,
            'recipient': recipient, 'bridge': bridge, 'mapped_state_sha256': bridge_record['mapped_state_sha256'],
            'original_recipe': original, 'recipe': config, 'recipe_sha256': json_hash(config),
            'common_step_policy': {'steps_by_seed': {str(seed): step for seed, step in COMMON_STEPS.items()},
                                   'step': step, 'original_transfer': ('reuse_selected_ema' if matching_transfer
                                                                       else 'matching_checkpoint_unavailable'),
                                   'intervention': 'fit_from_verified_bridge',
                                   'cm': ('reuse_selected_ema' if plan['checkpoints']['CM']['metadata']['best_step']
                                          == step else 'matching_checkpoint_unavailable')}}


def load_plan(slot, seed, digest):
    path = slot / 'plan.json'
    plan = read_json(base.verify_file(path, digest))
    if (plan['schema'] != SCHEMA or plan['seed'] != seed or plan['fold'] != 0
            or plan['common_step_policy']['step'] != COMMON_STEPS[seed]
            or plan['recipe'] != arm_config(plan['original_recipe']['config'])
            or json_hash(plan['recipe']) != plan['recipe_sha256']):
        raise ValueError('invalid learning-rate plan or recipe parity')
    if plan['primary_screen']:
        check_binding(plan['primary_screen'])
    return plan


def native_fit(runtime, directory, plan):
    use_runtime(runtime)
    import native_training as native

    request = read_json(check_binding(plan['recipient'])['path'])
    worker = {**request, 'identity': {**base.transfer.transfer_identity(request, base.REGIME), 'arm': ARM}}
    bridge = check_binding(plan['bridge'])
    initial = {}
    model_class = native.JointPerturbDiff

    class InitializedModel(model_class):
        def __init__(self, inputs, config, directory_out):
            super().__init__(inputs, config, directory_out)
            record = base.transfer.load_initialization(self, Path(bridge['path']).parent, bridge['sha256'], request)
            if record['mapped_state_sha256'] != plan['mapped_state_sha256']:
                raise ValueError('fit initialization differs from the verified Transfer state')
            dose = self.get_parameter(base.transfer.DOSE_ROWS)
            if not dose.requires_grad:
                raise ValueError('mapped dose embeddings must be trainable')
            initial.update(mapped_state_sha256=record['mapped_state_sha256'], mapped_dose_rows_trainable=True)

    load_recipe = native.load_recipe

    def resolve(runtime, path, request):
        original = load_recipe(runtime, path, request)
        if original != plan['original_recipe']['config']:
            raise ValueError('loaded native recipe differs from the original Transfer recipe')
        return arm_config(original)

    with tempfile.TemporaryDirectory(dir=directory.parent) as temp:
        recipe = Path(temp) / 'recipe.json'
        write_json(recipe, plan['original_recipe']['recipe'])
        with patch.multiple(native, JointPerturbDiff=InitializedModel, load_recipe=resolve):
            try:
                fit_diagnostic(native, runtime, worker, recipe, directory, plan['common_step_policy']['step'])
            except InterruptedError:
                pass
    return initial


def checkpoints(directory, plan):
    run = read_json(directory / 'run.json')
    identity = {**read_json(check_binding(plan['recipient'])['path'])['identity'], 'arm': ARM,
                'regime': base.REGIME, 'seed': plan['seed'], 'stage': 'formal'}
    if (run['identity'] != identity or run['recipe'] != plan['recipe'] or run['monitor'] != MONITOR
            or run['recipe_sha256'] != plan['recipe_sha256']):
        raise ValueError('learning-rate fit identity, recipe or selection monitor changed')
    result = {}
    native = directory / 'native-selection.json'
    if native.exists():
        selection = read_json(native)
        metadata = {key: run[key] for key in base.CHECKPOINT_KEYS}
        if (any(selection[key] != value for key, value in metadata.items())
                or selection['history_sha256'] != json_hash(run['history'])
                or selection['stop_step'] != run['history'][-1]['step']
                or selection['monitor'] != MONITOR
                or selection['stop_reason'] != ('validation_early_stop' if run['history'][-1]['early_stop']
                                                else 'max_steps')):
            raise ValueError('native selection changed after the sealed stop')
        verify_best_checkpoint(run, ARM)
        if (run['status'] == 'completed' and (run['reload']['validation_loss'] != run['best_validation_loss']
                                             or run['reload']['prediction_max_abs_error'] != 0)):
            raise ValueError('native loss-selected checkpoint did not reload exactly')
        result['native'] = {'checkpoint': bind(directory / NATIVE_CHECKPOINT, metadata['checkpoint_sha256']),
                            'metadata': metadata, 'selection': bind(native)}
    diagnostic = directory / 'diagnostic.json'
    if diagnostic.exists():
        meta = read_json(diagnostic)
        step = plan['common_step_policy']['step']
        if (meta['best_step'] != step or meta['identity'] != identity
                or meta['recipe_sha256'] != plan['recipe_sha256'] or run['optimizer_steps'] < step):
            raise ValueError('diagnostic checkpoint differs from the frozen common step')
        result['common'] = {'checkpoint': bind(directory / common_checkpoint(step), meta['checkpoint_sha256']),
                            'metadata': meta, 'selection': bind(diagnostic)}
        if ('native' in result and result['native']['metadata']['best_step'] == step
                and result['native']['metadata']['model_state_sha256'] != meta['model_state_sha256']):
            raise ValueError('same-step native and diagnostic EMA states differ')
    return result


def fit(runtime, seed, plan, digest):
    directory = checkpoint_root(runtime, seed)
    directory.parent.mkdir(parents=True, exist_ok=True)
    with slot_lock(directory.parent / f'.{directory.name}-lock'):
        if directory.exists():
            verify_fit(directory, plan, digest)
            return
        attempt = Path(tempfile.mkdtemp(prefix=f'{directory.name}-attempt-', dir=directory.parent))
        staging = attempt / 'ready'
        print(f'Learning-rate fit attempt evidence: {staging}', flush=True)
        initial = native_fit(runtime, staging, plan)
        selected = checkpoints(staging, plan)
        run = read_json(staging / 'run.json')
        write_json(staging / 'fit.json', {'plan_sha256': digest, 'initialization': initial,
                   'status': run['status'], 'reason': run['stop_reason'],
                   'files': {path.name: sha256(path) for path in sorted(staging.iterdir())},
                   'available_modes': sorted(selected)})
        staging.rename(directory)
        attempt.rmdir()


def verify_fit(directory, plan, digest):
    manifest = read_json(directory / 'fit.json')
    if (manifest['plan_sha256'] != digest
            or (manifest['available_modes'] and manifest['initialization'] != {
                'mapped_state_sha256': plan['mapped_state_sha256'], 'mapped_dose_rows_trainable': True})):
        raise ValueError('fit did not consume the frozen learning-rate plan and original bridge')
    for name, expected in manifest['files'].items():
        bind(directory / name, expected)
    selected = checkpoints(directory, plan)
    if sorted(selected) != manifest['available_modes']:
        raise ValueError('fit checkpoint inventory changed')
    return selected


def prediction_mode(mode, selected):
    if mode == 'common' and 'native' in selected:
        common, native = (selected[key]['metadata'] for key in ('common', 'native'))
        if (common['best_step'], common['model_state_sha256']) == (native['best_step'], native['model_state_sha256']):
            return 'native'
    return mode


def predict(runtime, slot, plan, digest):
    use_runtime(runtime)
    directory = checkpoint_root(runtime, plan['seed'])
    selected = verify_fit(directory, plan, digest)
    baseline = read_json(check_binding(plan['baseline'])['path'])
    for mode, item in selected.items():
        if prediction_mode(mode, selected) == mode:
            base.predict_checkpoint(baseline['request'], plan['recipe'], baseline['inference'], item['checkpoint'],
                                    item['metadata'], slot / 'predictions' / mode,
                                    {'plan_sha256': digest, 'arm': ARM, 'mode': mode})


def contrast(cm, original, intervention):
    metrics = {}
    for name in METRICS:
        gap = paired_gain(cm[name], original[name]) if cm is not None else None
        gain = paired_gain(intervention[name], original[name])
        residual = paired_gain(cm[name], intervention[name]) if cm is not None else None
        recovery = gain / gap if gain is not None and gap is not None and gap > 0 else None
        metrics[name] = {key: None if value is None else float(value) for key, value in
                         (('cm_gap', gap), ('gain', gain), ('residual_gap', residual), ('recovery', recovery))}
        if gain is not None:
            metrics[name]['gain_decimal'] = str(gain)
    return {'dose_um': original['dose_um'], **metrics}


def score(runtime, slot, plan, digest):
    seed = plan['seed']
    baseline_binding = check_binding(plan['baseline'])
    baseline_root = Path(baseline_binding['path']).parent
    baseline = read_json(baseline_binding['path'])
    request = baseline['request']
    directory = checkpoint_root(runtime, seed)
    selected = verify_fit(directory, plan, digest)
    values, evidence = {}, {}
    for arm, item in baseline['checkpoints'].items():
        destination = baseline_root / arm
        values[arm] = base.verify_prediction(request, item['checkpoint'], item['metadata'], destination,
                                             {'plan_sha256': baseline_binding['sha256'], 'arm': arm})
        evidence[arm] = bind(destination / 'prediction.json')
    for mode in selected:
        actual = prediction_mode(mode, selected)
        item = selected[actual]
        destination = slot / 'predictions' / actual
        values[mode] = base.verify_prediction(request, item['checkpoint'], item['metadata'], destination,
                                              {'plan_sha256': digest, 'arm': ARM, 'mode': actual})
        evidence[mode] = bind(destination / 'prediction.json')

    full = read_json(check_binding(plan['recipient'])['path'])
    if base.inference_request(full) != request:
        raise ValueError('validation truth request differs from the prediction request')
    truth = bind(full['predictor_files']['validation_expression'], full['file_sha256']['validation_expression'])
    observed = base.load_array(truth['path'], (baseline['case_count'], request['n_genes']))
    grouped, summaries = {}, {}
    for arm, delta in values.items():
        grouped[arm], summaries[arm] = base.score_prediction(request, observed, delta)
    comparisons = {}
    for mode in ('common', 'native'):
        if mode == 'common' and plan['common_step_policy']['original_transfer'] != 'reuse_selected_ema':
            comparisons[mode] = {'status': 'incomplete', 'reason': 'original Transfer common-step EMA unavailable',
                                 'paired': []}
            continue
        if mode not in selected:
            comparisons[mode] = {'status': 'incomplete', 'reason': f'{mode} intervention checkpoint unavailable',
                                 'paired': []}
            continue
        matched = mode == 'native' or plan['common_step_policy']['cm'] == 'reuse_selected_ema'
        cms = summaries['CM'] if matched else [None] * len(summaries['Transfer'])
        comparisons[mode] = {
            'status': 'completed', 'cm_status': 'available' if matched else 'matching_checkpoint_unavailable',
            'paired': [contrast(cm, original, treatment) for cm, original, treatment in
                       zip(cms, summaries['Transfer'], summaries[mode], strict=True)],
        }
    destination = slot / 'scores'
    with slot_lock(slot / '.score-lock'), tempfile.TemporaryDirectory(dir=slot) as temp:
        staging = Path(temp) / 'ready'
        write_csv(staging / 'per_group.csv', ({'arm_or_mode': arm, **row}
                                             for arm, group in grouped.items() for row in group))
        write_csv(staging / 'per_dose.csv', ({'arm_or_mode': arm, **row}
                                            for arm, summary in summaries.items() for row in summary))
        report = {'schema': SCHEMA, 'seed': seed, 'fold': 0, 'plan_sha256': digest, 'truth': truth,
                  'common_step_policy': plan['common_step_policy'], 'fit': bind(directory / 'fit.json'),
                  'baseline_checkpoints': baseline['checkpoints'], 'intervention_checkpoints': selected,
                  'predictions': evidence, 'per_dose': summaries, 'comparisons': comparisons,
                  'files': {path.name: sha256(path) for path in sorted(staging.iterdir())}}
        write_json(staging / 'scores.json', report)
        if destination.exists():
            for path in staging.iterdir():
                bind(destination / path.name, sha256(path))
        else:
            staging.rename(destination)


def primary_gain(report):
    rows = report['comparisons']['common']['paired']
    return next((row['mean_ap100'].get('gain_decimal') for row in rows if row['dose_um'] == 'pooled'), None)


def decide(root, unavailable=None):
    unavailable = unavailable or {}
    if set(unavailable) - set(SEEDS) or any(not reason.strip() for reason in unavailable.values()):
        raise ValueError('unavailable evidence requires a screening seed and a terminal reason')
    gains, evidence = {}, {}
    for seed in SEEDS:
        slot = seed_root(root, seed)
        path = slot / 'scores/scores.json'
        if seed in unavailable:
            if path.exists():
                raise ValueError(f'seed {seed} already has a score seal')
            gains[seed] = None
            continue
        if not path.exists():
            continue
        report = read_json(path)
        if (report['schema'] != SCHEMA or report['seed'] != seed
                or report['common_step_policy']['step'] != COMMON_STEPS[seed]):
            raise ValueError('unexpected matched-step score record')
        plan = load_plan(slot, seed, report['plan_sha256'])
        if report['common_step_policy'] != plan['common_step_policy']:
            raise ValueError('score and plan common-step policies differ')
        for item in [report['fit'], report['truth'], *report['predictions'].values()]:
            check_binding(item)
        for name, digest in report['files'].items():
            bind(path.parent / name, digest)
        gains[seed] = primary_gain(report)
        evidence[str(seed)] = bind(path)
    decision = screen_decision(gains)
    result = {'schema': SCHEMA, 'fold': 0, 'gate_mode': 'common', **decision,
              'common_steps': {str(seed): step for seed, step in COMMON_STEPS.items()},
              'gains': {str(seed): gain for seed, gain in gains.items()}, 'scores': evidence}
    if unavailable:
        result['unavailable'] = {str(seed): reason for seed, reason in unavailable.items()}
    if not decision['next_seeds']:
        publish(root / 'decision.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('verify', 'fit', 'predict', 'score', 'decision'))
    parser.add_argument('root', type=Path)
    parser.add_argument('--runtime', type=Path, default=base.REPO)
    parser.add_argument('--seed', type=int, choices=SEEDS, default=17)
    parser.add_argument('--baseline-plan', type=Path)
    parser.add_argument('--baseline-plan-sha256')
    parser.add_argument('--plan-sha256')
    parser.add_argument('--unavailable', nargs=2, metavar=('SEED', 'REASON'), action='append', default=[])
    args = parser.parse_args()
    if args.unavailable and args.phase != 'decision':
        parser.error('--unavailable is only accepted by decision')
    root, runtime = args.root.resolve(), args.runtime.resolve()
    slot = seed_root(root, args.seed)
    if args.phase == 'decision':
        unavailable = {int(seed): reason for seed, reason in args.unavailable}
        if len(unavailable) != len(args.unavailable):
            parser.error('duplicate unavailable screening seed')
        print(decide(root, unavailable), flush=True)
        return
    if args.phase == 'verify':
        if not args.baseline_plan or not args.baseline_plan_sha256:
            parser.error('verify requires --baseline-plan and --baseline-plan-sha256')
        plan = freeze(root, args.seed, args.baseline_plan.resolve(), args.baseline_plan_sha256)
        print(f"Learning-rate plan SHA-256 {publish(slot / 'plan.json', plan)['sha256']}", flush=True)
        return
    if not args.plan_sha256:
        parser.error(f'{args.phase} requires --plan-sha256')
    plan = load_plan(slot, args.seed, args.plan_sha256)
    if args.phase == 'fit':
        fit(runtime, args.seed, plan, args.plan_sha256)
    elif args.phase == 'predict':
        predict(runtime, slot, plan, args.plan_sha256)
    else:
        score(runtime, slot, plan, args.plan_sha256)


if __name__ == '__main__':
    main()
