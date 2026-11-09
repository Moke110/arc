#!/usr/bin/env python3
"""Fit and score one Step 06 intervention against same-seed original Transfer."""
import argparse
from pathlib import Path
from unittest.mock import patch

import paired_reference as reference
from artifacts import SEEDS, read_json, sha256, write_json
from dose_reset import load_dose_reset
from formal import completed_fit, fit_slot, slot_lock, use_runtime, verify_fit_seal
from formal_matrix import bind, check_binding
from learning_rate import publish
from tahoe_bridge import (CHECKPOINT_FIELDS, DOSE_ROWS, load_initialization,
                          request_digest, state_digest)
from training_inputs import json_hash

SCHEMA = 'step06-paired-target-v1'
PAIR_SCHEMA = 'step06-paired-differences-v1'
CHECKPOINT_NUMBERS = {'random': '04', 'all-half': '06'}
ORIGINAL_BASE_LR = 2e-4
HALF_BASE_LR = 1e-4
METRICS = ('mean_ap100', 'pcc')


def reference_scores(binding, plan, digest):
    path = Path(check_binding(binding)['path'])
    scores = read_json(path)
    expected = {
        'schema': reference.SCORE_SCHEMA, 'plan_sha256': digest,
        'protocol_sha256': plan['protocol_sha256'], 'seed': plan['seed'],
        'regime': plan['regime'], 'fold': plan['fold'], 'arm': reference.ARM,
        'selection': plan['checkpoints'][reference.ARM], 'truth': plan['truth'],
    }
    if (any(scores[key] != value for key, value in expected.items())
            or {row['dose_um']: row['group_count'] for row in scores['per_dose']} != plan['dose_group_counts']):
        raise ValueError('reference scores differ from the same-seed frozen validation protocol')
    return scores


def resolve_config(arm, original):
    if arm == 'all-half':
        if original['learning_rate'] != ORIGINAL_BASE_LR:
            raise ValueError('all-half requires the original Transfer base learning rate of 2e-4')
        return {**original, 'learning_rate': HALF_BASE_LR}
    return original


def freeze(root, arm, seed, reference_root, reference_digest):
    reference_plan = reference.load_plan(reference_root, seed, reference_digest)
    source = reference_plan['source']
    recipient = read_json(check_binding(reference_plan['recipient_request'])['path'])
    bridge = read_json(check_binding(source['bridge'])['path'])
    source_checkpoint = source['checkpoint']
    if (reference_plan['config']['seed'] != seed
            or reference.base.inference_request(recipient) != reference_plan['request']
            or bridge['recipient_identity'] != recipient['identity']
            or bridge['recipient_request_sha256'] != request_digest(recipient)
            or bridge['source']['path'] != source_checkpoint['path']
            or bridge['source']['checkpoint_sha256'] != source_checkpoint['sha256']
            or bridge['source']['identity']['seed'] != seed
            or any(bridge['source'][key] != source['metadata'][key] for key in CHECKPOINT_FIELDS)):
        raise ValueError('original Tahoe bridge differs from the same-seed source or recipient')
    check_binding(reference_plan['target_recipe'])
    scores = bind(reference_root / 'scores/scores.json')
    reference_scores(scores, reference_plan, reference_digest)
    slot = root.parents[2] / 'checkpoints' / f'{CHECKPOINT_NUMBERS[arm]}-paired-{arm}-fold-0-seed-{seed}'
    plan = {
        'schema': SCHEMA, 'arm': arm, 'seed': seed, 'fold': 0, 'regime': reference.REGIME,
        'reference': bind(reference_root / 'plan.json', reference_digest), 'reference_scores': scores,
        'protocol_sha256': reference_plan['protocol_sha256'], 'source': source,
        'config': resolve_config(arm, reference_plan['config']), 'checkpoint_slot': str(slot),
    }
    return publish(root / 'plan.json', plan)


def load_plan(root, arm, seed, digest):
    plan = read_json(bind(root / 'plan.json', digest)['path'])
    if (plan['schema'], plan['arm'], plan['seed'], plan['fold'], plan['regime']) != (
            SCHEMA, arm, seed, 0, reference.REGIME):
        raise ValueError('target plan does not match the requested arm, split or seed')
    pinned = plan['reference']
    original = reference.load_plan(Path(pinned['path']).parent, seed, pinned['sha256'])
    if (plan['protocol_sha256'] != original['protocol_sha256']
            or plan['config'] != resolve_config(arm, original['config']) or plan['source'] != original['source']):
        raise ValueError('target plan differs from the frozen Transfer protocol, recipe or source')
    return plan, original


def fit_identity(plan, original):
    return {**original['request']['identity'], 'arm': plan['arm'], 'regime': plan['regime'],
            'seed': plan['seed'], 'stage': plan['config']['stage']}


def verify_training(arm, training):
    if (training['optimizer_state_entries'], training['scheduler_start_step'], training['ema_start_step']) != (0, 0, 0):
        raise ValueError('target adaptation requires a fresh optimizer, scheduler and EMA')
    if arm == 'random' and not training.get('positive_dose_gradient_observed'):
        raise ValueError('no training gradient reached the restored positive-dose rows')
    if arm == 'all-half':
        groups = training['optimizer_groups']
        if not groups or any(group['base_lr'] != HALF_BASE_LR for group in groups):
            raise ValueError('all-half requires a base learning rate of 1e-4 for every optimizer group')


def selected_fit(plan, original, digest):
    directory = Path(plan['checkpoint_slot']) / 'fit'
    verify_fit_seal(directory, digest)
    run = completed_fit(directory, fit_identity(plan, original), plan['config'])
    groups = list(dict.fromkeys(row['group_id'] for row in reference.base.case_rows(original['request'])))
    batching = run['validation_batching']
    verify_training(plan['arm'], run['training_initial_state'])
    if (batching['group_ids'] != groups or batching['group_ids_sha256'] != json_hash(groups)
            or run['source'] != plan['source'] or run['optimizer_steps'] != run['history'][-1]['step']
            or run['stop_reason'] not in ('max_steps', 'validation_early_stop')):
        raise ValueError('target fit differs from the pinned validation groups, source or stopping rule')
    return {'checkpoint': bind(directory / 'best.pt', run['checkpoint_sha256']),
            'history': bind(directory / 'run.json'), 'fit_seal': bind(directory / 'fit.json'),
            'metadata': {key: run[key] for key in reference.base.CHECKPOINT_KEYS}}


def initialize(model, plan, recipient, training):
    bridge = plan['source']['bridge']
    directory = Path(bridge['path']).parent
    if plan['arm'] == 'all-half':
        return load_initialization(model, directory, bridge['sha256'], recipient)
    _, evidence = load_dose_reset(model, directory, bridge['sha256'], recipient)
    parameter = model.get_parameter(DOSE_ROWS)

    def observe(gradient):
        if gradient[1:].count_nonzero().item():
            training['positive_dose_gradient_observed'] = True
            hook.remove()

    hook = parameter.register_hook(observe)
    return evidence


def model_class(native, plan, initialize_model=None):
    class PairedModel(native.JointPerturbDiff):
        def __init__(self, inputs, config, directory):
            super().__init__(inputs, config, directory)
            if initialize_model is not None:
                initialize_model(self)

    return PairedModel


def fit(runtime, plan, original, digest):
    use_runtime(runtime)
    import native_training as native

    slot = Path(plan['checkpoint_slot'])
    with slot_lock(slot):
        if (slot / 'fit').exists():
            return selected_fit(plan, original, digest)

        recipient = read_json(check_binding(original['recipient_request'])['path'])
        if reference.base.inference_request(recipient) != original['request']:
            raise ValueError('target fit request differs from the frozen predictor inputs')
        for key in (*reference.base.INFERENCE_FILES, 'train_expression', 'validation_expression'):
            bind(recipient['predictor_files'][key], recipient['file_sha256'][key])
        identity = fit_identity(plan, original)
        worker = {**recipient, 'identity': {**recipient['identity'], 'arm': plan['arm'], 'regime': plan['regime']}}
        recipe = Path(check_binding(original['target_recipe'])['path'])
        if native.load_recipe(runtime, recipe, worker) != original['config']:
            raise ValueError('target recipe differs from the frozen Transfer resolution')

        evidence, training = {}, {}
        parameter_names = {}
        native_fit, native_optimizer = native.fit, native.get_optimizer
        native_scheduler = native.get_cosine_schedule_with_warmup

        def init_model(model):
            evidence.update(initialize(model, plan, recipient, training))
            parameter_names.update({id(parameter): name for name, parameter in model.named_parameters()})

        def optimizer(groups, config):
            result = native_optimizer(groups, config)
            training['optimizer_state_entries'] = len(result.state)
            training['optimizer_groups'] = [{
                'base_lr': group['lr'], 'weight_decay': group['weight_decay'],
                'parameters': [parameter_names[id(parameter)] for parameter in group['params']],
            } for group in result.param_groups]
            return result

        def schedule(*args, **kwargs):
            result = native_scheduler(*args, **kwargs)
            training['scheduler_start_step'] = result.last_epoch
            return result

        class FreshEMA(native.EMAWeightAveraging):
            def setup(self, trainer, model, stage):
                super().setup(trainer, model, stage)
                if (self._average_model.n_averaged.item() or trainer.global_step
                        or state_digest(self._average_model.module.state_dict()) != state_digest(model.state_dict())):
                    raise ValueError('fresh EMA must start at the initialized model state')
                training['ema_start_step'] = trainer.global_step

        def fit_with_initialization(runtime, request, recipe, directory):
            with patch.multiple(native, JointPerturbDiff=model_class(native, plan, init_model),
                                get_optimizer=optimizer, get_cosine_schedule_with_warmup=schedule,
                                EMAWeightAveraging=FreshEMA,
                                load_recipe=lambda *_: plan['config']):
                native_fit(runtime, request, recipe, directory)
            verify_training(plan['arm'], training)
            history = directory / 'run.json'
            write_json(history, {**read_json(history), 'source': plan['source'], 'initialization': evidence,
                                 'training_initial_state': training})

        with patch.object(native, 'fit', fit_with_initialization):
            fit_slot(slot, digest, runtime, worker, identity, recipe, plan['config'])
        return selected_fit(plan, original, digest)


def paired_differences(scores, baseline):
    original = {row['dose_um']: row for row in baseline['per_dose']}
    rows = []
    for row in scores['per_dose']:
        base = original[row['dose_um']]
        paired = {**row, 'reference_step': base['step']}
        for metric in METRICS:
            paired[f'reference_{metric}'] = base[metric]
            paired[f'{metric}_difference'] = (row[metric] - base[metric]
                                              if row[metric] is not None and base[metric] is not None else None)
        rows.append(paired)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('verify', 'fit', 'predict', 'score'))
    parser.add_argument('root', type=Path)
    parser.add_argument('--runtime', type=Path, default=reference.base.REPO)
    parser.add_argument('--arm', choices=CHECKPOINT_NUMBERS, required=True)
    parser.add_argument('--seed', type=int, choices=SEEDS, required=True)
    parser.add_argument('--reference', type=Path)
    parser.add_argument('--reference-sha256')
    parser.add_argument('--plan-sha256')
    args = parser.parse_args(argv)
    root, runtime = args.root.resolve(), args.runtime.resolve()
    if args.phase == 'verify':
        if not args.reference_sha256:
            parser.error('verify requires --reference-sha256')
        upstream = args.reference or root.parent.parent / 'reference' / f'seed-{args.seed}'
        result = freeze(root, args.arm, args.seed, upstream.resolve(), args.reference_sha256)
    else:
        digest = args.plan_sha256
        if not digest:
            parser.error(f'{args.phase} requires --plan-sha256')
        plan, original = load_plan(root, args.arm, args.seed, digest)
        if args.phase == 'fit':
            result = fit(runtime, plan, original, digest)['fit_seal']
        else:
            selected = selected_fit(plan, original, digest)
            lineage = {'plan_sha256': digest, 'protocol_sha256': plan['protocol_sha256'], 'arm': args.arm}
            prediction = root / 'prediction'
            if args.phase == 'predict':
                use_runtime(runtime)
                reference.base.predict_checkpoint(original['request'], plan['config'], original['inference'],
                                                  selected['checkpoint'], selected['metadata'], prediction, lineage)
                result = bind(prediction / 'prediction.json')
            else:
                baseline = reference_scores(plan['reference_scores'], original, plan['reference']['sha256'])
                score = reference.score_paired(root, original, digest, arm=args.arm, selected=selected,
                                               prediction_dir=prediction, lineage=lineage)
                result = publish(root / 'paired.json', {
                    'schema': PAIR_SCHEMA, 'plan_sha256': digest, 'protocol_sha256': plan['protocol_sha256'],
                    'arm': args.arm, 'seed': args.seed, 'scores': score, 'reference': plan['reference_scores'],
                    'per_dose': paired_differences(read_json(score['path']), baseline),
                })
    print(f'{result["path"]} SHA-256 {result["sha256"]}', flush=True)


if __name__ == '__main__':
    main()
