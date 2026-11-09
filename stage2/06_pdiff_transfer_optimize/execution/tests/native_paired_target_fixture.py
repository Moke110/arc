#!/usr/bin/env python3
"""Fit, select and pair target interventions through the ADR 0049 runner on CPU."""
import argparse
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))

import paired_reference as ref
import paired_target as cli
from artifacts import read_json, sha256, write_json
from formal_matrix import bind
from native_replication_fixture import source_checkpoint, target_request
from tahoe_bridge import load_initialization, seal_initialization, state_digest
from training_inputs import json_hash


def fixture(runtime, root, arm='random'):
    ref.use_runtime(runtime)
    import native_training as native

    request = target_request(root, native.np)
    request['identity']['fold'] = 0
    config = {'stage': 'formal', 'seed': 17, 'fold': 0, 'threads': 2, 'device': 'cpu',
              'expression_scale': 1, 'chemical_dimension': 16, 'batch_groups': 1,
              'population_size': 2, 'learning_rate': .0002, 'weight_decay': .001,
              'sampling_start_time': 1, 'warmup_steps': 1, 'schedule_training_steps': 6,
              'ema_decay': .9, 'ema_update_steps': 1, 'max_steps': 6,
              'validation_interval': 1, 'minimum_steps': 7, 'patience_steps': 7,
              'max_training_seconds': 120, 'min_delta': 0}
    recipient, recipe = root / 'request.json', root / 'recipe.json'
    write_json(recipient, request)
    write_json(recipe, config)
    source = source_checkpoint(native, root, request, config)
    bridge = root / 'bridge.json'
    native.seed_randomness(17)
    bridge_sha = seal_initialization(native.JointPerturbDiff(native.PredictorInputs(request), config, root),
                                     root, source, request)
    reference_root = root / 'output/paired-transfer/reference/seed-17'
    target_root = root / 'output/paired-transfer' / arm / 'seed-17'
    reference_fit = reference_root / 'fit'

    original_names, original_groups = {}, []
    original_rates, target_rates, optimizer_updates = [], [], []

    class Original(native.JointPerturbDiff):
        def __init__(self, inputs, config, directory):
            super().__init__(inputs, config, directory)
            load_initialization(self, root, bridge_sha, request)
            if arm == 'all-half':
                original_names.update({id(parameter): name for name, parameter in self.named_parameters()
                                       if parameter.requires_grad})

    native_optimizer = native.get_optimizer

    def original_optimizer(groups, config):
        optimizer = native_optimizer(groups, config)
        if arm == 'all-half':
            original_groups.extend({
                'base_lr': group['lr'], 'weight_decay': group['weight_decay'],
                'parameters': [original_names[id(parameter)] for parameter in group['params']],
            } for group in optimizer.param_groups)
            optimizer.register_step_post_hook(
                lambda opt, _args, _kwargs: original_rates.append([group['lr'] for group in opt.param_groups])
            )
        return optimizer

    def group_digests(optimizer):
        return [state_digest({str(i): parameter for i, parameter in enumerate(group['params'])})
                for group in optimizer.param_groups]

    def target_optimizer(groups, config):
        optimizer = native_optimizer(groups, config)
        if arm == 'all-half':
            before = group_digests(optimizer)

            def observe_step(optimizer, _args, _kwargs):
                rates = [group['lr'] for group in optimizer.param_groups]
                target_rates.append(rates)
                if not optimizer_updates and all(rate > 0 for rate in rates):
                    optimizer_updates.extend(current != initial for current, initial
                                             in zip(group_digests(optimizer), before, strict=True))
            optimizer.register_step_post_hook(observe_step)
        return optimizer

    cases = ref.base.case_rows(ref.base.inference_request(request))
    frozen = {'request': ref.base.inference_request(request), 'fold': 0, 'regime': ref.REGIME,
              'inference': {'seed_offset': 200000, 'max_groups': 1, 'max_profiles': 2},
              'case_count': 2, 'group_count': 1, 'dose_group_counts': {'1.0': 1},
              'truth': bind(request['predictor_files']['validation_expression']),
              'ordered_case_group_sha256': json_hash(ref.base.case_order(cases))}
    original_request = {**request, 'identity': {**request['identity'], 'arm': ref.ARM, 'regime': ref.REGIME}}
    with (patch.multiple(native.torch.cuda, reset_peak_memory_stats=lambda device: None,
                         synchronize=lambda device: None),
          patch.object(native, 'load_recipe', return_value=config)):
        # Controlled author losses exercise distinct minima and strict ties while
        # the trainer still performs actual optimizer, EMA and prediction steps.
        with (patch.object(native, 'JointPerturbDiff', Original),
              patch.object(native, 'get_optimizer', original_optimizer),
              patch.object(native, 'validate', side_effect=[(loss, [loss], [1])
                                                            for loss in (.5, .7, .6, .8, .6, .7, .5)])):
            native.fit(runtime, original_request, recipe, reference_fit)
        run = read_json(reference_fit / 'run.json')
        selected = {'checkpoint': bind(reference_fit / 'best.pt'), 'history': bind(reference_fit / 'run.json'),
                    'metadata': {key: run[key] for key in ref.base.CHECKPOINT_KEYS}}
        plan = {**frozen, 'schema': ref.SCHEMA, 'seed': 17, 'protocol_sha256': json_hash(frozen),
                'config': config, 'recipient_request': bind(recipient), 'target_recipe': bind(recipe),
                'checkpoints': {ref.ARM: selected}, 'source': {
                    'bridge': bind(bridge), 'checkpoint': {'path': source['path'],
                                                           'sha256': source['checkpoint_sha256']},
                    'metadata': {**source}}, 'release': {'recipient_split_sha256': 'fixture'}}
        write_json(reference_root / 'plan.json', plan)
        reference_sha = sha256(reference_root / 'plan.json')
        common = ['--seed', '17', '--runtime', str(runtime)]
        for phase in ('predict', 'score'):
            ref.main([phase, str(reference_root), *common, '--plan-sha256', reference_sha])
        args = [str(target_root), '--arm', arm, *common]
        cli.main(['verify', *args, '--reference-sha256', reference_sha])
        args.extend(['--plan-sha256', sha256(target_root / 'plan.json')])
        with (patch.object(native, 'validate', side_effect=[(loss, [loss], [1])
                                                            for loss in (.7, .5, .5, .8, .7, .6, .5)]),
              patch.object(native, 'get_optimizer', target_optimizer)):
            cli.main(['fit', *args])
        for phase in ('fit', 'predict', 'predict', 'score', 'score'):
            cli.main([phase, *args])
    paired = read_json(target_root / 'paired.json')
    scores = read_json(target_root / 'scores/scores.json')
    assert scores['optimizer_steps'] == 6
    assert scores['selection']['metadata']['best_step'] == 2
    fit = read_json(scores['selection']['history']['path'])
    if arm == 'random':
        reset = fit['initialization']
        assert reset['reset_doses_um'] == [.2, 1., 10.]
        assert reset['native_positive_rows_sha256'] == reset['restored_positive_rows_sha256']
        assert reset['mapped_state_sha256'] != reset['reset_state_sha256']
    else:
        target_plan = read_json(target_root / 'plan.json')
        assert target_plan['config'] == {**config, 'learning_rate': 1e-4}
        assert Path(target_plan['checkpoint_slot']).name == '06-paired-all-half-fold-0-seed-17'
        assert fit['initialization']['mapped_state_sha256'] == read_json(bridge)['mapped_state_sha256']
        assert len(original_groups) == 2
        assert all(group['base_lr'] == 2e-4 for group in original_groups)
        groups = fit['training_initial_state']['optimizer_groups']
        assert groups == [{**group, 'base_lr': 1e-4} for group in original_groups]
        checkpoint = native.torch.load(scores['selection']['checkpoint']['path'], map_location='cpu',
                                       weights_only=True, mmap=True)
        assert checkpoint['recipe'] == target_plan['config']
        assert len(original_rates) == 6
        assert target_rates == [[rate / 2 for rate in rates] for rates in original_rates]
        assert optimizer_updates == [True, True]
    row = paired['per_dose'][0]
    assert (row['step'], row['reference_step']) == (2, 1)
    assert row['pcc_difference'] == row['pcc'] - row['reference_pcc']
    print(f'PASS {arm} fit, strict loss selection, selected EMA replay and same-seed pairing')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime', required=True, type=Path)
    parser.add_argument('--arm', choices=('random', 'all-half'), default='random')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as temp:
        fixture(args.runtime.resolve(), Path(temp), args.arm)
