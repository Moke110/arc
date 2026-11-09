#!/usr/bin/env python3
"""Fit, select and pair the random and log dose arms through the ADR 0049 runner on CPU."""
import argparse
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))

import paired_log_source as log_source
import paired_reference as ref
import paired_target as cli
from artifacts import read_json, sha256, write_json
from formal_matrix import bind
from native_log_dose_fixture import source_release
from native_replication_fixture import source_checkpoint, target_request
from tahoe_bridge import load_initialization, seal_initialization
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
    target_root = root / f'output/paired-transfer/{arm}/seed-17'
    reference_fit = reference_root / 'fit'

    class Original(native.JointPerturbDiff):
        def __init__(self, inputs, config, directory):
            super().__init__(inputs, config, directory)
            load_initialization(self, root, bridge_sha, request)

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
        if arm == 'log':
            source_config = {**config, 'max_steps': 2}
            plan['source'] = source_release(native, root / 'tahoe-release', source_config,
                                            request['identity']['gene_sha256'])
        write_json(reference_root / 'plan.json', plan)
        reference_sha = sha256(reference_root / 'plan.json')
        common = ['--seed', '17', '--runtime', str(runtime)]
        for phase in ('predict', 'score'):
            ref.main([phase, str(reference_root), *common, '--plan-sha256', reference_sha])
        if arm == 'log':
            with (patch.object(native, 'load_recipe', return_value=source_config),
                  patch.object(log_source, 'verify_request', return_value=None)):
                log_source.main(['fit', str(target_root), *common, '--reference', str(reference_root),
                                 '--plan-sha256', reference_sha])
        args = [str(target_root), '--arm', arm, *common]
        cli.main(['verify', *args, '--reference-sha256', reference_sha])
        target_plan = target_root / 'plan.json'
        if arm == 'log':
            verified = read_json(target_plan)
            wrong = {**verified, 'source': {**verified['source'], 'source_plan_sha256': 'other-source'}}
            write_json(target_plan, wrong)
            with unittest.TestCase().assertRaisesRegex(ValueError, 'differs from the frozen Transfer'):
                cli.main(['fit', *args, '--plan-sha256', sha256(target_plan)])
            write_json(target_plan, verified)
        args.extend(['--plan-sha256', sha256(target_plan)])
        with patch.object(native, 'validate', side_effect=[(loss, [loss], [1])
                                                            for loss in (.7, .5, .5, .8, .7, .6, .5)]):
            cli.main(['fit', *args])
        for phase in ('fit', 'predict', 'predict', 'score', 'score'):
            cli.main([phase, *args])
    paired = read_json(target_root / 'paired.json')
    scores = read_json(target_root / 'scores/scores.json')
    assert scores['optimizer_steps'] == 6
    assert scores['selection']['metadata']['best_step'] == 2
    fit = read_json(scores['selection']['history']['path'])
    initialization = fit['initialization']
    if arm == 'log':
        source = log_source.selected_source(target_root, plan, reference_sha)
        assert fit['source'] == source
        assert initialization['source_selection'] == source['metadata']
        assert initialization['recipient_doses_um'] == [.2, 1., 10.]
        assert abs(initialization['recipient_log_features'][-1] - 2.30103) < 1e-6
        assert 'native.cov_encoder.dose_encoder.control' in {
            parameter for group in fit['training_initial_state']['optimizer_groups']
            for parameter in group['parameters']
        }
    else:
        assert initialization['reset_doses_um'] == [.2, 1., 10.]
        assert initialization['native_positive_rows_sha256'] == initialization['restored_positive_rows_sha256']
        assert initialization['mapped_state_sha256'] != initialization['reset_state_sha256']
    row = paired['per_dose'][0]
    assert (row['step'], row['reference_step']) == (2, 1)
    assert row['pcc_difference'] == row['pcc'] - row['reference_pcc']
    print(f'PASS {arm} fit, strict loss selection, selected EMA replay and same-seed pairing')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime', required=True, type=Path)
    parser.add_argument('--arm', choices=('random', 'log'), default='random')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as temp:
        fixture(args.runtime.resolve(), Path(temp), args.arm)
