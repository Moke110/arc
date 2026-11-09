#!/usr/bin/env python3
"""Bounded fold 1 release, native selection, independent reload and no-run fixture."""
import argparse
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))

import confirmation as cli
from artifacts import SEEDS, read_json, rows, sha256, write_csv, write_json
from formal_matrix import bind
from training_inputs import json_hash, request_for


def score_screens(output, dose_gains):
    for seed, value in dose_gains.items():
        slot = cli.dose_reset.run_root(output, seed)
        score = slot / 'scores'
        score.mkdir(parents=True, exist_ok=True)
        write_json(slot / 'plan.json', {})
        files = {}
        for name in ('per_group.csv', 'per_dose.csv'):
            (score / name).write_text('fixture\n')
            files[name] = sha256(score / name)
        write_json(score / 'scores.json', {
            'schema': cli.dose_reset.SCORE_SCHEMA, 'fold': 0, 'seed': seed, 'comparison_mode': 'loss-selected',
            'plan_sha256': sha256(slot / 'plan.json'), 'files': files,
            'per_dose': [{'dose_um': 'pooled', 'arm': 'Transfer', 'mean_ap100': 0.1},
                         {'dose_um': 'pooled', 'arm': cli.dose_reset.ARM, 'mean_ap100': value}],
        })
    cli.dose_reset.write_decision(output, {})
    root = output / 'learning-rate'
    slot = cli.learning_rate.seed_root(root, 17)
    score = slot / 'scores'
    score.mkdir(parents=True, exist_ok=True)
    original = {'learning_rate': 0.0002}
    config = cli.learning_rate.arm_config(original)
    policy = {'step': cli.learning_rate.COMMON_STEPS[17]}
    plan = {'schema': cli.learning_rate.SCHEMA, 'seed': 17, 'fold': 0, 'common_step_policy': policy,
            'primary_screen': None, 'recipe': config, 'original_recipe': {'config': original},
            'recipe_sha256': json_hash(config)}
    write_json(slot / 'plan.json', plan)
    write_json(score / 'scores.json', {
        'schema': cli.learning_rate.SCHEMA, 'seed': 17, 'fold': 0,
        'common_step_policy': policy, 'plan_sha256': sha256(slot / 'plan.json'),
        'fit': bind(slot / 'plan.json'), 'truth': bind(slot / 'plan.json'), 'files': {}, 'predictions': {},
        'comparisons': {'common': {'paired': []}},
    })
    assert cli.learning_rate.decide(root)['status'] == 'incomplete'
    return root


def source_context(root, config, np):
    paths = cli.cold.PATHS
    source = root / 'tahoe'
    (source / 'fold_0').mkdir(parents=True)
    (source / paths['genes']).write_text('symbol\nA\nB\nC\n')
    membership = [
        ('train_keep', 'train', 'allowed-train'), ('train_drop', 'train', 'tenv-validation'),
        ('val_drop', 'validation', 'tenv-test'), ('val_keep', 'validation', 'allowed-val'),
        ('test_keep', 'test', 'tenv-validation'),
    ]
    write_csv(source / 'folds.csv', [{'drug_key': drug, 'role': role, 'parent_connectivity': parent, 'fold': 0}
                                    for drug, role, parent in membership])
    write_csv(source / paths['drugs'], [{'drug_key': drug, 'morgan_index': i, 'parent_connectivity': parent}
                                       for i, (drug, _, parent) in enumerate(membership)])
    np.save(source / paths['controls'], np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]], dtype='float32'))
    np.save(source / paths['morgan'], np.eye(len(membership), 2048, dtype='float32'))
    for role in cli.cold.ROLES:
        cases = [{'case_id': f'{drug}-{i}', 'drug_key': drug, 'dose_um': '1', 'control_index': i,
                  'group_id': drug} for drug, r, _ in membership if role == r for i in range(2)]
        write_csv(source / paths[f'{role}_cases'], cases)
        if role != 'test':
            expression = np.arange(len(cases) * 3, dtype='float32').reshape(len(cases), 3) / 16
            np.save(source / paths[f'{role}_expression'], expression)
    files = {key: sha256(source / name) for key, name in paths.items()}
    split = sha256(source / 'folds.csv')
    manifest = source / 'manifest.json'
    write_json(manifest, {'schema': 'tahoe-input-v1',
                          'files': {**{name: files[key] for key, name in paths.items()}, 'folds.csv': split}})
    request = request_for(source, {'dataset': 'tahoe', 'fold': 0, 'gene_sha256': files['genes'],
                                   'split_sha256': split, 'source_input_sha256': sha256(manifest),
                                   'predictor_files': paths, 'file_sha256': files, 'target_kind': 'delta',
                                   'n_genes': 3})
    target = {'request': {'identity': {'dataset': 'tenv', 'fold': 1, 'gene_sha256': files['genes']}},
              'roles': {role: {'parent_connectivity': parents} for role, parents in
                        [('train', ['allowed-train']), ('validation', ['tenv-validation']), ('test', ['tenv-test'])]}}
    raw = root / 'recipe.json'
    write_json(raw, {'seed': 17, 'stage': 'formal', 'fold': 0})
    producer = {'request': request, 'recipes': {'17': {**bind(raw), 'resolved_sha256': json_hash(config)}},
                **{key: bind(manifest) for key in ('producer_plan', 'producer_request', 'fixed_config')}}
    return target, producer, source


def expect_error(fn, message):
    try:
        fn()
    except ValueError as error:
        assert message in str(error), str(error)
    else:
        raise AssertionError(f'expected error containing {message}')


def fixture(runtime, root):
    cli.use_runtime(runtime)
    import native_training as native

    skipped = root / 'skipped'
    skipped_lr = score_screens(skipped, {17: 0.11, 42: 0.1, 2026: 0.11})
    with patch.object(cli, 'recipient', side_effect=AssertionError('skipped plan opened a release')):
        digest = cli.freeze(runtime, skipped, skipped_lr)['sha256']
        plan = cli.load_plan(runtime, skipped, digest)
        assert plan['status'] == 'skipped' and plan['screens']['dose_reset']['status'] == 'inconsistent'
        assert plan['screens']['learning_rate']['status'] == 'incomplete'
        expect_error(lambda: cli.prepare_release(skipped, plan, digest), 'require a promoted intervention')
        assert not (skipped / cli.RELEASE_PATH).exists()
        assert not any(cli.source_path(skipped, seed).exists() for seed in SEEDS)
        cli.source_path(skipped, 17).mkdir(parents=True)
        expect_error(lambda: cli.load_plan(runtime, skipped, digest), 'skipped confirmation contains')

    output = root / 'promoted'
    lr_root = score_screens(output, {17: 0.11, 42: 0.101, 2026: 0.11})
    config = {'stage': 'formal', 'seed': 17, 'fold': 0, 'threads': 2, 'device': 'cpu',
              'expression_scale': 1, 'chemical_dimension': 16, 'batch_groups': 1,
              'population_size': 2, 'learning_rate': 0.0002, 'weight_decay': 0,
              'sampling_start_time': 1, 'warmup_steps': 1, 'schedule_training_steps': 4,
              'ema_decay': 0.9, 'ema_update_steps': 1, 'max_steps': 2,
              'validation_interval': 1, 'minimum_steps': 2, 'patience_steps': 1,
              'max_training_seconds': 120, 'min_delta': 1e6}
    target, producer, original = source_context(root, config, native.np)
    with patch.object(cli, 'recipient', return_value=target), patch.object(cli, 'source_predecessor', return_value=producer):
        digest = cli.freeze(runtime, output, lr_root)['sha256']
    plan = cli.load_plan(runtime, output, digest)
    assert plan['status'] == 'ready' and plan['promoted_arms'] == ['dose_reset']
    assert plan['screens']['dose_reset']['gains']['17'] == '0.01'
    assert plan['seeds'] == list(SEEDS)
    with patch.object(cli.cold, 'N_GENES', 3):
        request = cli.prepare_release(output, plan, digest)
    assert cli.verify_release(output, plan, digest) == request
    for role, parent in [('train', 'train_drop'), ('validation', 'val_drop')]:
        key = f'{role}_cases'
        path = Path(request['predictor_files'][key])
        snapshot = path.read_bytes()
        forbidden = next(row for row in rows(original / cli.cold.PATHS[key]) if row['drug_key'] == parent)
        write_csv(path, [*rows(path), forbidden])
        expect_error(lambda: cli.verify_release(output, plan, digest), f'forbidden parent in Tahoe {role}')
        path.write_bytes(snapshot)
    assert [row['drug_key'] for row in rows(request['predictor_files']['validation_cases'])] == ['val_keep'] * 2

    admission = root / 'preflight.json'
    command = [sys.executable, cli.__file__, 'fit', str(output), '--seed', '17',
               '--plan-sha256', digest, '--admission', str(admission)]
    write_json(admission, {'sample': {'fixture': True}, 'decisions': {'fixture': {'admit': True}},
                           'jobs': [{'id': 'fixture', 'kind': 'training', 'memory_gib': 18, 'command': command}]})
    with patch.multiple(native.torch.cuda, reset_peak_memory_stats=lambda device: None,
                        synchronize=lambda device: None), patch.object(native, 'load_recipe', return_value=config):
        fit = cli.fit_source(runtime, output, plan, digest, 17, admission)
    directory = cli.source_path(output, 17)
    assert cli.fit_source(runtime, output, plan, digest, 17, None) == fit
    assert not (directory / 'source.json').exists()
    source = cli.reload_source(runtime, output, plan, digest, 17)
    assert source == cli.load_source(runtime, output, digest, 17) == fit['source']
    history = read_json(directory / 'run.json')['history']
    best = min(history, key=lambda event: event['validation_loss'])
    assert (source['best_step'], source['best_validation_loss']) == (best['step'], best['validation_loss'])
    assert Path(source['path']).name == f"01-ema-selected-step-{best['step']:05}.pt"
    assert sha256(source['path']) == source['checkpoint_sha256']
    assert read_json(directory / 'replay.json')['prediction_max_abs_error'] == 0
    print('PASS fold 1 source membership, exact +0.010 promotion, no-run, native fit and independent replay')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime', type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as temp:
        fixture(args.runtime.resolve(), Path(temp))
