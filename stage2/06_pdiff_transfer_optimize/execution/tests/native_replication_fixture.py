#!/usr/bin/env python3
"""Fit and independently replay the promoted fold 1 target arms on a tiny CPU release."""
import argparse
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))

import replication as cli
from artifacts import read_json, sha256, write_csv, write_json
from formal_matrix import bind
from tahoe_bridge import dose_layout
from training_inputs import json_hash


def target_request(root, np):
    axis = root / 'genes.csv'
    axis.write_text('symbol\nA\nB\nC\n')
    files = {'genes': axis.name}
    roles = ('train', 'validation', 'test')
    for i, (role, dose) in enumerate(zip(roles, (0.2, 1, 10))):
        files[f'{role}_cases'] = f'{role}.csv'
        write_csv(root / files[f'{role}_cases'], [
            {'case_id': f'{role}-{well}', 'group_id': role, 'drug_key': role,
             'dose_um': dose, 'control_index': well} for well in range(2)])
        if role != 'test':
            files[f'{role}_expression'] = f'{role}.npy'
            np.save(root / files[f'{role}_expression'],
                    np.array([[.25, -.25, .5], [.5, -.5, .25]], dtype='float32') / (i + 1))
    files.update(controls='controls.npy', morgan='morgan.npy', drugs='drugs.csv')
    np.save(root / files['controls'], np.array([[.1, .2, .3], [.4, .5, .6]], dtype='float32'))
    np.save(root / files['morgan'], np.eye(3, 2048, dtype='float32'))
    write_csv(root / files['drugs'], [{'drug_key': role, 'morgan_index': i} for i, role in enumerate(roles)])
    hashes = {key: sha256(root / path) for key, path in files.items()}
    manifest = root / 'manifest.json'
    write_json(manifest, {'files': {f'fold_1/evaluator/{name}': 'unopened' for name in
                                   ('groups.csv', 'truth.npy', 'positive.npy', 'negative.npy')}})
    return {'identity': {'dataset': 'tenv', 'fold': 1, 'gene_sha256': hashes['genes'],
                         'source_input_sha256': sha256(manifest)},
            'input_root': str(root), 'input': {'predictor_files': files},
            'n_genes': 3, 'target_kind': 'delta', 'file_sha256': hashes,
            'predictor_files': {key: str(root / path) for key, path in files.items()}}


def source_checkpoint(native, root, request, config):
    inputs = native.PredictorInputs(request)
    source_cases = {role: [{**case, 'dose_um': dose} for case in inputs.cases[role]]
                    for role, dose in zip(inputs.cases, (0.05, 0.5, 5))}
    native.torch.set_num_threads(config['threads'])
    native.seed_randomness(23)
    source_inputs = SimpleNamespace(n_genes=inputs.n_genes, cases=source_cases,
                                    morgan=inputs.morgan, drug_index=inputs.drug_index)
    model = native.JointPerturbDiff(source_inputs, config, root)
    path = root / 'tahoe.pt'
    selection = {'identity': {'dataset': 'tahoe', 'seed': config['seed']}, 'recipe_sha256': json_hash(config),
                 'best_step': 1, 'best_validation_loss': .5, 'model_state_sha256': native.digest(model)}
    native.torch.save({**selection, 'format': native.FORMAT, 'model_state': model.state_dict()}, path)
    doses, indices = dose_layout(model.keys)
    return {**selection, 'path': str(path), 'checkpoint_sha256': sha256(path),
            'doses_um': doses, 'dose_indices': indices}


def fixture(runtime, output):
    cli.use_runtime(runtime)
    import native_training as native

    seed = 17
    config = {'stage': 'formal', 'seed': seed, 'fold': 1, 'threads': 2, 'device': 'cpu',
              'expression_scale': 1, 'chemical_dimension': 16, 'batch_groups': 1,
              'population_size': 2, 'learning_rate': .0002, 'weight_decay': 0,
              'sampling_start_time': 1, 'warmup_steps': 1, 'schedule_training_steps': 4,
              'ema_decay': .9, 'ema_update_steps': 1, 'max_steps': 3,
              'validation_interval': 1, 'minimum_steps': 2, 'patience_steps': 1,
              'max_training_seconds': 120, 'min_delta': 1e6}
    request = target_request(output, native.np)
    inference = cli.baseline.inference_request(request, role='test')
    source = source_checkpoint(native, output, request, config)
    source_path = cli.confirmation.source_path(output, seed) / 'source.json'
    write_json(source_path, source)
    recipient, recipe = output / 'request.json', output / 'recipe.json'
    write_json(recipient, request)
    write_json(recipe, config)
    upstream = {'status': 'ready', 'seeds': [seed], 'promoted_arms': list(cli.PRIMARY_MODES),
                'screens': {}, 'recipient': {'request': request, 'inference': inference,
                    'producer_request': bind(recipient), 'fixed_config': bind(recipe),
                    'roles': {'validation': {'group_order_sha256': json_hash(['validation'])}}},
                'comparison': {'common_steps': {str(seed): 3}, 'score': 'directional mean AP@100',
                               'prediction': {'seed_offset': 300000, 'max_groups': 1, 'max_profiles': 2}}}
    cm_root = output / 'cm'
    with (patch.multiple(native.torch.cuda, reset_peak_memory_stats=lambda device: None,
                         synchronize=lambda device: None),
          patch.object(native, 'load_recipe', return_value=config)):
        native.fit(runtime, request, recipe, cm_root)
        run = read_json(cm_root / 'run.json')
        replay = cm_root / 'replay.json'
        write_json(replay, run['reload'])
        cm = {'config': config, 'metadata': {key: run[key] for key in cli.baseline.CHECKPOINT_KEYS},
              'recipe': bind(recipe), 'request': bind(recipient), 'history': bind(cm_root / 'run.json'),
              'replay': bind(replay), 'checkpoint': bind(cm_root / 'best.pt')}
        with (patch.object(cli.confirmation, 'load_plan', return_value=upstream),
              patch.object(cli.confirmation, 'load_source', return_value=source),
              patch.object(cli, 'historical_scoring', return_value={'step04_checkpoints': bind(recipient)}),
              patch.object(cli, 'matched_cm', return_value=cm)):
            digest = cli.freeze(runtime, output, 'fixture-confirmation')['sha256']
            plan = cli.load_plan(runtime, output, digest)
            for arm in (cli.ORIGINAL, *cli.PRIMARY_MODES):
                admission = output / f'{arm}-admission.json'
                command = [sys.executable, cli.__file__, 'fit', str(output), '--seed', str(seed),
                           '--arm', arm, '--plan-sha256', digest, '--admission', str(admission)]
                write_json(admission, {'sample': {'fixture': True}, 'decisions': {'fixture': {'admit': True}},
                                       'jobs': [{'id': 'fixture', 'kind': 'training',
                                                 'memory_gib': 18, 'command': command}]})
                fit = cli.fit(runtime, output, plan, digest, seed, arm, admission)
                assert cli.fit(runtime, output, plan, digest, seed, arm, None) == fit
                assert fit['status'] == 'completed' and fit['reason'] == 'validation_early_stop'
                history = read_json(fit['history']['path'])
                assert [event['step'] for event in history['history']] == [1, 2]
                assert history['optimizer_steps'] == (2 if arm == 'dose_reset' else 3)

    choices, fits = cli.comparisons_for(output, plan, digest, seed)
    assert set(choices) == {'CM/native', 'Transfer/native', 'Transfer/common', 'dose_reset/native',
                            'learning_rate/native', 'learning_rate/common'}
    assert fits['dose_reset']['initialization']['dose_reset']['reset_doses_um'] == [0.2, 1., 10.]
    assert plan['seeds'][str(seed)]['arms']['learning_rate']['config'] == {**config, 'learning_rate': .0001}
    cli.predict(runtime, output, plan, digest, seed)
    for key, selected in choices.items():
        if key != 'CM/native':
            checkpoint = Path(selected['checkpoint']['path'])
            arm, mode = key.split('/')
            best = selected['metadata']['best_step']
            if mode == 'native':
                assert best in (1, 2)
                name = f'01-ema-selected-step-{best:05}.pt'
            else:
                assert best == 3
                name = '02-common-step-3.pt'
            assert checkpoint == cli.fit_binding(output, seed, arm) / name
            assert checkpoint.parents[1].name == '10-fold-1-target-seed-17'
        values = cli.baseline.verify_prediction(
            inference, selected['checkpoint'], selected['metadata'],
            cli.seed_slot(output, seed) / 'predictions' / key, cli.prediction_lineage(digest, seed, key, selected))
        assert values.shape == (2, 3) and native.np.isfinite(values).all()
    print('PASS promoted fold 1 target fits, numbered checkpoints and target-free native replay')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime', type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as temp:
        fixture(args.runtime.resolve(), Path(temp))
