#!/usr/bin/env python3
"""Run target-free validation replay through the pinned PerturbDiff interpreter."""
import argparse
from pathlib import Path
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))

from baseline import inference_request, predict_checkpoint, score_prediction, verify_prediction
from artifacts import sha256, write_csv
from formal import use_runtime
from training_inputs import json_hash


def replay(runtime, root):
    use_runtime(runtime)
    from native_training import FORMAT, JointPerturbDiff, digest, np, seed_randomness, torch
    from types import SimpleNamespace

    axis = root / 'genes.csv'
    axis.write_text('symbol\nA\nB\nC\n')
    files = {'genes': str(axis)}
    inventories = {}
    for index, role in enumerate(('train', 'validation', 'test')):
        drug = role
        count = index + 1
        inventories[role] = [
            {'case_id': f'{role}-{well}', 'drug_key': drug, 'dose_um': '1',
             'control_index': str(well % 2), 'group_id': role} for well in range(count)
        ]
        files[f'{role}_cases'] = str(root / f'{role}.csv')
        write_csv(files[f'{role}_cases'], inventories[role])

    files.update({key: str(root / name) for key, name in
                  (('controls', 'controls.npy'), ('morgan', 'morgan.npy'), ('drugs', 'drugs.csv'))})
    controls = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32)
    morgan = np.eye(3, 2048, dtype=np.float32)
    np.save(files['controls'], controls)
    np.save(files['morgan'], morgan)
    drug_index = {role: index for index, role in enumerate(inventories)}
    write_csv(files['drugs'], [{'drug_key': role, 'morgan_index': index} for role, index in drug_index.items()])
    identity = {'dataset': 'tenv', 'fold': 0, 'gene_sha256': sha256(axis)}
    request = inference_request({'identity': identity, 'n_genes': 3, 'target_kind': 'delta',
                                 'file_sha256': {key: sha256(path) for key, path in files.items()},
                                 'predictor_files': {**files, 'validation_expression': str(root / 'absent.npy')}})
    config = {'seed': 17, 'threads': 2, 'device': 'cpu', 'expression_scale': 1,
              'chemical_dimension': 16, 'batch_groups': 2, 'learning_rate': 0.0002,
              'weight_decay': 0, 'sampling_start_time': 1}
    settings = {'seed_offset': 300000, 'max_groups': 2, 'max_profiles': 2}
    torch.set_num_threads(2)
    seed_randomness(config['seed'])
    inputs = SimpleNamespace(cases=inventories, n_genes=3, morgan=morgan, drug_index=drug_index)
    model = JointPerturbDiff(inputs, config, root)
    state = {'format': FORMAT, 'identity': {**identity, 'seed': 17, 'stage': 'fixture'},
             'recipe_sha256': json_hash(config), 'best_step': 1, 'best_validation_loss': 0.5,
             'model_state_sha256': digest(model), 'model_state': model.state_dict()}
    checkpoint = root / 'ema.pt'
    torch.save(state, checkpoint)
    metadata = {key: value for key, value in state.items() if key not in ('format', 'model_state')}
    metadata['checkpoint_sha256'] = sha256(checkpoint)
    binding = {'path': str(checkpoint), 'sha256': sha256(checkpoint)}
    lineage = {'arm': 'fixture'}
    for name in ('first', 'second'):
        predict_checkpoint(request, config, settings, binding, metadata, root / name, lineage)
    assert sha256(root / 'first/prediction.npy') == sha256(root / 'second/prediction.npy')
    wrong_identity = {**metadata, 'identity': {**metadata['identity'], 'dataset': 'tahoe'}}
    try:
        predict_checkpoint(request, config, settings, binding, wrong_identity, root / 'mismatch', lineage)
    except ValueError as error:
        assert 'case/gene identity' in str(error)
    else:
        raise AssertionError('mismatched checkpoint identity was accepted')
    values = verify_prediction(request, binding, metadata, root / 'first', lineage)
    truth = np.array([[1, -1, 0], [1, -1, 0]], dtype=np.float32)
    group, summary = score_prediction(request, truth, values)
    assert group[0]['n_wells'] == 2 and summary[0]['group_count'] == 1
    assert summary[0]['mean_ap100'] is not None
    print('PASS target-free native validation replay and scoring fixture', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime', type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as temp:
        replay(args.runtime.resolve(), Path(temp))
