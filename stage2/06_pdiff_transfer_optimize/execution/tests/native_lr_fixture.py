#!/usr/bin/env python3
"""Bounded CPU exercise of the native trainer's selection and diagnostic snapshot."""
import argparse
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))

from baseline import write_csv, write_json, read_json, sha256
from formal import use_runtime
from native_lr import fit_diagnostic


def fixture(runtime, root):
    use_runtime(runtime)
    import native_training as native

    files = {}
    for role, count in [('train', 2), ('validation', 2), ('test', 1)]:
        cases = root / f'{role}.csv'
        write_csv(cases, [{'case_id': f'{role}-{i}', 'group_id': role, 'drug_key': role,
                           'dose_um': '1', 'control_index': i % 2} for i in range(count)])
        files[f'{role}_cases'] = str(cases)
        if role != 'test':
            targets = root / f'{role}.npy'
            native.np.save(targets, native.np.array([[.25, -.25, .5], [.5, -.5, .25]], dtype='float32'))
            files[f'{role}_expression'] = str(targets)
    for name, values in [('controls', [[.1, .2, .3], [.4, .5, .6]]),
                         ('morgan', native.np.eye(3, 2048))]:
        path = root / f'{name}.npy'
        native.np.save(path, native.np.array(values, dtype='float32'))
        files[name] = str(path)
    drugs = root / 'drugs.csv'
    write_csv(drugs, [{'drug_key': role, 'morgan_index': i}
                      for i, role in enumerate(('train', 'validation', 'test'))])
    files['drugs'] = str(drugs)
    request = {'identity': {'dataset': 'tenv', 'fold': 0}, 'target_kind': 'delta',
               'n_genes': 3, 'predictor_files': files}
    config = {'stage': 'formal', 'seed': 17, 'fold': 0, 'threads': 2, 'device': 'cpu',
              'expression_scale': 1, 'chemical_dimension': 16, 'batch_groups': 1,
              'population_size': 2, 'learning_rate': 0.0001, 'weight_decay': 0,
              'sampling_start_time': 1, 'warmup_steps': 1, 'schedule_training_steps': 4,
              'ema_decay': 0.9, 'ema_update_steps': 1, 'max_steps': 4,
              'validation_interval': 1, 'minimum_steps': 2, 'patience_steps': 1,
              'max_training_seconds': 120, 'min_delta': 1e6}
    recipe = root / 'recipe.json'
    write_json(recipe, {})
    with patch.multiple(native.torch.cuda, reset_peak_memory_stats=lambda device: None,
                        synchronize=lambda device: None), patch.object(native, 'load_recipe', return_value=config):
        native.fit(runtime, request, recipe, root / 'early_reference')
        fit_diagnostic(native, runtime, request, recipe, root / 'early', 3)
    reference = read_json(root / 'early_reference/run.json')
    early = read_json(root / 'early/run.json')
    seal = read_json(root / 'early/native-selection.json')
    assert (seal['stop_step'], seal['stop_reason'], early['optimizer_steps']) == (2, 'validation_early_stop', 3)
    assert early['best_step'] == reference['best_step']
    assert early['model_state_sha256'] == reference['model_state_sha256']
    assert early['reload'] == reference['reload']
    assert [event['step'] for event in early['history']] == [1, 2]
    assert early['history'][-1]['early_stop']
    diagnostic = read_json(root / 'early/diagnostic.json')
    assert diagnostic['best_step'] == 3
    assert diagnostic['checkpoint_sha256'] == sha256(root / 'early/02-common-step-3.pt')

    config = {**config, 'minimum_steps': 5, 'min_delta': 0}
    with patch.multiple(native.torch.cuda, reset_peak_memory_stats=lambda device: None,
                        synchronize=lambda device: None), patch.object(native, 'load_recipe', return_value=config):
        native.fit(runtime, request, recipe, root / 'late_reference')
        fit_diagnostic(native, runtime, request, recipe, root / 'late', 1)
    reference = read_json(root / 'late_reference/run.json')
    late = read_json(root / 'late/run.json')
    seal = read_json(root / 'late/native-selection.json')
    assert late['best_step'] == reference['best_step'] > 1
    assert late['model_state_sha256'] == reference['model_state_sha256']
    assert late['reload'] == reference['reload']
    assert (seal['stop_step'], seal['stop_reason'], late['optimizer_steps']) == (4, 'max_steps', 4)
    print('PASS early stop before diagnostic and native selection after diagnostic', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime', type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as temp:
        fixture(args.runtime.resolve(), Path(temp))
