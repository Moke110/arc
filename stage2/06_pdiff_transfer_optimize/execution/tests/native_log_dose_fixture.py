#!/usr/bin/env python3
"""Exercise continuous dose through native covariates and selected-state transfer."""
import argparse
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))

from baseline import sha256
from formal import use_runtime


def fixture(runtime, root):
    use_runtime(runtime)
    from log_dose import install_encoder, load_initialization
    from native_training import FORMAT, JointPerturbDiff, digest, np, seed_randomness, torch
    from tahoe_bridge import CHEMICAL_LOOKUP, DOSE_LOOKUP, dose_layout
    from training_inputs import json_hash

    def model(doses, fingerprint_index):
        cases = [{'drug_key': 'drug', 'dose_um': str(dose)} for dose in doses]
        morgan = np.eye(2, 2048, dtype=np.float32)
        inputs = SimpleNamespace(n_genes=3, cases={'train': cases},
                                 morgan=morgan, drug_index={'drug': fingerprint_index})
        config = {'chemical_dimension': 16, 'learning_rate': 0.0002, 'weight_decay': 0,
                  'expression_scale': 1, 'batch_groups': 2}
        return install_encoder(JointPerturbDiff(inputs, config, root))

    def doses_through_native(model, values):
        encoder = model.native.cov_encoder
        cov_pert = torch.tensor([[model.lookup[('drug', dose)] if dose else len(model.keys)]
                                 for dose in values])
        observed = []
        with encoder.dose_encoder.register_forward_hook(lambda _mod, _args, out: observed.append(out)):
            model.native._encode_covariates({'cov_pert': cov_pert, 'cov_celltype': torch.zeros_like(cov_pert),
                                             'cov_batch': torch.zeros_like(cov_pert)})
        return observed[0][:, 0]

    torch.set_num_threads(2)
    seed_randomness(17)
    source = model((0.05, 0.5, 1, 5), 0)
    dose_encoder = source.native.cov_encoder.dose_encoder
    with torch.no_grad():
        dose_encoder.projection.weight.fill_(2)
        dose_encoder.projection.bias.fill_(3)
        dose_encoder.control.fill_(7)
    measured = doses_through_native(source, [0, 0.05, 0.5, 5])
    torch.testing.assert_close(measured[:, 0], torch.tensor([7., 3., 5., 7.]))
    measured[0].sum().backward()
    assert torch.all(dose_encoder.control.grad == 1)
    assert torch.count_nonzero(dose_encoder.projection.weight.grad) == 0

    checkpoint = root / 'source.pt'
    recipe = {'seed': 17, 'stage': 'formal'}
    selected = {'identity': {'dataset': 'tahoe', 'seed': 17, 'fold': 0, 'stage': 'formal',
                             'arm': 'log', 'regime': 'cold-cm-drug'},
                'recipe_sha256': json_hash(recipe), 'best_step': 2, 'best_validation_loss': 0.5,
                'model_state_sha256': digest(source)}
    torch.save({**selected, 'format': FORMAT, 'recipe': recipe,
                'model_state': source.state_dict()}, checkpoint)
    metadata = {**selected, 'checkpoint_sha256': sha256(checkpoint)}
    source_doses, source_indices = dose_layout(source.keys)
    binding = {'metadata': metadata, 'checkpoint': {'path': str(checkpoint), 'sha256': sha256(checkpoint)},
               'doses_um': source_doses, 'dose_indices': source_indices}
    recipient = model((0.2, 1, 10), 1)
    request = {'identity': {'dataset': 'tenv', 'fold': 0}}
    original_lookups = {key: recipient.state_dict()[key].clone() for key in (CHEMICAL_LOOKUP, DOSE_LOOKUP)}
    load_initialization(recipient, binding, request)
    assert all(torch.equal(recipient.state_dict()[key], original_lookups[key]) for key in original_lookups)
    assert all(torch.equal(recipient.state_dict()[key], value)
               for key, value in source.state_dict().items() if key not in original_lookups)
    measured = doses_through_native(recipient, [0, 0.2, 1, 10])
    torch.testing.assert_close(measured[:, 0], torch.tensor([7., 4.20412, 5.60206, 7.60206]))
    reloaded = model((0.2, 1, 10), 1)
    reloaded.load_state_dict(recipient.state_dict(), strict=True)
    assert torch.equal(doses_through_native(reloaded, [0, 0.2, 1, 10]), measured)
    print('PASS native log covariates, source transfer, recipient lookups and exact reload', flush=True)



def fit_source_fixture(runtime, root):
    from unittest.mock import patch

    import native_training as native
    import paired_log_source as worker
    from artifacts import read_json, write_csv, write_json
    from formal_matrix import bind
    from paired_reference import REGIME, SCHEMA
    from training_inputs import json_hash

    root.mkdir()
    files, roles = {}, {}
    for role, dose in [('train', 0.05), ('validation', 0.5), ('test', 5)]:
        cases = root / f'{role}.csv'
        write_csv(cases, [{'case_id': f'{role}-{i}', 'group_id': role, 'drug_key': role,
                           'dose_um': dose, 'control_index': i} for i in range(2)])
        files[f'{role}_cases'] = str(cases)
        roles[role] = {'parent_connectivity': [role], 'case_count': 2, 'group_count': 1,
                       'cases_sha256': sha256(cases)}
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
    write_csv(drugs, [{'drug_key': role, 'parent_connectivity': role, 'morgan_index': i}
                      for i, role in enumerate(roles)])
    files['drugs'] = str(drugs)
    request_identity = {'dataset': 'tahoe', 'fold': 0, 'gene_sha256': 'fixture'}
    request = {'identity': request_identity, 'target_kind': 'delta', 'n_genes': 3, 'predictor_files': files,
               'file_sha256': {key: sha256(path) for key, path in files.items()}}
    config = {'stage': 'formal', 'seed': 17, 'fold': 0, 'threads': 2, 'device': 'cpu',
              'expression_scale': 1, 'chemical_dimension': 16, 'batch_groups': 1,
              'population_size': 2, 'learning_rate': 0.0002, 'weight_decay': 0,
              'sampling_start_time': 1, 'warmup_steps': 1, 'schedule_training_steps': 4,
              'ema_decay': 0.9, 'ema_update_steps': 1, 'max_steps': 4,
              'validation_interval': 1, 'minimum_steps': 2, 'patience_steps': 1,
              'max_training_seconds': 120, 'min_delta': 1e6}
    selected = {'identity': {**request_identity, 'seed': 17, 'stage': 'formal'},
                'recipe_sha256': json_hash(config), 'best_step': 1, 'best_validation_loss': 0.5,
                'checkpoint_sha256': 'original-checkpoint', 'model_state_sha256': 'original-state'}
    original = {'roles': roles, 'metadata': selected, 'source_sha256': 'original-release'}
    for key, value in [('request', request), ('recipe', {}),
                       ('history', {**selected, 'recipe': config})]:
        path = root / f'{key}.json'
        write_json(path, value)
        original[key] = bind(path)
    reference = root / 'reference'
    reference.mkdir()
    plan = {'schema': SCHEMA, 'fold': 0, 'seed': 17, 'regime': REGIME, 'source': original,
            'request': {'role': 'validation'}, 'protocol_sha256': 'fixture'}
    write_json(reference / 'plan.json', plan)
    digest = sha256(reference / 'plan.json')
    target = root / 'output' / 'paired-transfer' / 'log' / 'seed-17'
    args = [str(target), '--runtime', str(runtime), '--reference', str(reference),
            '--seed', '17', '--plan-sha256', digest]
    with patch.multiple(native.torch.cuda, reset_peak_memory_stats=lambda device: None,
                        synchronize=lambda device: None), patch.object(native, 'load_recipe', return_value=config), \
            patch.object(worker, 'verify_request', return_value=None):
        worker.main(['fit', *args])
    worker.main(['verify', *args])
    evidence = worker.selected_source(target, plan, digest)
    assert Path(evidence['fit_seal']['path']).parent.parent == worker.source_slot(target, 17)
    run = read_json(evidence['history']['path'])
    assert (run['optimizer_steps'], run['stop_reason']) == (2, 'validation_early_stop')
    assert evidence['metadata']['best_step'] == min(run['history'], key=lambda event: event['validation_loss'])['step']
    assert run['reload']['prediction_max_abs_error'] == 0
    assert evidence['metadata']['model_state_sha256'] != selected['model_state_sha256']
    assert evidence['doses_um'] == [0.05, 0.5, 5.]
    print('PASS bounded log Tahoe source fit, author-loss selection and exact prediction reload', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime', type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as temp:
        fixture(args.runtime.resolve(), Path(temp))
        fit_source_fixture(args.runtime.resolve(), Path(temp) / 'source-fit')
