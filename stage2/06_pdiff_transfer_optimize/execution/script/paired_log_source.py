#!/usr/bin/env python3
"""Fit or verify the ADR 0049 log-dose Tahoe source for one paired seed."""
import argparse
from pathlib import Path
from unittest.mock import patch

import baseline as base
from artifacts import SEEDS, read_json
from formal import completed_fit, fit_slot, slot_lock, use_runtime, verify_fit_seal
from formal_matrix import bind, check_binding
from paired_reference import REGIME, TRANSFER, load_plan
from tahoe_bridge import dose_layout
from training_inputs import json_hash, verify_request

SCHEMA = 'step06-paired-log-source-v1'


def source_inputs(plan, reference_digest):
    original = plan['source']
    request = read_json(check_binding(original['request'])['path'])
    recipe = Path(check_binding(original['recipe'])['path'])
    history = read_json(check_binding(original['history'])['path'])
    config = history['recipe']
    identity = {**request['identity'], 'seed': plan['seed'], 'stage': 'formal'}
    if ((identity['dataset'], identity['fold'], plan['regime']) != ('tahoe', 0, REGIME)
            or original['metadata']['identity'] != identity
            or any(history[key] != original['metadata'][key] for key in base.CHECKPOINT_KEYS)
            or json_hash(config) != original['metadata']['recipe_sha256']):
        raise ValueError('original Tahoe history, recipe, or request differs from the paired reference')
    identity = {**identity, 'arm': 'log', 'regime': REGIME}
    for role, evidence in original['roles'].items():
        key = f'{role}_cases'
        bind(request['predictor_files'][key], evidence['cases_sha256'])
    bind(request['predictor_files']['drugs'], request['file_sha256']['drugs'])
    cases = TRANSFER.tahoe.role_cases(request)
    parents = TRANSFER.tahoe.parent_blocks(request, cases)
    groups = {role: list(dict.fromkeys(row['group_id'] for row in members))
              for role, members in cases.items()}
    roles = {role: {'parent_connectivity': parents[role], 'case_count': len(members),
                    'group_count': len(groups[role]), 'cases_sha256': request['file_sha256'][f'{role}_cases']}
             for role, members in cases.items()}
    if roles != original['roles']:
        raise ValueError('log source case roles differ from the original Tahoe exposure')
    doses, indices = dose_layout(TRANSFER.tahoe.key_layout(cases))
    protocol = {
        'schema': SCHEMA, 'reference_plan_sha256': reference_digest,
        'protocol_sha256': plan['protocol_sha256'], 'identity': identity,
        'original_source_release_sha256': original['source_sha256'],
        'original_request': original['request'], 'original_recipe': original['recipe'],
        'original_history': original['history'], 'recipe_sha256': json_hash(config),
        'roles': roles, 'doses_um': doses, 'dose_indices': indices,
        'validation_groups_sha256': json_hash(groups['validation']),
    }
    return request, recipe, config, protocol


def selected_source(root, plan, reference_digest):
    """Verify the source loss minimum and bind its selected checkpoint for tenv."""
    _, _, config, protocol = source_inputs(plan, reference_digest)
    directory = root / 'source' / 'fit'
    digest = json_hash(protocol)
    verify_fit_seal(directory, digest)
    run = completed_fit(directory, protocol['identity'], config)
    last = run['history'][-1]
    batching = run['validation_batching']
    if (run['optimizer_steps'] != last['step']
            or (run['stop_reason'] != 'validation_early_stop' if last['early_stop']
                else run['stop_reason'] != 'max_steps' or last['step'] != config['max_steps'])
            or batching['group_ids_sha256'] != protocol['validation_groups_sha256']
            or json_hash(batching['group_ids']) != protocol['validation_groups_sha256']):
        raise ValueError('log source stopping or validation groups differ from the frozen protocol')
    return {
        **protocol, 'source_plan_sha256': digest, 'fit_seal': bind(directory / 'fit.json'),
        'checkpoint': bind(directory / 'best.pt', run['checkpoint_sha256']),
        'history': bind(directory / 'run.json'),
        'metadata': {key: run[key] for key in base.CHECKPOINT_KEYS},
        'optimizer_steps': run['optimizer_steps'],
    }


def fit_source(runtime, root, plan, reference_digest):
    use_runtime(runtime)
    import native_training as native
    from log_dose import install_encoder

    request, recipe, config, protocol = source_inputs(plan, reference_digest)
    verify_request(request, runtime)
    if native.load_recipe(runtime, recipe, request) != config:
        raise ValueError('resolved log source recipe differs from the original Tahoe recipe')
    worker_request = {**request, 'identity': {**request['identity'], 'arm': 'log', 'regime': REGIME}}

    class LogTahoe(native.JointPerturbDiff):
        def __init__(self, inputs, recipe, directory):
            super().__init__(inputs, recipe, directory)
            install_encoder(self)

    slot = root / 'source'
    with slot_lock(slot), patch.object(native, 'JointPerturbDiff', LogTahoe):
        fit_slot(slot, json_hash(protocol), runtime, worker_request, protocol['identity'], recipe, config)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('fit', 'verify'))
    parser.add_argument('root', type=Path)
    parser.add_argument('--runtime', type=Path, default=base.REPO)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--seed', type=int, choices=SEEDS, required=True)
    parser.add_argument('--plan-sha256', required=True)
    args = parser.parse_args(argv)
    runtime, root = args.runtime.resolve(), args.root.resolve()
    plan = load_plan(args.reference.resolve(), args.seed, args.plan_sha256)
    if args.phase == 'fit':
        fit_source(runtime, root, plan, args.plan_sha256)
    evidence = selected_source(root, plan, args.plan_sha256)
    result = evidence['fit_seal']
    print(f'{result["path"]} SHA-256 {result["sha256"]}', flush=True)


if __name__ == '__main__':
    main()
