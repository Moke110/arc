#!/usr/bin/env python3
"""Verify frozen Tahoe EMA states and score zero-shot predictions on tenv."""
import argparse
from contextlib import contextmanager
import csv
import hashlib
from pathlib import Path
import subprocess
import tempfile
import time

from warm_tenv import REPO, STEP04, N_GENES, ROLES, expression_sources, grouped_truth, release, verified_plan
from cold_tenv import pinned, unique_rows, verify_split
from evaluate import (GROUP_FIELDS, RECALL_GRID, RUN_FIELDS, SEEDS, evaluator_inputs, fold_inputs,
                      grouped_cases, load_array, load_schema, np, read_json, rows, run_scores,
                      sha256, summaries, verify_file, write_csv, write_json)
from formal import CHECKPOINT_FIELDS, completed_fit, slot_lock, use_runtime, verify_fit_seal
from tahoe_bridge import (FORMAT, LIMITATION, dose_layout, load_initialization, request_digest,
                          seal_initialization)
from training_inputs import json_hash, verify_request

SOURCE_SCHEMA = 'step05-tahoe-source-v1'
PREDICTION_SCHEMA = 'step05-tahoe-heldout-prediction-v1'
REGIMES = ('warm-start', 'cold-cm-drug')
STEP05 = Path('stage2/05_pdiff_transfer')
TAHOE_PRODUCER = 'd2ae285feb4138c0716105919aac9ca06c6ebc3e'


def historical_record(dataset):
    return REPO / STEP04 / 'execution/formal_record' / dataset


def role_cases(request):
    return {role: rows(request['predictor_files'][f'{role}_cases']) for role in ROLES}


def key_layout(cases):
    return sorted({(row['drug_key'], float(row['dose_um'])) for members in cases.values() for row in members})


def parent_blocks(request, cases):
    parents = {row['drug_key']: row['parent_connectivity']
               for row in rows(request['predictor_files']['drugs'])}
    return {role: sorted({parents[row['drug_key']] for row in members}) for role, members in cases.items()}


def source_slot(formal, record, plan_sha, entry, request, seed, outcome, validation_groups):
    recipe_entry = entry['recipes'][str(seed)]
    recipe = read_json(pinned(formal, recipe_entry['path'], recipe_entry['sha256']))
    config_sha = entry['recipe_sha256_by_seed'][str(seed)]
    slot = formal / 'slots' / f'fold-0-seed-{seed}'
    identity = {**request['identity'], 'seed': seed, 'stage': 'formal'}
    if ((outcome['status'], outcome['dataset'], outcome['seed']) != ('completed', 'tahoe', seed)
            or outcome['prediction'] != str((slot / 'prediction/prediction.json').relative_to(formal))):
        raise ValueError('unexpected historical Tahoe prediction outcome')
    compact = record / 'slots/fold-0' / f'seed-{seed}'
    manifest = read_json(pinned(slot / 'prediction', 'prediction.json', outcome['prediction_manifest_sha256']))
    verify_file(compact / 'prediction.json', outcome['prediction_manifest_sha256'])
    verify_fit_seal(slot / 'fit', plan_sha)
    verify_file(compact / 'run.json', manifest['run_sha256'])
    pinned(slot / 'fit', 'run.json', manifest['run_sha256'])
    run = read_json(compact / 'run.json')
    completed_fit(slot / 'fit', identity, run['recipe'])
    replay = read_json(pinned(slot / 'prediction', 'replay.json', manifest['replay_sha256']))
    verify_file(compact / 'replay.json', manifest['replay_sha256'])
    selected = {key: run[key] for key in CHECKPOINT_FIELDS}
    if (any(run['recipe'][key] != value for key, value in recipe.items())
            or run['recipe_sha256'] != config_sha
            or run['validation_batching']['group_ids'] != validation_groups
            or run['validation_batching']['group_ids_sha256'] != json_hash(validation_groups)
            or manifest['plan_sha256'] != plan_sha
            or any(manifest[key] != value for key, value in selected.items())
            or replay != {**run['reload'], **selected, 'plan_sha256': plan_sha}):
        raise ValueError('Tahoe checkpoint, recipe, validation membership or replay identity mismatch')
    checkpoint = pinned(slot / 'fit', 'best.pt', run['checkpoint_sha256'])
    if not (slot / 'prediction/best.pt').samefile(checkpoint):
        pinned(slot / 'prediction', 'best.pt', run['checkpoint_sha256'])
    return {**selected, 'path': str(checkpoint),
            'producer': {'plan_sha256': plan_sha, 'request_sha256': entry['request_sha256'],
                         'raw_recipe_sha256': recipe_entry['sha256'], 'run_sha256': manifest['run_sha256'],
                         'replay_sha256': manifest['replay_sha256'],
                         'prediction_manifest_sha256': outcome['prediction_manifest_sha256']}}


def verify_source(runtime, formal, output):
    if output.exists() or output.is_symlink():
        raise ValueError('source verification requires a fresh output path')
    record = historical_record('tahoe')
    plan_sha = sha256(record / 'plan.json')
    plan = load_schema(pinned(formal, 'plan.json', plan_sha), 'step04-evaluation-plan-v1')
    outcomes_sha = sha256(record / 'outcomes.json')
    outcomes = load_schema(pinned(formal, 'outcomes.json', outcomes_sha), 'step04-evaluation-outcomes-v1')
    if plan['dataset'] != 'tahoe' or outcomes['plan_sha256'] != plan_sha:
        raise ValueError('expected frozen Step 04 Tahoe evidence')
    pinned(formal, 'selection.json', plan['selection_sha256'])
    # The historical temporary checkout is gone. Its Step 04 files match the
    # producer revision, while its shared adapter and environment still match live files.
    for filename, digest in plan['launch_files'].items():
        parts = Path(filename).parts
        prefix = 'stage2' if 'stage2' in parts else 'pert_model'
        relative = Path(*parts[parts.index(prefix):])
        if relative.is_relative_to(STEP04):
            original = subprocess.check_output(['git', 'show', f'{TAHOE_PRODUCER}:{relative}'], cwd=runtime)
            if hashlib.sha256(original).hexdigest() != digest:
                raise ValueError(f'historical Tahoe producer checksum mismatch: {relative}')
        else:
            verify_file(runtime / relative, digest)
    entry = plan['folds'][0]
    request = read_json(pinned(formal, entry['request'], entry['request_sha256']))
    verify_request(request, runtime)
    if (request['identity']['dataset'], request['identity']['fold'], request['n_genes']) != ('tahoe', 0, N_GENES):
        raise ValueError('expected the full-axis Tahoe fold 0 request')
    cases = role_cases(request)
    validation_groups = list(dict.fromkeys(row['group_id'] for row in cases['validation']))
    doses, dose_indices = dose_layout(key_layout(cases))
    if doses != sorted({float(row['dose_um']) for row in cases['train']}):
        raise ValueError('the full Tahoe dose vocabulary must be present in training')
    slots = unique_rows([row for row in outcomes['outcomes'] if row['fold'] == 0], 'seed', 'Tahoe outcome')
    if set(slots) != set(SEEDS):
        raise ValueError('expected three historical Tahoe checkpoints')
    checkpoints = {}
    for seed in SEEDS:
        checkpoint = source_slot(formal, record, plan_sha, entry, request, seed, slots[seed], validation_groups)
        checkpoints[str(seed)] = {**checkpoint, 'doses_um': doses, 'dose_indices': dose_indices}
    parents = parent_blocks(request, cases)
    write_json(output, {'schema': SOURCE_SCHEMA, 'plan_sha256': plan_sha, 'outcomes_sha256': outcomes_sha,
                        'producer_revision': TAHOE_PRODUCER, 'request_identity': request['identity'],
                        'request_file_sha256': entry['request_sha256'],
                        'roles': {role: {'parent_connectivity': parents[role],
                                         'case_count': len(cases[role]),
                                         'group_count': len(grouped_cases(cases[role])[0]),
                                         'cases_sha256': request['file_sha256'][f'{role}_cases']}
                                  for role in ROLES},
                        'checkpoints': checkpoints})
    print(f'Verified Tahoe fold 0 source SHA-256 {sha256(output)}', flush=True)


def recipient(runtime, root, regime, digest):
    base = read_json(historical_record('tenv') / 'plan.json')
    fold, warm_source = None, None
    if regime == 'warm-start':
        plan, request = verified_plan(root, digest, runtime)
        inference = plan['inference']
        warm_source = Path(plan['source_root'])
    else:
        entry = base['folds'][0]
        record = read_json(REPO / STEP05 / 'execution/cold_record/replay.json')
        if digest != record['local_artifact_sha256']['cold-tenv-ticket76/split.json']:
            raise ValueError('cold tenv split differs from the verified replay')
        split = read_json(pinned(root, 'split.json', digest))
        request = read_json(pinned(root, 'predictor_request.json', entry['request_sha256']))
        verify_request(request, runtime)
        # The cold replay persists a byte-identical Step 04 request at a different name.
        fold = fold_inputs(root, {**entry, 'request': 'predictor_request.json'}, 'formal')
        if split != verify_split(request, fold):
            raise ValueError('unexpected cold tenv membership')
        inference = base['inference']
    return request, inference, fold, warm_source


def overlap(source, request):
    parents = parent_blocks(request, role_cases(request))
    return {source_role: {target_role: sorted(set(source['roles'][source_role]['parent_connectivity'])
                                              .intersection(parents[target_role]))
                          for target_role in ROLES}
            for source_role in ROLES}


def verify_pair(source, request):
    identity = source['request_identity']
    if (set(source['checkpoints']) != set(map(str, SEEDS))
            or (identity['dataset'], identity['fold']) != ('tahoe', 0)
            or (request['identity']['dataset'], request['identity']['fold']) != ('tenv', 0)
            or request['identity']['gene_sha256'] != identity['gene_sha256']
            or request['n_genes'] != N_GENES):
        raise ValueError('expected matching full-axis Tahoe and tenv fold 0 inputs')


def sample_cases(model, inputs, config, inference, directory, device):
    from native_training import groups_for, torch, validation_random_stream

    cases = inputs.cases['test']
    buckets = {}
    for group in groups_for(cases):
        buckets.setdefault(len(group), []).append(group)
    seen = np.zeros(len(cases), dtype=bool)
    predictions = np.lib.format.open_memmap(directory / 'prediction.npy', mode='w+', dtype='float32',
                                           shape=(len(cases), inputs.n_genes))
    inference_seed = config['seed'] + inference['seed_offset']
    generator = torch.Generator(device=device).manual_seed(inference_seed)
    with torch.no_grad(), validation_random_stream(inference_seed, device):
        for size, groups in buckets.items():
            batch_size = min(inference['max_groups'], max(1, inference['max_profiles'] // size))
            for start in range(0, len(groups), batch_size):
                chosen = groups[start:start + batch_size]
                batch, indices = model.batch(inputs, 'test', chosen, device)
                if 'pert_emb' in batch or not np.array_equal(indices, chosen):
                    raise ValueError('prediction changed frozen case membership')
                noise = torch.randn(batch['cont_emb'].shape, device=device, generator=generator)
                delta = (model.sample(batch, config, noise) - batch['cont_emb']) * config['expression_scale']
                values = delta.cpu().numpy()
                indices = indices.reshape(-1)
                if (values.shape != (len(chosen), size, inputs.n_genes) or values.dtype != np.float32
                        or not np.isfinite(values).all() or seen[indices].any()):
                    raise ValueError('invalid or repeated tenv prediction rows')
                predictions[indices] = values.reshape(-1, inputs.n_genes)
                seen[indices] = True
    if not seen.all():
        raise ValueError('prediction did not cover every frozen tenv case')
    predictions.flush()
    write_json(directory / 'cases.json', [[row['case_id'], row['group_id']] for row in cases])
    return {'case_count': len(cases), 'group_count': sum(map(len, buckets.values()))}


def predict(runtime, source, source_sha, root, regime, recipient_sha, seed, destination):
    request, inference, _, _ = recipient(runtime, root, regime, recipient_sha)
    verify_pair(source, request)
    predict_recipient(runtime, source, source_sha, request, inference, regime, recipient_sha, seed, destination)


def predict_recipient(runtime, source, source_sha, request, inference, regime, recipient_sha, seed, destination):
    checkpoint = source['checkpoints'][str(seed)]
    use_runtime(runtime)
    from native_training import JointPerturbDiff, PredictorInputs, resources, seed_randomness, torch
    from warm_inputs import WarmInputs

    # Sampling depends on the fixed recipe and seed, never on target outcomes.
    formal_recipe = read_json(Path(checkpoint['path']).parent / 'run.json')['recipe']
    if json_hash(formal_recipe) != checkpoint['recipe_sha256']:
        raise ValueError('source inference recipe differs from verified training recipe')
    destination.parent.mkdir(parents=True, exist_ok=True)
    with slot_lock(destination.parent / f'.{destination.name}-lock'):
        if destination.exists() or destination.is_symlink():
            verify_prediction(destination, source, source_sha, request, regime, recipient_sha, seed)
            return
        staging = Path(tempfile.mkdtemp(prefix=f'{destination.name}-attempt-', dir=destination.parent))
        print(f'Prediction attempt evidence: {staging}', flush=True)
        started = time.monotonic()
        torch.set_num_threads(formal_recipe['threads'])
        seed_randomness(seed)
        device = torch.device(formal_recipe['device'])
        torch.cuda.reset_peak_memory_stats(device)
        inputs = (WarmInputs if regime == 'warm-start' else PredictorInputs)(request)
        model = JointPerturbDiff(inputs, formal_recipe, staging)
        bridge_sha = seal_initialization(model, staging, checkpoint, request)
        bridge = load_initialization(model, staging, bridge_sha, request)
        model.to(device)
        model.native.to(device)
        model.eval()
        counts = sample_cases(model, inputs, formal_recipe, inference, staging, device)
        write_json(staging / 'prediction.json', {
            'schema': PREDICTION_SCHEMA, 'arm': 'Tahoe', 'regime': regime, 'seed': seed, **counts,
            'source_sha256': source_sha, 'recipient_sha256': recipient_sha,
            'recipient_identity': request['identity'], 'recipient_request_sha256': request_digest(request),
            'source_parent_overlap': overlap(source, request),
            'bridge_sha256': bridge_sha, 'mapped_state_sha256': bridge['mapped_state_sha256'],
            'dose_limitation': LIMITATION, 'cases_sha256': sha256(staging / 'cases.json'),
            'prediction_sha256': sha256(staging / 'prediction.npy'),
            'wall_seconds': time.monotonic() - started, 'resources': resources(device),
        })
        staging.rename(destination)


def verify_prediction(directory, source, source_sha, request, regime, recipient_sha, seed, *, load_values=True):
    manifest = load_schema(directory / 'prediction.json', PREDICTION_SCHEMA)
    bridge = load_schema(pinned(directory, 'bridge.json', manifest['bridge_sha256']), FORMAT)
    cases = rows(request['predictor_files']['test_cases'])
    groups, _ = grouped_cases(cases)
    expected = {'arm': 'Tahoe', 'regime': regime, 'seed': seed, 'case_count': len(cases),
                'group_count': len(groups), 'source_sha256': source_sha, 'recipient_sha256': recipient_sha,
                'recipient_identity': request['identity'], 'recipient_request_sha256': request_digest(request),
                'source_parent_overlap': overlap(source, request), 'dose_limitation': LIMITATION}
    if (any(manifest[key] != value for key, value in expected.items())
            or any(bridge[key] != expected[key] for key in
                   ('recipient_identity', 'recipient_request_sha256', 'dose_limitation'))
            or bridge['source'] != source['checkpoints'][str(seed)]
            or bridge['mapped_state_sha256'] != manifest['mapped_state_sha256']):
        raise ValueError('prediction evidence differs from frozen Tahoe and tenv identities')
    order = read_json(pinned(directory, 'cases.json', manifest['cases_sha256']))
    if order != [[row['case_id'], row['group_id']] for row in cases]:
        raise ValueError('prediction cases differ from frozen tenv case order')
    if not load_values:
        return manifest, None
    values = load_array(pinned(directory, 'prediction.npy', manifest['prediction_sha256']),
                        (len(cases), request['n_genes']))
    if any(not np.isfinite(values[start:start + 32]).all() for start in range(0, len(cases), 32)):
        raise ValueError('nonfinite Tahoe zero-shot prediction')
    return manifest, values


@contextmanager
def score_truth(runtime, request, fold, warm_source, scratch):
    if fold is not None:
        groups, truth, labels, _ = evaluator_inputs(fold)
        yield groups, truth, labels
        return
    cases = rows(request['predictor_files']['test_cases'])
    groups, _ = grouped_cases(cases)
    _, manifest = release(runtime, warm_source)
    locations = expression_sources(warm_source, manifest['files'], cases)
    grouped_truth(scratch, groups, cases, locations)
    shape = (len(groups), request['n_genes'])
    yield (rows(scratch / 'groups.csv'), load_array(scratch / 'truth.npy', shape),
           [load_array(scratch / f'{direction}.npy', shape, dtype='bool') for direction in ('positive', 'negative')])


def score(runtime, source, source_sha, root, regime, recipient_sha, runs_root, output, seeds):
    request, _, fold, warm_source = recipient(runtime, root, regime, recipient_sha)
    verify_pair(source, request)
    verified = {seed: verify_prediction(runs_root / f'seed-{seed}', source, source_sha,
                                        request, regime, recipient_sha, seed) for seed in seeds}
    write_scores(runtime, request, fold, warm_source, output, 'Tahoe',
                 {seed: values for seed, (_, values) in verified.items()}, {
                     'schema': 'step05-tahoe-evaluation-v1', 'regime': regime,
                     'dose_limitation': LIMITATION, 'source_sha256': source_sha,
                     'recipient_sha256': recipient_sha, 'scored_seeds': list(seeds),
                     'complete': tuple(seeds) == SEEDS, 'recipient_identity': request['identity'],
                     'source_parent_overlap': overlap(source, request),
                     'prediction_manifests': {
                         str(seed): sha256(runs_root / f'seed-{seed}/prediction.json') for seed in seeds},
                 })


def write_scores(runtime, request, fold, warm_source, output, arm, predictions, evidence):
    if output.exists() or output.is_symlink():
        raise ValueError('score requires a fresh output directory')
    indices, _ = grouped_cases(rows(request['predictor_files']['test_cases']))
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as temporary:
        scratch = Path(temporary)
        pending = scratch / 'scores'
        pending.mkdir()
        curves, run_rows = [], []
        with score_truth(runtime, request, fold, warm_source, scratch) as (groups, truth, labels):
            with (pending / 'per_group.csv').open('w', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=GROUP_FIELDS)
                writer.writeheader()
                for seed, values in predictions.items():
                    slot = {'dataset': 'tenv', 'fold': 0, 'model': arm, 'seed': seed}
                    run_rows.append(run_scores(writer, slot, 'completed', '', groups, truth, labels,
                                               indices, values, curves))
        write_csv(pending / 'per_run.csv', run_rows, RUN_FIELDS)
        write_csv(pending / 'summary.csv', summaries(run_rows))
        write_json(pending / 'tenv/pr_curves.json', {
            'schema': 'step04-pr-curves-v1', 'dataset': 'tenv', 'recall': RECALL_GRID.tolist(),
            'per_run_sha256': sha256(pending / 'per_run.csv'), 'runs': curves, 'dose_limitation': LIMITATION,
        })
        write_json(pending / 'evaluation.json', {
            **evidence,
            'files': {str(path.relative_to(pending)): sha256(path)
                      for path in sorted(pending.rglob('*')) if path.is_file()},
        })
        pending.rename(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime', default=REPO, type=Path)
    commands = parser.add_subparsers(dest='command', required=True)
    source_parser = commands.add_parser('verify-source')
    source_parser.add_argument('formal', type=Path)
    source_parser.add_argument('output', type=Path)
    for name in ('predict', 'score'):
        worker = commands.add_parser(name)
        worker.add_argument('source', type=Path)
        worker.add_argument('recipient', type=Path)
        worker.add_argument('regime', choices=REGIMES)
        if name == 'predict':
            worker.add_argument('seed', type=int, choices=SEEDS)
        else:
            worker.add_argument('runs', type=Path)
            worker.add_argument('--seed', type=int, choices=SEEDS)
        worker.add_argument('output', type=Path)
        worker.add_argument('--source-sha256', required=True)
        worker.add_argument('--recipient-sha256', required=True)
    args = parser.parse_args()
    runtime, output = args.runtime.resolve(), args.output.absolute()
    if args.command == 'verify-source':
        verify_source(runtime, args.formal.resolve(), output)
        return
    source = load_schema(verify_file(args.source, args.source_sha256), SOURCE_SCHEMA)
    inputs = (runtime, source, args.source_sha256, args.recipient.resolve(), args.regime, args.recipient_sha256)
    if args.command == 'predict':
        predict(*inputs, args.seed, output)
    else:
        score(*inputs, args.runs.resolve(), output, SEEDS if args.seed is None else (args.seed,))


if __name__ == '__main__':
    main()
