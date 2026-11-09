#!/usr/bin/env python3
"""Verify and replay the selected Step 05 baselines on fold 0 validation."""
import argparse
from collections import Counter
from pathlib import Path
import sys
import tempfile
import time

sys.dont_write_bytecode = True
REPO = next(path for path in Path(__file__).resolve().parents if (path / '.git').exists())
sys.path.insert(0, str(REPO / 'stage2/05_pdiff_transfer/execution/script'))

import transfer
from artifacts import SEEDS, read_json, rows, sha256, write_csv, write_json
from evaluate import grouped_cases, load_array, mean_delta, np, verify_file
from formal import completed_fit, slot_lock, use_runtime, verify_fit_seal
from formal_matrix import bind, check_binding
from metrics import directional_curve, score_group, truth_labels
from training_inputs import json_hash

SCHEMA = 'step06-validation-baseline-v1'
PREDICTION_SCHEMA = 'step06-sealed-prediction-v1'
CHECKPOINT_KEYS = ('identity', 'recipe_sha256', 'checkpoint_sha256', 'best_step',
                   'best_validation_loss', 'model_state_sha256')
INFERENCE_FILES = ('genes', 'controls', 'morgan', 'drugs', 'train_cases', 'validation_cases', 'test_cases')
ARMS = ('CM', 'Transfer')
REGIME = 'cold-drug'


def inference_request(request, role='validation'):
    if role not in ('validation', 'test') or request['target_kind'] != 'delta':
        raise ValueError('expected held-out delta prediction')
    files = request['predictor_files']
    return {'role': role, 'identity': request['identity'], 'n_genes': request['n_genes'],
            'predictor_files': {key: files[key] for key in INFERENCE_FILES},
            'file_sha256': {key: request['file_sha256'][key] for key in INFERENCE_FILES}}


def verify_inference(request):
    files, hashes = request['predictor_files'], request['file_sha256']
    if (files.keys() != set(INFERENCE_FILES) or hashes.keys() != files.keys()
            or hashes['genes'] != request['identity']['gene_sha256']
            or request['role'] not in ('validation', 'test')):
        raise ValueError('unexpected target-free predictor allowlist or gene axis')
    for key in INFERENCE_FILES:
        check_binding({'path': files[key], 'sha256': hashes[key]})
    return case_rows(request)


def case_rows(request):
    return rows(request['predictor_files'][f"{request['role']}_cases"])


def case_order(cases):
    return [[row['case_id'], row['group_id']] for row in cases]


def selection(directory, record, identity, config):
    """Bind a published native-loss selection to the checkpoint used in replay."""
    for name in ('checkpoint', 'history'):
        check_binding(record[name])
    if (Path(record['checkpoint']['path']) != directory / 'best.pt'
            or Path(record['history']['path']) != directory / 'run.json'):
        raise ValueError('published checkpoint and history point to a different fit')
    run = completed_fit(directory, identity, config)
    if (any(run[key] != record[key] for key in CHECKPOINT_KEYS if key != 'checkpoint_sha256')
            or run['checkpoint_sha256'] != record['checkpoint']['sha256']):
        raise ValueError('published EMA selection differs from the sealed fit')
    return {key: run[key] for key in CHECKPOINT_KEYS}


def verify_baseline(runtime, seed):
    step = runtime / transfer.STEP05
    record = transfer.RECORD / 'formal_record'
    seal = read_json(record / 'run.json')
    publication = read_json(verify_file(record / 'publication.json',
                                        seal['publication']['manifest']['sha256']))
    index_path = verify_file(record / 'checkpoints.json', publication['files']['checkpoints.json'])
    index = read_json(index_path)
    matrix_binding = seal['formal']['plan']
    matrix_path = verify_file(step / matrix_binding['path'], matrix_binding['sha256'])
    matrix = read_json(matrix_path)
    worker = matrix['worker_plans']['transfer']
    worker_path = Path(check_binding(worker)['path'])
    root, digest = worker_path.parent, worker['sha256']
    plan, context, _, _ = transfer.verify(runtime, root, digest, REGIME)
    request, _, _ = context[REGIME]
    target = plan['targets'][REGIME]
    inference = inference_request(request)
    cases = verify_inference(inference)
    groups, _ = grouped_cases(cases)
    doses = Counter(float(cases[members[0]]['dose_um']) for members in groups.values())
    if request['n_genes'] != 18164 or doses != {0.2: 46, 1.0: 46, 10.0: 46}:
        raise ValueError('unexpected fold 0 validation axis or dose coverage')

    inputs = transfer.cold_inputs(runtime)
    source, source_sha = transfer.cold_source(runtime, root, plan, digest, seed, inputs)
    source_slot = index['slots'][f'{REGIME}-tahoe-{seed}']
    source_record = index['checkpoints'][source_slot['source']]
    source_checkpoint = source['checkpoints'][str(seed)]
    source_dir = Path(source_checkpoint['path']).parent
    check_binding(source_slot['source_fit_seal'])
    if Path(source_slot['source_fit_seal']['path']) != source_dir / 'fit.json':
        raise ValueError('published source seal points to a different fit')
    verify_fit_seal(source_dir, digest)
    source_selected = selection(source_dir, source_record, source_checkpoint['identity'],
                                read_json(source_dir / 'run.json')['recipe'])
    if source_selected['checkpoint_sha256'] != source_checkpoint['checkpoint_sha256']:
        raise ValueError('published source differs from the globally cold Tahoe EMA')
    bridge_dir, bridge_manifest, _ = transfer.verify_tahoe(
        runtime, root, REGIME, seed, request, target, source, source_sha, stage='formal', load_values=False)
    bridge = bind(bridge_dir / 'bridge.json', bridge_manifest['bridge_sha256'])

    item = plan['active_recipes'][str(seed)]
    lineage = transfer.transfer_lineage(source, source_sha, target, item, bridge_manifest, seed)
    config = item['config']
    checkpoints = {}
    for arm in ARMS:
        slot = index['slots'][f'{REGIME}-{arm.lower()}-{seed}']
        selected = index['checkpoints'][slot['target']]
        identity = request['identity']
        directory = Path(selected['checkpoint']['path']).parent
        if arm == 'Transfer':
            identity = transfer.transfer_identity(request, REGIME)
            target_slot = root / 'slots' / REGIME / f'seed-{seed}'
            check_binding(slot['fit_seal'])
            if (directory != target_slot / 'fit' or Path(slot['fit_seal']['path']) != directory / 'fit.json'
                    or slot['initialization'] != lineage):
                raise ValueError('published Transfer seal or initialization differs from the original bridge')
            transfer.verify_transfer_fit(target_slot, digest, {**identity, 'seed': seed, 'stage': 'formal'},
                                         item, lineage)
        else:
            cm = matrix['reuse_cm'][str(seed)]
            if (selected['checkpoint']['sha256'] != cm['checkpoint_sha256']
                    or directory != Path(cm['checkpoint_directory'])
                    or read_json(check_binding(selected['request'])['path']) != request):
                raise ValueError('CM checkpoint differs from the matched frozen recipient')
            check_binding({'path': str(directory / 'replay.json'), 'sha256': cm['replay_sha256']})
        metadata = selection(directory, selected, {**identity, 'seed': seed, 'stage': 'formal'}, config)
        batching = read_json(directory / 'run.json')['validation_batching']
        if batching['group_ids'] != list(groups) or batching['group_ids_sha256'] != json_hash(list(groups)):
            raise ValueError(f'{arm} native selection used a different validation group order')
        if seed == 17 and metadata['best_step'] != 6500:
            raise ValueError('seed 17 baseline must reuse the selected step 6500 checkpoint')
        checkpoints[arm] = {'checkpoint': selected['checkpoint'], 'history': selected['history'],
                            'metadata': metadata}
    return {'schema': SCHEMA, 'seed': seed, 'request': inference, 'config': config,
            'inference': target['inference'], 'ordered_case_group_sha256': json_hash(case_order(cases)),
            'case_count': len(cases), 'group_count': len(groups),
            'dose_group_counts': {str(dose): doses[dose] for dose in sorted(doses)},
            'checkpoints': checkpoints, 'source': {'checkpoint': source_record['checkpoint'],
                'history': source_record['history'], 'fit_seal': source_slot['source_fit_seal'],
                'source_sha256': source_sha, 'roles': source['roles'], 'bridge': bridge},
            'release': {'step05_matrix': bind(matrix_path, matrix_binding['sha256']),
                        'step05_checkpoints': bind(index_path, publication['files']['checkpoints.json']),
                        'step05_transfer': worker, 'recipient_split_sha256': target['split_sha256']}}


def predict_checkpoint(request, config, inference, checkpoint, metadata, destination, lineage):
    """Run the Step 05 bounded sampler on a held-out role with no target input."""
    from native_training import (FORMAT, JointPerturbDiff, PredictorInputs, digest, resources,
                                 seed_randomness, torch)

    cases = verify_inference(request)
    checkpoint = check_binding(checkpoint)
    if (checkpoint['sha256'] != metadata['checkpoint_sha256']
            or metadata['recipe_sha256'] != json_hash(config)
            or metadata['identity']['seed'] != config['seed']
            or any(metadata['identity'][key] != request['identity'][key] for key in request['identity'])):
        raise ValueError('checkpoint binding, recipe or case/gene identity differs from selected metadata')

    class HeldOutInputs(PredictorInputs):
        def __init__(self):
            files = request['predictor_files']
            self.n_genes = request['n_genes']
            # The model builds its vocabulary from all roles. The Step 05
            # sampler uses its prediction role under the name test.
            self.cases = {f'vocabulary-{role}': rows(files[f'{role}_cases'])
                          for role in transfer.cold.ROLES if role != request['role']}
            self.cases['test'] = cases
            self.controls = np.load(files['controls'], mmap_mode='r', allow_pickle=False)
            self.morgan = np.load(files['morgan'], mmap_mode='r', allow_pickle=False)
            self.drug_index = {row['drug_key']: int(row['morgan_index']) for row in rows(files['drugs'])}
            self.delta = {}

    destination.parent.mkdir(parents=True, exist_ok=True)
    with slot_lock(destination.parent / f'.{destination.name}-lock'):
        if destination.exists():
            verify_prediction(request, checkpoint, metadata, destination, lineage)
            return
        staging = Path(tempfile.mkdtemp(prefix=f'{destination.name}-attempt-', dir=destination.parent))
        print(f'Prediction attempt evidence: {staging}', flush=True)
        started = time.monotonic()
        torch.set_num_threads(config['threads'])
        seed_randomness(config['seed'])
        device = torch.device(config['device'])
        if device.type == 'cuda':
            torch.cuda.reset_peak_memory_stats(device)
        inputs = HeldOutInputs()
        model = JointPerturbDiff(inputs, config, staging)
        state = torch.load(checkpoint['path'], map_location='cpu', weights_only=True)
        if state['format'] != FORMAT or any(state[key] != metadata[key]
                                             for key in CHECKPOINT_KEYS if key != 'checkpoint_sha256'):
            raise ValueError('selected EMA checkpoint metadata mismatch')
        model.load_state_dict(state['model_state'], strict=True)
        del state
        if digest(model) != metadata['model_state_sha256']:
            raise ValueError('selected EMA state digest mismatch')
        model.to(device)
        model.native.to(device)
        model.eval()
        counts = transfer.tahoe.sample_cases(model, inputs, config, inference, staging, device)
        write_json(staging / 'prediction.json', {
            'schema': PREDICTION_SCHEMA, 'lineage': lineage, 'checkpoint': checkpoint,
            'metadata': metadata, 'request_sha256': json_hash(request), **counts,
            'ordered_case_group_sha256': json_hash(case_order(cases)),
            'cases_sha256': sha256(staging / 'cases.json'),
            'prediction_sha256': sha256(staging / 'prediction.npy'),
            'wall_seconds': time.monotonic() - started,
            **({'resources': resources(device)} if device.type == 'cuda' else {}),
        })
        staging.rename(destination)


def verify_prediction(request, checkpoint, metadata, destination, lineage):
    cases = verify_inference(request)
    groups, _ = grouped_cases(cases)
    manifest = read_json(destination / 'prediction.json')
    expected = {'schema': PREDICTION_SCHEMA, 'lineage': lineage, 'checkpoint': checkpoint,
                'metadata': metadata, 'request_sha256': json_hash(request), 'case_count': len(cases),
                'group_count': len(groups), 'ordered_case_group_sha256': json_hash(case_order(cases))}
    order = read_json(verify_file(destination / 'cases.json', manifest['cases_sha256']))
    if any(manifest[key] != value for key, value in expected.items()) or order != case_order(cases):
        raise ValueError('prediction lineage or ordered validation case membership differs')
    return load_array(verify_file(destination / 'prediction.npy', manifest['prediction_sha256']),
                      (len(cases), request['n_genes']))


def score_prediction(request, observed, predicted, *, include_pcc=False):
    cases = case_rows(request)
    groups, _ = grouped_cases(cases)
    shape = (len(cases), request['n_genes'])
    if observed.shape != shape or predicted.shape != shape:
        raise ValueError('prediction or truth differs from the frozen case and gene axes')
    group_rows = []
    for group, members in groups.items():
        case = cases[members[0]]
        if any(cases[index][key] != case[key] for index in members for key in ('drug_key', 'dose_um')):
            raise ValueError('a validation group spans multiple drug–dose conditions')
        truth = mean_delta(observed, members).astype('float32')
        delta = mean_delta(predicted, members)
        positive, negative = truth_labels(truth)
        row = {'group_id': group, 'drug_key': case['drug_key'], 'dose_um': float(case['dose_um']),
               'n_wells': len(members)}
        if include_pcc:
            values, reasons, _ = score_group(truth, delta, (positive, negative))
            for field, metric in (('positive_ap100', 'AUPRC_positive@100'),
                                  ('negative_ap100', 'AUPRC_negative@100'), ('pcc', 'PCC')):
                reason = field.replace('_ap100', '') + '_reason'
                row[field] = values[metric] if metric not in reasons else None
                row[reason] = reasons.get(metric, '')
        else:
            for name, labels, values in (('positive', positive, delta), ('negative', negative, -delta)):
                ap, _, _ = directional_curve(labels, values)
                row[f'{name}_ap100'] = float(ap) if np.isfinite(ap) else None
                row[f'{name}_reason'] = '' if np.isfinite(ap) else 'empty_direction'
        group_rows.append(row)
    return group_rows, summarize_groups(group_rows)


def summarize_groups(group_rows, *, metric_named_support=False):
    summaries = []
    strata = [('pooled', group_rows)]
    for dose in sorted({row['dose_um'] for row in group_rows}):
        strata.append((str(dose), [row for row in group_rows if row['dose_um'] == dose]))
    for dose, members in strata:
        summary = {'dose_um': dose, 'group_count': len(members)}
        for direction in ('positive', 'negative'):
            key = f'{direction}_ap100'
            valid = [row[key] for row in members if row[key] is not None]
            summary[key] = float(np.mean(valid)) if valid else None
            support = key if metric_named_support else direction
            summary[f'{support}_valid_groups'] = len(valid)
            summary[f'{support}_undefined_groups'] = len(members) - len(valid)
        values = [summary[key] for key in ('positive_ap100', 'negative_ap100')]
        summary['mean_ap100'] = None if None in values else float(np.mean(values))
        summaries.append(summary)
    return summaries


def score_baseline(runtime, root, plan, digest):
    request = plan['request']
    values, evidence = {}, {}
    for arm in ARMS:
        selected = plan['checkpoints'][arm]
        check_binding(selected['checkpoint'])
        check_binding(selected['history'])
        destination = root / arm
        values[arm] = verify_prediction(request, selected['checkpoint'], selected['metadata'],
                                        destination, {'plan_sha256': digest, 'arm': arm})
        evidence[arm] = bind(destination / 'prediction.json')

    # Join validation targets from the pinned Step 04 request only after both
    # predictions have been verified. The prediction worker never opens them.
    historical = transfer.tahoe.historical_record('tenv')
    entry = read_json(historical / 'plan.json')['folds'][0]
    predecessor = runtime / transfer.STEP05 / 'execution/output/cold-tenv-ticket76'
    full = read_json(verify_file(predecessor / 'predictor_request.json', entry['request_sha256']))
    if inference_request(full) != request:
        raise ValueError('validation truth request does not match the prediction request')
    truth = bind(full['predictor_files']['validation_expression'], full['file_sha256']['validation_expression'])
    observed = load_array(truth['path'], (plan['case_count'], request['n_genes']))
    by_arm, group_rows, summary_rows = {}, [], []
    for arm in ARMS:
        group, summary = score_prediction(request, observed, values[arm])
        labels = {'seed': plan['seed'], 'arm': arm, 'step': plan['checkpoints'][arm]['metadata']['best_step']}
        group_rows.extend({**labels, **row} for row in group)
        summary_rows.extend({**labels, **row} for row in summary)
        by_arm[arm] = summary
    paired = []
    for cm, original in zip(by_arm['CM'], by_arm['Transfer'], strict=True):
        gaps = {}
        for metric in ('positive_ap100', 'negative_ap100', 'mean_ap100'):
            gaps[metric] = None if None in (cm[metric], original[metric]) else cm[metric] - original[metric]
        paired.append({'dose_um': cm['dose_um'], 'cm_minus_transfer': gaps})
    destination = root / 'scores'
    with slot_lock(root / '.score-lock'), tempfile.TemporaryDirectory(dir=root) as temp:
        staging = Path(temp) / 'ready'
        write_csv(staging / 'per_group.csv', group_rows)
        write_csv(staging / 'per_dose.csv', summary_rows)
        report = {'schema': 'step06-validation-scores-v1', 'seed': plan['seed'], 'fold': 0,
                  'plan_sha256': digest, 'truth': truth, 'predictions': evidence,
                  'checkpoints': plan['checkpoints'], 'per_dose': summary_rows, 'paired': paired,
                  'files': {path.name: sha256(path) for path in sorted(staging.iterdir())}}
        write_json(staging / 'scores.json', report)
        if destination.exists():
            for path in staging.iterdir():
                verify_file(destination / path.name, sha256(path))
        else:
            staging.rename(destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('verify', 'predict', 'score'))
    parser.add_argument('root', type=Path)
    parser.add_argument('--runtime', type=Path, default=REPO)
    parser.add_argument('--seed', type=int, choices=SEEDS, default=17)
    parser.add_argument('--plan-sha256')
    args = parser.parse_args()
    root, runtime = args.root.resolve(), args.runtime.resolve()
    path = root / 'plan.json'
    if args.phase == 'verify':
        plan = verify_baseline(runtime, args.seed)
        root.mkdir(parents=True, exist_ok=True)
        with slot_lock(root / '.plan-lock'), tempfile.TemporaryDirectory(dir=root) as temp:
            pending = Path(temp) / 'plan.json'
            write_json(pending, plan)
            if path.exists():
                verify_file(path, sha256(pending))
            else:
                pending.rename(path)
        print(f'Validation baseline plan SHA-256 {sha256(path)}', flush=True)
        return
    plan = read_json(verify_file(path, args.plan_sha256))
    if plan['schema'] != SCHEMA or plan['seed'] != args.seed:
        raise ValueError('unexpected baseline plan schema or seed')
    if args.phase == 'predict':
        use_runtime(runtime)
        for arm in ARMS:
            selected = plan['checkpoints'][arm]
            predict_checkpoint(plan['request'], plan['config'], plan['inference'], selected['checkpoint'],
                               selected['metadata'], root / arm, {'plan_sha256': args.plan_sha256, 'arm': arm})
    else:
        score_baseline(runtime, root, plan, args.plan_sha256)


if __name__ == '__main__':
    main()
