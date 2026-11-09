#!/usr/bin/env python3
"""Replay the original Step 05 cold-cm-drug Transfer at its selected EMA."""
import argparse
from collections import Counter
from pathlib import Path
import tempfile

import baseline as base
from artifacts import SEEDS, read_json, sha256, write_csv, write_json
from evaluate import grouped_cases, load_array, np, verify_file
from formal import slot_lock, use_runtime, verify_fit_seal
from formal_matrix import bind, check_binding
from learning_rate import publish
from training_inputs import json_hash

SCHEMA = 'step06-paired-reference-v1'
SCORE_SCHEMA = 'step06-paired-scores-v1'
REGIME = 'cold-cm-drug'
ARM = 'Transfer'
DOSES = (0.2, 1.0, 10.0)
TRANSFER = base.transfer


def score_prediction(request, observed, predicted):
    groups, summaries = base.score_prediction(request, observed, predicted, include_pcc=True)
    return groups, summarize_doses(groups, summaries)


def summarize_doses(groups, summaries):
    doses = []
    for summary in summaries[1:]:
        members = [row for row in groups if row['dose_um'] == float(summary['dose_um'])]
        valid = [row['pcc'] for row in members if row['pcc'] is not None]
        doses.append({
            **summary,
            'pcc': float(np.mean(valid)) if valid else None,
            'pcc_valid_groups': len(valid),
            'pcc_undefined_groups': len(members) - len(valid),
            **{f'{metric}_undefined_reasons': dict(Counter(row[f'{metric}_reason'] for row in members
                                                          if row[f'{metric}_reason']))
               for metric in ('positive', 'negative', 'pcc')},
        })
    return doses


def verify_baseline(runtime, seed):
    step = runtime / TRANSFER.STEP05
    record = TRANSFER.RECORD / 'formal_record'
    seal = read_json(record / 'run.json')
    publication = read_json(verify_file(record / 'publication.json',
                                        seal['publication']['manifest']['sha256']))
    index = bind(record / 'checkpoints.json', publication['files']['checkpoints.json'])
    checkpoints = read_json(index['path'])
    matrix_binding = seal['formal']['plan']
    matrix = bind(step / matrix_binding['path'], matrix_binding['sha256'])
    worker = check_binding(read_json(matrix['path'])['worker_plans']['transfer'])
    root, digest = Path(worker['path']).parent, worker['sha256']
    plan, context, source, source_sha = TRANSFER.verify(runtime, root, digest, REGIME)
    request, _, _ = context[REGIME]
    target = plan['targets'][REGIME]
    inference = base.inference_request(request)
    cases = base.verify_inference(inference)
    groups, _ = grouped_cases(cases)
    doses = Counter(float(cases[members[0]]['dose_um']) for members in groups.values())
    if doses != dict.fromkeys(DOSES, 46):
        raise ValueError('unexpected cold-cm-drug fold 0 validation dose coverage')
    slot = checkpoints['slots'][f'{REGIME}-transfer-{seed}']
    if (slot['regime'], slot['arm'], slot['fold'], slot['seed'], slot['status']) != (
            REGIME, ARM, 0, seed, 'completed'):
        raise ValueError('expected the original completed Transfer fit for this seed')
    source_record = checkpoints['checkpoints'][slot['source']]
    source_checkpoint = source['checkpoints'][str(seed)]
    source_dir = Path(source_checkpoint['path']).parent
    source_manifest = bind(step / 'execution/output/tahoe-zero-shot-ticket78/source.json', source_sha)
    source_request = check_binding(source_record['request'])
    raw_source = read_json(source_request['path'])
    if (source_request['sha256'] != source['request_file_sha256']
            or raw_source['identity'] != source['request_identity']):
        raise ValueError('published Tahoe request differs from the frozen source')
    source_cases = TRANSFER.tahoe.role_cases(raw_source)
    source_parents = TRANSFER.tahoe.parent_blocks(raw_source, source_cases)
    for role, members in source_cases.items():
        case_key = f'{role}_cases'
        digest_cases = raw_source['file_sha256'][case_key]
        bind(raw_source['predictor_files'][case_key], digest_cases)
        if source['roles'][role] != {
                'parent_connectivity': source_parents[role], 'case_count': len(members),
                'group_count': len(grouped_cases(members)[0]), 'cases_sha256': digest_cases}:
            raise ValueError('Tahoe fitting membership differs from the Step 05 release')
    bind(raw_source['predictor_files']['drugs'], raw_source['file_sha256']['drugs'])
    verify_fit_seal(source_dir, source['plan_sha256'])
    source_run = read_json(source_dir / 'run.json')
    source_batching = source_run['validation_batching']
    source_groups = list(grouped_cases(source_cases['validation'])[0])
    source_selected = base.selection(source_dir, source_record, source_checkpoint['identity'], source_run['recipe'])
    if (source_selected != {key: source_checkpoint[key] for key in base.CHECKPOINT_KEYS}
            or source_record['history']['sha256'] != source_checkpoint['producer']['run_sha256']
            or source_batching['group_ids'] != source_groups
            or source_batching['group_ids_sha256'] != json_hash(source_groups)):
        raise ValueError('published Tahoe EMA differs from the frozen source selection')

    bridge_dir, bridge_manifest, _ = TRANSFER.verify_tahoe(
        runtime, root, REGIME, seed, request, target, source, source_sha, stage='formal', load_values=False)
    bridge = bind(bridge_dir / 'bridge.json', bridge_manifest['bridge_sha256'])
    item = plan['active_recipes'][str(seed)]
    lineage = TRANSFER.transfer_lineage(source, source_sha, target, item, bridge_manifest, seed)
    selected = checkpoints['checkpoints'][slot['target']]
    recipient = check_binding(selected['request'])
    if read_json(recipient['path']) != request:
        raise ValueError('published Transfer request differs from the frozen recipient')
    identity = {**TRANSFER.transfer_identity(request, REGIME), 'seed': seed, 'stage': 'formal'}
    target_slot = root / 'slots' / REGIME / f'seed-{seed}'
    directory = target_slot / 'fit'
    fit_seal = check_binding(slot['fit_seal'])
    if (Path(fit_seal['path']) != directory / 'fit.json' or slot['initialization'] != lineage):
        raise ValueError('published Transfer seal or initialization differs from the original bridge')
    run = TRANSFER.verify_transfer_fit(target_slot, digest, identity, item, lineage)
    metadata = base.selection(directory, selected, identity, item['config'])
    batching = run['validation_batching']
    if batching['group_ids'] != list(groups) or batching['group_ids_sha256'] != json_hash(list(groups)):
        raise ValueError('Transfer native selection used a different validation group order')
    truth = bind(request['predictor_files']['validation_expression'],
                 request['file_sha256']['validation_expression'])
    frozen = {
        'fold': 0, 'regime': REGIME, 'request': inference, 'inference': target['inference'],
        'ordered_case_group_sha256': json_hash(base.case_order(cases)),
        'case_count': len(cases), 'group_count': len(groups),
        'dose_group_counts': {str(dose): doses[dose] for dose in DOSES}, 'truth': truth,
        'source_exposure': {
            'policy': 'tenv parent-connectivity holdout with original Tahoe membership',
            'source_manifest': source_manifest, 'source_request': source_request, 'roles': source['roles'],
            'parent_overlap': bridge_manifest['source_parent_overlap'],
        },
    }
    return {
        'schema': SCHEMA, 'seed': seed, **frozen, 'protocol_sha256': json_hash(frozen),
        'config': item['config'], 'recipient_request': recipient,
        'target_recipe': check_binding(selected['recipe']),
        'checkpoints': {ARM: {'checkpoint': selected['checkpoint'], 'history': selected['history'],
                              'metadata': metadata, 'fit_seal': fit_seal}},
        'source': {'checkpoint': source_record['checkpoint'], 'history': source_record['history'],
                   'metadata': source_selected, 'fit_seal': bind(source_dir / 'fit.json'),
                   'request': source_request, 'recipe': check_binding(source_record['recipe']),
                   'source_sha256': source_sha, 'roles': source['roles'], 'bridge': bridge,
                   'parent_overlap': bridge_manifest['source_parent_overlap']},
        'release': {'step05_matrix': matrix, 'step05_checkpoints': index, 'step05_transfer': worker,
                    'recipient_split_sha256': target['split_sha256']},
    }


def score_paired(root, plan, digest, *, arm, selected, prediction_dir, lineage):
    """Publish the common validation score artifact for one selected EMA."""
    request = plan['request']
    cases = base.verify_inference(request)
    groups, _ = grouped_cases(cases)
    if (plan['case_count'] != len(cases) or plan['group_count'] != len(groups)
            or plan['ordered_case_group_sha256'] != json_hash(base.case_order(cases))):
        raise ValueError('score case membership differs from the frozen reference')
    checkpoint = check_binding(selected['checkpoint'])
    history = check_binding(selected['history'])
    run = read_json(history['path'])
    record = {**selected['metadata'], 'checkpoint': checkpoint, 'history': history}
    metadata = base.selection(Path(checkpoint['path']).parent, record, record['identity'], run['recipe'])
    if (metadata != selected['metadata'] or metadata['identity']['seed'] != plan['seed']
            or any(metadata['identity'][key] != value for key, value in request['identity'].items())):
        raise ValueError('selected EMA lineage differs from the paired validation protocol')
    predicted = base.verify_prediction(request, checkpoint, metadata, prediction_dir, lineage)
    truth = check_binding(plan['truth'])
    observed = load_array(Path(truth['path']), (len(cases), request['n_genes']))
    group_rows, dose_rows = score_prediction(request, observed, predicted)
    if {row['dose_um']: row['group_count'] for row in dose_rows} != plan['dose_group_counts']:
        raise ValueError('score dose membership differs from the frozen reference')
    labels = {'seed': plan['seed'], 'arm': arm, 'step': metadata['best_step']}
    summaries = [{**labels, **row} for row in dose_rows]
    destination = root / 'scores'
    with slot_lock(root / '.score-lock'), tempfile.TemporaryDirectory(dir=root) as temp:
        staging = Path(temp) / 'ready'
        write_csv(staging / 'per_group.csv', [{**labels, **row} for row in group_rows])
        write_csv(staging / 'per_dose.csv', summaries)
        report = {
            'schema': SCORE_SCHEMA, 'seed': plan['seed'], 'fold': plan['fold'], 'regime': plan['regime'],
            'arm': arm, 'plan_sha256': digest, 'protocol_sha256': plan['protocol_sha256'],
            'truth': truth, 'prediction': bind(prediction_dir / 'prediction.json'),
            'selection': selected, 'optimizer_steps': run['optimizer_steps'], 'per_dose': summaries,
            'files': {path.name: sha256(path) for path in sorted(staging.iterdir())},
        }
        write_json(staging / 'scores.json', report)
        if destination.exists():
            for path in staging.iterdir():
                verify_file(destination / path.name, sha256(path))
        else:
            staging.rename(destination)
    return bind(destination / 'scores.json')


def load_plan(root, seed, digest):
    plan = read_json(verify_file(root / 'plan.json', digest))
    if ((plan['schema'], plan['seed'], plan['fold'], plan['regime'], plan['request']['role'])
            != (SCHEMA, seed, 0, REGIME, 'validation')):
        raise ValueError('expected cold-cm-drug fold 0 reference plan for this seed')
    return plan


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('verify', 'predict', 'score'))
    parser.add_argument('root', type=Path)
    parser.add_argument('--runtime', type=Path, default=base.REPO)
    parser.add_argument('--seed', type=int, choices=SEEDS, required=True)
    parser.add_argument('--plan-sha256')
    args = parser.parse_args(argv)
    root, runtime = args.root.resolve(), args.runtime.resolve()
    if args.phase == 'verify':
        plan = verify_baseline(runtime, args.seed)
        result = publish(root / 'plan.json', plan)
    else:
        plan = load_plan(root, args.seed, args.plan_sha256)
        selected = plan['checkpoints'][ARM]
        prediction_dir = root / ARM
        lineage = {'plan_sha256': args.plan_sha256, 'protocol_sha256': plan['protocol_sha256'], 'arm': ARM}
        if args.phase == 'predict':
            use_runtime(runtime)
            base.predict_checkpoint(plan['request'], plan['config'], plan['inference'], selected['checkpoint'],
                                    selected['metadata'], prediction_dir, lineage)
            result = bind(prediction_dir / 'prediction.json')
        else:
            result = score_paired(root, plan, args.plan_sha256, arm=ARM, selected=selected,
                                  prediction_dir=prediction_dir, lineage=lineage)
    print(f'{result["path"]} SHA-256 {result["sha256"]}', flush=True)


if __name__ == '__main__':
    main()
