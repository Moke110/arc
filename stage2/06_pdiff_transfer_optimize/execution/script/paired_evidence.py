#!/usr/bin/env python3
"""Verify sealed Step 06 paired Transfer evidence for all three seeds."""
import argparse
import math
from pathlib import Path

import paired_reference as reference
import paired_target as target
from artifacts import SEEDS, read_json, rows
from formal import verify_fit_seal
from formal_matrix import bind, check_binding
from learning_rate import publish
from training_inputs import json_hash

SCHEMA = 'step06-paired-evidence-v1'
GROUP_METRICS = {
    'positive_ap100': ('positive_reason', {'empty_direction'}),
    'negative_ap100': ('negative_reason', {'empty_direction'}),
    'pcc': ('pcc_reason', {'constant_truth', 'constant_prediction'}),
}


def verify_group_scores(path, plan, labels):
    cases = reference.base.case_rows(plan['request'])
    groups, _ = reference.grouped_cases(cases)
    records = rows(path)
    if (len(cases) != plan['case_count'] or len(groups) != plan['group_count']
            or json_hash(reference.base.case_order(cases)) != plan['ordered_case_group_sha256']
            or [row['group_id'] for row in records] != list(groups)):
        raise ValueError('group score membership differs from the frozen protocol')
    for row in records:
        members = groups[row['group_id']]
        case = cases[members[0]]
        if (any(row[key] != str(value) for key, value in labels.items())
                or row['drug_key'] != case['drug_key'] or float(row['dose_um']) != float(case['dose_um'])
                or int(row['n_wells']) != len(members)):
            raise ValueError('group score identity differs from the frozen protocol or EMA')
        row['dose_um'] = float(row['dose_um'])
        for field, (reason, allowed) in GROUP_METRICS.items():
            value = float(row[field]) if row[field] else None
            if ((value is None and row[reason] not in allowed)
                    or (value is not None and (row[reason] or not math.isfinite(value)))):
                raise ValueError('group score has invalid metric support')
            row[field] = value
    return reference.summarize_doses(records, reference.base.summarize_groups(records))


def verify_scores(directory, prediction_dir, plan, plan_digest, arm, selection, scores_binding):
    if bind(directory / 'scores.json', scores_binding['sha256']) != scores_binding:
        raise ValueError('score seal points to a different run')
    report = read_json(scores_binding['path'])
    lineage = {'plan_sha256': plan_digest, 'protocol_sha256': plan['protocol_sha256'], 'arm': arm}
    run = read_json(selection['history']['path'])
    expected = {
        **lineage, 'schema': reference.SCORE_SCHEMA, 'seed': plan['seed'], 'fold': plan['fold'],
        'regime': plan['regime'], 'selection': selection, 'truth': plan['truth'],
        'optimizer_steps': run['optimizer_steps'],
    }
    if any(report[key] != value for key, value in expected.items()):
        raise ValueError('score lineage differs from the frozen protocol or selected fit')
    prediction = bind(prediction_dir / 'prediction.json', report['prediction']['sha256'])
    if prediction != report['prediction']:
        raise ValueError('score seal points to a different prediction')
    reference.base.verify_prediction(plan['request'], selection['checkpoint'], selection['metadata'],
                                     prediction_dir, lineage)
    check_binding(plan['truth'])
    files = {name: bind(directory / name, report['files'][name]) for name in ('per_group.csv', 'per_dose.csv')}
    labels = {'seed': plan['seed'], 'arm': arm, 'step': selection['metadata']['best_step']}
    summaries = [{**labels, **row} for row in verify_group_scores(files['per_group.csv']['path'], plan, labels)]
    csv_summaries = [{key: '' if value is None else str(value) for key, value in row.items()}
                     for row in summaries]
    if (summaries != report['per_dose'] or csv_summaries != rows(files['per_dose.csv']['path'])
            or {row['dose_um']: row['group_count'] for row in summaries} != plan['dose_group_counts']
            or tuple(float(row['dose_um']) for row in summaries) != reference.DOSES):
        raise ValueError('dose summaries differ from the sealed group scores or support')
    return report, files


def verify_cell(root, arm, seed):
    """Read one intervention and the same-seed reference with their artifact bindings."""
    upstream = root / 'reference' / f'seed-{seed}'
    cell = root / arm / f'seed-{seed}'
    plan_binding = bind(cell / 'plan.json')
    pair_binding = bind(cell / 'paired.json')
    pair = read_json(pair_binding['path'])
    plan_digest = plan_binding['sha256']
    plan, original = target.load_plan(cell, arm, seed, plan_digest)
    if (plan['reference'] != bind(upstream / 'plan.json', plan['reference']['sha256'])
            or any(pair[key] != value for key, value in {
                'schema': target.PAIR_SCHEMA, 'arm': arm, 'seed': seed, 'plan_sha256': plan_digest,
                'protocol_sha256': plan['protocol_sha256'], 'reference': plan['reference_scores'],
            }.items())):
        raise ValueError('paired score lineage differs from the same-seed reference or intervention')

    baseline = original['checkpoints'][reference.ARM]
    directory = Path(check_binding(baseline['fit_seal'])['path']).parent
    verify_fit_seal(directory, read_json(directory / 'fit.json')['plan_sha256'])
    identity = target.fit_identity({'arm': reference.ARM, 'regime': plan['regime'],
                                    'seed': seed, 'config': original['config']}, original)
    if baseline['metadata'] != reference.base.selection(directory, {**baseline, **baseline['metadata']},
                                                         identity, original['config']):
        raise ValueError('original Transfer selection differs from the same-seed reference')
    original_report, original_files = verify_scores(
        upstream / 'scores', upstream / reference.ARM, original, plan['reference']['sha256'],
        reference.ARM, baseline, pair['reference'])

    selection = target.selected_fit(plan, original, plan_digest)
    score_report, files = verify_scores(cell / 'scores', cell / 'prediction', original, plan_digest,
                                        arm, selection, pair['scores'])
    if pair['per_dose'] != target.paired_differences(score_report, original_report):
        raise ValueError('paired dose differences differ from the same-seed scores')
    return {
        'seed': seed, 'arm': arm, 'protocol_sha256': plan['protocol_sha256'],
        'plan': plan_binding, 'paired': pair_binding, 'source': plan['source'], 'target': selection,
        'optimizer_steps': score_report['optimizer_steps'], 'scores': pair['scores'],
        'prediction': score_report['prediction'], 'score_files': files, 'per_dose': pair['per_dose'],
        'reference': {
            'plan': plan['reference'], 'target': baseline, 'scores': pair['reference'],
            'optimizer_steps': original_report['optimizer_steps'], 'prediction': original_report['prediction'],
            'score_files': original_files, 'per_dose': original_report['per_dose'],
        },
    }


def verify_matrix(root, arm):
    cells = [verify_cell(root, arm, seed) for seed in SEEDS]
    protocol = cells[0]['protocol_sha256']
    if any(cell['protocol_sha256'] != protocol for cell in cells):
        raise ValueError('paired seeds do not share the frozen validation protocol')
    return {'schema': SCHEMA, 'arm': arm, 'regime': reference.REGIME, 'fold': 0,
            'protocol_sha256': protocol, 'cells': cells}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path, help='paired-transfer output directory')
    parser.add_argument('--arm', choices=target.CHECKPOINT_NUMBERS, default='log')
    parser.add_argument('--manifest', type=Path, help='defaults to <root>/<arm>/matrix.json')
    args = parser.parse_args(argv)
    root = args.root.resolve()
    destination = args.manifest or root / args.arm / 'matrix.json'
    binding = publish(destination.resolve(), verify_matrix(root, args.arm))
    print(f'{binding["path"]} SHA-256 {binding["sha256"]}', flush=True)


if __name__ == '__main__':
    main()
