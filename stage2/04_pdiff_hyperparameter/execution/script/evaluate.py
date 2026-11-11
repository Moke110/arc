#!/usr/bin/env python3
"""Verify and score frozen Step 04 held-out prediction artifacts."""
import argparse
import csv
import math
import os
from pathlib import Path
import re
import sys
import tempfile
from collections import defaultdict

for name in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(name, '2')
sys.dont_write_bytecode = True
REPO = next(path for path in Path(__file__).resolve().parents if (path / '.git').exists())
sys.path.insert(0, str(REPO / 'stage2/03_modes_tenv/execution/script'))

import numpy as np

from artifacts import SEEDS, read_json, rows, sha256, write_csv, write_json
from metrics import METRICS, RECALL_GRID, score_group, truth_labels
from training_inputs import INPUT_FIELDS, STEP03, json_hash, request_for, tenv_paths
from training_policy import MONITOR, verify_best_checkpoint

DATASETS = ('tahoe', 'tenv')
REFERENCES = ('no_change', 'delta_mean')
MODEL = 'perturbdiff'
PR_MODELS = (MODEL, 'delta_mean')
PR_FILENAME = 'pr_curves.json'
PR_SCHEMA = 'step04-pr-curves-v1'
SLOT_FIELDS = ('dataset', 'fold', 'seed', 'model')
GROUP_FIELDS = (*SLOT_FIELDS, 'group_id', 'n_wells', *METRICS,
                *(f'{metric}_reason' for metric in METRICS))
RUN_FIELDS = (*SLOT_FIELDS, 'status', 'reason', 'group_count', *METRICS,
              *(f'{metric}_valid_groups' for metric in METRICS))
TERMINAL_FAILURES = {'failed', 'timeout', 'incomplete', 'interrupted', 'cancelled'}


def verify_file(path, expected):
    if not valid_hash(expected) or not path.is_file() or sha256(path) != expected:
        raise ValueError(f'file checksum mismatch: {path}')
    return path


def valid_hash(value):
    return isinstance(value, str) and re.fullmatch('[0-9a-f]{64}', value) is not None


def relative_path(root, name):
    path = root / name
    if Path(name).is_absolute() or '..' in Path(name).parts or path.resolve() != path:
        raise ValueError(f'unsafe artifact path: {name}')
    return path


def load_schema(path, schema):
    document = read_json(path)
    if document['schema'] != schema:
        raise ValueError(f'expected {schema}: {path}')
    return document


def grouped_cases(cases):
    indices = defaultdict(list)
    case_ids = set()
    for index, case in enumerate(cases):
        case_id, group_id = case['case_id'], case['group_id']
        if not case_id or not group_id or case_id in case_ids:
            raise ValueError('empty or duplicate case identity')
        case_ids.add(case_id)
        indices[group_id].append(index)
    if not indices:
        raise ValueError('empty case inventory')
    return indices, case_ids


def fold_inputs(root, entry, stage):
    request_path = verify_file(relative_path(root, entry['request']), entry['request_sha256'])
    request = read_json(request_path)
    consumer, source = request['input'], Path(request['input_root'])
    identity = request['identity']
    dataset, fold, n_genes = identity['dataset'], identity['fold'], request['n_genes']
    if (consumer.keys() != INPUT_FIELDS or consumer['schema'] != 'step04-input-v1'
            or request['target_kind'] != 'delta' or request != request_for(source, consumer)
            or dataset not in DATASETS or type(fold) is not int or fold not in range(3)
            or type(n_genes) is not int or n_genes < 1 or (stage == 'formal' and n_genes != 18164)
            or source.resolve() != source):
        raise ValueError('invalid evaluator input identity')
    if stage == 'formal' and dataset == 'tenv':
        pin = read_json(REPO / STEP03 / 'execution/input_record/verification.json')['manifest_sha256']
        if identity['source_input_sha256'] != pin:
            raise ValueError('tenv release differs from the frozen Step 03 identity')
    manifest_path = verify_file(relative_path(source, 'manifest.json'), identity['source_input_sha256'])
    release = load_schema(manifest_path, f'{dataset}-input-v1')
    files = release['files']
    paths = tenv_paths(fold)
    required = ('genes', 'test_cases', 'train_cases', 'train_expression')
    if (release['target_kind'] != 'delta' or release['gene_count'] != n_genes
            or files['genes.csv'] != identity['gene_sha256']
            or files['folds.csv'] != identity['split_sha256']
            or any(consumer['predictor_files'][key] != paths[key]
                   or request['file_sha256'][key] != files[paths[key]] for key in required)):
        raise ValueError('evaluator input paths or release hashes differ from predictor request')
    cases_path = verify_file(relative_path(source, paths['test_cases']), files[paths['test_cases']])
    cases = rows(cases_path)
    groups, case_ids = grouped_cases(cases)
    recipes = entry['recipe_sha256_by_seed']
    seeds = tuple(seed for seed in SEEDS if str(seed) in recipes)
    if (not seeds or len(seeds) != len(recipes) or not all(map(valid_hash, recipes.values()))
            or (stage == 'formal' and seeds != SEEDS)):
        raise ValueError('invalid expected seed or recipe inventory')
    return {'identity': identity, 'root': source, 'files': files, 'paths': paths,
            'cases': cases, 'groups': groups, 'case_ids': case_ids, 'n_genes': n_genes,
            'recipes': {seed: recipes[str(seed)] for seed in seeds}}


def load_array(path, shape, dtype='float32'):
    array = np.load(path, mmap_mode='r', allow_pickle=False)
    if array.shape != shape or array.dtype != dtype:
        raise ValueError(f'invalid shape or dtype: {path}')
    return array


def prediction(root, outcome, fold, stage, producer_plan_sha=None):
    path = verify_file(relative_path(root, outcome['prediction']), outcome['prediction_manifest_sha256'])
    manifest = load_schema(path, 'step04-heldout-prediction-v1')
    seed = outcome['seed']
    identity = {**fold['identity'], 'seed': seed, 'stage': stage}
    recipe_hash = fold['recipes'][seed]
    run = read_json(verify_file(relative_path(path.parent, 'run.json'), manifest['run_sha256']))
    checkpoint_keys = ('identity', 'recipe_sha256', 'checkpoint_sha256', 'best_step',
                       'best_validation_loss', 'model_state_sha256')
    if (run['status'] != 'completed' or run['monitor'] != MONITOR
            or run['identity'] != identity or run['recipe_sha256'] != recipe_hash
            or json_hash(run['recipe']) != recipe_hash
            or any(run['recipe'][key] != identity[key] for key in ('seed', 'fold', 'stage'))
            or any(run[key] != manifest[key] for key in checkpoint_keys)
            or not valid_hash(run['model_state_sha256'])):
        raise ValueError('prediction checkpoint identity or frozen recipe mismatch')
    best = verify_best_checkpoint(run, 'prediction')
    verify_file(relative_path(path.parent, 'best.pt'), run['checkpoint_sha256'])
    if producer_plan_sha is not None:
        replay = read_json(verify_file(relative_path(path.parent, 'replay.json'), manifest['replay_sha256']))
        if (manifest['plan_sha256'] != producer_plan_sha
                or replay != {**run['reload'], **{key: run[key] for key in checkpoint_keys},
                              'plan_sha256': producer_plan_sha}
                or replay['validation_loss'] != best['validation_loss']
                or replay['prediction_max_abs_error'] != 0):
            raise ValueError('prediction replay does not match the frozen plan and EMA selection')
    cases = read_json(verify_file(relative_path(path.parent, 'cases.json'), manifest['cases_sha256']))
    if cases != [[case['case_id'], case['group_id']] for case in fold['cases']]:
        raise ValueError('prediction case and group identities differ from frozen order')
    array = verify_file(relative_path(path.parent, 'prediction.npy'), manifest['prediction_sha256'])
    return load_array(array, (len(cases), fold['n_genes']))


def mean_delta(array, indices):
    """Reduce one group in bounded slices, rejecting nonfinite rows."""
    total = np.zeros(array.shape[1], dtype=np.float64)
    for start in range(0, len(indices), 32):
        batch = array[indices[start:start + 32]]
        if not np.isfinite(batch).all():
            raise ValueError('nonfinite delta rows')
        total += batch.sum(axis=0, dtype=np.float64)
    return total / len(indices)


def evaluator_inputs(fold):
    root, files, paths, n_genes = (fold[key] for key in ('root', 'files', 'paths', 'n_genes'))

    def pinned(relative):
        return verify_file(relative_path(root, relative), files[relative])

    prefix = f"fold_{fold['identity']['fold']}/evaluator/"
    groups = rows(pinned(prefix + 'groups.csv'))
    test_groups, test_cases = fold['groups'], fold['cases']
    group_ids = [group['group_id'] for group in groups]
    if len(group_ids) != len(test_groups) or set(group_ids) != set(test_groups):
        raise ValueError('frozen evaluator group identities differ from test cases')
    # Tahoe's frozen group_id binds the cell line. Predictor case rows omit that column.
    for group in groups:
        members = test_groups[group['group_id']]
        if (int(group['n_wells']) != len(members)
                or any(group[key] != test_cases[index][key] for index in members for key in ('drug_key', 'dose_um'))):
            raise ValueError('frozen evaluator group membership differs from test cases')
    shape = (len(groups), n_genes)
    truth = load_array(pinned(prefix + 'truth.npy'), shape)
    labels = [load_array(pinned(prefix + name), shape, dtype='bool')
              for name in ('positive.npy', 'negative.npy')]
    for index, observed in enumerate(truth):
        if not np.isfinite(observed).all() or any(
                not np.array_equal(label[index], expected)
                for label, expected in zip(labels, truth_labels(observed))):
            raise ValueError('invalid frozen directional truth')
    training = rows(pinned(paths['train_cases']))
    train_groups, train_ids = grouped_cases(training)
    if train_ids & fold['case_ids'] or train_groups.keys() & test_groups.keys():
        raise ValueError('training and held-out cases or groups overlap')
    train_delta = load_array(pinned(paths['train_expression']), (len(training), n_genes))
    delta_mean = np.zeros(n_genes, dtype=np.float64)
    for members in train_groups.values():
        delta_mean += mean_delta(train_delta, members)
    delta_mean /= len(train_groups)
    return groups, truth, labels, {'no_change': np.zeros(n_genes), 'delta_mean': delta_mean}


def run_scores(writer, slot, status, reason, groups, truth, labels, indices, delta, pr_runs=None):
    raise ValueError('Old score releases are retired; use stage2/evaluation/evaluate.py with the frozen inventory')


def summaries(runs):
    buckets = defaultdict(list)
    for run in runs:
        buckets[run['dataset'], run['model']].append(run)
    for (dataset, model), slots in buckets.items():
        completed = sum(row['status'] == 'completed' for row in slots)
        row = {'dataset': dataset, 'model': model, 'expected_runs': len(slots),
               'completed_runs': completed, 'complete': completed == len(slots)}
        for metric in METRICS:
            valid = [slot[metric] for slot in slots if slot[metric] is not None]
            row[f'{metric}_mean'] = float(np.mean(valid)) if valid else None
            row[f'{metric}_sample_sd'] = float(np.std(valid, ddof=1)) if len(valid) > 1 else None
            row[f'{metric}_valid_runs'] = len(valid)
        yield row


def write_results(output, folds, outcomes, predictions, audit, include_pr):
    runs = []
    pr_runs = defaultdict(list)
    with (output / 'per_group.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=GROUP_FIELDS)
        writer.writeheader()
        for (dataset, fold_number), fold in sorted(folds.items()):
            groups, truth, labels, references = evaluator_inputs(fold)
            directory = output / 'references' / f'{dataset}_fold_{fold_number}'
            write_json(directory / 'group_ids.json', [group['group_id'] for group in groups])
            for name, delta in references.items():
                np.save(directory / f'{name}.npy', delta)
            for model, seed in [(MODEL, seed) for seed in fold['recipes']] + [
                    (reference, None) for reference in REFERENCES]:
                key = (dataset, fold_number, seed)
                outcome = outcomes.get(key, {'status': 'missing', 'reason': 'no_outcome'})
                status, reason = (outcome['status'], outcome.get('reason', '')) if seed else ('completed', '')
                delta = predictions.get(key) if seed else references[model]
                slot = {'dataset': dataset, 'fold': fold_number, 'model': model, 'seed': seed}
                curves = pr_runs[dataset] if include_pr and model in PR_MODELS else None
                runs.append(run_scores(writer, slot, status, reason, groups, truth, labels,
                                       fold['groups'], delta, curves))
    write_csv(output / 'per_run.csv', runs, RUN_FIELDS)
    write_csv(output / 'summary.csv', summaries(runs))
    if include_pr:
        for dataset, curves in pr_runs.items():
            write_json(output / dataset / PR_FILENAME, {
                'schema': PR_SCHEMA, 'dataset': dataset, 'recall': RECALL_GRID.tolist(),
                'plan_sha256': audit['plan_sha256'], 'outcomes_sha256': audit['outcomes_sha256'],
                'per_run_sha256': sha256(output / 'per_run.csv'), 'runs': curves,
            })
    audit['complete'] = all(run['status'] == 'completed' for run in runs)
    audit['files'] = {str(path.relative_to(output)): sha256(path)
                      for path in sorted(output.rglob('*')) if path.is_file()}
    write_json(output / 'evaluation.json', audit)


def evaluate(plan_path, outcomes_path, output, expected_plan_sha=None, include_pr=False):
    plan_sha = sha256(plan_path)
    if expected_plan_sha is not None and expected_plan_sha != plan_sha:
        raise ValueError('frozen plan hash mismatch')
    plan = load_schema(plan_path, 'step04-evaluation-plan-v1')
    stage = plan['stage']
    if stage not in ('fixture', 'formal'):
        raise ValueError('expected fixture or formal inventory')
    if stage == 'formal' and expected_plan_sha is None:
        raise ValueError('formal evaluation requires a frozen plan hash')
    scope = DATASETS
    if 'dataset' in plan:
        if stage != 'formal' or plan['dataset'] not in DATASETS:
            raise ValueError('dataset-scoped evaluation requires a valid formal dataset')
        scope = (plan['dataset'],)
    producer_plan_sha = plan_sha if stage == 'formal' or 'launch_files' in plan else None
    inventory = load_schema(outcomes_path, 'step04-evaluation-outcomes-v1')
    if producer_plan_sha is not None and inventory.get('plan_sha256') != producer_plan_sha:
        raise ValueError('outcomes do not match the frozen producer plan')
    folds = {}
    for entry in plan['folds']:
        fold = fold_inputs(plan_path.parent, entry, stage)
        identity = fold['identity']
        key = (identity['dataset'], identity['fold'])
        if key in folds:
            raise ValueError('duplicate expected dataset/fold')
        folds[key] = fold
    if (not folds or (stage == 'formal' and folds.keys() != {(d, f) for d in scope for f in range(3)})
            or len({fold['identity']['gene_sha256'] for fold in folds.values()}) != 1):
        raise ValueError('expected dataset/fold inventory or shared gene axis mismatch')
    for dataset in DATASETS:
        releases = {fold['identity']['source_input_sha256'] for key, fold in folds.items() if key[0] == dataset}
        if len(releases) > 1:
            raise ValueError(f'multiple source releases for {dataset}')
    expected = {(dataset, fold, seed) for (dataset, fold), inputs in folds.items()
                for seed in inputs['recipes']}
    outcomes, predictions = {}, {}
    for entry in inventory['outcomes']:
        key = (entry['dataset'], entry['fold'], entry['seed'])
        if type(entry['fold']) is not int or type(entry['seed']) is not int or key not in expected or key in outcomes:
            raise ValueError('unexpected or duplicate outcome slot')
        outcomes[key] = entry
        if entry['status'] == 'completed':
            predictions[key] = prediction(outcomes_path.parent, entry, folds[key[:2]], stage, producer_plan_sha)
        elif entry['status'] not in TERMINAL_FAILURES or not entry.get('reason') or 'prediction' in entry:
            raise ValueError('outcome requires a completed prediction or explicit terminal failure')
    if output.exists() or output.is_symlink():
        raise ValueError('evaluation requires a fresh output directory')
    output.parent.mkdir(parents=True, exist_ok=True)
    audit = {'schema': 'step04-evaluation-v1', 'stage': stage,
             'plan_sha256': plan_sha, 'outcomes_sha256': sha256(outcomes_path),
             'datasets': sorted({dataset for dataset, _ in folds}),
             'tenv_test_memberships_previously_scored_in_step03': stage == 'formal' and 'tenv' in scope}
    with tempfile.TemporaryDirectory(dir=output.parent) as temporary:
        pending = Path(temporary) / 'evaluation'
        pending.mkdir()
        write_results(pending, folds, outcomes, predictions, audit, include_pr)
        pending.rename(output)


def main():
    sys.path.insert(0, str(REPO))
    from stage2.evaluation.evaluate import main as shared_main
    sys.argv.extend(['--step', '4'])
    shared_main()


if __name__ == '__main__':
    main()
