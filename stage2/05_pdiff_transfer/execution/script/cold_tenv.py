#!/usr/bin/env python3
"""Verify and score the Step 04 tenv fold shared by both Step 05 cold regimes."""
import argparse
from pathlib import Path
import shutil
import sys
import tempfile

sys.dont_write_bytecode = True
REPO = next(path for path in Path(__file__).resolve().parents if (path / '.git').exists())
STEP04 = REPO / 'stage2/04_pdiff_hyperparameter'
RECORD = STEP04 / 'execution/formal_record/tenv'
sys.path.insert(0, str(STEP04 / 'execution/script'))

from evaluate import (SEEDS, fold_inputs, grouped_cases, load_array, load_schema, prediction,
                      read_json, relative_path, rows, sha256, verify_file, write_json, write_results)
from training_inputs import verify_request

REGIMES = ('cold-cm-drug', 'cold-drug')
ROLES = ('train', 'validation', 'test')
CHECKPOINT_FIELDS = ('identity', 'recipe_sha256', 'checkpoint_sha256', 'model_state_sha256',
                     'best_step', 'best_validation_loss', 'run_sha256', 'replay_sha256')


def pinned(root, name, digest):
    return verify_file(relative_path(root, name), digest)


def unique_rows(entries, key, subject):
    by_key = {entry[key]: entry for entry in entries}
    if len(by_key) != len(entries) or not all(by_key):
        raise ValueError(f'empty or duplicate {subject} identity')
    return by_key


def verify_split(request, fold):
    root, files, n_genes = (fold[key] for key in ('root', 'files', 'n_genes'))

    def source_rows(name):
        return rows(pinned(root, name, files[name]))

    genes = unique_rows(rows(root / 'genes.csv'), 'symbol', 'gene')
    if len(genes) != n_genes:
        raise ValueError('unexpected ordered gene axis length')
    control_rows = source_rows('controls.csv')
    load_array(root / 'controls.npy', (len(control_rows), n_genes))
    if ([int(row['control_index']) for row in control_rows] != list(range(len(control_rows)))
            or any(int(row['n_wells']) <= 0 for row in control_rows)):
        raise ValueError('invalid matched-control inventory')

    drugs = unique_rows(rows(root / 'drugs.csv'), 'drug_key', 'drug')
    memberships = unique_rows([row for row in rows(root / 'folds.csv') if row['fold'] == '0'],
                              'drug_key', 'fold membership')
    if (memberships.keys() != drugs.keys()
            or any(row['role'] not in ROLES or not row['parent_connectivity']
                   or row['parent_connectivity'] != drugs[drug]['parent_connectivity']
                   for drug, row in memberships.items())):
        raise ValueError('fold memberships differ from the drug parent-connectivity inventory')
    wells = source_rows('wells.csv')
    unique_rows(wells, 'case_id', 'well')
    roles = {}
    assigned_parents = set()
    for role in ROLES:
        cases = rows(root / fold['paths'][f'{role}_cases'])
        expected = [well for well in wells if memberships[well['drug_key']]['role'] == role]
        parents = {row['parent_connectivity'] for row in memberships.values() if row['role'] == role}
        if cases != expected or any(not 0 <= int(case['control_index']) < len(control_rows) for case in cases):
            raise ValueError(f'{role} case order, membership or matched controls differ from the frozen release')
        if parents & assigned_parents:
            raise ValueError(f'parent-connectivity overlap at {role}')
        assigned_parents.update(parents)
        groups, _ = grouped_cases(cases)
        roles[role] = {'case_count': len(cases), 'group_count': len(groups),
                       'parent_connectivity': sorted(parents)}
    return {'schema': 'step05-cold-tenv-split-v1', **{key: request[key] for key in
            ('identity', 'input_root', 'predictor_files', 'file_sha256', 'n_genes', 'target_kind')},
            'release_metadata_sha256': {name: files[name] for name in ('folds.csv', 'wells.csv', 'controls.csv')},
            'roles': roles}


def verify_reused_cm(formal, release, index, entry, fold, outcome, plan_sha):
    seed = outcome['seed']
    recipe = entry['recipes'][str(seed)]
    pinned(formal, recipe['path'], recipe['sha256'])
    if (index['directory'] != f'fold-0/seed-{seed}'
            or index['identity'] != {**fold['identity'], 'seed': seed, 'stage': 'formal'}
            or index['recipe_sha256'] != fold['recipes'][seed]
            or index['raw_recipe_sha256'] != recipe['sha256']
            or index['policy_sha256'] != entry['policy_sha256']
            or index['request_sha256'] != entry['request_sha256']):
        raise ValueError(f'CM checkpoint, seed or recipe index mismatch for seed {seed}')
    slot = relative_path(release, index['directory'])
    for name, field in (('run.json', 'run_sha256'), ('replay.json', 'replay_sha256'),
                        ('request.json', 'request_sha256'), ('recipe.json', 'raw_recipe_sha256')):
        pinned(slot, name, index[field])
    manifest_path = pinned(formal, outcome['prediction'], outcome['prediction_manifest_sha256'])
    manifest = read_json(manifest_path)
    if (any(manifest[key] != index[key] for key in CHECKPOINT_FIELDS)
            or manifest['case_count'] != len(fold['cases'])
            or manifest['group_count'] != len(fold['groups'])):
        raise ValueError(f'CM prediction differs from the indexed checkpoint or case inventory for seed {seed}')
    values = prediction(formal, outcome, fold, 'formal', plan_sha)
    checkpoint = relative_path(slot, 'best.pt')
    if not checkpoint.samefile(manifest_path.parent / 'best.pt'):
        verify_file(checkpoint, index['checkpoint_sha256'])
    return values, {**index, 'checkpoint_directory': str(slot), 'prediction_manifest': str(manifest_path),
                    'prediction_manifest_sha256': outcome['prediction_manifest_sha256'],
                    'prediction_sha256': manifest['prediction_sha256'], 'cases_sha256': manifest['cases_sha256']}


def compare_scores(output, evaluation):
    files = evaluation['files']
    expected = rows(pinned(RECORD, 'per_run.csv', files['per_run.csv']))
    if rows(output / 'per_run.csv') != [row for row in expected if row['fold'] == '0']:
        raise ValueError('fold 0 CM or reference scores differ from the Step 04 record')
    for path in sorted((output / 'references/tenv_fold_0').iterdir()):
        verify_file(path, files[str(path.relative_to(output))])
    reference_pr = read_json(pinned(RECORD, 'pr_curves.json', files.get('tenv/pr_curves.json')
                                   or sha256(RECORD / 'pr_curves.json')))
    replay_pr = read_json(output / 'tenv/pr_curves.json')
    if (replay_pr['recall'] != reference_pr['recall']
            or replay_pr['runs'] != [row for row in reference_pr['runs'] if row['fold'] == 0]):
        raise ValueError('fold 0 directional PR curves differ from the Step 04 record')


def replay(formal, release, output):
    if output.exists() or output.is_symlink():
        raise ValueError('cold comparison requires a fresh output directory')
    plan_sha = sha256(RECORD / 'plan.json')
    outcomes_sha = sha256(RECORD / 'outcomes.json')
    plan = load_schema(pinned(formal, 'plan.json', plan_sha), 'step04-evaluation-plan-v1')
    inventory = load_schema(pinned(formal, 'outcomes.json', outcomes_sha), 'step04-evaluation-outcomes-v1')
    evaluation_path = RECORD / 'evaluation.json'
    evaluation = load_schema(evaluation_path, 'step04-evaluation-v1')
    index = load_schema(pinned(release, 'checkpoints.json', sha256(RECORD / 'checkpoints.json')),
                        'step04-tenv-checkpoints-v1')
    if ((plan['dataset'], plan['stage']) != ('tenv', 'formal')
            or inventory['plan_sha256'] != plan_sha or not evaluation['complete']
            or evaluation['datasets'] != ['tenv'] or evaluation['stage'] != 'formal'
            or any(doc['plan_sha256'] != plan_sha or doc['outcomes_sha256'] != outcomes_sha
                   for doc in (index, evaluation))
            or index['evaluation_sha256'] != sha256(evaluation_path)
            or index['selection_sha256'] != plan['selection_sha256']):
        raise ValueError('Step 04 formal evidence identities do not agree')
    pinned(formal, 'selection.json', plan['selection_sha256'])
    entry = plan['folds'][0]
    request_path = pinned(formal, entry['request'], entry['request_sha256'])
    request = read_json(request_path)
    if (request['identity']['dataset'], request['identity']['fold']) != ('tenv', 0):
        raise ValueError('expected the frozen Step 04 tenv fold 0')
    verify_request(request, REPO)
    fold = fold_inputs(formal, entry, 'formal')
    split = verify_split(request, fold)
    outcomes = unique_rows([row for row in inventory['outcomes'] if row['fold'] == 0], 'seed', 'outcome')
    checkpoints = {row['identity']['seed']: row for row in index['checkpoints'] if row['identity']['fold'] == 0}
    if (outcomes.keys() != set(SEEDS) or checkpoints.keys() != set(SEEDS)
            or sum(row['identity']['fold'] == 0 for row in index['checkpoints']) != len(SEEDS)
            or any((row['dataset'], row['status']) != ('tenv', 'completed') for row in outcomes.values())):
        raise ValueError('CM reuse requires all three completed tenv fold 0 seeds')
    cm, predictions, scored_outcomes = {}, {}, {}
    for seed in SEEDS:
        slot = ('tenv', 0, seed)
        predictions[slot], cm[str(seed)] = verify_reused_cm(
            formal, release, checkpoints[seed], entry, fold, outcomes[seed], plan_sha)
        scored_outcomes[slot] = outcomes[seed]

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as temporary:
        pending = Path(temporary) / 'comparison'
        pending.mkdir()
        # Only the evaluator opens held-out expression, after all prediction seals pass.
        write_results(pending, {('tenv', 0): fold}, scored_outcomes, predictions,
                      {'schema': 'step04-evaluation-v1', 'stage': 'formal', 'datasets': ['tenv'],
                       'plan_sha256': plan_sha, 'outcomes_sha256': outcomes_sha}, include_pr=True)
        compare_scores(pending, evaluation)
        write_json(pending / 'split.json', split)
        shutil.copyfile(request_path, pending / 'predictor_request.json')
        shared = {'split': 'split.json', 'scores': 'per_run.csv', 'references': 'references/tenv_fold_0'}
        write_json(pending / 'comparison.json', {
            'schema': 'step05-cold-tenv-comparison-v1', 'regimes': dict.fromkeys(REGIMES, shared),
            'cm': cm, 'step04_plan_sha256': plan_sha, 'step04_outcomes_sha256': outcomes_sha,
            'step04_checkpoint_index_sha256': sha256(RECORD / 'checkpoints.json'),
            'step04_evaluation_sha256': sha256(evaluation_path), 'step04_score_and_pr_replay': 'exact',
            'files': {str(path.relative_to(pending)): sha256(path)
                      for path in sorted(pending.rglob('*')) if path.is_file()},
        })
        pending.rename(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('formal', type=Path, help='Step 04 local tenv formal run root')
    parser.add_argument('checkpoints', type=Path, help='Step 04 local tenv checkpoint release')
    parser.add_argument('output', type=Path, help='fresh comparison destination')
    args = parser.parse_args()
    try:
        replay(args.formal.resolve(), args.checkpoints.resolve(), args.output.absolute())
    except (KeyError, ValueError, TypeError, IndexError, OSError, EOFError) as error:
        parser.exit(1, f'cold tenv comparison rejected: {error}\n')


if __name__ == '__main__':
    main()
