#!/usr/bin/env python3
"""Freeze, fit, predict and score the Step 05 pooled Mix arm."""
import argparse
import csv
from pathlib import Path
import sys
import tempfile

sys.dont_write_bytecode = True

import cold_tahoe
import warm_tenv
from artifacts import SEEDS, read_json, rows, sha256, write_csv, write_json
from evaluate import (GROUP_FIELDS, PR_SCHEMA, RECALL_GRID, RUN_FIELDS, evaluator_inputs, fold_inputs,
                      grouped_cases, load_array, np, prediction, run_scores, summaries)
from formal import completed_fit, fit_slot, pinned, predict_slot, slot_lock, use_runtime, verify_fit_seal
from training_inputs import json_hash, verify_files, verify_request

REPO = cold_tahoe.REPO
STEP04, STEP05 = cold_tahoe.STEP04, cold_tahoe.STEP05
REGIMES = ('warm-start', 'cold-cm-drug', 'cold-drug')
SCHEMA = 'step05-mix-plan-v1'
TAHOE_FILES = ('train_cases', 'train_expression', 'controls', 'morgan', 'drugs', 'genes')


def historical_inputs(runtime, regime, cache):
    if 'base' not in cache:
        record = REPO / STEP04 / 'execution/formal_record/tenv/plan.json'
        formal = runtime / STEP04 / 'execution/output/tenv-formal-ticket72'
        plan = read_json(pinned(formal, 'plan.json', sha256(record)))
        entry = plan['folds'][0]
        fixed, = (value for path, value in plan['launch_files'].items()
                  if path.endswith('execution/models/perturbdiff/v2.config.json'))
        pinned(runtime, 'stage2/03_modes_tenv/execution/models/perturbdiff/v2.config.json', fixed)
        cache['base'] = formal, plan, entry
    formal, plan, entry = cache['base']
    execution = runtime / STEP05 / 'execution/output'
    split = execution / 'cold-tenv-ticket76/split.json'
    if regime == 'warm-start':
        record = read_json(REPO / STEP05 / 'execution/warm_record/replay.json')
        root = execution / 'warm-tenv-ticket77'
        digest = record['local_artifact_sha256'][f'{root.name}/plan.json']
        warm_plan, tenv = warm_tenv.verified_plan(root, digest, runtime)
        cache['warm_plan'] = warm_plan
    else:
        if 'tenv' not in cache:
            tenv = read_json(pinned(formal, entry['request'], entry['request_sha256']))
            verify_request(tenv, runtime)
            record = read_json(REPO / STEP05 / 'execution/cold_record/replay.json')
            key = f'{split.parent.name}/{split.name}'
            pinned(split.parent, split.name, record['local_artifact_sha256'][key])
            if tenv['identity'] != record['split_identity']:
                raise ValueError('cold target differs from the verified Step 05 split')
            cache['tenv'] = tenv
        tenv = cache['tenv']

    tahoe_formal = runtime / STEP04 / 'execution/output/tahoe-formal-ticket73'
    if regime == 'cold-drug':
        root = execution / 'cold-tahoe-ticket79'
        record = read_json(REPO / STEP05 / 'execution/cold_tahoe_record/qualification.json')
        pinned(root, 'release.json', record['local_artifact_sha256'][f'{root.name}/release.json'])
        tahoe = cold_tahoe.verify_release(runtime, tahoe_formal, split, root)
    else:
        if 'tahoe' not in cache:
            root, consumer, _, _, _ = cold_tahoe.frozen_inputs(runtime, tahoe_formal, split)
            verify_files(root, {key: consumer['predictor_files'][key] for key in TAHOE_FILES},
                         consumer['file_sha256'])
            cache['tahoe'] = cold_tahoe.request_for(root, consumer)
        tahoe = cache['tahoe']
    return plan, entry, tenv, tahoe


def pool_request(regime, tenv, tahoe):
    n_genes = tenv['n_genes']
    if (tahoe['n_genes'] != n_genes or tenv['identity']['gene_sha256'] != tahoe['identity']['gene_sha256']):
        raise ValueError('pooled sources must share the ordered gene axis')
    sources = {'tenv': tenv['predictor_files'],
               'tahoe': {key: tahoe['predictor_files'][key] for key in TAHOE_FILES}}
    inventory, dose_union = {}, set()
    for name, files in sources.items():
        cases = rows(files['train_cases'])
        groups, ids = grouped_cases(cases)
        doses = {float(case['dose_um']) for case in cases}
        if not all(np.isfinite(dose) and dose > 0 for dose in doses):
            raise ValueError('invalid training dose')
        for members in groups.values():
            if len({(cases[i]['drug_key'], float(cases[i]['dose_um'])) for i in members}) != 1:
                raise ValueError('pooled source group spans multiple drug doses')
        dose_union.update(doses)
        inventory[name] = {'train_cases': len(ids), 'train_groups': len(groups),
                           'ordered_namespaced_groups_sha256': json_hash([f'{name}:{g}' for g in groups]),
                           'dose_um': sorted(doses)}

    for role in ('validation', 'test'):
        cases = rows(sources['tenv'][f'{role}_cases'])
        if not {float(case['dose_um']) for case in cases} <= dose_union:
            raise ValueError('tenv holdout introduces a dose outside the Mix training vocabulary')
    payload = {
        'schema': 'step05-mix-request-v1', 'n_genes': n_genes, **sources,
        'dose_um': sorted(dose_union), 'train_inventory': inventory,
        'source_identities': {'tenv': tenv['identity'], 'tahoe': tahoe['identity']},
        'predictor_sha256': {
            'tenv': tenv['file_sha256'], 'tahoe': {key: tahoe['file_sha256'][key] for key in TAHOE_FILES},
        },
    }
    return {**payload, 'identity': {**tenv['identity'], 'input_sha256': json_hash(payload),
                                    'regime': regime, 'arm': 'Mix'}}


def manifest(runtime, stage, cache, regimes=REGIMES):
    requests = {}
    for regime in regimes:
        base, entry, tenv, tahoe = historical_inputs(runtime, regime, cache)
        requests[regime] = pool_request(regime, tenv, tahoe)
        cache[regime] = tenv
    formal = cache['base'][0]
    recipes = {}
    for seed in SEEDS:
        item = entry['recipes'][str(seed)]
        raw = read_json(pinned(formal, item['path'], item['sha256']))
        recipes[str(seed)] = {'sha256': item['sha256'], 'recipe': raw,
                              'config_sha256': entry['recipe_sha256_by_seed'][str(seed)]}
    active = recipes if stage == 'formal' else {
        '17': {'recipe': {**recipes['17']['recipe'], **warm_tenv.SMOKE}}}
    return {
        'schema': SCHEMA, 'stage': stage, 'formal_recipes': recipes, 'active_recipes': active,
        'inference': base['inference'],
        'requests': requests,
        'formal_slots': [{**request['identity'], 'seed': seed, 'stage': 'formal',
                          'recipe_sha256': recipes[str(seed)]['config_sha256']}
                         for request in requests.values() for seed in SEEDS],
    }


def prepare(runtime, root, stage):
    if root.exists() or root.is_symlink():
        raise ValueError('prepare requires a fresh output directory')
    plan = manifest(runtime, stage, {})
    root.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=root.parent, prefix='.mix-') as temporary:
        pending = Path(temporary) / 'mix'
        write_json(pending / 'plan.json', plan)
        pending.rename(root)
    print(f'Frozen Mix plan SHA-256 {sha256(root / "plan.json")}')


def verify(runtime, root, digest, regime=None):
    plan = read_json(pinned(root, 'plan.json', digest))
    if plan['schema'] != SCHEMA or plan['stage'] not in ('smoke', 'formal'):
        raise ValueError('unexpected Mix plan schema or stage')
    if set(plan['requests']) != set(REGIMES):
        raise ValueError('Mix plan must contain all three regimes')
    cache = {}
    expected = manifest(runtime, plan['stage'], cache, (regime,) if regime else REGIMES)
    selected = dict(plan)
    if regime:
        selected['requests'] = {regime: plan['requests'][regime]}
        selected['formal_slots'] = [slot for slot in plan['formal_slots'] if slot['regime'] == regime]
    if selected != expected:
        raise ValueError('Mix plan differs from the frozen source, split, or recipe inventory')
    return plan, cache


def slot_for(root, regime, seed):
    return root / 'slots' / regime / f'seed-{seed}'


def worker(runtime, root, digest, plan, regime, seed, command):
    use_runtime(runtime)
    import native_training
    from mix_inputs import MixInputs

    native_training.PredictorInputs = MixInputs
    request = plan['requests'][regime]
    item = plan['active_recipes'][str(seed)]
    identity = {**request['identity'], 'seed': seed, 'stage': plan['stage']}
    slot = slot_for(root, regime, seed)
    with slot_lock(slot), tempfile.TemporaryDirectory(dir=slot) as temporary:
        recipe = Path(temporary) / 'recipe.json'
        write_json(recipe, item['recipe'])
        config = native_training.load_recipe(runtime, recipe, request)
        if 'config_sha256' in item and json_hash(config) != item['config_sha256']:
            raise ValueError('Mix recipe differs from the Step 04 CM recipe')
        if command == 'fit':
            fit_slot(slot, digest, runtime, request, identity, recipe, config)
        else:
            predict_slot(slot, plan, digest, request, identity, config)


def heldout_predictions(root, digest, plan, regime, fold):
    outcomes, predictions = [], {}
    for seed_text, item in plan['active_recipes'].items():
        seed = int(seed_text)
        slot = slot_for(root, regime, seed)
        verify_fit_seal(slot / 'fit', digest)
        run = read_json(slot / 'fit/run.json')
        formal = plan['formal_recipes'][seed_text]
        if (any(run['recipe'][key] != value for key, value in item['recipe'].items())
                or json_hash({**run['recipe'], **formal['recipe']}) != formal['config_sha256']):
            raise ValueError('Mix fit did not use the selected CM recipe')
        identity = {**fold['identity'], 'seed': seed, 'stage': plan['stage']}
        completed_fit(slot / 'fit', identity, run['recipe'])
        path = slot / 'prediction/prediction.json'
        seal = read_json(path)
        if (seal['run_sha256'] != sha256(slot / 'fit/run.json') or seal['case_count'] != len(fold['cases'])
                or seal['group_count'] != len(fold['groups'])):
            raise ValueError('Mix held-out prediction does not match the sealed fit and frozen cases')
        fold['recipes'][seed] = item.get('config_sha256', run['recipe_sha256'])
        outcome = {'regime': regime, 'fold': 0, 'arm': 'Mix', 'seed': seed, 'status': 'completed',
                   'prediction': str(path.relative_to(root)), 'prediction_manifest_sha256': sha256(path)}
        values = prediction(root, outcome, fold, plan['stage'], digest)
        if any(not np.isfinite(values[start:start + 32]).all() for start in range(0, len(values), 32)):
            raise ValueError('nonfinite Mix held-out prediction')
        predictions[seed] = values
        outcomes.append(outcome)
    return outcomes, predictions


def score(root, digest, plan, cache, regime, output):
    if output.exists() or output.is_symlink():
        raise ValueError('score requires a fresh output directory')
    request, tenv = plan['requests'][regime], cache[regime]
    cases = rows(tenv['predictor_files']['test_cases'])
    groups, ids = grouped_cases(cases)
    fold = {'identity': request['identity'], 'cases': cases, 'groups': groups,
            'case_ids': ids, 'n_genes': request['n_genes'], 'recipes': {}}
    outcomes, predictions = heldout_predictions(root, digest, plan, regime, fold)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent, prefix='.mix-score-') as temporary:
        work = Path(temporary)
        # The evaluator resolves treated expression only after all predictions pass.
        if regime == 'warm-start':
            warm = cache['warm_plan']
            source = Path(warm['source_root'])
            release = read_json(source / 'manifest.json')
            locations = warm_tenv.expression_sources(source, release['files'], cases)
            warm_tenv.grouped_truth(work, groups, cases, locations)
            inventory = rows(work / 'groups.csv')
            shape = (len(inventory), fold['n_genes'])
            truth = load_array(work / 'truth.npy', shape)
            labels = [load_array(work / f'{direction}.npy', shape, 'bool')
                      for direction in ('positive', 'negative')]
            evaluator_sha = {path.name: sha256(path) for path in sorted(work.iterdir())}
        else:
            formal, _, entry = cache['base']
            frozen = fold_inputs(formal, entry, 'formal')
            inventory, truth, labels, _ = evaluator_inputs(frozen)
            evaluator_sha = {Path(name).name: value for name, value in frozen['files'].items()
                             if name.startswith('fold_0/evaluator/')}
        pending = work / 'scores'
        pending.mkdir()
        runs, curves = [], []
        with (pending / 'per_group.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=GROUP_FIELDS)
            writer.writeheader()
            for seed, values in predictions.items():
                slot = {'dataset': 'tenv', 'fold': 0, 'model': 'Mix', 'seed': seed}
                runs.append(run_scores(writer, slot, 'completed', '', inventory, truth, labels,
                                       groups, values, curves))
        write_csv(pending / 'per_run.csv', runs, RUN_FIELDS)
        write_csv(pending / 'summary.csv', summaries(runs))
        outcome_sha = json_hash(outcomes)
        write_json(pending / 'tenv/pr_curves.json', {
            'schema': PR_SCHEMA, 'dataset': 'tenv', 'recall': RECALL_GRID.tolist(),
            'plan_sha256': digest, 'outcomes_sha256': outcome_sha,
            'per_run_sha256': sha256(pending / 'per_run.csv'), 'runs': curves,
        })
        write_json(pending / 'evaluation.json', {
            'schema': 'step04-evaluation-v1', 'stage': plan['stage'], 'regime': regime, 'datasets': ['tenv'],
            'plan_sha256': digest, 'outcomes_sha256': outcome_sha, 'outcomes': outcomes,
            'evaluator_files': evaluator_sha, 'complete': True,
            'files': {str(path.relative_to(pending)): sha256(path)
                      for path in sorted(pending.rglob('*')) if path.is_file()},
        })
        pending.rename(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime', type=Path, default=REPO)
    commands = parser.add_subparsers(dest='command', required=True)
    build = commands.add_parser('prepare')
    build.add_argument('root', type=Path)
    build.add_argument('--smoke', action='store_true')
    for name in ('verify', 'fit', 'predict', 'score'):
        command = commands.add_parser(name)
        command.add_argument('root', type=Path)
        command.add_argument('--plan-sha256', required=True)
        if name != 'verify':
            command.add_argument('regime', choices=REGIMES)
        if name in ('fit', 'predict'):
            command.add_argument('seed', type=int, choices=SEEDS)
        if name == 'score':
            command.add_argument('output', type=Path)
    args = parser.parse_args()
    runtime, root = args.runtime.resolve(), args.root.absolute()
    try:
        if args.command == 'prepare':
            prepare(runtime, root, 'smoke' if args.smoke else 'formal')
            return
        regime = getattr(args, 'regime', None)
        plan, cache = verify(runtime, root, args.plan_sha256, regime)
        if args.command == 'score':
            score(root, args.plan_sha256, plan, cache, regime, args.output.absolute())
        elif args.command != 'verify':
            worker(runtime, root, args.plan_sha256, plan, regime, args.seed, args.command)
    except (ArithmeticError, KeyError, ValueError, TypeError, IndexError, OSError, EOFError) as error:
        parser.exit(1, f'Mix {args.command} rejected: {error}\n')


if __name__ == '__main__':
    main()
