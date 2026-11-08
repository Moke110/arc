#!/usr/bin/env python3
"""Prepare, verify, fit, predict and score the Step 05 warm-start tenv input."""
import argparse
from collections import defaultdict
from decimal import Decimal
from errno import EXDEV
import os
from pathlib import Path
from random import Random
import shutil
import sys
import tempfile

sys.dont_write_bytecode = True
REPO = next(path for path in Path(__file__).resolve().parents if (path / '.git').exists())
STEP04 = Path('stage2/04_pdiff_hyperparameter')
sys.path.insert(0, str(REPO / STEP04 / 'execution/script'))

from cold_tenv import pinned, unique_rows
from evaluate import (SEEDS, grouped_cases, load_array, load_schema, np, prediction,
                      read_json, rows, sha256, truth_labels, write_json, write_results)
from artifacts import write_csv
from training_inputs import STEP03, predictor_paths, request_for

ROLES = ('train', 'validation', 'test')
SPLIT_SEED = 20261105
N_GENES = 18164
SCHEMA = 'step05-warm-tenv-plan-v1'
SMOKE = {'stage': 'smoke', 'max_steps': 12, 'validation_interval': 6, 'minimum_steps': 12,
         'patience_steps': 12, 'max_training_seconds': 900}


def release(runtime, source):
    record = runtime / STEP03 / 'execution/input_record/verification.json'
    digest = read_json(record)['manifest_sha256']
    manifest = load_schema(pinned(source, 'manifest.json', digest), 'tenv-input-v1')
    if (manifest['gene_count'], manifest['target_kind']) != (N_GENES, 'delta'):
        raise ValueError('unexpected Step 03 tenv release contract')
    return digest, manifest


def membership(source, manifest):
    files = manifest['files']
    drugs = unique_rows(rows(pinned(source, 'drugs.csv', files['drugs.csv'])), 'drug_key', 'drug')
    wells = rows(pinned(source, 'wells.csv', files['wells.csv']))
    unique_rows(wells, 'case_id', 'well')
    group_keys = {}
    combinations = defaultdict(set)
    keyed = []
    for well in wells:
        parent = drugs[well['drug_key']]['parent_connectivity']
        dose = Decimal(well['dose_um'])
        if not parent or not dose.is_finite() or dose <= 0:
            raise ValueError('invalid frozen parent connectivity or dose')
        key = (parent, dose)
        group = (well['drug_key'], dose)
        if group_keys.setdefault(well['group_id'], group) != group:
            raise ValueError('frozen group spans multiple drug–dose combinations')
        keyed.append((key, well))
        combinations[parent].add(dose)

    # The frozen tenv inventory has the same three doses for every parent.
    # Shuffling parents and cycling shuffled doses balances test counts exactly to within one.
    parents = sorted(combinations)
    doses = sorted(set().union(*combinations.values()))
    if len(doses) < 3 or any(combinations[parent] != set(doses) for parent in parents):
        raise ValueError('warm-start requires the complete common tenv dose grid')
    rng = Random(SPLIT_SEED)
    rng.shuffle(parents)
    rng.shuffle(doses)
    test = {(parent, doses[index % len(doses)]) for index, parent in enumerate(parents)}
    candidates = [(parent, dose) for parent in sorted(parents) for dose in sorted(doses)
                  if (parent, dose) not in test]
    rng.shuffle(candidates)
    validation = {}
    target = round(0.15 * len(candidates))
    for parent, dose in candidates:
        if len(validation) == target:
            break
        validation.setdefault(parent, dose)
    val = set(validation.items())
    train = set(candidates) - val
    if (len(val) != target or not val or {parent for parent, _ in train} != set(parents)
            or {dose for _, dose in train} != set(doses)):
        raise ValueError('warm-start validation count or training coverage failed')
    by_role = {role: [] for role in ROLES}
    assignment = {key: role for role, keys in zip(ROLES, (train, val, test)) for key in keys}
    for key, well in keyed:
        by_role[assignment[key]].append(well)
    split = {'schema': 'step05-warm-tenv-split-v1', 'seed': SPLIT_SEED, 'fold': 0,
             'regime': 'warm-start',
             'combinations': [{'parent_connectivity': parent, 'dose_um': format(dose.normalize(), 'f'),
                               'role': assignment[parent, dose]} for parent, dose in sorted(assignment)],
             'roles': {role: {'case_count': len(cases), 'group_count': len(grouped_cases(cases)[0])}
                       for role, cases in by_role.items()}}
    return split, by_role


def expression_sources(source, files, selected):
    """Resolve only requested cases to frozen per-well deltas in legacy train/validation arrays."""
    missing = unique_rows(selected, 'case_id', 'selected case')
    locations = {}
    for fold in range(3):
        for role in ('train', 'validation'):
            path = f'fold_{fold}/{role}.csv'
            cases = rows(pinned(source, path, files[path]))
            matches = []
            for index, case in enumerate(cases):
                case_id = case['case_id']
                if case_id in missing:
                    if missing.pop(case_id) != case:
                        raise ValueError('legacy expression row differs from canonical well identity')
                    matches.append((case_id, index))
            if matches:
                expression = f'fold_{fold}/{role}_expression.npy'
                array = load_array(pinned(source, expression, files[expression]), (len(cases), N_GENES))
                locations.update({case_id: (array, index) for case_id, index in matches})
            if not missing:
                return locations
    raise ValueError('missing per-well expression in the Step 03 release')


def make_request(root, manifest, source_digest, split_sha, hashes):
    return request_for(root, {'schema': 'step04-input-v1', 'dataset': 'tenv', 'fold': 0,
                              'target_kind': 'delta', 'n_genes': N_GENES,
                              'source_input_sha256': source_digest,
                              'corpus_sha256': manifest.get('source_sha256'),
                              'preprocessing': manifest.get('expression'),
                              'split_sha256': split_sha, 'gene_sha256': hashes['genes'],
                              'predictor_files': predictor_paths(0), 'file_sha256': hashes})


def prepare(runtime, formal, output, smoke):
    if output.exists() or output.is_symlink():
        raise ValueError('prepare requires a fresh output directory')
    record = runtime / STEP04 / 'execution/formal_record/tenv/plan.json'
    base_digest = sha256(record)
    base = load_schema(pinned(formal, 'plan.json', base_digest), 'step04-evaluation-plan-v1')
    if (base['dataset'], base['stage']) != ('tenv', 'formal'):
        raise ValueError('warm-start requires the frozen formal tenv plan')
    entry = base['folds'][0]
    base_request = read_json(pinned(formal, entry['request'], entry['request_sha256']))
    source = Path(base_request['input_root']).resolve()
    source_digest, manifest = release(runtime, source)
    split, by_role = membership(source, manifest)
    locations = expression_sources(source, manifest['files'], by_role['train'] + by_role['validation'])
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as temp:
        pending = Path(temp) / 'warm'
        pending.mkdir()
        paths = predictor_paths(0)
        for key in ('genes', 'drugs', 'morgan', 'controls'):
            name = paths[key]
            shutil.copyfile(pinned(source, name, manifest['files'][name]), pending / name)
        for role, cases in by_role.items():
            write_csv(pending / paths[f'{role}_cases'], cases)
            if role != 'test':
                array = np.lib.format.open_memmap(pending / paths[f'{role}_expression'], mode='w+',
                                                  dtype='float32', shape=(len(cases), N_GENES))
                for index, case in enumerate(cases):
                    origin, source_index = locations[case['case_id']]
                    array[index] = origin[source_index]
                array.flush()
                del array
        write_json(pending / 'split.json', split)
        split_sha = sha256(pending / 'split.json')
        hashes = {key: sha256(pending / path) for key, path in paths.items()}
        write_json(pending / 'request.json', make_request(output, manifest, source_digest, split_sha, hashes))

        recipes = {}
        for seed in map(str, SEEDS):
            item = entry['recipes'][seed]
            path = f'recipes/seed-{seed}.json'
            (pending / path).parent.mkdir(exist_ok=True)
            shutil.copyfile(pinned(formal, item['path'], item['sha256']), pending / path)
            recipes[seed] = {'path': path, 'sha256': item['sha256'],
                             'config_sha256': entry['recipe_sha256_by_seed'][seed]}
        active = recipes
        if smoke:
            smoke_recipe = {**read_json(pending / recipes['17']['path']), **SMOKE}
            write_json(pending / 'smoke.json', smoke_recipe)
            active = {'17': {'path': 'smoke.json', 'sha256': sha256(pending / 'smoke.json')}}
        write_json(pending / 'plan.json', {
            'schema': SCHEMA, 'stage': 'smoke' if smoke else 'formal', 'regime': 'warm-start',
            'source_root': str(source), 'source_sha256': source_digest, 'split_sha256': split_sha,
            'request_sha256': sha256(pending / 'request.json'), 'step04_plan_sha256': base_digest,
            'formal_recipes': recipes, 'active_recipes': active, 'inference': base['inference'],
        })
        pending.rename(output)
    print(f'Frozen warm-start plan SHA-256 {sha256(output / "plan.json")}')


def verified_plan(root, digest, runtime=REPO):
    plan = load_schema(pinned(root, 'plan.json', digest), SCHEMA)
    record = runtime / STEP04 / 'execution/formal_record/tenv'
    base = load_schema(pinned(record, 'plan.json', plan['step04_plan_sha256']), 'step04-evaluation-plan-v1')
    base_recipes = base['folds'][0]
    recipes = plan['formal_recipes']
    if set(recipes) != set(map(str, SEEDS)) or plan['inference'] != base['inference']:
        raise ValueError('unexpected warm-start recipes or inference')
    for seed, item in recipes.items():
        pinned(root, item['path'], base_recipes['recipes'][seed]['sha256'])
        if (item['sha256'] != base_recipes['recipes'][seed]['sha256']
                or item['config_sha256'] != base_recipes['recipe_sha256_by_seed'][seed]):
            raise ValueError('warm-start formal recipe differs from Step 04')
    active = plan['active_recipes']
    if plan['stage'] == 'smoke':
        if set(active) != {'17'} or read_json(pinned(root, active['17']['path'], active['17']['sha256'])) != {
                **read_json(root / recipes['17']['path']), **SMOKE}:
            raise ValueError('unexpected bounded warm-start recipe')
    elif plan['stage'] != 'formal' or active != recipes:
        raise ValueError('unexpected warm-start stage or recipe inventory')

    source = Path(plan['source_root'])
    source_digest, manifest = release(runtime, source)
    if source_digest != plan['source_sha256']:
        raise ValueError('warm-start source release differs')
    split, expected = membership(source, manifest)
    if read_json(pinned(root, 'split.json', plan['split_sha256'])) != split:
        raise ValueError('warm-start split differs from deterministic parent–dose membership')
    request = read_json(pinned(root, 'request.json', plan['request_sha256']))
    paths, hashes = predictor_paths(0), request['file_sha256']
    if (hashes.keys() != paths.keys()
            or request != make_request(root, manifest, source_digest, plan['split_sha256'], hashes)
            or any(hashes[key] != manifest['files'][paths[key]] for key in ('genes', 'drugs', 'morgan', 'controls'))):
        raise ValueError('warm-start predictor request or file allowlist differs')
    for key, path in paths.items():
        pinned(root, path, hashes[key])
    for role, cases in expected.items():
        if rows(root / paths[f'{role}_cases']) != cases:
            raise ValueError(f'{role} cases differ from frozen parent–dose membership')
    return plan, request


def worker(runtime, root, digest, phase, seed):
    plan, request = verified_plan(root, digest, runtime)
    item = plan['active_recipes'][str(seed)]
    from formal import fit_slot, predict_slot, slot_lock, use_runtime

    use_runtime(runtime)
    import native_training
    from warm_inputs import WarmInputs

    # Use the same native fit and held-out workers with the warm-start input policy.
    native_training.PredictorInputs = WarmInputs
    recipe = pinned(root, item['path'], item['sha256'])
    config = native_training.load_recipe(runtime, recipe, request)
    if (config['stage'], config['seed']) != (plan['stage'], seed):
        raise ValueError('warm-start stage or seed mismatch')
    if 'config_sha256' in item and native_training.json_hash(config) != item['config_sha256']:
        raise ValueError('warm-start configuration differs from Step 04')
    identity = {**request['identity'], 'seed': seed, 'stage': plan['stage']}
    slot = root / 'slots' / f'seed-{seed}'
    with slot_lock(slot):
        if phase == 'fit':
            fit_slot(slot, digest, runtime, request, identity, recipe, config)
        else:
            predict_slot(slot, plan, digest, request, identity, config)


def score(runtime, root, digest, output):
    from formal import verify_fit_seal
    from training_inputs import json_hash

    plan, request = verified_plan(root, digest, runtime)
    if output.exists() or output.is_symlink():
        raise ValueError('score requires a fresh output directory')
    paths = request['input']['predictor_files']
    cases = rows(root / paths['test_cases'])
    groups, case_ids = grouped_cases(cases)
    fold = {'identity': request['identity'], 'cases': cases, 'groups': groups, 'case_ids': case_ids,
            'n_genes': N_GENES, 'recipes': {}}
    outcomes, predictions = {}, {}
    for seed_text, item in plan['active_recipes'].items():
        seed = int(seed_text)
        path = root / 'slots' / f'seed-{seed}' / 'prediction' / 'prediction.json'
        # Pin the independent completed fit as well as the files in the prediction seal.
        manifest = load_schema(path, 'step04-heldout-prediction-v1')
        fit = path.parent.parent / 'fit'
        verify_fit_seal(fit, digest)
        run = read_json(fit / 'run.json')
        if manifest['run_sha256'] != sha256(fit / 'run.json'):
            raise ValueError('warm-start prediction differs from the sealed fit')
        raw = read_json(root / item['path'])
        if any(run['recipe'][key] != value for key, value in raw.items()):
            raise ValueError('warm-start prediction used a different recipe')
        fold['recipes'][seed] = item.get('config_sha256', run['recipe_sha256'])
        slot = ('tenv', 0, seed)
        outcomes[slot] = {'seed': seed, 'status': 'completed', 'prediction': str(path.relative_to(root)),
                          'prediction_manifest_sha256': sha256(path)}
        values = prediction(root, outcomes[slot], fold, plan['stage'], digest)
        if (manifest['case_count'] != len(cases) or manifest['group_count'] != len(groups)
                or any(not np.isfinite(values[start:start + 32]).all()
                       for start in range(0, len(cases), 32))):
            raise ValueError('warm-start prediction coverage or values differ from the frozen cases')
        predictions[slot] = values
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as temp:
        scratch = Path(temp)
        evaluator = scratch / 'fold_0/evaluator'
        evaluator.mkdir(parents=True)
        # Resolve held-out treated rows only inside the evaluator and after prediction verification.
        source = Path(plan['source_root'])
        _, source_release = release(runtime, source)
        locations = expression_sources(source, source_release['files'], cases)
        grouped_truth(evaluator, groups, cases, locations)
        for key in ('train_cases', 'train_expression'):
            target = scratch / paths[key]
            try:
                os.link(root / paths[key], target)
            except OSError as error:
                if error.errno != EXDEV:
                    raise
                shutil.copyfile(root / paths[key], target)
        fold.update(root=scratch, paths=paths, files={
            str(path.relative_to(scratch)): sha256(path)
            for path in sorted(scratch.rglob('*')) if path.is_file()})
        pending = scratch / 'scores'
        pending.mkdir()
        write_results(pending, {('tenv', 0): fold}, outcomes, predictions, {
            'schema': 'step04-evaluation-v1', 'stage': plan['stage'], 'datasets': ['tenv'],
            'plan_sha256': digest, 'outcomes_sha256': json_hash(list(outcomes.values())),
            'regime': 'warm-start', 'outcomes': list(outcomes.values()),
            'evaluator_files': {name: sha for name, sha in fold['files'].items() if '/evaluator/' in name},
        }, include_pr=True)
        pending.rename(output)


def grouped_truth(output, groups, cases, locations):
    inventory = []
    # Group means and labels are persisted with the evaluator's exact axis and dtypes.
    truth, positive, negative = [
        np.lib.format.open_memmap(output / f'{name}.npy', mode='w+', dtype=dtype,
                                  shape=(len(groups), N_GENES))
        for name, dtype in (('truth', 'float32'), ('positive', 'bool'), ('negative', 'bool'))
    ]
    for index, (group, members) in enumerate(groups.items()):
        case = cases[members[0]]
        inventory.append({'group_id': group, 'drug_key': case['drug_key'], 'dose_um': case['dose_um'],
                          'n_wells': len(members)})
        total = np.zeros(N_GENES, dtype=np.float64)
        for member in members:
            array, row = locations[cases[member]['case_id']]
            observed = array[row]
            if not np.isfinite(observed).all():
                raise ValueError('nonfinite warm-start held-out expression')
            total += observed
        truth[index] = total / len(members)
        positive[index], negative[index] = truth_labels(truth[index])
    for array in (truth, positive, negative):
        array.flush()
    write_csv(output / 'groups.csv', inventory)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime', type=Path, default=REPO)
    commands = parser.add_subparsers(dest='command', required=True)
    build = commands.add_parser('prepare')
    build.add_argument('formal', type=Path)
    build.add_argument('root', type=Path)
    build.add_argument('--smoke', action='store_true')
    for command in ('verify', 'fit', 'predict', 'score'):
        phase = commands.add_parser(command)
        phase.add_argument('root', type=Path)
        phase.add_argument('--plan-sha256', required=True)
        if command in ('fit', 'predict'):
            phase.add_argument('seed', type=int, choices=SEEDS)
        if command == 'score':
            phase.add_argument('output', type=Path)
    args = parser.parse_args()
    runtime, root = args.runtime.resolve(), args.root.absolute()
    try:
        if args.command == 'prepare':
            prepare(runtime, args.formal.resolve(), root, args.smoke)
        elif args.command == 'verify':
            verified_plan(root, args.plan_sha256, runtime)
        elif args.command == 'score':
            score(runtime, root, args.plan_sha256, args.output.absolute())
        else:
            worker(runtime, root, args.plan_sha256, args.command, args.seed)
    except (ArithmeticError, KeyError, ValueError, TypeError, IndexError, OSError, EOFError) as error:
        parser.exit(1, f'warm-start tenv {args.command} rejected: {error}\n')


if __name__ == '__main__':
    main()
