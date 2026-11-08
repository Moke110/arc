#!/usr/bin/env python3
"""Build, verify and qualify the globally cold Tahoe fold 0 source release."""
import argparse
import os
from pathlib import Path
import shutil
import sys
import tempfile

sys.dont_write_bytecode = True
REPO = next(path for path in Path(__file__).resolve().parents if (path / '.git').exists())
STEP03 = Path('stage2/03_modes_tenv')
STEP04 = Path('stage2/04_pdiff_hyperparameter')
STEP05 = Path('stage2/05_pdiff_transfer')
sys.path[:0] = [str(REPO / step / 'execution/script') for step in (STEP04, STEP03)]
for name in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(name, '2')

import numpy as np

from artifacts import read_json, rows, SEEDS, sha256, write_csv, write_json
from training_inputs import predictor_paths, request_for, verify_files

PATHS = predictor_paths(0)
ROLES = ('train', 'validation', 'test')
FITTED = ROLES[:2]
EXPRESSION = {f'{role}_expression' for role in FITTED}
CASE_FILES = {f'{role}_cases' for role in ROLES}
SHARED = PATHS.keys() - EXPRESSION - CASE_FILES
FORMAL_RECORD = REPO / STEP04 / 'execution/formal_record/tahoe/plan.json'
TENV_RECORD = REPO / STEP05 / 'execution/cold_record/replay.json'
SMOKE_RECIPE = REPO / STEP04 / 'execution/tahoe_smoke_record/recipe.json'
SCHEMA = 'step05-cold-tahoe-source-v1'
N_GENES = 18164
ROW_BATCH = 128


def pin(root, name, expected):
    verify_files(root, {'artifact': name}, {'artifact': expected})
    return root / name


def frozen_inputs(runtime, formal, split):
    plan_sha = sha256(FORMAL_RECORD)
    plan = read_json(pin(formal, 'plan.json', plan_sha))
    fold = plan['folds'][0]
    source_request = read_json(pin(formal, fold['request'], fold['request_sha256']))
    source_root = Path(source_request['input_root'])
    source = source_request['input']
    if (source_request != request_for(source_root, source)
            or (source['dataset'], source['fold'], source['n_genes'], source['target_kind'])
            != ('tahoe', 0, N_GENES, 'delta') or source['predictor_files'] != PATHS):
        raise ValueError('unexpected Tahoe source request')
    source_manifest = read_json(pin(source_root, 'manifest.json', source['source_input_sha256']))
    if (source_manifest['schema'] != 'tahoe-input-v1'
            or any(source_manifest['files'][name] != source['file_sha256'][key]
                   for key, name in PATHS.items())
            or source_manifest['files']['folds.csv'] != source['split_sha256']):
        raise ValueError('Tahoe release differs from frozen source request')
    paths = {key: PATHS[key] for key in CASE_FILES | {'drugs'}}
    verify_files(source_root, {**paths, 'folds': 'folds.csv'},
                 {**source['file_sha256'], 'folds': source['split_sha256']})

    tenv_record = read_json(TENV_RECORD)
    split_sha = tenv_record['local_artifact_sha256']['cold-tenv-ticket76/split.json']
    split = read_json(pin(split.parent, split.name, split_sha))
    if (split['schema'] != 'step05-cold-tenv-split-v1'
            or split['identity'] != tenv_record['split_identity']
            or split['identity']['gene_sha256'] != source['gene_sha256']):
        raise ValueError('unexpected cold tenv split identity')
    target = {role: set(split['roles'][role]['parent_connectivity']) for role in ROLES}
    fixed = STEP03 / 'execution/models/perturbdiff/v2.config.json'
    fixed_sha, = (value for path, value in plan['launch_files'].items() if path.endswith(str(fixed)))
    pin(runtime, str(fixed), fixed_sha)
    recipes = fold['recipes']
    if set(recipes) != set(map(str, SEEDS)):
        raise ValueError('Tahoe fold 0 does not contain the three formal recipes')
    recipe_paths = {seed: pin(formal, entry['path'], entry['sha256']) for seed, entry in recipes.items()}
    provenance = {'tahoe_plan_sha256': plan_sha, 'tahoe_request_sha256': fold['request_sha256'],
                  'tahoe_input_identity': source_request['identity'], 'tenv_split_sha256': split_sha,
                  'tenv_split_identity': split['identity'], 'fixed_config_sha256': fixed_sha,
                  'recipes': {seed: {'path': f'recipes/seed-{seed}.json', 'sha256': entry['sha256'],
                                     'resolved_sha256': fold['recipe_sha256_by_seed'][seed]}
                              for seed, entry in recipes.items()}}
    return source_root, source, target, recipe_paths, provenance


def parent_set(cases, parents):
    return {parents[case['drug_key']] for case in cases}


def role_counts(cases, parents):
    return {'wells': len(cases), 'groups': len({case['group_id'] for case in cases}),
            'drugs': len({case['drug_key'] for case in cases}),
            'parent_connectivity_blocks': len(parent_set(cases, parents))}


def keep_indices(role, cases, parents, forbidden):
    return [index for index, case in enumerate(cases)
            if role == 'test' or parents[case['drug_key']] not in forbidden]


def verify_roles(root, original, parents, target):
    """Verify the filtered release against ordered source case membership."""
    forbidden = target['validation'] | target['test']
    counts = {}
    for role in ROLES:
        before = original[role]
        retained = rows(root / PATHS[f'{role}_cases'])
        retained_parents = parent_set(retained, parents)
        if role in FITTED and retained_parents & forbidden:
            raise ValueError(f'forbidden parent in Tahoe {role}')
        keep = keep_indices(role, before, parents, forbidden)
        if not retained or retained != [before[index] for index in keep]:
            raise ValueError(f'Tahoe {role} membership, order or control pairing changed')
        removed = [case for case in before if role in FITTED and parents[case['drug_key']] in forbidden]
        before_parents = parent_set(before, parents)
        counts[role] = {
            'original': role_counts(before, parents), 'retained': role_counts(retained, parents),
            'excluded': role_counts(removed, parents), 'excluded_parents': sorted(parent_set(removed, parents)),
            'original_target_overlap': {name: len(before_parents & blocks) for name, blocks in target.items()},
            'retained_target_overlap': {name: len(retained_parents & blocks) for name, blocks in target.items()},
        }
    return counts


def membership(source_root, forbidden):
    cases = {role: rows(source_root / PATHS[f'{role}_cases']) for role in ROLES}
    parents = {row['drug_key']: row['parent_connectivity'] for row in rows(source_root / PATHS['drugs'])}
    fold = [row for row in rows(source_root / 'folds.csv') if row['fold'] == '0']
    retained = []
    for row in fold:
        parent = parents[row['drug_key']]
        if parent != row['parent_connectivity']:
            raise ValueError('Tahoe parent connectivity differs between drugs and folds')
        if row['role'] == 'test' or parent not in forbidden:
            retained.append(row)
    return cases, parents, retained


def write_expression(source, destination, expected_rows, keep):
    # Stream the destination so only the source mapping and a small row batch
    # contribute to process RSS. The source column order and float32 bytes survive.
    matrix = np.load(source, mmap_mode='r', allow_pickle=False)
    if matrix.shape != (expected_rows, N_GENES) or matrix.dtype != np.dtype('float32'):
        raise ValueError('Tahoe source expression shape or dtype mismatch')
    with destination.open('wb') as stream:
        np.lib.format.write_array_header_2_0(stream, {
            'descr': np.lib.format.dtype_to_descr(matrix.dtype),
            'fortran_order': False, 'shape': (len(keep), N_GENES),
        })
        for start in range(0, len(keep), ROW_BATCH):
            matrix[keep[start:start + ROW_BATCH]].tofile(stream)


def manifest_for(provenance, forbidden, counts, hashes, split_sha):
    return {'schema': SCHEMA, **provenance, 'forbidden_parents': sorted(forbidden),
            'roles': counts, 'file_sha256': hashes, 'folds_sha256': split_sha}


def filtered_request(root, source, manifest, manifest_sha):
    consumer = {**source, 'source_input_sha256': manifest_sha,
                'split_sha256': manifest['folds_sha256'], 'file_sha256': manifest['file_sha256']}
    return request_for(root, consumer)


def build_release(runtime, formal, split, output):
    source_root, source, target, recipes, provenance = frozen_inputs(runtime, formal, split)
    forbidden = target['validation'] | target['test']
    original, parents, fold = membership(source_root, forbidden)
    verify_files(source_root, {key: PATHS[key] for key in SHARED | EXPRESSION}, source['file_sha256'])
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(output)
    with tempfile.TemporaryDirectory(dir=output.parent, prefix='.cold-tahoe-') as pending:
        root = Path(pending)
        (root / 'fold_0').mkdir()
        (root / 'recipes').mkdir()
        for key in SHARED | {'test_cases'}:
            shutil.copyfile(source_root / PATHS[key], root / PATHS[key])
        for role in FITTED:
            cases = original[role]
            keep = keep_indices(role, cases, parents, forbidden)
            write_csv(root / PATHS[f'{role}_cases'], [cases[index] for index in keep], tuple(cases[0]))
            path = PATHS[f'{role}_expression']
            write_expression(source_root / path, root / path, len(cases), keep)
        write_csv(root / 'folds.csv', fold)
        for seed, path in recipes.items():
            shutil.copyfile(path, root / provenance['recipes'][seed]['path'])
        counts = verify_roles(root, original, parents, target)
        hashes = {key: sha256(root / name) for key, name in PATHS.items()}
        manifest = manifest_for(provenance, forbidden, counts, hashes, sha256(root / 'folds.csv'))
        write_json(root / 'release.json', manifest)
        request = filtered_request(output, source, manifest, sha256(root / 'release.json'))
        write_json(root / 'request.json', request)
        root.rename(output)


def verify_release(runtime, formal, split, root):
    source_root, source, target, _, provenance = frozen_inputs(runtime, formal, split)
    forbidden = target['validation'] | target['test']
    original, parents, fold = membership(source_root, forbidden)
    counts = verify_roles(root, original, parents, target)
    if rows(root / 'folds.csv') != fold:
        raise ValueError('Tahoe filtered fold membership differs from the frozen rule')
    manifest = read_json(root / 'release.json')
    hashes = manifest['file_sha256']
    if (hashes.keys() != PATHS.keys() or any(hashes[key] != source['file_sha256'][key] for key in SHARED)
            or manifest != manifest_for(provenance, forbidden, counts, hashes, sha256(root / 'folds.csv'))):
        raise ValueError('Tahoe filtered release manifest differs from the frozen rule')
    verify_files(root, PATHS, hashes)
    for entry in provenance['recipes'].values():
        pin(root, entry['path'], entry['sha256'])
    request = filtered_request(root, source, manifest, sha256(root / 'release.json'))
    if request != read_json(root / 'request.json'):
        raise ValueError('Tahoe filtered request differs from the verified release')
    return request


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime-root', type=Path, default=REPO)
    parser.add_argument('--source-formal', type=Path)
    parser.add_argument('--tenv-split', type=Path)
    commands = parser.add_subparsers(dest='command', required=True)
    for name in ('build', 'verify', 'smoke-fit', 'formal-fit', 'reload'):
        command = commands.add_parser(name)
        command.add_argument('release', type=Path)
        if name in ('smoke-fit', 'formal-fit', 'reload'):
            command.add_argument('run', type=Path)
        if name == 'formal-fit':
            command.add_argument('--seed', type=int, choices=SEEDS, required=True)
        if name == 'reload':
            command.add_argument('replay', type=Path)
    args = parser.parse_args()
    runtime = args.runtime_root.resolve()
    sys.path.insert(0, str(runtime / STEP03 / 'execution/script'))
    formal = (args.source_formal or runtime / STEP04 / 'execution/output/tahoe-formal-ticket73').resolve()
    split = (args.tenv_split or runtime / STEP05 / 'execution/output/cold-tenv-ticket76/split.json').absolute()
    root = args.release.resolve()
    if args.command == 'build':
        build_release(runtime, formal, split, root)
        return
    request = verify_release(runtime, formal, split, root)
    if args.command == 'verify':
        return
    from native_training import fit, reload
    if args.command == 'reload':
        reload(request, args.run, args.replay)
    else:
        recipe = SMOKE_RECIPE if args.command == 'smoke-fit' else root / f'recipes/seed-{args.seed}.json'
        fit(runtime, request, recipe, args.run)


if __name__ == '__main__':
    main()
