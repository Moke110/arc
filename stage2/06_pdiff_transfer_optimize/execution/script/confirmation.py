#!/usr/bin/env python3
"""Freeze the fold 1 confirmation decision and release its shared Tahoe EMAs."""
import argparse
import mmap
import os
from pathlib import Path
import shutil
import tempfile

import baseline
import cold_tahoe as cold
import dose_reset
import learning_rate
from artifacts import SEEDS, read_json, rows, sha256, write_csv, write_json
from cold_tenv import unique_rows
from evaluate import grouped_cases, load_array
from formal import CHECKPOINT_FIELDS, completed_fit, pinned, read_plan, slot_lock, use_runtime
from formal_matrix import bind, check_binding
from screen import finite_decimal, screen_decision
from tahoe_bridge import dose_layout
from tahoe_zero_shot import key_layout, parent_blocks, role_cases
from training_inputs import json_hash, request_for, tenv_consumer, verify_files

SCHEMA = 'step06-confirmation-plan-v1'
SOURCE_SCHEMA = 'step06-confirmation-tahoe-release-v1'
FIT_SCHEMA = 'step06-confirmation-tahoe-fit-v1'
FOLD = 1
TAHOE_FOLD = 0
FIXED = cold.STEP03 / 'execution/models/perturbdiff/v2.config.json'
PLAN_PATH = Path('confirmation/plan.json')
RELEASE_PATH = Path('confirmation/tahoe-fold-1-release')


def verified_screen(output, root, arm):
    binding = bind(root / 'decision.json')
    report = read_json(binding['path'])
    if arm == 'dose_reset':
        if ((report['schema'], report['fold'], report['arm'], report['comparison_mode'])
                != ('step06-screen-decision-v1', 0, dose_reset.ARM, 'loss-selected')):
            raise ValueError('unexpected positive-dose screen decision')
        published = report['overall_gains']
        unavailable = report['unavailable']
    else:
        if ((report['schema'], report['fold'], report['gate_mode'], report['common_steps'])
                != (learning_rate.SCHEMA, 0, 'common',
                    {str(seed): step for seed, step in learning_rate.COMMON_STEPS.items()})):
            raise ValueError('unexpected learning-rate screen decision')
        published = report['gains']
        unavailable = report.get('unavailable', {})
    if (set(unavailable) & set(report['scores'])
            or set(published) != set(report['scores']) | set(unavailable)
            or any(not isinstance(reason, str) or not reason.strip() for reason in unavailable.values())):
        raise ValueError(f'invalid terminal seed evidence for {arm}')

    def score_slot(seed):
        if arm == 'dose_reset':
            return dose_reset.run_root(output, seed)
        return learning_rate.seed_root(root, seed)

    if any((score_slot(int(seed)) / 'scores/scores.json').exists() for seed in unavailable):
        raise ValueError(f'{arm} has a score seal for an explicitly unavailable seed')
    gains = {int(seed): None for seed in unavailable}
    incomplete = dict(unavailable)
    for key, score in report['scores'].items():
        seed = int(key)
        slot = score_slot(seed)
        path = slot / 'scores/scores.json'
        if score != bind(path, score['sha256']):
            raise ValueError(f'{arm} score is outside its frozen screening slot')
        if arm == 'dose_reset':
            gains[seed] = dose_reset.score_gain(output, seed)
            if gains[seed] is None:
                incomplete[key] = 'overall paired directional mean AP@100 is undefined'
        else:
            evidence = read_json(path)
            plan = learning_rate.load_plan(slot, seed, evidence['plan_sha256'])
            if ((evidence['schema'], evidence['fold'], evidence['seed'], evidence['common_step_policy'])
                    != (learning_rate.SCHEMA, 0, seed, plan['common_step_policy'])):
                raise ValueError('learning-rate score differs from its frozen matched-step policy')
            for item in [evidence['fit'], evidence['truth'], *evidence['predictions'].values()]:
                check_binding(item)
            for name, digest in evidence['files'].items():
                pinned(path.parent, name, digest)
            gains[seed] = learning_rate.primary_gain(evidence)
            if gains[seed] is None:
                incomplete[key] = evidence['comparisons']['common'].get(
                    'reason', 'matched-step paired directional mean AP@100 is undefined')
    for seed, gain in gains.items():
        recorded = published[str(seed)]
        expected = None if gain is None else (float(gain) if arm == 'dose_reset' else finite_decimal(gain))
        actual = recorded if arm == 'dose_reset' or recorded is None else finite_decimal(recorded)
        if expected != actual:
            raise ValueError(f'{arm} published gain differs from seed {seed} evidence')
    decision = screen_decision(gains)
    if decision['next_seeds'] or any(report[key] != value for key, value in decision.items()):
        raise ValueError(f'{arm} requires a verified terminal screen decision')
    return {**decision, 'decision': binding, 'scores': report['scores'], 'incomplete_seeds': incomplete,
            'gains': {str(seed): None if gain is None else str(gain) for seed, gain in gains.items()}}


def predecessor(runtime, dataset, fold):
    stage = runtime / cold.STEP04 / 'execution'
    number = {'tenv': 72, 'tahoe': 73}[dataset]
    root = stage / f'output/{dataset}-formal-ticket{number}'
    record = bind(stage / f'formal_record/{dataset}/plan.json')
    path = root / 'plan.json'
    plan = read_plan(path, record['sha256'])
    if plan['dataset'] != dataset:
        raise ValueError(f'unexpected {dataset} predecessor plan')
    entry = plan['folds'][fold]
    request_path = pinned(root, entry['request'], entry['request_sha256'])
    request = read_json(request_path)
    consumer = request['input']
    source = Path(request['input_root'])
    if (source != source.resolve() or request != request_for(source, consumer)
            or (consumer['dataset'], consumer['fold'], consumer['target_kind']) != (dataset, fold, 'delta')):
        raise ValueError(f'unexpected {dataset} fold {fold} request')
    fixed_sha, = (value for name, value in plan['launch_files'].items() if name.endswith(str(FIXED)))
    fixed = bind(runtime / FIXED, fixed_sha)
    recipes = {}
    for seed in map(str, SEEDS):
        item = entry['recipes'][seed]
        recipe = bind(pinned(root, item['path'], item['sha256']))
        recipes[seed] = {**recipe, 'resolved_sha256': entry['recipe_sha256_by_seed'][seed]}
    return {'producer_plan': record, 'producer_request': bind(request_path), 'fixed_config': fixed,
            'request': request, 'recipes': recipes}


def recipient(runtime):
    target = predecessor(runtime, 'tenv', FOLD)
    request = target['request']
    source = Path(request['input_root'])
    if request['input'] != tenv_consumer(runtime, source, FOLD):
        raise ValueError('recipient differs from the frozen tenv fold 1 request')
    inference = baseline.inference_request(request, role='test')
    baseline.verify_inference(inference)
    manifest = read_json(source / 'manifest.json')
    files = manifest['files']
    metadata = {name: bind(pinned(source, name, files[name]))
                for name in ('wells.csv', 'folds.csv', 'controls.csv')}
    controls = rows(source / 'controls.csv')
    if ([int(row['control_index']) for row in controls] != list(range(len(controls)))
            or any(int(row['n_wells']) <= 0 for row in controls)):
        raise ValueError('invalid fold 1 control inventory')
    load_array(source / 'controls.npy', (len(controls), request['n_genes']))
    genes = unique_rows(rows(source / 'genes.csv'), 'symbol', 'gene')
    if len(genes) != request['n_genes']:
        raise ValueError('unexpected recipient gene axis length')
    drugs = unique_rows(rows(source / 'drugs.csv'), 'drug_key', 'tenv drug')
    memberships = unique_rows([row for row in rows(source / 'folds.csv') if row['fold'] == str(FOLD)],
                              'drug_key', 'tenv fold 1 membership')
    if (memberships.keys() != drugs.keys()
            or any(row['role'] not in cold.ROLES or not row['parent_connectivity']
                   or row['parent_connectivity'] != drugs[key]['parent_connectivity']
                   for key, row in memberships.items())):
        raise ValueError('fold 1 memberships differ from tenv drug parent connectivity')
    wells = rows(source / 'wells.csv')
    unique_rows(wells, 'case_id', 'tenv well')
    roles, assigned = {}, set()
    for role, cases in role_cases(request).items():
        expected = [well for well in wells if memberships[well['drug_key']]['role'] == role]
        parents = {row['parent_connectivity'] for row in memberships.values() if row['role'] == role}
        if (cases != expected or assigned & parents
                or any(not 0 <= int(case['control_index']) < len(controls) for case in cases)):
            raise ValueError(f'tenv fold 1 {role} case, control or parent membership mismatch')
        assigned.update(parents)
        groups, _ = grouped_cases(cases)
        roles[role] = {'parent_connectivity': sorted(parents), 'case_count': len(cases),
                       'group_count': len(groups), 'group_order_sha256': json_hash(list(groups)),
                       'case_order_sha256': json_hash(baseline.case_order(cases))}
    return {**target, 'inference': inference, 'metadata': metadata, 'roles': roles}


def source_predecessor(runtime, target):
    producer = predecessor(runtime, 'tahoe', TAHOE_FOLD)
    request = producer['request']
    source = request['input']
    root = Path(request['input_root'])
    manifest = read_json(pinned(root, 'manifest.json', source['source_input_sha256']))
    if (source['predictor_files'] != cold.PATHS or source['n_genes'] != target['request']['n_genes']
            or source['gene_sha256'] != target['request']['identity']['gene_sha256']
            or manifest['schema'] != 'tahoe-input-v1'
            or any(manifest['files'][name] != source['file_sha256'][key] for key, name in cold.PATHS.items())
            or manifest['files']['folds.csv'] != source['split_sha256']):
        raise ValueError('Tahoe predecessor differs from the pinned source or recipient gene axis')
    return producer


def plan_for(runtime, output, learning_rate_root):
    screens = {'dose_reset': verified_screen(output, output / 'dose_reset', 'dose_reset'),
               'learning_rate': verified_screen(output, learning_rate_root, 'learning_rate')}
    promoted = [arm for arm, screen in screens.items() if screen['status'] == 'promoted']
    plan = {'schema': SCHEMA, 'fold': FOLD, 'status': 'ready' if promoted else 'skipped',
            'seeds': list(SEEDS) if promoted else [], 'promoted_arms': promoted, 'screens': screens,
            'runtime': str(runtime), 'output': str(output)}
    if not promoted:
        return {**plan, 'reason': 'neither fold 0 intervention qualified for confirmation'}
    target = recipient(runtime)
    source = source_predecessor(runtime, target)
    steps = {}
    if 'learning_rate' in promoted:
        for seed, step in learning_rate.COMMON_STEPS.items():
            config = read_json(target['recipes'][str(seed)]['path'])
            if step > config['max_steps'] or (step % config['validation_interval'] and step != config['max_steps']):
                raise ValueError(f'fold 1 common step is outside the native cadence for seed {seed}')
            steps[str(seed)] = step
    forbidden = set(target['roles']['validation']['parent_connectivity'])
    forbidden.update(target['roles']['test']['parent_connectivity'])
    return {**plan, 'recipient': target, 'source_predecessor': source,
            'source_policy': {'base_fold': TAHOE_FOLD, 'target_fold': FOLD, 'fit_roles': list(cold.FITTED),
                              'forbidden_parents': sorted(forbidden)},
            'comparison': {'dose_reset': 'loss-selected' if 'dose_reset' in promoted else None,
                           'common_steps': steps, 'learning_rate': 'equal-step' if steps else None,
                           'prediction': {'seed_offset': 200000, 'max_groups': 8, 'max_profiles': 24},
                           'score': 'paired directional mean AP@100 over ordered test groups'}}


def freeze(runtime, output, learning_rate_root):
    plan = plan_for(runtime, output, learning_rate_root)
    verify_no_run(output, plan)
    return learning_rate.publish(output / PLAN_PATH, plan)


def verify_no_run(output, plan):
    if plan['status'] == 'skipped' and (
            (output / RELEASE_PATH).exists() or any(source_path(output, seed).exists() for seed in SEEDS)):
        raise ValueError('skipped confirmation contains Tahoe preparation or fit artifacts')


def load_plan(runtime, output, digest):
    plan = read_json(pinned(output, PLAN_PATH, digest))
    promoted = [arm for arm, screen in plan['screens'].items() if screen['status'] == 'promoted']
    if (plan['schema'] != SCHEMA or plan['fold'] != FOLD or plan['runtime'] != str(runtime)
            or plan['output'] != str(output) or plan['promoted_arms'] != promoted
            or plan['status'] != ('ready' if promoted else 'skipped')
            or plan['seeds'] != (list(SEEDS) if promoted else [])
            or set(plan['screens']) != {'dose_reset', 'learning_rate'}):
        raise ValueError('unexpected fold 1 confirmation plan')
    for screen in plan['screens'].values():
        check_binding(screen['decision'])
        gains = {int(seed): gain for seed, gain in screen['gains'].items()}
        decision = screen_decision(gains)
        if decision['next_seeds'] or any(screen[key] != value for key, value in decision.items()):
            raise ValueError('frozen screen decision differs from its paired effects')
        for binding in screen['scores'].values():
            check_binding(binding)
    verify_no_run(output, plan)
    return plan


def membership(plan):
    if plan['status'] != 'ready':
        raise ValueError('Tahoe preparation and fitting require a promoted intervention')
    request = plan['source_predecessor']['request']
    root = Path(request['input_root'])
    source = request['input']
    for key in ('producer_plan', 'producer_request', 'fixed_config'):
        check_binding(plan['source_predecessor'][key])
    pinned(root, 'manifest.json', source['source_input_sha256'])
    verify_files(root, {**{key: cold.PATHS[key] for key in cold.CASE_FILES | {'drugs'}}, 'folds': 'folds.csv'},
                 {**source['file_sha256'], 'folds': source['split_sha256']})
    target = {role: set(record['parent_connectivity']) for role, record in plan['recipient']['roles'].items()}
    forbidden = target['validation'] | target['test']
    if plan['source_policy'] != {'base_fold': TAHOE_FOLD, 'target_fold': FOLD, 'fit_roles': list(cold.FITTED),
                                'forbidden_parents': sorted(forbidden)}:
        raise ValueError('Tahoe exclusion policy differs from the frozen confirmation recipient')
    cases, parents, fold = cold.membership(root, forbidden)
    return root, source, target, forbidden, cases, parents, fold


def release_manifest(root, plan, digest, source, target, forbidden, cases, parents, fold,
                     *, expression_hashes=None):
    roles = cold.verify_roles(root, cases, parents, target)
    if rows(root / 'folds.csv') != fold:
        raise ValueError('Tahoe fold membership differs from the frozen exclusion')
    expression_hashes = expression_hashes or {}
    files = {key: expression_hashes[key] if key in expression_hashes else sha256(root / path)
             for key, path in cold.PATHS.items()}
    if any(files[key] != source['file_sha256'][key] for key in cold.SHARED):
        raise ValueError('Tahoe controls, Morgan fingerprints or gene order changed')
    return {'schema': SOURCE_SCHEMA, 'plan_sha256': digest, 'tahoe_fold': TAHOE_FOLD, 'tenv_fold': FOLD,
            'recipient_identity': plan['recipient']['request']['identity'],
            'tahoe_predecessor_identity': plan['source_predecessor']['request']['identity'],
            'roles': roles, 'forbidden_parents': sorted(forbidden),
            'file_sha256': files, 'folds_sha256': sha256(root / 'folds.csv')}


def prepare_release(output, plan, digest):
    original, source, target, forbidden, cases, parents, fold = membership(plan)
    root = output / RELEASE_PATH
    with slot_lock(root.parent / '.tahoe-release-lock'):
        if root.exists():
            return verify_release(output, plan, digest)
        verify_files(original, cold.PATHS, source['file_sha256'])
        with tempfile.TemporaryDirectory(dir=root.parent, prefix='.tahoe-release-') as temp:
            staging = Path(temp)
            (staging / f'fold_{TAHOE_FOLD}').mkdir()
            for key in cold.SHARED | {'test_cases'}:
                shutil.copyfile(original / cold.PATHS[key], staging / cold.PATHS[key])
            for role in cold.FITTED:
                keep = cold.keep_indices(role, cases[role], parents, forbidden)
                write_csv(staging / cold.PATHS[f'{role}_cases'], [cases[role][index] for index in keep],
                          tuple(cases[role][0]))
                path = cold.PATHS[f'{role}_expression']
                cold.write_expression(original / path, staging / path, len(cases[role]), keep)
            write_csv(staging / 'folds.csv', fold)
            manifest = release_manifest(staging, plan, digest, source, target, forbidden, cases, parents, fold)
            write_json(staging / 'release.json', manifest)
            request = cold.filtered_request(root, source, manifest, sha256(staging / 'release.json'))
            write_json(staging / 'request.json', request)
            staging.rename(root)
    return request


def verify_release(output, plan, digest, *, check_expression=True):
    original, source, target, forbidden, cases, parents, fold = membership(plan)
    root = output / RELEASE_PATH
    path = root / 'release.json'
    manifest = read_json(path)
    unchecked = {} if check_expression else {key: manifest['file_sha256'][key] for key in cold.EXPRESSION}
    expected = release_manifest(root, plan, digest, source, target, forbidden, cases, parents, fold,
                                expression_hashes=unchecked)
    if manifest != expected:
        raise ValueError('Tahoe release manifest differs from the frozen plan or files')
    request = cold.filtered_request(root, source, manifest, sha256(path))
    if read_json(root / 'request.json') != request:
        raise ValueError('Tahoe request differs from its sealed release')
    if not check_expression:
        return request
    for role in cold.FITTED:
        key = f'{role}_expression'
        rel = cold.PATHS[key]
        keep = cold.keep_indices(role, cases[role], parents, forbidden)
        before = load_array(pinned(original, rel, source['file_sha256'][key]), (len(cases[role]), source['n_genes']))
        after = load_array(root / rel, (len(keep), source['n_genes']))
        for offset in range(0, len(keep), cold.ROW_BATCH):
            stop = offset + cold.ROW_BATCH
            if before[keep[offset:stop]].tobytes() != after[offset:stop].tobytes():
                raise ValueError(f'Tahoe {role} expression differs from ordered retained cases')
            before._mmap.madvise(mmap.MADV_DONTNEED)
            after._mmap.madvise(mmap.MADV_DONTNEED)
    return request


def source_path(output, seed):
    ordinal = 7 + SEEDS.index(seed)
    return output / 'checkpoints' / f'{ordinal:02}-tahoe-source-fold-1-seed-{seed}'


def source_metadata(output, plan, digest, seed, request, directory):
    producer = plan['source_predecessor']
    recipe = producer['recipes'][str(seed)]
    check_binding(recipe)
    identity = {**request['identity'], 'stage': 'formal', 'seed': seed}
    config = read_json(directory / 'run.json')['recipe']
    run = completed_fit(directory, identity, config)
    cases = role_cases(request)
    groups = list(grouped_cases(cases['validation'])[0])
    if (run['recipe_sha256'] != recipe['resolved_sha256']
            or any(config[key] != value for key, value in read_json(recipe['path']).items())
            or run['validation_batching']['group_ids'] != groups
            or run['validation_batching']['group_ids_sha256'] != json_hash(groups)):
        raise ValueError('Tahoe fit differs from the frozen recipe or validation group order')
    checkpoint = pinned(directory, f"01-ema-selected-step-{run['best_step']:05}.pt", run['checkpoint_sha256'])
    doses, indices = dose_layout(key_layout(cases))
    if doses != sorted({float(case['dose_um']) for case in cases['train']}):
        raise ValueError('Tahoe train role lacks part of the source dose vocabulary')
    source = {**{key: run[key] for key in CHECKPOINT_FIELDS}, 'path': str(source_path(output, seed) / checkpoint.name),
              'doses_um': doses, 'dose_indices': indices}
    admission = read_json(directory / 'admission.json')
    return {'schema': FIT_SCHEMA, 'plan_sha256': digest, 'source': source,
            'release_sha256': request['identity']['source_input_sha256'],
            'request_sha256': sha256(output / RELEASE_PATH / 'request.json'),
            'recipient_identity': plan['recipient']['request']['identity'],
            'parent_connectivity': parent_blocks(request, cases), 'recipe': recipe,
            'run_sha256': sha256(directory / 'run.json'), 'resources': run.get('resources'),
            'wall_seconds': run['wall_seconds'], 'admission': admission,
            'admission_sha256': sha256(directory / 'admission.json')}


def checked_fit(output, plan, digest, seed, request):
    directory = source_path(output, seed)
    expected = source_metadata(output, plan, digest, seed, request, directory)
    if read_json(directory / 'fit.json') != expected:
        raise ValueError('Tahoe fit seal differs from its frozen lineage')
    return directory, expected


def verify_admission(path, output, digest, seed):
    return training_admission(path, output, digest, seed, script=__file__, memory_gib=18,
                              scope='Tahoe source')


def training_admission(path, output, digest, seed, *, script, memory_gib, scope, extra_flags=None):
    if path is None:
        raise ValueError(f'{scope} fitting requires scheduler admission evidence')
    preflight = read_json(path)
    if not preflight['sample']:
        raise ValueError('scheduler preflight has no resource sample')
    expected_flags = {'--seed': str(seed), '--plan-sha256': digest, '--admission': str(path),
                      **(extra_flags or {})}
    for job in preflight['jobs']:
        command = job['command']
        if (job['kind'] == 'training' and job['memory_gib'] == memory_gib
                and preflight['decisions'][job['id']]['admit'] is True
                and len(command) > 3 and Path(command[1]).resolve() == Path(script).resolve()
                and command[2:4] == ['fit', str(output)]
                and all(flag in command and command[command.index(flag) + 1:command.index(flag) + 2] == [value]
                        for flag, value in expected_flags.items())):
            return preflight
    raise ValueError(f'no admitted {memory_gib} GiB training job matches this frozen {scope}')


def fit_source(runtime, output, plan, digest, seed, admission):
    if seed not in plan['seeds']:
        raise ValueError('no fold 1 Tahoe source was authorized for this seed')
    request = verify_release(output, plan, digest)
    directory = source_path(output, seed)
    with slot_lock(directory.parent / f'.{directory.name}-lock'):
        if directory.exists():
            return checked_fit(output, plan, digest, seed, request)[1]
        preflight = verify_admission(admission, output, digest, seed)
        use_runtime(runtime)
        from native_training import fit, load_recipe

        item = plan['source_predecessor']['recipes'][str(seed)]
        recipe = Path(check_binding(item)['path'])
        config = load_recipe(runtime, recipe, request)
        if json_hash(config) != item['resolved_sha256']:
            raise ValueError('resolved Tahoe recipe changed')
        attempt = Path(tempfile.mkdtemp(dir=directory.parent, prefix=f'.{directory.name}-attempt-'))
        staging = attempt / 'run'
        print(f'Tahoe fit attempt evidence: {staging}', flush=True)
        fit(runtime, request, recipe, staging)
        run = completed_fit(staging, {**request['identity'], 'stage': 'formal', 'seed': seed}, config)
        (staging / 'best.pt').rename(staging / f"01-ema-selected-step-{run['best_step']:05}.pt")
        write_json(staging / 'admission.json', preflight)
        seal = source_metadata(output, plan, digest, seed, request, staging)
        write_json(staging / 'fit.json', seal)
        staging.rename(directory)
        attempt.rmdir()
        return seal


def verify_replay(directory, fit):
    run = read_json(directory / 'run.json')
    keys = ('identity', 'best_step', 'checkpoint_sha256')
    expected = {key: run[key] for key in keys} | run['reload']
    path = directory / 'replay.json'
    if read_json(path) != expected:
        raise ValueError('independent Tahoe replay differs from the selected native EMA')
    return {'fit_sha256': sha256(directory / 'fit.json'), 'replay_sha256': sha256(path),
            'source': fit['source']}


def reload_source(runtime, output, plan, digest, seed):
    if seed not in plan['seeds']:
        raise ValueError('no fold 1 Tahoe source was authorized for this seed')
    request = verify_release(output, plan, digest)
    directory, fit = checked_fit(output, plan, digest, seed, request)
    with slot_lock(directory.parent / f'.{directory.name}-lock'):
        seal = directory / 'source.json'
        if not seal.exists():
            use_runtime(runtime)
            from native_training import reload

            with tempfile.TemporaryDirectory(dir=directory, prefix='.replay-') as temp:
                staging = Path(temp)
                os.link(directory / 'run.json', staging / 'run.json')
                os.link(fit['source']['path'], staging / 'best.pt')
                reload(request, staging, staging / 'replay.json')
                (staging / 'replay.json').rename(directory / 'replay.json')
            learning_rate.publish(seal, verify_replay(directory, fit))
        return load_source(runtime, output, digest, seed)


def load_source(runtime, output, digest, seed):
    """Return the independently replayed EMA binding accepted by the Tahoe bridge."""
    plan = load_plan(runtime, output, digest)
    if seed not in plan['seeds']:
        raise ValueError('no fold 1 Tahoe source was authorized for this seed')
    request = verify_release(output, plan, digest, check_expression=False)
    directory, fit = checked_fit(output, plan, digest, seed, request)
    replay = verify_replay(directory, fit)
    if read_json(directory / 'source.json') != replay:
        raise ValueError('Tahoe source seal differs from its independent replay')
    return replay['source']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('freeze', 'verify', 'prepare', 'fit', 'reload', 'verify-source'))
    parser.add_argument('output', type=Path)
    parser.add_argument('--runtime', type=Path, default=baseline.REPO)
    parser.add_argument('--learning-rate-root', type=Path)
    parser.add_argument('--plan-sha256')
    parser.add_argument('--seed', type=int, choices=SEEDS)
    parser.add_argument('--admission', type=Path)
    args = parser.parse_args()
    runtime, output = args.runtime.resolve(), args.output.resolve()
    if args.phase == 'freeze':
        learning_rate_root = (args.learning_rate_root or output / 'learning-rate').resolve()
        print(freeze(runtime, output, learning_rate_root)['sha256'], flush=True)
        return
    if not args.plan_sha256:
        parser.error(f'{args.phase} requires --plan-sha256')
    digest = args.plan_sha256
    plan = load_plan(runtime, output, digest)
    if args.phase == 'verify':
        print(plan['status'], flush=True)
    elif args.phase == 'prepare':
        request = prepare_release(output, plan, digest)
        print(request['identity']['source_input_sha256'], flush=True)
    else:
        seed = args.seed
        if seed not in plan['seeds']:
            parser.error('this phase requires an authorized confirmation seed')
        if args.phase == 'fit':
            result = fit_source(runtime, output, plan, digest, seed, args.admission)
            print(result['source']['checkpoint_sha256'], flush=True)
        elif args.phase == 'reload':
            print(reload_source(runtime, output, plan, digest, seed)['checkpoint_sha256'], flush=True)
        else:
            print(load_source(runtime, output, digest, seed)['checkpoint_sha256'], flush=True)


if __name__ == '__main__':
    main()
