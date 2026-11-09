#!/usr/bin/env python3
"""Fit and test promoted Transfer interventions on the frozen cold-drug fold 1."""
import argparse
from pathlib import Path
import tempfile
from unittest.mock import patch

import baseline
import confirmation
import dose_reset
import learning_rate as lr
from artifacts import SEEDS, read_json, rows, write_csv, write_json
from evaluate import grouped_cases, load_array, mean_delta, np
from formal import completed_fit, pinned, slot_lock, use_runtime
from formal_matrix import bind, check_binding
from native_lr import NATIVE_CHECKPOINT, common_checkpoint, fit_diagnostic
from tahoe_bridge import load_initialization, request_digest, seal_initialization
from training_inputs import json_hash, verify_files
from training_policy import MONITOR, verify_best_checkpoint

FOLD = 1
PLAN = Path('replication/plan.json')
RESULT = Path('replication/result.json')
PLAN_SCHEMA = 'step06-replication-plan-v1'
RESULT_SCHEMA = 'step06-replication-result-v1'
SCORE_SCHEMA = 'step06-replication-scores-v1'
ORIGINAL = 'Transfer'
PRIMARY_MODES = {'dose_reset': 'native', 'learning_rate': 'common'}


def seed_slot(output, seed):
    return output / 'replication' / f'fold-1-seed-{seed}'


def checkpoint_slot(output, seed):
    return output / 'checkpoints' / f'{10 + SEEDS.index(seed):02}-fold-1-target-seed-{seed}'


def ensure_no_run(output):
    if any(seed_slot(output, seed).exists() or checkpoint_slot(output, seed).exists() for seed in SEEDS):
        raise ValueError('skipped replication contains fold 1 fit, prediction or score artifacts')


def terminal_dispositions(entries, allowed):
    if set(entries) - set(allowed):
        raise ValueError('unavailable evidence is outside the frozen replication')
    dispositions = {}
    for key, (status, reason) in entries.items():
        if status not in ('failed', 'incomplete', 'skipped') or not reason.strip():
            raise ValueError('unavailable evidence requires a terminal status and reason')
        dispositions[key] = {'status': status, 'reason': reason}
    return dispositions


def freeze(runtime, output, confirmation_sha, source_unavailable=None):
    upstream = confirmation.load_plan(runtime, output, confirmation_sha)
    sources = terminal_dispositions(source_unavailable or {}, map(str, upstream['seeds']))
    plan = {'schema': PLAN_SCHEMA, 'fold': FOLD, 'runtime': str(runtime), 'output': str(output),
            'confirmation_sha256': confirmation_sha, 'promoted_arms': upstream['promoted_arms'],
            'screens': upstream['screens'], 'status': upstream['status'], 'seeds': {}}
    if upstream['status'] == 'skipped':
        ensure_no_run(output)
        return lr.publish(output / PLAN, plan)
    target = upstream['recipient']
    request = target['request']
    cases = baseline.verify_inference(target['inference'])
    test_groups, _ = grouped_cases(cases)
    plan.update(inference=target['inference'], prediction=upstream['comparison']['prediction'],
                ordered_case_group_sha256=json_hash(baseline.case_order(cases)),
                test_group_order_sha256=json_hash(list(test_groups)),
                validation_group_order_sha256=target["roles"]["validation"]["group_order_sha256"],
                score=upstream['comparison']['score'], recipient=target['producer_request'],
                fixed_config=target['fixed_config'], truth=truth_inventory(request),
                historical_scoring=historical_scoring(runtime, target, cases))
    cm_index = plan['historical_scoring']['step04_checkpoints']
    for seed in upstream['seeds']:
        key = str(seed)
        source_path = confirmation.source_path(output, seed) / 'source.json'
        if key in sources:
            if source_path.exists():
                raise ValueError(f'Tahoe source seed {seed} already has a source seal')
            plan['seeds'][key] = sources[key]
            continue
        source = confirmation.load_source(runtime, output, confirmation_sha, seed)
        source_seal = bind(source_path)
        cm = matched_cm(runtime, target, seed, cm_index)
        common_step = upstream['comparison']['common_steps'].get(key)
        config = cm['config']
        plan['seeds'][key] = {'status': 'planned', 'reason': 'awaiting paired fold 1 target fits',
                              'source': source, 'source_seal': source_seal, 'cm': cm, 'common_step': common_step,
                              'cm_common_available': common_step == cm['metadata']['best_step'],
                              'arms': {arm: {'config': lr.arm_config(config) if arm == 'learning_rate' else config,
                                             'identity': {**cm['metadata']['identity'], 'arm': arm,
                                                          'regime': baseline.REGIME}}
                                       for arm in (ORIGINAL, *upstream['promoted_arms'])}}
    return lr.publish(output / PLAN, plan)


def load_plan(runtime, output, digest):
    plan = read_json(pinned(output, PLAN, digest))
    upstream = confirmation.load_plan(runtime, output, plan['confirmation_sha256'])
    if ((plan['schema'], plan['fold'], plan['runtime'], plan['output'])
            != (PLAN_SCHEMA, FOLD, str(runtime), str(output))
            or any(plan[key] != upstream[key] for key in ('status', 'screens', 'promoted_arms'))
            or list(plan['seeds']) != list(map(str, upstream['seeds']))):
        raise ValueError('replication plan differs from frozen confirmation')
    if plan['status'] == 'skipped':
        ensure_no_run(output)
    return plan


def finalize(output, plan, digest, unavailable):
    report = {'schema': RESULT_SCHEMA, 'fold': FOLD, 'plan_sha256': digest, 'screens': plan['screens'],
              'promoted_arms': plan['promoted_arms'], 'arms': {},
              'historical_scoring': plan.get('historical_scoring')}
    if plan['status'] == 'skipped':
        if unavailable:
            raise ValueError('skipped replication has no authorized seeds')
        ensure_no_run(output)
        report.update(status='skipped', reason='neither fold 0 intervention qualified for confirmation')
    else:
        report.update(summarize(output, plan, digest, unavailable))
    lr.publish(output / RESULT, report)
    return report


def truth_inventory(request):
    root = Path(request['input_root'])
    files = read_json(pinned(root, 'manifest.json', request['identity']['source_input_sha256']))['files']
    prefix = f'fold_{FOLD}/evaluator'
    # These values are copied from the sealed release without opening evaluator files.
    return {name: {'path': str(root / prefix / name), 'sha256': files[f'{prefix}/{name}']}
            for name in ('groups.csv', 'truth.npy', 'positive.npy', 'negative.npy')}


def historical_scoring(runtime, target, cases):
    step03 = runtime / confirmation.cold.STEP03 / 'execution/final_record'
    step04 = runtime / confirmation.cold.STEP04 / 'execution/formal_record/tenv'
    step05 = runtime / baseline.transfer.STEP05 / 'execution/formal_record'
    index = bind(step04 / 'checkpoints.json')
    release = read_json(index['path'])
    evaluation = bind(step04 / 'evaluation.json', release['evaluation_sha256'])
    manifest = read_json(evaluation['path'])
    if (manifest['tenv_test_memberships_previously_scored_in_step03'] is not True
            or manifest['plan_sha256'] != target['producer_plan']['sha256']):
        raise ValueError('historical scoring does not match the fold 1 producer')
    memberships = rows(check_binding(target['metadata']['folds.csv'])['path'])
    fold0 = {row['drug_key']: row for row in memberships if row['fold'] == '0'}
    if len(fold0) != len([row for row in memberships if row['fold'] == '0']):
        raise ValueError('duplicate fold 0 drug membership')
    parent_roles = {}
    for row in fold0.values():
        parent, role = row['parent_connectivity'], row['role']
        if role not in confirmation.cold.ROLES or parent_roles.setdefault(parent, role) != role:
            raise ValueError('inconsistent fold 0 parent membership')
    test_parents = target['roles']['test']['parent_connectivity']
    overlap = {}
    for role in confirmation.cold.ROLES:
        members = [case for case in cases if fold0[case['drug_key']]['role'] == role]
        overlap[role] = {'parent_count': sum(parent_roles[parent] == role for parent in test_parents),
                         'drug_count': len({case['drug_key'] for case in members}),
                         'well_count': len(members),
                         'group_count': len({case['group_id'] for case in members})}
    return {'step03_evaluation': bind(step03 / 'evaluation.json'), 'step04_evaluation': evaluation,
            'step04_checkpoints': index, 'step04_per_run': bind(step04 / 'per_run.csv',
                                                               manifest['files']['per_run.csv']),
            'step05_fold0_context': bind(step05 / 'paired_ap100.json'),
            'fold0_overlap': overlap, 'fold1_test_previously_scored': True,
            'interpretation': 'cross-fold replication on previously scored cold-drug memberships'}


def matched_cm(runtime, target, seed, index_binding):
    release = read_json(check_binding(index_binding)['path'])
    if (release['schema'] != 'step04-tenv-checkpoints-v1'
            or release['plan_sha256'] != target['producer_plan']['sha256']):
        raise ValueError('CM index differs from the fold 1 Step 04 plan')
    request, recipe = target['request'], target['recipes'][str(seed)]
    identity = {**request['identity'], 'seed': seed, 'stage': 'formal'}
    record, = (item for item in release['checkpoints'] if item['identity'] == identity)
    relative = f'fold-1/seed-{seed}'
    plan_entry = read_json(check_binding(target['producer_plan'])['path'])['folds'][FOLD]
    if ((record['directory'], record['recipe_sha256'], record['raw_recipe_sha256'],
         record['request_sha256'], record['policy_sha256'])
            != (relative, recipe['resolved_sha256'], recipe['sha256'],
                target['producer_request']['sha256'], plan_entry['policy_sha256'])):
        raise ValueError('CM identity or recipe differs from the frozen recipient')
    root = runtime / confirmation.cold.STEP04 / 'output/checkpoints/tenv' / relative
    bindings = {key: bind(root / name, record[field]) for key, name, field in (
        ('checkpoint', 'best.pt', 'checkpoint_sha256'), ('history', 'run.json', 'run_sha256'),
        ('replay', 'replay.json', 'replay_sha256'), ('recipe', 'recipe.json', 'raw_recipe_sha256'),
        ('request', 'request.json', 'request_sha256'))}
    run = read_json(bindings['history']['path'])
    config = run['recipe']
    selected = baseline.selection(root, {**record, **bindings}, identity, config)
    replay = read_json(bindings['replay']['path'])
    expected = {key: run[key] for key in ('identity', 'best_step', 'checkpoint_sha256')} | run['reload']
    if (replay != expected or read_json(bindings['request']['path']) != request
            or run['validation_batching']['group_ids_sha256']
            != target['roles']['validation']['group_order_sha256']):
        raise ValueError('CM replay or selection used different fold 1 inputs')
    return {**bindings, 'metadata': selected, 'config': config}


def fit_binding(output, seed, arm):
    arms = (ORIGINAL, *PRIMARY_MODES)
    return checkpoint_slot(output, seed) / f'{arms.index(arm) + 1:02}-{arm.replace("_", "-")}'


def fit_selected(output, plan, digest, seed, arm):
    slot = fit_binding(output, seed, arm)
    item = plan['seeds'][str(seed)]
    seal = read_json(slot / 'fit.json')
    run = read_json(check_binding(seal['history'])['path'])
    bridge = read_json(check_binding(seal['bridge'])['path'])
    config, identity = (item['arms'][arm][key] for key in ('config', 'identity'))
    request = read_json(check_binding(plan['recipient'])['path'])
    if (seal['plan_sha256'] != digest or seal['seed'] != seed or seal['arm'] != arm
            or seal['history']['path'] != str(slot / 'run.json')
            or seal['bridge']['path'] != str(slot / 'bridge.json')
            or bridge['source'] != item['source']
            or bridge['recipient_request_sha256'] != request_digest(request)
            or run['identity'] != identity or run['recipe'] != config or run['recipe_sha256'] != json_hash(config)
            or run['status'] != seal['status'] or run['stop_reason'] != seal['reason']
            or run['validation_batching']['group_ids_sha256'] != plan['validation_group_order_sha256']):
        raise ValueError('target fit differs from frozen source, recipient or recipe')
    initial = seal['initialization']
    if (initial['mapped_state_sha256'] != bridge['mapped_state_sha256']
            or (arm == 'dose_reset') != ('dose_reset' in initial)
            or (arm == 'dose_reset' and (initial['dose_reset']['reset_doses_um'] != dose_reset.DOSES
                or initial['dose_reset']['native_positive_rows_sha256']
                != initial['dose_reset']['restored_positive_rows_sha256']
                or initial['dose_reset']['torch_cpu_rng_preserved'] is not True))):
        raise ValueError('target initialization differs from the paired bridge policy')
    expected = fit_checkpoints(slot, slot, run, item, arm)
    if seal['checkpoints'] != expected:
        raise ValueError('target fit checkpoint inventory changed')
    return seal


def fit_checkpoints(directory, destination, run, item, arm):
    config, identity = (item['arms'][arm][key] for key in ('config', 'identity'))

    def selected(name, metadata):
        checkpoint = bind(directory / name, metadata['checkpoint_sha256'])
        return {'checkpoint': {**checkpoint, 'path': str(destination / name)}, 'metadata': metadata}

    expected = {}
    native_seal = directory / 'native-selection.json'
    if run['status'] == 'completed' or native_seal.is_file():
        if run['status'] == 'completed':
            completed_fit(directory, identity, config)
        metadata = {key: run[key] for key in baseline.CHECKPOINT_KEYS}
        if native_seal.is_file():
            seal = read_json(native_seal)
            last = run['history'][-1]
            if (any(seal[key] != value for key, value in metadata.items())
                    or seal['history_sha256'] != json_hash(run['history']) or seal['stop_step'] != last['step']
                    or seal['monitor'] != MONITOR or run['monitor'] != MONITOR
                    or seal['stop_reason'] != ('validation_early_stop' if last['early_stop'] else 'max_steps')):
                raise ValueError('native selection changed after its sealed stop')
            verify_best_checkpoint(run, arm)
        expected['native'] = selected(native_name(run['best_step']), metadata)
    diagnostic = directory / 'diagnostic.json'
    if diagnostic.is_file():
        meta = read_json(diagnostic)
        step = item['common_step']
        if (arm == 'dose_reset' or step is None or run['optimizer_steps'] < step
                or (meta['identity'], meta['recipe_sha256'], meta['best_step']) != (identity, json_hash(config), step)
                or (run.get('best_step') == step and run['model_state_sha256'] != meta['model_state_sha256'])):
            raise ValueError('diagnostic EMA differs from its frozen common step or same-step native state')
        expected['common'] = selected(common_checkpoint(step), meta)
    return expected


def native_name(step):
    return f'01-ema-selected-step-{step:05}.pt'


def admission_for(path, output, digest, seed, arm):
    return confirmation.training_admission(path, output, digest, seed, script=__file__, memory_gib=10,
                                           scope='target arm', extra_flags={'--arm': arm})


def fit(runtime, output, plan, digest, seed, arm, admission):
    item = plan['seeds'].get(str(seed), {})
    if arm not in item.get('arms', {}):
        raise ValueError('no promoted replication fit was authorized for this seed and arm')
    slot = fit_binding(output, seed, arm)
    slot.parent.mkdir(parents=True, exist_ok=True)
    with slot_lock(slot.parent / f'.{slot.name}-lock'):
        if arm in arm_dispositions(output, plan, digest, seed):
            raise ValueError(f'target arm {arm} already has a terminal disposition')
        if slot.exists():
            return fit_selected(output, plan, digest, seed, arm)
        preflight = admission_for(admission, output, digest, seed, arm)
        check_binding(item['source_seal'])
        source = confirmation.load_source(runtime, output, plan['confirmation_sha256'], seed)
        if source != item['source']:
            raise ValueError('shared Tahoe source changed after replication freeze')
        use_runtime(runtime)
        import native_training as native

        request = read_json(check_binding(plan['recipient'])['path'])
        check_binding(plan['fixed_config'])
        verify_files(Path(request['input_root']), request['input']['predictor_files'], request['file_sha256'])
        config, identity = (item['arms'][arm][key] for key in ('config', 'identity'))
        worker = {**request, 'identity': {key: value for key, value in identity.items()
                                         if key not in ('seed', 'stage')}}
        recipe = Path(check_binding(item['cm']['recipe'])['path'])
        original = item['cm']['config']
        if native.load_recipe(runtime, recipe, worker) != original:
            raise ValueError('loaded recipe differs from the frozen original target recipe')
        attempt = Path(tempfile.mkdtemp(prefix=f'.{slot.name}-attempt-', dir=slot.parent))
        staging = attempt / 'run'
        print(f'Target fit attempt evidence: {staging}', flush=True)
        initial = {}

        class InitializedModel(native.JointPerturbDiff):
            def __init__(self, inputs, active_config, directory):
                super().__init__(inputs, active_config, directory)
                bridge_sha = seal_initialization(self, directory, source, request)
                if arm == 'dose_reset':
                    bridge, initial['dose_reset'] = dose_reset.load_dose_reset(self, directory, bridge_sha, request)
                else:
                    bridge = load_initialization(self, directory, bridge_sha, request)
                if not self.get_parameter(baseline.transfer.DOSE_ROWS).requires_grad:
                    raise ValueError('mapped dose embeddings must be trainable')
                initial['mapped_state_sha256'] = bridge['mapped_state_sha256']

        def resolve(_runtime, path, current):
            if path != recipe or current != worker:
                raise ValueError('native worker changed its target recipe or request')
            return config

        step = item['common_step'] if arm != 'dose_reset' else None
        try:
            with patch.multiple(native, JointPerturbDiff=InitializedModel, load_recipe=resolve):
                if step is None:
                    native.fit(runtime, worker, recipe, staging)
                else:
                    fit_diagnostic(native, runtime, worker, recipe, staging, step)
        except InterruptedError:
            pass
        run = read_json(staging / 'run.json')
        native_file = staging / ('best.pt' if step is None else NATIVE_CHECKPOINT)
        if native_file.exists():
            native_file.rename(staging / native_name(run['best_step']))
        selected = fit_checkpoints(staging, slot, run, item, arm)
        seal = {'plan_sha256': digest, 'seed': seed, 'arm': arm, 'status': run['status'],
                'reason': run['stop_reason'], 'checkpoints': selected, 'initialization': initial,
                'admission': preflight}
        for key, name in (('history', 'run.json'), ('bridge', 'bridge.json')):
            seal[key] = {**bind(staging / name), 'path': str(slot / name)}
        write_json(staging / 'fit.json', seal)
        staging.rename(slot)
        attempt.rmdir()
        return fit_selected(output, plan, digest, seed, arm)


def arm_disposition_path(output, seed):
    return seed_slot(output, seed) / 'arm_dispositions.json'


def arm_dispositions(output, plan, digest, seed):
    path = arm_disposition_path(output, seed)
    if not path.exists():
        return {}
    record = read_json(path)
    if (record['plan_sha256'], record['seed']) != (digest, seed):
        raise ValueError('terminal target dispositions differ from the frozen replication')
    entries = {arm: (item['status'], item['reason']) for arm, item in record['arms'].items()}
    return terminal_dispositions(entries, plan['seeds'][str(seed)]['arms'])


def comparisons_for(output, plan, digest, seed, unavailable):
    item = plan['seeds'][str(seed)]
    check_binding(item['source_seal'])
    fits = {}
    for arm in item['arms']:
        if arm in unavailable:
            if (fit_binding(output, seed, arm) / 'fit.json').exists():
                raise ValueError(f'target arm {arm} already has a fit seal')
        else:
            fits[arm] = fit_selected(output, plan, digest, seed, arm)
    if len({fit['bridge']['sha256'] for fit in fits.values()}) > 1:
        raise ValueError('paired fits did not start from the same mapped Tahoe source')
    cm = item['cm']
    for key in ('history', 'replay', 'recipe', 'request', 'checkpoint'):
        check_binding(cm[key])
    choices = {'CM/native': {**cm, 'fit': cm['history']}}
    for arm, fit in fits.items():
        for mode, selected in fit['checkpoints'].items():
            choices[f'{arm}/{mode}'] = {**selected, 'config': item['arms'][arm]['config'],
                                        'fit': bind(fit_binding(output, seed, arm) / 'fit.json')}
    return choices, fits


def prediction_lineage(digest, seed, key, selected):
    return {'plan_sha256': digest, 'seed': seed, 'arm_mode': key, 'fit': selected['fit']}


def predict(runtime, output, plan, digest, seed, unavailable=None):
    dispositions = arm_dispositions(output, plan, digest, seed)
    if unavailable:
        declared = terminal_dispositions(unavailable, plan['seeds'][str(seed)]['arms'])
        if dispositions and dispositions != declared:
            raise ValueError('terminal target dispositions are already sealed')
        dispositions = declared
    choices, _ = comparisons_for(output, plan, digest, seed, dispositions)
    if unavailable:
        if (seed_slot(output, seed) / 'scores.json').exists():
            raise ValueError(f'seed {seed} already has a score seal')
        lr.publish(arm_disposition_path(output, seed),
                   {'plan_sha256': digest, 'seed': seed, 'arms': dispositions})
    use_runtime(runtime)
    for key, selected in choices.items():
        baseline.predict_checkpoint(plan['inference'], selected['config'], plan['prediction'],
                                    selected['checkpoint'], selected['metadata'],
                                    seed_slot(output, seed) / 'predictions' / key,
                                    prediction_lineage(digest, seed, key, selected))


def score_groups(cases, groups, observed, labels, predicted):
    case_groups, _ = grouped_cases(cases)
    if ([row['group_id'] for row in groups] != list(case_groups)
            or any(array.shape != (len(groups), predicted.shape[1]) for array in (observed, *labels))
            or predicted.shape[0] != len(cases)):
        raise ValueError('evaluator axes or groups differ from the frozen fold 1 prediction')
    group_rows = []
    for index, group in enumerate(groups):
        members = case_groups[group['group_id']]
        if (int(group['n_wells']) != len(members)
                or any(group[key] != cases[member][key] for member in members for key in ('drug_key', 'dose_um'))):
            raise ValueError('evaluator group membership differs from the frozen test cases')
        truth = observed[index]
        if not np.isfinite(truth).all() or any(not np.array_equal(label[index], expected)
                                              for label, expected in zip(labels, baseline.truth_labels(truth))):
            raise ValueError('invalid frozen directional truth')
        delta = mean_delta(predicted, members)
        row = {'group_id': group['group_id'], 'drug_key': group['drug_key'],
               'dose_um': float(group['dose_um']), 'n_wells': len(members)}
        for name, label, values in zip(('positive', 'negative'), labels, (delta, -delta)):
            ap, _, _ = baseline.directional_curve(label[index], values)
            row[f'{name}_ap100'] = float(ap) if np.isfinite(ap) else None
            row[f'{name}_reason'] = '' if np.isfinite(ap) else 'empty_direction'
        group_rows.append(row)
    return group_rows, baseline.summarize_groups(group_rows, metric_named_support=True)


def compare(summaries, arm, mode, cm_common_available, fits=None, unavailable=None):
    left, right = (summaries.get(f'{name}/{mode}') for name in (ORIGINAL, arm))
    cm = summaries['CM/native'] if mode == 'native' or cm_common_available else None
    cm_reason = None if cm is not None else 'matching CM EMA unavailable at the frozen optimizer step'
    if left is None or right is None:
        reasons = []
        for name, rows in ((ORIGINAL, left), (arm, right)):
            if rows is None:
                reason = (unavailable or {}).get(name) or (fits or {}).get(name)
                if reason and reason['reason']:
                    reasons.append(f"{name} ({reason['status']}): {reason['reason']}")
        return {'status': 'incomplete', 'reason': '; '.join(reasons) or
                'paired target checkpoint or prediction unavailable', 'cm_reason': cm_reason, 'contrasts': []}
    rows_by_dose = []
    for index, (original, intervention) in enumerate(zip(left, right, strict=True)):
        comparator = cm[index] if cm is not None else None
        stratum = (original['dose_um'], original['group_count'])
        if any((row['dose_um'], row['group_count']) != stratum
               for row in (intervention, comparator) if row is not None):
            raise ValueError('paired comparison did not score identical dose groups')
        rows_by_dose.append({**lr.contrast(comparator, original, intervention), 'group_count': stratum[1]})
    complete = rows_by_dose[0]['mean_ap100']['gain'] is not None
    return {'status': 'completed' if complete else 'incomplete',
            'reason': None if complete else 'pooled paired directional mean AP@100 undefined',
            'cm_reason': cm_reason, 'contrasts': rows_by_dose}


def score(output, plan, digest, seed):
    item = plan['seeds'][str(seed)]
    unavailable = arm_dispositions(output, plan, digest, seed)
    choices, fits = comparisons_for(output, plan, digest, seed, unavailable)
    slot = seed_slot(output, seed)
    values, evidence = {}, {}
    for key, selected in choices.items():
        path = slot / 'predictions' / key
        values[key] = baseline.verify_prediction(plan['inference'], selected['checkpoint'], selected['metadata'],
                                                 path, prediction_lineage(digest, seed, key, selected))
        evidence[key] = bind(path / 'prediction.json')
    cases = baseline.verify_inference(plan['inference'])
    groups, _ = grouped_cases(cases)
    if (json_hash(baseline.case_order(cases)) != plan['ordered_case_group_sha256']
            or json_hash(list(groups)) != plan['test_group_order_sha256']):
        raise ValueError('test cases differ from the frozen ordering')
    # Test targets are opened exclusively after all available predictions verify.
    truth = {key: check_binding(value) for key, value in plan['truth'].items()}
    groups = rows(truth['groups.csv']['path'])
    shape = (len(groups), plan['inference']['n_genes'])
    observed = load_array(Path(truth['truth.npy']['path']), shape)
    labels = [load_array(Path(truth[f'{name}.npy']['path']), shape, dtype='bool')
              for name in ('positive', 'negative')]
    results = {key: score_groups(cases, groups, observed, labels, predicted)
               for key, predicted in values.items()}
    summaries = {key: value[1] for key, value in results.items()}
    comparisons = {arm: {mode: compare(summaries, arm, mode, item['cm_common_available'], fits, unavailable)
                         for mode in (('native', 'common') if arm == 'learning_rate' else ('native',))}
                   for arm in plan['promoted_arms']}
    report = {'schema': SCORE_SCHEMA, 'fold': FOLD, 'seed': seed, 'plan_sha256': digest,
              'common_step': item['common_step'], 'predictions': evidence, 'truth': truth,
              'fits': {arm: {'binding': bind(fit_binding(output, seed, arm) / 'fit.json'),
                             'status': fit['status'], 'reason': fit['reason'],
                             'steps': {mode: selected['metadata']['best_step']
                                       for mode, selected in fit['checkpoints'].items()}}
                       for arm, fit in fits.items()},
              'cm_native_step': item['cm']['metadata']['best_step'], 'summaries': summaries,
              'comparisons': comparisons}
    if unavailable:
        report['arm_dispositions'] = bind(arm_disposition_path(output, seed))
    with tempfile.TemporaryDirectory(dir=slot, prefix='.scores-') as temp:
        pending = Path(temp) / 'per_group.csv'
        rows_out = [{'arm_mode': key, **row} for key, value in results.items() for row in value[0]]
        write_csv(pending, rows_out)
        with slot_lock(slot / '.score-lock'):
            path = slot / 'per_group.csv'
            if path.exists():
                bind(path, bind(pending)['sha256'])
            else:
                pending.rename(path)
            report['per_group'] = bind(path)
            lr.publish(slot / 'scores.json', report)
    return report


def aggregate(per_seed):
    paired = [(seed, row['contrasts'][0]) for seed, row in per_seed.items() if row['contrasts']]
    metrics = {}
    for metric in lr.METRICS:
        valid = [(seed, row[metric]) for seed, row in paired if row[metric]['gain'] is not None]
        matched = [(seed, row) for seed, row in valid if row['cm_gap'] is not None]
        gain = float(np.mean([row['gain'] for _, row in valid])) if valid else None
        average = {key: float(np.mean([row[key] for _, row in matched])) if matched else None
                   for key in ('cm_gap', 'gain', 'residual_gap')}
        gap = average['cm_gap']
        metrics[metric] = {'gain': gain, 'paired_seeds': [int(seed) for seed, _ in valid],
                           'matched_cm_seeds': [int(seed) for seed, _ in matched],
                           'cm_gap': gap, 'gain_on_cm_matched_seeds': average['gain'],
                           'residual_gap': average['residual_gap'],
                           'recovery': average['gain'] / gap if gap is not None and gap > 0 else None,
                           'planned_seed_count': len(per_seed)}
    return metrics


def summarize(output, plan, digest, unavailable):
    terminal = terminal_dispositions(unavailable, plan['seeds'])
    score_evidence, scores, dispositions = {}, {}, {}
    for seed, item in plan['seeds'].items():
        path = seed_slot(output, int(seed)) / 'scores.json'
        disposition = terminal.get(seed)
        if item['status'] != 'planned':
            frozen = {key: item[key] for key in ('status', 'reason')}
            if disposition is not None and disposition != frozen:
                raise ValueError(f'seed {seed} disposition differs from the frozen Tahoe source disposition')
            disposition = frozen
        if disposition is not None:
            if path.exists():
                raise ValueError(f'seed {seed} already has a score seal')
            dispositions[seed] = disposition
            continue
        if not path.exists():
            raise ValueError(f'seed {seed} is still planned and requires scores or a terminal disposition')
        binding = bind(path)
        record = read_json(binding['path'])
        if ((record['schema'], record['fold'], record['seed'], record['plan_sha256'])
                != (SCORE_SCHEMA, FOLD, int(seed), digest)):
            raise ValueError('seed score differs from the frozen replication plan')
        for evidence in (record['per_group'], *record['truth'].values(), *record['predictions'].values(),
                         *(fit['binding'] for fit in record['fits'].values())):
            check_binding(evidence)
        if 'arm_dispositions' in record:
            expected = bind(arm_disposition_path(output, int(seed)))
            if record['arm_dispositions'] != expected:
                raise ValueError('terminal target disposition seal changed')
        score_evidence[seed], scores[seed] = binding, record
        missing = [arm for arm in plan['promoted_arms']
                   if record['comparisons'][arm][PRIMARY_MODES[arm]]['status'] != 'completed']
        dispositions[seed] = {'status': 'incomplete' if missing else 'completed',
                              'reason': f"missing primary comparison for {', '.join(missing)}" if missing else None}
    arms = {}
    for arm in plan['promoted_arms']:
        primary = PRIMARY_MODES[arm]
        by_mode = {}
        for mode in (('native', 'common') if primary == 'common' else ('native',)):
            rows_by_seed = {seed: (scores[seed]['comparisons'][arm][mode] if seed in scores
                                   else {**dispositions[seed], 'cm_reason': None, 'contrasts': []})
                            for seed in plan['seeds']}
            by_mode[mode] = {'per_seed': rows_by_seed, 'pooled': aggregate(rows_by_seed)}
        metric = by_mode[primary]['pooled']['mean_ap100']
        complete = metric['paired_seeds'] == list(SEEDS)
        arms[arm] = {'primary_mode': primary, 'comparisons': by_mode,
                     'status': 'completed' if complete else 'incomplete',
                     'reason': None if complete else 'one or more frozen seeds lack a paired primary score',
                     'replicated': metric['gain'] > 0 if complete else None}
    states = {row['status'] for row in dispositions.values()}
    complete = all(arm['status'] == 'completed' for arm in arms.values())
    status = 'completed' if complete else 'failed' if states == {'failed'} else 'incomplete'
    return {'status': status, 'reason': None if complete else 'one or more primary comparisons are unavailable',
            'arms': arms, 'scores': score_evidence, 'seed_dispositions': dispositions}


def cli_dispositions(parser, entries):
    unavailable = {key: (status, reason) for key, status, reason in entries}
    if len(unavailable) != len(entries):
        parser.error('duplicate unavailable entry')
    return unavailable


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('freeze', 'status', 'fit', 'predict', 'score', 'finalize'))
    parser.add_argument('output', type=Path)
    parser.add_argument('--runtime', type=Path, default=baseline.REPO)
    parser.add_argument('--confirmation-sha256')
    parser.add_argument('--plan-sha256')
    parser.add_argument('--seed', type=int, choices=SEEDS)
    parser.add_argument('--arm', choices=(ORIGINAL, *PRIMARY_MODES))
    parser.add_argument('--admission', type=Path)
    parser.add_argument('--unavailable', nargs=3, metavar=('SEED', 'STATUS', 'REASON'),
                        action='append', default=[])
    parser.add_argument('--source-unavailable', nargs=3, metavar=('SEED', 'STATUS', 'REASON'),
                        action='append', default=[])
    parser.add_argument('--arm-unavailable', nargs=3, metavar=('ARM', 'STATUS', 'REASON'),
                        action='append', default=[])
    args = parser.parse_args()
    if args.arm_unavailable and args.phase != 'predict':
        parser.error('--arm-unavailable is only accepted by predict')
    if args.source_unavailable and args.phase != 'freeze':
        parser.error('--source-unavailable is only accepted by freeze')
    if args.unavailable and args.phase != 'finalize':
        parser.error('--unavailable is only accepted by finalize')
    runtime, output = args.runtime.resolve(), args.output.resolve()
    if args.phase == 'freeze':
        unavailable = cli_dispositions(parser, args.source_unavailable)
        binding = freeze(runtime, output, args.confirmation_sha256, unavailable)
        print(f"Frozen replication plan SHA-256 {binding['sha256']}")
        return
    digest = args.plan_sha256
    plan = load_plan(runtime, output, digest)
    seed = args.seed
    if args.phase in ('fit', 'predict', 'score'):
        if plan['seeds'].get(str(seed), {}).get('status') != 'planned':
            parser.error('phase requires a frozen seed and a promoted intervention')
        if args.phase == 'fit':
            fit(runtime, output, plan, digest, seed, args.arm, args.admission)
        elif args.phase == 'predict':
            unavailable = cli_dispositions(parser, args.arm_unavailable)
            predict(runtime, output, plan, digest, seed, unavailable)
        else:
            score(output, plan, digest, seed)
    elif args.phase == 'finalize':
        unavailable = cli_dispositions(parser, args.unavailable)
        print(finalize(output, plan, digest, unavailable)['status'])
    else:
        result = output / RESULT
        if result.exists():
            print(read_json(result)['status'])
        else:
            print('skipped' if plan['status'] == 'skipped' else 'planned')


if __name__ == '__main__':
    main()
