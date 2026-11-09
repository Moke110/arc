#!/usr/bin/env python3
"""Fit and evaluate the independent Step 06 positive-dose reset."""
import argparse
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

sys.dont_write_bytecode = True

import baseline
import transfer
from artifacts import SEEDS, read_json, sha256, write_csv, write_json
from evaluate import load_array
from formal import slot_lock, verify_fit_seal, completed_fit
from formal_matrix import bind, check_binding
from screen import SCREEN_SEED, paired_gain, screen_decision
from tahoe_bridge import (CHEMICAL_LOOKUP, DOSE_LOOKUP, DOSE_ROWS, dose_layout,
                          load_initialization, request_digest, state_digest)
from training_inputs import json_hash

STEP = Path('stage2/06_pdiff_transfer_optimize')
PLAN_SCHEMA = 'step06-dose-reset-plan-v1'
SCORE_SCHEMA = 'step06-dose-reset-scores-v1'
ARM = 'dose_reset'
DOSES = [0.2, 1.0, 10.0]


def load_dose_reset(model, directory, expected_sha256, recipient):
    """Apply the sealed bridge and restore the same construction's native dose rows."""
    import torch

    doses, _ = dose_layout(model.keys)
    parameter = model.get_parameter(DOSE_ROWS)
    if doses != DOSES or len(parameter) != len(DOSES) + 1 or not parameter.requires_grad:
        raise ValueError('expected three trainable positive tenv dose rows and a control row')

    native_rows = parameter.detach()[1:].clone()
    rng = torch.get_rng_state()
    bridge = load_initialization(model, directory, expected_sha256, recipient)
    mapped = model.state_dict()
    mapped_rows = mapped[DOSE_ROWS].clone()
    with torch.no_grad():
        parameter[1:].copy_(native_rows)

    restored = model.state_dict()
    reset_rows = restored[DOSE_ROWS]
    replay_rows = torch.cat((reset_rows[:1], mapped_rows[1:]))
    replay_sha = state_digest({**restored, DOSE_ROWS: replay_rows})
    lookup_sha = state_digest({key: restored[key] for key in (CHEMICAL_LOOKUP, DOSE_LOOKUP)})
    if (replay_sha != bridge['mapped_state_sha256']
            or lookup_sha != bridge['recipient_lookup_sha256']
            or not torch.equal(reset_rows[1:], native_rows)
            or not torch.equal(rng, torch.get_rng_state())):
        raise ValueError('dose reset changed the original bridge outside the permitted rows or consumed RNG')

    return bridge, {
        'reset_doses_um': doses,
        'mapped_state_sha256': replay_sha,
        'reset_state_sha256': state_digest(restored),
        'recipient_lookup_sha256': lookup_sha,
        'retained_zero_row_sha256': state_digest({DOSE_ROWS: reset_rows[:1]}),
        'native_positive_rows_sha256': state_digest({DOSE_ROWS: native_rows}),
        'restored_positive_rows_sha256': state_digest({DOSE_ROWS: reset_rows[1:]}),
        'torch_cpu_rng_preserved': True,
    }


def baseline_root(output, seed):
    return output / 'baseline' / f'fold-0-seed-{seed}'


def run_root(output, seed):
    return output / 'dose_reset' / f'fold-0-seed-{seed}'


def checkpoint_root(output, seed):
    return output / 'checkpoints' / f'01-dose-reset-fold-0-seed-{seed}'


def freeze_plan(output, seed, baseline_sha):
    if seed != SCREEN_SEED and screen_decision({17: score_gain(output, SCREEN_SEED)})['status'] != 'pending_replication':
        raise ValueError('additional seeds require a sealed seed 17 gain of at least +0.010')
    base = baseline_root(output, seed)
    pinned = bind(base / 'plan.json', baseline_sha)
    plan = read_json(pinned['path'])
    if plan['schema'] != baseline.SCHEMA or plan['seed'] != seed or plan['request']['identity']['fold'] != 0:
        raise ValueError('expected verified fold 0 validation baseline for this seed')
    history = read_json(check_binding(plan['checkpoints']['Transfer']['history'])['path'])
    lineage = history['initialization']
    source = plan['source']
    bridge = read_json(check_binding(source['bridge'])['path'])
    if (lineage['bridge_sha256'] != source['bridge']['sha256']
            or lineage['mapped_state_sha256'] != bridge['mapped_state_sha256']
            or lineage['source_checkpoint_sha256'] != source['checkpoint']['sha256']
            or bridge['source']['path'] != source['checkpoint']['path']
            or bridge['source']['checkpoint_sha256'] != source['checkpoint']['sha256']
            or lineage['source_sha256'] != source['source_sha256']
            or lineage['recipe_sha256'] != json_hash(plan['config'])):
        raise ValueError('baseline Transfer recipe, source and bridge identities differ')
    frozen = {'schema': PLAN_SCHEMA, 'seed': seed, 'baseline': pinned,
              'initialization': {**lineage, 'intervention': ARM, 'reset_doses_um': DOSES}}
    path = run_root(output, seed) / 'plan.json'
    with slot_lock(path.parent / '.plan-lock'), tempfile.TemporaryDirectory(dir=path.parent) as tmp:
        pending = Path(tmp) / 'plan.json'
        write_json(pending, frozen)
        if path.exists():
            bind(path, sha256(pending))
        else:
            pending.rename(path)
    print(f'Frozen dose reset plan SHA-256 {sha256(path)}', flush=True)


def selected_fit(destination, plan, base, digest):
    verify_fit_seal(destination, digest)
    identity = {**base['checkpoints']['Transfer']['metadata']['identity'], 'arm': ARM}
    run = completed_fit(destination, identity, base['config'])
    initialization = dict(run['initialization'])
    parity = initialization.pop('parity')
    if (initialization != plan['initialization']
            or parity['reset_doses_um'] != DOSES
            or parity['mapped_state_sha256'] != initialization['mapped_state_sha256']
            or parity['native_positive_rows_sha256'] != parity['restored_positive_rows_sha256']
            or parity['torch_cpu_rng_preserved'] is not True
            or run['training_initial_state'] != {
                'optimizer_state_entries': 0, 'scheduler_start_step': 0, 'ema_start_step': 0,
                'mapped_dose_rows_trainable': True, 'mapped_dose_gradient_observed': True}):
        raise ValueError('fit did not verify the permitted dose reset and fresh training state')
    name = f"01-ema-selected-step-{run['best_step']:05d}.pt"
    return {'checkpoint': bind(destination / name, run['checkpoint_sha256']),
            'history': bind(destination / 'run.json'), 'fit_seal': bind(destination / 'fit.json'),
            'metadata': {key: run[key] for key in baseline.CHECKPOINT_KEYS}}


def full_recipient(runtime, base):
    root = runtime / transfer.STEP05 / 'execution/output/cold-tenv-ticket76'
    entry = read_json(transfer.tahoe.historical_record('tenv') / 'plan.json')['folds'][0]
    request = read_json(bind(root / 'predictor_request.json', entry['request_sha256'])['path'])
    if baseline.inference_request(request) != base['request']:
        raise ValueError('frozen recipient request differs from the validation replay')
    return request


def fit_reset(runtime, output, plan, base, digest):
    from formal import use_runtime
    use_runtime(runtime)
    import native_training as native

    seed = plan['seed']
    destination = checkpoint_root(output, seed)
    with slot_lock(run_root(output, seed) / '.fit-lock'):
        if destination.exists():
            return selected_fit(destination, plan, base, digest)
        request = full_recipient(runtime, base)
        lineage = dict(plan['initialization'])
        bridge = check_binding(base['source']['bridge'])
        if request_digest(request) != lineage['recipient_request_sha256']:
            raise ValueError('full tenv recipient differs from the original Transfer bridge')
        worker_request = {**request, 'identity': {**transfer.transfer_identity(request, baseline.REGIME),
                                                   'arm': ARM}}
        config = base['config']
        identity = {**worker_request['identity'], 'seed': seed, 'stage': config['stage']}
        item = transfer.recipes(runtime, 'formal')[1][str(seed)]
        if item['config'] != config or item['recipe_sha256'] != lineage['recipe_sha256']:
            raise ValueError('dose reset must use the original tenv recipe')
        destination.parent.mkdir(parents=True, exist_ok=True)
        attempt = Path(tempfile.mkdtemp(prefix=f'{destination.name}-attempt-', dir=destination.parent))

        def initialize(model, directory, bridge_sha, recipient):
            bridge_record, parity = load_dose_reset(model, directory, bridge_sha, recipient)
            lineage['parity'] = parity
            return bridge_record

        with patch.object(transfer, 'load_initialization', initialize):
            transfer.fit_transfer(attempt, digest, runtime, worker_request, request, identity,
                                  item, Path(bridge['path']).parent, lineage, native)
        staged = attempt / 'fit'
        run = read_json(staged / 'run.json')
        (staged / 'best.pt').rename(staged / f"01-ema-selected-step-{run['best_step']:05d}.pt")
        selected_fit(staged, plan, base, digest)
        staged.rename(destination)
        attempt.rmdir()
        return selected_fit(destination, plan, base, digest)


def paired_summaries(arms):
    def differences(left, right):
        results = {}
        for metric in ('positive_ap100', 'negative_ap100', 'mean_ap100'):
            value = paired_gain(left[metric], right[metric])
            results[metric] = None if value is None else float(value)
        return results

    paired = []
    for reset, original, cm in zip(arms[ARM], arms['Transfer'], arms['CM'], strict=True):
        if any((row['dose_um'], row['group_count']) != (original['dose_um'], original['group_count'])
               for row in (reset, cm)):
            raise ValueError('paired dose strata differ from original Transfer')
        paired.append({'dose_um': original['dose_um'], 'group_count': original['group_count'],
                       'gain': differences(reset, original), 'cm_minus_transfer': differences(cm, original),
                       'cm_minus_intervention': differences(cm, reset)})
    return paired


def predict_reset(runtime, output, plan, base, digest):
    from formal import use_runtime
    use_runtime(runtime)
    seed = plan['seed']
    selected = selected_fit(checkpoint_root(output, seed), plan, base, digest)
    baseline.predict_checkpoint(base['request'], base['config'], base['inference'],
                                selected['checkpoint'], selected['metadata'],
                                run_root(output, seed) / 'prediction',
                                {'plan_sha256': digest, 'arm': ARM})


def score_reset(runtime, output, plan, base, digest):
    seed = plan['seed']
    root = run_root(output, seed)
    base_root = baseline_root(output, seed)
    selected = selected_fit(checkpoint_root(output, seed), plan, base, digest)
    checkpoints = {**base['checkpoints'], ARM: selected}
    values, evidence = {}, {}
    for arm, item in checkpoints.items():
        check_binding(item['checkpoint'])
        check_binding(item['history'])
        directory = root / 'prediction' if arm == ARM else base_root / arm
        arm_digest = digest if arm == ARM else plan['baseline']['sha256']
        values[arm] = baseline.verify_prediction(base['request'], item['checkpoint'], item['metadata'],
                                                 directory, {'plan_sha256': arm_digest, 'arm': arm})
        evidence[arm] = bind(directory / 'prediction.json')

    request = full_recipient(runtime, base)
    truth = bind(request['predictor_files']['validation_expression'],
                 request['file_sha256']['validation_expression'])
    observed = load_array(truth['path'], (base['case_count'], base['request']['n_genes']))
    group_rows, summary_rows, by_arm = [], [], {}
    for arm, prediction in values.items():
        groups, summary = baseline.score_prediction(base['request'], observed, prediction)
        labels = {'seed': seed, 'arm': arm, 'step': checkpoints[arm]['metadata']['best_step']}
        group_rows.extend({**labels, **row} for row in groups)
        summary_rows.extend({**labels, **row} for row in summary)
        by_arm[arm] = summary
    paired = paired_summaries(by_arm)
    pooled = paired[0]['gain']['mean_ap100']
    with slot_lock(root / '.score-lock'), tempfile.TemporaryDirectory(dir=root) as tmp:
        staging = Path(tmp) / 'scores'
        write_csv(staging / 'per_group.csv', group_rows)
        write_csv(staging / 'per_dose.csv', summary_rows)
        report = {'schema': SCORE_SCHEMA, 'fold': 0, 'seed': seed, 'comparison_mode': 'loss-selected',
                  'plan_sha256': digest, 'baseline': plan['baseline'], 'truth': truth,
                  'predictions': evidence, 'checkpoints': checkpoints,
                  'per_dose': summary_rows, 'paired': paired, 'overall_gain': pooled,
                  'files': {path.name: sha256(path) for path in sorted(staging.iterdir())}}
        write_json(staging / 'scores.json', report)
        destination = root / 'scores'
        if destination.exists():
            for path in staging.iterdir():
                bind(destination / path.name, sha256(path))
        else:
            staging.rename(destination)


def score_gain(output, seed):
    root = run_root(output, seed)
    path = root / 'scores/scores.json'
    report = read_json(path)
    if ((report['schema'], report['fold'], report['seed'], report['comparison_mode'])
            != (SCORE_SCHEMA, 0, seed, 'loss-selected')
            or report['plan_sha256'] != sha256(root / 'plan.json')):
        raise ValueError('dose reset score is detached from its frozen plan')
    for file in ('per_group.csv', 'per_dose.csv'):
        bind(path.parent / file, report['files'][file])
    pooled = {row['arm']: row['mean_ap100'] for row in report['per_dose'] if row['dose_um'] == 'pooled'}
    return paired_gain(pooled[ARM], pooled['Transfer'])


def write_decision(output, unavailable):
    gains, scores = {}, {}
    for seed in SEEDS:
        score = run_root(output, seed) / 'scores/scores.json'
        if score.exists():
            if seed in unavailable:
                raise ValueError(f'seed {seed} already has a score')
            gains[seed] = score_gain(output, seed)
            scores[seed] = bind(score)
        elif seed in unavailable:
            gains[seed] = None
    result = screen_decision(gains)
    root = output / 'dose_reset'
    with slot_lock(root / '.decision-lock'), tempfile.TemporaryDirectory(dir=root) as tmp:
        staging = Path(tmp) / 'decision.json'
        write_json(staging, {'schema': 'step06-screen-decision-v1', 'arm': ARM, 'fold': 0,
                             'comparison_mode': 'loss-selected', **result, 'scores': scores,
                             'overall_gains': {seed: None if gain is None else float(gain)
                                               for seed, gain in gains.items()},
                             'unavailable': unavailable})
        staging.replace(root / 'decision.json')
    print(result['status'], result['reason'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('verify', 'fit', 'predict', 'score', 'decision'))
    parser.add_argument('--runtime', type=Path, default=baseline.REPO)
    parser.add_argument('--seed', type=int, choices=SEEDS, default=SCREEN_SEED)
    parser.add_argument('--baseline-plan-sha256')
    parser.add_argument('--plan-sha256')
    parser.add_argument('--unavailable', nargs=2, metavar=('SEED', 'REASON'), action='append', default=[])
    args = parser.parse_args()
    runtime = args.runtime.resolve()
    output = runtime / STEP / 'output'
    if args.phase == 'decision':
        unavailable = {int(seed): reason for seed, reason in args.unavailable}
        if set(unavailable) - set(SEEDS) or any(not reason.strip() for reason in unavailable.values()):
            parser.error('--unavailable requires a screening seed and a terminal reason')
        write_decision(output, unavailable)
    elif args.phase == 'verify':
        if not args.baseline_plan_sha256:
            parser.error('verify requires --baseline-plan-sha256')
        freeze_plan(output, args.seed, args.baseline_plan_sha256)
    else:
        if not args.plan_sha256:
            parser.error(f'{args.phase} requires --plan-sha256')
        pinned = bind(run_root(output, args.seed) / 'plan.json', args.plan_sha256)
        plan = read_json(pinned['path'])
        if plan['schema'] != PLAN_SCHEMA or plan['seed'] != args.seed:
            raise ValueError('dose reset plan does not match the requested seed')
        base = read_json(check_binding(plan['baseline'])['path'])
        action = {'fit': fit_reset, 'predict': predict_reset, 'score': score_reset}[args.phase]
        action(runtime, output, plan, base, args.plan_sha256)


if __name__ == '__main__':
    main()
