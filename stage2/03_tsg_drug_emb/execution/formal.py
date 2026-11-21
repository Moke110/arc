"""Register, run, resume and verify the 21 tenv cold-drug fits."""
import argparse
from collections import Counter
import fcntl
from pathlib import Path
import re
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

from fit import (ARMS, DOSE, MAX_EPOCHS, PARTITIONS, PATIENCE, SPLIT_SEED, STEP,
                 TRAINING_SEED, WIDTHS, Inputs, RepresentationInputs, dose,
                 phase_directory, read_json, runtime_identity, score_test, seal,
                 sha256, summary, verify_phase, write_json)
from tenv5k.protocol import mean_prediction


HERE = Path(__file__).resolve()
FIT = HERE.with_name('fit.py')
CAP = STEP.parent / '06_tsg_transfer/execution/capped_phase.sh'
# Sealed plans predate the transfer step's relocation to Step 06.
PRE_MOVE_DRIVER_SHA256 = '1159d78c887350c12d43e3485977dc10d66167a54c0d2e4786bd8cda30a3646d'
FOLDS = range(3)
FIT_MATRIX = [(fold, arm) for fold in FOLDS for arm in ARMS]
ACTIONS = {'train': 'fit', 'replay-validation': 'validation-replay',
           'evaluate': 'test', 'replay-test': 'test-replay'}
RESOURCE_FILES = frozenset(('command.txt', 'preflight.txt', 'resources-10s.log',
                            'scope.txt', 'stderr.txt', 'stdout.txt', 'terminal.txt', 'time.txt'))


def validate_fit_grid(records):
    if Counter((fit['fold'], fit['arm']) for fit in records) != Counter(FIT_MATRIX):
        raise ValueError('expected exactly seven independent fits in each of the three folds')
    for fold in FOLDS:
        group = [record for record in records if record['fold'] == fold]
        for key in ('order', 'transform', 'downstream_sha256', 'rng_sha256'):
            if any(record[key] != group[0][key] for record in group):
                raise ValueError(f'fold {fold} differs by representation arm: {key}')


def verified_inputs(path):
    approval = read_json(DOSE / 'log/247-axis-local-smoke.json')
    if sha256(path.parent / 'manifest.json') != approval['pair_sha256']:
        raise ValueError('tenv inputs are outside the approved prepared pair')
    pair = read_json(path.parent / 'manifest.json')
    if sha256(path / 'manifest.json') != pair['sources']['tenv']:
        raise ValueError('tenv prepared input digest mismatch')
    return Inputs(path)


def plan(args):
    inputs_digest = sha256(args.inputs / 'manifest.json')
    representation_digest = sha256(args.representations / 'manifest.json')
    smoke = read_json(args.smoke)
    group = smoke['fits']
    if (smoke['schema'] != 'tsg-drug-emb-bounded-smoke-v1' or smoke['status'] != 'PASS'
            or smoke['formal_fits'] != 0 or smoke['train_seed'] != TRAINING_SEED
            or smoke['device'] != 'cuda'
            or smoke['inputs_sha256'] != inputs_digest
            or smoke['representations_sha256'] != representation_digest
            or Counter((fit['arm'], fit['width']) for fit in group) != Counter(WIDTHS.items())
            or any(any(fit[key] != group[0][key] for fit in group)
                   for key in ('downstream_sha256', 'rng_sha256', 'split_order', 'transform'))):
        raise ValueError('formal run requires an accepted smoke on the same frozen inputs and seven arms')
    return dict(schema='tsg-drug-emb-formal-plan-v1', scenario='cold-drug', device='cuda',
                epochs=MAX_EPOCHS, patience=PATIENCE, seed=TRAINING_SEED, split_seed=SPLIT_SEED,
                inputs_sha256=inputs_digest, representations_sha256=representation_digest,
                fit_sha256=sha256(FIT), driver_sha256=sha256(HERE), source_sha256=dose.SOURCE_DIGEST,
                smoke_sha256=sha256(args.smoke), runs=[dict(fold=fold, arm=arm) for fold, arm in FIT_MATRIX])


def atomic_json(path, value):
    with tempfile.TemporaryDirectory(prefix='.write-receipt-', dir=path.parent) as temporary:
        candidate = Path(temporary) / path.name
        write_json(candidate, value)
        candidate.replace(path)


def registered(path, expected):
    if path.exists():
        if read_json(path) != expected:
            raise ValueError(f'previous registration or receipt differs: {path}')
    else:
        atomic_json(path, expected)


def coordinate(root, fold, arm):
    return root / 'runs' / f'fold-{fold}' / arm


def resources(root, identifier):
    attempts = root / 'resources/phases'
    successful = []
    for folder in sorted(attempts.glob(f'{identifier}-*')):
        terminal = folder / 'terminal.txt'
        if terminal.is_file() and terminal.read_text().endswith('status=0\n'):
            successful.append(folder)
    if len(successful) != 1:
        raise ValueError(f'{identifier}: expected one successful capped attempt')
    attempt, = successful
    entries = {path.name for path in attempt.iterdir() if path.is_file()}
    timing = (attempt / 'time.txt').read_text()
    memory = re.search(r'Maximum resident set size \(kbytes\): (\d+)', timing)
    elapsed = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)', timing)
    scope = (attempt / 'scope.txt').read_text()
    if (entries != RESOURCE_FILES or not (attempt / 'resources-10s.log').stat().st_size
            or 'Exit status: 0' not in timing or memory is None or elapsed is None
            or 'MemoryMax=10737418240' not in scope or 'MemorySwapMax=8589934592' not in scope):
        raise ValueError(f'{identifier}: missing or incomplete resource evidence')
    return dict(path=attempt.relative_to(root).as_posix(), elapsed=elapsed.group(1),
                peak_rss_kib=int(memory.group(1)),
                hashes={name: sha256(attempt / name) for name in sorted(entries)})


def ensure_phase(root, phase_path, phase, identifier, command):
    if not phase_path.exists():
        print(f'{identifier}: running under workstation caps', flush=True)
        result = subprocess.run(['bash', str(CAP), str(root / 'resources'), identifier, *command],
                                text=True, capture_output=True)
        if result.returncode:
            raise RuntimeError(f'{identifier} failed ({result.returncode}); resources: '
                               f'{result.stdout.strip()} {result.stderr.strip()}')
    receipt = verify_phase(phase_path, phase)
    resources(root, identifier)
    return receipt


def validate_fit(path, fold, arm, expected, rows, transform):
    fit = verify_phase(path / 'fit', 'fit')
    config = fit['config']
    checks = dict(arm=arm, fold=fold, width=WIDTHS[arm], seed=expected['seed'],
                  split_seed=expected['split_seed'], scenario='cold-drug', device='cuda',
                  inputs_sha256=expected['inputs_sha256'], max_epochs=MAX_EPOCHS, patience=PATIENCE,
                  representations_sha256=expected['representations_sha256'], order=rows, smoke_rows=None,
                  native_source_sha256=expected['source_sha256'])
    architecture = dict(features_dim=WIDTHS[arm], features_embed_dim=[400], method='film',
                        film_dimensions=[1, 32, 800], learning_rate=1e-4, weight_decay=1e-5,
                        batch_size=64, molecule_feature=arm, transform=transform)
    if (any(config[key] != value for key, value in checks.items())
            or any(config['architecture'][key] != value for key, value in architecture.items())
            or config['runtime']['training'] != expected['fit_sha256']
            or config['runtime'] != dict(runtime_identity(), training=expected['fit_sha256'],
                                         representations=sha256(STEP / 'output/api/tsg_drug_emb/representations.py'))):
        raise ValueError(f'formal fit configuration mismatch: {path}')
    init = fit['initialization']
    if init['tensors']['feat_embeddings.0.weight']['shape'] != [400, WIDTHS[arm]]:
        raise ValueError(f'first affine shape mismatch: {path}')
    curve = pd.read_csv(path / 'fit/learning_curve.csv', float_precision='round_trip')
    best = curve.loc[curve.validation_mse.idxmin()]
    n = len(curve)
    full, leftover = divmod(len(rows['train']), 64)
    steps = full + int(leftover > 1)
    if (fit['epochs'] != n or not 1 <= n <= MAX_EPOCHS
            or curve.epoch.tolist() != list(range(1, n + 1))
            or curve.optimizer_steps.tolist() != [steps * epoch for epoch in range(1, n + 1)]
            or fit['optimizer_steps'] != n * steps
            or not np.isfinite(curve.to_numpy()).all()
            or not np.all(curve.lr == 1e-4)
            or best.epoch != fit['best_epoch'] or best.validation_mse != fit['validation_mse']
            or (n != MAX_EPOCHS and n - int(best.epoch) != PATIENCE)):
        raise ValueError(f'invalid selection, training budget or optimizer step count: {path}')
    return fit


def select_fold(root, fold, expected, orders, transform):
    group = []
    for arm in ARMS:
        path = coordinate(root, fold, arm)
        fit = validate_fit(path, fold, arm, expected, orders, transform)
        replay = verify_phase(path / 'validation-replay', 'validation-replay')
        fit_sha = sha256(path / 'fit/manifest.json')
        if (replay['fit_sha256'] != fit_sha or replay['process_token'] == fit['process_token']
                or replay['prediction_sha256'] != fit['outputs']['validation_prediction.npy']
                or replay['validation_mse'] != fit['validation_mse']
                or replay['rows'] != len(orders['validation'])):
            raise ValueError(f'validation replay does not bind selected fit: {path}')
        group.append(dict(arm=arm, best_epoch=fit['best_epoch'], validation_mse=fit['validation_mse'],
                          fit_sha256=fit_sha, replay_sha256=sha256(path / 'validation-replay/manifest.json'),
                          downstream_sha256=fit['initialization']['downstream_sha256'],
                          rng_sha256=fit['initialization']['rng_sha256']))
    for field in ('downstream_sha256', 'rng_sha256'):
        if any(fit[field] != group[0][field] for fit in group):
            raise ValueError(f'fold {fold}: initial downstream architecture, tensors or RNG differ')
    return dict(schema='tsg-drug-emb-fold-validation-selection-v1', fold=fold,
                rule='first minimum inner validation MSE', fits=group)


def score_rows(directory, data, rows):
    prediction = np.load(directory / 'test_prediction.npy', mmap_mode='r')
    if (prediction.shape != (len(rows), 5000) or prediction.dtype != np.float32
            or not np.isfinite(prediction).all()):
        raise ValueError(f'test prediction shape or precision mismatch: {directory}')
    metrics = pd.read_csv(directory / 'test_per_sample.csv', float_precision='round_trip')
    reference = data.obs.iloc[rows]
    if (metrics.profile_id.tolist() != reference.profile_id.tolist()
            or metrics.drug_parent.tolist() != reference.drug_parent_inchikey.tolist()
            or metrics.cell_line.tolist() != reference.cell_line.tolist()
            or metrics.dose_um.tolist() != reference.dose_um.tolist()
            or metrics.nonfinite.any()
            or summary(metrics) != read_json(directory / 'test_metrics.json')):
        raise ValueError(f'test score profile identity/order/validity mismatch: {directory}')
    return len(metrics)


def mean_cell(args):
    data = verified_inputs(args.inputs)
    train = data.rows('cold-drug', args.fold, 'train')
    test = data.rows('cold-drug', args.fold, 'test')
    _, first, inverse = np.unique(data.obs.iloc[test].cell_line, return_index=True, return_inverse=True)
    reference = mean_prediction(data.expression, data.obs, data.controls, train, test[first],
                                'Mean_cell', 'cold-drug')
    prediction = reference[inverse]
    scores = score_test(data, test, prediction)
    expected = dict(arm='Mean_cell', fold=args.fold, rows=len(test), train_order=train.tolist(),
                    test_order=test.tolist(), inputs_sha256=sha256(args.inputs / 'manifest.json'),
                    runtime=runtime_identity(), rule='recomputed from exact fitting rows and matched controls')
    with phase_directory(args.destination / 'mean-cell', f'fold-{args.fold}') as directory:
        np.save(directory / 'test_prediction.npy', prediction)
        scores.to_csv(directory / 'test_per_sample.csv', index=False)
        write_json(directory / 'test_metrics.json', summary(scores))
        seal(directory, 'test', expected)
    print(f'Mean_cell fold {args.fold}: recomputed canonical reference', flush=True)


def common_arguments(args):
    return ['--inputs', str(args.inputs), '--representations', str(args.representations),
            '--source', str(args.source)]


def fold_orders(data, fold):
    return {name: data.rows('cold-drug', fold, name).tolist() for name in PARTITIONS}


def fold_transform(data, orders):
    return dose.fit_transform('film', data.doses(np.asarray(orders['train'])))


def run_matrix(args, data, expected):
    root = args.destination
    for fold in FOLDS:
        orders = fold_orders(data, fold)
        transform = fold_transform(data, orders)
        for actions in (('train', 'replay-validation'), ('evaluate', 'replay-test')):
            if actions[0] == 'evaluate':
                registered(root / f'selection-fold-{fold}.json',
                           select_fold(root, fold, expected, orders, transform))
            for arm in ARMS:
                path = coordinate(root, fold, arm)
                command = [sys.executable, str(FIT), *common_arguments(args), '--device', 'cuda',
                           '--arm', arm, '--fold', str(fold), '--run', str(path)]
                for action in actions:
                    phase = ACTIONS[action]
                    ensure_phase(root, path / phase, phase, f'fold-{fold}-{arm.lower()}-{action}',
                                 [*command, action])
        ensure_phase(root, root / 'mean-cell' / f'fold-{fold}', 'test', f'fold-{fold}-mean-cell',
                     [sys.executable, str(HERE), 'mean-cell', '--inputs', str(args.inputs), '--destination',
                      str(root), '--fold', str(fold)])


def verify_matrix(args, data, expected):
    root = args.destination
    records, baselines = [], []
    observed = {(path.parent.name, path.name)
                for path in (root / 'runs').glob('*/*') if path.is_dir()}
    if observed != {(f'fold-{fold}', arm) for fold, arm in FIT_MATRIX}:
        raise ValueError('formal run directory coverage mismatch')
    for fold in FOLDS:
        orders = fold_orders(data, fold)
        transform = fold_transform(data, orders)
        selection = select_fold(root, fold, expected, orders, transform)
        selection_path = root / f'selection-fold-{fold}.json'
        if read_json(selection_path) != selection:
            raise ValueError(f'fold {fold}: sealed validation choices differ')
        for selected in selection['fits']:
            arm = selected['arm']
            path = coordinate(root, fold, arm)
            fit = read_json(path / 'fit/manifest.json')
            validation = read_json(path / 'validation-replay/manifest.json')
            test = verify_phase(path / 'test', 'test')
            replay = verify_phase(path / 'test-replay', 'test-replay')
            if (len({phase['process_token'] for phase in (fit, validation, test, replay)}) != 4
                    or test['fit_sha256'] != selected['fit_sha256']
                    or test['validation_replay_sha256'] != selected['replay_sha256']
                    or test['rows'] != len(orders['test'])
                    or replay['fit_sha256'] != selected['fit_sha256']
                    or replay['validation_replay_sha256'] != selected['replay_sha256']
                    or replay['rows'] != test['rows']
                    or replay['test_sha256'] != sha256(path / 'test/manifest.json')
                    or replay['prediction_sha256'] != test['outputs']['test_prediction.npy']):
                raise ValueError(f'test replay does not match evaluation and selected fit: {path}')
            score_rows(path / 'test', data, orders['test'])
            records.append(dict(selected, fold=fold, order=orders, transform=transform,
                                epochs=fit['epochs'], optimizer_steps=fit['optimizer_steps'],
                                selection_sha256=sha256(selection_path),
                                test_sha256=sha256(path / 'test/manifest.json'),
                                test_replay_sha256=sha256(path / 'test-replay/manifest.json'),
                                resources={action: resources(root, f'fold-{fold}-{arm.lower()}-{action}')
                                           for action in ACTIONS}))
        reference = root / 'mean-cell' / f'fold-{fold}'
        receipt = verify_phase(reference, 'test')
        if (receipt['arm'] != 'Mean_cell' or receipt['fold'] != fold
                or receipt['inputs_sha256'] != expected['inputs_sha256']
                or receipt['train_order'] != orders['train'] or receipt['test_order'] != orders['test']
                or receipt['runtime'] != runtime_identity()
                or receipt['rows'] != score_rows(reference, data, orders['test'])):
            raise ValueError(f'fold {fold}: baseline compute contract differs')
        baselines.append(dict(fold=fold, manifest_sha256=sha256(reference / 'manifest.json'),
                              resources=resources(root, f'fold-{fold}-mean-cell')))
    validate_fit_grid(records)
    return dict(schema='tsg-drug-emb-formal-acceptance-v1', status='PASS', formal_fits=21,
                plan_sha256=sha256(root / 'formal-plan.json'), fits=records, mean_cell=baselines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    for command in ('run', 'verify'):
        group = commands.add_parser(command)
        for name in ('inputs', 'representations', 'source', 'smoke', 'destination'):
            group.add_argument('--' + name, type=Path, required=True)
    group = commands.add_parser('mean-cell')
    group.add_argument('--inputs', type=Path, required=True)
    group.add_argument('--destination', type=Path, required=True)
    group.add_argument('--fold', type=int, choices=FOLDS, required=True)
    args = parser.parse_args()
    for name, value in vars(args).items():
        if isinstance(value, Path):
            setattr(args, name, value.resolve())
    if args.command == 'mean-cell':
        mean_cell(args)
        return
    expected = plan(args)
    data = verified_inputs(args.inputs)
    RepresentationInputs(data, args.representations)
    root = args.destination
    plan_path = root / 'formal-plan.json'
    if args.command == 'verify':
        recorded = read_json(plan_path)
        if recorded.get('driver_sha256') == PRE_MOVE_DRIVER_SHA256:
            expected['driver_sha256'] = PRE_MOVE_DRIVER_SHA256
        if recorded != expected or verify_matrix(args, data, expected) != read_json(root / 'acceptance.json'):
            raise ValueError('formal acceptance differs from independent verification')
    else:
        root.mkdir(parents=True, exist_ok=True)
        with (root / '.orchestration.lock').open('a') as lease:
            fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
            registered(plan_path, expected)
            run_matrix(args, data, expected)
            registered(root / 'acceptance.json', verify_matrix(args, data, expected))
    print('PASS: 21 matched formal fits, exact fresh-process replay and 3 recomputed Mean_cell folds', flush=True)


if __name__ == '__main__':
    main()
