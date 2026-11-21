"""Run and replay the tenv KPGT_MLP and Mean_cell controls."""

import argparse
import os
from pathlib import Path
import secrets
import subprocess
import sys
from tempfile import TemporaryDirectory
import time

import numpy as np
import pandas as pd
import torch

from formal import (
    BASELINES,
    EXECUTION,
    STEP,
    candidate_lock,
    memory_preflight,
    planned_runs,
    read_json,
    write_json,
)
from tenv5k import model_mlp
from tenv5k.learning import EarlyStopping
from tsg_kpgt.inputs import Inputs, WIDTH, sha256
from tsg_kpgt.learner import environment, runtime_identity

sys.path.append(str(EXECUTION / 'axis_local'))
from baselines import mean_cell, output_hashes, score_rows, validation_mse

MEAN_CELL = 'Mean_cell'
MLP = 'KPGT_MLP'
SCHEMA = 'tsg-dose-tenv-kpgt-baseline-v1'
SEED = 42
ARCHITECTURE = [WIDTH, 128, 5000]
CONFIG = dict(
    architecture=ARCHITECTURE,
    batch_size=64,
    compute_dtype='float32',
    learning_rate=1e-4,
    weight_decay=1e-5,
    loss='MSE',
    scheduler=None,
    max_epochs=1000,
    patience=50,
    chemical_input_only=True,
)
PHASES = ('fit', 'replay')


def configuration(args, origin):
    source = Path(sys.modules[mean_cell.__module__].__file__)
    manifest = sha256(args.inputs / 'manifest.json')
    if manifest != sha256(STEP / 'log/274-kpgt-smoke/input-manifest.json'):
        raise ValueError('baseline inputs differ from accepted #274 generation')
    return dict(
        schema=SCHEMA,
        **origin,
        input_manifest_sha256=manifest,
        runtime=runtime_identity(),
        runner_sha256=sha256(Path(__file__)),
        reference_sha256=sha256(source),
        model_mlp_sha256=model_mlp.SOURCE_DIGEST,
        environment=environment(args.device),
        device=args.device,
        training=CONFIG if origin['method'] == MLP else None,
    )


def build_model(device):
    torch.manual_seed(SEED)
    return torch.nn.Sequential(
        torch.nn.Linear(WIDTH, 128), torch.nn.ReLU(), torch.nn.Linear(128, 5000)
    ).to(device=device, dtype=torch.float32)


def fit_mlp(path, data, config):
    device = config['device']
    model = build_model(device)
    initial = {
        key: value.detach().cpu().clone() for key, value in model.state_dict().items()
    }
    optimizer = model_mlp.make_optimizer(model)
    train = data.rows(config['scenario'], config['fold'], 'train')
    validation = data.rows(config['scenario'], config['fold'], 'validation')
    settings = config['training']
    stopper = EarlyStopping(
        patience=settings['patience'], max_epochs=settings['max_epochs']
    )
    shuffle = np.random.default_rng(SEED)
    curves = []
    best_weights = None
    while not stopper.stopped:
        epoch = len(curves) + 1
        loss = model_mlp.train_epoch(model, data, train, optimizer, shuffle, device)
        prediction = model_mlp.sampled_prediction(
            model, data, validation, seed=SEED, device=device
        )
        mse = validation_mse(prediction, data.expression, validation)
        curves.append(dict(epoch=epoch, train_loss=loss, validation_mse=mse))
        if stopper.observe(mse, epoch):
            best_weights = {
                k: v.detach().cpu().clone() for k, v in model.state_dict().items()
            }
        if epoch % 50 == 0 or stopper.stopped:
            print(
                f'{config["scenario"]} fold {config["fold"]}: epoch {epoch}', flush=True
            )
    model.load_state_dict(best_weights)
    if all(torch.equal(value, initial[key]) for key, value in best_weights.items()):
        raise ValueError('KPGT_MLP best weights did not update')
    restored = model_mlp.sampled_prediction(
        model, data, validation, seed=SEED, device=device
    )
    if validation_mse(restored, data.expression, validation) != stopper.best:
        raise ValueError('restored KPGT_MLP validation MSE differs from the best epoch')
    torch.save(best_weights, path / 'checkpoint.pt')
    pd.DataFrame(curves).to_csv(path / 'learning_curve.csv', index=False)
    return model, dict(
        epochs=len(curves),
        best_epoch=stopper.best_epoch,
        validation_mse=stopper.best,
        weights_updated=True,
    )


def partitions(config):
    return ('validation', 'test') if config['method'] == MLP else ('test',)


def predict(data, config, partition, model):
    scenario, fold = config['scenario'], config['fold']
    rows = data.rows(scenario, fold, partition)
    if config['method'] == MEAN_CELL:
        result = mean_cell(data, data.rows(scenario, fold, 'train'), rows)
    else:
        result = model_mlp.sampled_prediction(
            model, data, rows, seed=SEED, device=config['device']
        )
    if result.dtype != np.float32 or not np.isfinite(result).all():
        raise ValueError('baseline produced a non-finite or non-float32 prediction')
    return result


def fit(path, data, config):
    from tsg_round1.protocol import summary

    model, evidence = (
        fit_mlp(path, data, config) if config['method'] == MLP else (None, {})
    )
    for partition in partitions(config):
        result = predict(data, config, partition, model)
        np.save(path / f'{partition}_prediction.npy', result)
        scores = score_rows(data, config, partition, result)
        scores.to_csv(path / f'{partition}_per_sample.csv', index=False)
        write_json(path / f'{partition}_metrics.json', summary(scores))
    write_json(
        path / 'manifest.json',
        dict(
            config=config,
            fit=evidence,
            outputs=output_hashes(path),
            process_token=secrets.token_hex(32),
        ),
    )


def verify_fit(path, config):
    receipt = read_json(path / 'manifest.json')
    expected = {
        f'{partition}_{artifact}'
        for partition in partitions(config)
        for artifact in ('prediction.npy', 'per_sample.csv', 'metrics.json')
    }
    if config['method'] == MLP:
        expected.update(('checkpoint.pt', 'learning_curve.csv'))
    if (
        receipt['config'] != config
        or set(receipt['outputs']) != expected
        or receipt['outputs'] != output_hashes(path)
        or len(receipt['process_token']) != 64
    ):
        raise ValueError('baseline fit config or output identity mismatch')
    return receipt


def replay(path, data, config):
    receipt = verify_fit(path, config)
    model = None
    if config['method'] == MLP:
        model = build_model(config['device'])
        model.load_state_dict(
            torch.load(path / 'checkpoint.pt', weights_only=True, map_location='cpu')
        )
    for partition in partitions(config):
        result = predict(data, config, partition, model)
        expected = np.load(path / f'{partition}_prediction.npy', mmap_mode='r')
        if not np.array_equal(result, expected):
            raise ValueError(f'{partition} baseline prediction failed exact replay')
        scores = score_rows(data, config, partition, result)
        old = pd.read_csv(
            path / f'{partition}_per_sample.csv',
            dtype={'profile_id': str, 'drug_parent': str, 'cell_line': str},
            float_precision='round_trip',
        )
        pd.testing.assert_frame_equal(scores, old, check_dtype=False, check_exact=True)
    return dict(
        schema=SCHEMA,
        partitions=list(partitions(config)),
        status='exact fresh-process replay',
        process_token=secrets.token_hex(32),
        fit_process_token=receipt['process_token'],
        manifest_sha256=sha256(path / 'manifest.json'),
    )


def paths(root, origin):
    coordinate = Path(origin['scenario']) / f'fold-{origin["fold"]}' / origin['method']
    return root / 'runs' / coordinate, root / 'resources' / coordinate


def run_phase(args, origin, phase):
    data = Inputs(args.inputs)
    config = configuration(args, origin)
    path, _ = paths(args.candidate, origin)
    if args.device.startswith('cuda'):
        torch.use_deterministic_algorithms(True)
        torch.backends.cuda.matmul.allow_tf32 = False
    if phase == 'replay':
        receipt = replay(path, data, config)
        write_json(path.parent / f'{origin["method"]}-replay.json', receipt)
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(prefix='.baseline-', dir=path.parent) as temporary:
            stage = Path(temporary)
            fit(stage, data, config)
            stage.rename(path)


def invoke(args, origin, phase, lock_fd):
    available = memory_preflight()
    _, resources = paths(args.candidate, origin)
    attempt = resources / f'{phase}-{time.time_ns()}'
    attempt.mkdir(parents=True)
    started = time.time()
    write_json(
        attempt / 'phase.json',
        dict(
            origin=origin,
            phase=phase,
            available_ram_kib=available,
            started_utc_unix_seconds=started,
        ),
    )
    command = [
        '/usr/bin/time',
        '-v',
        '-o',
        str(attempt / 'time.txt'),
        sys.executable,
        str(Path(__file__).resolve()),
        '--inputs',
        str(args.inputs),
        '--candidate',
        str(args.candidate),
        '--device',
        args.device,
        '--unit',
        phase,
        '--scenario',
        origin['scenario'],
        '--fold',
        str(origin['fold']),
        '--method',
        origin['method'],
    ]
    print(
        f'{origin["scenario"]} fold {origin["fold"]} {origin["method"]} {phase}',
        flush=True,
    )
    with (attempt / 'subprocess.log').open('w') as output:
        result = subprocess.run(
            command, stdout=output, stderr=subprocess.STDOUT, pass_fds=(lock_fd,)
        )
    write_json(
        attempt / 'result.json',
        dict(
            returncode=result.returncode,
            elapsed_seconds=time.time() - started,
            time_sha256=sha256(attempt / 'time.txt'),
        ),
    )
    result.check_returncode()


def verify_replay(path, origin, fit):
    receipt = read_json(path.parent / f'{origin["method"]}-replay.json')
    if (
        receipt['schema'] != SCHEMA
        or receipt['status'] != 'exact fresh-process replay'
        or receipt['partitions'] != list(partitions(fit['config']))
        or receipt['fit_process_token'] != fit['process_token']
        or receipt['process_token'] == fit['process_token']
        or len(receipt['process_token']) != 64
        or receipt['manifest_sha256'] != sha256(path / 'manifest.json')
    ):
        raise ValueError('baseline replay is unbound or not in a fresh process')
    return receipt


def drive(args, origins):
    Inputs(args.inputs)
    with candidate_lock(args.candidate) as lock_fd:
        for origin in origins:
            path, _ = paths(args.candidate, origin)
            if not path.exists():
                invoke(args, origin, 'fit', lock_fd)
            receipt = verify_fit(path, configuration(args, origin))
            if not (path.parent / f'{origin["method"]}-replay.json').exists():
                invoke(args, origin, 'replay', lock_fd)
            verify_replay(path, origin, receipt)


def main():
    os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--device', default='cuda')
    parser.add_argument('--unit', choices=PHASES, help=argparse.SUPPRESS)
    parser.add_argument('--scenario', choices=('cold-drug', 'cold-drug-dose'))
    parser.add_argument('--fold', type=int, choices=range(3))
    parser.add_argument('--method', choices=BASELINES)
    args = parser.parse_args()
    origins = planned_runs()['baselines']
    if args.unit:
        if None in (args.scenario, args.fold, args.method):
            parser.error('--unit requires a scenario, fold and method')
        origin = next(
            row
            for row in origins
            if (row['scenario'], row['fold'], row['method'])
            == (args.scenario, args.fold, args.method)
        )
        run_phase(args, origin, args.unit)
    else:
        if any(value is not None for value in (args.scenario, args.fold, args.method)):
            parser.error('scenario, fold and method are only accepted by --unit')
        drive(args, origins)


if __name__ == '__main__':
    main()
