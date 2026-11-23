"""Fit and independently replay Stage 03 comparison baselines on frozen rows."""
import argparse
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from formal import SCENARIOS, candidate_lock, memory_preflight, read_json, write_json
from tsg_round1.inputs import Inputs, SOURCE_POPULATIONS, sha256
from tsg_round1.learner import environment
from tsg_round1.protocol import score_prediction, summary
from tenv5k import model_mlp
from tenv5k.learning import EarlyStopping

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP.parent / 'util/output/api'))
from stage2_util.ecfp4 import ECFP4_MLP, open_ecfp4, publish_ecfp4

MEAN_CELL = 'Mean_cell'
METHODS = (MEAN_CELL, ECFP4_MLP)
SEED = 42
SCHEMA = 'tsg-dose-round1-baseline-v1'
PARTITIONS = ('validation', 'test')


def chunks(rows, batch_size):
    for start in range(0, len(rows), batch_size):
        yield rows[start:start + batch_size]


def mean_cell(data, fit_rows, test_rows):
    """Pool fitted treated expression and unique original vehicles per cell line."""
    fit_cells = data.obs.iloc[fit_rows].cell_line.to_numpy()
    test_cells = data.obs.iloc[test_rows].cell_line.to_numpy()
    result = np.empty((len(test_rows), data.expression.shape[1]), dtype=np.float32)
    for cell in sorted(set(test_cells)):
        treated = fit_rows[fit_cells == cell]
        if not len(treated):
            raise ValueError(f'missing training cell coverage: {cell}')
        vehicles = sorted({control for row in treated for control in data.controls[row]})
        rows = np.concatenate([treated, np.asarray(vehicles, dtype=int)])
        total = np.zeros(data.expression.shape[1], dtype=np.float64)
        for batch in chunks(rows, 256):
            total += data.expression[batch].sum(axis=0, dtype=np.float64)
        result[test_cells == cell] = total / len(rows)
    return result


class FingerprintInputs:
    def __init__(self, data, candidate):
        self.data = data
        self.expression = data.expression
        self.ecfp4 = open_ecfp4(candidate, source_manifest=data.directory / 'manifest.json')

    def fingerprints(self, rows):
        parents = self.data.obs.iloc[rows].drug_parent_inchikey.map(self.data.parent_lookup)
        return np.asarray(self.ecfp4[parents.to_numpy()], dtype=np.float32)


def run_path(root, source, scenario, fold, method):
    return root / 'runs' / source / scenario / str(fold) / method


def frozen_pair(directory):
    manifest = directory / 'manifest.json'
    evidence = read_json(STEP / 'log/247-axis-local-smoke.json')
    pair = read_json(manifest)
    if (sha256(manifest) != evidence['pair_sha256']
            or pair['sources'] != {source: sha256(directory / source / 'manifest.json')
                                  for source in SOURCE_POPULATIONS}):
        raise ValueError('baseline inputs differ from accepted Stage 03 pair')
    return sha256(manifest)


def configuration(args, data, fingerprints):
    import ecfp4
    import tenv5k.learning
    import tsg_round1.inputs
    import tsg_round1.protocol

    paths = [Path(__file__), Path(ecfp4.__file__), Path(model_mlp.__file__),
             Path(tenv5k.learning.__file__), Path(tsg_round1.inputs.__file__),
             Path(tsg_round1.protocol.__file__)]
    return dict(schema=SCHEMA, source=data.source, scenario=args.scenario, fold=args.fold,
                method=args.method, seed=SEED, device=args.device,
                input_manifest_sha256=sha256(data.directory / 'manifest.json'),
                fingerprint_manifest_sha256=sha256(fingerprints / 'manifest.json'),
                source_pair_sha256=frozen_pair(args.inputs),
                runtime={path.name: sha256(path) for path in paths},
                environment=environment(args.device), max_epochs=1000, patience=50,
                architecture=model_mlp.CONFIG['architecture'], batch_size=64,
                compute_dtype='float32', learning_rate=1e-4, weight_decay=1e-5,
                loss='MSE', scheduler=None, fingerprint_only=True)


def fit_rows(data, config, partition):
    return data.rows(config['scenario'], config['fold'], partition)


def mlp_prediction(model, features, rows, device):
    return model_mlp.sampled_prediction(model, features, rows, seed=SEED, device=device)


def validation_mse(prediction, expression, rows):
    return float(np.mean(np.square(prediction - expression[rows])))


def fit_mlp(candidate, data, features, config):
    import torch

    device = config['device']
    model = model_mlp.build_model(None, SEED, device)
    initial_weights = {key: value.detach().cpu().clone() for key, value in model.state_dict().items()}
    optimizer = model_mlp.make_optimizer(model)
    training = fit_rows(data, config, 'train')
    validation = fit_rows(data, config, 'validation')
    stopping = EarlyStopping(config['patience'], config['max_epochs'])
    shuffle = np.random.default_rng(SEED)
    curves = []
    best_weights = None
    while not stopping.stopped:
        epoch = len(curves) + 1
        loss = model_mlp.train_epoch(model, features, training, optimizer, shuffle, device)
        prediction = mlp_prediction(model, features, validation, device)
        mse = validation_mse(prediction, data.expression, validation)
        curves.append(dict(epoch=epoch, train_loss=loss, validation_mse=mse))
        if stopping.observe(mse, epoch):
            best_weights = {key: value.detach().cpu().clone() for key, value in model.state_dict().items()}
        if epoch % 50 == 0 or stopping.stopped:
            print(f'{data.source} {config["scenario"]} fold {config["fold"]}: epoch {epoch}', flush=True)
    model.load_state_dict(best_weights)
    if all(torch.equal(weight, initial_weights[key]) for key, weight in best_weights.items()):
        raise ValueError('ECFP4_MLP best weights did not update')
    restored = mlp_prediction(model, features, validation, device)
    if validation_mse(restored, data.expression, validation) != stopping.best:
        raise ValueError('restored ECFP4_MLP differs from best validation MSE')
    torch.save(best_weights, candidate / 'checkpoint.pt')
    pd.DataFrame(curves).to_csv(candidate / 'learning_curve.csv', index=False)
    return model, dict(epochs=len(curves), best_epoch=stopping.best_epoch,
                       validation_mse=stopping.best)


def predict_rows(data, features, config, partition, model):
    rows = fit_rows(data, config, partition)
    if config['method'] == MEAN_CELL:
        return mean_cell(data, fit_rows(data, config, 'train'), rows)
    return mlp_prediction(model, features, rows, config['device'])


def score_rows(data, config, partition, prediction):
    rows = fit_rows(data, config, partition)
    return pd.concat([score_prediction(data.obs.iloc[rows[batch]], data.expression[rows[batch]],
                                      prediction[batch], data.control[rows[batch]],
                                      prediction_ids=data.obs.index[rows[batch]])
                      for batch in chunks(np.arange(len(rows)), 256)], ignore_index=True)


def output_hashes(candidate):
    return {path.name: sha256(path) for path in sorted(candidate.iterdir())
            if path.is_file() and path.name != 'manifest.json'}


def execute_fit(candidate, data, features, config):
    model, fit = (fit_mlp(candidate, data, features, config) if config['method'] == ECFP4_MLP
                  else (None, {}))
    partitions = PARTITIONS if model is not None else ('test',)
    for partition in partitions:
        prediction = predict_rows(data, features, config, partition, model)
        if prediction.dtype != np.float32 or not np.isfinite(prediction).all():
            raise ValueError('non-finite or non-float32 baseline prediction')
        np.save(candidate / f'{partition}_prediction.npy', prediction)
        scores = score_rows(data, config, partition, prediction)
        scores.to_csv(candidate / f'{partition}_per_sample.csv', index=False)
        write_json(candidate / f'{partition}_metrics.json', summary(scores))
    write_json(candidate / 'manifest.json', dict(config=config, fit=fit, outputs=output_hashes(candidate)))


def replay(candidate, data, features, config):
    import torch

    receipt = read_json(candidate / 'manifest.json')
    if receipt['config'] != config or receipt['outputs'] != output_hashes(candidate):
        raise ValueError('baseline fit config or output digest mismatch')
    model = None
    partitions = ('test',)
    if config['method'] == ECFP4_MLP:
        model = model_mlp.build_model(None, SEED, config['device'])
        model.load_state_dict(torch.load(candidate / 'checkpoint.pt', weights_only=True))
        partitions = PARTITIONS
    for partition in partitions:
        result = predict_rows(data, features, config, partition, model)
        expected = np.load(candidate / f'{partition}_prediction.npy', mmap_mode='r')
        if not np.array_equal(result, expected):
            raise ValueError(f'{partition} baseline prediction failed exact replay')
        scored = score_rows(data, config, partition, result)
        old = pd.read_csv(candidate / f'{partition}_per_sample.csv',
                          dtype={'profile_id': str, 'drug_parent': str, 'cell_line': str},
                          float_precision='round_trip')
        pd.testing.assert_frame_equal(scored, old, check_dtype=False, check_exact=True)
    return dict(schema=SCHEMA, status='exact fresh-process replay', partitions=list(partitions))


def unit(args):
    import torch

    if args.device.startswith('cuda'):
        torch.use_deterministic_algorithms(True)
        torch.backends.cuda.matmul.allow_tf32 = False
    data = Inputs(args.inputs / args.source)
    fingerprints = args.candidate / 'ecfp4' / args.source
    features = FingerprintInputs(data, fingerprints)
    config = configuration(args, data, fingerprints)
    target = run_path(args.candidate, args.source, args.scenario, args.fold, args.method)
    if args.unit == 'replay':
        result = replay(target, data, features, config)
        result['manifest_sha256'] = sha256(target / 'manifest.json')
        write_json(target.parent / f'{args.method}-replay.json', result)
    else:
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            raise ValueError(f'baseline run already exists: {target}')
        with TemporaryDirectory(prefix='.baseline-', dir=target.parent) as temporary:
            path = Path(temporary)
            execute_fit(path, data, features, config)
            path.rename(target)


def timed_unit(args, source, scenario, fold, method, phase, lock_fd):
    memory_preflight()
    target = run_path(args.candidate, source, scenario, fold, method)
    resources = args.candidate / 'resources' / target.relative_to(args.candidate / 'runs')
    resources.mkdir(parents=True, exist_ok=True)
    command = ['/usr/bin/time', '-v', '-o', str(resources / f'{phase}.time'),
               sys.executable, str(Path(__file__).resolve()), '--inputs', str(args.inputs),
               '--candidate', str(args.candidate), '--device', args.device, '--unit', phase,
               '--source', source, '--scenario', scenario, '--fold', str(fold), '--method', method]
    subprocess.run(command, check=True, pass_fds=(lock_fd,))


def drive(args):
    frozen_pair(args.inputs)
    with candidate_lock(args.candidate) as lock_fd:
        for source in SOURCE_POPULATIONS:
            fingerprints = args.candidate / 'ecfp4' / source
            if not fingerprints.exists():
                publish_ecfp4(args.inputs / source, fingerprints)
            for scenario in SCENARIOS:
                for fold in range(3):
                    for method in METHODS:
                        print(f'{source} {scenario} fold {fold} {method}', flush=True)
                        target = run_path(args.candidate, source, scenario, fold, method)
                        if not target.exists():
                            timed_unit(args, source, scenario, fold, method, 'fit', lock_fd)
                        timed_unit(args, source, scenario, fold, method, 'replay', lock_fd)


def main():
    os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--device', default='cuda')
    parser.add_argument('--unit', choices=('fit', 'replay'), help=argparse.SUPPRESS)
    parser.add_argument('--source', choices=SOURCE_POPULATIONS)
    parser.add_argument('--scenario', choices=SCENARIOS)
    parser.add_argument('--fold', type=int, choices=range(3))
    parser.add_argument('--method', choices=METHODS)
    args = parser.parse_args()
    if args.unit:
        if None in (args.source, args.scenario, args.fold, args.method):
            parser.error('an individual fit or replay requires source, scenario, fold and method')
        unit(args)
    else:
        drive(args)


if __name__ == '__main__':
    main()
