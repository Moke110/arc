"""Fit and independently replay a registered tenv cold-drug representation arm."""
import argparse
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import secrets
import sys
import tempfile
import time
from types import MethodType

import numpy as np
import pandas as pd
import torch
from torch import nn

STEP = Path(__file__).resolve().parents[1]
DOSE = STEP.parent / '04_tsg_dose'

sys.path[:0] = [str(STEP / 'output/api'), str(DOSE / 'output/api')]
from tsg_drug_emb.representations import ARMS, WIDTHS, RepresentationInputs, sha256
from tsg_dose import model as dose
from tsg_round1.inputs import Inputs
from tsg_round1.learner import environment, predict, runtime_identity, smoke_subset
from tsg_round1.protocol import SPLIT_SEED, TRAINING_SEED, score_prediction, summary
from tenv5k.learning import EarlyStopping


PROCESS_TOKEN = secrets.token_hex(16)
PARTITIONS = ('train', 'validation', 'test')
FIRST_LAYER = frozenset(('feat_embeddings.0.weight', 'feat_embeddings.0.bias'))
VALIDATION_SEED = 22617
TEST_SEED = TRAINING_SEED + 100000
PATIENCE = 50
MAX_EPOCHS = 1000
PHASE_FILES = {
    'fit': {'checkpoint.pt', 'learning_curve.csv', 'validation_prediction.npy'},
    'validation-replay': set(),
    'test': {'test_prediction.npy', 'test_per_sample.csv', 'test_metrics.json'},
    'test-replay': set(),
}


def packed_film_forward(model, x1, features):
    """Keep the accepted FiLM forward order with a width-independent dose slot."""
    width = model.feat_embeddings[0].in_features
    if features.ndim != 2 or features.shape[1] != width + 1:
        raise ValueError('dose transport requires the native-width representation plus a scalar')
    chemical_features, scalar = features[:, :width], features[:, width:]
    z1, mu1, logvar1 = model.encode_x1(x1)
    x1_rec = model.decode_x1(z1)
    chemical = model.feat_embeddings(chemical_features)
    gamma, beta = model.film(scalar).chunk(2, dim=1)
    conditioned = torch.cat([z1, (1 + gamma) * chemical + beta], dim=1)
    mu_pred = model.mu_z2Fz1(conditioned)
    logvar_pred = model.logvar_z2Fz1(conditioned)
    z2_pred = model.sample_latent(mu_pred, logvar_pred)
    return x1_rec, mu1, logvar1, model.decode_x2(z2_pred), mu_pred, logvar_pred, z2_pred


class RepresentationAdapter(dose.DoseAdapter):
    def __init__(self, source, transform, width):
        super().__init__(source, 'film', transform)
        self.width = width

    def build(self, seed, device):
        network = super().build(seed, device)
        first = network.feat_embeddings[0]
        if self.width != first.in_features:
            # Constructor draws for this replacement cannot advance the shared
            # CPU or CUDA streams used by native latent sampling and dropout.
            with torch.random.fork_rng(devices=[]):
                layer = nn.Linear(self.width, first.out_features, device='cpu', dtype=torch.float32)
            network.feat_embeddings[0] = layer.to(device)
        network.forward = MethodType(packed_film_forward, network)
        return network


class ArmInputs:
    """Serve representations by treated row while preserving frozen data access."""
    def __init__(self, data, representations, arm):
        self.data = data
        self.representations = representations
        self.arm = arm

    def __getattr__(self, name):
        return getattr(self.data, name)

    def fingerprints(self, rows):
        return self.representations.vectors(self.arm, rows)


def canonical_digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def initialization(model):
    tensors = {}
    for name, tensor in model.state_dict().items():
        array = tensor.detach().cpu().contiguous().numpy()
        tensors[name] = dict(shape=list(array.shape), dtype=str(array.dtype),
                             sha256=hashlib.sha256(array.tobytes()).hexdigest())
    device = next(model.parameters()).device
    rng = {'cpu': hashlib.sha256(torch.get_rng_state().numpy().tobytes()).hexdigest()}
    if device.type == 'cuda':
        rng['cuda'] = hashlib.sha256(torch.cuda.get_rng_state(device).numpy().tobytes()).hexdigest()
    return dict(tensors=tensors, rng_sha256=rng, downstream_sha256=canonical_digest({
        name: value for name, value in tensors.items() if name not in FIRST_LAYER
    }))


def read_json(path):
    return json.loads(path.read_text())


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


@contextmanager
def phase_directory(run, phase):
    run.mkdir(parents=True, exist_ok=True)
    target = run / phase
    if target.exists():
        raise FileExistsError(f'phase already sealed: {target}')
    with tempfile.TemporaryDirectory(prefix=f'.{phase}-', dir=run) as temporary:
        staging = Path(temporary) / phase
        staging.mkdir()
        yield staging
        staging.rename(target)


def seal(directory, phase, evidence):
    files = {path.name: sha256(path) for path in sorted(directory.iterdir())}
    receipt = dict(schema='tsg-drug-emb-phase-v1', phase=phase, process_token=PROCESS_TOKEN,
                   outputs=files, **evidence)
    write_json(directory / 'manifest.json', receipt)
    return receipt


def verify_phase(directory, phase):
    receipt = read_json(directory / 'manifest.json')
    paths = {path.name: path for path in directory.iterdir() if path.name != 'manifest.json'}
    if (receipt.get('schema') != 'tsg-drug-emb-phase-v1' or receipt.get('phase') != phase
            or set(paths) != set(receipt['outputs']) or set(paths) != PHASE_FILES[phase]):
        raise ValueError(f'invalid phase inventory: {directory}')
    for name, path in paths.items():
        if path.is_symlink() or not path.is_file() or sha256(path) != receipt['outputs'][name]:
            raise ValueError(f'digest mismatch: {path}')
    return receipt


def load_arm(args):
    approved = read_json(DOSE / 'log/247-axis-local-smoke.json')
    if sha256(args.inputs.parent / 'manifest.json') != approved['pair_sha256']:
        raise ValueError('inputs are outside the accepted axis-local pair')
    pair = read_json(args.inputs.parent / 'manifest.json')
    input_digest = sha256(args.inputs / 'manifest.json')
    if input_digest != pair['sources']['tenv']:
        raise ValueError('prepared tenv input digest differs from the accepted pair')

    data = Inputs(args.inputs)
    representations = RepresentationInputs(data, args.representations)
    rows = {partition: data.rows('cold-drug', args.fold, partition) for partition in PARTITIONS}
    parents = {partition: set(data.obs.iloc[part].drug_parent_inchikey)
               for partition, part in rows.items()}
    if sum(map(len, parents.values())) != len(set.union(*parents.values())):
        raise ValueError('parent leakage between fold partitions')
    if args.smoke_rows is not None:
        rows = {partition: smoke_subset(data, part, args.smoke_rows)
                for partition, part in rows.items()}
    transform = dose.fit_transform('film', data.doses(rows['train']))
    runtime = dict(runtime_identity(), training=sha256(Path(__file__)),
                   representations=sha256(STEP / 'output/api/tsg_drug_emb/representations.py'))
    architecture = dose.configuration('film', transform)
    architecture.update(features_dim=representations.width(args.arm), molecule_feature=args.arm,
                        method_identity='tsg-drug-emb-film-v1')
    architecture['adaptations'] = [
        adaptation for adaptation in architecture['adaptations']
        if adaptation != 'shared Morgan2048 ECFP4'
    ] + ['native-width representation with no normalization']
    config = dict(arm=args.arm, width=representations.width(args.arm), scenario='cold-drug',
                  fold=args.fold, seed=TRAINING_SEED, split_seed=SPLIT_SEED, device=args.device,
                  inputs_sha256=input_digest, pair_sha256=approved['pair_sha256'],
                  representations_sha256=sha256(args.representations / 'manifest.json'),
                  order={key: value.tolist() for key, value in rows.items()},
                  max_epochs=args.smoke_epochs or MAX_EPOCHS, smoke_rows=args.smoke_rows,
                  patience=PATIENCE, validation_seed=VALIDATION_SEED, test_seed=TEST_SEED,
                  runtime=runtime, native_source_sha256=dose.SOURCE_DIGEST,
                  architecture=architecture,
                  environment=dict(environment(args.device), torch=str(torch.__version__)),
                  backend=dose.configure_backend(args.device))
    adapter = RepresentationAdapter(args.source, transform, config['width'])
    return ArmInputs(data, representations, args.arm), rows, config, adapter


def train(args, data, rows, config, adapter):
    start = time.monotonic()
    training, validation = rows['train'], rows['validation']
    model = adapter.build(TRAINING_SEED, args.device)
    initial = initialization(model)
    optimizer = adapter.optimizer(model)
    stopper = EarlyStopping(patience=PATIENCE, max_epochs=config['max_epochs'])
    shuffle = np.random.default_rng(TRAINING_SEED)
    history = []
    steps = 0

    def count_step(_optimizer, _args, _kwargs):
        nonlocal steps
        steps += 1

    with phase_directory(args.run, 'fit') as directory:
        optimizer.register_step_post_hook(count_step)
        target = data.expression[validation]
        for epoch in range(1, config['max_epochs'] + 1):
            loss = adapter.train_epoch(model, data, training, optimizer, shuffle)
            prediction = predict(adapter, model, data, validation, VALIDATION_SEED)
            mse = float(np.mean(np.square(prediction - target)))
            if not np.isfinite([loss, mse]).all():
                raise ValueError(f'nonfinite fit at epoch {epoch}')
            if stopper.observe(mse, epoch):
                checkpoint = {key: tensor.detach().to('cpu', copy=True)
                              for key, tensor in model.state_dict().items()}
                best_prediction = prediction
            history.append(dict(epoch=epoch, train_loss=loss, validation_mse=mse,
                                best_mse=stopper.best, stale=stopper.stale,
                                lr=optimizer.param_groups[0]['lr'], optimizer_steps=steps))
            if stopper.stopped:
                break

        if steps < 1 or initialization(model)['tensors'] == initial['tensors']:
            raise ValueError('training did not update model state')
        model.load_state_dict(checkpoint)
        restored = predict(adapter, model, data, validation, VALIDATION_SEED)
        if not np.array_equal(best_prediction, restored):
            raise ValueError('best checkpoint does not restore the selected validation prediction')
        state = dict(config=config, state_dict=checkpoint, best_epoch=stopper.best_epoch,
                     validation_mse=stopper.best)
        torch.save(state, directory / 'checkpoint.pt')
        np.save(directory / 'validation_prediction.npy', restored)
        pd.DataFrame(history).to_csv(directory / 'learning_curve.csv', index=False)
        return seal(directory, 'fit', dict(config=config, initialization=initial,
                    best_epoch=stopper.best_epoch, validation_mse=stopper.best, epochs=len(history),
                    optimizer_steps=steps, elapsed_seconds=time.monotonic() - start))


def load_checkpoint(args, config, adapter, fit):
    if fit['config'] != config:
        raise ValueError('current inputs, runtime or fold differ from the sealed fit')
    history = pd.read_csv(args.run / 'fit/learning_curve.csv', float_precision='round_trip')
    best = history.loc[history.validation_mse.idxmin()]
    if (history.epoch.tolist() != list(range(1, fit['epochs'] + 1))
            or not history.optimizer_steps.diff().iloc[1:].gt(0).all()
            or fit['best_epoch'] != int(best.epoch) or fit['validation_mse'] != float(best.validation_mse)
            or fit['optimizer_steps'] != int(history.optimizer_steps.iloc[-1])
            or len(history) != fit['epochs']):
        raise ValueError('fit selection or step count differs from the validation history')
    saved = torch.load(args.run / 'fit/checkpoint.pt', map_location='cpu', weights_only=True)
    if (saved['config'] != config or saved['best_epoch'] != fit['best_epoch']
            or saved['validation_mse'] != fit['validation_mse']):
        raise ValueError('checkpoint does not match the sealed validation selection')
    model = adapter.build(TRAINING_SEED, args.device)
    if initialization(model) != fit['initialization']:
        raise ValueError('reconstructed initial architecture or tensors differ')
    model.load_state_dict(saved['state_dict'], strict=True)
    return model


def check_replay_tokens(*receipts):
    tokens = {receipt['process_token'] for receipt in receipts}
    if len(tokens) != len(receipts) or PROCESS_TOKEN in tokens:
        raise ValueError('fit, validation replay, test evaluation and test replay require fresh processes')


def score_test(data, rows, prediction):
    return score_prediction(data.obs.iloc[rows], data.expression[rows], prediction,
                            data.control[rows], prediction_ids=data.obs.index[rows])


def advance(args, data, rows, config, adapter):
    fit = verify_phase(args.run / 'fit', 'fit')
    check_replay_tokens(fit)
    model = load_checkpoint(args, config, adapter, fit)
    fit_digest = sha256(args.run / 'fit/manifest.json')
    validation = rows['validation']
    if args.phase == 'replay-validation':
        actual = predict(adapter, model, data, validation, VALIDATION_SEED)
        expected = np.load(args.run / 'fit/validation_prediction.npy', mmap_mode='r')
        if (not np.array_equal(actual, expected)
                or float(np.mean(np.square(actual - data.expression[validation]))) != fit['validation_mse']):
            raise ValueError('fresh-process validation replay differs')
        with phase_directory(args.run, 'validation-replay') as directory:
            return seal(directory, 'validation-replay', dict(fit_sha256=fit_digest,
                        rows=len(validation), prediction_sha256=fit['outputs']['validation_prediction.npy'],
                        validation_mse=fit['validation_mse']))

    replay = verify_phase(args.run / 'validation-replay', 'validation-replay')
    check_replay_tokens(fit, replay)
    if (replay['fit_sha256'] != fit_digest
            or replay['prediction_sha256'] != fit['outputs']['validation_prediction.npy']
            or replay['rows'] != len(validation) or replay['validation_mse'] != fit['validation_mse']):
        raise ValueError('validation replay does not bind the selected checkpoint')
    validation_digest = sha256(args.run / 'validation-replay/manifest.json')
    test_rows = rows['test']
    evidence = dict(fit_sha256=fit_digest, validation_replay_sha256=validation_digest, rows=len(test_rows))
    if args.phase == 'evaluate':
        prediction = predict(adapter, model, data, test_rows, TEST_SEED)
        scores = score_test(data, test_rows, prediction)
        with phase_directory(args.run, 'test') as directory:
            np.save(directory / 'test_prediction.npy', prediction)
            scores.to_csv(directory / 'test_per_sample.csv', index=False)
            write_json(directory / 'test_metrics.json', summary(scores))
            return seal(directory, 'test', evidence)

    test = verify_phase(args.run / 'test', 'test')
    check_replay_tokens(fit, replay, test)
    if any(test[key] != value for key, value in evidence.items()):
        raise ValueError('test evidence refers to different fit, validation or rows')
    prediction = predict(adapter, model, data, test_rows, TEST_SEED)
    expected = np.load(args.run / 'test/test_prediction.npy', mmap_mode='r')
    scores = score_test(data, test_rows, prediction)
    if (not np.array_equal(prediction, expected)
            or scores.to_csv(index=False).encode() != (args.run / 'test/test_per_sample.csv').read_bytes()
            or summary(scores) != read_json(args.run / 'test/test_metrics.json')):
        raise ValueError('fresh-process test prediction or metric replay differs')
    with phase_directory(args.run, 'test-replay') as directory:
        return seal(directory, 'test-replay', dict(evidence,
                    test_sha256=sha256(args.run / 'test/manifest.json'),
                    prediction_sha256=test['outputs']['test_prediction.npy']))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('train', 'replay-validation', 'evaluate', 'replay-test'))
    for argument in ('inputs', 'representations', 'source', 'run'):
        parser.add_argument('--' + argument, type=Path, required=True)
    parser.add_argument('--arm', choices=ARMS, required=True)
    parser.add_argument('--fold', type=int, choices=range(3), required=True)
    parser.add_argument('--device', choices=('cpu', 'cuda'), required=True)
    parser.add_argument('--smoke-rows', type=int)
    parser.add_argument('--smoke-epochs', type=int)
    args = parser.parse_args()
    if ((args.smoke_rows is None) != (args.smoke_epochs is None)
            or (args.smoke_rows is not None and
                (not 3 <= args.smoke_rows <= 64 or not 1 <= args.smoke_epochs <= 5))):
        parser.error('bounded smoke requires 3–64 rows and 1–5 epochs')
    os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
    torch.use_deterministic_algorithms(True)
    torch.set_num_threads(2)
    data, rows, config, adapter = load_arm(args)
    receipt = (train if args.phase == 'train' else advance)(args, data, rows, config, adapter)
    print(args.arm, args.fold, receipt['phase'], 'PASS', flush=True)


if __name__ == '__main__':
    main()
