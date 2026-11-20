"""Native TranSiGen fitting and exact independent validation replay."""
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import secrets
import time

import numpy as np
import pandas as pd
import torch

from . import STEP_03
from tsg_round1 import model as models
from tsg_round1.inputs import Inputs, sha256
from tsg_round1.learner import environment, runtime_identity as upstream_runtime
from tsg_round1.protocol import TRAINING_SEED
from tenv5k.learning import EarlyStopping, fixed_rng

METHOD = 'latent-embedding'
LEARNING_RATES = {'tahoe': 1e-4, 'tenv-only': 1e-5, 'transfer': 1e-5}
VALIDATION_SEED = 22617
PREDICTION_SEED = 100042
PROCESS_TOKEN = secrets.token_hex(32)
FIT_SCHEMA = 'stage2-tsg-transfer-fit-v1'


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(f'.{path.name}.{secrets.token_hex(8)}.tmp')
    try:
        temporary.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + '\n')
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def load_inputs(path, expected_source):
    """Read only the accepted ADR 0078 frozen source via the Step 03 API."""
    path = Path(path)
    pair_path = path.parent / 'manifest.json'
    evidence = read_json(STEP_03 / 'log/247-axis-local-smoke.json')
    if sha256(pair_path) != evidence['pair_sha256']:
        raise ValueError('input pair differs from accepted Step 03 identity')
    pair = read_json(pair_path)
    if sha256(path / 'manifest.json') != pair['sources'][expected_source]:
        raise ValueError('source manifest differs from the frozen input pair')
    data = Inputs(path)
    if data.source != expected_source:
        raise ValueError('source does not match the requested experiment')
    return data


def predict(native, model, data, rows, seed=PREDICTION_SEED):
    with fixed_rng(seed, next(model.parameters()).device):
        result = native.predict(model, data, rows)
    if result.dtype != np.float32 or not np.isfinite(result).all():
        raise ValueError('native prediction must be finite float32')
    return result


def _validation(native, model, data, rows):
    """Evaluate all validation rows in frozen order, keeping only one batch in RAM."""
    predictions = hashlib.sha256()
    squared_error = 0.
    values = 0
    size = native.batch_size
    with fixed_rng(VALIDATION_SEED, next(model.parameters()).device):
        for offset in range(0, len(rows), size):
            batch = rows[offset:offset + size]
            predicted = native.predict(model, data, batch)
            truth = data.expression[batch]
            if predicted.dtype != np.float32 or predicted.shape != truth.shape:
                raise ValueError('invalid native validation prediction shape or dtype')
            predictions.update(predicted.tobytes())
            squared_error += float(np.square(predicted.astype(np.float64) - truth).sum())
            values += truth.size
    mse = squared_error / values
    if not math.isfinite(mse):
        raise ValueError('nonfinite validation treated-expression MSE')
    return mse, predictions.hexdigest()


def _runtime():
    return dict(fitting=sha256(Path(__file__)), upstream=upstream_runtime())


def _configuration(data, train_rows, validation_rows, kind, fold, device, smoke_epochs, origin):
    train_rows = np.asarray(train_rows, dtype=int)
    validation_rows = np.asarray(validation_rows, dtype=int)
    for rows in (train_rows, validation_rows):
        if (rows.ndim != 1 or not len(rows) or np.unique(rows).size != len(rows)
                or (rows < 0).any() or (rows >= len(data.expression)).any()):
            raise ValueError('fit row identities are empty, duplicate or out of range')
    if len(train_rows) < 2 or np.intersect1d(train_rows, validation_rows).size:
        raise ValueError('invalid disjoint fitting partitions')
    if kind not in LEARNING_RATES or fold not in range(3):
        raise ValueError('unregistered fit identity')
    if kind == 'tahoe' and data.source != 'tahoe_bulk' or kind != 'tahoe' and data.source != 'tenv':
        raise ValueError('fit kind does not match input source')
    if smoke_epochs is not None and not 1 <= smoke_epochs < 1000:
        raise ValueError('smoke epoch cap must be between 1 and 999')
    if (kind == 'transfer') != bool(origin):
        raise ValueError('only transfer fits require a locked source origin')
    model_environment = environment(device)
    model_environment['torch'] = str(model_environment['torch'])
    transform = models.fit_transform(METHOD, data.doses(train_rows))
    return dict(
        schema=FIT_SCHEMA, source=data.source, kind=kind, fold=fold, device=str(device),
        seed=TRAINING_SEED, validation_seed=VALIDATION_SEED, origin=origin,
        train_rows=train_rows.tolist(), validation_rows=validation_rows.tolist(),
        smoke_epochs=smoke_epochs, transform=transform, max_epochs=smoke_epochs or 1000,
        learning_rate=LEARNING_RATES[kind], patience=50, batch_size=64,
        inputs_sha256=sha256(data.directory / 'manifest.json'),
        model_configuration=models.configuration(METHOD, transform),
        native_source_sha256=models.SOURCE_DIGEST, runtime=_runtime(),
        environment=model_environment,
    )


def _backend(device):
    os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
    torch.use_deterministic_algorithms(True)
    return models.configure_backend(device)


def train(data, source, candidate, *, train_rows, validation_rows, fold, kind, device,
          smoke_epochs=None, initialize=None, origin=None):
    """Fit one immutable attempt and restore the earliest strict validation minimum."""
    candidate = Path(candidate)
    train_rows = np.asarray(train_rows, dtype=int)
    validation_rows = np.asarray(validation_rows, dtype=int)
    config = _configuration(data, train_rows, validation_rows, kind, fold, device, smoke_epochs, origin)
    if (kind == 'transfer') != (initialize is not None):
        raise ValueError('only transfer fits require winner initialization')
    candidate.mkdir(parents=True)
    started = time.monotonic()
    config['backend'] = _backend(device)
    native = models.adapter(source, METHOD, config['transform'])
    model = native.build(TRAINING_SEED, device)
    initialization = initialize(model, config['transform']) if initialize is not None else None
    optimizer = native.optimizer(model)
    for group in optimizer.param_groups:
        group['lr'] = config['learning_rate']
    if optimizer.state:
        raise ValueError('fit requires fresh Adam state')
    initial_weight = next(model.parameters()).detach().cpu().clone()
    stopper = EarlyStopping(max_epochs=config['max_epochs'], patience=config['patience'])
    shuffle = np.random.default_rng(TRAINING_SEED)
    history = []
    checkpoint = candidate / 'checkpoint.pt'
    best_digest = None
    steps = math.ceil(len(train_rows) / native.batch_size) - (len(train_rows) % native.batch_size == 1)
    try:
        for epoch in range(1, config['max_epochs'] + 1):
            loss = native.train_epoch(model, data, train_rows, optimizer, shuffle)
            if not math.isfinite(loss):
                raise ValueError('nonfinite native loss')
            mse, digest = _validation(native, model, data, validation_rows)
            improved = stopper.observe(mse, epoch)
            history.append(dict(epoch=epoch, train_loss=loss, validation_mse=mse,
                                improved=improved, training_steps=epoch * steps))
            if improved:
                torch.save(dict(state_dict=model.state_dict(), config=config, best_epoch=epoch,
                                validation_mse=mse), checkpoint)
                best_digest = digest
            if stopper.stopped:
                break
        if torch.equal(initial_weight, next(model.parameters()).detach().cpu()):
            raise ValueError('training did not update the first native weight')
        saved = torch.load(checkpoint, map_location='cpu', weights_only=True)
        model.load_state_dict(saved['state_dict'])
        mse, digest = _validation(native, model, data, validation_rows)
        if (mse, digest) != (stopper.best, best_digest):
            raise ValueError('restored validation checkpoint differs from its best epoch')
        pd.DataFrame(history).to_csv(candidate / 'learning_curve.csv', index=False)
        receipt = dict(
            config, best_epoch=stopper.best_epoch, validation_mse=mse,
            validation_prediction_sha256=digest, epochs=len(history),
            training_steps=len(history) * steps, training_process_token=PROCESS_TOKEN,
            elapsed_seconds=time.monotonic() - started, peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            initialization=initialization, stop_reason='patience' if stopper.stale >= stopper.patience else 'max_epochs',
            files={name: sha256(candidate / name) for name in ('checkpoint.pt', 'learning_curve.csv')},
        )
        write_json(candidate / 'fit.json', receipt)
        return receipt
    except Exception as error:
        write_json(candidate / 'failure.json', dict(error=str(error), error_type=type(error).__name__))
        raise


def _verify_fit(data, candidate, train_rows, validation_rows):
    receipt = read_json(candidate / 'fit.json')
    expected = _configuration(data, train_rows, validation_rows, receipt['kind'], receipt['fold'],
                              receipt['device'], receipt['smoke_epochs'], receipt['origin'])
    if any(receipt.get(key) != value for key, value in expected.items()):
        raise ValueError('fit provenance or partition differs')
    expected['backend'] = receipt['backend']
    files = receipt['files']
    if set(files) != {'checkpoint.pt', 'learning_curve.csv'}:
        raise ValueError('fit artifact coverage differs')
    for name, digest in files.items():
        if sha256(candidate / name) != digest:
            raise ValueError(f'fit artifact digest differs: {name}')
    curve = pd.read_csv(candidate / 'learning_curve.csv', float_precision='round_trip')
    stopper = EarlyStopping(max_epochs=expected['max_epochs'], patience=expected['patience'])
    if not np.isfinite(curve[['train_loss', 'validation_mse']].to_numpy()).all():
        raise ValueError('nonfinite fit history')
    for epoch, row in enumerate(curve.itertuples(index=False), 1):
        if row.epoch != epoch or row.improved != stopper.observe(row.validation_mse, epoch):
            raise ValueError('fit history violates strict best checkpoint selection')
    if (not stopper.stopped or stopper.best_epoch != receipt['best_epoch']
            or stopper.best != receipt['validation_mse'] or len(curve) != receipt['epochs']
            or curve.training_steps.iloc[-1] != receipt['training_steps']):
        raise ValueError('fit receipt differs from learning history')
    return expected, receipt


def load_fit(data, source, candidate, *, train_rows, validation_rows):
    """Verify provenance before loading a float32 validation-best checkpoint."""
    candidate = Path(candidate)
    expected, receipt = _verify_fit(data, candidate, train_rows, validation_rows)
    if _backend(receipt['device']) != receipt['backend']:
        raise ValueError('backend differs from training')
    saved = torch.load(candidate / 'checkpoint.pt', map_location='cpu', weights_only=True)
    if (saved['config'] != expected or saved['best_epoch'] != receipt['best_epoch']
            or saved['validation_mse'] != receipt['validation_mse']):
        raise ValueError('checkpoint contract differs from fit receipt')
    if any(value.is_floating_point() and value.dtype != torch.float32
           for value in saved['state_dict'].values()):
        raise ValueError('checkpoint is not native float32')
    native = models.adapter(source, METHOD, receipt['transform'])
    model = native.build(TRAINING_SEED, receipt['device'])
    model.load_state_dict(saved['state_dict'])
    return native, model, receipt


def replay_fit(data, source, candidate, *, train_rows, validation_rows):
    candidate = Path(candidate)
    native, model, receipt = load_fit(data, source, candidate, train_rows=train_rows,
                                      validation_rows=validation_rows)
    if receipt['training_process_token'] == PROCESS_TOKEN:
        raise ValueError('validation replay requires a fresh process')
    mse, digest = _validation(native, model, data, np.asarray(validation_rows, dtype=int))
    if mse != receipt['validation_mse'] or digest != receipt['validation_prediction_sha256']:
        raise ValueError('independent validation replay differs')
    proof = dict(fit_sha256=sha256(candidate / 'fit.json'), validation_mse=mse,
                 validation_prediction_sha256=digest, replay_process_token=PROCESS_TOKEN)
    path = candidate / 'validation-replay.json'
    if path.exists():
        proof = read_json(path)
        verify_replay(data, candidate, train_rows=train_rows, validation_rows=validation_rows,
                      kind=receipt['kind'], fold=receipt['fold'], formal=receipt['smoke_epochs'] is None)
    else:
        write_json(path, proof)
    return proof


def verify_replay(data, candidate, *, train_rows, validation_rows, kind, fold, formal=True):
    candidate = Path(candidate)
    _, receipt = _verify_fit(data, candidate, train_rows, validation_rows)
    if receipt['kind'] != kind or receipt['fold'] != fold:
        raise ValueError('replayed fit identity differs')
    if formal and (receipt['smoke_epochs'] is not None or torch.device(receipt['device']).type != 'cuda'):
        raise ValueError('formal verification requires a full-budget CUDA fit')
    proof = read_json(candidate / 'validation-replay.json')
    if (proof['fit_sha256'] != sha256(candidate / 'fit.json')
            or not proof['replay_process_token']
            or proof['replay_process_token'] == receipt['training_process_token']
            or proof['validation_mse'] != receipt['validation_mse']
            or proof['validation_prediction_sha256'] != receipt['validation_prediction_sha256']):
        raise ValueError('independent validation replay proof differs')
    return receipt
