"""Shared candidate lifecycle: strict stopping, best restore and replay acceptance."""
import importlib
import json
import platform
import os
import time
import resource
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from .inputs import Inputs, sha256
from .attempt_identity import read_registration
from .learning import EarlyStopping
from .protocol import score_prediction


def load_adapter(name):
    if name not in ('prnet', 'transigen', 'mlp'): raise ValueError('unknown adapter')
    return importlib.import_module('tenv5k.model_' + name)


def replay(inputs, checkpoint, source, *, adapter, device=None):
    """Fresh-process-safe reload with verified input/source identity."""
    native = load_adapter(adapter)
    data = Inputs(inputs)
    saved = torch.load(checkpoint, map_location='cpu', weights_only=False)
    if saved['input_manifest_sha256'] != sha256(Path(inputs) / 'manifest.json'):
        raise ValueError('checkpoint input identity mismatch')
    if saved['config']['native_source_sha256'] != native.SOURCE_DIGEST:
        raise ValueError('checkpoint source identity mismatch')
    config = saved['config']
    if saved.get('execution_binding') != config.get('execution_binding'):
        raise ValueError('checkpoint attempt registration mismatch')
    if config['adapter'] != adapter:
        raise ValueError('checkpoint adapter identity mismatch')
    for key in ('method_identity', 'compute_dtype', 'initialization_dtype', 'prediction_dtype', 'noise_dtype'):
        if key in native.CONFIG and config.get(key) != native.CONFIG[key]:
            raise ValueError('checkpoint numerical method identity mismatch: ' + key)
    device = device or config['device']
    if str(device) != config['device']:
        raise ValueError('exact replay requires recorded device')
    os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
    torch.use_deterministic_algorithms(True)
    if hasattr(native, 'configure_backend'):
        backend = native.configure_backend(device)
        if str(device).startswith('cuda') and config.get('backend') != backend:
            raise ValueError('checkpoint CUDA backend identity is missing or differs; do not relabel old runs')
    if config.get('compute_dtype') is not None:
        expected_dtype = getattr(torch, config['compute_dtype'])
        if any(v.is_floating_point() and v.dtype != expected_dtype for v in saved['state_dict'].values()):
            raise ValueError('checkpoint state precision differs from numerical method identity')
    model = native.build_model(source, config['seed'], device)
    model.load_state_dict(saved['state_dict'])
    rows = data.rows(config['scenario'], config['fold'], 'test')
    if config['smoke_rows'] is not None: rows = rows[:config['smoke_rows']]
    return serialize_prediction(native.sampled_prediction(model, data, rows, seed=config['inference_seed'], device=device), config)


def serialize_prediction(prediction, config):
    """Apply the declared storage dtype and check the stored values themselves."""
    dtype = config.get('prediction_dtype')
    with np.errstate(over='ignore', invalid='ignore'):
        result = np.asarray(prediction, dtype=dtype)
    if not np.isfinite(result).all():
        raise ValueError('non-finite prediction after serialization conversion')
    return result


def train(inputs, candidate, source, *, adapter, scenario, fold, seed, device='cpu',
          smoke_rows=None, max_epochs=1000, registration=None):
    """Build an isolated candidate. Smoke identities cannot be formal runs."""
    native = load_adapter(adapter)
    candidate = Path(candidate)
    if candidate.exists() and any(candidate.iterdir()):
        raise ValueError('candidate must be empty')
    candidate.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    try:
        if scenario not in ('cold-drug', 'warm-start') or fold not in range(3) or seed not in (17, 42, 2026):
            raise ValueError('unregistered run identity')
        if smoke_rows is None and max_epochs != 1000:
            raise ValueError('formal runs require max_epochs=1000')
        data = Inputs(inputs)
        train_rows, valid_rows, test_rows = [data.rows(scenario, fold, p) for p in ('train', 'validation', 'test')]
        if smoke_rows is not None:
            if smoke_rows < 2: raise ValueError('smoke needs at least two rows')
            train_rows, valid_rows, test_rows = [r[:smoke_rows] for r in (train_rows, valid_rows, test_rows)]
        if min(map(len, (train_rows, valid_rows, test_rows))) < 1:
            raise ValueError('empty split')
        # Author execution settings must be applied and recorded before native
        # construction in both fitting and the fresh-process replay entry point.
        backend = native.configure_backend(device) if hasattr(native, 'configure_backend') else None
        config = {**native.CONFIG, 'adapter': adapter, 'scenario': scenario, 'fold': fold, 'seed': seed,
                  'device': str(device), 'validation_seed': data.config['early_stopping']['validation_seed'],
                  'inference_seed': seed + 100000, 'smoke_rows': smoke_rows,
                  'max_epochs': max_epochs, 'patience': 50, 'strict_lower': True,
                  'native_source_sha256': native.SOURCE_DIGEST,
                  'validation_order': valid_rows.tolist(), 'test_order': test_rows.tolist(),
                  'singleton_train_batch': native.CONFIG.get('singleton_train_batch', 'merge into preceding batch (BatchNorm needs >=2)'),
                  'adaptations': native.CONFIG.get('adaptations', ['shared split/input readers', 'fixed-validation treated MSE early stopping'])}
        binding = read_registration(registration, adapter=adapter, scenario=scenario, fold=fold, seed=seed, device=device)
        if binding is not None:
            config['execution_binding'] = binding
        if backend is not None:
            config['backend'] = backend
        (candidate / 'config.json').write_text(json.dumps(config, indent=2) + '\n')
        os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
        torch.use_deterministic_algorithms(True)
        model = native.build_model(source, seed, device)
        optimizer = native.make_optimizer(model)
        scheduler = native.make_scheduler(optimizer)
        stopper = EarlyStopping(max_epochs=max_epochs)
        shuffle = np.random.default_rng(seed)
        history = []
        initial = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
        for epoch in range(1, max_epochs + 1):
            loss_value = native.train_epoch(model, data, train_rows, optimizer, shuffle, device)
            validation = native.sampled_prediction(model, data, valid_rows, seed=config['validation_seed'], device=device)
            mse = float(np.mean(np.square(validation.astype(np.float32) - np.asarray(data.expression[valid_rows], dtype=np.float32))))
            improved = stopper.observe(mse, epoch)
            if scheduler is not None: scheduler.step(mse)
            history.append(dict(epoch=epoch, train_loss=loss_value, validation_mse=mse,
                                best_mse=stopper.best, stale=stopper.stale, lr=optimizer.param_groups[0]['lr']))
            if improved:
                torch.save(dict(state_dict=model.state_dict(), config=config, execution_binding=binding, best_epoch=epoch,
                                validation_mse=mse, input_manifest_sha256=sha256(Path(inputs) / 'manifest.json')),
                           candidate / 'checkpoint.pt')
            if stopper.stopped: break
        changed = any(not torch.equal(initial[k], v.detach().cpu()) for k, v in model.state_dict().items() if k.endswith('weight'))
        if not changed: raise ValueError('training did not update weights')
        best = torch.load(candidate / 'checkpoint.pt', map_location=device, weights_only=False)
        model.load_state_dict(best['state_dict'])
        prediction = native.sampled_prediction(model, data, test_rows, seed=config['inference_seed'], device=device)
        prediction = serialize_prediction(prediction, config)
        np.save(candidate / 'prediction.npy', prediction)
        data.obs.iloc[test_rows][['profile_id', 'drug_parent_inchikey']].to_csv(candidate / 'identity.csv', index=False)
        wells, parents = score_prediction(data.obs.iloc[test_rows], data.expression[test_rows], prediction,
                                         data.control[test_rows], prediction_ids=data.obs.index[test_rows])
        wells.to_csv(candidate / 'per_well.csv', index=False)
        parents.to_csv(candidate / 'per_drug.csv', index=False)
        pd.DataFrame(history).to_csv(candidate / 'learning_curve.csv', index=False)
        environment = dict(python=platform.python_version(), torch=torch.__version__, numpy=np.__version__,
                           device=str(device), cuda=torch.version.cuda)
        (candidate / 'environment.json').write_text(json.dumps(environment, indent=2) + '\n')
        terminal = dict(status='trained-awaiting-fresh-process-replay', smoke=smoke_rows is not None,
                        best_epoch=stopper.best_epoch, epochs=len(history), validation_mse=stopper.best,
                        elapsed_seconds=time.monotonic() - started, peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                        weights_updated=changed)
        (candidate / 'terminal.json').write_text(json.dumps(terminal, indent=2) + '\n')
        return terminal
    except Exception as error:
        (candidate / 'terminal.json').write_text(json.dumps(dict(status='failed', error=str(error)), indent=2) + '\n')
        raise


def accept_replay(inputs, candidate, source, *, adapter, registration=None):
    """No training: exact recorded-environment replay before candidate completion."""
    candidate = Path(candidate)
    terminal = json.loads((candidate / 'terminal.json').read_text())
    if terminal['status'] not in ('trained-awaiting-fresh-process-replay', 'complete'):
        raise ValueError('cannot accept failed candidate')
    config = json.loads((candidate / 'config.json').read_text())
    saved = torch.load(candidate / 'checkpoint.pt', map_location='cpu', weights_only=False)
    binding = read_registration(registration, adapter=adapter, scenario=config['scenario'], fold=config['fold'],
                                seed=config['seed'], device=config['device'])
    if binding is not None and config.get('execution_binding') != binding:
        raise ValueError('replay attempt registration mismatch')
    if saved.get('execution_binding') != config.get('execution_binding'):
        raise ValueError('checkpoint attempt registration mismatch')
    if config != saved['config']:
        raise ValueError('candidate/checkpoint configuration mismatch')
    environment = json.loads((candidate / 'environment.json').read_text())
    if environment['torch'] != torch.__version__ or environment['numpy'] != np.__version__ or environment['python'] != platform.python_version():
        raise ValueError('exact replay requires recorded software environment')
    prediction = replay(inputs, candidate / 'checkpoint.pt', source, adapter=adapter)
    stored = np.load(candidate / 'prediction.npy')
    if prediction.dtype != stored.dtype or not np.array_equal(prediction, stored):
        raise ValueError('fresh process prediction replay differs')
    terminal.update(status='complete', replay='exact native sampled prediction')
    (candidate / 'terminal.json').write_text(json.dumps(terminal, indent=2) + '\n')
    outputs = {str(p.relative_to(candidate)): sha256(p) for p in sorted(candidate.iterdir()) if p.is_file() and p.name != 'manifest.json'}
    (candidate / 'manifest.json').write_text(json.dumps(dict(schema='tenv5k-learner-v1',
        input_manifest_sha256=sha256(Path(inputs) / 'manifest.json'), execution_binding=config.get('execution_binding'), outputs=outputs), indent=2) + '\n')
    return terminal
