"""Fit one registered coordinate, then replay and score in a fresh process."""

import argparse
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import torch

from run import METHODS, STEP04, digest, load_json, save_json

sys.path.insert(0, str(STEP04 / 'output/api'))
from tsg_kpgt.inputs import Inputs
from tsg_kpgt.learner import environment, runtime_identity
from tsg_kpgt.model import configure_backend, fit_transform
from tsg_round1.learner import predict, smoke_subset
from tsg_round1.protocol import score_prediction, summary
from tenv5k.learning import EarlyStopping

from model import adapter, dimensions

VALIDATION_SEED = 22617
INFERENCE_SEED_OFFSET = 100000


def runtime(device):
    # TorchVersion subclasses str, but checkpoints must stay weights-only-loadable.
    environment_details = dict(environment(device), torch=str(torch.__version__))
    return dict(environment=environment_details, step04_runtime=runtime_identity())


def validation_mse(prediction, data, rows):
    return float(np.square(prediction - data.expression[rows]).mean())


def train(native, model, data, config, root):
    seed = config['origin']['seed']
    optimizer = native.optimizer(model)
    max_epochs = config['identity']['max_epochs']
    stopper = EarlyStopping(max_epochs=max_epochs)
    shuffle = np.random.default_rng(seed)
    training = np.asarray(config['train_order'])
    validation = np.asarray(config['validation_order'])
    history = []
    for epoch in range(1, max_epochs + 1):
        loss = native.train_epoch(model, data, training, optimizer, shuffle)
        prediction = predict(native, model, data, validation, VALIDATION_SEED)
        mse = validation_mse(prediction, data, validation)
        if not np.isfinite([loss, mse]).all():
            raise ValueError(f'Non-finite training evidence at epoch {epoch}')
        if stopper.observe(mse, epoch):
            torch.save(dict(config=config, state_dict=model.state_dict(),
                            best_epoch=epoch, validation_mse=mse), root / 'checkpoint.pt')
            np.save(root / 'validation_prediction.npy', prediction)
        history.append(dict(epoch=epoch, train_loss=loss, validation_mse=mse,
                            best_mse=stopper.best, stale=stopper.stale,
                            lr=optimizer.param_groups[0]['lr']))
        if stopper.stopped:
            break
    pd.DataFrame(history).to_csv(root / 'learning_curve.csv', index=False)
    save_json(root / 'terminal.json', dict(
        status='trained', training_pid=os.getpid(), epochs=len(history),
        validation_mse=stopper.best, best_epoch=stopper.best_epoch,
        parameter_count=sum(parameter.numel() for parameter in model.parameters()),
    ))


def evaluate(native, model, data, config, root):
    terminal = load_json(root / 'terminal.json')
    if terminal['status'] != 'trained' or terminal['training_pid'] == os.getpid():
        raise ValueError('Test scoring requires a checkpoint trained in a separate process')
    checkpoint = torch.load(root / 'checkpoint.pt', map_location='cpu', weights_only=True)
    if checkpoint['config'] != config or any(checkpoint[key] != terminal[key]
                                            for key in ('best_epoch', 'validation_mse')):
        raise ValueError('Validation checkpoint identity changed')
    model.load_state_dict(checkpoint['state_dict'])
    del checkpoint
    validation = predict(native, model, data, np.asarray(config['validation_order']), VALIDATION_SEED)
    if not np.array_equal(validation, np.load(root / 'validation_prediction.npy')):
        raise ValueError('Best checkpoint did not exactly replay the validation predictions')
    origin = config['origin']
    test = smoke_subset(data, data.rows(origin['scenario'], origin['fold'], 'test'),
                        config['identity']['smoke_rows'])
    prediction = predict(native, model, data, test, origin['seed'] + INFERENCE_SEED_OFFSET)
    samples = score_prediction(data.obs.iloc[test], data.expression[test], prediction,
                               data.control[test], prediction_ids=data.obs.index[test])
    samples.to_csv(root / 'test_per_sample.csv', index=False)
    save_json(root / 'test_metrics.json', summary(samples))
    save_json(root / 'terminal.json', dict(terminal, status='complete',
              validation_replay='exact', evaluation_pid=os.getpid(), test_samples=len(samples)))


def execute(inputs, source, root, phase):
    config = load_json(root / 'config.json')
    origin, identity = config['origin'], config['identity']
    if digest(inputs / 'manifest.json') != identity['input_manifest_sha256']:
        raise ValueError('Input manifest changed')
    data = Inputs(inputs)
    fold, scenario = origin['fold'], origin['scenario']
    training, validation = (smoke_subset(data, data.rows(scenario, fold, partition), identity['smoke_rows'])
                            for partition in ('train', 'validation'))
    if len(training) < 2 or not len(validation):
        raise ValueError('Insufficient fold rows for fitting or evaluation')
    method, width, device = origin['method'], origin['width'], identity['device']
    if runtime(device) != identity['runtime']:
        raise ValueError('Registered runtime or environment changed')
    transform = fit_transform(METHODS[method], data.doses(training))
    os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
    torch.use_deterministic_algorithms(True)
    details = dict(
        transform=transform, dimensions=dimensions(width), train_order=training.tolist(),
        validation_order=validation.tolist(), backend=configure_backend(device),
        validation_seed=VALIDATION_SEED, inference_seed=origin['seed'] + INFERENCE_SEED_OFFSET,
    )
    config_path = root / 'config.json'
    if phase == 'train':
        config.update(details)
        save_json(config_path, config)
    elif any(config[key] != value for key, value in details.items()):
        raise ValueError('Training runtime, partition or dose transform changed')
    native = adapter(source, width, method, transform)
    model = native.build(origin['seed'], device)
    function = train if phase == 'train' else evaluate
    function(native, model, data, config, root)


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--runtime':
        print(json.dumps(runtime(sys.argv[2])))
        sys.exit(0)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', type=Path)
    parser.add_argument('source', type=Path)
    parser.add_argument('root', type=Path)
    parser.add_argument('phase', choices=('train', 'test'))
    execute(**vars(parser.parse_args()))
