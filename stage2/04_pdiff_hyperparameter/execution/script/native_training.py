"""Random-init training, EMA checkpoint selection and native prediction replay."""
from importlib.util import find_spec
import math
from pathlib import Path
import random
import resource
import signal
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace


def verify_native_source():
    # The Step 03 adapter resolves its vendor tree relative to its own file.
    adapter = Path(find_spec('perturbdiff_v2_adapter').origin).resolve()
    checkout = next(path for path in adapter.parents if (path / '.git').exists())
    verifier = checkout / 'pert_model/PerturbDiff/script/verify.sh'
    try:
        subprocess.run(['bash', str(verifier)], check=True, stdout=subprocess.DEVNULL)
    except subprocess.CalledProcessError as error:
        raise ValueError('vendored PerturbDiff source checksum mismatch') from error


verify_native_source()

import numpy as np
import torch

from artifacts import read_json, sha256, write_json
from perturbdiff_v2_adapter import (
    EMAWeightAveraging, JointPerturbDiff, digest, get_cosine_schedule_with_warmup,
    get_optimizer, groups_for, save, validation_random_stream,
)
from predictor_inputs import PredictorInputs
from training_inputs import STEP03, json_hash
from training_policy import MONITOR, Selection, native_batch_mean, native_epoch_mean

FORMAT = 'step04-ema-prediction-v1'
POSITIVE_INTS = ('max_steps', 'validation_interval', 'minimum_steps', 'patience_steps')


def load_recipe(runtime, path, request):
    recipe = read_json(path)
    keys = {'stage', 'fold', 'seed', 'device', 'min_delta', 'max_training_seconds', *POSITIVE_INTS}
    if (set(recipe) != keys or recipe['stage'] not in ('smoke', 'pilot', 'formal')
            or recipe['fold'] != request['identity']['fold'] or recipe['seed'] not in (17, 42, 2026)
            or recipe['device'] != 'cuda'
            or any(type(recipe[key]) is not int or recipe[key] < 1 for key in POSITIVE_INTS)
            or type(recipe['min_delta']) not in (int, float) or not 0 <= recipe['min_delta'] < math.inf
            or type(recipe['max_training_seconds']) not in (int, float)
            or not 0 < recipe['max_training_seconds'] < math.inf):
        raise ValueError('invalid training recipe')
    fixed = read_json(runtime / STEP03 / 'execution/models/perturbdiff/v2.config.json')
    fixed_keys = ('batch_groups', 'population_size', 'chemical_dimension', 'learning_rate', 'weight_decay',
                  'warmup_steps', 'sampling_start_time', 'threads', 'expression_scale',
                  'schedule_training_steps', 'ema_decay', 'ema_update_steps')
    return {**{key: fixed[key] for key in fixed_keys}, **recipe}


def seed_randomness(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def objective(model, batch):
    losses = model.native._compute_loss(batch)
    return native_batch_mean(losses['loss'], losses['weights'])


def validate(model, inputs, groups, config, device, check_stop):
    values, sizes = [], []
    rng_seed = config['seed'] + 100000
    rng = np.random.default_rng(rng_seed)
    with torch.no_grad(), validation_random_stream(rng_seed, device):
        for offset in range(0, len(groups), config['batch_groups']):
            check_stop()
            part = groups[offset:offset + config['batch_groups']]
            batch, _ = model.batch(inputs, 'validation', part, device, rng, config['population_size'])
            values.append(objective(model, batch))
            sizes.append(len(part))
        loss = native_epoch_mean(values, torch.stack)
    return loss, [float(value) for value in values], sizes


def predict_probe(model, inputs, groups, config, device):
    """Exercise the native predictor with the first frozen validation group."""
    batch, indices = model.batch(inputs, 'validation', groups[:1], device)
    batch.pop('pert_emb', None)
    probe_seed = config['seed'] + 200000
    generator = torch.Generator(device=device).manual_seed(probe_seed)
    noise = torch.randn(batch['cont_emb'].shape, device=device, generator=generator)
    # DDIM draws from the global RNG even at eta=0. Isolate those draws from training.
    with torch.no_grad(), validation_random_stream(probe_seed, device):
        delta = (model.sample(batch, config, noise) - batch['cont_emb']) * config['expression_scale']
    if not torch.isfinite(delta).all():
        raise ValueError('nonfinite native prediction')
    array = delta.cpu().numpy().reshape(-1, inputs.n_genes)
    case_ids = [inputs.cases['validation'][int(index)]['case_id'] for index in indices.reshape(-1)]
    return array, case_ids


def checkpoint_reload(directory, metadata, model, inputs, groups, config, device, check_stop):
    path = directory / 'best.pt'
    check_stop()
    if sha256(path) != metadata['checkpoint_sha256']:
        raise ValueError('best checkpoint file hash mismatch')
    checkpoint = torch.load(path, map_location='cpu', weights_only=True)
    for key in ('identity', 'recipe_sha256', 'best_step', 'best_validation_loss', 'model_state_sha256'):
        if checkpoint[key] != metadata[key]:
            raise ValueError(f'best checkpoint metadata mismatch: {key}')
    if checkpoint['format'] != FORMAT:
        raise ValueError('unsupported checkpoint format')
    model.load_state_dict(checkpoint['model_state'], strict=True)
    if digest(model) != metadata['model_state_sha256']:
        raise ValueError('best checkpoint model state hash mismatch')
    model.eval()
    loss, _, _ = validate(model, inputs, groups, config, device, check_stop)
    predicted, case_ids = predict_probe(model, inputs, groups, config, device)
    reference = checkpoint['probe_delta'].numpy()
    if predicted.shape != reference.shape or case_ids != checkpoint['probe_case_ids']:
        raise ValueError('validation probe case/axis mismatch')
    error = float(np.max(np.abs(predicted - reference)))
    if error != 0 or loss != metadata['best_validation_loss']:
        raise ValueError(f'EMA reload differs from selection (prediction error {error}, loss {loss})')
    import hashlib
    return {'validation_loss': loss, 'prediction_sha256': hashlib.sha256(predicted.tobytes()).hexdigest(),
            'prediction_max_abs_error': error, 'case_ids': case_ids, 'prediction_shape': list(predicted.shape)}


def resources(device):
    return {'peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'peak_gpu_allocated_bytes': torch.cuda.max_memory_allocated(device),
            'peak_gpu_reserved_bytes': torch.cuda.max_memory_reserved(device)}


def fit(runtime, request, recipe_path, directory):
    config = load_recipe(runtime, recipe_path, request)
    directory.mkdir(parents=True)
    identity = {**request['identity'], 'seed': config['seed'], 'stage': config['stage']}
    metadata = {'identity': identity, 'recipe': config, 'recipe_sha256': json_hash(config),
                'monitor': MONITOR, 'status': 'running', 'history': [],
                'upstream_commit': '31be89cb8cce5dddef0782bac29d5c748ce55931',
                'environment': {'python': sys.version, 'torch': str(torch.__version__),
                                'numpy': np.__version__, 'cuda': torch.version.cuda}}
    write_json(directory / 'run.json', metadata)
    started = time.monotonic()
    deadline = started + config['max_training_seconds']
    stop_reason = None

    def handle_stop(signum, _frame):
        nonlocal stop_reason
        stop_reason = 'signal'

    def check_stop():
        nonlocal stop_reason
        if stop_reason == 'signal':
            raise InterruptedError('signal')
        if time.monotonic() >= deadline:
            stop_reason = 'training_budget'
            raise InterruptedError(stop_reason)

    signal.signal(signal.SIGTERM, handle_stop)
    signal.signal(signal.SIGINT, handle_stop)
    device = torch.device(config['device'])
    step, training_seconds = 0, 0.0
    try:
        torch.set_num_threads(config['threads'])
        seed_randomness(config['seed'])
        torch.cuda.reset_peak_memory_stats(device)
        inputs = PredictorInputs(request)
        groups = {role: groups_for(inputs.cases[role]) for role in ('train', 'validation')}
        ordered_group_ids = list(dict.fromkeys(case['group_id'] for case in inputs.cases['validation']))
        metadata['validation_batching'] = {
            'group_ids': ordered_group_ids, 'group_ids_sha256': json_hash(ordered_group_ids),
            'order': 'first appearance in frozen case CSV',
            'batch_groups': config['batch_groups'], 'population_size': config['population_size'],
            'random_seed': config['seed'] + 100000, 'remainder': 'keep and weight as one batch',
            'population_replacement': 'only if group has fewer wells than population_size',
        }
        model = JointPerturbDiff(inputs, config, directory).to(device)
        model.native.to(device)
        decay, no_decay = [], []
        for name, parameter in model.named_parameters():
            if parameter.requires_grad:
                (no_decay if 'bias' in name or 'norm' in name else decay).append(parameter)
        optimizer = get_optimizer([{'params': decay, 'weight_decay': config['weight_decay']},
                                   {'params': no_decay, 'weight_decay': 0.}], model.native.optimizer_cfg.optimizer)
        scheduler = get_cosine_schedule_with_warmup(optimizer, num_warmup_steps=config['warmup_steps'],
                                                     num_training_steps=config['schedule_training_steps'],
                                                     min_ratio=.1, plateau_ratio=.1)
        ema = EMAWeightAveraging(decay=config['ema_decay'], update_steps=config['ema_update_steps'])
        context = SimpleNamespace(global_step=0, current_epoch=0)
        ema.setup(context, model.native, 'fit')
        if model.native.device != next(ema._average_model.parameters()).device:
            raise ValueError('native EMA device mismatch')
        rng = np.random.default_rng(config['seed'])
        policy = Selection(config['minimum_steps'], config['patience_steps'], config['min_delta'])
        while step < config['max_steps'] and stop_reason is None:
            order = rng.permutation(len(groups['train']))
            for offset in range(0, len(order), config['batch_groups']):
                check_stop()
                step_started = time.monotonic()
                model.train()
                selected = [groups['train'][i] for i in order[offset:offset + config['batch_groups']]]
                batch, _ = model.batch(inputs, 'train', selected, device, rng, config['population_size'])
                optimizer.zero_grad(set_to_none=True)
                loss = objective(model, batch)
                if not torch.isfinite(loss):
                    raise ValueError('nonfinite native training objective')
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                optimizer.step()
                scheduler.step()
                step += 1
                context.global_step = step
                ema.on_train_batch_end(context, model.native, None, batch, offset // config['batch_groups'])
                if offset + config['batch_groups'] >= len(order):
                    ema.on_train_epoch_end(context, model.native)
                    context.current_epoch += 1
                torch.cuda.synchronize(device)
                training_seconds += time.monotonic() - step_started
                if step % config['validation_interval'] == 0 or step == config['max_steps']:
                    model.eval()
                    ema.on_validation_epoch_start(context, model.native)
                    try:
                        validation_start = time.monotonic()
                        value, losses, sizes = validate(model, inputs, groups['validation'], config, device, check_stop)
                        event = policy.observe(step, value)
                        event.update(training_loss=float(loss.detach()), batch_losses=losses, batch_sizes=sizes,
                                     training_seconds=training_seconds,
                                     validation_seconds=time.monotonic() - validation_start,
                                     checkpoint_seconds=0.)
                        if event['checkpoint_refreshed']:
                            checkpoint_started = time.monotonic()
                            predicted, case_ids = predict_probe(model, inputs, groups['validation'], config, device)
                            metadata.update(best_step=step, best_validation_loss=value, model_state_sha256=digest(model))
                            checkpoint = {key: metadata[key] for key in ('identity', 'recipe_sha256', 'best_step',
                                                                        'best_validation_loss', 'model_state_sha256')}
                            save(directory / 'best.pt', {**checkpoint, 'format': FORMAT, 'recipe': config,
                                                         'model_state': model.state_dict(),
                                                         'probe_case_ids': case_ids,
                                                         'probe_delta': torch.from_numpy(predicted)})
                            metadata['checkpoint_sha256'] = sha256(directory / 'best.pt')
                            event['checkpoint_seconds'] = time.monotonic() - checkpoint_started
                        metadata['history'].append(event)
                        write_json(directory / 'run.json', metadata)
                    finally:
                        ema.on_validation_epoch_end(context, model.native)
                    if event['early_stop']:
                        stop_reason = 'validation_early_stop'
                if step >= config['max_steps'] or stop_reason is not None:
                    break
        check_stop()
        del ema, optimizer, scheduler
        reload_started = time.monotonic()
        metadata['reload'] = checkpoint_reload(directory, metadata, model, inputs, groups['validation'],
                                               config, device, check_stop)
        check_stop()
        metadata['reload_seconds'] = time.monotonic() - reload_started
        metadata.update(status='completed', stop_reason=stop_reason or 'max_steps')
    except InterruptedError:
        metadata.update(status='incomplete', stop_reason=stop_reason)
        raise
    except Exception as error:
        metadata.update(status='failed', stop_reason='worker_error', error=str(error))
        raise
    finally:
        metadata.update(optimizer_steps=step, wall_seconds=time.monotonic() - started,
                        training_seconds=training_seconds)
        if torch.cuda.is_initialized():
            metadata['resources'] = resources(device)
        write_json(directory / 'run.json', metadata)


def reload(request, directory, output):
    metadata = read_json(directory / 'run.json')
    if (metadata['status'] != 'completed'
            or any(metadata['identity'][key] != value for key, value in request['identity'].items())
            or json_hash(metadata['recipe']) != metadata['recipe_sha256']):
        raise ValueError('completed run, request identity and recipe are required for replay')
    config = metadata['recipe']
    device = torch.device(config['device'])
    torch.set_num_threads(config['threads'])
    seed_randomness(config['seed'])
    inputs = PredictorInputs(request)
    groups = groups_for(inputs.cases['validation'])
    with tempfile.TemporaryDirectory() as temporary:
        model = JointPerturbDiff(inputs, config, Path(temporary)).to(device)
        model.native.to(device)
    evidence = checkpoint_reload(directory, metadata, model, inputs, groups, config, device, lambda: None)
    if evidence != metadata['reload']:
        raise ValueError('replay differs from completed run')
    write_json(output, {'identity': metadata['identity'], 'best_step': metadata['best_step'],
                        'checkpoint_sha256': metadata['checkpoint_sha256'], **evidence})
