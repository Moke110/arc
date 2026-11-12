"""Step-owned PRnet fit/predict command for the tenv comparison (#163).

The upstream predictor is imported unchanged from
``pert_model/PRnet/source/models/PRnet.py`` (class ``PGM``); the source tree is
never edited.  This command owns everything the retired adapter owned and keeps
it auditable in the step:

* a jointly trained ``Linear(2048, 1024)`` Morgan-2048 adapter, multiplied by
  ``log10(dose_um + 1)``, feeds the native ``PGM``;
* the objective is the author's Gaussian NLL on the expression target
  (``delta + same-well control``), with the native positive mean and
  ``softplus`` variance;
* training wells only drive the fit; validation wells only drive the
  epoch-based early stopping and checkpoint selection of ADR 0063;
* the selected checkpoint is reloaded from disk and the published prediction
  comes from that reloaded artifact, never from in-memory training state;
* a completed run publishes ``prediction.npy`` (all genes, finite), case
  identities, checkpoint, effective config, an identity record of hashes and
  the measured resource peak; a failed, interrupted or over-budget run writes
  an honest ``terminal.json`` and no prediction.

Run with the pinned isolated interpreter, never the shared 3.13 environment:

    /home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python \
        execution/script/prnet_fit_predict.py --request REQUEST --config CONFIG
"""

import os
import sys

sys.dont_write_bytecode = True

for _name in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
              'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_name, '2')

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import random
import resource
import signal
import subprocess
import tempfile
import time
import traceback

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F

# The pinned PRnet environment has numba 0.60 + scanpy 1.9.8 and cannot write a
# numba cache beside the package; without a cache directory ``import scanpy``
# (pulled in by the upstream module) fails while defining decorating helpers.
os.environ.setdefault(
    'NUMBA_CACHE_DIR', str(Path(tempfile.gettempdir()) / 'arc-stage2-prnet-numba'))
os.environ.setdefault(
    'MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'arc-stage2-prnet-matplotlib'))


def repo_root():
    for candidate in Path(__file__).resolve().parents:
        if (candidate / '.git').exists():
            return candidate
    raise RuntimeError('repository root not found')


def upstream_source():
    """The untouched upstream PRnet tree, from this checkout or the primary one.

    ADR 0030 keeps external ``source/`` trees out of git, so a linked worktree
    does not carry ``pert_model/PRnet/source``.  The primary checkout that owns
    the worktree's git directory does, and it is the same pinned tree.
    """
    local = REPO / 'pert_model/PRnet/source'
    if (local / 'models/PRnet.py').is_file():
        return local
    common = subprocess.run(
        ['git', '-C', str(REPO), 'rev-parse', '--git-common-dir'],
        capture_output=True, text=True, check=True).stdout.strip()
    common_dir = Path(common)
    if not common_dir.is_absolute():
        common_dir = REPO / common_dir
    candidate = common_dir.resolve().parent / 'pert_model/PRnet/source'
    if (candidate / 'models/PRnet.py').is_file():
        return candidate
    raise RuntimeError(f'upstream PRnet source not found beside {REPO}')


REPO = repo_root()
sys.path.insert(0, str(Path(__file__).resolve().parent))
from artifacts import read_json, read_rows, sha256, write_json  # noqa: E402
from early_stopping import STOP_REASONS, EarlyStopper  # noqa: E402

UPSTREAM = upstream_source()
sys.path.insert(0, str(UPSTREAM))
from models.PRnet import PGM  # noqa: E402

SCRIPT_DIR = Path(__file__).resolve().parent

CASE_FIELDS = ['case_id', 'drug_key', 'dose_um', 'group_id', 'control_group',
               'source_row', 'morgan_index']
ROLES = ('train', 'validation', 'test')
SCHEMA = 'stage2-02-models-tenv-inputs-v1'
SELECTION = ('minimum author Gaussian NLL on validation; '
             'fixed validation noise stream')
OUTPUT = 'native conditional Gaussian mean minus approved same-well control'


class RunFailure(RuntimeError):
    """A declared run precondition or behavioral gate failed."""


class OverBudget(RunFailure):
    """The run hit a declared execution ceiling before completion."""


class Interrupted(RunFailure):
    """The run received a termination request before completion."""


_STOP = False


def _request_stop(signum, frame):
    global _STOP
    _STOP = True


def check_stop():
    if _STOP:
        raise Interrupted('termination request received; no completed prediction')


def config_digest(config):
    return hashlib.sha256(
        json.dumps(config, sort_keys=True).encode()).hexdigest()


def _integer(config, name, low, high):
    value = config[name]
    if isinstance(value, bool) or not isinstance(value, int):
        raise RunFailure(f'config {name} must be an integer, found {value!r}')
    if not low <= value <= high:
        raise RunFailure(
            f'config {name} must be within [{low}, {high}], found {value}')
    return value


def _number(config, name, low, high, *, low_open=False, high_open=False):
    value = config[name]
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise RunFailure(f'config {name} must be a number, found {value!r}')
    value = float(value)
    if not np.isfinite(value):
        raise RunFailure(f'config {name} must be finite, found {value!r}')
    too_low = value <= low if low_open else value < low
    too_high = value >= high if high_open else value > high
    if too_low or too_high:
        raise RunFailure(
            f'config {name} is outside its declared range, found {value}')
    return value


def _layer_sizes(config, name):
    value = config[name]
    if (not isinstance(value, list) or not value
            or any(isinstance(item, bool) or not isinstance(item, int) or item < 1
                   for item in value)):
        raise RunFailure(f'config {name} must be a non-empty list of positive ints')
    return value


def load_config(path):
    config = read_json(path)
    required = ('max_epochs', 'batch_size', 'hidden_layer_sizes', 'z_dimension',
                'adaptor_layer_sizes', 'comb_dimension', 'drug_dimension',
                'dropout', 'noise_dimension', 'learning_rate', 'weight_decay',
                'scheduler_factor', 'scheduler_patience',
                'minimum_learning_rate', 'no_strict_decrease_epochs',
                'above_running_min_epochs', 'early_stopping',
                'max_training_seconds', 'inference_max_seconds', 'device',
                'threads')
    missing = [name for name in required if name not in config]
    if missing:
        raise RunFailure(f'config missing keys: {missing}')
    # ADR 0063: at most 100 epochs, both counters declared, 2 threads, 10 GiB.
    _integer(config, 'max_epochs', 1, 100)
    _integer(config, 'batch_size', 2, 4096)
    _integer(config, 'threads', 1, 2)
    _integer(config, 'z_dimension', 1, 4096)
    _integer(config, 'comb_dimension', 1, 4096)
    _integer(config, 'noise_dimension', 1, 4096)
    if _integer(config, 'drug_dimension', 1, 4096) != 1024:
        raise RunFailure('ADR 0063 adapter requires drug_dimension 1024')
    _layer_sizes(config, 'hidden_layer_sizes')
    _layer_sizes(config, 'adaptor_layer_sizes')
    _number(config, 'dropout', 0.0, 1.0, high_open=True)
    _number(config, 'learning_rate', 0.0, 1.0, low_open=True)
    _number(config, 'weight_decay', 0.0, 1.0)
    _number(config, 'scheduler_factor', 0.0, 1.0, low_open=True, high_open=True)
    _integer(config, 'scheduler_patience', 1, 1000)
    minimum = _number(config, 'minimum_learning_rate', 0.0, 1.0)
    if minimum >= float(config['learning_rate']):
        raise RunFailure('config minimum_learning_rate must be below learning_rate')
    _integer(config, 'no_strict_decrease_epochs', 1, 1000)
    _integer(config, 'above_running_min_epochs', 1, 1000)
    _number(config, 'max_training_seconds', 0.0, 86400.0)
    _number(config, 'inference_max_seconds', 0.0, 86400.0)
    early = config['early_stopping']
    expected = {
        'no_strict_decrease_epochs': config['no_strict_decrease_epochs'],
        'above_running_min_epochs': config['above_running_min_epochs'],
        'max_epochs': config['max_epochs'],
    }
    if early != expected:
        raise RunFailure(
            f'config early_stopping {early} disagrees with top-level recipe {expected}')
    if config['device'] not in ('auto', 'cpu', 'cuda'):
        raise RunFailure('config device must be auto, cpu or cuda')
    return config


def _array(path, shape, name):
    array = np.load(path, mmap_mode='r', allow_pickle=False)
    if array.dtype != np.float32:
        raise RunFailure(f'{name}: expected float32, found {array.dtype}')
    if array.shape != shape:
        raise RunFailure(f'{name}: expected shape {shape}, found {array.shape}')
    if not np.isfinite(np.asarray(array)).all():
        raise RunFailure(f'{name}: non-finite input value')
    return array


class Release:
    """The step-owned predictor input release, validated through its manifest."""

    def __init__(self, request):
        for key in ('identity', 'release', 'n_genes', 'target_kind',
                    'output_directory'):
            if key not in request:
                raise RunFailure(f'request missing key: {key}')
        if request['target_kind'] != 'delta':
            raise RunFailure('PRnet requires the frozen delta target contract')
        self.root = Path(request['release'])
        manifest_path = self.root / 'manifest.json'
        if not manifest_path.is_file():
            raise RunFailure(f'missing input manifest: {manifest_path}')
        self.manifest_path = manifest_path
        self.manifest = read_json(manifest_path)
        manifest = self.manifest
        if manifest.get('schema') != SCHEMA:
            raise RunFailure(f'unsupported input schema: {manifest.get("schema")}')
        if manifest.get('target_kind') != 'delta':
            raise RunFailure('input release target kind is not delta')
        if manifest.get('n_genes') != request['n_genes']:
            raise RunFailure(
                f'n_genes mismatch: request {request["n_genes"]} '
                f'manifest {manifest.get("n_genes")}')
        self.n_genes = int(request['n_genes'])
        if self.n_genes < 1:
            raise RunFailure('n_genes must be positive')
        expected_identity = {'model': request['identity'].get('model')}
        expected_identity.update(manifest.get('identity') or {})
        if request['identity'] != expected_identity:
            raise RunFailure('request identity differs from the release identity')
        self.identity = request['identity']
        seed = self.identity.get('seed')
        if isinstance(seed, bool) or not isinstance(seed, int):
            raise RunFailure('release identity seed must be an integer')

        files = manifest.get('files') or {}
        digests = manifest.get('sha256') or {}
        required = ['drugs', 'morgan']
        for role in ROLES:
            required += [f'{role}_cases', f'{role}_baseline']
        for role in ('train', 'validation'):
            required.append(f'{role}_delta')
        for name in required:
            if name not in files:
                raise RunFailure(f'input manifest missing file role: {name}')
            path = self.root / files[name]
            if not path.resolve().is_relative_to(self.root.resolve()):
                raise RunFailure(f'input path escapes the release: {name}')
            if not path.is_file():
                raise RunFailure(f'missing input file: {name}')
            if digests.get(name) != sha256(path):
                raise RunFailure(f'{name}: input digest mismatch')
        self.files = {name: self.root / files[name] for name in required}

        self.cases = {role: self._cases(role) for role in ROLES}
        self._reject_drug_leakage()
        self._drugs()
        self.arrays = {}
        for role in ROLES:
            self.arrays[f'{role}_baseline'] = _array(
                self.files[f'{role}_baseline'],
                (len(self.cases[role]), self.n_genes), f'{role}_baseline')
        for role in ('train', 'validation'):
            self.arrays[f'{role}_delta'] = _array(
                self.files[f'{role}_delta'],
                (len(self.cases[role]), self.n_genes), f'{role}_delta')

    def _cases(self, role):
        rows = read_rows(self.files[f'{role}_cases'])
        if not rows:
            raise RunFailure(f'empty {role} case table')
        for row in rows:
            if set(row) != set(CASE_FIELDS):
                raise RunFailure(f'{role} case schema differs from the allowlist')
            if not row['case_id'] or not row['drug_key'] or not row['group_id']:
                raise RunFailure(f'{role} case identity is incomplete')
            try:
                dose = float(row['dose_um'])
                source_row = int(row['source_row'])
                morgan_index = int(row['morgan_index'])
            except ValueError as error:
                raise RunFailure(f'{role} case has a malformed field: {error}') from error
            if not np.isfinite(dose) or dose <= 0:
                raise RunFailure(f'{role} case has an invalid treated dose')
            if source_row < 0 or morgan_index < 0:
                raise RunFailure(f'{role} case has a negative identity field')
            row['dose_um'] = dose
            row['source_row'] = source_row
            row['morgan_index'] = morgan_index
        if len({row['case_id'] for row in rows}) != len(rows):
            raise RunFailure(f'{role} case ids are not unique')
        return rows

    def _reject_drug_leakage(self):
        for left, right in (('train', 'validation'), ('train', 'test'),
                            ('validation', 'test')):
            overlap = ({row['drug_key'] for row in self.cases[left]}
                       & {row['drug_key'] for row in self.cases[right]})
            if overlap:
                raise RunFailure(
                    f'drug leakage between {left} and {right}: {sorted(overlap)[:3]}')

    def _drugs(self):
        rows = read_rows(self.files['drugs'])
        if not rows:
            raise RunFailure('empty drug table')
        morgan = np.load(self.files['morgan'], mmap_mode='r', allow_pickle=False)
        if morgan.dtype != np.uint8:
            raise RunFailure(f'morgan: expected uint8, found {morgan.dtype}')
        if morgan.ndim != 2 or morgan.shape[1] != 2048:
            raise RunFailure(f'morgan: expected (*, 2048), found {morgan.shape}')
        if len(rows) != morgan.shape[0]:
            raise RunFailure('drug table and Morgan row count differ')
        index = {}
        for row in rows:
            if set(row) != {'drug_key', 'morgan_index'}:
                raise RunFailure('drug table schema differs from the allowlist')
            key = row['drug_key']
            position = int(row['morgan_index'])
            if key in index or position in index.values():
                raise RunFailure('duplicate drug or Morgan identity')
            if position < 0 or position >= morgan.shape[0]:
                raise RunFailure('Morgan identity out of range')
            index[key] = position
        if sorted(index.values()) != list(range(morgan.shape[0])):
            raise RunFailure('Morgan indices are not a dense row prefix')
        if not np.isin(np.asarray(morgan), (0, 1)).all():
            raise RunFailure('Morgan fingerprints must be binary')
        self.drug_index = index
        self.morgan = morgan
        for role in ROLES:
            for row in self.cases[role]:
                if row['drug_key'] not in index:
                    raise RunFailure(f'unknown drug identity: {row["drug_key"]}')
                if row['morgan_index'] != index[row['drug_key']]:
                    raise RunFailure('case morgan_index disagrees with the drug table')

    def batch(self, role, indices, device):
        cases = [self.cases[role][int(index)] for index in indices]
        baseline = np.asarray(
            self.arrays[f'{role}_baseline'][indices], dtype='float32')
        fingerprint = np.asarray(
            self.morgan[[row['morgan_index'] for row in cases]], dtype='float32')
        dose = np.asarray([row['dose_um'] for row in cases], dtype='float32')
        tensors = [baseline, fingerprint, dose]
        if role in ('train', 'validation'):
            tensors.append(np.asarray(
                self.arrays[f'{role}_delta'][indices], dtype='float32') + baseline)
        if not all(np.isfinite(array).all() for array in tensors):
            raise RunFailure('non-finite predictor input')
        return [torch.from_numpy(array).to(device) for array in tensors]


def minibatches(count, batch_size, generator=None, training=False):
    indices = (generator.permutation(count) if generator is not None
               else np.arange(count))
    if training and count < 2:
        raise RunFailure('native PRnet BatchNorm needs at least two training cases')
    start = 0
    while start < count:
        end = min(start + batch_size, count)
        if training and count - end == 1:
            end = count
        yield indices[start:end]
        start = end


class JointPRnet(nn.Module):
    """The ADR 0063 seam: shared Morgan adapter plus the unmodified native PGM."""

    def __init__(self, n_genes, config):
        super().__init__()
        self.adapter = nn.Linear(2048, config['drug_dimension'])
        self.predictor = PGM(n_genes, config['comb_dimension'],
                             config['noise_dimension'],
                             config['hidden_layer_sizes'],
                             config['z_dimension'],
                             config['adaptor_layer_sizes'],
                             config['drug_dimension'], config['dropout'])
        # Upstream PRnetTrainer.weight_init overrides native BatchNorm defaults.
        for module in self.predictor.modules():
            if isinstance(module, nn.BatchNorm1d):
                nn.init.normal_(module.weight, 1., .02)
                if module.bias is not None:
                    nn.init.zeros_(module.bias)

    def forward(self, baseline, fingerprints, dose, noise):
        chemical = self.adapter(fingerprints) * torch.log10(dose[:, None] + 1.)
        return self.predictor(baseline, chemical, noise)


def device_for(config):
    if config['device'] == 'cpu':
        return torch.device('cpu')
    if config['device'] == 'cuda':
        if not torch.cuda.is_available():
            raise RunFailure('config requires CUDA but no device is available')
        return torch.device('cuda')
    return torch.device('cuda' if torch.cuda.is_available() else 'cpu')


def predict_deltas(model, release, config, device, seed):
    """Frozen test prediction: one batch of (indices, delta) at a time.

    Both the published prediction and the fresh-process replay consume this one
    generator, so a replay from ``best.pt`` is exact by construction.
    """
    n_genes = release.n_genes
    generator = torch.Generator(device=device).manual_seed(seed + 200000)
    with torch.no_grad():
        for indices in minibatches(len(release.cases['test']), config['batch_size']):
            baseline, fingerprint, dose = release.batch('test', indices, device)
            noise = torch.randn(len(indices), config['noise_dimension'],
                                generator=generator, device=device)
            mean = model(baseline, fingerprint, dose, noise)[:, :n_genes]
            yield indices, (mean - baseline).cpu().numpy()


def replay_prediction(request, config, directory):
    """Recompute the published prediction from ``best.pt`` in a fresh process."""
    directory = Path(directory)
    release = Release(request)
    torch.set_num_threads(int(config['threads']))
    device = device_for(config)
    config_sha = config_digest(config)
    if not (directory / 'best.pt').is_file():
        raise RunFailure(f'missing checkpoint for replay: {directory / "best.pt"}')
    checkpoint = torch.load(directory / 'best.pt', map_location=device,
                            weights_only=False)
    if checkpoint['identity'] != release.identity:
        raise RunFailure('stale selected checkpoint identity')
    if checkpoint['config_sha256'] != config_sha:
        raise RunFailure('stale selected checkpoint config')
    model = JointPRnet(release.n_genes, config).to(device)
    model.load_state_dict(checkpoint['model_state'])
    model.eval()
    partial = directory / 'replay.npy.partial'
    replay = np.lib.format.open_memmap(
        partial, mode='w+', dtype='float32',
        shape=(len(release.cases['test']), release.n_genes))
    for indices, delta in predict_deltas(model, release, config, device,
                                         release.identity['seed']):
        if not np.isfinite(delta).all():
            raise RunFailure('nonfinite replay prediction')
        replay[indices] = delta
    replay.flush()
    del replay
    partial.replace(directory / 'replay.npy')
    summary = {
        'replay_sha256': sha256(directory / 'replay.npy'),
        'shape': [len(release.cases['test']), release.n_genes],
        'best_epoch': checkpoint['epoch'],
        'config_sha256': config_sha,
        'device': str(device),
    }
    published = directory / 'prediction.npy'
    if published.is_file():
        summary['matches_published'] = bool(np.array_equal(
            np.load(published), np.load(directory / 'replay.npy')))
    return summary


def parameter_vector(module):
    return torch.cat([p.detach().cpu().reshape(-1) for p in module.parameters()])


def save_torch(path, value):
    temporary = path.with_suffix(path.suffix + '.partial')
    torch.save(value, temporary)
    temporary.replace(path)


def environment(device):
    return {
        'python': sys.version,
        'executable': sys.executable,
        'torch': torch.__version__,
        'numpy': np.__version__,
        'cuda': torch.version.cuda,
        'cuda_available': torch.cuda.is_available(),
        'device': str(device),
        'threads': torch.get_num_threads(),
        'numba_cache_dir': os.environ.get('NUMBA_CACHE_DIR'),
        'upstream_source': str(UPSTREAM),
        'upstream_prnet_module_sha256': sha256(UPSTREAM / 'models/PRnet.py'),
    }


def code_hashes():
    names = ('prnet_fit_predict.py', 'early_stopping.py', 'artifacts.py')
    return {name: sha256(SCRIPT_DIR / name) for name in names}


def fit_predict(request, config, release, directory, state):
    """Fit, select, reload, predict, and publish; raises on any declared failure."""
    seed = release.identity['seed']
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
        torch.cuda.reset_peak_memory_stats()
    device = device_for(config)
    torch.set_num_threads(int(config['threads']))
    n_genes = release.n_genes
    config_sha = config_digest(config)
    directory = Path(directory)

    state['environment'] = environment(device)
    state['device'] = str(device)
    state['started'] = time.monotonic()
    state['phase'] = 'training'
    write_json(directory / 'progress.json', {'phase': 'training', 'epoch': 0})

    model = JointPRnet(n_genes, config).to(device)
    initial_adapter = parameter_vector(model.adapter)
    initial_predictor = parameter_vector(model.predictor)
    optimizer = torch.optim.Adam(model.parameters(),
                                 lr=config['learning_rate'],
                                 weight_decay=config['weight_decay'])
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode='min', factor=config['scheduler_factor'],
        patience=config['scheduler_patience'],
        min_lr=config['minimum_learning_rate'])
    objective = nn.GaussianNLLLoss()
    stopper = EarlyStopper(
        max_epochs=config['max_epochs'],
        no_strict_decrease_epochs=config['no_strict_decrease_epochs'],
        above_running_min_epochs=config['above_running_min_epochs'])
    rng = np.random.default_rng(seed)
    history = []
    training_start = time.monotonic()
    best = float('inf')
    best_epoch = None
    stop_reason = 'max_epochs'

    def snapshot(epoch):
        return {
            'identity': release.identity,
            'config': config,
            'config_sha256': config_sha,
            'model_state': model.state_dict(),
            'epoch': epoch,
            'best_validation_loss': best,
            'training_seconds': time.monotonic() - training_start,
        }

    for epoch in range(1, config['max_epochs'] + 1):
        check_stop()
        model.train()
        total = 0.
        count = 0
        for indices in minibatches(len(release.cases['train']),
                                   config['batch_size'], rng, training=True):
            check_stop()
            if time.monotonic() - training_start >= config['max_training_seconds']:
                raise OverBudget(
                    f'training budget of {config["max_training_seconds"]} s '
                    f'exceeded in epoch {epoch}')
            baseline, fingerprint, dose, target = release.batch(
                'train', indices, device)
            optimizer.zero_grad(set_to_none=True)
            output = model(baseline, fingerprint, dose,
                           torch.randn(len(indices), config['noise_dimension'],
                                       device=device))
            mean, variance = output.chunk(2, dim=1)
            loss = objective(mean, target, F.softplus(variance))
            if not torch.isfinite(loss):
                raise RunFailure(
                    f'nonfinite author training objective in epoch {epoch}')
            loss.backward()
            optimizer.step()
            total += float(loss.detach()) * len(indices)
            count += len(indices)
            write_json(directory / 'progress.json',
                       {'phase': 'training', 'epoch': epoch})

        check_stop()
        model.eval()
        validation_total = 0.
        validation_count = 0
        valid_generator = torch.Generator(device=device).manual_seed(seed + 100000)
        with torch.no_grad():
            for indices in minibatches(len(release.cases['validation']),
                                       config['batch_size']):
                check_stop()
                if time.monotonic() - training_start >= config['max_training_seconds']:
                    raise OverBudget(
                        f'training budget of {config["max_training_seconds"]} s '
                        f'exceeded during validation of epoch {epoch}')
                baseline, fingerprint, dose, target = release.batch(
                    'validation', indices, device)
                noise = torch.randn(len(indices), config['noise_dimension'],
                                    generator=valid_generator, device=device)
                mean, variance = model(baseline, fingerprint, dose, noise).chunk(2, dim=1)
                validation_total += float(
                    objective(mean, target, F.softplus(variance))) * len(indices)
                validation_count += len(indices)
        validation = validation_total / validation_count
        if not np.isfinite(validation):
            raise RunFailure(f'nonfinite author validation objective in epoch {epoch}')
        decision = stopper.update(epoch, validation)
        improved = decision['best_epoch'] == epoch
        scheduler.step(validation)
        history.append({
            'epoch': epoch,
            'training_loss': total / count,
            'validation_loss': validation,
            'learning_rate': optimizer.param_groups[0]['lr'],
            'improved': improved,
        })
        if improved:
            best = validation
            best_epoch = epoch
            save_torch(directory / 'best.pt', snapshot(epoch))
        write_json(directory / 'progress.json', {
            'phase': 'training', 'epoch': epoch, 'best_epoch': best_epoch,
            'best_validation_loss': best})
        if decision['stop']:
            stop_reason = decision['reason']
            break
    training_seconds = time.monotonic() - training_start
    if best_epoch is None or not (directory / 'best.pt').is_file():
        raise RunFailure('fit produced no validation checkpoint')
    stopping = stopper.record()
    state['early_stopping'] = stopping
    if stop_reason not in STOP_REASONS:
        raise RunFailure(f'unexpected early-stopping reason: {stop_reason}')
    if not stopping['dominance_observed']:
        raise RunFailure(
            'the recorded early-stopping dominance accounting is inverted')
    # Real selection accounting: the saved checkpoint must be an argmin of the
    # observed validation-loss history, not merely "the last improved epoch".
    minimum_loss = min(entry['validation_loss'] for entry in history)
    minimizers = [entry['epoch'] for entry in history
                  if entry['validation_loss'] == minimum_loss]
    if best_epoch != min(minimizers) or best != minimum_loss:
        raise RunFailure(
            f'selected checkpoint epoch {best_epoch} loss {best!r} is not the '
            f'lowest-validation-loss epoch {min(minimizers)} loss {minimum_loss!r}')

    # Reload the selected checkpoint from disk; the published prediction must come
    # from this reloaded artifact, not from in-memory training state.
    inference_start = time.monotonic()
    state['phase'] = 'reload'
    write_json(directory / 'progress.json', {
        'phase': 'reload', 'training_seconds': training_seconds,
        'best_epoch': best_epoch, 'stop_reason': stop_reason})
    checkpoint = torch.load(directory / 'best.pt', map_location=device,
                            weights_only=False)
    if checkpoint['identity'] != release.identity:
        raise RunFailure('stale selected checkpoint identity')
    if checkpoint['config_sha256'] != config_sha:
        raise RunFailure('stale selected checkpoint config')
    if checkpoint['epoch'] != best_epoch:
        raise RunFailure(
            f'saved checkpoint epoch {checkpoint["epoch"]} is not the selected '
            f'epoch {best_epoch}')
    model.load_state_dict(checkpoint['model_state'])
    model.eval()
    probe_count = min(2, len(release.cases['test']))
    probe = release.batch('test', np.arange(probe_count), device)
    probe_noise = torch.zeros(probe_count, config['noise_dimension'], device=device)
    with torch.no_grad():
        original = model(*probe, probe_noise).detach()
        changed = model(torch.zeros_like(probe[0]), probe[1], probe[2], probe_noise)
        baseline_sensitivity = float(torch.linalg.vector_norm(
            original[:, :n_genes] - changed[:, :n_genes]))
        adapter_change = float(torch.linalg.vector_norm(
            parameter_vector(model.adapter) - initial_adapter))
        predictor_change = float(torch.linalg.vector_norm(
            parameter_vector(model.predictor) - initial_predictor))
        replayed_model = JointPRnet(n_genes, config).to(device)
        replayed_model.load_state_dict(torch.load(
            directory / 'best.pt', map_location=device,
            weights_only=False)['model_state'])
        replayed_model.eval()
        replayed = replayed_model(*probe, probe_noise)
        reload_error = float(torch.max(torch.abs(original - replayed)))
    if adapter_change <= 0 or predictor_change <= 0:
        raise RunFailure('joint adapter/predictor training gate failed')
    if not baseline_sensitivity > 0:
        raise RunFailure('baseline-expression sensitivity gate failed')
    if reload_error != 0.:
        raise RunFailure(f'checkpoint reload is not exact: {reload_error}')
    behavior = {
        'adapter_parameter_change_l2': adapter_change,
        'predictor_parameter_change_l2': predictor_change,
        'baseline_sensitivity_l2': baseline_sensitivity,
        'reload_max_abs_error': reload_error,
    }
    trainable_parameters = sum(p.numel() for p in model.parameters())
    predictor_parameters = sum(p.numel() for p in model.predictor.parameters())
    del replayed_model, probe, probe_noise, original, changed

    # ``model`` now holds the reloaded disk state; every published value below
    # is produced by it, never by pre-reload training state.
    state['phase'] = 'inference'
    write_json(directory / 'progress.json', {
        'phase': 'inference', 'training_seconds': training_seconds,
        'inference_seconds': time.monotonic() - inference_start})
    partial = directory / 'prediction.npy.partial'
    prediction = np.lib.format.open_memmap(
        partial, mode='w+', dtype='float32',
        shape=(len(release.cases['test']), n_genes))
    for indices, delta in predict_deltas(model, release, config, device, seed):
        check_stop()
        if time.monotonic() - inference_start >= config['inference_max_seconds']:
            raise OverBudget(
                f'inference budget of {config["inference_max_seconds"]} s exceeded')
        if not np.isfinite(delta).all():
            raise RunFailure('nonfinite full-axis prediction')
        prediction[indices] = delta
    prediction.flush()
    del prediction
    partial.replace(directory / 'prediction.npy')
    inference_seconds = time.monotonic() - inference_start
    case_ids = [row['case_id'] for row in release.cases['test']]
    write_json(directory / 'case_ids.json', case_ids)

    peak_rss = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    result = {
        'identity': release.identity,
        'author_objective': 'GaussianNLLLoss',
        'n_genes': n_genes,
        'n_cases': len(case_ids),
        'target_kind': 'delta',
        'selection': SELECTION,
        'output': OUTPUT,
        'best_epoch': best_epoch,
        'best_validation_loss': best,
        'history': history,
        'stop_reason': stop_reason,
        'early_stopping': stopping,
        'training_seconds': training_seconds,
        'inference_seconds': inference_seconds,
        'trainable_parameters': trainable_parameters,
        'adapter_parameters': 2048 * 1024 + 1024,
        'predictor_parameters': predictor_parameters,
        'behavior': behavior,
        'config': config,
        'config_sha256': config_sha,
        'device': str(device),
        'environment': state['environment'],
        'peak_rss_kbytes': peak_rss,
        'peak_gpu_bytes': int(torch.cuda.max_memory_allocated())
        if torch.cuda.is_available() else 0,
        'blobs': {'best.pt': sha256(directory / 'best.pt')},
    }
    write_json(directory / 'checkpoint.json', result)
    identity = {
        'identity': release.identity,
        'n_genes': n_genes,
        'n_cases': len(case_ids),
        'release_manifest_sha256': sha256(release.manifest_path),
        'release_files_sha256': {
            name: sha256(path) for name, path in release.files.items()},
        'config_sha256': sha256(directory / 'config.json'),
        'config_recipe_sha256': config_sha,
        'code_sha256': code_hashes(),
        'checkpoint_sha256': sha256(directory / 'best.pt'),
        'prediction_sha256': sha256(directory / 'prediction.npy'),
        'case_ids_sha256': sha256(directory / 'case_ids.json'),
        'environment': state['environment'],
        'wall_seconds': time.monotonic() - state['started'],
        'peak_rss_kbytes': peak_rss,
    }
    write_json(directory / 'identity.json', identity)
    state['phase'] = 'complete'
    write_json(directory / 'progress.json', {
        'phase': 'complete', 'training_seconds': training_seconds,
        'inference_seconds': inference_seconds})
    return result


def write_terminal(directory, status, phase, state, error=None):
    directory = Path(directory)
    if status != 'complete':
        for leftover in ('prediction.npy', 'prediction.npy.partial'):
            path = directory / leftover
            if path.exists():
                path.unlink()
    record = {
        'status': status,
        'phase': phase,
        'error': error,
        'prediction': None if status != 'complete' else 'prediction.npy',
        'wall_seconds': time.monotonic() - state.get('started', time.monotonic()),
        'peak_rss_kbytes': int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
        'finished_utc': datetime.now(timezone.utc).isoformat(),
        'early_stopping': state.get('early_stopping'),
    }
    write_json(directory / 'terminal.json', record)
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--request', required=True)
    parser.add_argument('--config',
                        default=str(SCRIPT_DIR.parent / 'models/prnet/config.json'))
    parser.add_argument(
        '--replay', action='store_true',
        help='recompute the published prediction from best.pt in this fresh '
             'process (writes replay.npy; never rewrites run artifacts)')
    args = parser.parse_args(argv)
    signal.signal(signal.SIGTERM, _request_stop)
    signal.signal(signal.SIGINT, _request_stop)
    request = read_json(args.request)
    directory = Path(request['output_directory'])
    if args.replay:
        try:
            summary = replay_prediction(request, load_config(args.config), directory)
            print(json.dumps(summary, sort_keys=True), flush=True)
            return 0
        except Exception as error:  # noqa: BLE001 - replay is a read-only check
            traceback.print_exc()
            print(f'PRnet replay failed: {error}', file=sys.stderr)
            return 1
    directory.mkdir(parents=True, exist_ok=True)
    for name in ('prediction.npy', 'prediction.npy.partial', 'terminal.json',
                 'checkpoint.json', 'identity.json', 'case_ids.json', 'best.pt'):
        stale = directory / name
        if stale.exists():
            stale.unlink()
    state = {'started': time.monotonic(), 'phase': 'validate'}
    try:
        config = load_config(args.config)
        write_json(directory / 'config.json', config)
        release = Release(request)
        result = fit_predict(request, config, release, directory, state)
        write_terminal(directory, 'complete', 'complete', state)
        print(json.dumps({
            'best_epoch': result['best_epoch'],
            'best_validation_loss': result['best_validation_loss'],
            'stop_reason': result['stop_reason'],
            'training_seconds': result['training_seconds'],
            'inference_seconds': result['inference_seconds'],
            'behavior': result['behavior'],
        }), flush=True)
        return 0
    except OverBudget as error:
        write_terminal(directory, 'over_budget', state.get('phase'), state,
                       str(error))
        print(f'PRnet run over budget: {error}', file=sys.stderr)
    except Interrupted as error:
        write_terminal(directory, 'interrupted', state.get('phase'), state,
                       str(error))
        print(f'PRnet run interrupted: {error}', file=sys.stderr)
    except Exception as error:  # noqa: BLE001 - every failure gets a terminal record
        traceback.print_exc()
        write_terminal(directory, 'failed', state.get('phase'), state,
                       f'{type(error).__name__}: {error}')
        print(f'PRnet run failed: {error}', file=sys.stderr)
    return 1


if __name__ == '__main__':
    sys.exit(main())
