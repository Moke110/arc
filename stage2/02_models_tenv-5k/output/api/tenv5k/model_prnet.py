"""Native Sci-Plex PRnet: random weights, GaussianNLL and Normal draws.

Only chemical interface changes: shared binary Morgan2048 replaces upstream
feature-Morgan1024. The upstream constant log10(1+dose_um) multiplier remains.
Source files are imported unchanged and their digests are checkpoint identity.
"""
import importlib.util
import json
import math
from pathlib import Path
import platform
import time
import resource
import numpy as np
import pandas as pd
import torch
from torch.nn import functional as F
from .inputs import Inputs, sha256
from .learning import EarlyStopping, fixed_rng
from .protocol import score_prediction

UPSTREAM = 'f19174bde3ed2633f54c7831799cc38c4ffc7a0d'
SOURCE_DIGEST = '024f39a07e2887fea45a42f60b07b3ba17dc7501f134285ff72798eacb8bea3e'
CONFIG = dict(method='PRnet', x_dimension=5000, hidden_layer_sizes=[128],
              z_dimension=64, adaptor_layer_sizes=[128], comb_dimension=64,
              comb_num=1, drug_dimension=2048, dr_rate=0.05, batch_size=512,
              lr=1e-3, weight_decay=1e-8, scheduler_factor=0.5,
              scheduler_patience=5, scheduler_min_lr=1e-8,
              loss='GaussianNLLLoss(mean,target,softplus(raw_variance))',
              inference='native Normal(clamp(mean),clamp(sqrt(softplus(var)))).sample',
              chemical_interface='shared Morgan2048 * log10(1+1uM)',
              initialization='native Linear N(0,.02), BatchNorm N(1,.02), zero bias',
              pretrained=False, upstream_commit=UPSTREAM,
              adaptations=['5000 released gene axis', '2048 shared Morgan replaces feature Morgan1024',
                           'shared readers and splits', 'fixed-validation treated MSE stopping'],
              upstream_files={'models/PRnet.py': SOURCE_DIGEST,
                 'train_sciplex.py': 'e1a99dd5d9c344e8dc6ee9273faf7e997d8eadf8347b790c8bc46add2f5366bd',
                 'trainer/PRnetTrainer.py': '0d0eb2d8e19ce284d56313667fc2ffb532050ede8b93b4a40bc5ffa6329e8038'})


def build_model(source, seed, device='cpu'):
    source = Path(source)
    path = source / 'models/PRnet.py'
    for name, expected in CONFIG['upstream_files'].items():
        if sha256(source / name) != expected:
            raise ValueError('PRnet native source digest mismatch: ' + name)
    torch.manual_seed(seed)
    spec = importlib.util.spec_from_file_location('tenv5k_native_prnet', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    model = module.PGM(5000, 64, 10, [128], 64, [128], 2048, .05)
    def initialize(layer):
        if isinstance(layer, torch.nn.Linear):
            torch.nn.init.normal_(layer.weight, 0, .02)
            if layer.bias is not None: torch.nn.init.zeros_(layer.bias)
        elif isinstance(layer, torch.nn.BatchNorm1d):
            torch.nn.init.normal_(layer.weight, 1, .02)
            if layer.bias is not None: torch.nn.init.zeros_(layer.bias)
    # Native trainer resets its run seed immediately before weight_init.
    torch.manual_seed(seed)
    model.apply(initialize)
    return model.to(device)


def batch(data, rows, device):
    control = torch.as_tensor(np.array(data.control[rows]), device=device)
    chemical = torch.as_tensor(data.fingerprints(rows), device=device) * math.log10(2)
    return control, chemical


def sampled_prediction(model, data, rows, *, seed, device='cpu', batch_size=512):
    """Predict treated expression with no access to any treated labels."""
    model.eval()
    values = []
    with fixed_rng(seed, device), torch.no_grad():
        for start in range(0, len(rows), batch_size):
            selected = rows[start:start + batch_size]
            control, chemical = batch(data, selected, device)
            noise = torch.randn(len(selected), 10).to(device)
            output = model(control, chemical, noise)
            mean, raw_var = output.chunk(2, dim=1)
            distribution = torch.distributions.Normal(mean.clamp(1e-3, 1e3), F.softplus(raw_var).sqrt().clamp(1e-3, 1e3))
            values.append(distribution.sample().cpu().numpy())
    prediction = np.concatenate(values).astype(np.float32)
    if not np.isfinite(prediction).all():
        raise ValueError('nonfinite native prediction')
    return prediction


def make_optimizer(model):
    return torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-8)


def make_scheduler(optimizer):
    return torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, factor=.5, patience=5, min_lr=1e-8)


def train_epoch(model, data, rows, optimizer, shuffle, device):
    model.train()
    order = shuffle.permutation(rows)
    batches = [order[start:start + 512] for start in range(0, len(order), 512)]
    if len(batches[-1]) == 1:
        if len(batches) == 1: raise ValueError('singleton train split')
        batches[-2] = np.concatenate([batches[-2], batches[-1]])
        batches.pop()
    loss_sum = 0.
    for selected in batches:
        control, chemical = batch(data, selected, device)
        target = torch.as_tensor(np.array(data.expression[selected]), device=device)
        optimizer.zero_grad()
        mean, raw_var = model(control, chemical, torch.randn(len(selected), 10).to(device)).chunk(2, dim=1)
        loss = F.gaussian_nll_loss(mean, target, F.softplus(raw_var))
        if not torch.isfinite(loss): raise ValueError('nonfinite training loss')
        loss.backward()
        optimizer.step()
        loss_sum += loss.item() * len(selected)
    return loss_sum / len(rows)
