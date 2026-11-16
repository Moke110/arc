"""XPert-defined chemical baseline with the project's preregistered configuration.

Official GitHub/Figshare did not supply complete baseline training code. This is
not a claim of reproducing an author implementation. Only Morgan2048 is input.
"""
import hashlib
from pathlib import Path
import numpy as np
import torch
from .learning import fixed_rng

CONFIG = dict(method='MLP_Morgan', batch_size=64, architecture=[2048, 128, 5000],
              hidden_activation='ReLU', output_activation='linear', loss='MSE',
              learning_rate=1e-3, weight_decay=1e-5, scheduler=None,
              baseline_input=False, dose_input=False, dropout=0,
              initialization_model='random', pretraining=False,
              singleton_train_batch='retain singleton (no BatchNorm)',
              configuration_origin='XPert-defined baseline; project preregistered configuration, not author-code reproduction',
              adaptations=['shared Morgan2048 and splits',
                           'shared fixed validation treated-expression MSE and best restore'])
SOURCE_DIGEST = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def build_model(source, seed, device):
    torch.manual_seed(seed)
    return torch.nn.Sequential(torch.nn.Linear(2048, 128), torch.nn.ReLU(),
                               torch.nn.Linear(128, 5000)).to(device)


def make_optimizer(model):
    return torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-5)


def make_scheduler(optimizer):
    return None


def train_epoch(model, data, rows, optimizer, shuffle, device):
    model.train()
    order = shuffle.permutation(rows)
    total = 0.0
    for start in range(0, len(order), 64):
        batch = order[start:start + 64]
        features = torch.as_tensor(data.fingerprints(batch), dtype=torch.float32, device=device)
        target = torch.tensor(np.asarray(data.expression[batch]), dtype=torch.float32, device=device)
        optimizer.zero_grad()
        loss = torch.nn.functional.mse_loss(model(features), target)
        if not torch.isfinite(loss):
            raise ValueError('non-finite MLP_Morgan loss')
        loss.backward()
        optimizer.step()
        total += loss.item() * len(batch)
    return total / len(order)


def sampled_prediction(model, data, rows, seed=42, device='cpu', batch_size=64):
    model.eval()
    result = []
    with fixed_rng(seed, device), torch.no_grad():
        for start in range(0, len(rows), batch_size):
            batch = rows[start:start + batch_size]
            features = torch.as_tensor(data.fingerprints(batch), dtype=torch.float32, device=device)
            result.append(model(features).cpu().numpy())
    return np.concatenate(result) if result else np.empty((0, 5000), np.float32)
