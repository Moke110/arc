"""Shared strict stopping and isolated fixed evaluation randomness."""
import math
from contextlib import contextmanager


class EarlyStopping:
    def __init__(self, patience=50, max_epochs=1000):
        if patience < 1 or not 1 <= max_epochs <= 1000:
            raise ValueError('invalid stopping budget')
        self.patience, self.max_epochs = patience, max_epochs
        self.best, self.best_epoch, self.stale = math.inf, None, 0
        self.stopped = False

    def observe(self, value, epoch):
        if self.stopped or not math.isfinite(value):
            raise ValueError('stopped or nonfinite validation')
        improved = value < self.best
        if improved:
            self.best, self.best_epoch, self.stale = value, epoch, 0
        else:
            self.stale += 1
        self.stopped = self.stale >= self.patience or epoch >= self.max_epochs
        return improved


@contextmanager
def fixed_rng(seed, device='cpu'):
    import torch
    devices = [] if str(device) == 'cpu' else [torch.device(device).index or 0]
    with torch.random.fork_rng(devices=devices):
        torch.manual_seed(seed)
        if devices:
            torch.cuda.manual_seed_all(seed)
        yield
