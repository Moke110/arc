"""The pinned PerturbDiff monitor and optimizer-step stopping policy."""
from dataclasses import dataclass
import math

MONITOR = 'validation_validation_loss_epoch'


def native_batch_mean(loss, weights):
    return (loss * weights).mean()


def native_epoch_mean(batch_scalars, stack):
    return float(stack(batch_scalars).mean())


@dataclass
class Selection:
    minimum_steps: int
    patience_steps: int
    min_delta: float
    best_loss: float = math.inf
    best_step: int | None = None
    reference_loss: float = math.inf
    reference_step: int = 0
    last_step: int = 0

    def observe(self, step, loss):
        if step <= self.last_step or not math.isfinite(loss):
            raise ValueError('validation requires increasing optimizer steps and finite loss')
        self.last_step = step
        strict_best = loss < self.best_loss
        if strict_best:
            self.best_loss, self.best_step = loss, step
        if self.reference_loss - loss > self.min_delta:
            self.reference_loss, self.reference_step = loss, step
        stop = step >= self.minimum_steps and step - self.reference_step >= self.patience_steps
        return {'step': step, 'validation_loss': loss, 'checkpoint_refreshed': strict_best,
                'patience_reference_step': self.reference_step,
                'patience_reference_loss': self.reference_loss, 'early_stop': stop}
