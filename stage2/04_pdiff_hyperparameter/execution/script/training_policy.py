"""The pinned PerturbDiff monitor and optimizer-step stopping policy."""
from dataclasses import dataclass
import math

MONITOR = 'validation_validation_loss_epoch'


def native_batch_mean(loss, weights):
    return (loss * weights).mean()


def native_epoch_mean(batch_scalars, stack):
    return float(stack(batch_scalars).mean())


def verify_best_checkpoint(run, subject):
    history = run['history']
    steps = [event['step'] for event in history]
    if (not steps or any(type(step) is not int or step < 1 for step in steps)
            or steps != sorted(set(steps))
            or not all(math.isfinite(event['validation_loss']) for event in history)):
        raise ValueError(f'invalid {subject} validation history')
    best = min(history, key=lambda event: event['validation_loss'])
    if (type(run['best_step']) is not int
            or (run['best_step'], run['best_validation_loss']) != (best['step'], best['validation_loss'])):
        raise ValueError(f'{subject} checkpoint differs from the strict validation minimum')
    return best


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
