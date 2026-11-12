"""Epoch-based early stopping for the PRnet entrant (ADR 0063).

Two independent stop conditions are implemented and recorded:

* ``no_strict_decrease`` - stop after 5 consecutive epochs in which the
  validation loss does not strictly decrease;
* ``above_running_min``   - stop after 20 consecutive epochs in which the
  validation loss is strictly above the running minimum.

The run never exceeds 100 epochs.  The lowest-validation-loss epoch is the
selected checkpoint, which the fit loop owns.

The 5-epoch rule dominates the 20-epoch rule: every epoch counted by the
20-epoch counter is also an epoch without a strict decrease, and the counters
are consecutive, so a 20-epoch run above the minimum always contains five
consecutive non-decreasing epochs.  The 5-epoch condition therefore fires
first, and ``record()`` writes that dominance instead of presenting the two
rules as equally likely to stop a run.  The 20-epoch rule is still implemented
independently (its counter is checked on its own), because the ADR requires
both conditions.
"""

RULES = {
    'no_strict_decrease_epochs': 5,
    'above_running_min_epochs': 20,
    'max_epochs': 100,
}

DOMINANCE = (
    'The 5-epoch rule dominates the 20-epoch rule: any 20 consecutive epochs '
    'above the running minimum are also 20 consecutive epochs without a strict '
    'decrease, so the 5-epoch condition fires first and the 20-epoch condition '
    'cannot be the first to stop a default-configuration run.'
)

STOP_REASONS = ('no_strict_decrease', 'above_running_min', 'max_epochs')


class EarlyStopper:
    """Stateful counter pair over one run's validation-loss trajectory."""

    def __init__(self, max_epochs=None, no_strict_decrease_epochs=None,
                 above_running_min_epochs=None):
        self.max_epochs = _positive('max_epochs', max_epochs, RULES['max_epochs'])
        self.no_strict_decrease_epochs = _positive(
            'no_strict_decrease_epochs', no_strict_decrease_epochs,
            RULES['no_strict_decrease_epochs'])
        self.above_running_min_epochs = _positive(
            'above_running_min_epochs', above_running_min_epochs,
            RULES['above_running_min_epochs'])
        self.best = float('inf')
        self.best_epoch = None
        self.no_strict_decrease = 0
        self.above_running_min = 0
        self.stop_reason = None

    def update(self, epoch, validation_loss):
        """Consume one epoch; return the stop decision and the counters seen."""
        if isinstance(epoch, bool) or not isinstance(epoch, int) or epoch < 1:
            raise ValueError(f'invalid epoch: {epoch!r}')
        value = float(validation_loss)
        if value != value or value in (float('inf'), float('-inf')):
            raise ValueError(f'nonfinite validation loss: {validation_loss!r}')
        previous = self.best
        if value < previous:
            self.best = value
            self.best_epoch = epoch
            self.no_strict_decrease = 0
            self.above_running_min = 0
        else:
            self.no_strict_decrease += 1
            if value > previous:
                self.above_running_min += 1
            else:  # a tie with the minimum is not "above" it
                self.above_running_min = 0
        reason = None
        if self.no_strict_decrease >= self.no_strict_decrease_epochs:
            reason = 'no_strict_decrease'
        elif self.above_running_min >= self.above_running_min_epochs:
            reason = 'above_running_min'
        elif epoch >= self.max_epochs:
            reason = 'max_epochs'
        if reason is not None:
            self.stop_reason = reason
        return {
            'epoch': epoch,
            'stop': reason is not None,
            'reason': reason,
            'best_validation_loss': self.best,
            'best_epoch': self.best_epoch,
            'no_strict_decrease': self.no_strict_decrease,
            'above_running_min': self.above_running_min,
        }

    def record(self):
        """Machine-readable rule declaration plus the observed counters."""
        return {
            'rules': dict(self.rules()),
            'dominance': DOMINANCE,
            'dominance_observed': self.above_running_min <= self.no_strict_decrease,
            'stop_reason': self.stop_reason,
            'observed': {
                'best_validation_loss': self.best,
                'best_epoch': self.best_epoch,
                'no_strict_decrease': self.no_strict_decrease,
                'above_running_min': self.above_running_min,
            },
        }

    def rules(self):
        return {
            'no_strict_decrease_epochs': self.no_strict_decrease_epochs,
            'above_running_min_epochs': self.above_running_min_epochs,
            'max_epochs': self.max_epochs,
        }


def _positive(name, value, default):
    if value is None:
        return default
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f'{name} must be a positive integer: {value!r}')
    return value
