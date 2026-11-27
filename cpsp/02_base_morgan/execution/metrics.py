"""Fixed three-direction metrics, shared by training and saved-prediction evaluation."""
import numpy as np

CLASS_ORDER = (-1, 0, 1)
METRIC_CONVENTION = {
    'class_order': list(CLASS_ORDER),
    'zero_division': 0,
    'probability_sum_atol': 1e-6,
    'log_loss_dtype': 'float64',
    'log_loss_epsilon': float(np.finfo(np.float64).eps),
    'log_loss_clip': '[epsilon, 1-epsilon], then row-normalize',
    'argmax_ties': 'first in class_order',
}


def validate_probabilities(probabilities):
    """Return float64 N×3 probabilities; reject invalid or empty inputs."""
    values = np.asarray(probabilities, dtype=np.float64)
    if values.ndim != 2 or values.shape[1] != 3 or not len(values):
        raise ValueError('probabilities must have nonempty shape (N, 3)')
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('probabilities must be finite and nonnegative')
    if not np.allclose(values.sum(axis=1), 1, atol=1e-6, rtol=0):
        raise ValueError('probabilities must sum to one within atol=1e-6')
    return values


def hard_predictions(probabilities):
    """Map argmax columns to directions, resolving ties in fixed class order."""
    return np.asarray(CLASS_ORDER)[validate_probabilities(probabilities).argmax(axis=1)]


def evaluate_predictions(y_true, probabilities, y_pred=None):
    """Return JSON-safe metrics, checking any supplied hard predictions against argmax."""
    values = validate_probabilities(probabilities)
    truth = np.asarray(y_true)
    if truth.shape != (len(values),) or not np.isin(truth, CLASS_ORDER).all():
        raise ValueError('truth must contain one direction in [-1, 0, 1] per row')
    predicted = np.asarray(CLASS_ORDER)[values.argmax(axis=1)]
    if y_pred is not None:
        supplied = np.asarray(y_pred)
        if supplied.shape != predicted.shape or not np.array_equal(supplied, predicted):
            raise ValueError('hard predictions must match probability argmax in class order')
    confusion = np.zeros((3, 3), dtype=np.int64)
    np.add.at(confusion, (truth.astype(int) + 1, predicted + 1), 1)
    correct = confusion.diagonal()
    precision = np.divide(correct, confusion.sum(axis=0), out=np.zeros(3), where=confusion.sum(axis=0) != 0)
    recall = np.divide(correct, confusion.sum(axis=1), out=np.zeros(3), where=confusion.sum(axis=1) != 0)
    f1 = np.divide(2 * precision * recall, precision + recall, out=np.zeros(3), where=precision + recall != 0)
    epsilon = METRIC_CONVENTION['log_loss_epsilon']
    safe = np.clip(values, epsilon, 1 - epsilon)
    safe /= safe.sum(axis=1, keepdims=True)
    loss = -np.log(safe[np.arange(len(truth)), truth.astype(int) + 1]).mean()
    return {'class_order': list(CLASS_ORDER), 'n': len(truth), 'macro_f1': float(f1.mean()),
            'precision': precision.tolist(), 'recall': recall.tolist(), 'f1': f1.tolist(),
            'confusion_matrix': confusion.tolist(), 'log_loss': float(loss)}
