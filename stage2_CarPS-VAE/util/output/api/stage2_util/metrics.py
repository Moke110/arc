"""Generic expression-vector metrics shared by Stage 2 studies."""

import math

import numpy as np
from scipy.stats import rankdata


def vector_metrics(truth, prediction):
    """Score paired vectors with the accepted constant-truth mask."""
    truth = np.asarray(truth, dtype=np.float64)
    prediction = np.asarray(prediction, dtype=np.float64)
    if (truth.ndim != 1 or truth.size == 0 or truth.shape != prediction.shape
            or not np.isfinite(truth).all() or not np.isfinite(prediction).all()):
        raise ValueError("expected paired nonempty, finite vectors")

    constant_truth = bool(np.ptp(truth) == 0)
    constant_prediction = bool(np.ptp(prediction) == 0)
    squared_error = (truth - prediction) ** 2
    denominator = float(np.sum((truth - truth.mean()) ** 2))
    if constant_truth:
        pcc = spearman = float("nan")
    elif constant_prediction:
        pcc = spearman = 0.0
    else:
        pcc = float(np.corrcoef(truth, prediction)[0, 1])
        spearman = float(np.corrcoef(rankdata(truth), rankdata(prediction))[0, 1])
    return {
        "pcc": pcc,
        "spearman": spearman,
        "rmse": math.sqrt(float(squared_error.mean())),
        "r2": 1 - float(squared_error.sum()) / denominator if denominator else float("nan"),
        "constant_truth": constant_truth,
        "constant_prediction": constant_prediction,
    }
