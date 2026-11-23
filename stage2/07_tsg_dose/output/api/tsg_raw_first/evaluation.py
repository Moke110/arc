"""Compute the dose study's sample-equal Delta metrics on verified profiles."""
import numpy as np
import pandas as pd

from stage2_util import metrics
from tsg_round1.protocol import METRICS, summary

MAX_ROWS = 128


def score_prediction(obs, truth, prediction, control, prediction_ids=None):
    """Score each row against the verified raw-first matched-control expression."""
    if prediction_ids is not None and list(map(str, prediction_ids)) != list(map(str, obs.index)):
        raise ValueError("prediction identity order mismatch")
    truth, prediction, control = map(np.asarray, (truth, prediction, control))
    if truth.ndim != 2 or truth.shape != prediction.shape or truth.shape != control.shape or len(obs) != len(truth):
        raise ValueError("prediction identity or shape mismatch")
    values = []
    for measured, predicted, matched in zip(truth, prediction, control, strict=True):
        delta = measured.astype(np.float64) - matched
        predicted_delta = predicted.astype(np.float64) - matched
        if np.isfinite(delta).all() and np.isfinite(predicted_delta).all():
            row = metrics.vector_metrics(delta, predicted_delta)
            row["nonfinite"] = False
        else:
            row = dict.fromkeys(METRICS, float("nan"))
            row.update(constant_truth=False, constant_prediction=False, nonfinite=True)
        values.append(row)
    identities = pd.DataFrame({"profile_id": obs.index.astype(str),
                               "drug_parent": obs.drug_parent_inchikey.to_numpy(),
                               "dose_um": obs.dose_um.to_numpy()})
    if "cell_line" in obs:
        identities["cell_line"] = obs.cell_line.to_numpy()
    return identities.join(pd.DataFrame(values))


def score(data, scenario, fold, partition, prediction, prediction_ids):
    """Validate partition and prediction order before reading bounded test windows."""
    rows = data.rows(scenario, fold, partition)
    if not len(rows) or prediction.shape != (len(rows), len(data.genes)):
        raise ValueError("prediction shapes or partition coverage differ")
    if list(map(str, prediction_ids)) != list(map(str, data.obs.index[rows])):
        raise ValueError("prediction identity order mismatch")
    outputs = []
    for start in range(0, len(rows), MAX_ROWS):
        stop = start + MAX_ROWS
        indices = rows[start:stop]
        outputs.append(score_prediction(data.obs.iloc[indices], data.expression[indices],
                                        prediction[start:stop], data.control[indices]))
    samples = pd.concat(outputs, ignore_index=True)
    return samples, summary(samples)
