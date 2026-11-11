"""Stage 2 reporting metrics; ranking labels remain only for frozen input replay."""
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0, str(next(parent for parent in Path(__file__).resolve().parents if (parent / '.git').exists())))
from stage2.evaluation.metrics import METRICS as RETRIEVAL_METRICS, retrieval

METRICS = (*RETRIEVAL_METRICS, 'PCC', 'Spearman', 'RMSE', 'R2')
RECALL_GRID = np.linspace(0, 1, 101)
CURVE_CONVENTION = 'Full-CM direction-matched significance, fixed-group pooled thresholds, attainable endpoints'


def directional_curve(*args, **kwargs):
    raise ValueError('Rank-defined directional AP is retired. Use stage2/evaluation/evaluate.py with full-CM significance labels.')


def truth_labels(truth):
    truth = np.asarray(truth)
    labels = []
    for sign in (1, -1):
        axis = sign * truth
        eligible = np.flatnonzero(axis > 0)
        selected = eligible[np.argsort(-axis[eligible], kind='stable')[:100]]
        label = np.zeros(len(truth), dtype=bool)
        label[selected] = True
        labels.append(label)
    return labels



def score_group(truth, prediction, labels=None):
    """Require signed DESeq2 truth; old ranking-label tuples are rejected."""
    if labels is None or isinstance(labels, (tuple, list)):
        raise ValueError('Scoring requires signed DESeq2 labels; use stage2/evaluation/evaluate.py')
    observed, prediction = np.asarray(truth), np.asarray(prediction)
    curves = retrieval(prediction, labels, np.arange(len(prediction)))
    values = {name: curve.ap if curve is not None else np.nan for name, curve in curves.items()}
    reasons = {name: 'zero_significant' for name, curve in curves.items() if curve is None}
    from scipy.stats import rankdata
    a, b = observed - observed.mean(), prediction - prediction.mean()
    variance = float(a @ a)
    if variance and b @ b:
        values['PCC'] = float((a @ b) / np.sqrt(variance * (b @ b)))
        values['Spearman'] = float(np.corrcoef(rankdata(observed), rankdata(prediction))[0, 1])
    else:
        values['PCC'] = values['Spearman'] = np.nan
        reasons['PCC'] = reasons['Spearman'] = 'constant_truth_or_prediction'
    error = float((prediction - observed) @ (prediction - observed))
    values['RMSE'] = float(np.sqrt(error / len(observed)))
    values['R2'] = 1 - error / variance if variance else np.nan
    return values, reasons, curves
