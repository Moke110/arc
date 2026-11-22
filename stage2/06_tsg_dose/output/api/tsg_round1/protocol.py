"""Source-specific parent partitions and sample-equal Delta metrics."""
import math

import numpy as np
import pandas as pd

from tsg_dose import shared  # Registers the accepted tenv5k scientific API.
from tenv5k.protocol import vector_metrics

SPLIT_SEED = 226
TRAINING_SEED = 42
SCENARIOS = ('cold-drug', 'cold-drug-dose')
METRICS = ('pcc', 'spearman', 'rmse', 'r2')


def make_splits(obs):
    """Assign each sample a scored, development or exclusion role in each fold."""
    parents = obs.drug_parent_inchikey.to_numpy(dtype=str)
    doses = obs.dose_um.to_numpy(dtype=float)
    unique = np.unique(parents)
    dose_levels = np.unique(doses)
    if obs.index.has_duplicates or obs.drug_parent_inchikey.isna().any():
        raise ValueError('missing or duplicate partition identities')
    if len(unique) < 3 or '' in unique:
        raise ValueError('at least three nonempty parent identities required')
    if len(dose_levels) != 3 or not np.isfinite(doses).all() or (doses <= 0).any():
        raise ValueError('exactly three positive finite doses required')

    rows = pd.DataFrame(dict(row=np.arange(len(obs)), profile_id=obs.index.astype(str),
                             drug_parent=parents, dose_um=doses))
    groups = np.array_split(np.random.default_rng(SPLIT_SEED).permutation(unique), 3)
    frames = []
    for scenario in SCENARIOS:
        for fold, held_parents in enumerate(groups):
            held_parent = np.isin(parents, held_parents)
            omitted_dose = None if scenario == 'cold-drug' else dose_levels[fold]
            held_dose = doses == omitted_dose
            development = ~held_parent & ~held_dose
            test = held_parent if omitted_dose is None else held_parent & held_dose
            units = np.unique(parents[development])
            inner_rng = np.random.default_rng(SPLIT_SEED + 1000 + fold)
            valid_parents = inner_rng.permutation(units)[:math.ceil(.1 * len(units))]
            validation = development & np.isin(parents, valid_parents)
            train = development & ~validation
            if not all(mask.any() for mask in (train, validation, test)):
                raise ValueError(f'empty partition: {scenario} fold {fold}')
            partitions = np.select(
                [train, validation, test, held_parent, held_dose],
                ['train', 'validation', 'test', 'excluded-held-parent', 'excluded-held-dose'],
                default='unassigned',
            )
            frames.append(rows.assign(scenario=scenario, fold=fold, partition=partitions))
    return pd.concat(frames, ignore_index=True)


def score_prediction(obs, truth, prediction, control, prediction_ids=None):
    """Compute Delta metrics independently on each sample and record invalid rows."""
    if prediction_ids is not None and list(map(str, prediction_ids)) != list(map(str, obs.index)):
        raise ValueError('prediction identity order mismatch')
    truth, prediction, control = map(np.asarray, (truth, prediction, control))
    if truth.ndim != 2 or truth.shape != prediction.shape or truth.shape != control.shape or len(obs) != len(truth):
        raise ValueError('prediction identity/shape mismatch')
    result = []
    for measured, predicted, matched_control in zip(truth, prediction, control, strict=True):
        delta = measured.astype(np.float64) - matched_control
        predicted_delta = predicted.astype(np.float64) - matched_control
        finite = np.isfinite(delta).all() and np.isfinite(predicted_delta).all()
        if finite:
            metrics = vector_metrics(delta, predicted_delta)
        else:
            metrics = dict.fromkeys(METRICS, float('nan'))
            metrics.update(constant_truth=False, constant_prediction=False)
        result.append(dict(metrics, nonfinite=not finite))
    identities = pd.DataFrame(dict(profile_id=obs.index.astype(str),
                                   drug_parent=obs.drug_parent_inchikey.to_numpy(),
                                   dose_um=obs.dose_um.to_numpy()))
    if 'cell_line' in obs:
        identities['cell_line'] = obs.cell_line.to_numpy()
    return identities.join(pd.DataFrame(result))


def summary(samples):
    """Pool all samples equally, recording unavailable metrics independently."""
    result = dict(n_samples=len(samples),
                  constant_truth_count=int(samples.constant_truth.sum()),
                  constant_prediction_count=int(samples.constant_prediction.sum()),
                  nonfinite_count=int(samples.nonfinite.sum()))
    for key in METRICS:
        finite = samples[key][np.isfinite(samples[key])]
        result[key] = float(finite.mean()) if len(finite) else None
        result[f'undefined_{key}_count'] = len(samples) - len(finite)
    return result
