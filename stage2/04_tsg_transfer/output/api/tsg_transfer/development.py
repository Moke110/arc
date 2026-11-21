"""Bounded development fits on disjoint parent groups."""
import hashlib
import math

import numpy as np
import torch

from . import fit
from tsg_round1 import model as models


def split_validation(data, train, validation):
    """Split frozen validation parents into selection and independent assessment."""
    train = np.asarray(train, dtype=int)
    validation = np.asarray(validation, dtype=int)
    parents = data.obs.drug_parent_inchikey.to_numpy()
    order = sorted(set(parents[validation]), key=lambda key: hashlib.sha256(key.encode()).digest())
    if len(order) < 2 or set(parents[train]).intersection(order):
        raise ValueError('expected at least two validation parents disjoint from training')
    selection = np.isin(parents[validation], order[::2])
    return dict(train=train, selection=validation[selection], assessment=validation[~selection])


def delta_scores(truth, prediction, control):
    """Return per-profile gene-wise Delta PCC and Delta RMSE."""
    truth, prediction, control = (np.asarray(item, dtype=np.float64) for item in (truth, prediction, control))
    if truth.ndim != 2 or truth.shape != prediction.shape or truth.shape != control.shape:
        raise ValueError('Delta score matrices must have the same two-dimensional shape')
    measured = truth - control
    predicted = prediction - control
    rmse = np.sqrt(np.square(measured - predicted).mean(axis=1))
    measured -= measured.mean(axis=1, keepdims=True)
    predicted -= predicted.mean(axis=1, keepdims=True)
    norm = np.linalg.norm(measured, axis=1) * np.linalg.norm(predicted, axis=1)
    pcc = np.divide(np.einsum('ij,ij->i', measured, predicted), norm,
                    out=np.full(len(truth), np.nan), where=norm > 0)
    return np.clip(pcc, -1, 1), rmse


def evaluate(native, model, data, rows):
    """Score rows at a fixed native sampling seed and retain profile identities."""
    prediction = fit.predict(native, model, data, rows, seed=fit.VALIDATION_SEED)
    pcc, rmse = delta_scores(data.expression[rows], prediction, data.control[rows])
    scores = data.obs.iloc[rows][['drug_parent_inchikey', 'dose_um']].copy()
    return scores.assign(pcc=pcc, rmse=rmse)


def summarize(scores):
    """Weight each defined parent equally after averaging its profiles."""
    by_parent = scores.groupby('drug_parent_inchikey')[['pcc', 'rmse']].mean()
    return dict(pcc=by_parent.pcc.mean(), rmse=by_parent.rmse.mean(),
                n_parents=len(by_parent), n_pcc_parents=int(by_parent.pcc.notna().sum()),
                n_profiles=len(scores), n_pcc_profiles=int(scores.pcc.notna().sum()))


def run_fit(data, source, roles, *, learning_rate, seed, epochs, initialize=None, device='cuda',
            assessment_evaluator=evaluate):
    """Train to a common budget, independently assess both validation-selected epochs."""
    train, selection, assessment = (np.asarray(roles[name], dtype=int)
                                    for name in ('train', 'selection', 'assessment'))
    parents = data.obs.drug_parent_inchikey.to_numpy()
    groups = [set(parents[rows]) for rows in (train, selection, assessment)]
    if not all(groups) or sum(map(len, groups)) != len(set.union(*groups)):
        raise ValueError('train, selection and assessment parents must be nonempty and disjoint')
    if not 1 <= epochs < 1000 or not math.isfinite(learning_rate) or learning_rate <= 0:
        raise ValueError('invalid bounded training budget or learning rate')

    fit._backend(device)
    transform = models.fit_transform(fit.METHOD, data.doses(train))
    native = models.adapter(source, fit.METHOD, transform)
    model = native.build(seed, device)
    if initialize is not None:
        initialize(model, transform)
    optimizer = native.optimizer(model)
    for group in optimizer.param_groups:
        group['lr'] = learning_rate

    shuffle = np.random.default_rng(seed)
    updates = math.ceil(len(train) / native.batch_size) - (len(train) % native.batch_size == 1)
    history, chosen = [], {}
    for epoch in range(1, epochs + 1):
        loss = native.train_epoch(model, data, train, optimizer, shuffle)
        scores = evaluate(native, model, data, selection)
        metrics = summarize(scores)
        row = dict(epoch=epoch, updates=epoch * updates, loss=loss, **metrics)
        if not all(math.isfinite(row[name]) for name in ('loss', 'pcc', 'rmse')):
            raise ValueError('development loss or selection score is undefined')
        # Equal-sized profiles make mean squared profile RMSE equal treated-expression MSE.
        row['mse'] = float(np.square(scores.rmse).mean())
        history.append(row)
        for criterion, value in (('mse', row['mse']), ('pcc', -row['pcc'])):
            if criterion not in chosen or value < chosen[criterion][0]:
                state = {name: tensor.detach().to('cpu', copy=True)
                         for name, tensor in model.state_dict().items()}
                chosen[criterion] = (value, row, state)

    result = {}
    for criterion, (_, selected, state) in chosen.items():
        model.load_state_dict(state)
        scores = assessment_evaluator(native, model, data, assessment)
        result[criterion] = dict(selected=selected, assessment=summarize(scores), scores=scores)
    return dict(history=history, selected=result, transform=transform)


def load_tahoe_winner(inputs, candidate, source):
    """Load and verify the frozen, source-only winner without changing fitting code."""
    from . import tahoe

    data = fit.load_inputs(inputs, 'tahoe_bulk')
    winner = tahoe.lock_winner(data, candidate)['winner']
    train, validation = tahoe.tahoe_rows(data, tahoe.prepare_tahoe(data, candidate), winner['fold'])
    _, model, _ = fit.load_fit(data, source, candidate / f"fold-{winner['fold']}",
                               train_rows=train, validation_rows=validation)
    return model, winner, tuple(data.obs.dose_um.unique())
