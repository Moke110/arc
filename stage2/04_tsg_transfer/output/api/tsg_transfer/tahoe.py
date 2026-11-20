"""Tahoe-only condition validation and strict lowest-MSE checkpoint selection."""
from collections import Counter

import numpy as np
import pandas as pd

from tsg_round1.protocol import SPLIT_SEED


def make_tahoe_splits(obs):
    """Assign each parent × cell line × dose group to at most one validation fold."""
    columns = ['drug_parent_inchikey', 'cell_line', 'dose_um']
    dose = obs.dose_um.to_numpy(dtype=float)
    if (obs.index.has_duplicates or obs[columns].isna().any().any()
            or obs[columns[:2]].eq('').any().any()
            or not np.isfinite(dose).all() or (dose <= 0).any()):
        raise ValueError('invalid Tahoe condition or profile identity')
    conditions = list(obs[columns].itertuples(index=False, name=None))
    groups = sorted(set(conditions))
    marginals = [Counter(values) for values in zip(*groups, strict=True)]
    support = [[counts.copy() for counts in marginals] for _ in range(3)]
    sizes = Counter(conditions)
    fold_sizes = [0, 0, 0]
    assignments = {}
    rng = np.random.default_rng(SPLIT_SEED)
    for offset in rng.permutation(len(groups)):
        key = groups[offset]
        eligible = [fold for fold in range(3)
                    if all(counts[value] > 1
                           for counts, value in zip(support[fold], key, strict=True))]
        if not eligible:
            continue
        fold = min(eligible, key=lambda f: (fold_sizes[f], f))
        assignments[key] = fold
        fold_sizes[fold] += sizes[key]
        for counts, value in zip(support[fold], key, strict=True):
            counts[value] -= 1
    if not all(fold_sizes):
        raise ValueError('cannot produce three nonempty Tahoe validation folds')
    assigned = np.array([assignments.get(key, -1) for key in conditions])
    rows = pd.DataFrame(dict(profile_id=obs.index.astype(str), row=np.arange(len(obs)),
                             drug_parent=obs.drug_parent_inchikey.to_numpy(),
                             cell_line=obs.cell_line.to_numpy(), dose_um=dose))
    rows['validation_coverage'] = np.where(assigned == -1, 'training-only', 'assigned')
    return pd.concat([rows.assign(fold=fold, partition=np.where(assigned == fold,
                                                                'validation', 'train'))
                      for fold in range(3)], ignore_index=True)


def prepare_tahoe(data, candidate):
    """Freeze Tahoe folds and account for conditions with no validation exposure."""
    from pathlib import Path
    from .fit import read_json, sha256, write_json

    candidate = Path(candidate)
    if data.source != 'tahoe_bulk':
        raise ValueError('Tahoe preparation requires Tahoe inputs')
    splits = make_tahoe_splits(data.obs)
    candidate.mkdir(parents=True, exist_ok=True)
    path = candidate / 'splits.csv'
    receipt_path = candidate / 'splits.json'
    if not path.exists() and not receipt_path.exists():
        splits.to_csv(path, index=False)
        training_only = splits.loc[(splits.fold == 0) & (splits.validation_coverage == 'training-only')]
        write_json(receipt_path, dict(
            schema='stage2-tsg-transfer-tahoe-splits-v1',
            input_manifest_sha256=sha256(data.directory / 'manifest.json'),
            split_seed=SPLIT_SEED, splits_sha256=sha256(path),
            role_counts=[dict(fold=int(fold), partition=partition, rows=int(count))
                         for (fold, partition), count in
                         splits.groupby(['fold', 'partition']).size().items()],
            training_only_profiles=len(training_only),
            training_only_conditions=len(training_only[
                ['drug_parent', 'cell_line', 'dose_um']].drop_duplicates()),
        ))
    else:
        receipt = read_json(receipt_path)
        if receipt['input_manifest_sha256'] != sha256(data.directory / 'manifest.json'):
            raise ValueError('Tahoe split input identity differs')
        if sha256(path) != receipt['splits_sha256'] or not splits.equals(pd.read_csv(path)):
            raise ValueError('Tahoe split content differs from reconstructed assignment')
    return splits


def tahoe_rows(data, splits, fold, *, smoke_rows=None):
    """Resolve one preregistered training/validation fold and optional smoke slice."""
    if fold not in range(3):
        raise ValueError('unregistered Tahoe fold')
    selected = splits.loc[splits.fold == fold]
    train, validation = (selected.loc[selected.partition == role, 'row'].to_numpy(dtype=int)
                         for role in ('train', 'validation'))
    if smoke_rows is not None:
        from tsg_round1.learner import smoke_subset
        if smoke_rows < 4:
            raise ValueError('Tahoe smoke requires at least four rows')
        return smoke_subset(data, train, smoke_rows), smoke_subset(data, validation, smoke_rows)
    return train, validation


def lock_winner(data, candidate, *, smoke_rows=None):
    """Revalidate all three Tahoe fold receipts and preserve one earliest-fold winner."""
    from pathlib import Path
    from .fit import read_json, sha256, verify_replay, write_json

    candidate = Path(candidate)
    splits = prepare_tahoe(data, candidate)
    candidates = []
    for fold in range(3):
        path = candidate / f'fold-{fold}'
        train, validation = tahoe_rows(data, splits, fold, smoke_rows=smoke_rows)
        fit = verify_replay(data, path, train_rows=train, validation_rows=validation,
                            kind='tahoe', fold=fold, formal=smoke_rows is None)
        candidates.append(dict(
            fold=fold, best_epoch=fit['best_epoch'], validation_mse=fit['validation_mse'],
            transform=fit['transform'], checkpoint_sha256=fit['files']['checkpoint.pt'],
            fit_sha256=sha256(path / 'fit.json'),
            replay_sha256=sha256(path / 'validation-replay.json'),
        ))
    winner = min(candidates, key=lambda fit: (fit['validation_mse'], fit['fold']))
    record = dict(schema='stage2-tsg-transfer-tahoe-winner-v1', smoke_rows=smoke_rows,
                  input_manifest_sha256=sha256(data.directory / 'manifest.json'),
                  split_receipt_sha256=sha256(candidate / 'splits.json'),
                  candidates=candidates, winner=winner)
    path = candidate / 'winner.json'
    if path.exists():
        if read_json(path) != record:
            raise ValueError('locked Tahoe winner differs from candidate receipts')
    else:
        write_json(path, record)
    return record
