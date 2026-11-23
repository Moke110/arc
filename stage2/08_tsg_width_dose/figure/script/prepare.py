"""Build plot tables from a verified Step 05 formal run matrix."""

import argparse
from collections import defaultdict
from itertools import product
import math
from pathlib import Path
from statistics import mean, pstdev, stdev
import sys

import pandas as pd

from plot_contract import BAR, CURVES, FOLDS, METHODS, METRICS, OVERALL, PER_DOSE
from plot_contract import WIDTHS, read_rows, write_rows

STAGE = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(STAGE / '07_tsg_dose/output/api'))
from tsg_round1.protocol import summary

IDENTITY_COLUMNS = ['profile_id', 'drug_parent', 'dose_um']
DOT_COLUMNS = [f'fold_{fold}' for fold in FOLDS]
MIN_DIAMETER, MAX_DIAMETER = 12, 30
ZERO_Z_DIAMETER = 21
POINTS_PER_Z = 4


def scaled_diameter(fraction):
    return MIN_DIAMETER + (MAX_DIAMETER - MIN_DIAMETER) * fraction


def fold_scores(path, reference):
    rows = pd.read_csv(
        path / 'test_per_sample.csv',
        dtype={'profile_id': str, 'drug_parent': str},
        float_precision='round_trip',
    )
    identities = rows[IDENTITY_COLUMNS]
    doses = rows.dose_um.unique()
    if (
        len(doses) != 1
        or not math.isfinite(doses[0])
        or doses[0] <= 0
        or identities.isna().any().any()
        or rows.profile_id.duplicated().any()
        or rows.nonfinite.any()
    ):
        raise ValueError(f'invalid held-dose scores: {path}')
    if reference is not None and not identities.equals(reference):
        raise ValueError(f'held-dose identities differ across variants: {path}')
    scores = summary(rows)
    return doses[0], identities, {metric: scores[metric] for metric in METRICS}


def curves(path, coordinate):
    samples = read_rows(path, 'learning_curve')
    if not samples:
        raise ValueError(f'empty learning curve: {path}')
    rows = []
    for epoch, sample in enumerate(samples, start=1):
        train_loss = float(sample['train_loss'])
        validation_mse = float(sample['validation_mse'])
        if int(sample['epoch']) != epoch or not all(map(math.isfinite, (train_loss, validation_mse))):
            raise ValueError(f'invalid learning curve: {path}')
        rows.append(dict(coordinate, epoch=epoch, train_loss=train_loss, validation_mse=validation_mse))
    best = min(rows, key=lambda row: row['validation_mse'])['epoch']
    return [dict(row, best_validation=int(row['epoch'] == best)) for row in rows]


def summarize(run_doses):
    by_variant = defaultdict(list)
    by_metric = defaultdict(list)
    for row in run_doses:
        metric = row['metric']
        by_variant[row['width'], row['method'], metric].append(row['mean'])
        if row['mean'] is not None:
            by_metric[metric].append(row['mean'])
    stats = {
        metric: (mean(values), pstdev(values), min(values), max(values))
        for metric, values in by_metric.items()
    }

    balloons, bars = [], []
    for (width, method, metric), dots in by_variant.items():
        defined = [value for value in dots if value is not None]
        average = mean(defined) if defined else None
        coordinate = dict(width=width, method=method, metric=metric, mean=average)
        bars.append(dict(
            coordinate,
            sd=stdev(defined) if len(defined) > 1 else None,
            **dict(zip(DOT_COLUMNS, dots, strict=True)),
        ))
        z_score = diameter = None
        if average is not None:
            center, spread, _, _ = stats[metric]
            z_score = (average - center) / spread if spread else 0.0
            diameter = min(MAX_DIAMETER, max(MIN_DIAMETER, ZERO_Z_DIAMETER + POINTS_PER_Z * z_score))
        balloons.append(dict(coordinate, z_score=z_score, diameter_pt=diameter))

    for row in run_doses:
        value = row['mean']
        if value is None:
            row['diameter_pt'] = None
            continue
        _, _, low, high = stats[row['metric']]
        fraction = (value - low) / (high - low) if high != low else 0.5
        row['diameter_pt'] = scaled_diameter(fraction)

    return balloons, bars


def prepare(candidate, output):
    """Require all 48 score and history files, and matching test rows per fold."""
    references = {}
    run_doses, learning = [], []
    for width, method, fold in product(WIDTHS, METHODS, FOLDS):
        path = candidate / 'runs' / width / method / f'fold-{fold}'
        dose, identities, scores = fold_scores(path, references.get(fold))
        references[fold] = identities
        coordinate = dict(width=width, method=method, fold=fold)
        run_doses.extend(
            dict(coordinate, dose_um=dose, metric=metric, mean=value)
            for metric, value in scores.items()
        )
        learning.extend(curves(path, coordinate))

    held_doses = [frame.dose_um.iloc[0] for frame in references.values()]
    if len(set(held_doses)) != len(FOLDS):
        raise ValueError('each fold must hold a different dose')
    balloons, bars = summarize(run_doses)
    output.mkdir(parents=True, exist_ok=True)
    for name, rows in ((OVERALL, balloons), (BAR, bars), (PER_DOSE, run_doses), (CURVES, learning)):
        write_rows(output, name, rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    prepare(args.candidate, args.output)


if __name__ == '__main__':
    main()
