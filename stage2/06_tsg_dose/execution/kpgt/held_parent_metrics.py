"""Aggregate accepted held-dose profile scores into one dot per drug parent."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

METRICS = ('pcc', 'spearman', 'rmse', 'r2')
SCENARIO = 'cold-drug-dose'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def aggregate(run, manifest_digest):
    manifest = run / 'manifest.json'
    scores = run / 'test_per_sample.csv'
    if digest(manifest) != manifest_digest:
        raise ValueError(f'Run manifest differs from acceptance: {run}')
    expected = json.loads(manifest.read_text())['outputs'][scores.name]
    if digest(scores) != expected:
        raise ValueError(f'Profile scores differ from run manifest: {run}')
    columns = ['drug_parent', *METRICS]
    frame = pd.read_csv(scores, usecols=columns, float_precision='round_trip')
    frame[list(METRICS)] = frame[list(METRICS)].replace([np.inf, -np.inf], np.nan)
    return frame.groupby('drug_parent', sort=True)[list(METRICS)].mean().reset_index()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-root', type=Path, required=True)
    parser.add_argument('--acceptance', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    receipt = json.loads(args.acceptance.read_text())
    if receipt['status'] != 'passed':
        raise ValueError('Scientific acceptance required')
    parent_rows = []
    for family, records, manifest_key in (
        ('formal/runs/tenv', receipt['origins'], 'run_manifest_sha256'),
        ('baselines/runs', receipt['baselines'], 'manifest_sha256'),
    ):
        for record in records:
            origin = record['origin']
            if origin['scenario'] != SCENARIO:
                continue
            fold, method = origin['fold'], origin['method']
            run = args.run_root / family / SCENARIO / f'fold-{fold}' / method
            parent_rows.append(aggregate(run, record[manifest_key]).assign(fold=fold, method=method))
    result = pd.concat(parent_rows, ignore_index=True)
    result[['fold', 'method', 'drug_parent', *METRICS]].to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
