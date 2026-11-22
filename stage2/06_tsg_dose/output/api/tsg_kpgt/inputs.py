"""tenv-only KPGT transport over the accepted axis-local input population."""

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

from tsg_round1.inputs import Inputs as FrozenInputs, sha256

STEP = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(STEP.parent / '05_tsg_drug_emb/output/api'))
from tsg_drug_emb.representations import PARENT_COLUMNS, WIDTHS, checked_vectors

INPUT_SCHEMA = 'tsg-dose-tenv-kpgt-inputs-v1'
WIDTH = WIDTHS['KPGT']
ACCEPTED_TENV = 'e08ea61360eab0b3b0d1a4762652d9070acda3cbe00347285c1a62e6266a427e'
ACCEPTED_REPRESENTATIONS = (
    'bc10c6af897b1b0d939a685aa198e03b22a3392e3891b721f70cc49f1f028678'
)


def align_vectors(expected, supplied, values):
    """Join clean, pre-canonicalized parent identities without transforming vectors."""
    expected = expected[PARENT_COLUMNS]
    supplied = supplied[PARENT_COLUMNS]
    key = PARENT_COLUMNS[0]
    for parents in (expected, supplied):
        if (
            parents.isna().any().any()
            or parents[key].duplicated().any()
            or any(
                not parents[column].str.strip().ne('').all()
                for column in PARENT_COLUMNS
            )
        ):
            raise ValueError('KPGT parent identity has missing or duplicate values')
    supplied = supplied.set_index(key)
    wanted = expected[key].tolist()
    if set(wanted) != set(supplied.index):
        missing = sorted(set(wanted) - set(supplied.index))
        extra = sorted(set(supplied.index) - set(wanted))
        raise ValueError(
            f'KPGT parent coverage mismatch: missing {missing}, extra {extra}'
        )
    aligned = supplied.loc[wanted, PARENT_COLUMNS[1]].tolist()
    mismatches = [
        key
        for key, a, b in zip(wanted, expected[PARENT_COLUMNS[1]], aligned, strict=True)
        if a != b
    ]
    if mismatches:
        raise ValueError(f'KPGT parent canonical SMILES mismatch: {mismatches}')
    values = checked_vectors('KPGT', values, len(supplied))
    return values[supplied.index.get_indexer(wanted)]


class Inputs(FrozenInputs):
    def __init__(self, directory):
        directory = Path(directory)
        successor = json.loads((directory / 'manifest.json').read_text())
        actual = {
            'tenv/manifest.json': sha256(directory / 'tenv/manifest.json'),
            'KPGT.npy': sha256(directory / 'KPGT.npy'),
            'representation.json': sha256(directory / 'representation.json'),
            'representation_parents.csv': sha256(
                directory / 'representation_parents.csv'
            ),
        }
        if (
            successor.get('schema') != INPUT_SCHEMA
            or successor.get('files') != actual
            or actual['tenv/manifest.json'] != ACCEPTED_TENV
            or actual['representation.json'] != ACCEPTED_REPRESENTATIONS
        ):
            raise ValueError('tenv KPGT manifest identity or coverage mismatch')
        super().__init__(directory / 'tenv')
        if self.source != 'tenv':
            raise ValueError('KPGT inputs require accepted tenv')
        self.directory = directory
        self.manifest = successor
        representation = self.read_json('representation.json')
        if (
            representation['inputs_sha256'] != ACCEPTED_TENV
            or representation['files']['KPGT.npy'] != actual['KPGT.npy']
            or representation['files']['parents.csv']
            != actual['representation_parents.csv']
        ):
            raise ValueError(
                'KPGT artifact or parent roster differs from accepted Step 03'
            )
        self.kpgt = align_vectors(
            self.parents,
            pd.read_csv(directory / 'representation_parents.csv'),
            np.load(directory / 'KPGT.npy', mmap_mode='r'),
        )

    def fingerprints(self, rows):
        positions = self.obs.iloc[rows].drug_parent_inchikey.map(self.parent_lookup)
        return self.kpgt[positions.to_numpy()]
