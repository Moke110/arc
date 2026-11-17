"""Verified, bounded shared inputs for every single-dose learner."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .protocol import PROTOCOL_SCHEMA, INPUT_SCHEMA, FOLDS, make_splits


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


class Inputs:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.manifest = json.loads((self.directory / 'manifest.json').read_text())
        for name, expected in self.manifest['outputs'].items():
            if Path(name).name != name or sha256(self.directory / name) != expected:
                raise ValueError('input digest mismatch: ' + name)
        self.config = json.loads((self.directory / 'config.json').read_text())
        if (self.manifest.get('schema') != INPUT_SCHEMA or
                self.config.get('schema') != PROTOCOL_SCHEMA or
                self.config.get('folds') != FOLDS or self.config.get('split_seed') != 226):
            raise ValueError('input protocol differs from approved three-fold comparison')
        self.obs = pd.read_csv(self.directory / 'population.csv').set_index('profile_id', drop=False)
        self.all_obs = pd.read_csv(self.directory / 'all_obs.csv').set_index('profile_id', drop=False)
        self.parents = pd.read_csv(self.directory / 'parents.csv')
        self.genes = pd.read_csv(self.directory / 'genes.csv')
        self.splits = pd.read_csv(self.directory / 'splits.csv')
        self.controls = {int(k): v for k, v in json.loads((self.directory / 'matched_controls.json').read_text()).items()}
        for name in ('expression', 'control', 'delta', 'morgan'):
            setattr(self, name, np.load(self.directory / (name + '.npy'), mmap_mode='r'))
        if self.obs.index.duplicated().any() or list(self.all_obs.index[:len(self.obs)]) != list(self.obs.index):
            raise ValueError('population identity mismatch')
        if self.expression.shape != self.control.shape or self.delta.shape != self.expression.shape or self.expression.shape != (len(self.all_obs), len(self.genes)):
            raise ValueError('expression shape mismatch')
        if self.morgan.shape != (len(self.parents), 2048):
            raise ValueError('Morgan identity shape mismatch')
        expected_splits = make_splits(self.obs, seed=self.config['split_seed'])
        if not self.splits.equals(expected_splits):
            raise ValueError('split protocol or complete OOF coverage mismatch')
        for row in self.splits.itertuples():
            if row.profile_id != self.obs.index[row.row] or row.drug_parent != self.obs.iloc[row.row].drug_parent_inchikey:
                raise ValueError('split identity mismatch')

    def rows(self, scenario, fold, partition):
        selected = self.splits
        return selected.loc[(selected.scenario == scenario) & (selected.fold == fold) & (selected.partition == partition), 'row'].to_numpy(dtype=int)

    def fingerprints(self, rows):
        lookup = dict(zip(self.parents.drug_parent_inchikey, range(len(self.parents))))
        return np.asarray(self.morgan[[lookup[p] for p in self.obs.iloc[rows].drug_parent_inchikey]], dtype=np.float32)
