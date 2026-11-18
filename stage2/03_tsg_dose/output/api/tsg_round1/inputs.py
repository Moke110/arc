"""Verified inputs for one source of the ADR0077 dual-source comparison."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from tsg_dose import shared
from tenv5k.inputs import sha256

from .protocol import make_splits

INPUT_SCHEMA = 'tsg-dose-round1-axis-local-inputs-v1'
SOURCE_POPULATIONS = {
    'tenv': (4206, 459, 1, (.2, 1., 10.)),
    'tahoe_bulk': (65318, 352, 50, (.05, .5, 5.)),
}
REQUIRED = frozenset(('config.json', 'population.csv', 'all_obs.csv', 'parents.csv',
                      'genes.csv', 'splits.csv', 'matched_controls.json',
                      'expression.npy', 'control.npy', 'morgan.npy'))


def validate_population(source, obs, axis_size):
    if source not in SOURCE_POPULATIONS:
        raise ValueError('unknown round1 input source')
    count, parents, lines, doses = SOURCE_POPULATIONS[source]
    observed = (len(obs), obs.drug_parent_inchikey.nunique(),
                obs.cell_line.nunique(), tuple(sorted(obs.dose_um.unique())))
    if observed != (count, parents, lines, doses) or axis_size != 5000:
        raise ValueError(f'{source}: unexpected formal treated population or gene axis {observed}')


class Inputs:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.manifest = self.read_json('manifest.json')
        if self.manifest.get('schema') != INPUT_SCHEMA or set(self.manifest.get('outputs', {})) != REQUIRED:
            raise ValueError('round1 axis-local input manifest schema or coverage mismatch')
        for name, digest in self.manifest['outputs'].items():
            if sha256(self.directory / name) != digest:
                raise ValueError(f'round1 input digest mismatch: {name}')
        self.config = self.read_json('config.json')
        self.source = self.config['source']
        generation = self.config['view_generation']
        if (generation != self.manifest.get('view_generation')
                or generation.get('schema') != 'stage2-view-v3'
                or generation.get('axis', {}).get('size') != 5000
                or generation.get('normalization', {}).get('source') != 'ADR 0077 five-k-only axis-local normalization'
                or set(generation.get('payload_sha256', {})) != set(SOURCE_POPULATIONS)):
            raise ValueError('round1 axis-local view generation mismatch')

        self.obs = pd.read_csv(self.directory / 'population.csv').set_index('profile_id', drop=False)
        self.all_obs = pd.read_csv(self.directory / 'all_obs.csv', low_memory=False).set_index(
            'profile_id', drop=False)
        self.genes = pd.read_csv(self.directory / 'genes.csv')
        self.parents = pd.read_csv(self.directory / 'parents.csv')
        self.splits = pd.read_csv(self.directory / 'splits.csv')
        self.controls = {int(key): value for key, value in self.read_json('matched_controls.json').items()}
        for name in ('expression', 'control', 'morgan'):
            setattr(self, name, np.load(self.directory / f'{name}.npy', mmap_mode='r'))
        validate_population(self.source, self.obs, len(self.genes))
        keys = ['profile_id', 'view_row', 'role', 'drug_parent_inchikey',
                'dose_um', 'cell_line', 'matched_control_group']
        if (self.all_obs.index.has_duplicates
                or not self.all_obs.iloc[:len(self.obs)][keys].equals(self.obs[keys])):
            raise ValueError('round1 profile identity mismatch')
        vehicles = self.all_obs.iloc[len(self.obs):]
        vehicle_groups = {str(group): (offsets + len(self.obs)).tolist()
                          for group, offsets in vehicles.groupby('matched_control_group').indices.items()}
        if (not vehicles.role.eq('vehicle').all()
                or set(self.controls) != set(range(len(self.obs)))
                or any(self.controls[i] != vehicle_groups.get(str(group))
                       for i, group in enumerate(self.obs.matched_control_group))):
            raise ValueError('round1 matched-control identity mismatch')
        shape = (len(self.all_obs), len(self.genes))
        for matrix in (self.expression, self.control):
            if matrix.shape != shape or matrix.dtype != np.float32:
                raise ValueError('round1 expression shape or dtype mismatch')
        if (self.morgan.shape != (len(self.parents), 2048) or self.morgan.dtype != np.uint8
                or self.parents.drug_parent_inchikey.duplicated().any()
                or set(self.parents.drug_parent_inchikey) != set(self.obs.drug_parent_inchikey)):
            raise ValueError('round1 Morgan identity mismatch')
        if not self.splits.equals(make_splits(self.obs)):
            raise ValueError('round1 partition identity or coverage mismatch')
        self.parent_lookup = dict(zip(self.parents.drug_parent_inchikey, range(len(self.parents))))

    def read_json(self, name):
        return json.loads((self.directory / name).read_text())

    def rows(self, scenario, fold, partition):
        s = self.splits
        return s.loc[(s.scenario == scenario) & (s.fold == fold)
                     & (s.partition == partition), 'row'].to_numpy(dtype=int)

    def fingerprints(self, rows):
        positions = self.obs.iloc[rows].drug_parent_inchikey.map(self.parent_lookup)
        return np.asarray(self.morgan[positions.to_numpy()], dtype=np.float32)

    def doses(self, rows):
        return self.obs.iloc[rows].dose_um.to_numpy(dtype=np.float32)
