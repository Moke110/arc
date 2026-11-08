"""Provide a single native sampler pool over separate, memory-mapped sources."""
import numpy as np
import torch

from artifacts import rows
from predictor_inputs import PredictorInputs


class SourceInputs(PredictorInputs):
    """Use the Step 03 batch assembly with only the roles admitted to Mix."""

    def __init__(self, files, roles, n_genes):
        self.cases = {role: rows(files[f'{role}_cases']) for role in roles}
        self.controls = np.load(files['controls'], mmap_mode='r', allow_pickle=False)
        self.morgan = np.load(files['morgan'], mmap_mode='r', allow_pickle=False)
        self.drug_index = {row['drug_key']: int(row['morgan_index']) for row in rows(files['drugs'])}
        self.delta = {role: np.load(files[f'{role}_expression'], mmap_mode='r', allow_pickle=False)
                      for role in roles if role != 'test'}
        if (self.controls.ndim != 2 or self.controls.shape[1] != n_genes
                or self.morgan.ndim != 2 or self.morgan.shape[1] != 2048
                or any(array.shape != (len(self.cases[role]), n_genes) or array.dtype != 'float32'
                       for role, array in self.delta.items())):
            raise ValueError('invalid pooled predictor array shapes')


class MixInputs:
    def __init__(self, request):
        self.n_genes = request['n_genes']
        self.tenv = SourceInputs(request['tenv'], ('train', 'validation', 'test'), self.n_genes)
        self.tahoe = SourceInputs(request['tahoe'], ('train',), self.n_genes)
        # Native conditioning uses drug keys to index a frozen Morgan table.
        # Prefixing Tahoe keys keeps each source's verified fingerprint intact.
        self.morgan = np.concatenate([self.tenv.morgan, self.tahoe.morgan])
        self.drug_index = {
            **self.tenv.drug_index,
            **{f'tahoe:{drug}': index + len(self.tenv.morgan)
               for drug, index in self.tahoe.drug_index.items()},
        }
        self.cases = {**self.tenv.cases, 'train': []}
        for name, source in (('tenv', self.tenv), ('tahoe', self.tahoe)):
            for case in source.cases['train']:
                self.cases['train'].append({
                    **case, 'case_id': f"{name}:{case['case_id']}", 'group_id': f"{name}:{case['group_id']}",
                    'drug_key': f"tahoe:{case['drug_key']}" if name == 'tahoe' else case['drug_key'],
                })

    def batch(self, role, indices, device):
        if role != 'train':
            return self.tenv.batch(role, indices, device)
        boundary = len(self.tenv.cases['train'])
        indices = np.asarray(indices)
        tenv_rows = indices < boundary
        if tenv_rows.all():
            return self.tenv.batch(role, indices, device)
        if not tenv_rows.any():
            return self.tahoe.batch(role, indices - boundary, device)
        tenv = self.tenv.batch(role, indices[tenv_rows], device)
        tahoe = self.tahoe.batch(role, indices[~tenv_rows] - boundary, device)
        positions = np.concatenate([np.flatnonzero(tenv_rows), np.flatnonzero(~tenv_rows)])
        order = torch.as_tensor(np.argsort(positions), device=device)
        return [torch.cat([target, source]).index_select(0, order) for target, source in zip(tenv, tahoe)]
