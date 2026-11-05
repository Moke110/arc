"""Shared safe predictor requests and bounded minibatches; no model imports."""
import numpy as np
import torch
from artifacts import rows

class PredictorInputs:
    """Open only the request allowlist, with bounded copied minibatches."""
    def __init__(self, request):
        if request.get('target_kind') != 'delta':
            raise ValueError('PRnet requires the frozen delta target contract')
        files = request['predictor_files']
        self.n_genes = request['n_genes']
        self.cases = {role:rows(files[role+'_cases']) for role in ('train','validation','test')}
        allowed = {'case_id','drug_key','dose_um','control_index','group_id'}
        for role, cases in self.cases.items():
            if not cases or any(set(case) != allowed for case in cases):
                raise ValueError('predictor case leakage/invalid schema')
            if len({case['case_id'] for case in cases}) != len(cases):
                raise ValueError('duplicate case identity')
        for left,right in [('train','validation'),('train','test'),('validation','test')]:
            if {r['drug_key'] for r in self.cases[left]} & {r['drug_key'] for r in self.cases[right]}:
                raise ValueError('drug leakage between roles')
        self.controls = np.load(files['controls'],mmap_mode='r',allow_pickle=False)
        self.morgan = np.load(files['morgan'],mmap_mode='r',allow_pickle=False)
        self.drug_index = {r['drug_key']:int(r['morgan_index']) for r in rows(files['drugs'])}
        self.delta = {role:np.load(files[role+'_expression'],mmap_mode='r',allow_pickle=False)
                      for role in ('train','validation')}
        if self.controls.ndim != 2 or self.controls.shape[1] != self.n_genes:
            raise ValueError('control gene shape mismatch')
        if self.morgan.ndim != 2 or self.morgan.shape[1] != 2048:
            raise ValueError('Morgan shape mismatch')
        for role,array in self.delta.items():
            if array.shape != (len(self.cases[role]),self.n_genes):
                raise ValueError('target shape mismatch')
        for cases in self.cases.values():
            for case in cases:
                if case['drug_key'] not in self.drug_index:
                    raise ValueError('unknown drug identity')
                if not 0 <= int(case['control_index']) < len(self.controls):
                    raise ValueError('unknown control identity')
                if not 0 <= self.drug_index[case['drug_key']] < len(self.morgan):
                    raise ValueError('unknown Morgan identity')
                if not np.isfinite(float(case['dose_um'])) or float(case['dose_um']) <= 0:
                    raise ValueError('invalid treated dose')

    def batch(self, role, indices, device):
        cases = [self.cases[role][int(i)] for i in indices]
        baseline = np.array(self.controls[[int(r['control_index']) for r in cases]],dtype='float32')
        fingerprint = np.array(self.morgan[[self.drug_index[r['drug_key']] for r in cases]],dtype='float32')
        dose = np.array([float(r['dose_um']) for r in cases],dtype='float32')
        arrays = [baseline,fingerprint,dose]
        if role in self.delta:
            arrays.append(np.array(self.delta[role][indices],dtype='float32') + baseline)
        if not all(np.isfinite(array).all() for array in arrays):
            raise ValueError('nonfinite predictor input')
        return [torch.from_numpy(array).to(device) for array in arrays]


def minibatches(count, batch_size, generator=None, training=False):
    indices = generator.permutation(count) if generator is not None else np.arange(count)
    if training and count < 2:
        raise ValueError('native PRnet BatchNorm needs at least two training cases')
    start = 0
    while start < count:
        end = min(start+batch_size,count)
        if training and count-end == 1:
            end = count
        yield indices[start:end]
        start = end

