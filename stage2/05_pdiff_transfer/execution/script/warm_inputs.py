"""Predictor input for the warm-start parent–dose partition."""
from collections import defaultdict

import numpy as np

from artifacts import rows
from predictor_inputs import PredictorInputs


class WarmInputs(PredictorInputs):
    """Load only predictor covariates and the two permitted expression roles."""

    def __init__(self, request):
        self.n_genes = request['n_genes']
        files = request['predictor_files']
        self.cases = {role: rows(files[f'{role}_cases']) for role in ('train', 'validation', 'test')}
        drugs = rows(files['drugs'])
        self.drug_index = {row['drug_key']: int(row['morgan_index']) for row in drugs}
        parents = {row['drug_key']: row['parent_connectivity'] for row in drugs}
        self.controls = np.load(files['controls'], mmap_mode='r', allow_pickle=False)
        self.morgan = np.load(files['morgan'], mmap_mode='r', allow_pickle=False)
        self.delta = {role: np.load(files[f'{role}_expression'], mmap_mode='r', allow_pickle=False)
                      for role in ('train', 'validation')}
        if (request['target_kind'] != 'delta'
                or self.controls.ndim != 2 or self.controls.shape[1] != self.n_genes
                or self.morgan.ndim != 2 or self.morgan.shape[1] != 2048
                or any(array.shape != (len(self.cases[role]), self.n_genes) or array.dtype != 'float32'
                       for role, array in self.delta.items())):
            raise ValueError('invalid warm-start predictor shape or target')
        roles = defaultdict(set)
        seen = set()
        for role, cases in self.cases.items():
            if not cases:
                raise ValueError('empty warm-start predictor role')
            for case in cases:
                if set(case) != {'case_id', 'drug_key', 'dose_um', 'control_index', 'group_id'}:
                    raise ValueError('unexpected warm-start predictor case fields')
                drug, dose = case['drug_key'], float(case['dose_um'])
                if (not case['case_id'] or case['case_id'] in seen or not np.isfinite(dose) or dose <= 0
                        or not 0 <= int(case['control_index']) < len(self.controls)
                        or not 0 <= self.drug_index[drug] < len(self.morgan)):
                    raise ValueError('invalid or repeated warm-start predictor case')
                seen.add(case['case_id'])
                roles[parents[drug], dose].add(role)
        if any(len(assigned) != 1 for assigned in roles.values()):
            raise ValueError('warm-start parent–dose leakage')
