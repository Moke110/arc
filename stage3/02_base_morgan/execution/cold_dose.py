"""Deterministic globally held-dose partitions, grouped by parent and dose."""
from collections import Counter
import math
import random

DOSES = [0.04, 0.2, 1.0, 5.0]


def cold_dose(rows):
    from workflow import SPLIT_SEED, CLASSES
    if sorted({r['dose_um'] for r in rows}) != DOSES:
        raise ValueError('cold_dose requires all four dose levels 0.04, 0.2, 1, 5 µM; restore dose coverage')
    combinations = sorted({(r['parent_id'], r['dose_um']) for r in rows})
    targets = sorted({(r['condition'], r['phenotype']) for r in rows})
    splits = []
    for fold, held in enumerate(DOSES):
        available = [key for key in combinations if key[1] != held]
        # Coverage is required independently for each target as well as globally.
        coverage = Counter((r['condition'], r['phenotype'], r['parent_id'])
                           for r in { (r['condition'],r['phenotype'],r['parent_id'],r['dose_um']):r
                                      for r in rows if r['dose_um'] != held }.values())
        required = {(r['condition'],r['phenotype'],r['parent_id']) for r in rows}
        if any(coverage[key] == 0 for key in required):
            raise ValueError(f'cold_dose fold {fold}: parent coverage impossible; add another dose for every test parent and target')
        members = {key: {(r['condition'],r['phenotype'],r['parent_id']) for r in rows
                         if (r['parent_id'],r['dose_um']) == key} for key in available}
        candidates = list(available)
        random.Random(SPLIT_SEED + fold + 1).shuffle(candidates)
        count = math.ceil(0.2 * len(available))
        validation = set()
        for key in candidates:
            if all(coverage[m] > 1 for m in members[key]):
                validation.add(key)
                coverage.subtract(members[key])
                if len(validation) == count:
                    break
        validation_targets = {(r['condition'], r['phenotype']) for r in rows
                              if (r['parent_id'], r['dose_um']) in validation}
        if len(validation) != count or validation_targets != set(targets):
            validation = _validation(available, members, required, targets, count, fold)
        mapping = {key:'test' if key[1] == held else 'validation' if key in validation else 'train'
                   for key in combinations}
        roles = [mapping[(r['parent_id'],r['dose_um'])] for r in rows]
        for target in targets:
            if {role for r,role in zip(rows,roles) if (r['condition'],r['phenotype'])==target} != {'train','validation','test'}:
                raise ValueError(f'cold_dose fold {fold}: cannot populate all roles for {target}; add parent-dose coverage')
        counts = {}
        for role in ['train','validation','test']:
            selected = [r for r,a in zip(rows,roles) if a==role]
            counts[role] = dict(measurements=len(selected),parents=len({r['parent_id'] for r in selected}),
                                combinations=len({(r['parent_id'],r['dose_um']) for r in selected}),
                                classes={str(c):sum(r['direction']==c for r in selected) for c in CLASSES})
        splits.append(dict(setting='cold_dose',fold=fold,held_dose_um=held,
                           dose_generalization='extrapolation' if fold in (0,3) else 'interpolation',
                           combination_roles=[dict(parent_id=p,dose_um=d,role=mapping[(p,d)]) for p,d in combinations],
                           roles=roles,counts=counts))
    return splits


def _validation(available, members, required, targets, count, fold):
    """Resolve sparse coverage exactly when the deterministic greedy path fails."""
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import coo_matrix
    rr, cc, vv, lower, upper = [], [], [], [], []
    def constraint(indices, low, high):
        k = len(lower)
        for i in indices:
            rr.append(k); cc.append(i); vv.append(1.)
        lower.append(low); upper.append(high)
    constraint(range(len(available)), count, count)
    for member in sorted(required):
        indices = [i for i,key in enumerate(available) if member in members[key]]
        constraint(indices, 0, len(indices)-1)
    for target in targets:
        indices = [i for i,key in enumerate(available) if any(m[:2] == target for m in members[key])]
        constraint(indices, 1, len(indices)-1)
    matrix = coo_matrix((vv,(rr,cc)),shape=(len(lower),len(available))).tocsc()
    rng = random.Random(20261105 + fold)
    result = milp(np.array([rng.random() for _ in available]),
                  integrality=np.ones(len(available)), bounds=Bounds(0,1),
                  constraints=LinearConstraint(matrix,lower,upper), options={'time_limit':60})
    if result.x is None:
        reason = 'impossible' if result.status == 2 else 'unresolved within solver budget'
        raise ValueError(f'cold_dose fold {fold}: {reason} validation parent coverage; add dose coverage; {result.message}')
    return {key for i,key in enumerate(available) if result.x[i] > .5}
