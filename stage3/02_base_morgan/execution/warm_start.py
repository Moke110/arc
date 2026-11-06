"""Deterministic warm-start allocation of indivisible parent–dose combinations."""
from collections import Counter
import math
import random
import warnings

SPLIT_SEED = 20261105
ROLES = ('train', 'validation', 'test')


def warm_start(rows):
    """Test each combination once, retaining both marginal inputs in training.

    Seeded constrained allocation handles the normal dense dose grid cheaply.
    A binary feasibility solver resolves sparse graphs that defeat greedy choices;
    it never drops records or relaxes coverage to make a partition fit.
    """
    combinations = sorted({(r['parent_id'], r['dose_um']) for r in rows})
    indices = {key: i for i, key in enumerate(combinations)}
    if not combinations:
        raise ValueError('warm_start requires measured parent–dose combinations')
    target_keys = sorted({(r['condition'], r['phenotype']) for r in rows})
    targets = [{indices[r['parent_id'], r['dose_um']] for r in rows
                if (r['condition'], r['phenotype']) == target} for target in target_keys]
    groups = []
    for target_key, target in zip(target_keys, targets):
        for axis, label in [(0, 'parent'), (1, 'dose')]:
            marginal = {}
            for i in sorted(target):
                marginal.setdefault(combinations[i][axis], set()).add(i)
            for value, members in sorted(marginal.items()):
                if len(members) < 2:
                    raise ValueError(f'warm_start impossible {label} coverage for {value} '
                                     f'in target {target_key}: only one combination; '
                                     'add another measured combination for this marginal input')
                groups.append(members)
    # Index the constraints once: validation selection should scale with the
    # groups containing a combination, rather than all parents in the release.
    memberships = [[] for _ in combinations]
    for group_index, group in enumerate(groups):
        for i in group:
            memberships[i].append(group_index)
    parents = {}
    for i, (parent, dose) in enumerate(combinations):
        parents.setdefault(parent, []).append(i)

    def valid(assignments):
        for assigned in assignments:
            for group in groups:
                if not any(assigned[i] == 'train' for i in group):
                    return False
            if any({assigned[i] for i in target} != set(ROLES) for target in targets):
                return False
        return True

    assignments = None
    for attempt in range(128):
        rng = random.Random(SPLIT_SEED + attempt)
        order = list(range(len(combinations)))
        rng.shuffle(order)
        # Stagger each parent's seeded dose order across the three folds. A
        # four-dose parent contributes at most two test combinations per fold,
        # leaving room for validation while retaining training coverage.
        outer = [set() for _ in range(3)]
        for offset, members in enumerate(parents.values()):
            shuffled = list(members)
            rng.shuffle(shuffled)
            for position, i in enumerate(shuffled):
                outer[(offset + position) % 3].add(i)
        candidate = []
        for fold in range(3):
            test = outer[fold]
            available = [i for i in order if i not in test]
            rng.shuffle(available)
            assigned = ['test' if i in test else 'train' for i in range(len(combinations))]
            remaining = [sum(assigned[i] == 'train' for i in group) for group in groups]
            needed = math.ceil(.2 * len(available))
            for i in available:
                if needed and all(remaining[g] > 1 for g in memberships[i]):
                    assigned[i] = 'validation'
                    needed -= 1
                    for g in memberships[i]:
                        remaining[g] -= 1
            if needed:
                break
            candidate.append(assigned)
        if len(candidate) == 3 and valid(candidate):
            assignments = candidate
            break
    if assignments is None:
        assignments = _solve(combinations, groups, targets)
    if not valid(assignments):
        raise ValueError('warm_start coverage audit failed; inspect parent/dose and target coverage')
    splits = []
    for fold, assigned in enumerate(assignments):
        roles = [assigned[indices[r['parent_id'], r['dose_um']]] for r in rows]
        counts = {}
        for role in ROLES:
            selected = [r for r, actual in zip(rows, roles) if role == actual]
            classes = Counter(r['direction'] for r in selected)
            counts[role] = dict(measurements=len(selected), parents=len({r['parent_id'] for r in selected}),
                                combinations=len({(r['parent_id'], r['dose_um']) for r in selected}),
                                classes={str(c): classes[c] for c in [-1, 0, 1]})
        splits.append(dict(setting='warm_start', fold=fold, roles=roles, counts=counts,
                           combination_roles=[dict(parent_id=p, dose_um=d, role=assigned[i])
                                              for i, (p, d) in enumerate(combinations)]))
    return splits


def _solve(combinations, groups, targets):
    """Exact sparse binary constraints, including ceil(20% non-test) validation."""
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import coo_matrix

    n = len(combinations)
    # For each fold, binary test and validation; training is their complement.
    def variable(fold, role, i):
        return (fold * 2 + role) * n + i
    rr, cc, vv, lower, upper = [], [], [], [], []

    def constraint(terms, lo, hi):
        row = len(lower)
        for column, value in terms:
            rr.append(row)
            cc.append(column)
            vv.append(value)
        lower.append(lo)
        upper.append(hi)

    for i in range(n):
        constraint([(variable(f, 0, i), 1) for f in range(3)], 1, 1)
    for fold in range(3):
        for i in range(n):
            constraint([(variable(fold, role, i), 1) for role in range(2)], 0, 1)
        # n <= 5*validation + test <= n+4 is exactly validation=ceil((n-test)/5).
        constraint([(variable(fold, role, i), 1 if role == 0 else 5)
                    for role in range(2) for i in range(n)], n, n + 4)
        for group in groups:
            constraint([(variable(fold, role, i), 1) for role in range(2) for i in sorted(group)],
                       0, len(group) - 1)
        for target in targets:
            for role in range(2):
                constraint([(variable(fold, role, i), 1) for i in sorted(target)], 1, len(target))
            constraint([(variable(fold, role, i), 1) for role in range(2) for i in sorted(target)],
                       0, len(target) - 1)
    matrix = coo_matrix((vv, (rr, cc)), shape=(len(lower), 6 * n)).tocsc()
    rng = random.Random(SPLIT_SEED)
    objective = np.array([rng.random() for _ in range(6 * n)])
    # scipy forwards HiGHS' thread limit; keep fallback inside the shared CPU budget.
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message="Unrecognized options detected.*", category=RuntimeWarning)
        result = milp(objective, integrality=np.ones(6 * n), bounds=Bounds(0, 1),
                      constraints=LinearConstraint(matrix, lower, upper),
                      options={'time_limit': 60, 'threads': 1})
    if result.x is None:
        reason = 'impossible' if result.status == 2 else 'unresolved within solver budget'
        raise ValueError(f'warm_start {reason} parent/dose/target coverage with 20% validation; '
                         'add measured combinations spanning sparse parents and doses; '
                         f'solver diagnostic: {result.message}')
    return [['test' if result.x[variable(fold, 0, i)] > .5 else
             'validation' if result.x[variable(fold, 1, i)] > .5 else 'train'
             for i in range(n)] for fold in range(3)]
