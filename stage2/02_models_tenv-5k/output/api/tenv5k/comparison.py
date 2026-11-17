"""Parent-paired comparison of complete out-of-fold predictions (ADR 0074)."""
import itertools
import numpy as np
import pandas as pd

METRICS = ('pcc', 'spearman', 'rmse', 'r2')
CANDIDATES = ('PRnet', 'TranSiGen')
SEEDS = {17, 42, 2026}


def assemble_oof(population, folds):
    """Place each held-out well exactly once, independent of fold ordering."""
    population = list(map(str, population))
    if len(set(population)) != len(population): raise ValueError('duplicate population')
    lookup = dict(zip(population, range(len(population))))
    seen = set(); result = None
    for identities, prediction in folds:
        identities = list(map(str, identities)); prediction = np.asarray(prediction)
        if prediction.ndim != 2 or len(identities) != len(prediction) or not np.isfinite(prediction).all():
            raise ValueError('invalid OOF prediction')
        if len(set(identities)) != len(identities) or seen.intersection(identities) or not set(identities) <= set(population):
            raise ValueError('duplicate or unknown OOF well')
        if result is None: result = np.empty((len(population), prediction.shape[1]), dtype=prediction.dtype)
        if result.shape[1] != prediction.shape[1]: raise ValueError('OOF gene dimension mismatch')
        result[[lookup[p] for p in identities]] = prediction
        seen.update(identities)
    if seen != set(population): raise ValueError('incomplete OOF population')
    return result


def seed_average(rows):
    """Average three seed metric values per drug, never prediction vectors."""
    records = []
    for (scenario, method, parent), group in rows.groupby(['scenario','method','drug_parent'], sort=True):
        deterministic = method in ('Mean','Mean_cell','Mean_drug')
        expected = {-1} if deterministic else SEEDS
        if set(group.seed) != expected or len(group) != len(expected): raise ValueError('missing/duplicate seed metrics')
        if group.constant_truth.nunique() != 1 or group.n_wells.nunique() != 1: raise ValueError('seed population differs')
        constant = bool(group.constant_truth.iloc[0])
        record = dict(scenario=scenario, method=method, drug_parent=parent, n_wells=int(group.n_wells.iloc[0]),
                      n_seeds=len(expected), constant_truth=constant,
                      constant_prediction_count=int(group.constant_prediction.sum()))
        for metric in METRICS:
            values = group[metric].to_numpy(float)
            undefined = constant and metric in ('pcc','spearman','r2')
            if undefined:
                if not np.isnan(values).all(): raise ValueError('inconsistent undefined truth metric')
                record[metric] = float('nan')
            else:
                if not np.isfinite(values).all(): raise ValueError('non-finite seed metric')
                record[metric] = float(values.mean())
        records.append(record)
    return pd.DataFrame(records)


def paired_comparisons(parents, repetitions=10000, seed=2262026):
    """Shared drug resampling indices; percentile two-sided 95% intervals."""
    result=[]
    for scenario, rows in parents.groupby('scenario', sort=True):
        methods = sorted(rows.method.unique())
        populations = {m: rows[rows.method==m].set_index('drug_parent').sort_index() for m in methods}
        first=populations[methods[0]]
        for table in populations.values():
            if table.index.duplicated().any() or not table.index.equals(first.index) or not np.array_equal(table.constant_truth,first.constant_truth):
                raise ValueError('methods do not share parent population/truth mask')
        for metric in METRICS:
            mask = ~first.constant_truth.to_numpy(bool) if metric in ('pcc','spearman','r2') else np.ones(len(first),bool)
            n=int(mask.sum())
            if not n: raise ValueError('no defined parent metrics')
            values={m:t[metric].to_numpy(float)[mask] for m,t in populations.items()}
            if any(not np.isfinite(v).all() for v in values.values()): raise ValueError('non-finite paired metric')
            pairs=[(a,b) for a in CANDIDATES if a in methods for b in methods if b!=a]
            draws={pair:[] for pair in pairs}
            rng=np.random.default_rng(seed)
            for start in range(0,repetitions,128):
                indices=rng.integers(n,size=(min(128,repetitions-start),n))
                for pair in pairs:
                    draws[pair].append((values[pair[0]]-values[pair[1]])[indices].mean(axis=1))
            for a,b in pairs:
                low,high=np.quantile(np.concatenate(draws[(a,b)]),[.025,.975])
                result.append(dict(scenario=scenario,metric=metric,candidate=a,reference=b,
                    estimate=float((values[a]-values[b]).mean()),low=float(low),high=float(high),
                    n_parents=n,n_excluded_constant=int((~mask).sum()),repetitions=repetitions,statistical_seed=seed,
                    direction='candidate minus reference',zero=0.))
    return pd.DataFrame(result)


def select_winner(pairs):
    cold=pairs[(pairs.scenario=='cold-drug')&(pairs.metric=='pcc')]
    winners=[]
    for candidate in CANDIDATES:
        rows=cold[cold.candidate==candidate]
        required={'PRnet','TranSiGen','Mean','Mean_cell','MLP_Morgan'}-{candidate}
        if set(rows.reference)!=required: raise ValueError('incomplete primary comparisons')
        if (rows.low>0).all(): winners.append(candidate)
    if len(winners)>1: raise ValueError('contradictory winners')
    return dict(winner=winners[0] if winners else None, retained=winners or list(CANDIDATES),
                criterion='cold-drug parent Δ PCC; all paired 95% lower bounds strictly above zero')
