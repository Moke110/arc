"""Validation-only choice and parent-block inference for dose-input comparison."""
import itertools
import numpy as np
import pandas as pd
from .protocol import METHODS,SEEDS
from . import shared
from tenv5k.comparison import assemble_oof
METRICS=('pcc','spearman','rmse','r2')


def choose_method(rows):
    """Mean three seed metrics; exact ties use parameter count then frozen order."""
    if len(rows)!=21 or set(zip(rows.method,rows.seed))!={(m,s) for m in METHODS for s in SEEDS}:raise ValueError('choice requires21 unique method/seed scores')
    candidates=[]
    for method in METHODS:
        frame=rows[rows.method==method]
        if frame.parameter_count.nunique()!=1 or not np.isfinite(frame.pcc).all():raise ValueError('invalid validation scores/parameters')
        candidates.append((-float(frame.pcc.mean()),int(frame.parameter_count.iloc[0]),METHODS.index(method),method))
    return min(candidates)[-1]


def seed_average(rows):
    """Average complete seed metrics per parent, retaining dose/constant counts."""
    result=[]
    keys=['method','drug_parent']
    for identity,frame in rows.groupby(keys,sort=True):
        if len(frame)!=3 or set(frame.seed)!=set(SEEDS):raise ValueError('incomplete/duplicate seed metrics')
        record=dict(zip(keys,identity),n_seeds=3)
        for name in ['n_doses','n_wells','truth_constant_count','undefined_pcc_count','undefined_r2_count']:
            if name in frame:
                if frame[name].nunique()!=1:raise ValueError('seed populations differ')
                record[name]=int(frame[name].iloc[0])
        if 'prediction_constant_count' in frame:record['prediction_constant_count']=int(frame.prediction_constant_count.sum())
        for metric in METRICS:
            values=frame[metric].to_numpy(float)
            if np.isinf(values).any() or (np.isnan(values).any() and not np.isnan(values).all()):raise ValueError('inconsistent undefined seed metric')
            record[metric]=float(values.mean())
        result.append(record)
    return pd.DataFrame(result)


def paired_comparisons(parents,repetitions=10000,seed=2262026):
    """Resample parents carrying every dose and seed together, choices fixed."""
    if repetitions<1:raise ValueError('positive bootstrap count required')
    methods=sorted(parents.method.unique())
    tables={m:parents[parents.method==m].set_index('drug_parent').sort_index() for m in methods}
    if not methods:raise ValueError('empty comparison')
    first=tables[methods[0]]
    for frame in tables.values():
        if frame.index.duplicated().any() or not frame.index.equals(first.index):raise ValueError('unpaired parent populations')
    pairs=[('selected','no-dose')] if set(methods)=={'selected','no-dose'} else list(itertools.permutations(methods,2))
    result=[]
    for metric in METRICS:
        arrays={m:t[metric].to_numpy(float) for m,t in tables.items()}
        if any(np.isinf(a).any() for a in arrays.values()):raise ValueError('infinite paired metric')
        mask=np.isfinite(arrays[methods[0]])
        if any(not np.array_equal(np.isfinite(a),mask) for a in arrays.values()):raise ValueError('methods disagree on undefined parents')
        n=int(mask.sum())
        if not n:raise ValueError('no defined parents for '+metric)
        differences={pair:(arrays[pair[0]]-arrays[pair[1]])[mask] for pair in pairs}
        draws={pair:[] for pair in pairs};rng=np.random.default_rng(seed)
        for start in range(0,repetitions,128):
            index=rng.integers(n,size=(min(128,repetitions-start),n))
            for pair,values in differences.items():draws[pair].append(values[index].mean(axis=1))
        for pair,values in differences.items():
            low,high=np.quantile(np.concatenate(draws[pair]),[.025,.975])
            result.append(dict(candidate=pair[0],reference=pair[1],metric=metric,estimate=float(values.mean()),low=float(low),high=float(high),n_parents=n,n_excluded_undefined=len(first)-n,repetitions=repetitions,statistical_seed=seed,direction='candidate minus reference',zero=0.,inference='fixed-selection conditional parent-paired bootstrap; does not repeat method selection',role='primary' if pair==('selected','no-dose') else 'descriptive'))
    return pd.DataFrame(result)
