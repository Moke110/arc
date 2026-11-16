"""Shared single-dose split, baseline and scoring boundaries (ADR 0074)."""
from __future__ import annotations
import math
import numpy as np
import pandas as pd
from scipy.stats import rankdata


def mean_prediction(x, obs, controls, train, test, method, scenario):
    """Fit solely on training treated rows and unique matched original controls."""
    train=np.asarray(train,dtype=int); test=np.asarray(test,dtype=int)
    if set(train)&set(test):
        raise ValueError('training and prediction rows overlap')
    if method=='Mean_drug':
        if scenario=='cold-drug':
            return None
        result=[]
        for row in test:
            same=train[obs.iloc[train].drug_parent_inchikey.to_numpy()==obs.iloc[row].drug_parent_inchikey]
            if not len(same):
                raise ValueError('warm-start drug lacks training coverage')
            result.append(np.asarray(x[same],dtype=np.float64).mean(axis=0))
        return np.asarray(result,dtype=np.float32)
    if method not in ('Mean','Mean_cell'):
        raise ValueError('unknown Mean method')
    result=[]
    cached={}
    for row in test:
        key = str(obs.iloc[row].cell_line) if method == 'Mean_cell' else 'all'
        if key in cached:
            result.append(cached[key])
            continue
        selected=train
        if method=='Mean_cell':
            selected=train[obs.iloc[train].cell_line.to_numpy()==obs.iloc[row].cell_line]
        if not len(selected):
            raise ValueError('no training cell context')
        control_rows=sorted({int(c) for t in selected for c in controls[int(t)]})
        fit=np.concatenate([selected,np.asarray(control_rows,dtype=int)])
        if not np.isfinite(x[fit]).all():
            raise ValueError('non-finite fitting expression')
        cached[key] = np.asarray(x[fit],dtype=np.float64).mean(axis=0)
        result.append(cached[key])
    return np.asarray(result,dtype=np.float32)


def make_splits(obs, seed=226, folds=5):
    """Frozen integer rule: ceil(10% outer-training units); balanced outer folds.

    Cold units are sorted drug parents, shuffled once. Warm units are wells:
    each parent's independently shuffled wells are round-robin assigned with a
    random starting fold. Validation greedily selects shuffled wells while
    retaining one training well per parent; impossible coverage fails explicitly.
    """
    if obs.index.duplicated().any() or len(obs)<folds:
        raise ValueError('unique wells and sufficient population required')
    parents=obs.drug_parent_inchikey.astype(str).to_numpy()
    if any(not p for p in parents):
        raise ValueError('missing parent')
    unique=np.asarray(sorted(set(parents)))
    if len(unique)<folds:
        raise ValueError('too few parents')
    result=[]
    for scenario in ('cold-drug','warm-start'):
        rng=np.random.default_rng(seed)
        outer=np.empty(len(obs),dtype=int)
        if scenario=='cold-drug':
            groups=np.array_split(rng.permutation(unique),folds)
            for fold,group in enumerate(groups): outer[np.isin(parents,group)]=fold
        else:
            for parent in unique:
                rows=np.flatnonzero(parents==parent)
                if len(rows)<2: raise ValueError('warm-start requires at least two wells per parent')
                rows=rng.permutation(rows); offset=int(rng.integers(folds))
                outer[rows]=(np.arange(len(rows))+offset)%folds
        for fold in range(folds):
            test=np.flatnonzero(outer==fold); remain=np.flatnonzero(outer!=fold)
            vrng=np.random.default_rng(seed+1000+fold)
            if scenario=='cold-drug':
                units=np.asarray(sorted(set(parents[remain])))
                n=math.ceil(len(units)*0.1)
                val_parents=vrng.permutation(units)[:n]
                val=remain[np.isin(parents[remain],val_parents)]
            else:
                n=math.ceil(len(remain)*0.1); counts={p:int((parents[remain]==p).sum()) for p in unique}
                val=[]
                for row in vrng.permutation(remain):
                    if counts[parents[row]]>1:
                        val.append(int(row)); counts[parents[row]]-=1
                    if len(val)==n: break
                if len(val)!=n: raise ValueError('warm-start validation coverage impossible')
                val=np.asarray(val,dtype=int)
            train=np.setdiff1d(remain,val)
            if not len(test) or not len(val) or not len(train): raise ValueError('empty partition')
            if scenario=='warm-start' and not set(parents[np.r_[test,val]])<=set(parents[train]):
                raise ValueError('warm-start coverage impossible')
            for role,rows in [('train',train),('validation',val),('test',test)]:
                result.extend({'scenario':scenario,'fold':fold,'partition':role,'row':int(r),'profile_id':str(obs.index[r]),'drug_parent':parents[r]} for r in sorted(rows))
    return pd.DataFrame(result)


def vector_metrics(truth,prediction):
    """Constants share a truth mask; prediction constants score zero correlation."""
    truth=np.asarray(truth,dtype=np.float64); prediction=np.asarray(prediction,dtype=np.float64)
    if truth.shape!=prediction.shape or not np.isfinite(truth).all() or not np.isfinite(prediction).all():
        raise ValueError('shape mismatch or non-finite prediction/truth')
    constant_truth=bool(np.ptp(truth)==0); constant_prediction=bool(np.ptp(prediction)==0)
    mse=float(np.mean((truth-prediction)**2)); denom=float(np.sum((truth-truth.mean())**2))
    if constant_truth: pcc=spearman=float('nan')
    elif constant_prediction: pcc=spearman=0.0
    else:
        pcc=float(np.corrcoef(truth,prediction)[0,1]);spearman=float(np.corrcoef(rankdata(truth),rankdata(prediction))[0,1])
    return {'pcc':pcc,'spearman':spearman,'rmse':math.sqrt(mse),'r2':1-float(np.sum((truth-prediction)**2))/denom if denom else float('nan'),'constant_truth':constant_truth,'constant_prediction':constant_prediction}


def score_prediction(obs, truth, prediction, control, prediction_ids=None):
    """Return per-well and parent-level Δ metrics plus expression auxiliaries."""
    if prediction_ids is not None and list(map(str, prediction_ids)) != list(map(str, obs.index)):
        raise ValueError('prediction identity order mismatch')
    truth=np.asarray(truth); prediction=np.asarray(prediction);control=np.asarray(control)
    if truth.shape!=prediction.shape or control.shape!=truth.shape or len(obs)!=len(truth):
        raise ValueError('prediction identity/shape mismatch')
    rows=[]
    for i in range(len(obs)):
        rows.append({'profile_id':str(obs.index[i]),'drug_parent':str(obs.iloc[i].drug_parent_inchikey),**vector_metrics(truth[i]-control[i],prediction[i]-control[i]),**{'expression_'+k:v for k,v in vector_metrics(truth[i],prediction[i]).items()}})
    groups=[]
    for parent in sorted(set(obs.drug_parent_inchikey)):
        mask=obs.drug_parent_inchikey.to_numpy()==parent
        groups.append({'drug_parent':parent,'n_wells':int(mask.sum()),**vector_metrics((truth[mask]-control[mask]).mean(axis=0),(prediction[mask]-control[mask]).mean(axis=0))})
    return pd.DataFrame(rows),pd.DataFrame(groups)
