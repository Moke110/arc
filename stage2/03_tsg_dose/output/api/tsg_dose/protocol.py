"""Frozen all-dose identities, group partitions and repeat-first Δ scoring."""
import math
import numpy as np
import pandas as pd
from . import shared
from tenv5k.protocol import vector_metrics
PROTOCOL_SCHEMA='tsg-dose-protocol-v1'
INPUT_SCHEMA='tsg-dose-inputs-v1'
FOLDS=3
DOSES=(0.2,1.,10.)
SEEDS=(17,42,2026)
METHODS=('no-dose','scaling','input-numeric','latent-numeric','input-embedding','latent-embedding','film')


def make_splits(obs,seed=226):
    """Freeze partitions from sorted parent/profile/dose identities only.

    Cold-drug uses the accepted02 parent shuffle/inner seeds. Final uses a
    separate seed2264000 parent shuffle. Held-dose validation greedily selects
    seeded parent–dose groups while retaining the parent's other training dose.
    """
    if obs.index.duplicated().any():raise ValueError('unique profile identities required')
    parents=obs.drug_parent_inchikey.astype(str).to_numpy()
    doses=obs.dose_um.to_numpy(dtype=float)
    if not np.isfinite(doses).all() or not set(doses)<=set(DOSES):raise ValueError('invalid dose population')
    unique=np.array(sorted(set(parents)))
    if len(unique)<3 or any(not p or p=='nan' for p in unique):raise ValueError('too few or missing parents')
    rows=np.arange(len(obs));result=[]
    def append(scenario,fold,train,valid,test):
        for partition, selected in [('train',train),('validation',valid),('test',test)]:
            result.extend(dict(scenario=scenario,fold=fold,partition=partition,row=int(r),profile_id=str(obs.index[r]),drug_parent=parents[r],dose_um=float(doses[r])) for r in sorted(selected))
    outer=np.empty(len(obs),dtype=int)
    for f,group in enumerate(np.array_split(np.random.default_rng(seed).permutation(unique),3)):
        outer[np.isin(parents,group)]=f
    for f in range(3):
        remain=rows[outer!=f]; units=np.array(sorted(set(parents[remain])))
        valid_parents=np.random.default_rng(seed+1000+f).permutation(units)[:math.ceil(.1*len(units))]
        valid=remain[np.isin(parents[remain],valid_parents)]
        train=np.setdiff1d(remain,valid)
        if not len(train) or not len(valid):raise ValueError('empty cold partition')
        append('cold-drug',f,train,valid,rows[outer==f])
    for f,dose in enumerate(DOSES):
        test=rows[doses==dose];remain=rows[doses!=dose]
        groups=sorted({(parents[r],float(doses[r])) for r in remain})
        target=math.ceil(.1*len(groups));selected=[];used=set()
        for index in np.random.default_rng(seed+2000+f).permutation(len(groups)):
            group=groups[index]
            if group[0] not in used:
                selected.append(group);used.add(group[0])
            if len(selected)==target:break
        if len(selected)!=target:raise ValueError('held-dose validation coverage impossible')
        selected=set(selected)
        valid=np.array([r for r in remain if (parents[r],float(doses[r])) in selected],dtype=int)
        train=np.setdiff1d(remain,valid)
        if not len(test) or not len(valid) or not set(parents[np.r_[test,valid]])<=set(parents[train]):raise ValueError('held-dose training-parent coverage impossible')
        append('held-dose',f,train,valid,test)
    valid_parents=np.random.default_rng(seed+4000).permutation(unique)[:math.ceil(.1*len(unique))]
    valid=rows[np.isin(parents,valid_parents)];train=np.setdiff1d(rows,valid)
    append('final',0,train,valid,[])
    return pd.DataFrame(result)


def score_prediction(obs,truth,prediction,control,prediction_ids=None):
    """Return well and parent–dose metrics; never average different doses' vectors."""
    if prediction_ids is not None and list(map(str,prediction_ids))!=list(map(str,obs.index)):raise ValueError('prediction identity order mismatch')
    truth=np.asarray(truth);prediction=np.asarray(prediction);control=np.asarray(control)
    if truth.shape!=prediction.shape or truth.shape!=control.shape or len(obs)!=len(truth):raise ValueError('prediction identity/shape mismatch')
    if not np.isfinite(control).all():raise ValueError('non-finite control')
    delta=truth-control;pred=prediction-control
    wells=[];groups=[]
    for i in range(len(obs)):
        wells.append(dict(profile_id=str(obs.index[i]),drug_parent=str(obs.iloc[i].drug_parent_inchikey),dose_um=float(obs.iloc[i].dose_um),**vector_metrics(delta[i],pred[i]),**{'expression_'+k:v for k,v in vector_metrics(truth[i],prediction[i]).items()}))
    for (parent,dose),indices in obs.groupby(['drug_parent_inchikey','dose_um']).indices.items():
        groups.append(dict(drug_parent=str(parent),dose_um=float(dose),n_wells=len(indices),**vector_metrics(delta[indices].mean(axis=0),pred[indices].mean(axis=0))))
    return pd.DataFrame(wells),pd.DataFrame(groups)


def parent_metrics(groups):
    """Equal dose weights per parent, with explicit undefined/constant counts."""
    rows=[]
    for parent, frame in groups.groupby('drug_parent',sort=True):
        rows.append(dict(drug_parent=parent,n_doses=len(frame),n_wells=int(frame.n_wells.sum()),**{k:float(frame[k].mean()) for k in ['pcc','spearman','rmse','r2']},truth_constant_count=int(frame.constant_truth.sum()),prediction_constant_count=int(frame.constant_prediction.sum()),undefined_pcc_count=int(frame.pcc.isna().sum()),undefined_r2_count=int(frame.r2.isna().sum())))
    return pd.DataFrame(rows)


def summary(groups):
    parents=parent_metrics(groups)
    return {k:float(parents[k].mean()) for k in ['pcc','spearman','rmse','r2']}
