#!/usr/bin/env python3
"""Score indexed well predictions against only an authorized tenv partition."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from bases.hvg import sha256
from bases.metrics import score_groups
from bases.plots import render_retrieval
from bases.significance import PartitionTruth

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--predictions',type=Path,required=True,help='npy rows match samples CSV; gene columns match genes CSV')
    parser.add_argument('--samples',type=Path,required=True,help='sample_id,drug_parent,dose_um in prediction row order')
    parser.add_argument('--genes',type=Path,required=True,help='CSV gene_index column in prediction column order')
    parser.add_argument('--truth',type=Path,required=True)
    parser.add_argument('--assignments',type=Path,required=True,help='Complete group-isolated sample_id,role partition CSV')
    parser.add_argument('--fold',type=int)
    parser.add_argument('--role',choices=['train','validation','test'],required=True)
    parser.add_argument('--label',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    assignments=pd.read_csv(args.assignments,keep_default_na=False,dtype={'sample_id':str})
    if 'fold' in assignments:
        if args.fold is None:
            parser.error('--fold is required for a multi-fold assignment file')
        assignments=assignments[assignments.fold.eq(args.fold)]
    samples=pd.read_csv(args.samples,keep_default_na=False,dtype={'sample_id':str,'drug_parent':str})
    expected=assignments.loc[assignments.role.eq(args.role), 'sample_id']
    if (assignments.sample_id.duplicated().any() or samples.sample_id.duplicated().any()
            or expected.empty or set(samples.sample_id) != set(expected)):
        raise ValueError('prediction samples must exactly cover the requested evaluator role once')
    truth=PartitionTruth(args.truth,assignments,args.role)
    predictions=np.load(args.predictions,mmap_mode='r',allow_pickle=False)
    columns=pd.read_csv(args.genes).gene_index.to_numpy()
    result=score_groups(predictions,samples,truth,columns)
    args.output.mkdir(parents=True,exist_ok=False)
    result['groups'].to_csv(args.output/'groups.csv',index=False)
    render_retrieval([{'label':args.label,'summary':result['summary'],'curves':result['curves']}],
                     args.output/'pr.png')
    report={'schema':'stage2-indexed-prediction-evaluation-v1','source':'tenv',
            'role':args.role,'fold':args.fold,'label':args.label,'summary':result['summary'],
            'inputs':{name:{'path':str(path.resolve()),'sha256':sha256(path)} for name,path in (
                ('predictions',args.predictions),('samples',args.samples),('genes',args.genes),
                ('assignments',args.assignments),('truth_manifest',args.truth/'manifest.json'))},
            'groups_sha256':sha256(args.output/'groups.csv')}
    (args.output/'report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(report,indent=2,sort_keys=True))

if __name__=='__main__':
    main()
