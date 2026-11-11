#!/usr/bin/env python3
"""Score indexed well predictions against only an authorized tenv partition."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases.hvg import sha256
from bases.metrics import score_groups
from bases.plots import render_retrieval
from bases.significance import PartitionTruth

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--predictions',type=Path,required=True,help='npy rows match samples CSV; gene columns match genes CSV')
    parser.add_argument('--samples',type=Path,required=True,help='sample_id,drug_parent,dose_um in prediction row order')
    parser.add_argument('--genes',type=Path,required=True,help='CSV gene_index column in prediction column order')
    truth_source=parser.add_mutually_exclusive_group(required=True)
    truth_source.add_argument('--truth',type=Path,help='Original accepted truth release')
    truth_source.add_argument('--successor',type=Path,help='Accepted exact successor envelope')
    parser.add_argument('--training',type=Path,help='Successor training release directory; required with --successor')
    parser.add_argument('--diagnostic',action='store_true',help='Mark measurement-derived QA; never model performance')
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
    if args.successor:
        if args.training is None:
            parser.error('--training is required with --successor')
        from bases.truth_successor import SuccessorTruth
        truth=SuccessorTruth(args.successor,args.training,assignments,args.role)
    else:
        if args.training is not None:
            parser.error('--training requires --successor')
        truth=PartitionTruth(args.truth,assignments,args.role)
    predictions=np.load(args.predictions,mmap_mode='r',allow_pickle=False)
    columns=pd.read_csv(args.genes).gene_index.to_numpy()
    result=score_groups(predictions,samples,truth,columns)
    args.output.mkdir(parents=True,exist_ok=False)
    result['groups'].to_csv(args.output/'groups.csv',index=False)
    label=('Diagnostic: '+args.label) if args.diagnostic else args.label
    render_retrieval([{'label':label,'summary':result['summary'],'curves':result['curves']}],
                     args.output/'pr.png',diagnostic=args.diagnostic)
    report={'schema':'stage2-indexed-prediction-evaluation-v1','source':'tenv',
            'role':args.role,'fold':args.fold,'label':label,'summary':result['summary'],
            'diagnostic':args.diagnostic,'interpretation':('Measurement-derived QA only; no model evaluated' if args.diagnostic else 'Prediction evaluation'),
            'inputs':{name:{'path':str(path.resolve()),'sha256':sha256(path)} for name,path in (
                ('predictions',args.predictions),('samples',args.samples),('genes',args.genes),
                ('assignments',args.assignments),('truth_manifest',truth.root/'manifest.json'))},
            'groups_sha256':sha256(args.output/'groups.csv')}
    if args.successor:
        report['inputs']['successor_manifest']={'path':str((args.successor/'manifest.json').resolve()),
            'sha256':truth.successor_manifest_sha256}
        report['inputs']['training_manifest']={'path':str((args.training/'manifest.json').resolve()),
            'sha256':sha256(args.training/'manifest.json')}
    (args.output/'report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(report,indent=2,sort_keys=True))

if __name__=='__main__':
    main()
