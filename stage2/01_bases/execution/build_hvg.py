#!/usr/bin/env python3
"""Publish or replay standalone descriptive and train-only tenv HVG rankings."""
import argparse
import json
from pathlib import Path
import pandas as pd
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases.hvg import build_hvg, verify_hvg

STEP = Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--training',type=Path,default=STEP/'output/training-corpora/tenv.h5ad')
    parser.add_argument('--output',type=Path,default=STEP/'output/tenv-hvg')
    parser.add_argument('--cohorts',type=Path,help='Explicit complete sample_id/fold/role CSV; otherwise new deterministic cold-parent folds')
    parser.add_argument('--folds',type=int,default=3)
    parser.add_argument('--seed',default='bases-cold-parent-v1')
    parser.add_argument('--verify',action='store_true')
    parser.add_argument('--record',type=Path)
    args=parser.parse_args()
    result=verify_hvg(args.training,args.output) if args.verify else build_hvg(
        args.training,args.output,cohorts=pd.read_csv(args.cohorts,keep_default_na=False,dtype={'sample_id':str,'drug_parent':str}) if args.cohorts else None,
        folds=args.folds,seed=args.seed)
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.record:
        args.record.parent.mkdir(parents=True,exist_ok=True);args.record.write_text(text)
    print(text,end='')

if __name__=='__main__':
    main()
