#!/usr/bin/env python3
"""Publish or replay standalone descriptive and train-only tenv HVG rankings.

``--training`` names the live raw corpus payload the current pipeline builds
(``execution/build_corpus.py`` produces ``output/corpus/tenv.h5ad`` on the frozen
18,164-gene axis). The pre-0066 ``output/training-corpora/tenv.h5ad`` is retired
and is never read. The builder still consumes a training-style h5ad (obs carrying
``sample_id``, ``drug_parent`` and ``dose_um``) and the accepted
``output/tenv-hvg/`` bytes stay pinned, so deriving a training-style projection
from the corpus is the documented deferred re-point (ADR 0066 *Deferred*): the
default only stops naming a retired payload.
"""
import argparse
import json
from pathlib import Path
import pandas as pd
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases.hvg import build_hvg, verify_hvg

STEP = Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--training',type=Path,default=STEP/'output/corpus/tenv.h5ad')
    parser.add_argument('--output',type=Path,default=STEP/'output/tenv-hvg')
    parser.add_argument('--cohorts',type=Path,help='Explicit complete sample_id/fold/role CSV; otherwise new deterministic cold-parent folds')
    parser.add_argument('--folds',type=int,default=3)
    parser.add_argument('--seed',default='bases-cold-parent-v1')
    parser.add_argument('--verify',action='store_true')
    parser.add_argument('--record',type=Path)
    args=parser.parse_args()
    if not args.training.is_file():
        parser.error(f'--training payload not found: {args.training} '
                     '(publish the live corpus first: python stage2/01_bases/execution/build_corpus.py publish)')
    result=verify_hvg(args.training,args.output) if args.verify else build_hvg(
        args.training,args.output,cohorts=pd.read_csv(args.cohorts,keep_default_na=False,dtype={'sample_id':str,'drug_parent':str}) if args.cohorts else None,
        folds=args.folds,seed=args.seed)
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.record:
        args.record.parent.mkdir(parents=True,exist_ok=True);args.record.write_text(text)
    print(text,end='')

if __name__=='__main__':
    main()
