#!/usr/bin/env python3
"""Replay real HVG isolation and export measurement-derived diagnostic QA only."""
import argparse
import json
from pathlib import Path
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases.hvg_qa import verify_hvg_metrics

STEP = Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--training',type=Path,default=STEP/'output/training-corpora')
    parser.add_argument('--hvg',type=Path,default=STEP/'output/tenv-hvg')
    parser.add_argument('--successor',type=Path,default=STEP/'output/tenv-truth-successor')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--record',type=Path)
    args=parser.parse_args()
    result=verify_hvg_metrics(args.training,args.hvg,args.successor,args.output,
        progress=lambda event:print(json.dumps(event,sort_keys=True),flush=True))
    value=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.record:
        args.record.parent.mkdir(parents=True,exist_ok=True)
        args.record.write_text(value)
    print(value,end='')

if __name__=='__main__':
    main()
