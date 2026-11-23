"""Validation-first dose-aware native fitting with exact fresh-process replay."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'output/api'))
from tsg_dose.learner import train,accept_replay,evaluate
from tsg_dose.protocol import METHODS

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path,required=True);parser.add_argument('--candidate',type=Path,required=True)
    parser.add_argument('--source',type=Path,required=True);parser.add_argument('--method',choices=METHODS,default='no-dose')
    parser.add_argument('--scenario',choices=['cold-drug','held-dose','final'],default='cold-drug')
    parser.add_argument('--fold',type=int,default=0);parser.add_argument('--seed',type=int,default=17);parser.add_argument('--device',default='cpu')
    parser.add_argument('--smoke-rows',type=int);parser.add_argument('--max-epochs',type=int,default=1000)
    parser.add_argument('--registration',type=Path);parser.add_argument('--selection',type=Path)
    modes=parser.add_mutually_exclusive_group();modes.add_argument('--replay-only',action='store_true');modes.add_argument('--evaluate-only',action='store_true')
    args=parser.parse_args()
    if args.replay_only:print(json.dumps(accept_replay(args.inputs,args.candidate,args.source)))
    else:
        if args.evaluate_only:evaluate(args.inputs,args.candidate,args.source,args.selection)
        else:
            registration=json.loads(args.registration.read_text()) if args.registration else None
            train(args.inputs,args.candidate,args.source,method=args.method,scenario=args.scenario,fold=args.fold,seed=args.seed,device=args.device,smoke_rows=args.smoke_rows,max_epochs=args.max_epochs,registration=registration)
        subprocess.run([sys.executable,str(Path(__file__).resolve()),'--inputs',str(args.inputs),'--candidate',str(args.candidate),'--source',str(args.source),'--replay-only'],check=True)
