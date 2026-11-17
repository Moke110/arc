"""Train a native TranSiGen candidate, then replay it in a fresh process."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'output/api'))
from tenv5k.learner import train, accept_replay

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True, help='immutable TranSiGen source root')
    parser.add_argument('--scenario', choices=['cold-drug', 'warm-start'], default='cold-drug')
    parser.add_argument('--fold', type=int, default=0)
    parser.add_argument('--seed', type=int, default=17)
    parser.add_argument('--device', default='cpu')
    parser.add_argument('--smoke-rows', type=int)
    parser.add_argument('--max-epochs', type=int, default=1000)
    parser.add_argument('--registration', type=Path, help='matrix execution-attempt registration')
    parser.add_argument('--replay-only', action='store_true')
    args = parser.parse_args()
    if args.replay_only:
        print(json.dumps(accept_replay(args.inputs, args.candidate, args.source, adapter='transigen', registration=args.registration)))
    else:
        train(args.inputs, args.candidate, args.source, adapter='transigen', scenario=args.scenario, fold=args.fold,
              seed=args.seed, device=args.device, smoke_rows=args.smoke_rows, max_epochs=args.max_epochs, registration=args.registration)
        subprocess.run([sys.executable, str(Path(__file__).resolve()), '--inputs', str(args.inputs),
                        '--candidate', str(args.candidate), '--source', str(args.source), '--replay-only'] +
                       (['--registration', str(args.registration)] if args.registration is not None else []), check=True)
