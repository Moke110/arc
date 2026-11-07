#!/usr/bin/env python3
"""Step 04 PerturbDiff training artifact CLI."""
import argparse
import os
from pathlib import Path
import sys

from training_policy import MONITOR, Selection, native_batch_mean, native_epoch_mean

REPO = next(path for path in Path(__file__).resolve().parents if (path / '.git').exists())
STEP03 = Path('stage2/03_modes_tenv')


def replay_monitor(trace, output):
    import numpy as np
    from artifacts import read_json, write_json

    trace = read_json(trace)
    policy = Selection(trace['minimum_steps'], trace['patience_steps'], trace['min_delta'])
    history = []
    reason = 'trace_exhausted'
    for event in trace['events']:
        batches = [native_batch_mean(np.asarray(batch['loss'], dtype='float32'),
                                     np.asarray(batch['weights'], dtype='float32'))
                   for batch in event['batches']]
        history.append(policy.observe(event['step'], native_epoch_mean(batches, np.stack)))
        if history[-1]['early_stop']:
            reason = 'validation_early_stop'
            break
    write_json(output, {'monitor': MONITOR, 'best_step': policy.best_step,
                        'stop_reason': reason, 'history': history})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime-root', type=Path, default=REPO,
                        help='checkout with the isolated environment and pinned vendor tree')
    commands = parser.add_subparsers(dest='command', required=True)
    trace = commands.add_parser('replay-monitor')
    trace.add_argument('trace', type=Path)
    trace.add_argument('output', type=Path)
    for name in ('prepare-tenv', 'prepare-tahoe'):
        prepare = commands.add_parser(name)
        prepare.add_argument('input_root', type=Path)
        prepare.add_argument('request', type=Path)
        prepare.add_argument('--fold', type=int, required=True)
    for name in ('build-tahoe', 'verify-tahoe'):
        commands.add_parser(name).add_argument('input_root', type=Path)
    verify = commands.add_parser('verify-request')
    verify.add_argument('request', type=Path)
    fit = commands.add_parser('fit')
    fit.add_argument('request', type=Path)
    fit.add_argument('recipe', type=Path)
    fit.add_argument('output', type=Path)
    reload = commands.add_parser('reload')
    reload.add_argument('request', type=Path)
    reload.add_argument('run', type=Path)
    reload.add_argument('output', type=Path)
    args = parser.parse_args()
    for name in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ.setdefault(name, '2')
    sys.dont_write_bytecode = True
    runtime = args.runtime_root.resolve()
    sys.path.insert(0, str(REPO / STEP03 / 'execution/script'))
    sys.path.insert(0, str(runtime / STEP03 / 'execution/script'))
    if args.command == 'replay-monitor':
        replay_monitor(args.trace, args.output)
    elif args.command == 'prepare-tenv':
        from training_inputs import prepare_tenv
        prepare_tenv(runtime, args.input_root, args.request, args.fold)
    elif args.command in ('build-tahoe', 'verify-tahoe', 'prepare-tahoe'):
        from tahoe_inputs import build_tahoe, verify_tahoe, prepare_tahoe
        if args.command == 'build-tahoe':
            build_tahoe(runtime, args.input_root)
        elif args.command == 'verify-tahoe':
            verify_tahoe(runtime, args.input_root)
        else:
            prepare_tahoe(args.input_root, args.request, args.fold)
    else:
        from artifacts import read_json
        from training_inputs import verify_request
        request = read_json(args.request)
        verify_request(request)
        if args.command == 'fit':
            from native_training import fit
            fit(runtime, request, args.recipe, args.output)
        elif args.command == 'reload':
            from native_training import reload
            reload(request, args.run, args.output)


if __name__ == '__main__':
    main()
