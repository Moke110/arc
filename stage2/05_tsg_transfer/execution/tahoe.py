"""Prepare, fit, independently replay and lock the Tahoe source checkpoint."""
import argparse
import json
from pathlib import Path
import sys

STEP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STEP / 'output/api'))

from tsg_transfer.fit import load_inputs, replay_fit, train
from tsg_transfer.tahoe import lock_winner, prepare_tahoe, tahoe_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'fit', 'replay', 'select'))
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--native-source', type=Path)
    parser.add_argument('--fold', type=int, choices=range(3))
    parser.add_argument('--device', default='cuda')
    parser.add_argument('--smoke-epochs', type=int)
    parser.add_argument('--smoke-rows', type=int)
    args = parser.parse_args()
    if args.action in ('fit', 'replay') and (args.native_source is None or args.fold is None):
        parser.error('fit and replay require --native-source and --fold')
    if args.action == 'fit' and (args.smoke_epochs is None) != (args.smoke_rows is None):
        parser.error('a bounded fit requires both --smoke-epochs and --smoke-rows')
    if args.action == 'fit' and args.smoke_rows is None and not args.device.startswith('cuda'):
        parser.error('formal fit requires CUDA')
    data = load_inputs(args.inputs, 'tahoe_bulk')
    if args.action == 'select':
        report = lock_winner(data, args.candidate, smoke_rows=args.smoke_rows)
    else:
        splits = prepare_tahoe(data, args.candidate)
        if args.action == 'prepare':
            report = splits.groupby(['fold', 'partition']).size().unstack().to_dict()
        else:
            train_rows, validation_rows = tahoe_rows(data, splits, args.fold, smoke_rows=args.smoke_rows)
            attempt = args.candidate / f'fold-{args.fold}'
            partitions = dict(train_rows=train_rows, validation_rows=validation_rows)
            if args.action == 'fit':
                fit = train(data, args.native_source, attempt, **partitions, kind='tahoe',
                            fold=args.fold, device=args.device, smoke_epochs=args.smoke_epochs)
                report = {key: fit[key] for key in ('epochs', 'best_epoch', 'training_steps',
                                                   'validation_mse', 'elapsed_seconds', 'peak_rss_kib')}
            else:
                report = replay_fit(data, args.native_source, attempt, **partitions)
    print(json.dumps(report, sort_keys=True))


if __name__ == '__main__':
    main()
