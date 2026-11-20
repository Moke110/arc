"""Prepare, train and independently replay Step 04 tenv-only controls."""
import argparse
import json
from pathlib import Path
import sys

STEP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STEP / 'output/api'))

from tsg_transfer import fit, tenv
from tsg_round1.learner import smoke_subset


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'fit', 'replay', 'evaluate', 'replay-test'))
    parser.add_argument('--inputs', type=Path, required=True, help='Accepted Step 03 tenv input directory')
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--source', type=Path, help='Pinned native TranSiGen source checkout')
    parser.add_argument('--fold', type=int, choices=tenv.FOLDS)
    parser.add_argument('--device', default='cuda')
    parser.add_argument('--smoke-epochs', type=int)
    args = parser.parse_args()
    if args.action != 'prepare' and (args.source is None or args.fold is None):
        parser.error(f'{args.action} requires --source and --fold')
    if args.smoke_epochs is not None and args.action != 'fit':
        parser.error('--smoke-epochs is only accepted for fit')

    data = fit.load_inputs(args.inputs, 'tenv')
    split = args.candidate / 'tenv-split'
    if args.action == 'prepare':
        print(tenv.freeze_split(data, split))
        return

    rows = tenv.load_split(data, split, args.fold)
    run = args.candidate / 'tenv-only' / f'fold-{args.fold}'
    fit_dir = run / 'fit'
    smoke_epochs = (args.smoke_epochs if args.action == 'fit'
                    else fit.read_json(fit_dir / 'fit.json')['smoke_epochs'])
    train_rows = rows['train']
    validation_rows = rows['validation']
    if smoke_epochs is not None:
        train_rows = smoke_subset(data, train_rows, 128)
        validation_rows = smoke_subset(data, validation_rows, 64)
    partition = dict(train_rows=train_rows, validation_rows=validation_rows)
    if args.action == 'fit':
        result = fit.train(data, args.source, fit_dir, **partition, fold=args.fold,
                           kind='tenv-only', device=args.device, smoke_epochs=args.smoke_epochs)
    elif args.action == 'replay':
        result = fit.replay_fit(data, args.source, fit_dir, **partition)
    else:
        native, model, receipt = fit.load_fit(data, args.source, fit_dir, **partition)
        fit.verify_replay(data, fit_dir, **partition, kind='tenv-only', fold=args.fold,
                          formal=receipt['smoke_epochs'] is None)
        test = rows['test']
        prediction = fit.predict(native, model, data, test)
        run_test = tenv.export_test if args.action == 'evaluate' else tenv.replay_test
        result = run_test(data, split, args.fold, run / 'test', prediction, path='tenv-only',
                          lineage=tenv.control_lineage(fit_dir), prediction_ids=data.obs.index[test])
    print(json.dumps(result, default=str))


if __name__ == '__main__':
    main()
