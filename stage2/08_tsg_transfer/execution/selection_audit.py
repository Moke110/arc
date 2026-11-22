"""Audit preserved validation evidence and compare bounded development selection rules."""
import argparse
from pathlib import Path
import sys

import pandas as pd

STEP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STEP / 'output/api'))

from tsg_transfer import fit, tenv, transfer
import development


def audit(data, experiment):
    """Read only frozen splits and MSE-selected fit evidence."""
    report = []
    for fold in tenv.FOLDS:
        roles = tenv.load_split(data, experiment / 'tenv-split', fold)
        for arm in ('tenv-only', 'transfer'):
            path = experiment / arm / f'fold-{fold}' / 'fit'
            receipt = fit.verify_replay(data, path, train_rows=roles['train'],
                                        validation_rows=roles['validation'], kind=arm, fold=fold)
            curve = pd.read_csv(path / 'learning_curve.csv', float_precision='round_trip')
            selected = curve.iloc[receipt['best_epoch'] - 1]
            report.append(dict(fold=fold, arm=arm, seed=receipt['seed'],
                               n_validation_profiles=len(roles['validation']),
                               n_validation_parents=int(
                                   data.obs.iloc[roles['validation']].drug_parent_inchikey.nunique()),
                               first_mse=float(curve.validation_mse.iloc[0]),
                               best_mse=receipt['validation_mse'], best_epoch=receipt['best_epoch'],
                               updates_at_best=int(selected.training_steps), epochs=receipt['epochs'],
                               recoverable_epochs=[receipt['best_epoch']], delta_pcc_trace=False,
                               fit_sha256=fit.sha256(path / 'fit.json')))
    return report


def run(args):
    data = fit.load_inputs(args.inputs, 'tenv')
    audit_report = audit(data, args.experiment)
    if not args.epochs:
        return dict(audit=audit_report)

    roles = tenv.load_split(data, args.experiment / 'tenv-split', args.fold)
    split = development.split_validation(data, roles['train'], roles['validation'])
    source_fit, winner, source_doses = development.load_tahoe_winner(
        args.tahoe_inputs, args.experiment / 'tahoe', args.source,
    )

    def initialize(recipient, transform):
        return transfer.initialize_from_tahoe(source_fit, recipient, winner['transform'], transform,
                                             source_doses + tenv.DOSES, held_dose=tenv.DOSES[args.fold])

    results = []
    for seed in args.seeds:
        for lr in args.learning_rates:
            for arm, init in (('tenv-only', None), ('transfer', initialize)):
                print(f'seed={seed} learning_rate={lr} arm={arm}', file=sys.stderr, flush=True)
                result = development.run_fit(data, args.source, split, learning_rate=lr,
                                             seed=seed, epochs=args.epochs, initialize=init)
                selected = {name: {key: value[key] for key in ('selected', 'assessment')}
                            for name, value in result['selected'].items()}
                results.append(dict(seed=seed, learning_rate=lr, arm=arm, selected=selected,
                                    history=result['history']))
    return dict(audit=audit_report, development=dict(
        fold=args.fold, held_dose_um=tenv.DOSES[args.fold], epochs=args.epochs,
        split_sha256=fit.sha256(args.experiment / 'tenv-split' / 'manifest.json'),
        winner_sha256=winner['checkpoint_sha256'], learning_rates=args.learning_rates,
        seeds=args.seeds, identities=development.role_identities(data, split), results=results,
    ))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('inputs', 'experiment', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    for name in ('tahoe-inputs', 'source'):
        parser.add_argument('--' + name, type=Path)
    parser.add_argument('--fold', type=int, choices=tenv.FOLDS, default=0)
    parser.add_argument('--epochs', type=int, default=0)
    parser.add_argument('--seeds', nargs='+', type=int, default=[42, 17])
    parser.add_argument('--learning-rates', nargs='+', type=float, default=[1e-5, 1e-4])
    args = parser.parse_args()
    if args.epochs and (args.tahoe_inputs is None or args.source is None):
        parser.error('--epochs requires --tahoe-inputs and --source')
    if args.output.exists():
        parser.error('--output must be new')
    fit.write_json(args.output, run(args))


if __name__ == '__main__':
    main()
