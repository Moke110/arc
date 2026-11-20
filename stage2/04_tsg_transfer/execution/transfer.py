"""Check the dose handoff, fit transfer and replay both Tahoe-derived paths."""
import argparse
import json
from pathlib import Path
import sys

STEP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STEP / 'output/api'))

from tsg_transfer import fit, tahoe, tenv, transfer
from tsg_round1 import model as models
from tsg_round1.learner import smoke_subset
from tsg_round1.protocol import TRAINING_SEED


def run(args):
    source = fit.load_inputs(args.tahoe_inputs, 'tahoe_bulk')
    target = fit.load_inputs(args.tenv_inputs, 'tenv')
    if not source.genes.equals(target.genes):
        raise ValueError('Tahoe and tenv ordered G_5000 axes differ')

    winner_path = args.tahoe_candidate / 'winner.json'
    sealed = fit.read_json(winner_path)
    smoke_rows = sealed['smoke_rows']
    if (smoke_rows is not None) != args.smoke:
        raise ValueError('locked Tahoe winner does not match the requested smoke or formal run')
    lock = tahoe.lock_winner(source, args.tahoe_candidate, smoke_rows=smoke_rows)
    winner = lock['winner']
    source_fit = args.tahoe_candidate / f"fold-{winner['fold']}"
    source_origin = dict(
        winner_fold=winner['fold'], winner_sha256=fit.sha256(winner_path),
        winner_fit_sha256=winner['fit_sha256'], checkpoint_sha256=winner['checkpoint_sha256'],
    )

    split = args.tenv_candidate / 'tenv-split'
    rows = tenv.load_split(target, split, args.fold)
    control = args.tenv_candidate / 'tenv-only' / f'fold-{args.fold}'
    control_fit = control / 'fit'
    control_receipt = fit.read_json(control_fit / 'fit.json')
    smoke_epochs = control_receipt['smoke_epochs']
    if (smoke_epochs is not None) != args.smoke:
        raise ValueError('tenv-only control does not match the requested smoke or formal run')
    train_rows = rows['train']
    validation_rows = rows['validation']
    if args.smoke:
        train_rows = smoke_subset(target, train_rows, 128)
        validation_rows = smoke_subset(target, validation_rows, 64)
    partitions = dict(train_rows=train_rows, validation_rows=validation_rows)
    fit.verify_replay(target, control_fit, **partitions, kind='tenv-only', fold=args.fold,
                      formal=not args.smoke)
    target_origin = dict(tenv.control_lineage(control_fit),
                         split_manifest_sha256=fit.sha256(split / 'manifest.json'))

    origin = dict(source_origin, **target_origin, handoff_code_sha256=fit.sha256(transfer.__file__))
    transfer_fit = args.candidate / 'transfer' / f'fold-{args.fold}' / 'fit'
    needs_source_model = args.action in ('check', 'fit') or (
        args.path == 'tahoe-only' and args.action != 'verify'
    )
    if needs_source_model:
        source_rows = tahoe.tahoe_rows(
            source, tahoe.prepare_tahoe(source, args.tahoe_candidate),
            winner['fold'], smoke_rows=smoke_rows,
        )
        native, model, _ = fit.load_fit(source, args.native_source, source_fit,
                                        train_rows=source_rows[0], validation_rows=source_rows[1])

    if args.action in ('check', 'fit'):
        if args.device != control_receipt['device']:
            raise ValueError('transfer and tenv-only control require matching devices')
        physical_doses = tuple(source.obs.dose_um.unique()) + tenv.DOSES

        def initialize(recipient, standardizer):
            if standardizer != control_receipt['transform']:
                raise ValueError('tenv transfer dose transform differs from the tenv-only control')
            return transfer.initialize_from_tahoe(
                model, recipient, winner['transform'], standardizer, physical_doses,
                held_dose=tenv.DOSES[args.fold],
            )

        if args.action == 'check':
            standardizer = control_receipt['transform']
            recipient = models.adapter(args.native_source, fit.METHOD, standardizer).build(
                TRAINING_SEED, args.device,
            )
            return dict(origin=origin, bridge=initialize(recipient, standardizer))
        receipt = fit.train(target, args.native_source, transfer_fit, **partitions,
                            fold=args.fold, kind='transfer', device=args.device,
                            smoke_epochs=smoke_epochs, initialize=initialize, origin=origin)
        return {field: receipt[field] for field in (
            'initialization', 'epochs', 'best_epoch', 'training_steps',
            'validation_mse', 'elapsed_seconds', 'peak_rss_kib',
        )}

    if args.path == 'transfer':
        receipt = fit.read_json(transfer_fit / 'fit.json')
        if (receipt['origin'] != origin or receipt['transform'] != control_receipt['transform']
                or receipt['device'] != control_receipt['device']
                or receipt['smoke_epochs'] != smoke_epochs or receipt['kind'] != 'transfer'
                or receipt['fold'] != args.fold):
            raise ValueError('transfer checkpoint does not match the locked winner and control')
        if args.action == 'replay':
            return fit.replay_fit(target, args.native_source, transfer_fit, **partitions)
        fit.verify_replay(target, transfer_fit, **partitions, kind='transfer', fold=args.fold,
                          formal=not args.smoke)
        transfer.verify_bridge(receipt['initialization'], winner['transform'],
                               control_receipt['transform'],
                               tuple(source.obs.dose_um.unique()) + tenv.DOSES,
                               held_dose=tenv.DOSES[args.fold])
        if args.action != 'verify':
            native, model, _ = fit.load_fit(target, args.native_source, transfer_fit, **partitions)
        lineage = dict(origin, transfer_fit_sha256=fit.sha256(transfer_fit / 'fit.json'),
                       transfer_replay_sha256=fit.sha256(transfer_fit / 'validation-replay.json'))
    else:
        lineage = dict(source_origin, **target_origin)

    test_rows = rows['test']
    control_test = tenv.verify_test(target, split, args.fold, control / 'test',
                                    path='tenv-only', lineage=tenv.control_lineage(control_fit))
    output = args.candidate / args.path / f'fold-{args.fold}' / 'test'
    if args.action == 'verify':
        result = tenv.verify_test(target, split, args.fold, output, path=args.path, lineage=lineage)
    else:
        prediction = fit.predict(native, model, target, test_rows)
        operation = tenv.export_test if args.action == 'evaluate' else tenv.replay_test
        result = operation(target, split, args.fold, output, prediction, path=args.path,
                           lineage=lineage, prediction_ids=target.obs.index[test_rows])
    if fit.sha256(output / 'identity.csv') != control_test['files']['identity.csv']:
        raise ValueError('prediction identity differs from the tenv-only test path')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('check', 'fit', 'replay', 'evaluate', 'replay-test', 'verify'))
    for name in ('tahoe-inputs', 'tenv-inputs', 'tahoe-candidate', 'tenv-candidate',
                 'candidate', 'native-source'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--fold', type=int, choices=tenv.FOLDS, required=True)
    parser.add_argument('--path', choices=('transfer', 'tahoe-only'), default='transfer')
    parser.add_argument('--device', default='cuda')
    parser.add_argument('--smoke', action='store_true')
    args = parser.parse_args()
    if args.path == 'tahoe-only' and args.action not in ('evaluate', 'replay-test', 'verify'):
        parser.error('Tahoe-only direct prediction has no tenv fit or dose handoff')
    print(json.dumps(run(args), default=str, sort_keys=True))


if __name__ == '__main__':
    main()
