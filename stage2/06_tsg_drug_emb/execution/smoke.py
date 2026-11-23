"""Run every representation through bounded native fit, validation and test replay."""
import argparse
from pathlib import Path
import subprocess
import sys

from fit import ARMS, WIDTHS, sha256, verify_phase, write_json


def run(args):
    args.destination.mkdir(parents=True)
    fits = []
    command = [sys.executable, str(Path(__file__).with_name('fit.py')),
               '--inputs', str(args.inputs), '--representations', str(args.representations),
               '--source', str(args.source), '--fold', str(args.fold), '--device', args.device,
               '--smoke-rows', '12', '--smoke-epochs', '2']
    for arm in ARMS:
        candidate = args.destination / arm
        with (args.destination / f'{arm}.log').open('w') as log:
            for action in ('train', 'replay-validation', 'evaluate', 'replay-test'):
                subprocess.run(command + ['--arm', arm, '--run', str(candidate), action],
                               stdout=log, stderr=subprocess.STDOUT, check=True)
        phases = {phase: verify_phase(candidate / phase, phase)
                  for phase in ('fit', 'validation-replay', 'test', 'test-replay')}
        fit = phases['fit']
        config = fit['config']
        if (config['width'] != WIDTHS[arm] or config['arm'] != arm
                or config['seed'] != 42 or config['fold'] != args.fold
                or fit['epochs'] != 2 or fit['optimizer_steps'] != 2
                or len({receipt['process_token'] for receipt in phases.values()}) != 4
                or {len(rows) for rows in config['order'].values()} != {12}):
            raise ValueError(f'invalid bounded origin or replay for {arm}')
        tensor = fit['initialization']['tensors']['feat_embeddings.0.weight']
        if tensor['shape'] != [400, WIDTHS[arm]]:
            raise ValueError(f'first affine input/output shape differs for {arm}')
        fits.append(dict(arm=arm, width=config['width'], history_epochs=fit['epochs'],
                         optimizer_steps=fit['optimizer_steps'], best_epoch=fit['best_epoch'],
                         validation_mse=fit['validation_mse'],
                         downstream_sha256=fit['initialization']['downstream_sha256'],
                         rng_sha256=fit['initialization']['rng_sha256'],
                         split_order=config['order'], transform=config['architecture']['transform'],
                         phase_sha256={phase: sha256(candidate / phase / 'manifest.json')
                                       for phase in phases}))
        print(f'{arm}: bounded independent fit and exact validation/test replay', flush=True)
    for field in ('downstream_sha256', 'rng_sha256', 'split_order', 'transform'):
        if any(fit[field] != fits[0][field] for fit in fits):
            raise ValueError(f'representation arms differ outside the first drug layer: {field}')
    receipt = dict(schema='tsg-drug-emb-bounded-smoke-v1', formal_fits=0,
                   fold=args.fold, device=args.device, train_seed=42,
                   representations_sha256=sha256(args.representations / 'manifest.json'),
                   inputs_sha256=sha256(args.inputs / 'manifest.json'),
                   fits=fits, status='PASS')
    write_json(args.destination / 'acceptance.json', receipt)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ('inputs', 'representations', 'source', 'destination'):
        parser.add_argument('--' + argument, type=Path, required=True)
    parser.add_argument('--device', choices=('cpu', 'cuda'), required=True)
    parser.add_argument('--fold', type=int, choices=range(3), default=0)
    run(parser.parse_args())


if __name__ == '__main__':
    main()
