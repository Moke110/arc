"""Verify the complete formal matrix and prepare its Stage 2 Step 05 report."""
import argparse
import csv
from math import isfinite
from pathlib import Path
from statistics import mean

PATHS = ('tenv-only', 'tahoe-only', 'transfer')
DOSES = (0.2, 1.0, 10.0)
DIFFERENCES = ('transfer_minus_tenv_only', 'transfer_minus_tahoe_only')


def comparison(receipts):
    """Build per-fold and equal-fold-mean report rows for the nine paired evaluations."""
    indexed = {(receipt['fold'], receipt['path']): receipt for receipt in receipts}
    if len(receipts) != 9 or indexed.keys() != {
        (fold, path) for fold in range(3) for path in PATHS
    }:
        raise ValueError('report requires exactly nine fold-by-path test evaluations')
    scores = []
    differences = []
    for fold, dose in enumerate(DOSES):
        selected = [indexed[fold, path] for path in PATHS]
        if len({receipt['files']['identity.csv'] for receipt in selected}) != 1:
            raise ValueError('report paths do not share ordered test identities')
        current = []
        for path, receipt in zip(PATHS, selected, strict=True):
            metrics = receipt['metrics']
            score, count = metrics['pcc'], metrics['n_samples']
            undefined = metrics['undefined_pcc_count']
            if score is None or not isfinite(score) or not -1 <= score <= 1 or not 0 <= undefined < count:
                raise ValueError('report requires finite Delta PCC for every test path and fold')
            current.append(dict(fold=fold, held_dose_um=dose, path=path, pcc=score,
                                n_samples=count, undefined_pcc_count=undefined))
        if len({row['n_samples'] for row in current}) != 1:
            raise ValueError('report path test counts differ')
        tenv, tahoe, transfer = (row['pcc'] for row in current)
        scores.extend(current)
        differences.append(dict(fold=fold, held_dose_um=dose,
                                transfer_minus_tenv_only=transfer - tenv,
                                transfer_minus_tahoe_only=transfer - tahoe))
    for path in PATHS:
        selected = [row for row in scores if row['path'] == path]
        scores.append(dict(fold='mean', held_dose_um='', path=path,
                           pcc=mean(row['pcc'] for row in selected),
                           n_samples=sum(row['n_samples'] for row in selected),
                           undefined_pcc_count=sum(row['undefined_pcc_count'] for row in selected)))
    differences.append(dict(fold='mean', held_dose_um='',
                            **{key: mean(row[key] for row in differences) for key in DIFFERENCES}))
    return scores, differences


def verify_matrix(inputs, native_source, candidate):
    """Verify source selection, both derived paths and all nine fit/test receipts."""
    import sys
    from types import SimpleNamespace

    sys.path.insert(0, str(STEP / 'output/api'))
    from tsg_transfer import fit, tahoe, tenv
    from tenv5k.model_transigen import SCRIPT_HASH, SOURCE_HASH
    from transfer import run_with_inputs as verify_derived

    if (fit.sha256(native_source / 'src/model.py') != SOURCE_HASH
            or fit.sha256(native_source / 'src/train_TranSiGen.py') != SCRIPT_HASH):
        raise ValueError('native TranSiGen source differs from the pinned checkout')
    source = fit.load_inputs(inputs / 'tahoe_bulk', 'tahoe_bulk')
    target = fit.load_inputs(inputs / 'tenv', 'tenv')
    tahoe_candidate = candidate / 'tahoe'
    split = candidate / 'tenv-split'
    winner_path = tahoe_candidate / 'winner.json'
    fit.read_json(winner_path)
    winner = tahoe.lock_winner(source, tahoe_candidate)
    test_receipts = []
    fits = []
    tests = []

    for fold in tenv.FOLDS:
        rows = tenv.load_split(target, split, fold)
        control = candidate / 'tenv-only' / f'fold-{fold}'
        control_fit = control / 'fit'
        fit.verify_replay(target, control_fit, train_rows=rows['train'],
                          validation_rows=rows['validation'], kind='tenv-only', fold=fold)
        test_receipts.append(tenv.verify_test(
            target, split, fold, control / 'test', path='tenv-only',
            lineage=tenv.control_lineage(control_fit),
        ))
        for path in ('tahoe-only', 'transfer'):
            test_receipts.append(verify_derived(SimpleNamespace(
                action='verify', path=path, fold=fold, smoke=False, native_source=native_source,
                tahoe_inputs=inputs / 'tahoe_bulk', tenv_inputs=inputs / 'tenv',
                tahoe_candidate=tahoe_candidate, tenv_candidate=candidate, candidate=candidate,
            ), source, target))
        for kind, directory in (
            ('tahoe', tahoe_candidate / f'fold-{fold}'),
            ('tenv-only', control_fit),
            ('transfer', candidate / 'transfer' / f'fold-{fold}' / 'fit'),
        ):
            receipt = fit.read_json(directory / 'fit.json')
            fits.append(dict(
                kind=kind, fold=fold, fit_sha256=fit.sha256(directory / 'fit.json'),
                replay_sha256=fit.sha256(directory / 'validation-replay.json'),
                artifact_sha256=receipt['files'],
                **{key: receipt[key] for key in (
                    'epochs', 'best_epoch', 'training_steps', 'validation_mse', 'stop_reason',
                    'transform', 'elapsed_seconds', 'peak_rss_kib', 'environment', 'initialization',
                )},
            ))
        for path in PATHS:
            directory = candidate / path / f'fold-{fold}' / 'test'
            receipt = fit.read_json(directory / 'manifest.json')
            tests.append(dict(
                path=path, fold=fold, test_sha256=fit.sha256(directory / 'manifest.json'),
                replay_sha256=fit.sha256(directory / 'replay.json'),
                artifact_sha256=receipt['files'], lineage=receipt['lineage'], metrics=receipt['metrics'],
            ))

    evidence = dict(
        pair_sha256=fit.sha256(inputs / 'manifest.json'),
        tenv_split=fit.read_json(split / 'manifest.json'),
        tenv_split_sha256=fit.sha256(split / 'manifest.json'),
        tahoe_split=fit.read_json(tahoe_candidate / 'splits.json'),
        tahoe_winner=winner,
        tahoe_winner_sha256=fit.sha256(winner_path),
        fits=fits, tests=tests,
    )
    return comparison(test_receipts), evidence


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path, scores, differences, evidence):
    score_table = dict(((row['fold'], row['path']), row) for row in scores)
    winner = evidence['tahoe_winner']['winner']
    lines = [
        '# Accepted Tahoe-to-tenv comparison',
        '',
        f"Tahoe validation selected fold {winner['fold']} at epoch {winner['best_epoch']} "
        f"(treated-expression MSE {winner['validation_mse']:.8f}).",
        '',
        '| tenv fold | Held dose (µM) | Test profiles | tenv-only | Tahoe-only | Transfer | '
        'Transfer − tenv-only | Transfer − Tahoe-only |',
        '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |',
    ]
    for row in differences:
        fold = row['fold']
        label = 'Equal-fold mean' if fold == 'mean' else str(fold)
        dose = row['held_dose_um'] if fold != 'mean' else '—'
        selected = [score_table[fold, model_path] for model_path in PATHS]
        line = [label, str(dose), str(selected[0]['n_samples'])]
        line.extend(f"{score['pcc']:.6f}" for score in selected)
        line.extend(f'{row[key]:+.6f}' for key in DIFFERENCES)
        lines.append('| ' + ' | '.join(line) + ' |')
    lines += [
        '',
        'Scores are sample-equal Delta PCC within each fold. The final row weights folds equally. '
        'Its profile count is the total across the three folds. '
        'The [plot table](../figure/data/transfer_scores.csv) records undefined PCC counts by path and fold.',
        '',
        'Transfer − tenv-only measures the effect of Tahoe initialization under the fixed tenv training protocol. '
        'Transfer − Tahoe-only measures the effect of tenv adaptation from the locked Tahoe checkpoint.',
        '',
        'The [acceptance receipt](acceptance.json) binds all nine fits, nine test evaluations and fresh-process replays.',
        '',
    ]
    path.write_text('\n'.join(lines))


STEP = Path(__file__).resolve().parents[1]


def build_report(inputs, native_source, candidate, phases, report):
    """Create the complete report in an isolated, empty destination."""
    import subprocess
    import sys

    from report_resources import verify_resources

    (scores, differences), evidence = verify_matrix(inputs, native_source, candidate)
    resources = verify_resources(phases, evidence, candidate)
    data = report / 'figure/data'
    data.mkdir(parents=True)
    log = report / 'log'
    log.mkdir()
    write_csv(data / 'transfer_scores.csv', scores)
    write_csv(data / 'transfer_differences.csv', differences)
    subprocess.run([sys.executable, str(STEP / 'figure/script/render.py'),
                    '--data', str(data), '--output', str(report / 'figure/transfer.png')], check=True)
    from tsg_transfer import fit

    write_summary(log / 'summary.md', scores, differences, evidence)
    fit.write_json(log / 'acceptance.json', dict(
        evidence, schema='stage2-tsg-transfer-acceptance-v1', resources=resources,
        report_sha256={str(path.relative_to(report)): fit.sha256(path)
                       for path in sorted(report.rglob('*')) if path.is_file()},
    ))



def verify_report_files(staged, report):
    """Require the exact generated report tree and contents in its archive."""
    expected = {path.relative_to(staged) for path in staged.rglob('*')}
    actual = {path.relative_to(report) for path in report.rglob('*')}
    if actual != expected:
        raise ValueError('published report inventory differs')
    for relative in expected:
        source, destination = staged / relative, report / relative
        if (destination.is_symlink() or destination.is_file() != source.is_file()
                or (source.is_file() and source.read_bytes() != destination.read_bytes())):
            raise ValueError(f'published report differs: {relative}')

def main():
    import tempfile

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'verify'))
    for name in ('inputs', 'native-source', 'candidate', 'phases', 'report'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    if args.action == 'prepare':
        if args.report.exists():
            parser.error('prepare requires a fresh report destination')
        args.report.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=args.report.parent) as temporary:
        staged = Path(temporary) / 'report'
        build_report(args.inputs, args.native_source, args.candidate, args.phases, staged)
        if args.action == 'prepare':
            staged.rename(args.report)
        else:
            verify_report_files(staged, args.report)
    print(f'{args.action}: verified complete formal Stage 2 Step 05 report')


if __name__ == '__main__':
    main()
