"""Dispose of selected accepted Step 04 weights or verify the prediction archive."""
import argparse
import json
from pathlib import Path

from report_formal import STAGE, STEP, report, run_paths
from cpst_util.checkpoint_disposition import WEIGHT, dispose_checkpoints, verify_prediction_archive
from cpst_util.derivation.release import check_outputs, read_json, sha256

ACCEPTED = STEP / "log/399-tenv-v1-formal.json"
ACCEPTED_SHA256 = "e352942a1d45dd0759569259c21162419093a9419d760c374e3c3b7796fcda00"
FIT_SOURCE = "4ae0d64195f4d983b144326592608016e18057ad"
DISPOSITION = "checkpoint-disposition-v1.json"


def _accepted(args):
    if sha256(ACCEPTED) != ACCEPTED_SHA256:
        raise ValueError("Step 04 full acceptance receipt digest differs")
    accepted = read_json(ACCEPTED)
    if (accepted["result"] != "PASS" or accepted["fit_source_commit"] != FIT_SOURCE
            or accepted["matrix_manifest_sha256"] != sha256(args.matrix / "manifest.json")
            or accepted["report_manifest_sha256"] != sha256(args.report / "manifest.json")):
        raise ValueError("Step 04 archive is limited to the accepted installed rerun")
    return accepted


def _authorized_runs(matrix, accepted, names):
    accepted_fits = {(fit['axis'], fit['fold'], fit['seed']): fit for fit in accepted['matrix_fits']}
    paths = run_paths(matrix)
    if len(accepted['matrix_fits']) != len(paths) or set(accepted_fits) != set(paths):
        raise ValueError('Step 04 full acceptance does not bind the complete run roster')
    authorized = {(path / 'manifest.json').relative_to(matrix).as_posix(): accepted_fits[key]['run_manifest_sha256']
                  for key, path in paths.items()}
    if not names or len(names) != len(set(names)) or set(names) - authorized.keys():
        raise ValueError('name a unique subset of accepted run manifests relative to the matrix')
    return {name: authorized[name] for name in names}


def _retained_output_check(matrix, runs):
    def check(path, manifest):
        manifest_path = path / 'manifest.json'
        pin = runs.get(manifest_path.relative_to(matrix).as_posix())
        if pin is None:
            return check_outputs(path, manifest)
        if sha256(manifest_path) != pin or WEIGHT not in manifest['outputs']:
            raise ValueError(f'archived run manifest differs: {path}')
        weight = path / WEIGHT
        if weight.exists() or weight.is_symlink():
            raise ValueError(f'archived checkpoint remains at the run path: {weight}')
        retained = {name: digest for name, digest in manifest['outputs'].items() if name != WEIGHT}
        return check_outputs(path, {'outputs': retained})

    return check


def _verified_report(args, accepted, *, output_check=None):
    verified = report(args, archive_check=output_check)
    for key in ('matrix_manifest_sha256', 'report_manifest_sha256', 'csv_only_figure_sha256'):
        if verified[key] != accepted[key]:
            raise ValueError(f'prediction archive differs from the full acceptance: {key}')
    return verified


def _dispose_roster(args, accepted, receipt, pin):
    leaf = Path(receipt)
    if not leaf.name or leaf.name != receipt or leaf.name == '..':
        raise ValueError('disposition receipt must be a file at the matrix root')
    path = args.matrix / leaf
    if path.is_symlink() or sha256(path) != pin:
        raise ValueError('disposition receipt digest differs')
    record = read_json(path)
    runs = _authorized_runs(args.matrix, accepted, [item['manifest'] for item in record['weights']])
    verify_prediction_archive(args.matrix, receipt=receipt, receipt_sha256=pin, acceptance=ACCEPTED,
                              acceptance_sha256=ACCEPTED_SHA256, runs=runs)
    return runs


def verify(args):
    accepted = _accepted(args)
    receipt = getattr(args, 'receipt', DISPOSITION)
    runs = _dispose_roster(args, accepted, receipt, args.disposition_sha256)
    verified = _verified_report(args, accepted, output_check=_retained_output_check(args.matrix, runs))
    return {**verified, 'verification': 'predictions-and-figures',
            'disposition_sha256': args.disposition_sha256, 'disposed_weights': len(runs)}


def dispose(args):
    accepted = _accepted(args)
    runs = _authorized_runs(args.matrix, accepted, args.run_manifest)
    receipt = getattr(args, 'receipt', DISPOSITION)
    if Path(receipt).name != receipt or not receipt or receipt == '..':
        raise ValueError('disposition receipt must be a file at the matrix root')
    previous_receipt = getattr(args, 'previous_receipt', None)
    previous_pin = getattr(args, 'previous_disposition_sha256', None)
    if bool(previous_receipt) != bool(previous_pin):
        raise ValueError('previous disposition requires both its receipt and digest')
    previous = None
    prior_runs = {}
    if previous_receipt:
        prior_runs = _dispose_roster(args, accepted, previous_receipt, previous_pin)
        previous = previous_receipt, previous_pin, prior_runs

    def verify_before():
        checker = _retained_output_check(args.matrix, prior_runs) if prior_runs else None
        _verified_report(args, accepted, output_check=checker)

    def verify_retained(_entries):
        _verified_report(args, accepted, output_check=_retained_output_check(args.matrix, prior_runs | runs))

    pin = dispose_checkpoints(args.matrix, receipt=receipt, acceptance=ACCEPTED,
                              acceptance_sha256=ACCEPTED_SHA256, runs=runs, previous=previous,
                              verify_full=verify_before, verify_retained=verify_retained)
    return {'result': 'PASS', 'verification': 'predictions-and-figures',
            'disposition_sha256': pin, 'disposed_weights': len(runs)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='action', required=True)
    for action in (dispose, verify):
        command = commands.add_parser(action.__name__)
        command.add_argument('--matrix', type=Path, required=True)
        command.add_argument('--receipt', default=DISPOSITION)
        command.add_argument('--matrix-audit', type=Path, required=True)
        command.add_argument('--report', type=Path, required=True)
        command.add_argument('--figures', type=Path, required=True)
        command.add_argument('--logs', type=Path, required=True)
        command.add_argument('--step03', type=Path, default=STAGE / '03_data_train')
        command.add_argument('--raw', type=Path, default=STAGE / '01_data_raw/output')
        command.set_defaults(handler=action, command='verify')
    commands.choices['dispose'].add_argument('--run-manifest', action='append', required=True)
    commands.choices['dispose'].add_argument('--previous-receipt')
    commands.choices['dispose'].add_argument('--previous-disposition-sha256')
    commands.choices['verify'].add_argument('--disposition-sha256', required=True)
    args = parser.parse_args()
    args.matrix = args.matrix.resolve(strict=True)
    print(json.dumps(args.handler(args), sort_keys=True))


if __name__ == '__main__':
    main()
