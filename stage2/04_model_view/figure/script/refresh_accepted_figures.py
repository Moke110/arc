"""Verify a frozen accepted matrix, redraw its CSVs, and transactionally republish.

Only figure hashes change. Training, inference, and scientific scoring never run.
Interrupted candidate changes are restored before another refresh attempt.
"""
import argparse
import copy
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

STEP = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('frozen_matrix', STEP / 'execution/run_matrix.py')
matrix = importlib.util.module_from_spec(spec)
spec.loader.exec_module(matrix)


def restore(candidate, journal):
    state = json.loads(journal.read_text())
    if not state.get('committed'):
        backup = Path(state['backup'])
        shutil.rmtree(candidate / 'figures', ignore_errors=True)
        shutil.copytree(backup / 'figures', candidate / 'figures')
        shutil.copy2(backup / 'acceptance.json', candidate / 'acceptance.json')
    journal.unlink()


def refresh(args):
    candidate = args.candidate
    journal = candidate / 'figure-refresh-journal.json'
    if journal.exists():
        restore(candidate, journal)
    acceptance_path = candidate / 'acceptance.json'
    original_sha = matrix.digest(acceptance_path)
    accepted = json.loads(acceptance_path.read_text())
    frozen = json.loads((candidate / 'plan.json').read_text())
    if (accepted['schema'] != matrix.ACCEPTANCE_SCHEMA or accepted['runs'] != 54
            or frozen != matrix.frozen_plan(args)
            or accepted['plan_sha256'] != matrix.digest(candidate / 'plan.json')):
        raise ValueError('accepted matrix identity differs')
    for name, expected in accepted['files'].items():
        if Path(name).is_absolute() or '..' in Path(name).parts or matrix.digest(candidate / name) != expected:
            raise ValueError('accepted artifact changed: ' + name)
    for run in matrix.plan():
        matrix.verify_execution(matrix.run_path(candidate, run), run, frozen)
    science = {name: value for name, value in accepted['files'].items() if name.startswith('comparison/')}
    renderer = STEP / 'figure/script/render_figures.py'
    with tempfile.TemporaryDirectory(prefix='figure-refresh-', dir=candidate) as temporary:
        temporary = Path(temporary)
        for name in ('first', 'replay'):
            subprocess.run([args.python, str(renderer), '--figure-data', str(candidate / 'comparison'),
                            '--out-dir', str(temporary / name)], check=True)
        first = {p.name: matrix.digest(p) for p in (temporary / 'first').iterdir()}
        replay = {p.name: matrix.digest(p) for p in (temporary / 'replay').iterdir()}
        if first != replay or set(first) != {'learning_curve.png', 'metrics.png', 'paired_differences.png'}:
            raise ValueError('CSV-only exact figure replay differs')
        if any(matrix.digest(candidate / name) != value for name, value in science.items()):
            raise ValueError('scientific report changed during rendering')
        updated = copy.deepcopy(accepted)
        old_png = {name: value for name, value in accepted['files'].items() if name.startswith('figures/')}
        new_png = {'figures/' + name: value for name, value in first.items()}
        updated['files'].update(new_png)
        provenance = dict(reason='Separate native training loss and validation MSE onto labeled independent y axes',
                          previous_acceptance_sha256=original_sha, renderer_sha256=matrix.digest(renderer),
                          previous_figures_sha256=old_png, figures_sha256=new_png,
                          unchanged_scientific_files_sha256=science,
                          verified_fresh_execution_receipts=54, csv_only_exact_replay=True,
                          frozen_plan_sha256=accepted['plan_sha256'])
        updated['figure_refresh'] = provenance
        backup = candidate / 'figure-refresh-history' / original_sha
        if backup.exists():
            if matrix.digest(backup / 'acceptance.json') != original_sha:
                raise ValueError('figure refresh backup identity differs')
        else:
            backup.mkdir(parents=True)
            shutil.copytree(candidate / 'figures', backup / 'figures')
            shutil.copy2(acceptance_path, backup / 'acceptance.json')
        pending = temporary / 'acceptance.json'
        pending.write_text(json.dumps(updated, indent=2) + '\n')
        journal.write_text(json.dumps(dict(backup=str(backup), committed=False)) + '\n')
        try:
            shutil.rmtree(candidate / 'figures')
            os.rename(temporary / 'first', candidate / 'figures')
            os.replace(pending, acceptance_path)
            # Original immutable-runtime verifier and transactional installer.
            matrix.publish(args)
        except BaseException:
            restore(candidate, journal)
            raise
        journal.write_text(json.dumps(dict(backup=str(backup), committed=True)) + '\n')
        journal.unlink()
    (STEP / 'log/231-figure-refresh.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(json.dumps(provenance))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('candidate', 'inputs', 'means', 'prnet-source', 'transigen-source'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--device', default='cuda')
    parser.add_argument('--python', required=True)
    refresh(parser.parse_args())
