"""Plan, run and verify the tenv TranSiGen width-by-dose comparison."""

import argparse
import fcntl
import hashlib
from itertools import product
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
STEP04 = HERE.parents[1] / '07_tsg_dose'
ACCEPTED = STEP04 / 'log/274-kpgt-smoke/input-manifest.json'
WIDTHS = {'original': 2, '1.5x': 3, '2x': 4, '2.5x': 5}
METHODS = {
    'no_dose': 'no-dose',
    'latent_numeric': 'latent-numeric',
    'latent_embedding': 'latent-embedding',
    'FiLM': 'film',
}
SCENARIO = 'cold-drug-dose'
SEED = 42
SCHEMA = 'tsg-width-dose-fit-v1'
FILES = {'config.json', 'checkpoint.pt', 'learning_curve.csv', 'validation_prediction.npy',
         'test_per_sample.csv', 'test_metrics.json', 'terminal.json'}


def digest(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def load_json(path):
    return json.loads(path.read_text())


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def plan():
    return [dict(width=width, method=method, fold=fold, scenario=SCENARIO,
                 source='tenv', seed=SEED)
            for width, method, fold in product(WIDTHS, METHODS, range(3))]


def run_path(root, origin):
    return root / 'runs' / origin['width'] / origin['method'] / f"fold-{origin['fold']}"


def verify_run(path, origin, identity):
    expected = load_json(path / 'manifest.json')
    if {entry.name for entry in path.iterdir()} != FILES | {'manifest.json'}:
        raise ValueError(f'Unexpected file inventory in {path}')
    actual = dict(schema=SCHEMA, origin=origin, identity=identity,
                  files={name: digest(path / name) for name in FILES})
    if expected != actual:
        raise ValueError(f'Run identity or digest mismatch in {path}')


def verify(root):
    with (root / '.driver.lock').open('rb') as lock:
        fcntl.flock(lock, fcntl.LOCK_SH | fcntl.LOCK_NB)
        identity = load_json(root / 'plan.json')
        if identity['smoke_rows'] is not None or identity['max_epochs'] != 1000:
            raise ValueError('A formal report requires full runs')
        if identity['input_manifest_sha256'] != digest(ACCEPTED):
            raise ValueError('A formal report requires the accepted Step 04 KPGT freeze')
        origins = plan()
        paths = {run_path(root, origin) for origin in origins}
        actual = set((root / 'runs').glob('*/*/*'))
        if actual != paths:
            raise ValueError('The candidate does not contain the exact 48-run matrix')
        for origin in origins:
            verify_run(run_path(root, origin), origin, identity)
    return origins


def run(args):
    input_digest = digest(args.inputs / 'manifest.json')
    if input_digest != digest(ACCEPTED):
        raise ValueError('Inputs differ from the accepted Step 04 KPGT freeze')
    if args.smoke_rows is None and args.max_epochs != 1000:
        raise ValueError('Changing max epochs requires a bounded smoke run')
    worker = [sys.executable, str(HERE / 'learner.py')]
    identity = dict(
        input_manifest_sha256=input_digest,
        code={name: digest(HERE / name) for name in ('run.py', 'learner.py', 'model.py')},
        source={name: digest(args.source / 'src' / name)
                for name in ('model.py', 'train_TranSiGen.py', 'utils.py')},
        device=args.device, smoke_rows=args.smoke_rows, max_epochs=args.max_epochs,
    )
    args.candidate.mkdir(parents=True, exist_ok=True)
    with (args.candidate / '.driver.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        identity['runtime'] = json.loads(subprocess.check_output(
            [*worker, '--runtime', args.device], pass_fds=(lock.fileno(),),
        ))
        plan_path = args.candidate / 'plan.json'
        if plan_path.exists():
            if load_json(plan_path) != identity:
                raise ValueError('Candidate identity changed. Use a fresh candidate.')
        else:
            if (args.candidate / 'runs').exists():
                raise ValueError('Existing runs have no registered plan')
            with TemporaryDirectory(prefix='width-dose-plan-', dir=args.candidate) as work:
                temporary_plan = Path(work) / 'plan.json'
                save_json(temporary_plan, identity)
                temporary_plan.rename(plan_path)
        for origin in plan():
            if any(getattr(args, key) not in (None, origin[key])
                   for key in ('width', 'method', 'fold')):
                continue
            destination = run_path(args.candidate, origin)
            if destination.exists():
                verify_run(destination, origin, identity)
                continue
            destination.parent.mkdir(parents=True, exist_ok=True)
            print(json.dumps(origin), flush=True)
            with TemporaryDirectory(prefix='width-dose-', dir=args.candidate) as work:
                staging = Path(work) / 'run'
                staging.mkdir()
                save_json(staging / 'config.json', dict(origin=origin, identity=identity))
                command = [*worker, str(args.inputs), str(args.source), str(staging)]
                for phase in ('train', 'test'):
                    subprocess.run([*command, phase], check=True, pass_fds=(lock.fileno(),))
                save_json(staging / 'manifest.json', dict(
                    schema=SCHEMA, origin=origin, identity=identity,
                    files={name: digest(staging / name) for name in FILES},
                ))
                staging.rename(destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    commands.add_parser('plan')
    verification = commands.add_parser('verify')
    execution = commands.add_parser('run')
    for command in (execution, verification):
        command.add_argument('--candidate', type=Path, required=True)
    execution.add_argument('--inputs', type=Path, required=True)
    execution.add_argument('--source', type=Path, required=True)
    execution.add_argument('--device', default='cuda')
    execution.add_argument('--width', choices=WIDTHS)
    execution.add_argument('--method', choices=METHODS)
    execution.add_argument('--fold', type=int, choices=range(3))
    execution.add_argument('--smoke-rows', type=int, choices=range(2, 65), metavar='2..64')
    execution.add_argument('--max-epochs', type=int, choices=range(1, 1001),
                           default=1000, metavar='1..1000')
    args = parser.parse_args()
    if args.command == 'run':
        run(args)
    else:
        print(json.dumps(plan() if args.command == 'plan' else verify(args.candidate), indent=2))


if __name__ == '__main__':
    main()
