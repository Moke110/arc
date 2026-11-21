"""Run and verify the frozen tenv KPGT TranSiGen comparison."""
import argparse
from collections import defaultdict
from contextlib import contextmanager
import fcntl
import json
from pathlib import Path
import subprocess
import sys
import time

STEP = Path(__file__).resolve().parents[2]
EXECUTION = STEP / 'execution'
sys.path.insert(0, str(STEP / 'output/api'))
from tsg_kpgt.inputs import Inputs, sha256
from tsg_kpgt.learner import environment, runtime_identity, verify_candidate, write_json
from tsg_kpgt.matrix import (
    check_registration, lock_selection, plan, read_json, register, registration_path,
    run_path, selection_evidence, selection_path, validate_origins,
)
from tsg_round1.model import SOURCE_DIGEST
from tsg_round1.protocol import SCENARIOS
from tenv5k.attempt_identity import attempt_binding

TRANSIGEN = EXECUTION / 'kpgt/transigen.py'
PLAN_SCHEMA = 'tsg-dose-tenv-kpgt-formal-plan-v1'
BASELINES = ('KPGT_MLP', 'Mean_cell')
EXPECTED_OUTPUTS = {
    'checkpoint.pt', 'config.json', 'environment.json', 'learning_curve.csv', 'registration.json',
    'terminal.json', 'test_metrics.json', 'test_per_sample.csv', 'test_prediction.npy',
    'validation_identity.csv', 'validation_manifest.json', 'validation_metrics.json',
    'validation_per_sample.csv', 'validation_prediction.npy',
}


def next_action(status, phase):
    transitions = {
        'validation': {None: 'train', 'validation-trained-awaiting-replay': 'replay',
                       'validation-complete': None, 'test-awaiting-replay': None, 'complete': None},
        'test': {'validation-complete': 'evaluate', 'test-awaiting-replay': 'replay', 'complete': None},
    }
    if status not in transitions[phase]:
        raise ValueError(f'cannot resume {phase} from {status}')
    return transitions[phase][status]


def planned_runs():
    return dict(
        transigen=plan(),
        baselines=[dict(source='tenv', scenario=scenario, fold=fold, method=method, seed=42)
                   for scenario in SCENARIOS for fold in range(3) for method in BASELINES],
    )


def frozen_plan(args):
    """Bind the accepted input generation, execution code and both scenarios."""
    input_digest = sha256(args.inputs / 'manifest.json')
    accepted = STEP / 'log/274-kpgt-smoke/input-manifest.json'
    if sha256(accepted) != input_digest:
        raise ValueError('registered generation differs from accepted #274 inputs')
    native_path = args.source / 'src/model.py'
    if sha256(native_path) != SOURCE_DIGEST:
        raise ValueError('unregistered native TranSiGen source')
    files = (Path(__file__), TRANSIGEN, EXECUTION / 'capped.sh', accepted,
             STEP / 'log/274-kpgt-smoke/acceptance.json', args.source / 'src/train_TranSiGen.py')
    return dict(schema=PLAN_SCHEMA, origins=plan(), inputs=str(args.inputs), device=args.device,
                input_manifest_sha256=input_digest, native_source=str(args.source),
                native_model_sha256=SOURCE_DIGEST, learner_runtime=runtime_identity(),
                environment=environment(args.device),
                code_evidence_sha256={path.name: sha256(path) for path in files},
                max_epochs=1000, patience=50, validation_seed=22617)


def check_plan(args, *, create=False):
    current = frozen_plan(args)
    path = args.candidate / 'formal_plan.json'
    if not path.exists() and create:
        if (args.candidate / 'runs').exists() or (args.candidate / 'registrations').exists():
            raise ValueError('missing sealed plan for existing origins')
        write_json(path, current)
    elif read_json(path) != current:
        raise ValueError('sealed formal plan, runtime or input identity changed')
    return current


@contextmanager
def candidate_lock(root):
    root.mkdir(parents=True, exist_ok=True)
    with (root / '.driver.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield lock.fileno()


def folds(scenario):
    groups = defaultdict(list)
    for origin in plan(scenario):
        groups[origin['source'], origin['scenario'], origin['fold']].append(origin)
    return [groups[key] for key in sorted(groups, key=lambda key: (key[2], key[0], key[1]))]


def status(run):
    if (run / 'terminal.json').is_file():
        return read_json(run / 'terminal.json')['status']
    return 'interrupted' if run.exists() else None


def memory_preflight():
    for line in Path('/proc/meminfo').read_text().splitlines():
        if line.startswith('MemAvailable:'):
            available_kib = int(line.split()[1])
            if available_kib < 12 * 1024 * 1024:
                raise RuntimeError('available memory below 12 GiB; stop before next formal phase')
            return available_kib
    raise RuntimeError('missing required MemAvailable preflight')


def invoke(args, origin, action, *, phase, lock_fd, selection=None):
    available = memory_preflight()
    root = args.candidate
    run = run_path(root, origin)
    registration = registration_path(root, origin)
    if action == 'train' and not registration.exists():
        register(root, args.inputs, origin, args.device)
    check_registration(read_json(registration), args.inputs, origin, args.device)
    command = [
        sys.executable, str(TRANSIGEN), '--inputs', str(args.inputs),
        '--source', str(args.source), '--candidate', str(run), '--device', args.device,
        '--method', origin['method'], '--scenario', origin['scenario'],
        '--fold', str(origin['fold']), '--seed', str(origin['seed']),
    ]
    phase_options = {
        'train': ['--registration', str(registration), '--defer-replay'],
        'evaluate': ['--evaluate-only', '--selection', str(selection), '--defer-replay'],
        'replay': ['--replay-only'],
    }
    attempt = root / 'resources' / run.relative_to(root / 'runs') / f'{phase}-{action}-{time.time_ns()}'
    attempt.mkdir(parents=True)
    started = time.time()
    write_json(attempt / 'phase.json', dict(origin=origin, action=action, phase=phase,
               available_ram_kib=available, started_utc_unix_seconds=started,
               registration_sha256=sha256(registration)))
    command = ['/usr/bin/time', '-v', '-o', str(attempt / 'time.txt'), *command, *phase_options[action]]
    print(json.dumps(dict(event='start', origin=origin, action=action, phase=phase), sort_keys=True), flush=True)
    with (attempt / 'subprocess.log').open('w') as output:
        result = subprocess.run(command, stdout=output, stderr=subprocess.STDOUT, check=False,
                                pass_fds=(lock_fd,))
    write_json(attempt / 'result.json', dict(returncode=result.returncode,
               elapsed_seconds=time.time() - started, time_sha256=sha256(attempt / 'time.txt')))
    result.check_returncode()
    print(json.dumps(dict(event='complete', origin=origin, action=action, phase=phase), sort_keys=True), flush=True)


def run_matrix(args):
    """Each fold trains all seven validations before unlocking any outer label."""
    with candidate_lock(args.candidate) as lock_fd:
        Inputs(args.inputs)
        check_plan(args, create=True)
        for origins in folds(args.scenario):
            for origin in origins:
                while action := next_action(status(run_path(args.candidate, origin)), 'validation'):
                    invoke(args, origin, action, phase='validation', lock_fd=lock_fd)
            first = origins[0]
            selection = lock_selection(args.candidate, args.inputs, first)
            for origin in origins:
                while action := next_action(status(run_path(args.candidate, origin)), 'test'):
                    invoke(args, origin, action, phase='test', lock_fd=lock_fd, selection=selection)


def verify(args, scenario):
    """Verify formal manifests, sealed validation choices and fresh replay proof."""
    frozen = check_plan(args)
    root = args.candidate
    expected = plan(scenario)
    found = []
    for origins in folds(scenario):
        first = origins[0]
        inputs = args.inputs
        choice_path = selection_path(root, first)
        choice = read_json(choice_path)
        if choice != selection_evidence(root, inputs, first):
            raise ValueError('sealed selection receipt changed')
        choice_digest = sha256(choice_path)
        for origin in origins:
            path = run_path(root, origin)
            registration = read_json(registration_path(root, origin))
            check_registration(registration, inputs, origin, frozen['device'])
            bound = attempt_binding(registration)
            manifest = verify_candidate(path)
            config = read_json(path / 'config.json')
            terminal = read_json(path / 'terminal.json')
            if (set(manifest['outputs']) != EXPECTED_OUTPUTS
                    or manifest['input_manifest_sha256'] != frozen['input_manifest_sha256']
                    or manifest['execution_binding'] != bound or config['execution_binding'] != bound
                    or read_json(path / 'registration.json') != registration
                    or config['runtime_identity'] != frozen['learner_runtime']
                    or config['native_source_sha256'] != frozen['native_model_sha256']
                    or read_json(path / 'environment.json') != frozen['environment']
                    or any(config[key] != frozen[key] for key in ('max_epochs', 'patience', 'validation_seed'))
                    or config['smoke_rows'] is not None
                    or terminal['status'] != 'complete' or terminal['smoke'] or not terminal['weights_updated']
                    or terminal['test_selection_sha256'] != choice_digest):
                raise ValueError('formal origin identity or terminal state changed')
            validation = read_json(path / 'validation_manifest.json')
            if any(manifest['outputs'].get(name) != digest for name, digest in validation['outputs'].items()):
                raise ValueError('sealed validation artifacts changed before test')
            tokens = [terminal['training_process_token'], terminal['test_evaluation_process_token']]
            for phase in ('validation', 'test'):
                replay = terminal[phase + '_replay']
                if replay['attempt_nonce'] != bound['attempt_nonce'] or not replay['process_token']:
                    raise ValueError('formal replay evidence mismatch')
                tokens.append(replay['process_token'])
            if any(not token for token in tokens) or len(set(tokens)) != len(tokens):
                raise ValueError('formal replay was not run in separate processes')
            found.append(dict(origin, terminal=terminal, run_manifest_sha256=sha256(path / 'manifest.json'),
                              selection_sha256=choice_digest, selected_method=choice['selected_method']))
    validate_origins(expected, [record['origin'] for record in found])
    for tree in ('runs', 'registrations'):
        top = root / tree
        if {item.name for item in top.iterdir()} != {'tenv'}:
            raise ValueError('foreign source in formal candidate')
        source = top / 'tenv'
        if {item.name for item in source.iterdir()} - set(SCENARIOS):
            raise ValueError('foreign scenario in formal candidate')
        for label in SCENARIOS if scenario is None else (scenario,):
            base = source / label
            canonical = {f'fold-{fold}' for fold in range(3)}
            if {item.name for item in base.iterdir()} != canonical:
                raise ValueError('formal fold inventory mismatch')
            for fold in range(3):
                planned = [item for item in expected if item['scenario'] == label and item['fold'] == fold]
                filenames = {item['method'] + ('.json' if tree == 'registrations' else '')
                             for item in planned}
                if {item.name for item in (base / f'fold-{fold}').iterdir()} != filenames:
                    raise ValueError('formal run inventory mismatch')
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_subparsers(dest='action', required=True)
    actions.add_parser('plan')
    for action in ('run', 'verify'):
        command = actions.add_parser(action)
        for name in ('inputs', 'source', 'candidate'):
            command.add_argument('--' + name, type=Path, required=True)
        command.add_argument('--scenario', choices=SCENARIOS)
        command.add_argument('--device', default='cuda')
    args = parser.parse_args()
    if args.action == 'plan':
        print(json.dumps(planned_runs(), indent=2))
        return
    for name in ('inputs', 'source', 'candidate'):
        setattr(args, name, getattr(args, name).resolve())
    if args.action == 'run':
        run_matrix(args)
    else:
        with candidate_lock(args.candidate):
            print(json.dumps(verify(args, args.scenario)))


if __name__ == '__main__':
    main()
