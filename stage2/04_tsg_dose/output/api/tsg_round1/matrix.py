"""Source-specific registered origins and validation gates."""
import json
import math
from pathlib import Path
import secrets

from tsg_dose.protocol import METHODS

from .inputs import SOURCE_POPULATIONS, sha256
from tenv5k.attempt_identity import attempt_binding
from .protocol import SCENARIOS, TRAINING_SEED


REGISTRATION_SCHEMA = 'tsg-dose-round1-axis-local-registration-v1'


def plan(scenario=None):
    if scenario is not None and scenario not in SCENARIOS:
        raise ValueError('unknown formal scenario')
    return [
        dict(source=source, scenario=name, fold=fold, method=method, seed=TRAINING_SEED)
        for source in SOURCE_POPULATIONS
        for name in SCENARIOS if scenario is None or name == scenario
        for fold in range(3)
        for method in METHODS
    ]


def select(scores, scenario):
    if scenario not in SCENARIOS or set(scores) != set(METHODS):
        raise ValueError('selection requires all seven unique validation origins')
    for score, parameters in scores.values():
        if (score is not None and (not math.isfinite(score) or abs(score) > 1 + 1e-12)
                or type(parameters) is not int or parameters < 1):
            raise ValueError('invalid validation score or parameter count')
    if scenario == 'cold-drug-dose':
        return None
    candidates = [(-score, parameters, METHODS.index(method), method)
                  for method, (score, parameters) in scores.items() if score is not None]
    if not candidates:
        raise ValueError('no defined validation PCC')
    return min(candidates)[-1]


def run_path(candidate, origin):
    return Path(candidate) / 'runs' / origin['source'] / origin['scenario'] / f"fold-{origin['fold']}" / origin['method']


def registration_path(candidate, origin):
    return Path(candidate) / 'registrations' / origin['source'] / origin['scenario'] / f"fold-{origin['fold']}" / (origin['method'] + '.json')


def selection_path(candidate, origin):
    return Path(candidate) / 'selection' / origin['source'] / origin['scenario'] / f"fold-{origin['fold']}.json"


def register(candidate, inputs, origin, device):
    """Write one immutable preregistration before starting a fresh origin."""
    if origin not in plan():
        raise ValueError('invalid formal run identity')
    path = registration_path(candidate, origin)
    path.parent.mkdir(parents=True, exist_ok=True)
    value = dict(schema=REGISTRATION_SCHEMA, run=origin, origin='formal', device=str(device),
                 attempt_nonce=secrets.token_hex(32),
                 input_manifest_sha256=sha256(Path(inputs) / origin['source'] / 'manifest.json'))
    with path.open('x') as handle:
        json.dump(value, handle, sort_keys=True, indent=2)
        handle.write('\n')
    return path


def check_registration(registration, inputs, origin, device):
    if (registration.get('schema') != REGISTRATION_SCHEMA
            or registration.get('origin') != 'formal'
            or registration.get('run') != origin or origin not in plan()
            or registration.get('input_manifest_sha256') != sha256(Path(inputs) / 'manifest.json')
            or registration.get('device') != str(device)
            or len(registration.get('attempt_nonce', '')) != 64):
        raise ValueError('formal registration identity mismatch')


def read_json(path):
    return json.loads(Path(path).read_text())


def validate_origins(expected, found):
    """Reject missing, duplicate and foreign origins before aggregation."""
    def key(origin):
        return tuple(origin[name] for name in ('source', 'scenario', 'fold', 'method', 'seed'))
    expected_keys = [key(origin) for origin in expected]
    found_keys = [key(origin) for origin in found]
    if len(found_keys) != len(expected_keys) or set(found_keys) != set(expected_keys):
        raise ValueError('missing, duplicate or unexpected formal origins')


def lock_selection(candidate, inputs, origin):
    """Record source/fold-local validation evidence before any outer test."""
    path = selection_path(candidate, origin)
    value = selection_evidence(candidate, inputs, origin, before_test=not path.exists())
    if path.exists():
        if read_json(path) != value:
            raise ValueError('sealed selection evidence changed')
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('x') as handle:
            json.dump(value, handle, sort_keys=True, indent=2, allow_nan=False)
            handle.write('\n')
    return path


def selection_evidence(candidate, inputs, origin, *, before_test=False):
    scores = {}
    evidence = {}
    input_digest = sha256(Path(inputs) / 'manifest.json')
    origins = [item for item in plan(origin['scenario'])
               if item['source'] == origin['source'] and item['fold'] == origin['fold']]
    expected_paths = {item['method'] for item in origins}
    fold_path = run_path(candidate, origin).parent
    if not expected_paths or {path.name for path in fold_path.iterdir()} != expected_paths:
        raise ValueError('selection requires exactly seven fold origin directories')
    allowed_statuses = {'validation-complete'}
    if not before_test:
        allowed_statuses.update(('test-awaiting-replay', 'complete'))
    for item in origins:
        method = item['method']
        run = run_path(candidate, item)
        registration = read_json(registration_path(candidate, item))
        config = read_json(run / 'config.json')
        terminal = read_json(run / 'terminal.json')
        check_registration(registration, inputs, item, config['device'])
        if (read_json(run / 'registration.json') != registration
                or config['execution_binding'] != attempt_binding(registration)
                or any(config[key] != value for key, value in item.items())
                or terminal['smoke'] or terminal['status'] not in allowed_statuses
                or terminal['validation_replay']['attempt_nonce'] != registration['attempt_nonce']
                or terminal['validation_replay']['process_token'] == terminal['training_process_token']):
            raise ValueError('validation origin registration or lifecycle mismatch')
        manifest_path = run / 'validation_manifest.json'
        manifest = read_json(manifest_path)
        required = {'checkpoint.pt', 'config.json', 'environment.json', 'learning_curve.csv',
                    'registration.json', 'validation_identity.csv', 'validation_metrics.json',
                    'validation_per_sample.csv', 'validation_prediction.npy'}
        if (manifest['schema'] != 'tsg-dose-round1-axis-local-validation-v1'
                or manifest['input_manifest_sha256'] != input_digest
                or set(manifest['outputs']) != required
                or manifest['parameter_count'] != terminal['parameter_count']):
            raise ValueError('validation evidence identity mismatch')
        for name in ('config.json', 'registration.json', 'validation_metrics.json'):
            if sha256(run / name) != manifest['outputs'][name]:
                raise ValueError(f'validation artifact changed: {name}')
        scores[method] = (read_json(run / 'validation_metrics.json')['pcc'],
                          manifest['parameter_count'])
        evidence[method] = dict(validation_manifest_sha256=sha256(manifest_path),
                                validation_pcc=scores[method][0], parameter_count=scores[method][1])
    return dict(schema='tsg-dose-round1-axis-local-selection-v1',
                source=origin['source'], scenario=origin['scenario'], fold=origin['fold'],
                input_manifest_sha256=input_digest, evidence=evidence,
                selected_method=select(scores, origin['scenario']))
