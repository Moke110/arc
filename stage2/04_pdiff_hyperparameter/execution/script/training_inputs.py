"""Frozen input identity and predictor allowlist for Step 04."""
import hashlib
import json
from pathlib import Path

from artifacts import read_json, rows, sha256, write_json

STEP03 = Path('stage2/03_modes_tenv')
REQUEST_SCHEMA = 'step04-predictor-v1'
INPUT_FIELDS = {'schema', 'dataset', 'fold', 'n_genes', 'target_kind', 'source_input_sha256',
                'corpus_sha256', 'preprocessing', 'split_sha256', 'gene_sha256',
                'predictor_files', 'file_sha256'}


def json_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def predictor_paths(fold):
    if type(fold) is not int or fold not in range(3):
        raise ValueError('unknown outer fold')
    prefix = f'fold_{fold}/'
    return {'genes': 'genes.csv', 'drugs': 'drugs.csv', 'morgan': 'morgan.npy',
            'controls': 'controls.npy', 'train_cases': prefix + 'train.csv',
            'validation_cases': prefix + 'validation.csv', 'test_cases': prefix + 'test_cases.csv',
            'train_expression': prefix + 'train_expression.npy',
            'validation_expression': prefix + 'validation_expression.npy'}


# Step 04 evaluation originally imported the tenv name.
tenv_paths = predictor_paths


def verify_files(root, paths, hashes):
    for key, relative in paths.items():
        path = root / relative
        if path.resolve() != path or not path.is_file() or sha256(path) != hashes[key]:
            raise ValueError(f'predictor file hash/path mismatch: {key}')


def prepare_tenv(runtime_root, input_root, output, fold):
    root = input_root.resolve()
    pin = read_json(runtime_root / STEP03 / 'execution/input_record/verification.json')['manifest_sha256']
    manifest_path = root / 'manifest.json'
    if manifest_path.resolve() != manifest_path or sha256(manifest_path) != pin:
        raise ValueError('tenv release differs from the frozen Step 03 identity')
    source = read_json(manifest_path)
    if (source['schema'], source['target_kind'], source['gene_count']) != ('tenv-input-v1', 'delta', 18164):
        raise ValueError('unexpected tenv release contract')
    paths = predictor_paths(fold)
    hashes = {key: source['files'][relative] for key, relative in paths.items()}
    split_hash = source['files']['folds.csv']
    verify_files(root, {**paths, 'split': 'folds.csv'}, {**hashes, 'split': split_hash})
    consumer = {'schema': 'step04-input-v1', 'dataset': 'tenv', 'fold': fold, 'n_genes': source['gene_count'],
                'target_kind': source['target_kind'], 'source_input_sha256': pin,
                'corpus_sha256': source.get('source_sha256'), 'preprocessing': source.get('expression'),
                'split_sha256': split_hash, 'gene_sha256': hashes['genes'],
                'predictor_files': paths, 'file_sha256': hashes}
    write_json(output.with_suffix('.input.json'), consumer)
    write_json(output, request_for(root, consumer))


def request_for(root, consumer):
    identity = {key: consumer[key] for key in ('dataset', 'fold', 'gene_sha256', 'split_sha256', 'source_input_sha256')}
    identity['input_sha256'] = json_hash(consumer)
    return {'schema': REQUEST_SCHEMA, 'identity': identity, 'input_root': str(root),
            'n_genes': consumer['n_genes'], 'target_kind': consumer['target_kind'],
            'predictor_files': {key: str(root / path) for key, path in consumer['predictor_files'].items()},
            'file_sha256': consumer['file_sha256'], 'input': consumer}


def verify_request(request):
    consumer = request['input']
    root = Path(request['input_root'])
    if (consumer.keys() != INPUT_FIELDS or consumer['schema'] != 'step04-input-v1'
            or consumer['target_kind'] != 'delta' or consumer['n_genes'] != 18164 or root != root.resolve()):
        raise ValueError('unexpected predictor request contract')
    if request != request_for(root, consumer):
        raise ValueError('request input identity does not match the consumer manifest')
    if consumer['dataset'] not in ('tenv', 'tahoe'):
        raise ValueError('unsupported predictor dataset')
    paths = predictor_paths(consumer['fold'])
    hashes = consumer['file_sha256']
    if (consumer['predictor_files'] != paths or hashes.keys() != paths.keys()
            or hashes['genes'] != consumer['gene_sha256']):
        raise ValueError('unexpected predictor file allowlist')
    verify_files(root, paths, hashes)
    if consumer['dataset'] == 'tahoe':
        verify_tahoe_roles(root, paths)


def verify_tahoe_roles(root, paths):
    drugs = {row['drug_key']: row['parent_connectivity'] for row in rows(root / paths['drugs'])}
    assigned = {}
    for role in ('train', 'validation', 'test'):
        for case in rows(root / paths[f'{role}_cases']):
            parent = drugs.get(case['drug_key'])
            if not parent:
                raise ValueError('Tahoe predictor case has no parent-connectivity block')
            if assigned.setdefault(parent, role) != role:
                raise ValueError('Tahoe parent-connectivity block crosses predictor roles')
