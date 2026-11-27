"""Read and publish Stage 2 ECFP4 fingerprints without changing sealed inputs."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil

import numpy as np

ECFP4_MLP = 'ECFP4_MLP'
SCHEMA = 'stage2-ecfp4-v1'
LEGACY_SCHEMAS = {
    'tenv5k-inputs-v2': 'includeChirality',
    'tsg-dose-round1-axis-local-inputs-v1': 'include_chirality',
}
REPRESENTATION = {
    'atom_invariants': 'connectivity',
    'binary': True,
    'bits': 2048,
    'include_chirality': True,
    'radius': 2,
}


def _sha256(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def _read_json(path):
    return json.loads(path.read_text())


def _verified(directory, manifest, name):
    path = directory / name
    if _sha256(path) != manifest['outputs'].get(name):
        raise ValueError(f'ECFP4 input digest mismatch: {name}')
    return path


def _read_contract(directory, manifest):
    schema = manifest.get('schema')
    if schema == SCHEMA:
        parameters = manifest['ecfp4']
        if (manifest.get('source_schema') not in LEGACY_SCHEMAS
                or not re.fullmatch(r'[0-9a-f]{64}', manifest.get('source_manifest_sha256', ''))):
            raise ValueError('invalid ECFP4 source identity')
        fingerprint_name = 'ecfp4.npy'
        if set(manifest['outputs']) != {fingerprint_name, 'parents.csv'}:
            raise ValueError('ECFP4 manifest coverage mismatch')
    elif schema in LEGACY_SCHEMAS:
        config = _read_json(_verified(directory, manifest, 'config.json'))['morgan']
        parameters = {
            'atom_invariants': 'connectivity',
            'binary': config.get('binary', True),
            'bits': config['bits'],
            'include_chirality': config[LEGACY_SCHEMAS[schema]],
            'radius': config['radius'],
            'rdkit': config['rdkit'],
        }
        fingerprint_name = 'morgan.npy'
    else:
        raise ValueError(f'unknown Stage 2 ECFP4 input schema: {schema}')

    if (any(parameters.get(key) != value for key, value in REPRESENTATION.items())
            or not isinstance(parameters.get('rdkit'), str) or not parameters['rdkit']):
        raise ValueError('Stage 2 ECFP4 representation mismatch')
    fingerprints = _verified(directory, manifest, fingerprint_name)
    parents_path = _verified(directory, manifest, 'parents.csv')
    with parents_path.open(newline='') as stream:
        parents = [row['drug_parent_inchikey'] for row in csv.DictReader(stream)]
    bits = np.load(fingerprints, mmap_mode='r')
    if (not parents or len(set(parents)) != len(parents) or not all(parents)
            or bits.dtype != np.uint8 or bits.shape != (len(parents), REPRESENTATION['bits'])
            or np.any(bits > 1)):
        raise ValueError('Stage 2 ECFP4 parent identity, shape or binary dtype mismatch')
    return bits, fingerprints, parents_path, parameters


def open_ecfp4(directory, *, source_manifest=None):
    """Verify the binary matrix and optionally bind it to the original input freeze."""
    directory = Path(directory)
    manifest = _read_json(directory / 'manifest.json')
    if source_manifest is not None:
        source_manifest = Path(source_manifest)
        if (manifest.get('schema') != SCHEMA
                or manifest.get('source_manifest_sha256') != _sha256(source_manifest)):
            raise ValueError('ECFP4 source manifest mismatch')
        original = _read_json(source_manifest)
        if (manifest['source_schema'] != original['schema']
                or manifest['outputs'].get('ecfp4.npy') != original['outputs'].get('morgan.npy')
                or manifest['outputs'].get('parents.csv') != original['outputs'].get('parents.csv')):
            raise ValueError('ECFP4 source fingerprint or parent order mismatch')
    return _read_contract(directory, manifest)[0]


def publish_ecfp4(input_directory, candidate_directory):
    """Prepare a canonical fingerprint contract from a verified frozen source."""
    source = Path(input_directory).resolve()
    candidate = Path(candidate_directory).resolve()
    if source == candidate or source in candidate.parents:
        raise ValueError('ECFP4 candidate must be outside the frozen input directory')
    if candidate.exists() and any(candidate.iterdir()):
        raise ValueError('ECFP4 publication requires a new empty candidate')
    manifest_path = source / 'manifest.json'
    manifest_sha256 = _sha256(manifest_path)
    manifest = _read_json(manifest_path)
    if manifest.get('schema') not in LEGACY_SCHEMAS:
        raise ValueError('ECFP4 publication requires a frozen Stage 2 input generation')
    _, fingerprints, parents, parameters = _read_contract(source, manifest)
    candidate.mkdir(parents=True, exist_ok=True)
    copied = {'ecfp4.npy': fingerprints, 'parents.csv': parents}
    outputs = {}
    for name, original in copied.items():
        shutil.copyfile(original, candidate / name)
        outputs[name] = manifest['outputs'][original.name]
        if _sha256(candidate / name) != outputs[name]:
            raise ValueError(f'ECFP4 copy digest mismatch: {name}')
    if _sha256(manifest_path) != manifest_sha256:
        raise ValueError('frozen input manifest changed during ECFP4 publication')
    (candidate / 'manifest.json').write_text(json.dumps({
        'schema': SCHEMA,
        'ecfp4': parameters,
        'source_manifest_sha256': manifest_sha256,
        'source_schema': manifest['schema'],
        'outputs': outputs,
    }, indent=2, sort_keys=True) + '\n')
    open_ecfp4(candidate, source_manifest=manifest_path)
    return candidate


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    args = parser.parse_args()
    publish_ecfp4(args.inputs, args.candidate)
