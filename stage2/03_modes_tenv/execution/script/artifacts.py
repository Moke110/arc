"""Persisted run identities and prediction verification shared by CLI/scheduler."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

SEEDS = (17, 42, 2026)
TERMINAL = {'completed', 'failed', 'timeout', 'excluded'}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.partial')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def rows(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, records, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    records = list(records)
    fields = fields or list(records[0])
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def run_directory(output_root, model, fold, seed):
    if not model or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in model):
        raise ValueError('invalid model identifier')
    if fold not in range(3) or seed not in SEEDS:
        raise ValueError('invalid fold/training seed')
    return Path(output_root) / 'runs' / model / f'fold_{fold}' / f'seed_{seed}'


def verify_input_files(input_root):
    root=Path(input_root).resolve()
    manifest=read_json(root/'manifest.json')
    if manifest.get('synthetic'):
        return manifest
    if manifest.get('gene_count')!=18164 or manifest.get('drug_count')!=460 or manifest.get('training_seeds')!=list(SEEDS):
        raise ValueError('formal input release requires 460 drugs, 18164 genes and fixed seeds')
    if not manifest.get('files'):
        raise ValueError('input release missing file hashes')
    for relative,digest in manifest['files'].items():
        path=(root/relative).resolve()
        if not path.is_relative_to(root) or sha256(path)!=digest:
            raise ValueError(f'input hash/path mismatch: {relative}')
    return manifest


def run_identity(input_root, model, fold, seed, command=None):
    root = Path(input_root).resolve()
    sources = [Path(__file__).parent / name for name in ('artifacts.py','benchmark.py','baselines.py','metrics.py','plots.py')]
    source_hash = hashlib.sha256(''.join(f'{p.name}:{sha256(p)}\n' for p in sources).encode()).hexdigest()
    # Command and environment identity are explicit and checked on resume.
    import platform
    import sys
    environment = {'python':sys.version, 'executable':sys.executable, 'numpy':np.__version__, 'platform':platform.platform()}
    config = {'command':command, 'model':model, 'target_kind':'delta', 'command_evidence':command_evidence(command)}
    return {'model':model, 'fold':int(fold), 'seed':int(seed),
            'input_sha256':sha256(root / 'manifest.json'),
            'split_sha256':sha256(root / 'folds.csv'), 'gene_sha256':sha256(root / 'genes.csv'),
            'source_sha256':source_hash,
            'environment_sha256':hashlib.sha256(json.dumps(environment,sort_keys=True).encode()).hexdigest(),
            'config_sha256':hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()}


def validate_prediction(run_dir, cases, n_genes):
    run_dir = Path(run_dir)
    expected = [case['case_id'] for case in cases]
    if len(expected) != len(set(expected)):
        raise ValueError('duplicate input case identity')
    identifiers = read_json(run_dir / 'case_ids.json')
    if identifiers != expected:
        raise ValueError('prediction case identities missing, duplicate or misordered')
    prediction = np.load(run_dir / 'prediction.npy', mmap_mode='r', allow_pickle=False)
    if prediction.shape != (len(cases), n_genes):
        raise ValueError('wrong prediction shape')
    for start in range(0,len(cases),32):
        if not np.isfinite(prediction[start:start+32]).all():
            raise ValueError('nonfinite prediction')
    return prediction


def command_evidence(command):
    """Pin external adapter plus declared upstream, environment and recipe files.

    A command script foo.py must have foo.identity.json beside it, containing
    nonempty source_files/environment_files/config_files lists. Relative paths
    resolve against that sidecar. Checkpoints declare every binary in blobs.
    """
    if command is None:
        return {}
    scripts=[Path(arg).resolve() for arg in command if Path(arg).is_file() and str(arg).endswith('.py')]
    if not scripts:
        raise ValueError('external command needs an adapter script identity sidecar')
    evidence={}
    for script in scripts:
        sidecar=script.with_suffix('.identity.json')
        if not sidecar.is_file():
            raise ValueError('external adapter identity sidecar missing')
        declaration=read_json(sidecar)
        evidence[str(script)]=sha256(script)
        evidence[str(sidecar)]=sha256(sidecar)
        for category in ('source_files','environment_files','config_files'):
            if not declaration.get(category): raise ValueError(f'identity sidecar missing {category}')
            for name in declaration[category]:
                path=Path(name)
                if not path.is_absolute(): path=sidecar.parent/path
                evidence[str(path.resolve())]=sha256(path)
    for argument in command:
        path=Path(argument)
        if path.is_file(): evidence[str(path.resolve())]=sha256(path)
    return evidence


def checkpoint_artifacts(directory, checkpoint, require_blobs=False):
    directory=Path(directory).resolve()
    blobs=checkpoint.get('blobs',{})
    if require_blobs and not blobs:
        raise ValueError('external checkpoint must declare checkpoint blobs')
    for name,digest in blobs.items():
        path=(directory/name).resolve()
        if not path.is_relative_to(directory) or sha256(path)!=digest:
            raise ValueError('checkpoint blob checksum/path mismatch')
    return {name:sha256(directory/name) for name in blobs}


def validate_run(input_root, output_root, model, fold, seed, command=None):
    root = Path(input_root)
    verify_input_files(root)
    directory = run_directory(output_root,model,fold,seed)
    try:
        run = read_json(directory / 'run.json')
        if run['status'] != 'completed' or run['identity'] != run_identity(root,model,fold,seed,command):
            raise ValueError('stale run identity/status')
        for name in ('prediction.npy','case_ids.json','checkpoint.json','request.json'):
            if run['artifacts'][name] != sha256(directory / name):
                raise ValueError(f'artifact checksum mismatch: {name}')
        validate_prediction(directory,rows(root / f'fold_{fold}' / 'test_cases.csv'),len(rows(root / 'genes.csv')))
        checkpoint = read_json(directory / 'checkpoint.json')
        checkpoint_artifacts(directory,checkpoint,require_blobs=command is not None)
        if checkpoint.get('identity') != run['identity']:
            raise ValueError('checkpoint identity mismatch')
        return run
    except (OSError,KeyError,TypeError,json.JSONDecodeError) as exc:
        raise ValueError(f'invalid completed run: {exc}') from exc
