"""Freeze predictor-safe identities and read-only prediction bindings for rescoring."""
from pathlib import Path
import hashlib
import json
from functools import lru_cache

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / '.git').exists())
ROOT = REPO / 'stage2'
STEPS = {n: next(ROOT.glob(f'{n:02d}_*')) for n in (3, 4, 5, 6, 9)}


@lru_cache(maxsize=1024)
def _digest(path, signature):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def sha(path):
    path = Path(path)
    stat = path.stat()
    signature = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
    return _digest(str(path), signature)


def read(path):
    return json.loads(Path(path).read_text())


def binding(path, expected=None):
    path = Path(path).absolute()
    if path.resolve() != path or not path.is_relative_to(REPO):
        raise ValueError(f'unsafe artifact: {path}')
    digest = sha(path)
    if expected and digest != expected:
        raise ValueError(f'changed frozen artifact: {path}')
    return {'path': str(path.relative_to(REPO)), 'sha256': digest}


def iter_bindings(value):
    if isinstance(value, dict):
        if 'path' in value and 'sha256' in value:
            yield value
        else:
            for item in value.values():
                yield from iter_bindings(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from iter_bindings(item)


def load_inventory(path):
    inventory = read(path)
    fields = ('step', 'dataset', 'regime', 'role', 'model', 'fold', 'seed')
    identities = [tuple(run[key] for key in fields) for run in inventory['runs']]
    expected = [tuple(value) for value in inventory['expected_identities']]
    if (inventory.get('schema') != 'stage2-retrieval-inventory-v1'
            or len(identities) != len(set(identities)) or len(expected) != len(set(expected))
            or set(identities) != set(expected)
            or any(run['dataset'] != 'tenv' for run in inventory['runs'])):
        raise ValueError('incomplete or invalid frozen tenv prediction inventory')
    return inventory


def verify(path=None):
    path = path or Path(__file__).with_name('inventory.json')
    inventory = load_inventory(path)
    seen = set()
    for run in inventory['runs']:
        key = tuple(run[k] for k in ('step', 'dataset', 'regime', 'role', 'model', 'fold', 'seed'))
        if key in seen:
            raise ValueError(f'duplicate evaluation identity: {key}')
        seen.add(key)
        for value in iter_bindings(run):
            binding(REPO / value['path'], value['sha256'])
    return inventory


if __name__ == '__main__':
    inventory = verify()
    print(f"Verified {len(inventory['runs'])} frozen prediction/reference identities")
