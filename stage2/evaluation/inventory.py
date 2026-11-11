"""Freeze predictor-safe identities and read-only prediction bindings for rescoring."""
from pathlib import Path
import hashlib
import json
import pandas as pd
from functools import lru_cache

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / 'stage2'
SEEDS = (17, 42, 2026)
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


def verify(path=None):
    path = path or Path(__file__).with_name('inventory.json')
    inventory = read(path)
    seen = set()
    for run in inventory['runs']:
        key = tuple(run[k] for k in ('step', 'dataset', 'regime', 'role', 'model', 'fold', 'seed'))
        if key in seen:
            raise ValueError(f'duplicate evaluation identity: {key}')
        seen.add(key)
        for value in run.values():
            if isinstance(value, dict) and 'path' in value and 'sha256' in value:
                binding(REPO / value['path'], value['sha256'])
    return inventory


if __name__ == '__main__':
    inventory = verify()
    print(f"Verified {len(inventory['runs'])} frozen prediction/reference identities")
