"""Authenticate frozen Step 03 runs across the reporting-only migration."""
import hashlib
import json
from pathlib import Path
import subprocess

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / '.git').exists())


def accepts_frozen_run(path, previous, current):
    if previous == current:
        return True
    if {k: v for k, v in previous.items() if k != 'source_sha256'} != {
            k: v for k, v in current.items() if k != 'source_sha256'}:
        return False
    def committed_json(relative):
        physical = (REPO / relative).read_bytes()
        committed = subprocess.check_output(['git', '-C', str(REPO), 'show', f'dev:{relative}'])
        if physical != committed:
            raise ValueError('frozen compatibility evidence differs from committed dev bytes')
        return json.loads(committed)
    release = committed_json('stage2/evaluation/predictor-compatibility.json')
    if previous['source_sha256'] != release['previous_source_sha256']:
        return False
    for relative, digest in release['runtime_files'].items():
        if hashlib.sha256((REPO / relative).read_bytes()).hexdigest() != digest:
            return False
    inventory = committed_json('stage2/evaluation/inventory.json')
    path = Path(path).resolve()
    for item in inventory['runs']:
        original = item.get('run', {})
        if item['step'] == 3 and REPO / original.get('path', '') == path:
            return hashlib.sha256(path.read_bytes()).hexdigest() == original['sha256']
    return False
