"""Independently rescore twice, then atomically publish one Stage 2 release."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import uuid

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from stage2.evaluation.inventory import STEPS, read, sha


def verify_existing(destination):
    if not destination.exists():
        return
    manifest = read(destination / 'manifest.json')
    expected = {*manifest['files'], 'manifest.json'}
    if {p.name for p in destination.iterdir()} != expected or any(not p.is_file() for p in destination.iterdir()):
        raise ValueError('publication directory contains files outside the generated release')
    for name, digest in manifest['files'].items():
        if sha(destination / name) != digest:
            raise ValueError(f'published artifact changed: {name}')


def publish(step, inventory, labels_root=None, destination=None):
    destination = destination or STEPS[step] / 'output/retrieval'
    destination.parent.mkdir(parents=True, exist_ok=True)
    verify_existing(destination)
    with tempfile.TemporaryDirectory(dir=destination.parent, prefix='.retrieval-') as temp:
        temporary = Path(temp)
        for name in ('score', 'replay'):
            command = [sys.executable, str(Path(__file__).with_name('evaluate.py')), 'score', '--step', str(step),
                       '--inventory', str(inventory), '--destination', str(temporary / name)]
            if labels_root:
                command.extend(['--labels-root', str(labels_root)])
            subprocess.run(command, check=True)
        first, second = temporary / 'score', temporary / 'replay'
        if sha(first / 'manifest.json') != sha(second / 'manifest.json'):
            raise ValueError('independent process replay differs')
        for name, digest in read(first / 'manifest.json')['files'].items():
            if sha(second / name) != digest:
                raise ValueError(f'independent replay differs: {name}')
        previous = destination.with_name(f'.retrieval-previous-{uuid.uuid4().hex}')
        if destination.exists():
            destination.rename(previous)
        try:
            first.rename(destination)
        except BaseException:
            if previous.exists():
                previous.rename(destination)
            raise
        if previous.exists():
            shutil.rmtree(previous)
    print(f'Step {step:02d}: published three-metric tenv release after independent replay', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--step', type=int, choices=STEPS, action='append', help='repeat for multiple steps; default all')
    parser.add_argument('--inventory', type=Path, default=Path(__file__).with_name('inventory.json'))
    parser.add_argument('--labels-root', type=Path)
    args = parser.parse_args()
    for step in args.step or STEPS:
        publish(step, args.inventory, args.labels_root)


if __name__ == '__main__':
    main()
