"""Independently rescore twice, then atomically publish one Stage 2 release."""
import argparse
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import uuid

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / '.git').exists())
sys.path.insert(0, str(REPO))
from stage2.evaluation.inventory import STEPS, read, sha


def verify_existing(destination, figure_destination=None):
    if not destination.exists():
        return
    figure_destination = figure_destination or destination
    manifest = read(destination / 'manifest.json')
    expected = {*manifest['files'], 'manifest.json'}
    if figure_destination == destination:
        expected.update(manifest.get('figures', {}))
    if {p.name for p in destination.iterdir()} != expected or any(not p.is_file() for p in destination.iterdir()):
        raise ValueError('publication directory contains files outside the generated release')
    for directory, entries in ((destination, manifest['files']),
                               (figure_destination, manifest.get('figures', {}))):
        for name, digest in entries.items():
            if sha(directory / name) != digest:
                raise ValueError(f'published artifact changed: {name}')


def publish(step, inventory, labels_root=None, destination=None):
    destination = destination or STEPS[step] / 'output/retrieval'
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure_destination = STEPS[6] / 'output' if step == 6 and destination == STEPS[6] / 'output/retrieval' else destination
    verify_existing(destination, figure_destination)
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
        manifest = read(first / 'manifest.json')
        for name, digest in {**manifest['files'], **manifest['figures']}.items():
            if sha(first / name) != digest or sha(second / name) != digest:
                raise ValueError(f'independent replay differs: {name}')
        previous = destination.with_name(f'.retrieval-previous-{uuid.uuid4().hex}')
        old_figures, moved_figures = {}, []
        # Step 06's owning rule keeps result figures directly in output/.
        if figure_destination != destination:
            figure_destination.mkdir(parents=True, exist_ok=True)
            old = read(destination / 'manifest.json') if destination.exists() else {}
            for name in manifest['figures']:
                existing = figure_destination / name
                if existing.exists():
                    if name not in old.get('figures', {}):
                        raise ValueError(f'unowned result figure would be overwritten: {existing}')
                    backup = temporary / f'previous-{name}'
                    shutil.copy2(existing, backup)
                    old_figures[name] = backup
        exchanged = False
        try:
            if figure_destination != destination:
                for name in manifest['figures']:
                    (first / name).replace(figure_destination / name)
                    moved_figures.append(name)
            if destination.exists():
                destination.rename(previous)
            first.rename(destination)
            exchanged = True
        except BaseException:
            if not exchanged and previous.exists():
                previous.rename(destination)
            for name in moved_figures:
                target = figure_destination / name
                if name in old_figures:
                    old_figures[name].replace(target)
                else:
                    target.unlink(missing_ok=True)
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
