#!/usr/bin/env python3
"""Accept the independently replayed ADR0077 paired candidate; never rebuild axes.

Validate both sources and all input pins before moving live paths. Superseded
views are retained under the supplied recovery directory outside output. On an
installation error every moved path is restored. Recovery: move the new view-5k
out of output and move recovery/view-{2,3,4,5}k back; use the matching historical
reader from the sealed release when reconstructing that full-axis generation.
"""
import argparse
import json
from pathlib import Path
import shutil
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import public_api
from bases import release, view

STEP = Path('stage2/01_bases')
PROTECTED = ('gene-identity', 'raw_counts', 'corpus', 'tenv-hvg', 'disease-deg',
             'disease-fusion', 'gene-axis', 'fold_n_seed', 'palette')


def preservation(root):
    """Stream every protected scientific, frozen experiment and sealed-log byte."""
    paths = [root / STEP / 'output' / name for name in PROTECTED]
    paths += [root / STEP / name for name in ('figure', 'log')]
    for step in ('02_models_tenv-5k', '03_tsg_dose'):
        paths += [root / 'stage2' / step / name for name in
                  ('output', 'figure', 'log', 'execution/checkpoints')]
    records = {}
    for path in paths:
        if not path.exists():
            raise ValueError(f'missing preservation root: {path}')
        for file in sorted(path.rglob('*')) if path.is_dir() else [path]:
            if (file.is_file() and '__pycache__' not in file.parts and file.suffix != '.pyc'
                    and file.relative_to(root).as_posix() not in {
                        'stage2/01_bases/log/README.md', 'stage2/01_bases/log/acceptance.md'}):
                records[file.relative_to(root).as_posix()] = {
                    'bytes': file.stat().st_size, 'sha256': release.sha256(file)}
    return records


def validate(root, candidate):
    result = view.verify(candidate)
    manifest = release.read_json(candidate / 'manifest.json')
    for name, record in manifest['inputs'].items():
        if release.sha256(root / name) != record['sha256']:
            raise ValueError(f'candidate input generation differs: {name}')
    return result


def install(root, candidate, recovery):
    """Validate complete paired bytes, then replace with rollback on failure."""
    root, candidate, recovery = map(Path, (root, candidate, recovery))
    result = validate(root, candidate)
    output = root / STEP / 'output'
    if recovery.exists():
        raise ValueError('recovery destination must be fresh')
    recovery.mkdir(parents=True)
    staged = recovery / 'new-view-5k'
    # Finish the copy outside formal output before replacing any accepted path.
    shutil.copytree(candidate, staged)
    validate(root, staged)
    moved = []
    target = output / 'view-5k'
    try:
        for size in (2, 3, 4, 5):
            old = output / f'view-{size}k'
            if old.exists():
                old.rename(recovery / old.name)
                moved.append(old.name)
        staged.rename(target)
        validate(root, target)
    except BaseException:
        if target.exists():
            shutil.rmtree(target)
        for name in moved:
            (recovery / name).rename(output / name)
        raise
    return {'validation': result, 'retired_live_paths': [f'view-{i}k' for i in (2, 3, 4)],
            'recovery': recovery.relative_to(root).as_posix(),
            'manifest_sha256': release.sha256(target / 'manifest.json')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('snapshot', 'install', 'verify-preservation'))
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--candidate', type=Path)
    parser.add_argument('--recovery', type=Path)
    parser.add_argument('--record', type=Path, required=True)
    args = parser.parse_args()
    if args.action == 'snapshot':
        value = preservation(args.root)
    elif args.action == 'verify-preservation':
        before = release.read_json(args.record)
        for name, record in before.items():
            path = args.root / name
            if not path.is_file() or path.stat().st_size != record['bytes'] or release.sha256(path) != record['sha256']:
                raise ValueError(f'protected scientific byte changed: {name}')
        print(json.dumps({'status': 'PASS', 'files': len(before)}))
        return
    else:
        if args.candidate is None or args.recovery is None:
            parser.error('install requires --candidate and --recovery')
        value = install(args.root, args.candidate, args.recovery)
    release.write_json(args.record, value)

if __name__ == '__main__':
    main()
