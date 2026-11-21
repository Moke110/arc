"""Prepare a tenv-only input candidate from the accepted Steps 03 and 04 freezes."""

import argparse
import json
from pathlib import Path
import shutil
import sys
from tempfile import TemporaryDirectory

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP / 'output/api'))
from tsg_kpgt.inputs import (
    ACCEPTED_REPRESENTATIONS,
    ACCEPTED_TENV,
    INPUT_SCHEMA,
    Inputs,
    sha256,
)
from tsg_round1.inputs import REQUIRED


def prepare(tenv, representations, candidate):
    if candidate.exists():
        raise ValueError('KPGT inputs require a fresh candidate path')
    for source in (tenv, representations):
        if source.resolve() in candidate.resolve().parents:
            raise ValueError('KPGT candidate cannot reside inside a frozen source')
    if (
        sha256(tenv / 'manifest.json') != ACCEPTED_TENV
        or sha256(representations / 'manifest.json') != ACCEPTED_REPRESENTATIONS
    ):
        raise ValueError(
            'expected the accepted tenv input and Step 03 representation freezes'
        )
    candidate.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix='.kpgt-inputs-', dir=candidate.parent) as root:
        staging = Path(root) / 'inputs'
        (staging / 'tenv').mkdir(parents=True)
        for name in ('manifest.json', *sorted(REQUIRED)):
            shutil.copyfile(tenv / name, staging / 'tenv' / name)
        names = {
            'manifest.json': 'representation.json',
            'parents.csv': 'representation_parents.csv',
            'KPGT.npy': 'KPGT.npy',
        }
        for original, output in names.items():
            shutil.copyfile(representations / original, staging / output)
        files = {
            name: sha256(staging / name)
            for name in ('tenv/manifest.json', *names.values())
        }
        (staging / 'manifest.json').write_text(
            json.dumps(
                {'schema': INPUT_SCHEMA, 'files': files}, indent=2, sort_keys=True
            )
            + '\n'
        )
        Inputs(staging)
        staging.rename(candidate)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ('tenv', 'representations', 'candidate'):
        parser.add_argument('--' + argument, type=Path, required=True)
    args = parser.parse_args()
    prepare(args.tenv, args.representations, args.candidate)
