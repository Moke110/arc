#!/usr/bin/env python3
"""Build/replay a pre-analysis TOX admission candidate, without reading expression."""
import argparse
from pathlib import Path
import tempfile
import public_api
from bases.tox_admission import build_ledger, ToxAdmission
from bases.release import ReleaseError, sha256

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('publish', 'verify'))
    parser.add_argument('--repo', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    repo = args.repo or next(p for p in Path(__file__).resolve().parents if (p/'.git').exists())
    output = args.output or repo/'stage2/01_bases/execution/checkpoints/tox-admission'
    if args.command == 'publish':
        build_ledger(repo, output)
        ToxAdmission(output, repo)
    else:
        ToxAdmission(output, repo)
        with tempfile.TemporaryDirectory(prefix='tox-admission-replay-', dir=output.parent) as temporary:
            replay = Path(temporary)
            build_ledger(repo, replay)
            for name in ('admission.json', 'manifest.json'):
                if sha256(output/name) != sha256(replay/name):
                    raise ReleaseError('TOX admission replay differs: '+name)
    print('TOX admission '+args.command+' PASS: '+str(output))

if __name__ == '__main__':
    main()
