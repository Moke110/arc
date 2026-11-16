#!/usr/bin/env python3
"""Build/replay unified TOX only as a staged candidate (ADR 0072)."""
import argparse
from pathlib import Path
import public_api  # noqa: F401
from bases import tox


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('publish', 'verify'))
    parser.add_argument('--repo', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--workdir', type=Path)
    args = parser.parse_args()
    repo = args.repo or next(p for p in Path(__file__).resolve().parents if (p / '.git').exists())
    output = args.output or repo / 'stage2/01_bases/execution/checkpoints/unified-tox'
    if args.command == 'publish':
        catalog = tox.build_candidate(repo, output, args.workdir)
        print(f"TOX candidate publish PASS: {catalog['n_deg']} DEGs, {catalog['deg2000']} DEG2000")
    else:
        tox.verify_candidate(repo, output, args.workdir)
        print('TOX candidate verify PASS')


if __name__ == '__main__':
    main()
