#!/usr/bin/env python3
"""Publish or verify the accepted cold-drug assignment and fold/seed index tables.

The accepted tenv cold-drug split is re-homed once as
``output/fold_n_seed/`` (``cohorts.csv``, ``folds.csv``, ``seeds.csv`` and a
digest manifest); the published product no longer reads the HVG product.  The
build copies the assignment byte-identically, refuses a source that is not the
accepted digest, writes the two index tables deterministically, and then runs
the offline verification itself so a failed candidate is never reported as
accepted.

    python stage2/01_bases/execution/build_fold_n_seed.py \\
        --cohorts stage2/01_bases/output/tenv-hvg/cohorts.csv
    python stage2/01_bases/execution/build_fold_n_seed.py --verify

``--expected-cohorts-sha256`` exists for fixtures and candidate inspection; its
default is the accepted assignment digest.
"""
import argparse
import json
from pathlib import Path
import sys

import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases.fold_seed import (ACCEPTED_COHORTS_SHA256, FoldSeedError,
                             build_fold_seed, verify_fold_seed)

STEP = Path(__file__).resolve().parents[1]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cohorts', type=Path,
                        help='source assignment to re-home (required to build)')
    parser.add_argument('--output', type=Path, default=STEP / 'output/fold_n_seed')
    parser.add_argument('--expected-cohorts-sha256', default=ACCEPTED_COHORTS_SHA256,
                        help='accepted assignment digest (default: the accepted one)')
    parser.add_argument('--verify', action='store_true',
                        help='verify the published product instead of building it')
    parser.add_argument('--record', type=Path, help='write the JSON result here')
    args = parser.parse_args(argv)
    code = 0
    if args.verify:
        result = verify_fold_seed(args.output,
                                  expected_cohorts_sha256=args.expected_cohorts_sha256)
        code = 0 if result['result'] == 'PASS' else 1
    else:
        if args.cohorts is None:
            parser.error('--cohorts is required unless --verify is given')
        try:
            manifest = build_fold_seed(
                args.cohorts, args.output,
                accepted_sha256=args.expected_cohorts_sha256)
        except (FoldSeedError, OSError) as exc:
            result = {'result': 'FAIL', 'error': str(exc),
                      'output': str(args.output)}
            code = 1
        else:
            verification = verify_fold_seed(
                args.output, expected_cohorts_sha256=args.expected_cohorts_sha256)
            result = {'result': verification['result'], 'manifest': manifest,
                      'verification': verification}
            code = 0 if verification['result'] == 'PASS' else 1
    text = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.record:
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(text)
    print(text, end='')
    return code


if __name__ == '__main__':
    sys.exit(main())
