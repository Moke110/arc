#!/usr/bin/env python
"""Build/replay the bases-only Tahoe shared-axis raw-count release."""
import argparse
import json
from pathlib import Path

import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases.tahoe import build_tahoe, verify_tahoe
from bases.identity import AXIS_SHA256
import hashlib


def repository(start):
    for candidate in (start, *start.parents):
        if (candidate / '.git').exists():
            return candidate
    raise ValueError('repository root not found')


def main():
    root = repository(Path(__file__).resolve().parent)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('build', 'verify'))
    parser.add_argument('--source', type=Path, default=root / 'pert_data/Tahoe-100M/source')
    parser.add_argument('--source-manifest', type=Path, default=root / 'pert_data/Tahoe-100M/manifest.sha256')
    parser.add_argument('--axis', type=Path, default=root / 'stage2/01_bases/output/gene-identity/tables/genes.csv')
    parser.add_argument('--crosswalk', type=Path, default=root / 'stage2/01_bases/output/gene-identity/tables/gene_crosswalk.csv')
    parser.add_argument('--output', type=Path, default=root / 'stage2/01_bases/output/tahoe-raw-counts')
    parser.add_argument('--batch-size', type=int, default=1024)
    args = parser.parse_args()
    if args.batch_size <= 0:
        parser.error('--batch-size must be positive')
    # The shared resolver rejects ambiguous source identities and seals the
    # source-specific Ensembl/token evidence rather than guessing from symbols.
    if hashlib.sha256(args.axis.read_bytes()).hexdigest() != AXIS_SHA256:
        raise ValueError('frozen gene axis digest mismatch')
    from bases.raw import tahoe_axis_mapping
    mapping = tahoe_axis_mapping(args.axis, args.crosswalk, args.source / 'metadata/gene_metadata.parquet')
    method = build_tahoe if args.action == 'build' else verify_tahoe
    result = method(args.source, mapping, args.source_manifest, args.output, batch_size=args.batch_size)
    print(json.dumps({key: value for key, value in result.items() if key != 'inputs'}, sort_keys=True), flush=True)


if __name__ == '__main__':
    main()
