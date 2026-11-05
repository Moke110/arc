#!/usr/bin/env python3
"""Build or verify the complete first-wave phenotype release, entirely offline."""
import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import re
import sys
import tempfile

import corpus
import heliste_adapter
import sirenko_adapter
import viability_adapters

HERE = Path(__file__).resolve().parent
STEP = HERE.parent


def frozen_bundle():
    root = corpus.root_path()
    registry = json.loads((root/'pheno_data/registry.json').read_text())
    expected = {'sirenko_2013', 'heliste_2021', 'linn_2024', 'liu_2021', 'duan_2013'}
    entries = {row['dataset_id']: row for row in registry['items']}
    if not expected <= entries.keys():
        raise ValueError('First-wave source inventory incomplete')
    bundles = []
    if entries['sirenko_2013']['status'] == 'active':
        bundles.append(sirenko_adapter.extract(root)[0])
    if entries['heliste_2021']['status'] == 'active':
        bundle, audit = heliste_adapter.build_bundle(root)
        bundle['source_audits'] = {'heliste_2021': audit}
        bundles.append(bundle)
    viability = [key for key in ('linn_2024','liu_2021','duan_2013') if entries[key]['status'] == 'active']
    if viability:
        bundles.append(viability_adapters.build_sources(viability))
    combined = {'schema_version': 1, 'sources': [], 'records': [], 'source_audits': {}}
    for bundle in bundles:
        combined['sources'].extend(bundle['sources'])
        combined['records'].extend(bundle['records'])
        combined['source_audits'].update(bundle.get('source_audits', {}))
    for key in sorted(expected):
        entry = entries[key]
        if entry['status'] != 'active':
            reason = entry.get('reason')
            if not reason:
                raise ValueError('Non-active configured source requires recorded acquisition reason: '+key)
            combined['sources'].append({'source': key, 'item': 'pheno_data/'+entry['item'],
                'status': 'unavailable', 'test_only': False, 'assets': [], 'reason': reason})
            combined['source_audits'][key] = {'acquisition_status': entry['status'], 'reason': reason,
                                            'extraction_status': 'not_run'}
    combined['sources'].sort(key=lambda row: row['source'])
    for source in combined['sources']:
        if source['status'] != 'unavailable':
            item = root/source['item']
            provenance = item/'author_log.json'
            source['provenance'] = json.loads(provenance.read_text())
            source['provenance_sha256'] = hashlib.sha256(provenance.read_bytes()).hexdigest()
            source['source_manifest_sha256'] = hashlib.sha256((item/'manifest.sha256').read_bytes()).hexdigest()
    combined['records'].sort(key=lambda row: (row['source'],row['comparison_id']))
    combined['adapter_versions'] = {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        for name in ('release.py','sirenko_adapter.py','heliste_adapter.py','viability_adapters.py')}
    combined['adapter_software'] = {name: version(name) for name in ('xlrd','openpyxl')}
    return combined


def release_path(identity):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', identity):
        raise ValueError('Release ID must be a single versioned directory name')
    return STEP/'output/corpus'/identity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    bundle_parser = sub.add_parser('bundle', help='Write the frozen source-adapter input without publishing')
    bundle_parser.add_argument('--output', type=Path, required=True)
    for command in ('build','verify'):
        command_parser = sub.add_parser(command)
        command_parser.add_argument('--release-id', required=True)
    args = parser.parse_args()
    try:
        if args.command == 'bundle':
            bundle = frozen_bundle()
            args.output.parent.mkdir(parents=True, exist_ok=True)
            corpus.write_json(args.output,bundle)
            result = {'source_count': len(bundle['sources']), 'comparisons': len(bundle['records'])}
        elif args.command == 'verify':
            result = corpus.verify(release_path(args.release_id))
        else:
            bundle = frozen_bundle()
            with tempfile.TemporaryDirectory(prefix='stage3-frozen-input-') as directory:
                path = Path(directory)/'input.json'
                corpus.write_json(path,bundle)
                result = corpus.build(path,release_path(args.release_id))
        print(json.dumps(result,sort_keys=True))
    except (ValueError,KeyError,TypeError,OSError) as error:
        print('FAIL: '+str(error),file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
