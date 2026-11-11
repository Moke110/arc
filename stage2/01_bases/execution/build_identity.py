"""Publish/replay the small indexed identity release; no expression data is loaded."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from bases.identity import GeneIdentity, native_mapping_rows


def repo_root():
    return next(parent for parent in Path(__file__).resolve().parents if (parent / '.git').exists())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def publish(native, annotation, output, shared, verify=False):
    axis = native / 'tables/genes.csv'
    crosswalk = native / 'tables/gene_crosswalk.csv'
    mappings = native_mapping_rows(shared)
    identity = GeneIdentity.load(axis, crosswalk, annotation, native_mappings=mappings)
    all_indices = list(range(18164))
    if identity.indices(identity.labels(all_indices, 'hgnc_id')) != all_indices:
        raise ValueError('full-axis round trip failed')
    fam = identity.indices(['HGNC:27323'])[0]
    for source, identifier in [('tahoe','ENSG00000289731'), ('tenv','ENSG00000182230')]:
        if identity.indices([identifier], source=source) != [fam]:
            raise ValueError(f'FAM153B source identity differs: {source}')
    manifest = {
        'schema_version': 1, 'axis_size': 18164, 'gene_index': '0-based frozen row order',
        'inputs': {name: {'sha256': digest(path)} for name, path in
                   [('genes.csv',axis), ('gene_crosswalk.csv',crosswalk), ('hgnc_complete_set.txt',annotation), *[(f'{source}_mapping.csv', shared / f'{source}_mapping.csv') for source in ('tahoe','tenv')]]},
        'source_evidence_rows': len(identity.evidence),
        'fam153b_gene_index': fam,
        'checks': ['full-axis HGNC roundtrip', 'source-specific FAM153B Tahoe/tenv'],
    }
    import io
    stream = io.StringIO(newline='')
    fields = ['gene_index', *[key for key in identity.axis[0] if key != 'gene_index']]
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n')
    writer.writeheader()
    writer.writerows(identity.axis)
    payload = stream.getvalue()
    manifest['genes_sha256'] = hashlib.sha256(payload.encode()).hexdigest()
    serialized = json.dumps(manifest, indent=2, sort_keys=True) + '\n'
    if verify:
        if (output / 'genes.csv').read_text() != payload or (output / 'manifest.json').read_text() != serialized:
            raise ValueError('indexed identity release differs from pinned inputs')
    else:
        output.mkdir(parents=True, exist_ok=True)
        (output / 'genes.csv').write_text(payload)
        (output / 'manifest.json').write_text(serialized)
    return manifest


if __name__ == '__main__':
    root = repo_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', type=Path, default=root / 'stage2/02_dataset_corpus/output/corpus/v2')
    parser.add_argument('--annotation', type=Path, default=root / 'pert_data/Annotation/source/hgnc_complete_set.txt')
    parser.add_argument('--output', type=Path, default=root / 'stage2/01_bases/output/identity-v1')
    parser.add_argument('--shared', type=Path, default=root / 'stage2/02_dataset_corpus/output/corpus/shared-cm-v1')
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    print(json.dumps(publish(args.native, args.annotation, args.output, args.shared, args.verify), indent=2))
