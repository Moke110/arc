#!/usr/bin/env python
"""Publish/replay the per-disease level-1 fusion and ``DEG2000`` cut (ADR 0071).

The disease-half contract of the v2 gene axis: every disease category fuses its
source DEG sets at level 1 and keeps the first 2,000 genes of its DEG union.
ADR 0071 decisions 3 and 4; the pure definitions live in
``output/api/bases/disease_deg.py`` (``source_deg_rank``, ``level_one_rank``,
``fuse_disease``), so this driver only reads the accepted per-category products,
applies those definitions and publishes the artefacts.

Product layout (one directory per disease family)::

    output/disease-fusion/
        catalog.json                 # per disease: counts, ledger, exclusions
        manifest.json                # candidate-then-replace release manifest
        <directory>.deg2000.csv      # gene_index, disease_rank, disease_rank_sum

``catalog.json`` (``stage2-disease-fusion-v1``) is self-describing: it names the
frozen axis and its digest, the declared ``limit`` (2,000), the fusion
semantics, and one record per disease category::

    {
      "category": "nHF",                  # declared label, null if unpublished
      "directory": "nhf",                 # product directory under disease-deg/
      "published": true,                  # an accepted disease-deg product exists
      "tier": 1,                          # declared winning tier (1..3)
      "evidence_type": "cm_pseudobulk",
      "declared_source_count": 4,         # registry-declared sources
      "informative_source_count": 2,      # sources fused at level 1
      "n_deg": 1057,                      # DEG-union size (total DEG)
      "deg2000": 1057,                    # |DEG2000_d| == min(n_deg, 2000)
      "deg2000_file": "nhf.deg2000.csv",
      "sources": [                        # per-source ledger
        {"dataset_id": "gse270788__nicm", "tier": 1, "informative": true,
         "n_deg": 4, "n_deg_in_deg2000": 4, "evidence_type": "cm_pseudobulk",
         "flag_source": "collapsed_winning_contrast"}
      ],
      "excluded_sources": [               # published but non-informative
        {"dataset_id": "...", "tier": 1, "n_deg": 0, "reason": "..."}
      ],
      "not_published": [                  # declared but no accepted product yet
        {"dataset_id": "...", "reason": "declared by the source registry; no
         published disease-deg product in this checkout"}
      ],
      "tiers_present": [1],
      "reason": null                      # set only for an uninformative disease
    }

``<directory>.deg2000.csv`` is the per-gene ``DEG2000_d`` table: one row per
gene of the disease's DEG union whose disease rank is ``<= 2,000`` (all of them
when the union is smaller), ordered by ``disease_rank``. ``disease_rank_sum`` is
the level-1 ``Σ rank_norm`` that produced the rank. A disease with an empty union
still publishes a header-only file, so "no DEG" is representable rather than
absent.

The driver is deliberately roster-neutral: it fuses whatever sources the
``output/disease-deg/`` products declare and reports the rest of the source
registry (``disease_deg_sources.py``) as ``not_published``. The tier-order / stop
rule and the final roster are ticket #206's job, and the republish of the
per-category products (which fills the remaining sources) is ticket #8's; re-run
this driver after that republish. ``--sources``/``--input-root`` exist so the
hermetic tests can run the same driver over a fixture tree and declaration.

Usage (from the repository root, with the shared ``arc`` env's ``python``)::

    python stage2/01_bases/execution/build_disease_fusion.py publish
    python stage2/01_bases/execution/build_disease_fusion.py verify
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from importlib.machinery import SourceFileLoader
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve()
while not (REPO / '.git').exists():
    if REPO.parent == REPO:
        raise SystemExit('cannot locate the repository root')
    REPO = REPO.parent
EXECUTION = Path(__file__).resolve().parent
sys.path.insert(0, str(EXECUTION))
import public_api  # noqa: F401  (binds the single public implementation)
from bases import disease_deg, gene_axis, release  # noqa: E402

DEFAULT_OUTPUT = REPO / 'stage2/01_bases/output/disease-fusion'
DISEASE_DEG_RELATIVE = 'stage2/01_bases/output/disease-deg'
AXIS_RELATIVE = 'stage2/01_bases/output/gene-identity/genes.csv'
CHECKPOINTS = EXECUTION / 'checkpoints'

FUSION_SEMANTICS = {
    'deg': ('a gene is a DEG of a source when padj < 0.05 AND |log2FC| >= 2 in a '
            'declared contrast (ADR 0071 decision 1); each source ranks its DEGs '
            '1..n_DEG by descending effect with a gene_index tie-break and gives '
            'every other frozen-axis gene n_DEG + 1'),
    'level_one': ('disease_rank_sum = sum over the disease informative sources of '
                  'rank_norm = rank / n_DEG; the disease rank is the ascending '
                  'position of that sum with a gene_index tie-break'),
    'deg_union': ('the genes that are a DEG in at least one informative source; '
                  'its size is not padded and not constrained to 2,000'),
    'deg2000': (f'the DEG-union genes whose disease rank is <= '
                f'{disease_deg.DEG2000_LIMIT} (the whole union when it is '
                'smaller); a DEG union occupies the first ranks, so '
                f'|DEG2000| == min(|union|, {disease_deg.DEG2000_LIMIT})'),
    'non_informative': ('a source with n_DEG = 0 has no rank_norm (the ADR 0071 '
                        'denominator is zero) and a source whose measured effects '
                        'are all tied carries no ordering evidence; both are '
                        'excluded from the sum and recorded with their reason '
                        '(ADR 0068), so a NaN can never poison a disease sum'),
    'flag_derivation': ('a published v2 table carries the pass-any-contrast deg '
                        'flag (authoritative); a v1 table carries only the '
                        'collapsed winning-contrast rows, so its v2 flag is '
                        're-derived from those rows and the ledger records '
                        'flag_source = collapsed_winning_contrast'),
}


class Context:
    """The bound paths and source declaration one build runs against."""

    def __init__(self, output, input_root, sources, staging_root=None):
        self.output = Path(output)
        self.input_root = Path(input_root)
        self.sources = sources
        #: ADR 0062 keeps candidates outside ``output/``; the hermetic tests stage
        #: inside their own temporary root so the atomic rename stays on one
        #: filesystem.
        self.staging_root = Path(staging_root) if staging_root else CHECKPOINTS

    @property
    def disease_root(self):
        return self.input_root / DISEASE_DEG_RELATIVE

    def input(self, relative):
        return self.input_root / relative

    def relative(self, path):
        path = Path(path)
        try:
            return path.resolve().relative_to(self.input_root.resolve()).as_posix()
        except ValueError:
            raise release.ReleaseError(
                f'input outside the input root: {path}') from None


def sha256(path):
    return release.sha256(Path(path))


def load_sources(path=None):
    if path is None:
        import disease_deg_sources as module
        return module
    path = Path(path)
    loader = SourceFileLoader('disease_deg_sources_fixture', str(path))
    spec = importlib.util.spec_from_loader('disease_deg_sources_fixture', loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


# ---------------------------------------------------------------- input reads

def published_categories(context):
    """Every accepted disease-deg category, keyed by its product directory.

    The ``cardiotox/`` crosswalk directory carries no ``catalog.json`` and is not
    a disease category; every other directory with a catalog is one. The catalog
    is read and its ``catalog_sha256`` verified, and its digest must agree with
    the category manifest (the two-product contract of ADR 0062).
    """
    root = context.disease_root
    if not root.is_dir():
        raise release.ReleaseError(f'missing disease-deg product root: {root}')
    found = {}
    for directory in sorted(path for path in root.iterdir() if path.is_dir()):
        catalog_file = directory / 'catalog.json'
        if not catalog_file.is_file():
            continue  # the cardiotox/ crosswalk, not a disease product
        catalog = disease_deg.read_catalog(catalog_file)
        manifest_file = directory / 'manifest.json'
        if not manifest_file.is_file():
            raise release.ReleaseError(f'{directory.name}: no release manifest')
        manifest = release.read_json(manifest_file)
        if manifest.get('catalog_sha256') != catalog['catalog_sha256']:
            raise release.ReleaseError(
                f'{directory.name}: catalog digest differs from its manifest')
        found[directory.name] = {'catalog': catalog, 'manifest': manifest,
                                 'path': directory}
    if not found:
        raise release.ReleaseError('no published disease-deg category found')
    return found


def declared_sources(sources, directory):
    """The registry-declared dataset ids of one product directory (empty if none)."""
    for label, record in sources.CATEGORIES.items():
        if (record.get('directory') or label) == directory:
            return [dataset['dataset_id'] for dataset in sources.datasets(label)]
    return []


def registry_category(sources, directory):
    """The registry's category label and record for one product directory."""
    for label, record in sources.CATEGORIES.items():
        if (record.get('directory') or label) == directory:
            return label, record
    return None, None


def dataset_deg_table(entry, dataset):
    """Read one accepted compact DEG table, digest-checked against its manifest."""
    name = dataset.get('deg')
    declared = (entry['manifest'].get('outputs') or {}).get(name)
    digest = declared.get('sha256') if isinstance(declared, dict) else declared
    if not isinstance(digest, str) or len(digest) != 64:
        raise release.ReleaseError(
            f'{entry["path"].name}/{dataset["dataset_id"]}: deg table is not a '
            'declared output')
    path = entry['path'] / name
    if not path.is_file():
        raise release.ReleaseError(
            f'{entry["path"].name}/{dataset["dataset_id"]}: missing {name}')
    if sha256(path) != digest:
        raise release.ReleaseError(
            f'{entry["path"].name}/{dataset["dataset_id"]}: deg digest differs '
            'from the category manifest')
    return pd.read_csv(path, float_precision='round_trip')


def source_records_for_fusion(entry):
    """One fusion input per published dataset of a category (v2 DEG rank derived).

    The *published* ``informative`` flag is the declared ADR 0068 verdict and is
    used when the record carries it; the flag derived from the table is the
    fallback for a catalog that predates the field. A source whose published DEG
    count is zero is excluded by :func:`bases.disease_deg.fuse_disease` either way.
    """
    prepared = []
    for dataset in entry['catalog']['datasets']:
        derived = disease_deg.source_deg_rank(dataset_deg_table(entry, dataset))
        recorded = dataset.get('informativeness') or {}
        reason = recorded.get('reason') if isinstance(recorded, dict) else None
        prepared.append({
            'dataset_id': dataset['dataset_id'],
            'tier': dataset.get('tier'),
            'evidence_type': dataset.get('evidence_type'),
            'rank': derived['rank'],
            'n_deg': derived['n_deg'],
            'informative': bool(dataset.get('informative', derived['informative'])),
            'informativeness_reason': reason or derived['reason'],
            'flag_source': derived['flag_source'],
        })
    return prepared


# ---------------------------------------------------------------- computation

def disease_record(directory, entry, sources, axis_size):
    """The catalog record and fusion of one category.

    A declared category with no accepted disease-deg product stays in the catalog
    as ``published: false`` (reporting the roster gap) instead of being dropped.
    """
    declared = declared_sources(sources, directory)
    label, declared_record = registry_category(sources, directory)
    if entry is None:
        return {
            'category': label,
            'directory': directory,
            'published': False,
            'tier': None if declared_record is None else declared_record.get('tier'),
            'evidence_type': (None if declared_record is None
                              else declared_record.get('evidence_type')),
            'declared_source_count': len(declared),
            'informative_source_count': 0,
            'n_deg': None,
            'deg2000': None,
            'deg2000_file': None,
            'sources': [],
            'excluded_sources': [],
            'not_published': [
                {'dataset_id': dataset_id,
                 'reason': ('declared by the source registry; no accepted '
                            'disease-deg product in this checkout')}
                for dataset_id in declared],
            'tiers_present': [],
            'reason': ('declared but not published: the #8 republish must publish '
                       'its disease-deg products, then this driver is re-run'),
        }, None
    prepared = source_records_for_fusion(entry)
    fusion = disease_deg.fuse_disease(prepared, axis_size=axis_size)
    published_ids = [dataset['dataset_id'] for dataset in entry['catalog']['datasets']]
    not_published = [
        {'dataset_id': dataset_id,
         'reason': ('declared by the source registry; no accepted disease-deg '
                    'product in this checkout')}
        for dataset_id in declared if dataset_id not in set(published_ids)]
    evidence = {dataset['dataset_id']: dataset.get('evidence_type')
                for dataset in entry['catalog']['datasets']}
    flags = {source['dataset_id']: source['flag_source'] for source in prepared}
    ledger = []
    for source in fusion['sources']:
        row = dict(source)
        row['evidence_type'] = evidence[row['dataset_id']]
        row['flag_source'] = flags[row['dataset_id']]
        ledger.append(row)
    tiers_present = sorted({int(row['tier']) for row in ledger
                            if row['tier'] is not None})
    record = {
        'category': entry['catalog'].get('category'),
        'directory': directory,
        'published': True,
        'tier': entry['catalog'].get('tier'),
        'evidence_type': entry['catalog'].get('evidence_type'),
        'declared_source_count': len(set(declared) | set(published_ids)),
        'informative_source_count': fusion['n_sources_informative'],
        'n_deg': fusion['n_deg'],
        'deg2000': fusion['deg2000'],
        'deg2000_file': f'{directory}.deg2000.csv',
        'sources': ledger,
        'excluded_sources': fusion['excluded_sources'],
        'not_published': not_published,
        'tiers_present': tiers_present,
        'reason': fusion['reason'],
    }
    return record, fusion


def build(context, staging):
    """Write the candidate fusion catalog, manifest and per-disease DEG2000 tables."""
    staging = Path(staging)
    axis_path = context.input(AXIS_RELATIVE)
    axis_size = len(disease_deg.read_axis(axis_path))
    published = published_categories(context)
    declared_dirs = {(record.get('directory') or label)
                     for label, record in context.sources.CATEGORIES.items()}
    # Deterministic disease order: the axis roster order first, then any other
    # declared/published product directory sorted.
    directories = [name for name in gene_axis.CATEGORIES
                   if name in published or name in declared_dirs]
    directories += sorted((set(published) | declared_dirs) - set(directories))
    inputs, outputs, diseases = {}, {}, []
    for directory in directories:
        entry = published.get(directory)
        record, fusion = disease_record(directory, entry, context.sources, axis_size)
        diseases.append(record)
        if entry is None:
            continue
        rows = disease_deg.deg2000_rows(fusion)
        name = record['deg2000_file']
        rows.to_csv(staging / name, index=False, lineterminator='\n',
                    float_format='%.17g')
        outputs[name] = {'sha256': sha256(staging / name), 'rows': int(len(rows))}
        # Pin what the fusion read: the accepted category catalog, its manifest
        # and every compact DEG table it fused (repo-relative, digest-checked).
        for role, path in (('catalog', entry['path'] / 'catalog.json'),
                           ('release_manifest', entry['path'] / 'manifest.json')):
            inputs[context.relative(path)] = {
                'role': f'{directory}:{role}', 'sha256': sha256(path),
                'bytes': path.stat().st_size}
        for dataset in entry['catalog']['datasets']:
            path = entry['path'] / dataset['deg']
            inputs[context.relative(path)] = {
                'role': f'{directory}:{dataset["dataset_id"]}:deg',
                'sha256': sha256(path), 'bytes': path.stat().st_size}
    inputs[context.relative(axis_path)] = {
        'role': 'frozen_axis', 'sha256': sha256(axis_path),
        'bytes': axis_path.stat().st_size}
    catalog = {
        'schema': disease_deg.FUSION_SCHEMA,
        'generated_by': 'stage2/01_bases/execution/build_disease_fusion.py',
        'axis': {'path': AXIS_RELATIVE, 'sha256': sha256(axis_path),
                 'size': int(axis_size)},
        'rank_contract': disease_deg.RANK_SCHEMA,
        'rank_semantics': disease_deg.RANK_SEMANTICS,
        'limit': disease_deg.DEG2000_LIMIT,
        'semantics': dict(FUSION_SEMANTICS),
        'diseases': diseases,
    }
    disease_deg.validate_fusion_catalog(catalog)
    catalog_file = staging / 'catalog.json'
    catalog_sha = disease_deg.write_catalog(catalog_file, catalog)
    outputs['catalog.json'] = {'sha256': sha256(catalog_file)}
    manifest = {
        'schema': disease_deg.FUSION_RELEASE_SCHEMA,
        'recipe': 'execution/build_disease_fusion.py publish',
        'axis': {'path': AXIS_RELATIVE},
        'limit': disease_deg.DEG2000_LIMIT,
        'inputs': dict(sorted(inputs.items())),
        'outputs': dict(sorted(outputs.items())),
        'counts': {
            'diseases': len(diseases),
            'published': sum(1 for record in diseases if record['published']),
            'deg2000': sum(record['deg2000'] or 0 for record in diseases),
        },
        'catalog_sha256': catalog_sha,
    }
    release.check_relative_paths(manifest)
    release.write_json(staging / 'manifest.json', manifest)
    return manifest


def _verify_directory(path):
    """Structural + digest verification of a candidate or accepted fusion release."""
    path = Path(path)
    manifest = release.read_json(path / 'manifest.json')
    if manifest.get('schema') != disease_deg.FUSION_RELEASE_SCHEMA:
        raise release.ReleaseError('fusion manifest schema differs from the contract')
    release.check_manifest(path, manifest)
    catalog = disease_deg.read_catalog(path / 'catalog.json')
    disease_deg.validate_fusion_catalog(catalog)
    if manifest.get('catalog_sha256') != catalog['catalog_sha256']:
        raise release.ReleaseError('catalog_sha256 differs from the fusion manifest')
    if manifest.get('limit') != catalog['limit']:
        raise release.ReleaseError('manifest limit differs from the catalog')
    expected = release.output_digests(manifest)
    for record in catalog['diseases']:
        if not record['published']:
            continue
        name = record['deg2000_file']
        frame = pd.read_csv(path / name,
                            dtype={'gene_index': 'int64', 'disease_rank': 'int64'})
        if list(frame.columns) != list(disease_deg.DEG2000_COLUMNS):
            raise release.ReleaseError(f'{name}: unexpected DEG2000 columns')
        if len(frame) != record['deg2000']:
            raise release.ReleaseError(
                f'{name}: rows differ from the declared DEG2000 size')
        if record['deg2000']:
            if frame.disease_rank.tolist() != list(range(1, record['deg2000'] + 1)):
                raise release.ReleaseError(
                    f'{name}: rows are not in disease-rank order 1..deg2000')
            if not frame.gene_index.is_unique:
                raise release.ReleaseError(f'{name}: repeats a gene_index')
        elif len(frame):
            raise release.ReleaseError(f'{name}: an empty union publishes no row')
        if sha256(path / name) != expected[name]:
            raise release.ReleaseError(f'{name}: digest differs from the manifest')
    return manifest, catalog


def verify_fusion(path):
    manifest, catalog = _verify_directory(path)
    return {'result': 'PASS', 'diseases': len(catalog['diseases']),
            'published': sum(1 for record in catalog['diseases']
                             if record['published']),
            'manifest_sha256': sha256(Path(path) / 'manifest.json')}


def publish(context):
    release.publish(context.output, lambda staging: build(context, staging),
                    verify_fusion, root=context.staging_root)
    reader = disease_deg.DiseaseFusion(context.output)
    catalog = reader.catalog()
    accepted = reader.manifest()
    for record in catalog['diseases']:
        if record['published']:
            reader.deg2000(record['directory'])
    return {'result': 'PASS',
            'catalog_sha256': catalog['catalog_sha256'],
            'manifest_sha256': sha256(context.output / 'manifest.json'),
            'diseases': len(catalog['diseases']),
            'published': sum(1 for record in catalog['diseases']
                             if record['published']),
            'deg2000': accepted['counts']['deg2000']}


def replay(context):
    result = release.replay(context.output,
                            lambda staging: build(context, staging),
                            verify_fusion, root=context.staging_root)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('publish', 'verify'))
    parser.add_argument('--root', default=None,
                        help='fusion output root (default: the step output)')
    parser.add_argument('--input-root', default=None,
                        help='root the pinned input paths resolve against (tests)')
    parser.add_argument('--sources', default=None,
                        help='source declaration module path (tests)')
    args = parser.parse_args(argv)
    sources = load_sources(args.sources)
    context = Context(Path(args.root) if args.root else DEFAULT_OUTPUT,
                      Path(args.input_root) if args.input_root else REPO, sources)
    result = publish(context) if args.command == 'publish' else replay(context)
    print(json.dumps({'command': args.command, 'results': result},
                     indent=2, sort_keys=True))
    return 0


if __name__ == '__main__':
    sys.exit(main())
