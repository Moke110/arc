#!/usr/bin/env python
"""Publish/replay the per-category disease-DEG tables and within-dataset ranks.

ADR 0066 / ADR 0062. One directory per disease category under
``stage2/01_bases/output/disease-deg/``::

    <category>/catalog.json
    <category>/manifest.json
    <category>/<dataset>.deg.csv        # compact, one row per measured gene
    <category>/<dataset>.rank.csv       # frozen-axis rank frame (18,164 rows)

The driver is thin: it binds the pinned item-local DEG tables to the frozen axis
through :mod:`bases.disease_deg`, writes a candidate outside ``output/``
(``execution/checkpoints/``), verifies it, and only then replaces the accepted
category. ``verify`` replays the build from the pinned inputs and compares the
manifest and every output digest; the accepted release is never written.

The published catalog carries the ADR 0071 v2 fields (``n_deg``, ``deg2000``,
``tiering`` and the declared ``tier_sources`` ladder) once the driver runs: the
all-source tier order and the complete axis-resolved union are recorded by
:func:`bases.disease_deg.tier_ledger` over the sources this driver just measured,
so the record is never hand-written (ticket #206; products are republished by
#208).

Usage (from the repository root, with the shared ``arc`` env's ``python``)::

    python stage2/01_bases/execution/build_disease_deg.py publish --category nHF
    ... publish --all
    ... verify  --all
    ... ledger  --all        # print each category's tier ledger; writes no product

``--input-root``/``--sources`` exist so the hermetic tests can run the same
driver over a fixture tree and a fixture source declaration.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from importlib.machinery import SourceFileLoader
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve()
while not (REPO / '.git').exists():
    if REPO.parent == REPO:
        raise SystemExit('cannot locate the repository root')
    REPO = REPO.parent
EXECUTION = Path(__file__).resolve().parent
sys.path.insert(0, str(EXECUTION))
import public_api  # noqa: F401  (binds the single public implementation)
from bases import disease_deg, release  # noqa: E402

DEFAULT_OUTPUT = REPO / 'stage2/01_bases/output/disease-deg'
CHECKPOINTS = EXECUTION / 'checkpoints'
AXIS_RELATIVE = 'stage2/01_bases/output/gene-identity/genes.csv'
CROSSWALK_RELATIVE = 'stage2/01_bases/output/gene-identity/tables/gene_crosswalk.csv'


class Context:
    """The bound paths and source declaration one build runs against."""

    def __init__(self, output, input_root, sources):
        self.output = Path(output)
        self.input_root = Path(input_root)
        self.sources = sources

    def input(self, relative):
        return self.input_root / relative

    def relative(self, path):
        path = Path(path)
        try:
            return path.resolve().relative_to(self.input_root.resolve()).as_posix()
        except ValueError:
            raise release.ReleaseError(f'input outside the input root: {path}') from None


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


def item_input(context, dataset, relative):
    """Resolve an item-relative local path (no repo-relative literal in the source)."""
    path = context.input(dataset['source_item']) / relative
    if not path.is_file():
        raise release.ReleaseError(
            f"missing pinned item input {dataset['source_item']}/{relative}")
    return path


def source_inputs(context, dataset):
    """Pinned inputs of one dataset: tracked files by path, derived files by digest.

    Tracked inputs (item manifest, the item's tracked DEG summary, the axis and
    the crosswalk) are repo-relative. The item-local DEG/pseudobulk CSV artefacts
    are gitignored by their ADR 0030 item, so they are recorded as
    ``derived_inputs`` with an item-relative name and their sha256: the path is
    resolved at build time from ``source_item`` + the relative name, which keeps
    ``script/check_paths.sh`` green in a merged checkout that does not carry the
    per-worktree derived files, while every input is still digest-pinned.
    """
    inputs = {}
    item_manifest = context.input(Path(dataset['source_item']) / 'manifest.sha256')
    inputs[context.relative(item_manifest)] = sha256(item_manifest)
    summary = dataset.get('local_summary')
    if summary:
        path = item_input(context, dataset, summary)
        inputs[context.relative(path)] = sha256(path)
    derived = []
    for key in ('local_deg', 'cells_per_donor'):
        relative = dataset.get(key)
        if not relative:
            continue
        path = item_input(context, dataset, relative)
        derived.append({'relative': relative, 'sha256': sha256(path)})
    return inputs, derived


def read_category_sources(context, category, resolver):
    """Read, select, align and v2-rank every declared source of one category.

    One read path for both the product build and the tier ledger, so the
    declaration, the ladder and the published tables cannot drift. Each record
    carries the item-local table (canonical columns through the declaration's
    optional ``local_deg_columns`` alias), its axis-aligned rows, the alignment
    report and the v2 DEG-only table (``bases.deg_rank.dataset_deg_rank``, which
    carries the authoritative pass-any-contrast ``deg`` flag).
    """
    prepared = []
    for dataset in context.sources.datasets(category):
        local = disease_deg.read_local_deg(
            item_input(context, dataset, dataset['local_deg']),
            aliases=dataset.get('local_deg_columns'))
        local = disease_deg.contrast_selection(local, dataset['dataset_id'],
                                               dataset['contrasts'])
        aligned, report = disease_deg.align_source_table(local, resolver)
        v2_table, _v2_rank, v2_stats = disease_deg.dataset_deg_rank(
            aligned[['gene_index', 'log2fc', 'padj', 'contrast']])
        prepared.append({'dataset': dataset, 'local': local, 'aligned': aligned,
                         'alignment_report': report, 'v2_table': v2_table,
                         'v2_stats': v2_stats})
    return prepared


def source_deg_sets(prepared):
    """The measured v2 DEG set of every source, as :func:`tier_ledger` input.

    The set is the published-table convention of
    :func:`bases.disease_deg.source_deg_rank` applied to the v2 table (the
    authoritative ``deg`` flag), so it is the same set the republished product
    and its fusion will carry, resolved to the frozen axis by ``gene_index``.
    """
    measured = {}
    for item in prepared:
        derived = disease_deg.source_deg_rank(item['v2_table'])
        mask = np.asarray(derived['rank'], dtype=np.int64) <= int(derived['n_deg'])
        measured[item['dataset']['dataset_id']] = {
            'deg': np.flatnonzero(mask).astype(np.int64),
            'n_deg': int(derived['n_deg']),
            'informative': bool(derived['informative'])}
    return measured


def published_source_sets(context, category):
    """The measured v2 DEG set of every *published* source of one category.

    Fallback for a category whose sources this registry does not declare (the
    cardiotox classes, whose declaration lives in :mod:`bases.cardiotox` and whose
    accepted products are complete): read the digest-checked published tables and
    apply the same :func:`bases.disease_deg.source_deg_rank` convention.
    """
    root = disease_deg.category_root(context.output, category)
    catalog = json.loads((root / 'catalog.json').read_text())
    manifest = release.read_json(root / 'manifest.json')
    outputs = manifest.get('outputs') or {}
    datasets = catalog['datasets']
    measured = {}
    for record in datasets:
        name = record['deg']
        declared = outputs.get(name)
        expected = declared.get('sha256') if isinstance(declared, dict) else declared
        if not isinstance(expected, str) or len(expected) != 64:
            raise release.ReleaseError(
                f'{category}/{record["dataset_id"]}: {name} is not a declared output')
        path = root / name
        if sha256(path) != expected:
            raise release.ReleaseError(
                f'{category}/{record["dataset_id"]}: {name} digest differs from the '
                'category manifest')
        derived = disease_deg.source_deg_rank(
            pd.read_csv(path, float_precision='round_trip'))
        mask = np.asarray(derived['rank'], dtype=np.int64) <= int(derived['n_deg'])
        measured[record['dataset_id']] = {
            'deg': np.flatnonzero(mask).astype(np.int64),
            'n_deg': int(derived['n_deg']),
            'informative': bool(derived['informative'])}
    return datasets, measured


def category_ledger(context, category):
    """Compute one category's complete-source tier ledger from its declared sources.

    Registry-declared categories are measured from their item-local tables (the
    full declared roster, including products not yet published); a category whose
    declaration lives outside this registry (the cardiotox classes) is measured
    from its accepted products. Nothing is published: this is the audit/print
    path behind ``build_disease_deg.py ledger`` and the same computation the
    publish path writes into the catalog.
    """
    axis_path = context.input(AXIS_RELATIVE)
    axis = disease_deg.read_axis(axis_path)
    # Accept a registry label (nHF) or a product directory (nhf); a category whose
    # declaration lives outside this registry (the cardiotox classes) keeps its
    # product name.
    label = category
    for key, record in context.sources.CATEGORIES.items():
        if (record.get('directory') or key) == category:
            label = key
            break
    if any(record['category'] == label for record in context.sources.DATASETS):
        resolver = disease_deg.build_resolver(axis, disease_deg.read_crosswalk(
            context.input(CROSSWALK_RELATIVE)))
        prepared = read_category_sources(context, label, resolver)
        declared = context.sources.datasets(label)
        measured = source_deg_sets(prepared)
    else:
        declared, measured = published_source_sets(context, category)
    ledger = disease_deg.tier_ledger(label, declared, measured)
    return {'category': label,
            'directory': category,
            'tiering': ledger['tiering'],
            'tier_sources': ledger['tier_sources'],
            'n_deg': ledger['n_deg'],
            'deg2000': ledger['deg2000'],
            'sources': ledger['sources']}


def build_category(context, candidate, category):
    """Write one category's candidate directory; return its manifest."""
    candidate = Path(candidate)
    axis_path = context.input(AXIS_RELATIVE)
    crosswalk_path = context.input(CROSSWALK_RELATIVE)
    axis = disease_deg.read_axis(axis_path)
    crosswalk = disease_deg.read_crosswalk(crosswalk_path)
    resolver = disease_deg.build_resolver(axis, crosswalk)
    symbols = axis.symbol.to_numpy(dtype=object)
    category_record = context.sources.category_record(category)
    declared = context.sources.datasets(category)
    prepared = read_category_sources(context, category, resolver)
    records, outputs, inputs = [], {}, {}
    fusion_inputs = []
    for item in prepared:
        dataset = item['dataset']
        table, rank, stats, report = disease_deg.build_dataset_products(
            item['local'], resolver, symbols, dataset['dataset_id'])
        # The v2 source rank the tier ledger and the disease fusion consume; the
        # published table/rank above are still the current (pre-#208) product.
        derived = disease_deg.source_deg_rank(item['v2_table'])
        fusion_inputs.append({
            'dataset_id': dataset['dataset_id'], 'tier': dataset['tier'],
            'rank': derived['rank'], 'n_deg': derived['n_deg'],
            'informative': bool(derived['informative']),
            'informativeness_reason': derived['reason'],
        })
        disease_deg.validate_table(table)
        disease_deg.validate_rank(rank)
        deg_name = f"{dataset['dataset_id']}.deg.csv"
        rank_name = f"{dataset['dataset_id']}.rank.csv"
        # 17 significant digits round-trip every float64, so the published bytes
        # still encode the exact effect order the rank was derived from.
        table.to_csv(candidate / deg_name, index=False, float_format='%.17g')
        rank.to_csv(candidate / rank_name, index=False, float_format='%.17g')
        dataset_inputs, derived_inputs = source_inputs(context, dataset)
        inputs.update(dataset_inputs)
        method = dict(dataset['method'])
        method['alignment'] = 'symbol-first then Ensembl; unresolved rows reported'
        record = {
            'dataset_id': dataset['dataset_id'],
            'category': category,
            'source_item': dataset['source_item'],
            'arm_or_class': dataset['arm_or_class'],
            'species': dataset['species'],
            'model_system': dataset['model_system'],
            'evidence_type': dataset['evidence_type'],
            'tier': dataset['tier'],
            'genotype_status': dataset['genotype_status'],
            'time_window_basis': dataset['time_window_basis'],
            'exploratory': dataset['exploratory'],
            'shared_control_group': dataset['shared_control_group'],
            'n_case': dataset['n_case'],
            'n_control': dataset['n_control'],
            'contrasts': list(dataset['contrasts']),
            'n_measured': int(len(table)),
            'n_significant': int(table.significant.sum()),
            'n_axis_unmeasured': int(len(rank) - len(table)),
            # ADR 0068: a measured dataset whose effect scores are all tied
            # carries no ordering evidence. The table/rank are still published
            # unchanged; the fusion excludes the record (and a category with no
            # informative record entirely).
            'informative': bool(stats['informative']),
            'informativeness': {
                'rule': stats['informativeness_rule'],
                'n_distinct_effect': int(stats['n_distinct_effect']),
                'reason': stats['informativeness_reason'],
            },
            'deg': deg_name,
            'rank': rank_name,
            'method': method,
            'alignment': report,
            'inputs': dataset_inputs,
            'derived_inputs': derived_inputs,
        }
        # Optional declared evidence the canonical shape tolerates: the ADR 0066
        # isogenic-exception record, the multi-arm strata/pooling decision and the
        # ADR 0071 tier-2 non-genetic model basis / tier-3 whole-tissue basis.
        for optional in ('isogenic_exception', 'strata', 'non_genetic_model_basis',
                         'model_basis'):
            if optional in dataset:
                record[optional] = dataset[optional]
        records.append(record)
        for name in (deg_name, rank_name):
            path = candidate / name
            with path.open('rb') as handle:
                rows = sum(1 for _ in handle) - 1
            outputs[name] = {'sha256': sha256(path), 'rows': int(rows)}
    axis_sha = sha256(axis_path)
    inputs[context.relative(axis_path)] = axis_sha
    if crosswalk_path.is_file():
        inputs[context.relative(crosswalk_path)] = sha256(crosswalk_path)
    # ADR 0071 decision 2 (#206): the tier ladder, the cumulative axis-resolved
    # DEG union after every declared tier is computed from the declared
    # per-source tiers and the measured DEG sets, then published on the catalog.
    ledger = disease_deg.tier_ledger(category, declared, source_deg_sets(prepared))
    fusion = disease_deg.fuse_disease(fusion_inputs, axis_size=len(axis))
    if fusion['n_deg'] != ledger['n_deg']:
        raise release.ReleaseError(
            f'{category}: the tier ledger union ({ledger["n_deg"]}) differs from the '
            f'fusion union ({fusion["n_deg"]})')
    fusion_ledger = {row['dataset_id']: row for row in fusion['sources']}
    source_deg = {dataset_id: {'n_deg': int(row['n_deg']),
                               'n_deg_in_deg2000': int(row['n_deg_in_deg2000'])}
                  for dataset_id, row in fusion_ledger.items()}
    catalog = disease_deg.build_catalog(
        category, category_record['directory'], category_record['tier'],
        category_record['evidence_type'], axis, axis_sha, records,
        rejected=category_record['rejected'],
        isogenic_exception=category_record.get('isogenic_exception'),
        n_deg=fusion['n_deg'], deg2000=fusion['deg2000'],
        tiering=ledger['tiering'], tier_sources=ledger['tier_sources'],
        source_deg=source_deg)
    if 'source_review' in category_record:
        catalog['source_review'] = category_record['source_review']
    catalog_file = candidate / 'catalog.json'
    catalog_sha = disease_deg.write_catalog(catalog_file, catalog)
    outputs['catalog.json'] = {'sha256': sha256(catalog_file)}
    manifest = disease_deg.assemble_manifest(
        category,
        recipe=f"execution/build_disease_deg.py publish --category {category}",
        inputs=inputs, outputs=outputs,
        counts={'datasets': len(records),
                'n_measured': {record['dataset_id']: record['n_measured']
                               for record in records}},
        catalog_sha256=catalog_sha)
    release.write_json(candidate / 'manifest.json', manifest)
    return manifest


def _verify_directory(path):
    """Structural + digest verification of one category directory (candidate or accepted)."""
    path = Path(path)
    manifest = release.read_json(path / 'manifest.json')
    release.check_manifest(path, manifest)
    catalog = disease_deg.read_catalog(path / 'catalog.json')
    if catalog.get('schema') != disease_deg.CATALOG_SCHEMA:
        raise release.ReleaseError('catalog schema differs from the contract')
    if catalog['category'] != manifest['category']:
        raise release.ReleaseError('catalog/manifest category differ')
    if manifest['catalog_sha256'] != catalog['catalog_sha256']:
        raise release.ReleaseError('catalog_sha256 differs from the manifest')
    if not isinstance(catalog['datasets'], list) or not catalog['datasets']:
        raise release.ReleaseError('canonical catalog requires a non-empty dataset record list')
    expected = release.output_digests(manifest)
    for record in catalog['datasets']:
        dataset = record['dataset_id']
        if record['deg'] != f'{dataset}.deg.csv' or record['rank'] != f'{dataset}.rank.csv':
            raise release.ReleaseError(f'{dataset}: unexpected published filenames')
        table = disease_deg.validate_table(
            pd.read_csv(path / record['deg'], float_precision='round_trip'))
        rank = disease_deg.validate_rank(
            pd.read_csv(path / record['rank'], float_precision='round_trip'))
        if record['n_measured'] != len(table):
            raise release.ReleaseError(f'{dataset}: n_measured differs from the table')
        if record['n_axis_unmeasured'] != len(rank) - len(table):
            raise release.ReleaseError(f'{dataset}: n_axis_unmeasured differs from the rank')
        # A multi-contrast dataset collapses per gene by maximum effect, so the
        # published table carries the winning contrast label per gene: every
        # produced label must be a declared one and each declared label needs to
        # have won at least one gene (a declared-but-never-winning arm would mean
        # the declaration and the fit disagree).
        undeclared = sorted(set(table.contrast.astype(str)) - set(record['contrasts']))
        if undeclared:
            raise release.ReleaseError(f'{dataset}: undeclared contrasts {undeclared}')
        missing = sorted(set(record['contrasts']) - set(table.contrast.astype(str)))
        if missing:
            raise release.ReleaseError(f'{dataset}: declared contrasts {missing} won no gene')
        for name in (record['deg'], record['rank']):
            if sha256(path / name) != expected[name]:
                raise release.ReleaseError(f'{dataset}: {name} digest differs')
    if 'tiering' in catalog:
        # ADR 0073: a catalog that carries the complete tier ladder is validated
        # against the full v2 contract (the published products before #208 carry
        # only the v1 shape and skip this).
        try:
            disease_deg.validate_catalog(catalog)
        except ValueError as error:
            raise release.ReleaseError(f'{path.name}: {error}') from None
    return manifest, catalog


def verify_category(path):
    manifest, catalog = _verify_directory(path)
    return {'result': 'PASS', 'category': catalog['category'],
            'datasets': len(catalog['datasets']),
            'manifest_sha256': sha256(Path(path) / 'manifest.json')}


def verify_accepted(context, category):
    """Read the accepted release through the public reader after replacement."""
    directory = context.sources.category_record(category)['directory']
    path = disease_deg.category_root(context.output, directory)
    _manifest, catalog = _verify_directory(path)
    reader = disease_deg.DiseaseDeg(context.output)
    for record in catalog['datasets']:
        dataset = record['dataset_id']
        reader.table(category, dataset)
        reader.rank(category, dataset)
    return {'result': 'PASS', 'category': catalog['category'],
            'datasets': len(catalog['datasets']),
            'manifest_sha256': sha256(path / 'manifest.json')}


def publish_category(context, category):
    out = disease_deg.category_root(
        context.output, context.sources.category_record(category)['directory'])
    release.publish(out, lambda staging: build_category(context, staging, category),
                    verify_category, root=CHECKPOINTS)
    result = verify_accepted(context, category)
    accepted = release.read_json(out / 'manifest.json')
    result['outputs'] = {name: value['sha256']
                         for name, value in accepted['outputs'].items()}
    return result


def replay_category(context, category):
    out = disease_deg.category_root(
        context.output, context.sources.category_record(category)['directory'])
    result = release.replay(out, lambda staging: build_category(context, staging, category),
                            verify_category, root=CHECKPOINTS)
    result['category'] = category
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('publish', 'verify', 'ledger'))
    parser.add_argument('--category', action='append', default=[],
                        help='disease category (repeatable)')
    parser.add_argument('--all', action='store_true', help='every ticket-187 category')
    parser.add_argument('--root', default=None, help='output root (default: the step output)')
    parser.add_argument('--input-root', default=None,
                        help='root the pinned input paths resolve against (tests)')
    parser.add_argument('--sources', default=None,
                        help='source declaration module path (tests)')
    args = parser.parse_args(argv)
    sources = load_sources(args.sources)
    context = Context(Path(args.root) if args.root else DEFAULT_OUTPUT,
                      Path(args.input_root) if args.input_root else REPO, sources)
    categories = list(args.category)
    if args.all:
        if args.command == 'ledger':
            # Every category with a declaration or accepted products: the
            # registry's product directories plus the published categories whose
            # declaration lives elsewhere (the cardiotox classes). The crosswalk
            # directory is not a category.
            registry = {(record.get('directory') or label)
                        for label, record in sources.CATEGORIES.items()}
            categories = sorted(registry | set(
                disease_deg.DiseaseDeg(context.output).categories()))
        else:
            categories = sorted(sources.CATEGORIES)
    if not categories:
        parser.error('pass --category or --all')
    if args.command == 'ledger':
        results = [category_ledger(context, category) for category in categories]
        print(json.dumps({'command': 'ledger', 'results': results},
                         indent=2, sort_keys=True))
        return 0
    results = []
    for category in categories:
        if args.command == 'publish':
            results.append(publish_category(context, category))
        else:
            results.append(replay_category(context, category))
    print(json.dumps({'command': args.command, 'results': results},
                     indent=2, sort_keys=True))
    return 0


if __name__ == '__main__':
    sys.exit(main())
