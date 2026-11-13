#!/usr/bin/env python3
"""Build or replay the released Stage 2 gene axes ``G_2000`` / ``G_5000`` (ADR 0066).

The products are the two-level rank-sum fusion of the six cardiomyocyte-derived
disease-DEG categories plus the published tenv-HVG ranking: the master candidate
table, the two released axes, their manifest, the per-category composition table
and the rank-sum histogram evidence panel. No corpus, view, model or truth is
built; the disease-DEG products and the tenv-HVG ranking are read-only inputs.

Every scientific definition lives in ``output/api/bases/gene_axis.py`` (the single
public implementation, ADR 0062); this driver only reads the pinned inputs,
applies those definitions and publishes the products.

    python stage2/01_bases/execution/build_gene_axis.py
    python stage2/01_bases/execution/build_gene_axis.py --verify

Publication is candidate-then-replace (ADR 0062): a complete candidate is built
outside ``output/``, validated, and only then swapped in; a failed build leaves
the accepted release usable. ``--verify`` rebuilds from the pinned inputs into a
scratch candidate and compares the manifest and every output digest.
"""
import argparse
import hashlib
import json
import os
import platform
import shutil
import sys
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import pandas as pd

import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases import release
from bases.disease_deg import DiseaseDeg
from bases.gene_axis import (
    AXIS_SIZES, CATEGORIES, CATEGORY_LABELS, CANDIDATE_COLUMNS, FILL_RULE,
    N_FUSION_SOURCES, SCHEMA, axis_composition, axis_fill, axis_intersection,
    candidate_columns, candidate_table, category_rank, detect_inflection, fusion_rank_sum,
    inflection_stability, pn_distribution, rank_sum_histogram, winning_tier,
)
from bases.deg_rank import SCHEMA as RANK_SCHEMA, normalized_rank

#: The two accepted category product families that can contribute a level-1 source:
#: the #187/#188 within-dataset release and the #189 cardiotoxicity class release.
CATEGORY_MANIFEST_SCHEMAS = ('stage2-disease-deg-release-v1', 'stage2-cardiotox-class-v1')
CANDIDATE_NAME = 'candidate_rank_sum.csv'
COMPOSITION_NAME = 'axis_composition.csv'
FIGURE_NAME = 'gene_axis_rank_sum_histogram'
HISTOGRAM_BINS = 70
#: Inflection is diagnostic only; the declared decision range starts at N=100.
INFLECTION_MIN_N = 100
#: Values of rank_sum live in (0, 7]; the histogram covers the full possible range.
RANK_SUM_MAX = float(N_FUSION_SOURCES)
COMPOSITION_NOTE = ('pN is an ordinal cross-source replication count (0..7), '
                    'not a continuous quality score')

FIGURE_COLUMNS = ('record', 'panel', 'panel_title', 'title', 'axis', 'bin_left',
                  'bin_right', 'bin_center', 'count', 'share', 'interval_low',
                  'interval_high', 'median', 'inflection', 'inflection_n', 'note')
FIGURE_TITLE = ('Released Stage 2 gene axes · rank-sum distribution over all '
                'frozen-axis candidates with the G_2000 and G_5000 intervals')
COMPOSITION_COLUMNS = ('axis', 'n_axis', 'record', 'group', 'n_genes', 'share', 'note')


@dataclass(frozen=True)
class Config:
    """Every path the build reads or writes; all are repo-relative by default."""

    root: Path
    genes: Path
    tenv_hvg: Path
    disease_root: Path
    output: Path
    figure_data: Path
    axis_sizes: tuple = AXIS_SIZES


def repo_root():
    return next(parent for parent in Path(__file__).resolve().parents
                if (parent / '.git').exists())


def default_config(root=None, output=None, figure_data=None, axis_sizes=AXIS_SIZES):
    root = Path(root) if root else repo_root()
    step = root / 'stage2/01_bases'
    return Config(
        root=root,
        genes=step / 'output/gene-identity/genes.csv',
        tenv_hvg=step / 'output/tenv-hvg/descriptive.csv',
        disease_root=step / 'output/disease-deg',
        output=Path(output) if output else step / 'output/gene-axis',
        figure_data=Path(figure_data) if figure_data
        else step / 'figure/data' / f'{FIGURE_NAME}.csv',
        axis_sizes=tuple(int(size) for size in axis_sizes),
    )


def repo_relative(path, root):
    """Lexical repo-relative path; an output outside the repository fails closed."""
    lexical = Path(os.path.abspath(path))
    try:
        return lexical.relative_to(Path(os.path.abspath(root))).as_posix()
    except ValueError as error:
        raise ValueError(f'product path is outside the repository root: {path}') from error


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write_frame(frame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(frame.to_csv(index=False, lineterminator='\n',
                                 float_format='%.17g'))
    return path


def write_json(payload, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n')
    return path


# ---------------------------------------------------------------- input reads

def read_axis(path):
    """Read the frozen axis and require the contiguous ``0..n-1`` gene_index order."""
    axis = pd.read_csv(path, dtype={'symbol': str, 'hgnc_id': str,
                                    'ensembl_id': str, 'match_type': str})
    required = {'gene_index', 'symbol'}
    if not required <= set(axis.columns):
        raise ValueError(f'{path}: axis table requires {sorted(required)}')
    gene_index = np.asarray(axis.gene_index, dtype=np.int64)
    if len(gene_index) < 2 or not np.array_equal(
            gene_index, np.arange(len(gene_index), dtype=np.int64)):
        raise ValueError(f'{path}: axis gene_index is not the contiguous 0..n-1 order')
    if axis.symbol.isna().any() or axis.symbol.duplicated().any():
        raise ValueError(f'{path}: axis symbols must be unique and non-null')
    return axis


def _published_digest(path, manifest_path, getter):
    """Compare one input against its published manifest digest when one is present."""
    manifest_path = Path(manifest_path)
    if not manifest_path.is_file():
        return None
    manifest = json.loads(manifest_path.read_text())
    declared = getter(manifest)
    if not isinstance(declared, str) or len(declared) != 64:
        raise ValueError(f'{manifest_path}: no usable digest for {path}')
    actual = sha256(path)
    if actual != declared:
        raise ValueError(f'{path}: digest differs from its published manifest')
    return actual


def read_tenv(path, axis_size):
    """The published tenv-HVG dispersion rank: the seventh fusion source."""
    frame = pd.read_csv(path, dtype={'gene_index': np.int64})
    required = {'gene_index', 'rank'}
    if not required <= set(frame.columns):
        raise ValueError(f'{path}: tenv-HVG ranking requires {sorted(required)}')
    gene_index = np.asarray(frame.gene_index, dtype=np.int64)
    if len(frame) != axis_size or not np.array_equal(
            gene_index, np.arange(axis_size, dtype=np.int64)):
        raise ValueError(f'{path}: tenv-HVG ranking must cover the frozen axis')
    rank = np.asarray(frame['rank'], dtype=np.int64)
    if not np.array_equal(np.sort(rank), np.arange(1, axis_size + 1, dtype=np.int64)):
        raise ValueError(f'{path}: tenv-HVG rank is not a complete 1..{axis_size} permutation')
    _published_digest(path, Path(path).with_name('manifest.json'),
                      lambda manifest: manifest.get('populations', {})
                      .get('descriptive', {}).get('sha256'))
    return frame


def _require_digest_pinned(record, key):
    """A record's provenance must be digest-pinned; no unverified path is accepted."""
    if key not in record:
        raise ValueError(f'{record.get("dataset_id", "?")}: record lacks {key}')
    if key == 'inputs':
        value = record['inputs']
        if not isinstance(value, dict) or not value:
            raise ValueError(f'{record["dataset_id"]}: inputs are not a pinned mapping')
        for name, evidence in value.items():
            digest = evidence if isinstance(evidence, str) else (
                evidence.get('sha256') if isinstance(evidence, dict) else None)
            if not isinstance(digest, str) or len(digest) != 64:
                raise ValueError(f'{record["dataset_id"]}: input {name} has no sha256')
        return
    value = record[key]
    if not isinstance(value, list) or not value:
        raise ValueError(f'{record["dataset_id"]}: {key} are not a pinned list')
    for evidence in value:
        if (not isinstance(evidence, dict)
                or not str(evidence.get('relative', ''))
                or len(str(evidence.get('sha256', ''))) != 64):
            raise ValueError(f'{record["dataset_id"]}: {key} entry has no relative + sha256')


def validate_catalog_shape(category, catalog, axis_size):
    """Fail closed on any catalogue that does not match the canonical #187-#189 shape.

    ``datasets`` must be the list of full records (no ``dataset_records``/id-map
    fallback), every record must carry its own ``method``, a positive
    ``n_case``/``n_control``, and digest-pinned provenance; the category must
    carry ``catalog_sha256`` (verified by :class:`bases.disease_deg.DiseaseDeg`).
    A catalog declaring the within-dataset rank contract must pin its derived
    item-local inputs as ``derived_inputs``; the cardiotoxicity catalogs pin the
    upstream item bytes as ``inputs`` and declare no rank contract.
    """
    if not isinstance(catalog.get('datasets'), list) or not catalog['datasets']:
        raise ValueError(f'{category}: catalog datasets is not a non-empty list')
    if len(str(catalog.get('catalog_sha256', ''))) != 64:
        raise ValueError(f'{category}: catalog has no canonical catalog_sha256')
    if int(catalog.get('axis_size', axis_size)) != axis_size:
        raise ValueError(f'{category}: catalog axis_size differs from the frozen axis')
    rank_contract = catalog.get('rank_contract')
    for record in catalog['datasets']:
        dataset_id = record.get('dataset_id')
        if not dataset_id:
            raise ValueError(f'{category}: a dataset record has no dataset_id')
        for key in ('deg', 'rank', 'n_measured', 'n_significant'):
            if key not in record:
                raise ValueError(f'{category}/{dataset_id}: record lacks {key}')
        method = record.get('method')
        if not isinstance(method, dict) or not method:
            raise ValueError(f'{category}/{dataset_id}: record lacks its own method')
        for key in ('n_case', 'n_control'):
            value = record.get(key)
            if value is None or float(value) <= 0:
                raise ValueError(f'{category}/{dataset_id}: {key} must be positive')
        _require_digest_pinned(record, 'inputs')
        if rank_contract == RANK_SCHEMA:
            _require_digest_pinned(record, 'derived_inputs')


def read_categories(config, axis_size):
    """Read every category through ``bases.disease_deg.DiseaseDeg`` (digest-checked)."""
    reader = DiseaseDeg(config.disease_root)
    found = tuple(reader.categories())
    if len(found) != len(CATEGORIES) or set(found) != set(CATEGORIES):
        raise ValueError(f'expected the six declared categories {CATEGORIES}, found {found}')
    categories = {}
    for category in CATEGORIES:
        catalog = reader.catalog(category)
        validate_catalog_shape(category, catalog, axis_size)
        selection = winning_tier(catalog)
        manifest_path = config.disease_root / category / 'manifest.json'
        category_manifest = json.loads(manifest_path.read_text())
        if category_manifest.get('schema') not in CATEGORY_MANIFEST_SCHEMAS:
            raise ValueError(f'{manifest_path}: unsupported category release schema')
        if category_manifest.get('catalog_sha256') != catalog['catalog_sha256']:
            raise ValueError(f'{manifest_path}: catalog digest differs from the catalog')
        declared_outputs = category_manifest.get('outputs')
        if not isinstance(declared_outputs, dict) or not declared_outputs:
            raise ValueError(f'{manifest_path}: no declared outputs')
        ranks, measured_counts, datasets = [], [], []
        for record in selection['records']:
            rank_path = config.disease_root / category / record['rank']
            declared = declared_outputs.get(record['rank'])
            declared_sha = (declared.get('sha256') if isinstance(declared, dict)
                            else declared)
            if not isinstance(declared_sha, str) or len(declared_sha) != 64:
                raise ValueError(f'{category}/{record["dataset_id"]}: '
                                 'rank is not a declared output')
            if sha256(rank_path) != declared_sha:
                raise ValueError(f'{category}/{record["dataset_id"]}: '
                                 'rank digest differs from the category manifest')
            frame = pd.read_csv(rank_path, dtype={'gene_index': np.int64},
                                float_precision='round_trip')
            if list(frame.columns) != ['gene_index', 'rank', 'rank_norm']:
                raise ValueError(f'{category}/{record["dataset_id"]}: unexpected rank columns')
            gene_index = np.asarray(frame.gene_index, dtype=np.int64)
            if len(frame) != axis_size or not np.array_equal(
                    gene_index, np.arange(axis_size, dtype=np.int64)):
                raise ValueError(f'{category}/{record["dataset_id"]}: '
                                 'rank frame does not cover the frozen axis')
            n_measured = int(record['n_measured'])
            rank = np.asarray(frame['rank'], dtype=np.int64)
            expected = normalized_rank(rank, n_measured)
            # The producers write rank_norm with finite print precision; the fusion
            # itself recomputes rank / n_measured from the integer rank, never the
            # printed column. The tolerance only guards the published contract.
            if not np.allclose(np.asarray(frame.rank_norm, dtype=np.float64),
                               expected, rtol=1e-9, atol=1e-12):
                raise ValueError(f'{category}/{record["dataset_id"]}: '
                                 'published rank_norm differs from rank / n_measured')
            ranks.append(rank)
            measured_counts.append(n_measured)
            datasets.append({'dataset_id': record['dataset_id'],
                             'tier': int(record['tier']),
                             'evidence_type': record.get('evidence_type'),
                             'n_measured': n_measured,
                             'rank_sha256': declared_sha})
        categories[category] = {
            'catalog': catalog,
            'selection': selection,
            'ranks': ranks,
            'n_measured': measured_counts,
            'datasets': datasets,
            'catalog_sha256': catalog['catalog_sha256'],
            'manifest_sha256': sha256(manifest_path),
            'catalog_sha256_file': sha256(config.disease_root / category / 'catalog.json'),
        }
    return categories


# ---------------------------------------------------------------- computation

def rank_sum_curve(rank_sum):
    """The sorted ascending rank-sum curve: x = rank position, y = rank_sum."""
    rank_sum = np.asarray(rank_sum, dtype=np.float64)
    ordered = np.sort(rank_sum)
    return np.arange(1, len(ordered) + 1, dtype=np.int64), ordered


def compute(axis, tenv_frame, categories, axis_sizes):
    """Apply the ADR 0066 fusion and fill rule; returns every product payload."""
    axis_size = len(axis)
    gene_index = np.asarray(axis.gene_index, dtype=np.int64)
    tenv_rank = np.asarray(tenv_frame['rank'], dtype=np.int64)
    level_one = {}
    category_ranks, category_measured = [], []
    for category in CATEGORIES:
        entry = categories[category]
        result = category_rank(entry['ranks'], entry['n_measured'],
                               gene_index=gene_index, axis_size=axis_size)
        level_one[category] = result
        category_ranks.append(result['rank'])
        category_measured.append(result['measured'])
    fusion = fusion_rank_sum(category_ranks, category_measured, tenv_rank,
                             gene_index=gene_index, axis_size=axis_size,
                             sizes=axis_sizes)
    table = candidate_table(axis, tenv_rank, category_ranks, category_measured,
                            gene_index=gene_index, sizes=axis_sizes)
    axes = {}
    for size in axis_sizes:
        axes[int(size)] = axis_fill(tenv_rank, fusion['pn'][int(size)],
                                    fusion['rank_sum'], n=int(size),
                                    gene_index=gene_index)
    intervals = {}
    for size, genes in axes.items():
        values = fusion['rank_sum'][genes]
        intervals[size] = {'low': float(values.min()),
                           'high': float(values.max()),
                           'median': float(np.median(values))}
    histogram = rank_sum_histogram(fusion['rank_sum'], low=0.0,
                                   high=RANK_SUM_MAX, bins=HISTOGRAM_BINS)
    curve_x, curve_y = rank_sum_curve(fusion['rank_sum'])
    inflection = detect_inflection(curve_x, curve_y, min_n=INFLECTION_MIN_N)
    replay_curves = [(curve_x, curve_y)]
    for dropped in range(N_FUSION_SOURCES):
        sources = category_ranks + [tenv_rank]
        replay_sum = np.zeros(axis_size, dtype=np.float64)
        for index, rank in enumerate(sources):
            if index == dropped:
                continue
            replay_sum += rank / float(axis_size)
        replay_curves.append(rank_sum_curve(replay_sum))
    stability = inflection_stability(replay_curves)
    return {
        'level_one': level_one,
        'category_ranks': category_ranks,
        'category_measured': category_measured,
        'fusion': fusion,
        'table': table,
        'axes': axes,
        'intervals': intervals,
        'histogram': histogram,
        'inflection': inflection,
        'stability': stability,
    }


def composition_rows(source, axes, axis_sizes):
    """Per-category and pN composition of both released axes."""
    rows = []
    for size in axis_sizes:
        size = int(size)
        genes = source['axes'][size]
        composition = axis_composition(genes, source['category_measured'])
        for category in CATEGORIES:
            count = composition[category]
            rows.append({'axis': str(size), 'n_axis': size, 'record': 'composition',
                         'group': category, 'n_genes': count,
                         'share': count / size,
                         'note': 'axis genes the category source measures'})
        rows.append({'axis': str(size), 'n_axis': size, 'record': 'composition',
                     'group': 'tenv_hvg_top_half', 'n_genes': size // 2,
                     'share': 0.5,
                     'note': 'declared 0.5*N tenv-HVG half of the axis'})
        axis_counts = pn_distribution(source['table'][f'p{size}'].to_numpy()[genes])
        for value, count in enumerate(axis_counts.tolist()):
            rows.append({'axis': str(size), 'n_axis': size, 'record': 'pN',
                         'group': str(value), 'n_genes': count,
                         'share': count / size, 'note': COMPOSITION_NOTE})
        candidate_counts = pn_distribution(source['table'][f'p{size}'].to_numpy())
        for value, count in enumerate(candidate_counts.tolist()):
            rows.append({'axis': str(size), 'n_axis': size, 'record': 'candidate_pN',
                         'group': str(value),
                         'n_genes': count,
                         'share': count / len(source['table']),
                         'note': COMPOSITION_NOTE})
    return pd.DataFrame(rows, columns=list(COMPOSITION_COLUMNS))


def figure_rows(source, axis_sizes):
    """The complete plotted values of the rank-sum histogram evidence panel."""
    histogram = source['histogram']
    total = int(histogram['count'].sum())
    rows = [{'record': 'figure', 'title': FIGURE_TITLE},
            {'record': 'panel', 'panel': 'a',
             'panel_title': 'rank_sum over all candidates with both axis intervals'}]
    for row in histogram.itertuples(index=False):
        rows.append({'record': 'hist',
                     'bin_left': row.bin_left, 'bin_right': row.bin_right,
                     'bin_center': row.bin_center, 'count': int(row.count),
                     'share': int(row.count) / total})
    for size in axis_sizes:
        interval = source['intervals'][int(size)]
        rows.append({'record': 'interval', 'axis': str(int(size)),
                     'interval_low': interval['low'],
                     'interval_high': interval['high'],
                     'median': interval['median']})
    inflection = source['inflection']
    note = ('no inflection' if not inflection['inflection']
            else f'candidate N={inflection["n"]}')
    rows.append({'record': 'marker', 'inflection': inflection['inflection'],
                 'inflection_n': inflection['n'] if inflection['n'] is not None else '',
                 'note': note})
    return pd.DataFrame(rows, columns=list(FIGURE_COLUMNS))


# ---------------------------------------------------------------- release I/O

def input_entries(config, axis_size):
    """Every pinned input with its role and digest (repo-relative paths only)."""
    entries = {}
    for role, path in (('frozen_axis', config.genes), ('tenv_hvg_descriptive', config.tenv_hvg)):
        entries[repo_relative(path, config.root)] = {
            'role': role, 'sha256': sha256(path), 'bytes': path.stat().st_size}
    reader = DiseaseDeg(config.disease_root)
    for category in CATEGORIES:
        directory = config.disease_root / category
        for role, path in (('catalog', directory / 'catalog.json'),
                           ('release_manifest', directory / 'manifest.json')):
            entries[repo_relative(path, config.root)] = {
                'role': f'{category}:{role}', 'sha256': sha256(path),
                'bytes': path.stat().st_size}
        catalog = reader.catalog(category)
        selection = winning_tier(catalog)
        for record in selection['records']:
            rank_path = directory / record['rank']
            entries[repo_relative(rank_path, config.root)] = {
                'role': f'{category}:{record["dataset_id"]}:rank',
                'sha256': sha256(rank_path), 'bytes': rank_path.stat().st_size}
    return dict(sorted(entries.items()))


def _axis_manifest(config, axis, tenv_frame, categories, source, axis_sizes):
    axis_size = len(axis)
    fusion = source['fusion']
    table = source['table']
    composition = composition_rows(source, source['axes'], axis_sizes)
    axes_manifest, intersection = {}, None
    for size in axis_sizes:
        size = int(size)
        genes = source['axes'][size]
        composition_counts = axis_composition(genes, source['category_measured'])
        counts = pn_distribution(fusion['pn'][size])
        axes_manifest[str(size)] = {
            'n': int(len(genes)),
            'file': f'gene_axis_{size}.csv',
            'tenv_half': size // 2,
            'fill_half': size // 2,
            'rank_sum': source['intervals'][size],
            'tenv_rank': {'min': int(tenv_frame['rank'].to_numpy()[genes].min()),
                          'max': int(tenv_frame['rank'].to_numpy()[genes].max())},
            'pN_distribution': {str(value): int(count)
                                for value, count in enumerate(counts.tolist())},
            'composition': {category: int(count)
                            for category, count in composition_counts.items()},
        }
    left, right = axis_sizes[0], axis_sizes[1]
    intersection = axis_intersection(source['axes'][left], source['axes'][right])
    intersection = {'sizes': [int(left), int(right)], **intersection}
    category_manifest = {}
    for category in CATEGORIES:
        entry = categories[category]
        catalog = entry['catalog']
        category_manifest[category] = {
            'label': CATEGORY_LABELS[category],
            'product_dir': category,
            'catalog_sha256': entry['catalog_sha256'],
            'catalog_file_sha256': entry['catalog_sha256_file'],
            'release_manifest_sha256': entry['manifest_sha256'],
            'category_tier': catalog.get('tier'),
            'winning_tier': entry['selection']['winning_tier'],
            'winning_tier_rule': entry['selection']['rule'],
            'retained_by_exception': entry['selection']['retained_by_exception'],
            'excluded_tiers': entry['selection']['excluded_tiers'],
            'datasets': entry['datasets'],
            'n_measured_axis_genes': int(np.asarray(
                source['category_measured'][CATEGORIES.index(category)]).sum()),
        }
    return {
        'schema': SCHEMA,
        'generated_by': 'stage2/01_bases/execution/build_gene_axis.py',
        'release': {
            'axis_sizes': [int(size) for size in axis_sizes],
            'fill_rule': FILL_RULE,
            'fill_fraction': 0.5,
            'pN': ('number of the seven fusion sources whose top-N set contains the '
                   'gene; ordinal, not a continuous quality score'),
        },
        'environment': environment_facts(),
        'inputs': input_entries(config, axis_size),
        'fusion': {
            'n_sources': N_FUSION_SOURCES,
            'sources': list(CATEGORIES) + ['tenv_hvg'],
            'normalization': ('level-1 rank_norm = rank / n_measured per dataset; '
                              'level-2 each of the seven source ranks / frozen-axis size'),
            'rank_sum': 'ascending is better; ties on gene_index',
            'n_candidates': int(len(table)),
        },
        'tenv_hvg': {
            'input': repo_relative(config.tenv_hvg, config.root),
            'n_ranked': int(axis_size),
            'source': 'published descriptive normalized-dispersion rank (never recomputed)',
            'normalization': 'rank / frozen-axis size',
        },
        'categories': category_manifest,
        'counts': {
            'n_axis': int(axis_size),
            'n_candidates': int(len(table)),
            'rank_sum_min': float(fusion['rank_sum'].min()),
            'rank_sum_median': float(np.median(fusion['rank_sum'])),
            'rank_sum_max': float(fusion['rank_sum'].max()),
            'pN_distribution': {
                str(size): {str(value): int(count) for value, count in enumerate(
                    pn_distribution(fusion['pn'][int(size)]).tolist())}
                for size in axis_sizes},
            'n_sources_measured_distribution': {
                str(value): int(count) for value, count in enumerate(
                    np.bincount(np.asarray(table.n_sources_measured, dtype=np.int64),
                                minlength=N_FUSION_SOURCES + 1).tolist())},
        },
        'axes': axes_manifest,
        'intersection': intersection,
        'composition': {
            'file': COMPOSITION_NAME,
            'note': COMPOSITION_NOTE,
            'rows': int(len(composition)),
        },
        'diagnostics': {
            'inflection': {
                **{key: value for key, value in source['inflection'].items()},
                'role': 'diagnostic only; N is declared and never read off the curve',
                'stability': source['stability'],
            },
            'histogram': {
                'file': f'figure/data/{FIGURE_NAME}.csv',
                'low': 0.0, 'high': RANK_SUM_MAX, 'bins': HISTOGRAM_BINS,
                'axis_intervals': {str(size): source['intervals'][size]
                                   for size in axis_sizes},
            },
        },
    }


def environment_facts():
    return {
        'python': platform.python_version(),
        'numpy': np.__version__,
        'pandas': pd.__version__,
        'data_access': ('offline; reads only output/gene-identity, output/tenv-hvg and '
                        'output/disease-deg; never recomputes the tenv-HVG ranking'),
        'resource_contract': ('one heavy job at a time, 10 GiB RLIMIT_AS fallback, '
                              'OMP_NUM_THREADS=2; peak RSS recorded by the driver'),
    }


def product_paths(config, axis_sizes):
    """Repo-relative names of every product (manifest keys)."""
    output = repo_relative(config.output, config.root)
    figure = repo_relative(config.figure_data, config.root)
    names = {
        'candidate': f'{output}/{CANDIDATE_NAME}',
        'composition': f'{output}/{COMPOSITION_NAME}',
        'manifest': f'{output}/manifest.json',
        'figure': figure,
    }
    for size in axis_sizes:
        names[f'axis_{int(size)}'] = f'{output}/gene_axis_{int(size)}.csv'
    return names


def _build_payload(config, payload_root):
    """Build every product into ``payload_root`` (repo-relative mirror) and return the manifest."""
    paths = product_paths(config, config.axis_sizes)
    axis = read_axis(config.genes)
    axis_size = len(axis)
    tenv_frame = read_tenv(config.tenv_hvg, axis_size)
    categories = read_categories(config, axis_size)
    source = compute(axis, tenv_frame, categories, config.axis_sizes)

    outputs = {}
    candidate_path = payload_root / paths['candidate']
    write_frame(source['table'], candidate_path)
    outputs[paths['candidate']] = {'rows': int(len(source['table'])),
                                   'sha256': sha256(candidate_path),
                                   'columns': list(candidate_columns(config.axis_sizes))}
    for size in config.axis_sizes:
        size = int(size)
        genes = source['axes'][size]
        frame = pd.DataFrame({'gene_index': genes,
                              'symbol': axis.symbol.to_numpy()[genes]})
        path = payload_root / paths[f'axis_{size}']
        write_frame(frame, path)
        outputs[paths[f'axis_{size}']] = {'rows': int(len(frame)),
                                          'sha256': sha256(path)}
    composition = composition_rows(source, source['axes'], config.axis_sizes)
    composition_path = payload_root / paths['composition']
    write_frame(composition, composition_path)
    outputs[paths['composition']] = {'rows': int(len(composition)),
                                     'sha256': sha256(composition_path)}
    figure = figure_rows(source, config.axis_sizes)
    figure_path = payload_root / paths['figure']
    write_frame(figure, figure_path)
    outputs[paths['figure']] = {'rows': int(len(figure)),
                                'sha256': sha256(figure_path)}

    manifest = _axis_manifest(config, axis, tenv_frame, categories, source,
                             config.axis_sizes)
    manifest['outputs'] = dict(sorted(outputs.items()))
    manifest_path = payload_root / paths['manifest']
    write_json(manifest, manifest_path)
    validate_staged(config, source, payload_root, manifest, axis_size)
    return manifest


def validate_staged(config, source, payload_root, manifest, axis_size):
    """Structural acceptance of the candidate before it may replace the release."""
    table = source['table']
    if len(table) != axis_size:
        raise ValueError('candidate table does not have one row per frozen-axis gene')
    if list(table.columns) != list(candidate_columns(config.axis_sizes)):
        raise ValueError('candidate table columns differ from the contract')
    if table.gene_index.tolist() != list(range(axis_size)):
        raise ValueError('candidate table is not in frozen gene_index order')
    for size in config.axis_sizes:
        size = int(size)
        column = f'p{size}'
        if int(table[column].min()) < 0 or int(table[column].max()) > N_FUSION_SOURCES:
            raise ValueError(f'{column} is outside 0..{N_FUSION_SOURCES}')
        genes = source['axes'][size]
        if len(genes) != size:
            raise ValueError(f'G_{size} does not hold exactly {size} genes')
        if not np.array_equal(genes, np.sort(genes)):
            raise ValueError(f'G_{size} is not ordered by ascending gene_index')
        if len(np.unique(genes)) != size or genes.min() < 0 or genes.max() >= axis_size:
            raise ValueError(f'G_{size} is not a distinct subset of the frozen axis')
        expected = axis_fill(np.asarray(source['table'].tenv_hvg_rank, dtype=np.int64),
                             table[column].to_numpy(), table.rank_sum.to_numpy(),
                             n=size)
        if not np.array_equal(expected, genes):
            raise ValueError(f'G_{size} does not follow the declared fill rule')
    release.check_relative_paths(manifest)
    for name, record in manifest['outputs'].items():
        actual = sha256(payload_root / name)
        if actual != record['sha256']:
            raise ValueError(f'staged digest differs for {name}')
    release.check_outputs(payload_root, manifest)


def _stage_root(config):
    return Path(config.root) / 'stage2/01_bases/execution/checkpoints'


def publish(config):
    """Build, validate and atomically replace the accepted release."""
    paths = product_paths(config, config.axis_sizes)
    staging = release.candidate(config.output, root=_stage_root(config))
    try:
        payload = staging / 'payload'
        manifest = _build_payload(config, payload)
        # Candidate-then-replace: only this ticket's products are swapped in. The
        # pre-0066 gene-axis file set stays where it is until the removal ticket
        # (#192) deletes it after acceptance (ADR 0066 Removals / ADR 0062), so
        # the release directory is never renamed wholesale.
        output_root = payload / repo_relative(config.output, config.root)
        config.output.mkdir(parents=True, exist_ok=True)
        for key in ('candidate', 'composition'):
            os.replace(output_root / Path(paths[key]).name,
                       config.output / Path(paths[key]).name)
        for size in config.axis_sizes:
            name = Path(paths[f'axis_{int(size)}']).name
            os.replace(output_root / name, config.output / name)
        config.figure_data.parent.mkdir(parents=True, exist_ok=True)
        os.replace(payload / paths['figure'], config.figure_data)
        # The manifest is replaced last: it only announces products already in place.
        os.replace(output_root / 'manifest.json', config.output / 'manifest.json')
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    return summary(config, manifest)


def verify(config):
    """Rebuild from the pinned inputs and compare the manifest and every digest."""
    paths = product_paths(config, config.axis_sizes)
    accepted_path = config.output / 'manifest.json'
    if not accepted_path.is_file():
        raise ValueError(f'no accepted gene-axis manifest at {accepted_path}')
    accepted = json.loads(accepted_path.read_text())
    if accepted.get('schema') != SCHEMA:
        raise ValueError('unknown gene-axis release schema')
    release.check_manifest(config.root, accepted)
    for name, evidence in accepted['inputs'].items():
        if sha256(config.root / name) != evidence['sha256']:
            raise ValueError(f'pinned input digest differs: {name}')
    staging = release.candidate(config.output, root=_stage_root(config))
    try:
        payload = staging / 'payload'
        replayed = _build_payload(config, payload)
        if replayed != accepted:
            raise ValueError('replayed gene-axis manifest differs from the accepted release')
        for name, record in accepted['outputs'].items():
            if sha256(payload / name) != record['sha256']:
                raise ValueError(f'replayed digest differs: {name}')
        release.check_outputs(payload, accepted)
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    return {
        'result': 'PASS', 'schema': SCHEMA,
        'manifest_sha256': sha256(accepted_path),
        'inputs': len(accepted['inputs']),
        'outputs': len(accepted['outputs']),
        'axes': accepted['release']['axis_sizes'],
        'intersection': accepted['intersection']['n_intersection'],
        'inflection': accepted['diagnostics']['inflection']['inflection'],
    }


def summary(config, manifest):
    import resource
    return {
        'result': 'PASS',
        'schema': manifest['schema'],
        'n_candidates': manifest['counts']['n_candidates'],
        'axes': manifest['release']['axis_sizes'],
        'intersection': manifest['intersection']['n_intersection'],
        'inflection': manifest['diagnostics']['inflection']['inflection'],
        'manifest_sha256': sha256(config.output / 'manifest.json'),
        'outputs': {name: record['sha256'] for name, record in manifest['outputs'].items()},
        'peak_rss_kib': int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--root', type=Path)
    parser.add_argument('--genes', type=Path)
    parser.add_argument('--tenv-hvg', type=Path)
    parser.add_argument('--disease-root', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--figure-data', type=Path)
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--record', type=Path)
    args = parser.parse_args()
    config = default_config(args.root, args.output, args.figure_data)
    overrides = {'genes': args.genes, 'tenv_hvg': args.tenv_hvg,
                 'disease_root': args.disease_root}
    config = replace(config, **{key: value for key, value in overrides.items() if value})
    result = verify(config) if args.verify else publish(config)
    text = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.record:
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
