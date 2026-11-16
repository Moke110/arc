#!/usr/bin/env python3
"""Build or replay the released Stage 2 gene axes ``G_2000`` / ``G_5000`` (ADR 0066, ADR 0071).

The products are the cross-disease fusion of the published per-disease ``DEG2000``
sets plus the published tenv-HVG ranking: the master candidate table, the two
released axes, their manifest, the per-axis composition table, the source
coverage table, the two-panel rank-sum evidence figure values and the two-panel
source-coverage figure values. No corpus, view, model or truth is built; the
disease-fusion products (#205) and the tenv-HVG ranking are read-only inputs.

The disease side is **DEG-only** (ADR 0071): each disease contributes its
published ``DEG2000``, the axis disease half is ordered by ``f2000`` then
``rank_sum`` and any shortfall is filled from unused tenv-HVG genes (never from a
gene on disease evidence it does not have). The retired ``pN`` count is not read,
written or bounded here.

Every scientific definition lives in ``output/api/bases/gene_axis.py`` (the single
public implementation, ADR 0062); this driver only reads the pinned inputs,
applies those definitions and publishes the products.

    python stage2/01_bases/execution/build_gene_axis.py
    python stage2/01_bases/execution/build_gene_axis.py --verify

Publication is candidate-then-replace (ADR 0062): a complete candidate is built
outside ``output/``, validated, and only then swapped in; a failed build leaves
the accepted release usable. ``--verify`` rebuilds from the pinned inputs into a
scratch candidate and compares the manifest and every output digest. The
committed release predates the v2 rule until #208 republishes it, so ``--verify``
refuses it with a self-describing message instead of an ambiguous byte mismatch.
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
from bases.disease_deg import (DEG2000_COLUMNS, DEG2000_LIMIT, DiseaseFusion,
                               disease_eligibility)
from bases.gene_axis import (
    ALGORITHM, AXIS_SIZES, CATEGORIES, CATEGORY_LABELS, COMPOSITION_COLUMNS,
    COVERAGE_COLUMNS, F2000_RULE, FALLBACK_RULE, FILL_RULE, RANK_SUM_RULE,
    RELEASE_CATEGORIES, SCHEMA, axis_fill, axis_intersection, candidate_columns,
    candidate_table, cross_disease_rank_sum, deg2000_composition, deg2000_frame,
    detect_inflection, f2000_distribution, inflection_stability, rank_sum_histogram,
    source_topn_coverage,
)

CANDIDATE_NAME = 'candidate_rank_sum.csv'
COMPOSITION_NAME = 'axis_composition.csv'
COVERAGE_NAME = 'source_topn_coverage.csv'
FIGURE_NAME = 'gene_axis_rank_sum_histogram'
HISTOGRAM_BINS = 70
#: Inflection is diagnostic only; the declared decision range starts at N=100.
INFLECTION_MIN_N = 100
COMPOSITION_NOTE = (
    'share = n_genes / denominator. record=axis_fill reports the declared 0.5*N '
    'tenv-HVG half, the disease-ordered half and the tenv-HVG fallback; '
    'record=deg2000 reports |DEG2000_d ∩ G_N| over |DEG2000_d|; record=f2000 '
    'reports how many axis genes have each disease-replication count, over N.')
COVERAGE_NOTE = (
    'one row per axis and fusion group; numerator/denominator are explicit and '
    'rate = 100 * numerator / denominator. tenv_hvg is |top_N(tenv) ∩ G_N| / N; '
    'a disease is |DEG2000_d ∩ G_N| / |DEG2000_d| (ADR 0071 decision 8). A group '
    'excluded from the fusion keeps an empty row.')

#: The evidence-figure plotted values: one panel per released axis, the gray
#: candidate histogram and the count of that bin's candidates each axis selects.
FIGURE_COLUMNS = ('record', 'panel', 'panel_title', 'title', 'axis', 'bin_left',
                  'bin_right', 'bin_center', 'count', 'selected')
FIGURE_TITLE = ('Released Stage 2 gene axes · rank-sum distribution of the frozen-axis '
                "candidates with each axis's selected genes")

#: The source-coverage figure plotted values: one x group per fusion group and one
#: bar series per released axis. Panel ``a`` plots the absolute DEG2000 numerator
#: the axis contains (``count``) and panel ``b`` its rate against the group's own
#: denominator (``rate``); both values are carried on every bar so a panel
#: declares which column it draws. The plotted values are the coverage table's
#: numerator/rate columns, so the figure and the table cannot disagree.
SOURCE_COVERAGE_FIGURE_NAME = 'gene_axis_source_topn_coverage'
SOURCE_COVERAGE_FIGURE_COLUMNS = (
    'record', 'panel', 'panel_title', 'y_label', 'value',
    'series', 'series_label', 'source', 'source_label', 'count', 'rate')
SOURCE_COVERAGE_FIGURE_PANELS = (
    {'panel': 'a', 'title': 'Covered Gene Count', 'y_label': 'covered genes',
     'value': 'count'},
    # The rate's denominator is no longer N for every group: it is |DEG2000_d|
    # for a disease and N only for tenv-HVG (ADR 0071 decision 8), so the v1
    # "Top-N" wording misdescribes six of the seven groups.
    {'panel': 'b', 'title': 'Coverage of Each Group',
     'y_label': 'share of group covered (%)', 'value': 'rate'},
)
SOURCE_COVERAGE_FIGURE_LABELS = {'tenv_hvg': 'tenv-HVG', **CATEGORY_LABELS}


@dataclass(frozen=True)
class Config:
    """Every path the build reads or writes; all are repo-relative by default."""

    root: Path
    genes: Path
    tenv_hvg: Path
    fusion_root: Path
    output: Path
    figure_data: Path
    source_coverage_figure: Path
    axis_sizes: tuple = AXIS_SIZES


def repo_root():
    return next(parent for parent in Path(__file__).resolve().parents
                if (parent / '.git').exists())


def default_config(root=None, output=None, figure_data=None,
                   source_coverage_figure=None, axis_sizes=AXIS_SIZES):
    root = Path(root) if root else repo_root()
    step = root / 'stage2/01_bases'
    figure_data = (Path(figure_data) if figure_data
                   else step / 'figure/data' / f'{FIGURE_NAME}.csv')
    return Config(
        root=root,
        genes=step / 'output/gene-identity/genes.csv',
        tenv_hvg=step / 'output/tenv-hvg/descriptive.csv',
        fusion_root=step / 'output/disease-fusion',
        output=Path(output) if output else step / 'output/gene-axis',
        figure_data=figure_data,
        source_coverage_figure=(Path(source_coverage_figure) if source_coverage_figure
                                else figure_data.parent / f'{SOURCE_COVERAGE_FIGURE_NAME}.csv'),
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


def write_frame(frame, path, float_format='%.17g'):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(frame.to_csv(index=False, lineterminator='\n',
                                 float_format=float_format))
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
    """The published tenv-HVG dispersion rank: the axis's first-half source."""
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


def read_fusion(config, axis_size):
    """Read the published per-disease ``DEG2000`` products through ``DiseaseFusion``.

    ADR 0071 decisions 4/5: the axis disease side is built from each disease's
    published ``DEG2000`` -- the tiered, axis-resolved level-1 fusion product --
    never from a within-dataset measured-gene rank frame or from ADR 0066's
    minimum-tier dataset selection. Every released category must be published;
    a declared-but-unpublished category (``AF`` until #208) is tolerated, and a
    published category outside :data:`RELEASE_CATEGORIES` fails with the #208
    pointer instead of silently entering the axis. A published disease with an
    empty ``DEG2000`` is representable and is excluded from the fusion with its
    reason recorded (ADR 0068 has no ordering evidence to fuse).

    Returns ``{'catalog', 'manifest', 'catalog_file_sha256',
    'manifest_file_sha256', 'diseases', 'excluded'}``; ``diseases`` maps a
    released category to its complete frozen-axis ``DEG2000`` rank frame (the
    published ``disease_rank`` inside ``DEG2000_d``, ``n_deg + 1`` elsewhere), the
    union/DEG2000 sizes and the file digest.
    """
    reader = DiseaseFusion(config.fusion_root)
    catalog = reader.catalog()
    fusion_manifest = reader.manifest()
    if int(catalog['axis']['size']) != int(axis_size):
        raise ValueError('the disease-fusion catalog was built on a different frozen axis')
    if sha256(config.genes) != catalog['axis']['sha256']:
        raise ValueError(
            'the disease-fusion catalog pins a different frozen-axis digest than '
            f'{config.genes}')
    if int(catalog.get('limit', -1)) != int(DEG2000_LIMIT):
        raise ValueError(f'the disease-fusion DEG2000 limit must be {DEG2000_LIMIT}')
    by_directory = {record['directory']: record for record in catalog['diseases']}
    published = {directory for directory, record in by_directory.items()
                 if record['published']}
    unknown = sorted(published - set(CATEGORIES))
    if unknown:
        raise ValueError(
            f'undeclared disease-fusion product directories present: {unknown}; '
            f'the declared roster is {CATEGORIES}')
    missing = sorted(set(RELEASE_CATEGORIES) - published)
    if missing:
        raise ValueError(
            f'the disease-fusion catalog is missing the released diseases {missing}; '
            '#208 must publish their DEG2000 products before the axes can be rebuilt')
    unpublished = sorted(published - set(RELEASE_CATEGORIES))
    if unpublished:
        raise ValueError(
            f'the declared diseases {unpublished} have disease-fusion products but '
            'are not in gene_axis.RELEASE_CATEGORIES; #208 must publish them as one '
            'release (move them out of UNPUBLISHED_CATEGORIES) before the axes can '
            'be rebuilt')
    declared_outputs = fusion_manifest.get('outputs')
    if not isinstance(declared_outputs, dict) or not declared_outputs:
        raise ValueError(f'{config.fusion_root}/manifest.json declares no outputs')
    diseases, excluded = {}, []
    for category in RELEASE_CATEGORIES:
        record = by_directory[category]
        n_deg = int(record['n_deg'])
        deg2000 = int(record['deg2000'])
        admission = disease_eligibility(n_deg)
        # Old releases retain their original meaning until the coherent successor
        # is published; the explicit policy marks ADR 0072 candidate products.
        gated = catalog.get('eligibility_min_deg') is not None
        if (gated and not admission['eligible']) or n_deg <= 0 or deg2000 <= 0:
            excluded.append({'directory': category, 'label': CATEGORY_LABELS[category],
                             'n_deg': n_deg, 'deg2000': deg2000,
                             **admission})
            continue
        if deg2000 != min(n_deg, DEG2000_LIMIT):
            raise ValueError(
                f'{category}: published deg2000 must be min(n_deg, {DEG2000_LIMIT})')
        name = record['deg2000_file']
        declared = declared_outputs.get(name)
        digest = declared.get('sha256') if isinstance(declared, dict) else declared
        if not isinstance(digest, str) or len(digest) != 64:
            raise ValueError(f'{category}: {name} is not a declared fusion output')
        rows = reader.deg2000(category)  # digest-checked against the fusion manifest
        if list(rows.columns) != list(DEG2000_COLUMNS):
            raise ValueError(f'{category}: {name} has unexpected columns')
        if len(rows) != deg2000:
            raise ValueError(f'{category}: {name} holds {len(rows)} rows, expected {deg2000}')
        genes = np.asarray(rows.gene_index, dtype=np.int64)
        ranks = np.asarray(rows.disease_rank, dtype=np.int64)
        if genes.min() < 0 or genes.max() >= int(axis_size):
            raise ValueError(f'{category}: {name} holds a gene_index outside the frozen axis')
        if len(np.unique(genes)) != len(genes):
            raise ValueError(f'{category}: {name} repeats a gene_index')
        if not np.array_equal(np.sort(ranks),
                              np.arange(1, deg2000 + 1, dtype=np.int64)):
            raise ValueError(f'{category}: {name} ranks must be exactly 1..{deg2000}')
        frame = np.full(int(axis_size), n_deg + 1, dtype=np.int64)
        frame[genes] = ranks
        frame = deg2000_frame(frame, n_deg, int(axis_size), category)
        diseases[category] = {
            'record': record, 'n_deg': n_deg, 'deg2000': deg2000,
            'rank': frame, 'mask': frame <= deg2000,
            'deg2000_file': name, 'deg2000_sha256': digest,
        }
    if not diseases:
        raise ValueError(
            'no released disease has a published DEG2000; the axes cannot be built '
            '(every released category is empty)')
    return {
        'catalog': catalog,
        'manifest': fusion_manifest,
        'catalog_file_sha256': sha256(config.fusion_root / 'catalog.json'),
        'manifest_file_sha256': sha256(config.fusion_root / 'manifest.json'),
        'diseases': diseases,
        'excluded': excluded,
    }


# ---------------------------------------------------------------- computation

def rank_sum_curve(rank_sum):
    """The sorted ascending rank-sum curve: x = rank position, y = rank_sum."""
    rank_sum = np.asarray(rank_sum, dtype=np.float64)
    ordered = np.sort(rank_sum)
    return np.arange(1, len(ordered) + 1, dtype=np.int64), ordered


def compute(axis, tenv_frame, fusion, axis_sizes):
    """Apply the ADR 0071 cross-disease fusion and fill rule; returns every payload."""
    axis_size = len(axis)
    gene_index = np.asarray(axis.gene_index, dtype=np.int64)
    tenv_rank = np.asarray(tenv_frame['rank'], dtype=np.int64)
    active = [category for category in RELEASE_CATEGORIES
              if category in fusion['diseases']]
    excluded = [category for category in RELEASE_CATEGORIES
                if category not in fusion['diseases']]
    frames = [fusion['diseases'][category]['rank'] for category in active]
    n_deg = [fusion['diseases'][category]['n_deg'] for category in active]
    level_two = cross_disease_rank_sum(frames, n_deg, categories=active,
                                       gene_index=gene_index, axis_size=axis_size)
    table = candidate_table(axis, tenv_rank, frames, n_deg, categories=active,
                            gene_index=gene_index)
    axes, fill = {}, {}
    for size in axis_sizes:
        result = axis_fill(tenv_rank, level_two['f2000'], level_two['rank_sum'],
                           n=int(size), gene_index=gene_index,
                           n_diseases=len(active))
        axes[int(size)] = result['axis']
        fill[int(size)] = result
    intervals = {}
    for size, genes in axes.items():
        values = level_two['rank_sum'][genes]
        intervals[size] = {'low': float(values.min()),
                           'high': float(values.max()),
                           'median': float(np.median(values))}
    histogram_high = float(np.ceil(level_two['rank_sum'].max()))
    histogram = rank_sum_histogram(level_two['rank_sum'], low=0.0,
                                   high=histogram_high, bins=HISTOGRAM_BINS)
    curve_x, curve_y = rank_sum_curve(level_two['rank_sum'])
    inflection = detect_inflection(curve_x, curve_y, min_n=INFLECTION_MIN_N)
    replay_curves = [(curve_x, curve_y)]
    for dropped in range(len(active)):
        replay_sum = level_two['rank_sum'] - frames[dropped] / float(n_deg[dropped])
        replay_curves.append(rank_sum_curve(replay_sum))
    stability = inflection_stability(replay_curves)
    return {
        'active_categories': active,
        'excluded_categories': excluded,
        'frames': frames,
        'n_deg': n_deg,
        'level_two': level_two,
        'table': table,
        'axes': axes,
        'fill': fill,
        'intervals': intervals,
        'histogram_high': histogram_high,
        'histogram': histogram,
        'inflection': inflection,
        'stability': stability,
        'tenv_rank': tenv_rank,
    }


def _deg2000_masks(source):
    """One DEG2000 mask per **released** category, ``None`` for an excluded one."""
    return [source['level_two']['deg2000'][source['active_categories'].index(category)]
            if category in source['active_categories'] else None
            for category in RELEASE_CATEGORIES]


def _deg2000_sizes(source):
    """One DEG2000 denominator per **released** category, ``None`` for an excluded one."""
    return [int(source['level_two']['deg2000'][
                source['active_categories'].index(category)].sum())
            if category in source['active_categories'] else None
            for category in RELEASE_CATEGORIES]


def composition_rows(source, axis_sizes):
    """Per-axis fill, ``DEG2000`` composition and ``f2000`` distribution.

    ``record=axis_fill`` states the declared ``0.5*N`` tenv-HVG half, the
    disease-ordered half and the tenv-HVG fallback count; ``record=deg2000``
    states each disease's ``|DEG2000_d ∩ G_N|`` over ``|DEG2000_d|``;
    ``record=f2000`` states how many axis genes carry each disease count. A
    disease excluded from the fusion keeps a row with empty cells and an explicit
    note, so the table reads as "not fused", never as "fused but measured
    nothing".
    """
    active = list(source['active_categories'])
    masks = _deg2000_masks(source)
    sizes = _deg2000_sizes(source)
    rows = []
    for size in axis_sizes:
        size = int(size)
        genes = source['axes'][size]
        fill = source['fill'][size]
        half = size // 2
        rows.append({'axis': str(size), 'n_axis': size, 'record': 'axis_fill',
                     'group': 'tenv_hvg_first_half', 'n_genes': half,
                     'denominator': size, 'share': half / size,
                     'note': 'declared 0.5*N tenv-HVG half (tenv_hvg_rank <= N/2)'})
        rows.append({'axis': str(size), 'n_axis': size, 'record': 'axis_fill',
                     'group': 'deg2000_second_half',
                     'n_genes': fill['from_deg2000'], 'denominator': size,
                     'share': fill['from_deg2000'] / size,
                     'note': ('second-half genes selected from DEG2000 by f2000 '
                              'descending, rank_sum ascending, gene_index ascending')})
        rows.append({'axis': str(size), 'n_axis': size, 'record': 'axis_fill',
                     'group': 'tenv_hvg_fallback',
                     'n_genes': fill['fallback'], 'denominator': size,
                     'share': fill['fallback'] / size,
                     'note': ('second-half slots filled from unused tenv-HVG genes in '
                              'ascending tenv-HVG rank order (no disease evidence)')})
        composition = deg2000_composition(genes, masks, sizes,
                                          categories=RELEASE_CATEGORIES)
        for category in RELEASE_CATEGORIES:
            entry = composition[category]
            if entry is None:
                rows.append({'axis': str(size), 'n_axis': size, 'record': 'deg2000',
                             'group': category, 'n_genes': None,
                             'denominator': None, 'share': None,
                             'note': ('excluded from the fusion: the disease has no '
                                      'DEG2000, so no axis composition is defined')})
                continue
            rows.append({'axis': str(size), 'n_axis': size, 'record': 'deg2000',
                         'group': category, 'n_genes': entry['n_covered'],
                         'denominator': entry['n_deg2000'], 'share': entry['share'],
                         'note': ('|DEG2000_d ∩ G_N| / |DEG2000_d| for '
                                  f'{CATEGORY_LABELS[category]}')})
        counts = f2000_distribution(source['level_two']['f2000'][genes],
                                    n_diseases=len(active))
        for value, count in enumerate(counts.tolist()):
            rows.append({'axis': str(size), 'n_axis': size, 'record': 'f2000',
                         'group': str(value), 'n_genes': int(count),
                         'denominator': size, 'share': int(count) / size,
                         'note': F2000_RULE})
    frame = pd.DataFrame(rows, columns=list(COMPOSITION_COLUMNS))
    for column in ('n_genes', 'denominator'):
        frame[column] = pd.array(frame[column], dtype='Int64')
    return frame


def coverage_rows(source, axis_sizes):
    """Each fusion group's coverage of both released axes, with both numbers.

    The source ranks are the published tenv-HVG dispersion rank and each
    disease's published ``DEG2000`` mask; a non-fused disease contributes an
    empty row (ADR 0068).
    """
    return source_topn_coverage(source['axes'], source['tenv_rank'],
                                _deg2000_masks(source),
                                categories=RELEASE_CATEGORIES, sizes=axis_sizes)


def source_coverage_figure_rows(source, axis_sizes):
    """The complete plotted values of the two-panel source-coverage figure.

    One bar series per released axis (a bar's colour is its series) and one x
    group per fusion group. Every bar carries both the ``DEG2000`` numerator the
    axis contains and the rate against that group's own denominator, so each
    panel record declares which column it draws: panel ``a`` is the count and
    panel ``b`` the rate. An excluded group has no numerator and therefore no
    bar.
    """
    coverage = coverage_rows(source, axis_sizes)
    rows = [{'record': 'panel', 'panel': panel['panel'],
             'panel_title': panel['title'], 'y_label': panel['y_label'],
             'value': panel['value']}
            for panel in SOURCE_COVERAGE_FIGURE_PANELS]
    for record in coverage.to_dict('records'):
        size = int(record['N'])
        name = record['group']
        if pd.isna(record['numerator']):
            continue
        rows.append({'record': 'bar', 'series': f'G_{size}',
                     'series_label': f'G_{size} (N = {size:,})',
                     'source': name,
                     'source_label': SOURCE_COVERAGE_FIGURE_LABELS[name],
                     'count': int(record['numerator']),
                     'rate': float(record['rate'])})
    frame = pd.DataFrame(rows, columns=list(SOURCE_COVERAGE_FIGURE_COLUMNS))
    # The count is a gene count, so it is written as an integer (the panel rows
    # leave it empty); the rate keeps two decimals.
    frame['count'] = frame['count'].astype('Int64')
    return frame


def figure_rows(source, axis_sizes):
    """The complete plotted values of the two-panel rank-sum evidence figure.

    One panel per released axis. The gray bars are the frozen-axis candidates per
    ``rank_sum`` bin and are identical in every panel; ``selected`` says how many
    of that bin's candidates the panel's axis keeps (``sum(selected) == N``).
    """
    histogram = source['histogram']
    edges = np.append(histogram['bin_left'].to_numpy(dtype=float),
                      float(histogram['bin_right'].iloc[-1]))
    n_bins = len(histogram)
    rank_sum = source['level_two']['rank_sum']
    rows = [{'record': 'figure', 'title': FIGURE_TITLE}]
    for index, size in enumerate(axis_sizes):
        size = int(size)
        panel = chr(ord('a') + index)
        genes = source['axes'][size]
        bins = np.clip(np.searchsorted(edges, rank_sum[genes], side='right') - 1,
                       0, n_bins - 1)
        selected = np.bincount(bins, minlength=n_bins)
        rows.append({'record': 'panel', 'panel': panel, 'axis': str(size),
                     'panel_title': f'N = {size}'})
        for row, chosen in zip(histogram.itertuples(index=False), selected.tolist()):
            rows.append({'record': 'hist', 'panel': panel,
                         'bin_left': row.bin_left, 'bin_right': row.bin_right,
                         'bin_center': row.bin_center, 'count': int(row.count),
                         'selected': int(chosen)})
    return pd.DataFrame(rows, columns=list(FIGURE_COLUMNS))


# ---------------------------------------------------------------- release I/O

def input_entries(config, fusion):
    """Every pinned input with its role and digest (repo-relative paths only)."""
    entries = {}
    for role, path in (('frozen_axis', config.genes),
                       ('tenv_hvg_descriptive', config.tenv_hvg),
                       ('disease_fusion_catalog', config.fusion_root / 'catalog.json'),
                       ('disease_fusion_release_manifest',
                        config.fusion_root / 'manifest.json')):
        entries[repo_relative(path, config.root)] = {
            'role': role, 'sha256': sha256(path), 'bytes': path.stat().st_size}
    for category in RELEASE_CATEGORIES:
        record = fusion['diseases'].get(category)
        if record is None:
            continue  # an excluded disease contributes no payload
        path = config.fusion_root / record['deg2000_file']
        entries[repo_relative(path, config.root)] = {
            'role': f'{category}:deg2000', 'sha256': sha256(path),
            'bytes': path.stat().st_size}
    return dict(sorted(entries.items()))


def _axis_manifest(config, axis, tenv_frame, fusion, source, axis_sizes):
    axis_size = len(axis)
    fusion_result = source['level_two']
    table = source['table']
    active = list(source['active_categories'])
    n_diseases = len(active)
    composition = composition_rows(source, axis_sizes)
    coverage = coverage_rows(source, axis_sizes)
    axes_manifest = {}
    for size in axis_sizes:
        size = int(size)
        genes = source['axes'][size]
        fill = source['fill'][size]
        composition_counts = deg2000_composition(
            genes, _deg2000_masks(source), _deg2000_sizes(source),
            categories=RELEASE_CATEGORIES)
        counts = f2000_distribution(fusion_result['f2000'][genes],
                                    n_diseases=n_diseases)
        axes_manifest[str(size)] = {
            'n': int(len(genes)),
            'file': f'gene_axis_{size}.csv',
            'tenv_half': size // 2,
            'fill_half': size // 2,
            'second_half': {
                'from_deg2000': int(fill['from_deg2000']),
                'tenv_hvg_fallback': int(fill['fallback']),
                'n_deg2000_candidates': int(fill['n_deg2000_candidates']),
            },
            'rank_sum': source['intervals'][size],
            'tenv_rank': {'min': int(tenv_frame['rank'].to_numpy()[genes].min()),
                          'max': int(tenv_frame['rank'].to_numpy()[genes].max())},
            'f2000_distribution': {str(value): int(count)
                                   for value, count in enumerate(counts.tolist())},
            'deg2000_composition': {
                category: (None if composition_counts[category] is None else {
                    'n_covered': int(composition_counts[category]['n_covered']),
                    'n_deg2000': int(composition_counts[category]['n_deg2000']),
                    'share': float(composition_counts[category]['share'])})
                for category in RELEASE_CATEGORIES},
        }
    left, right = axis_sizes[0], axis_sizes[1]
    intersection = axis_intersection(source['axes'][left], source['axes'][right])
    intersection = {'sizes': [int(left), int(right)], **intersection}
    category_manifest = {}
    for category in RELEASE_CATEGORIES:
        record = fusion['diseases'].get(category)
        published = fusion['catalog']['diseases']
        catalog_record = next(item for item in published
                              if item['directory'] == category)
        entry = {
            'label': CATEGORY_LABELS[category],
            'product_dir': category,
            'published': bool(record is not None),
            'n_deg': (None if record is None else int(record['n_deg'])),
            'deg2000': (None if record is None else int(record['deg2000'])),
            'informative_source_count': int(catalog_record['informative_source_count']),
            'evidence_type': catalog_record.get('evidence_type'),
            'tier': catalog_record.get('tier'),
            'tiers_present': catalog_record.get('tiers_present'),
            'deg2000_file': (None if record is None else record['deg2000_file']),
            'deg2000_sha256': (None if record is None else record['deg2000_sha256']),
            'fusion_catalog_sha256': fusion['catalog'].get('catalog_sha256'),
            'sources': [dict(source_record) for source_record in
                        catalog_record.get('sources', [])],
            'excluded_sources': [dict(item) for item in
                                 catalog_record.get('excluded_sources', [])],
        }
        if record is None:
            entry['excluded_from_fusion'] = True
            entry['excluded_reason'] = (
                'the published disease has an empty DEG2000 (no informative '
                'source), so level-2 fusion and f2000 exclude it (ADR 0068)')
        category_manifest[category] = entry
    return {
        'schema': SCHEMA,
        'generated_by': 'stage2/01_bases/execution/build_gene_axis.py',
        'release': {
            'algorithm': ALGORITHM,
            'axis_sizes': [int(size) for size in axis_sizes],
            'fill_rule': FILL_RULE,
            'fill_fraction': 0.5,
            'f2000': F2000_RULE,
            'rank_sum': RANK_SUM_RULE,
            'fallback': FALLBACK_RULE,
            'fallback_second_half_slots': {
                str(int(size)): int(source['fill'][int(size)]['fallback'])
                for size in axis_sizes},
        },
        'environment': environment_facts(),
        'inputs': input_entries(config, fusion),
        'fusion': {
            'n_diseases': int(n_diseases),
            'diseases': active,
            'n_declared_diseases': len(CATEGORIES),
            'excluded_diseases': list(source['excluded_categories']),
            'rank_sum_sources': ('the active diseases only; tenv-HVG is not a '
                                 'rank_sum term (it supplies the axis first half)'),
            'rank_semantics': ('DEG2000 disease rank / DEG-union size (ADR 0071 '
                               'decision 5); a gene outside DEG2000_d contributes '
                               '(n_d + 1) / n_d'),
            'f2000_bounds': [0, int(n_diseases)],
            'n_candidates': int(len(table)),
        },
        'tenv_hvg': {
            'input': repo_relative(config.tenv_hvg, config.root),
            'n_ranked': int(axis_size),
            'source': 'published descriptive normalized-dispersion rank (never recomputed)',
            'role': ('first axis half (tenv_hvg_rank <= N/2) and the fallback order; '
                     'not a rank_sum term'),
        },
        'categories': category_manifest,
        'counts': {
            'n_axis': int(axis_size),
            'n_candidates': int(len(table)),
            'rank_sum_min': float(fusion_result['rank_sum'].min()),
            'rank_sum_median': float(np.median(fusion_result['rank_sum'])),
            'rank_sum_max': float(fusion_result['rank_sum'].max()),
            'f2000_distribution': {
                str(value): int(count) for value, count in enumerate(
                    f2000_distribution(fusion_result['f2000'],
                                       n_diseases=n_diseases).tolist())},
            'n_deg2000_candidates': int((fusion_result['f2000'] >= 1).sum()),
            'n_deg2000_overlap_rows': int(
                sum(int(record['deg2000']) for record in fusion['diseases'].values())
                - int((fusion_result['f2000'] >= 1).sum())),
        },
        'axes': axes_manifest,
        'intersection': intersection,
        'composition': {
            'file': COMPOSITION_NAME,
            'note': COMPOSITION_NOTE,
            'columns': list(COMPOSITION_COLUMNS),
            'rows': int(len(composition)),
        },
        'coverage': {
            'file': COVERAGE_NAME,
            'note': COVERAGE_NOTE,
            'columns': list(COVERAGE_COLUMNS),
            'rows': int(len(coverage)),
            'sources': {
                'tenv_hvg': 'published descriptive dispersion rank (never recomputed)',
                'categories': ("each disease's published DEG2000 set and its size; "
                               'numerator is |DEG2000_d ∩ G_N|'),
                'excluded': ('a disease with an empty DEG2000 has no coverage and '
                             'keeps an empty row (ADR 0068)'),
            },
            'figure': {
                'file': f'figure/data/{SOURCE_COVERAGE_FIGURE_NAME}.csv',
                'panels': [dict(panel) for panel in SOURCE_COVERAGE_FIGURE_PANELS],
                'value': ('count = |DEG2000_d ∩ G_N| (tenv-HVG: |top_N(tenv) ∩ G_N|); '
                          'rate = 100 * count / denominator'),
                'series': {f'G_{int(size)}': f'G_{int(size)} (N = {int(size):,})'
                           for size in axis_sizes},
                'sources': [SOURCE_COVERAGE_FIGURE_LABELS[name]
                            for name in ('tenv_hvg',) + tuple(source['active_categories'])],
            },
        },
        'diagnostics': {
            'inflection': {
                **{key: value for key, value in source['inflection'].items()},
                'role': 'diagnostic only; N is declared and never read off the curve',
                'stability': source['stability'],
            },
            'histogram': {
                'file': f'figure/data/{FIGURE_NAME}.csv',
                'low': 0.0, 'high': float(source['histogram_high']),
                'bins': HISTOGRAM_BINS,
                'panels': {str(int(size)): {
                    'panel': chr(ord('a') + index),
                    'panel_title': f'N = {int(size)}',
                    'axis': f'gene_axis_{int(size)}.csv',
                    'n_selected': int(size),
                    'n_candidates': int(len(table)),
                } for index, size in enumerate(axis_sizes)},
                'selection': ('the colored segment of each bar is the number of that '
                              "bin's candidates the panel's axis selects"),
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
                        'output/disease-fusion; never recomputes the tenv-HVG ranking '
                        'or the per-disease DEG2000 fusion'),
        'resource_contract': ('one heavy job at a time, 10 GiB RLIMIT_AS fallback, '
                              'OMP_NUM_THREADS=2; peak RSS recorded by the driver'),
    }


def product_paths(config, axis_sizes):
    """Repo-relative names of every product (manifest keys)."""
    output = repo_relative(config.output, config.root)
    figure = repo_relative(config.figure_data, config.root)
    source_coverage_figure = repo_relative(config.source_coverage_figure, config.root)
    names = {
        'candidate': f'{output}/{CANDIDATE_NAME}',
        'composition': f'{output}/{COMPOSITION_NAME}',
        'coverage': f'{output}/{COVERAGE_NAME}',
        'manifest': f'{output}/manifest.json',
        'figure': figure,
        'source_coverage_figure': source_coverage_figure,
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
    fusion = read_fusion(config, axis_size)
    source = compute(axis, tenv_frame, fusion, config.axis_sizes)

    outputs = {}
    candidate_path = payload_root / paths['candidate']
    write_frame(source['table'], candidate_path)
    outputs[paths['candidate']] = {'rows': int(len(source['table'])),
                                   'sha256': sha256(candidate_path),
                                   'columns': list(candidate_columns())}
    for size in config.axis_sizes:
        size = int(size)
        genes = source['axes'][size]
        frame = pd.DataFrame({'gene_index': genes,
                              'symbol': axis.symbol.to_numpy()[genes]})
        path = payload_root / paths[f'axis_{size}']
        write_frame(frame, path)
        outputs[paths[f'axis_{size}']] = {'rows': int(len(frame)),
                                          'sha256': sha256(path)}
    composition = composition_rows(source, config.axis_sizes)
    composition_path = payload_root / paths['composition']
    write_frame(composition, composition_path)
    outputs[paths['composition']] = {'rows': int(len(composition)),
                                     'sha256': sha256(composition_path),
                                     'columns': list(COMPOSITION_COLUMNS)}
    coverage = coverage_rows(source, config.axis_sizes)
    coverage_path = payload_root / paths['coverage']
    # Two decimals are the contract for the rate column: the cells are
    # percentages rounded by ``source_topn_coverage``, so printing 17 significant
    # digits would publish float representation noise (58.700000000000003).
    write_frame(coverage, coverage_path, float_format='%.2f')
    outputs[paths['coverage']] = {'rows': int(len(coverage)),
                                  'sha256': sha256(coverage_path),
                                  'columns': list(COVERAGE_COLUMNS)}
    figure = figure_rows(source, config.axis_sizes)
    figure_path = payload_root / paths['figure']
    write_frame(figure, figure_path)
    outputs[paths['figure']] = {'rows': int(len(figure)),
                                'sha256': sha256(figure_path)}
    source_coverage = source_coverage_figure_rows(source, config.axis_sizes)
    source_coverage_path = payload_root / paths['source_coverage_figure']
    # Two decimals are the contract for the rate column; the count column is
    # integer and is unaffected by the float format.
    write_frame(source_coverage, source_coverage_path, float_format='%.2f')
    outputs[paths['source_coverage_figure']] = {
        'rows': int(len(source_coverage)),
        'sha256': sha256(source_coverage_path),
        'columns': list(SOURCE_COVERAGE_FIGURE_COLUMNS)}

    manifest = _axis_manifest(config, axis, tenv_frame, fusion, source,
                              config.axis_sizes)
    manifest['outputs'] = dict(sorted(outputs.items()))
    manifest_path = payload_root / paths['manifest']
    write_json(manifest, manifest_path)
    validate_staged(config, source, payload_root, manifest, axis_size)
    return manifest


def validate_staged(config, source, payload_root, manifest, axis_size):
    """Structural acceptance of the candidate before it may replace the release."""
    table = source['table']
    active = list(source['active_categories'])
    n_diseases = len(active)
    if len(table) != axis_size:
        raise ValueError('candidate table does not have one row per frozen-axis gene')
    if list(table.columns) != list(candidate_columns()):
        raise ValueError('candidate table columns differ from the contract')
    if 'p2000' in table.columns or 'p5000' in table.columns:
        raise ValueError('the candidate table still carries a retired pN column')
    if table.gene_index.tolist() != list(range(axis_size)):
        raise ValueError('candidate table is not in frozen gene_index order')
    if int(table.f2000.min()) < 0 or int(table.f2000.max()) > n_diseases:
        raise ValueError(f'f2000 is outside 0..{n_diseases}')
    if not pd.api.types.is_integer_dtype(table.f2000):
        raise ValueError('f2000 must be an integer count')
    for category in RELEASE_CATEGORIES:
        column = table[f'{category}_rank']
        if category in active:
            if not column.notna().any():
                raise ValueError(f'{category}_rank is empty for a fused disease')
        elif column.notna().any():
            raise ValueError(f'the excluded {category}_rank column is not empty')
    for size in config.axis_sizes:
        size = int(size)
        genes = source['axes'][size]
        if len(genes) != size:
            raise ValueError(f'G_{size} does not hold exactly {size} genes')
        if not np.array_equal(genes, np.sort(genes)):
            raise ValueError(f'G_{size} is not ordered by ascending gene_index')
        if len(np.unique(genes)) != size or genes.min() < 0 or genes.max() >= axis_size:
            raise ValueError(f'G_{size} is not a distinct subset of the frozen axis')
        expected = axis_fill(np.asarray(source['table'].tenv_hvg_rank, dtype=np.int64),
                             table.f2000.to_numpy(),
                             table.rank_sum.to_numpy(),
                             n=size, n_diseases=n_diseases)
        if not np.array_equal(expected['axis'], genes):
            raise ValueError(f'G_{size} does not follow the declared fill rule')
        entry = manifest['axes'][str(size)]
        if entry['second_half']['from_deg2000'] != expected['from_deg2000']:
            raise ValueError(f'G_{size} reports the wrong DEG2000 second-half count')
        if entry['second_half']['tenv_hvg_fallback'] != expected['fallback']:
            raise ValueError(f'G_{size} reports the wrong tenv-HVG fallback count')
        if (entry['second_half']['from_deg2000']
                + entry['second_half']['tenv_hvg_fallback'] != size // 2):
            raise ValueError(f'G_{size} second-half counts do not fill 0.5*N')
        if manifest['release']['fallback_second_half_slots'][str(size)] != expected['fallback']:
            raise ValueError(f'G_{size} fallback count differs from the release record')
        # A gene with f2000 = 0 may only be in the tenv half or the fallback.
        chosen = np.setdiff1d(genes, np.flatnonzero(
            np.asarray(source['table'].tenv_hvg_rank, dtype=np.int64) <= size // 2))
        if (table.f2000.to_numpy()[chosen] == 0).sum() != expected['fallback']:
            raise ValueError(f'G_{size} admits a non-DEG gene outside the fallback')
    paths = product_paths(config, config.axis_sizes)
    coverage = pd.read_csv(payload_root / paths['coverage'])
    if list(coverage.columns) != list(COVERAGE_COLUMNS):
        raise ValueError('coverage table columns differ from the contract')
    if list(dict.fromkeys(coverage['N'].tolist())) != [
            int(size) for size in config.axis_sizes]:
        raise ValueError('coverage table rows are not the released axis sizes in order')
    if len(coverage) != len(config.axis_sizes) * (1 + len(RELEASE_CATEGORIES)):
        raise ValueError('coverage table does not hold one row per axis and group')
    fused = coverage.dropna(subset=['numerator'])
    if len(fused):
        if (fused['numerator'] > fused['denominator']).any():
            raise ValueError('a coverage numerator exceeds its denominator')
        if not fused['rate'].between(0.0, 100.0).all():
            raise ValueError('a coverage rate is outside 0..100 percent')
        expected_rate = (100.0 * fused['numerator'].astype(float)
                         / fused['denominator'].astype(float)).round(2)
        if (np.abs(fused['rate'].astype(float) - expected_rate) > 1e-9).any():
            raise ValueError('a coverage rate is not numerator / denominator')
    tenv_rows = coverage[coverage.group == 'tenv_hvg']
    if (tenv_rows['denominator'].astype(int) != tenv_rows['N'].astype(int)).any():
        raise ValueError('the tenv-HVG coverage denominator must be N')
    for category in RELEASE_CATEGORIES:
        rows = coverage[coverage.group == category]
        if len(rows) != len(config.axis_sizes):
            raise ValueError(f'coverage lacks a row per axis for {category}')
        if category in active:
            if rows['numerator'].isna().any():
                raise ValueError(f'coverage for the fused {category} is empty')
            sizes = _deg2000_sizes(source)
            denominator = sizes[RELEASE_CATEGORIES.index(category)]
            if (rows['denominator'].astype(int) != denominator).any():
                raise ValueError(f'coverage for {category} does not use |DEG2000_d|')
        elif rows['numerator'].notna().any():
            raise ValueError(f'coverage for the excluded {category} is not empty')
    composition = pd.read_csv(payload_root / paths['composition'])
    if list(composition.columns) != list(COMPOSITION_COLUMNS):
        raise ValueError('composition table columns differ from the contract')
    if set(composition.record) != {'axis_fill', 'deg2000', 'f2000'}:
        raise ValueError('composition table does not carry the contracted records')
    for size in config.axis_sizes:
        rows = composition[composition.axis.astype(str) == str(int(size))]
        fill = rows[rows.record == 'axis_fill'].set_index('group')
        second = source['fill'][int(size)]
        if int(fill.loc['deg2000_second_half', 'n_genes']) != second['from_deg2000']:
            raise ValueError('composition reports the wrong DEG2000 second-half count')
        if int(fill.loc['tenv_hvg_fallback', 'n_genes']) != second['fallback']:
            raise ValueError('composition reports the wrong tenv-HVG fallback count')
        if int(fill['n_genes'].sum()) != int(size):
            raise ValueError('the axis_fill composition does not cover the axis')
        counts = rows[rows.record == 'f2000']
        if int(counts['n_genes'].sum()) != int(size):
            raise ValueError('the f2000 composition does not cover the axis')
        deg2000 = rows[rows.record == 'deg2000']
        if len(deg2000) != len(RELEASE_CATEGORIES):
            raise ValueError('the composition lacks a DEG2000 row per released disease')
    figure = pd.read_csv(payload_root / paths['figure'])
    if list(figure.columns) != list(FIGURE_COLUMNS):
        raise ValueError('figure plotted-values columns differ from the contract')
    panels = figure[figure.record == 'panel']
    histogram = figure[figure.record == 'hist']
    sizes = {row.panel: int(row.axis) for row in panels.itertuples(index=False)}
    if sorted(sizes.values()) != [int(size) for size in config.axis_sizes]:
        raise ValueError('figure panels do not match the released axis sizes')
    if set(histogram.panel) != set(sizes):
        raise ValueError('a figure panel has no plotted candidate histogram')
    n_candidates = int(histogram[histogram.panel == panels.panel.iloc[0]]['count'].sum())
    if n_candidates != axis_size:
        raise ValueError('the figure histogram does not cover every candidate')
    for panel, size in sizes.items():
        rows = histogram[histogram.panel == panel]
        if len(rows) != HISTOGRAM_BINS:
            raise ValueError('a figure panel does not hold the declared bin count')
        if rows.selected.sum() != size:
            raise ValueError(f'figure panel {panel} does not select exactly N genes')
        if (rows.selected > rows['count']).any():
            raise ValueError('a figure bin selects more genes than it holds')
    coverage_figure = pd.read_csv(payload_root / paths['source_coverage_figure'])
    if list(coverage_figure.columns) != list(SOURCE_COVERAGE_FIGURE_COLUMNS):
        raise ValueError('source-coverage figure columns differ from the contract')
    panels = list(coverage_figure[coverage_figure.record == 'panel']
                  .itertuples(index=False))
    if len(panels) != len(SOURCE_COVERAGE_FIGURE_PANELS):
        raise ValueError('the source-coverage figure does not declare every panel')
    for expected, actual in zip(SOURCE_COVERAGE_FIGURE_PANELS, panels):
        if (actual.panel != expected['panel']
                or actual.panel_title != expected['title']
                or actual.y_label != expected['y_label']
                or actual.value != expected['value']):
            raise ValueError('a source-coverage panel differs from the contract')
    bars = coverage_figure[coverage_figure.record == 'bar']
    sources = ('tenv_hvg',) + tuple(active)
    axis_series = {f'G_{int(size)}': int(size) for size in config.axis_sizes}
    expected_bars = {(series, name) for series in axis_series for name in sources}
    if (set(zip(bars.series, bars.source)) != expected_bars
            or len(bars) != len(expected_bars)):
        raise ValueError('the source-coverage figure does not plot every bar once')
    if bars['count'].isna().any() or (bars['count'] < 0).any():
        raise ValueError('a source-coverage count is negative or missing')
    if (bars['count'] % 1 != 0).any():
        raise ValueError('a source-coverage count is not an integer')
    if bars['rate'].isna().any() or not bars['rate'].between(0.0, 100.0).all():
        raise ValueError('a source-coverage rate is outside 0..100 percent')
    published = {(int(record['N']), record['group']): record
                 for record in coverage.to_dict('records')}
    for row in bars.itertuples(index=False):
        record = published[(int(row.series.split('_')[1]), row.source)]
        if int(row.count) != int(record['numerator']):
            raise ValueError('a source-coverage bar differs from the published numerator')
        if abs(float(row.rate) - float(record['rate'])) > 1e-9:
            raise ValueError('a source-coverage bar differs from the published rate')
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
        # Candidate-then-replace: only this ticket's products are swapped in.
        output_root = payload / repo_relative(config.output, config.root)
        config.output.mkdir(parents=True, exist_ok=True)
        for key in ('candidate', 'composition', 'coverage'):
            os.replace(output_root / Path(paths[key]).name,
                       config.output / Path(paths[key]).name)
        for size in config.axis_sizes:
            name = Path(paths[f'axis_{int(size)}']).name
            os.replace(output_root / name, config.output / name)
        config.figure_data.parent.mkdir(parents=True, exist_ok=True)
        os.replace(payload / paths['figure'], config.figure_data)
        config.source_coverage_figure.parent.mkdir(parents=True, exist_ok=True)
        os.replace(payload / paths['source_coverage_figure'],
                   config.source_coverage_figure)
        # The manifest is replaced last: it only announces products already in place.
        os.replace(output_root / 'manifest.json', config.output / 'manifest.json')
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    return summary(config, manifest)


def verify(config):
    """Rebuild from the pinned inputs and compare the manifest and every digest.

    A release whose manifest predates the v2 ``f2000`` algorithm (its
    ``release.algorithm`` is absent or different, e.g. the committed v1 release
    with ``p2000``/``p5000`` candidate columns) is refused with a
    self-describing error before the rebuild: the products are stale by design
    until #208 republishes them, and an ambiguous digest mismatch would hide
    that. The verify path is never weakened: once the manifest carries the v2
    algorithm, the full replay comparison runs unchanged.
    """
    paths = product_paths(config, config.axis_sizes)
    accepted_path = config.output / 'manifest.json'
    if not accepted_path.is_file():
        raise ValueError(f'no accepted gene-axis manifest at {accepted_path}')
    accepted = json.loads(accepted_path.read_text())
    if accepted.get('schema') != SCHEMA:
        raise ValueError('unknown gene-axis release schema')
    algorithm = (accepted.get('release') or {}).get('algorithm')
    if algorithm != ALGORITHM:
        raise ValueError(
            'the committed gene-axis products predate the v2 f2000 algorithm: '
            f'release.algorithm is {algorithm!r}, expected {ALGORITHM!r}. The '
            'committed release was built under the retired v1 pN fill rule (its '
            'candidate table carries p2000/p5000 and its source-coverage table '
            'uses the measured-gene schema), so a replay cannot match by design. '
            '#208 must republish output/gene-axis/ with build_gene_axis.py before '
            '--verify can pass; this ticket does not touch the accepted products.')
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
        'result': 'PASS', 'schema': SCHEMA, 'algorithm': ALGORITHM,
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
        'algorithm': manifest['release']['algorithm'],
        'n_candidates': manifest['counts']['n_candidates'],
        'axes': manifest['release']['axis_sizes'],
        'intersection': manifest['intersection']['n_intersection'],
        'n_diseases': manifest['fusion']['n_diseases'],
        'excluded_diseases': manifest['fusion']['excluded_diseases'],
        'fallback_second_half_slots': manifest['release']['fallback_second_half_slots'],
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
    parser.add_argument('--fusion-root', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--figure-data', type=Path)
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--record', type=Path)
    args = parser.parse_args()
    config = default_config(args.root, args.output, args.figure_data)
    overrides = {'genes': args.genes, 'tenv_hvg': args.tenv_hvg,
                 'fusion_root': args.fusion_root}
    config = replace(config, **{key: value for key, value in overrides.items() if value})
    result = verify(config) if args.verify else publish(config)
    text = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.record:
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
