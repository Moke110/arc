"""Read accepted inputs and recompute HVG-DEG axis selection."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from cpst_util.derivation.disease_deg import (DEG1000_COLUMNS, DEG1000_LIMIT, DiseaseFusion,
                               disease_eligibility)
from axis_selection_kernel import (
    CATEGORIES, CATEGORY_LABELS, RELEASE_CATEGORIES, axis_fill, candidate_table,
    cross_disease_rank_sum, deg1000_frame, detect_inflection, inflection_stability,
    rank_sum_histogram,
)
from cpst_util.derivation.release import sha256

HISTOGRAM_BINS = 70
INFLECTION_MIN_N = 100


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
    """Read tier-first DEG1000 through the digest-checked fusion reader.

    Only diseases whose complete axis-resolved DEG union passes the 100-gene
    gate contribute. Every outside-cut gene receives n_d+1, including an actual
    DEG ranked 1001 or later; n_d remains the complete union size.
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
    if int(catalog.get('limit', -1)) != int(DEG1000_LIMIT):
        raise ValueError(f'the disease-fusion DEG1000 limit must be {DEG1000_LIMIT}')
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
            '#208 must publish their DEG1000 products before the axes can be rebuilt')
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
        deg1000 = int(record['deg1000'])
        admission = disease_eligibility(n_deg)
        if not admission['eligible'] or n_deg <= 0 or deg1000 <= 0:
            excluded.append({'directory': category, 'label': CATEGORY_LABELS[category],
                             'n_deg': n_deg, 'deg1000': deg1000,
                             **admission})
            continue
        if deg1000 != min(n_deg, DEG1000_LIMIT):
            raise ValueError(
                f'{category}: published deg1000 must be min(n_deg, {DEG1000_LIMIT})')
        name = record['deg1000_file']
        declared = declared_outputs.get(name)
        digest = declared.get('sha256') if isinstance(declared, dict) else declared
        if not isinstance(digest, str) or len(digest) != 64:
            raise ValueError(f'{category}: {name} is not a declared fusion output')
        rows = reader.deg1000(category)  # digest-checked against the fusion manifest
        if list(rows.columns) != list(DEG1000_COLUMNS):
            raise ValueError(f'{category}: {name} has unexpected columns')
        if len(rows) != deg1000:
            raise ValueError(f'{category}: {name} holds {len(rows)} rows, expected {deg1000}')
        genes = np.asarray(rows.gene_index, dtype=np.int64)
        ranks = np.asarray(rows.disease_rank, dtype=np.int64)
        if genes.min() < 0 or genes.max() >= int(axis_size):
            raise ValueError(f'{category}: {name} holds a gene_index outside the frozen axis')
        if len(np.unique(genes)) != len(genes):
            raise ValueError(f'{category}: {name} repeats a gene_index')
        if not np.array_equal(np.sort(ranks),
                              np.arange(1, deg1000 + 1, dtype=np.int64)):
            raise ValueError(f'{category}: {name} ranks must be exactly 1..{deg1000}')
        frame = np.full(int(axis_size), n_deg + 1, dtype=np.int64)
        frame[genes] = ranks
        frame = deg1000_frame(frame, n_deg, int(axis_size), category)
        diseases[category] = {
            'record': record, 'n_deg': n_deg, 'deg1000': deg1000,
            'rank': frame, 'mask': frame <= deg1000,
            'deg1000_file': name, 'deg1000_sha256': digest,
        }
    if not diseases:
        raise ValueError(
            'no released disease has a published DEG1000; the axes cannot be built '
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
        result = axis_fill(tenv_rank, level_two['f1000'], level_two['rank_sum'],
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
