"""Two-level rank-sum fusion and the released Stage 2 axes (ADR 0066).

This module owns every scientific definition of the released gene axis: the
winning-tier selection inside a disease category, the level-1 category rank, the
level-2 rank sum over the **fused** sources, the ``pN`` cross-source replication
count, the ``0.5*N + 0.5*N`` axis fill rule, the source top-``N`` coverage table,
the rank-sum histogram and the diagnostic Kneedle + second-difference inflection.
``execution/build_gene_axis.py`` reads the pinned inputs and publishes the
products; ``execution/tests/test_gene_axis.py`` exercises these definitions on
hermetic fixtures.

It is a pure function module (ADR 0062): no file is opened here except through
:class:`GeneAxis`, the public reader for the published ``output/gene-axis/``
release. Nothing here chooses an axis size: ``N`` is declared (2,000 and 5,000)
and the inflection diagnostic never selects it.

Vocabulary follows ADR 0066 and ``GLOSSARY.md`` -> *Stage 2 gene axis*:
``Source rank``, ``Class source``, ``Rank-sum fusion``, ``pN``, ``Axis fill
rule``, ``Inflection point``.

Index convention: every per-gene vector is indexed by the frozen ``gene_index``
(``0..axis_size-1``). A released axis is an ascending ordering of original
``gene_index`` values, never a renumbered ``0..N-1`` position vector.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

# The #186 contract: ``bases.gene_axis`` re-exports the DEG/rank functions now
# owned by ``bases.deg_rank`` (one live implementation, no second copy).
from .deg_rank import (
    ALIGNMENT_CATEGORIES,
    AXIS_SIZE,
    PADJ_FLOOR,
    PADJ_THRESHOLD,
    SymbolResolver,
    complete_rank,
    deg_effect,
    normalized_rank,
    per_gene_deg,
    rank_scores,
)

SCHEMA = 'stage2-gene-axis-release-v1'

# The five declared disease categories, in ADR 0066 roster order minus ``af``
# (removed by ADR 0069: its only dataset carried no ordering evidence). The
# product directory is the catalog directory under ``output/disease-deg/``; it is
# also the prefix of the candidate table's ``<category>_rank`` column.
CATEGORIES = ('nhf', 'icm', 'genetic_cardiomyopathy',
              'cardiotox_TOP2', 'cardiotox_other')
CATEGORY_LABELS = {'nhf': 'nHF', 'icm': 'ICM',
                   'genetic_cardiomyopathy': 'genetic_cardiomyopathy',
                   'cardiotox_TOP2': 'cardiotox_TOP2',
                   'cardiotox_other': 'cardiotox_other'}
CATEGORY_RANK_COLUMNS = tuple(f'{category}_rank' for category in CATEGORIES)

#: Declared released axis sizes (ADR 0066). N is chosen, never read off a curve.
AXIS_SIZES = (2000, 5000)
#: Half of an axis comes from tenv-HVG, half from cross-source replication.
FILL_FRACTION = 0.5
#: The declared maximum: five category class sources plus the published
#: tenv-HVG ranking. A category whose records are all non-informative is dropped
#: from the fusion (ADR 0068), so the *actual* source count is
#: :func:`n_fusion_sources` over the informative categories and is published in
#: the release manifest; ``pN`` is bounded by that count, not by this maximum.
N_FUSION_SOURCES = len(CATEGORIES) + 1


def n_fusion_sources(categories=None):
    """The number of sources actually fused: the given categories plus tenv-HVG."""
    categories = CATEGORIES if categories is None else tuple(categories)
    if not categories:
        return 1
    if len(set(categories)) != len(categories):
        raise ValueError('a fusion source list cannot repeat a category')
    unknown = set(categories) - set(CATEGORIES)
    if unknown:
        raise ValueError(f'unknown fusion categories: {sorted(unknown)}')
    return len(categories) + 1

def candidate_columns(sizes=AXIS_SIZES):
    """Contracted candidate-table columns for the requested axis sizes."""
    return (('symbol', 'gene_index', 'tenv_hvg_rank') + CATEGORY_RANK_COLUMNS
            + ('rank_sum', 'n_sources_measured')
            + tuple(f'p{int(size)}' for size in sizes))


CANDIDATE_COLUMNS = candidate_columns()
AXIS_COLUMNS = ('gene_index', 'symbol')
#: The source top-``N`` coverage table: one row per released axis, one column per
#: source (``tenv_hvg`` then the five categories in roster order). The header keeps
#: the literal ``N`` so a single table can carry every released axis.
COVERAGE_COLUMNS = ("N", "%tenv_hvg_gN") + tuple(
    f"%{category}_gN" for category in CATEGORIES)

#: The integer companion of :data:`COVERAGE_COLUMNS`: one row per released axis,
#: one column per source, holding how many of that source's top-``N`` genes the
#: axis contains. The columns are bare source names because the values are counts,
#: not percentages.
COUNT_COLUMNS = ("N", "tenv_hvg") + CATEGORIES

FILL_RULE = (
    'G_N = the top 0.5*N genes of tenv-HVG (tenv_hvg_rank <= N/2), then 0.5*N '
    'more from the remaining candidates ordered by pN descending, rank_sum '
    'ascending, gene_index ascending; the union is published ordered by '
    'ascending gene_index and is never renumbered.')


def winning_tier(catalog):
    """Select the datasets a category fuses, and report the rule that selected them.

    ADR 0066 *Source eligibility*: the highest available priority tier supplies
    the source and **all independent datasets at that tier** are kept. ADR 0066
    *Isogenic exception*: an in-vitro isogenic model is kept as its own evidence
    type even when a higher-tier human cohort exists for the same category, so a
    record whose ``isogenic_exception.applies`` is true is retained alongside the
    winning-tier records instead of being dropped by the tier rule.

    Returns ``{'winning_tier', 'datasets', 'records', 'retained_by_exception',
    'excluded_tiers', 'rule'}``. ``datasets`` is the catalog order, so the fused
    dataset order is the producer's declared order, not a re-sort.
    ``excluded_tiers`` names only tiers from which **no** record was retained: a
    tier whose records are all retained by the isogenic exception is not an
    excluded tier, and the manifest must not call it one (ADR 0068).
    """
    records = list(catalog.get('datasets') or ())
    if not records:
        raise ValueError('a disease category catalog needs at least one dataset record')
    tiers = []
    for record in records:
        tier = record.get('tier')
        if tier is None or int(tier) <= 0:
            raise ValueError(
                f'{record.get("dataset_id", "?")}: dataset record has no priority tier')
        tiers.append(int(tier))
    winning = min(tiers)
    retained = []
    retained_tiers = set()
    by_exception = []
    for record, tier in zip(records, tiers):
        isogenic = (record.get('isogenic_exception') or {}).get('applies') is True
        if tier == winning:
            retained.append(record)
            retained_tiers.add(tier)
        elif isogenic:
            retained.append(record)
            retained_tiers.add(tier)
            by_exception.append(record['dataset_id'])
    return {
        'winning_tier': winning,
        'datasets': [record['dataset_id'] for record in retained],
        'records': retained,
        'retained_by_exception': by_exception,
        'excluded_tiers': sorted(set(tiers) - retained_tiers),
        'rule': ('minimum catalog tier, plus records whose isogenic_exception.applies '
                 'is true (ADR 0066 Isogenic exception); excluded_tiers lists only '
                 'tiers with no retained record'),
    }


def category_rank(ranks, n_measured, gene_index=None, axis_size=None):
    """Level 1: a category's source rank from its datasets' within-dataset ranks.

    ADR 0066 *Candidate table and the two-level rank-sum fusion*::

        level-1 sum = sum(rank_norm)  over the datasets at the winning tier
        category rank = position in ascending level-1 sum, ties by gene_index

    ``rank_norm`` is recomputed from the raw rank and the dataset's
    ``n_measured`` with :func:`bases.deg_rank.normalized_rank` (every gene is
    ``rank/n_measured``; an unmeasured gene's rank is ``n_measured + 1``, so its
    value is ``(n_measured + 1)/n_measured``) so a category never trusts a
    hand-edited normalized value. Every dataset rank vector covers the frozen
    axis (``gene_index`` order) and is a complete ``1..n_measured+1`` assignment.

    Returns ``{'rank', 'sum', 'measured', 'n_datasets'}``: a 1-based
    ``1..axis_size`` permutation, the level-1 ``sum(rank_norm)``, and a boolean
    "this source measures the gene" mask (measured by at least one dataset).
    """
    ranks = [np.asarray(rank, dtype=np.int64) for rank in ranks]
    n_measured = [int(value) for value in n_measured]
    if not ranks:
        raise ValueError('a category needs at least one dataset rank')
    if len(ranks) != len(n_measured):
        raise ValueError('every dataset rank needs its n_measured')
    axis_size = int(axis_size if axis_size is not None else len(ranks[0]))
    gene_index = (np.arange(axis_size, dtype=np.int64) if gene_index is None
                  else np.asarray(gene_index, dtype=np.int64))
    if gene_index.shape != (axis_size,):
        raise ValueError('gene_index must cover the frozen axis')
    total = np.zeros(axis_size, dtype=np.float64)
    measured = np.zeros(axis_size, dtype=bool)
    for rank, n in zip(ranks, n_measured):
        if rank.shape != (axis_size,):
            raise ValueError('a dataset rank vector must cover the frozen axis')
        if n <= 0:
            raise ValueError('a dataset needs at least one measured gene')
        if rank.min() < 1 or rank.max() > n + 1:
            raise ValueError('a dataset rank must be in 1..n_measured+1')
        # A complete frozen-axis rank frame: measured genes get exactly 1..n and
        # every unmeasured gene gets n+1 (ADR 0066 / bases.deg_rank contract).
        if int((rank <= n).sum()) != n or np.any(rank > n + 1):
            raise ValueError('a dataset rank frame is not a complete 1..n_measured+1 assignment')
        if not np.array_equal(np.sort(rank[rank <= n]), np.arange(1, n + 1)):
            raise ValueError('measured ranks must be exactly 1..n_measured')
        measured |= rank <= n
        total += normalized_rank(rank, n)
    order = np.lexsort((gene_index, total))
    rank = np.empty(axis_size, dtype=np.int64)
    rank[order] = np.arange(1, axis_size + 1, dtype=np.int64)
    return {'rank': rank, 'sum': total, 'measured': measured, 'n_datasets': len(ranks)}


def _complete_rank(rank, axis_size, label):
    """Validate one source's frozen-axis rank vector is a complete 1..axis_size permutation."""
    rank = np.asarray(rank, dtype=np.int64)
    if rank.shape != (axis_size,):
        raise ValueError(f'{label} rank must cover the frozen axis')
    if not np.array_equal(np.sort(rank), np.arange(1, axis_size + 1, dtype=np.int64)):
        raise ValueError(f'{label} rank must be a complete 1..{axis_size} permutation')
    return rank


def fusion_rank_sum(category_ranks, categories_measured, tenv_rank,
                    categories=None, gene_index=None, axis_size=None,
                    sizes=AXIS_SIZES):
    """Level 2: ``rank_sum`` over the fused sources (ADR 0066, ADR 0068).

    ADR 0066 / spec D13: ``rank_sum = sum(rank_norm)`` over the category class
    sources plus the published tenv-HVG ranking, each normalized by the
    frozen-axis size; **ascending is better** and ties break on ``gene_index``.
    ``pN`` counts how many of the *fused* sources place the gene in their top-``N``
    set; ``n_sources_measured`` counts the fused sources that measure it
    (tenv-HVG ranks every frozen-axis gene, so it always contributes one).

    ``categories`` is the active (informative) category subset in :data:`CATEGORIES`
    order; a category whose every dataset is non-informative is passed by
    omission and is excluded from ``rank_sum``, ``pN`` and ``n_sources_measured``
    (ADR 0068). ``category_ranks`` are the level-1 category ranks (one complete
    permutation per active category), ``categories_measured`` the corresponding
    "source measures the gene" masks, ``tenv_rank`` the published tenv-HVG
    dispersion rank. Returns ``{'rank_sum', 'n_sources_measured', 'pn',
    'n_sources', 'categories'}`` where ``pn`` maps each requested ``N`` to its
    integer vector.
    """
    categories = tuple(CATEGORIES if categories is None else categories)
    n_sources = n_fusion_sources(categories)
    if len(category_ranks) != len(categories):
        raise ValueError(
            f'a level-2 fusion needs one rank per active category ({len(categories)})')
    if len(categories_measured) != len(categories):
        raise ValueError('every category rank needs its measured mask')
    tenv_rank = np.asarray(tenv_rank, dtype=np.int64)
    axis_size = int(axis_size if axis_size is not None else tenv_rank.size)
    gene_index = (np.arange(axis_size, dtype=np.int64) if gene_index is None
                  else np.asarray(gene_index, dtype=np.int64))
    if gene_index.shape != (axis_size,):
        raise ValueError('gene_index must cover the frozen axis')
    tenv_rank = _complete_rank(tenv_rank, axis_size, 'tenv-HVG')
    sources = [_complete_rank(rank, axis_size, category)
               for rank, category in zip(category_ranks, categories)]
    sources.append(tenv_rank)
    rank_sum = np.zeros(axis_size, dtype=np.float64)
    for rank in sources:
        rank_sum += rank / float(axis_size)
    measured = [np.asarray(mask, dtype=bool) for mask in categories_measured]
    for mask in measured:
        if mask.shape != (axis_size,):
            raise ValueError('a category measured mask must cover the frozen axis')
    n_sources_measured = 1 + np.sum(measured, axis=0).astype(np.int64)
    pn = {}
    for size in sizes:
        size = int(size)
        if size <= 0 or size > axis_size:
            raise ValueError('a pN size must be within the frozen axis')
        pn[size] = np.sum([rank <= size for rank in sources], axis=0).astype(np.int64)
    return {'rank_sum': rank_sum, 'n_sources_measured': n_sources_measured,
            'pn': pn, 'n_sources': n_sources, 'categories': categories}


def candidate_table(axis, tenv_rank, category_ranks, categories_measured,
                    categories=None, gene_index=None, sizes=AXIS_SIZES):
    """The master candidate table: one row per frozen-axis gene.

    ADR 0066: ``symbol``, ``gene_index``, ``tenv_hvg_rank``, one
    ``<category>_rank`` per category, ``rank_sum``, ``n_sources_measured``,
    ``p2000``, ``p5000``. Rows are in frozen-axis order; every vector is indexed
    by ``gene_index`` so no consumer can read a row position as a gene index.

    ``category_ranks`` is either a mapping ``{category: rank}`` or a sequence
    aligned with ``categories``. A category excluded from the fusion (ADR 0068:
    every dataset non-informative) keeps its contracted column but the column is
    empty (``NaN``) — it has no level-1 rank — and it does not contribute to
    ``rank_sum``/``pN``/``n_sources_measured``.
    """
    axis = pd.DataFrame(axis)
    required = {'gene_index', 'symbol'}
    if not required <= set(axis.columns):
        raise ValueError(f'axis table requires {sorted(required)}')
    axis_size = len(axis)
    gene_index = (np.arange(axis_size, dtype=np.int64) if gene_index is None
                  else np.asarray(gene_index, dtype=np.int64))
    if gene_index.shape != (axis_size,):
        raise ValueError('gene_index must cover the frozen axis')
    if not np.array_equal(np.asarray(axis.gene_index, dtype=np.int64), gene_index):
        raise ValueError('axis table must be in frozen gene_index order')
    if isinstance(category_ranks, dict):
        ranks_by_category = {str(key): np.asarray(value, dtype=np.int64)
                             for key, value in category_ranks.items()}
        active = tuple(CATEGORIES if categories is None else categories)
        if set(ranks_by_category) != set(active):
            raise ValueError('category rank mapping does not match the active categories')
        rank_list = [ranks_by_category[category] for category in active]
        measured_list = list(categories_measured)
    else:
        active = tuple(CATEGORIES if categories is None else categories)
        if len(category_ranks) != len(active):
            raise ValueError('a candidate table needs one rank per active category')
        ranks_by_category = {category: np.asarray(rank, dtype=np.int64)
                             for category, rank in zip(active, category_ranks)}
        rank_list = list(category_ranks)
        measured_list = list(categories_measured)
    fusion = fusion_rank_sum(rank_list, measured_list, tenv_rank,
                             categories=active, gene_index=gene_index,
                             axis_size=axis_size, sizes=sizes)
    columns = {'symbol': axis.symbol.to_numpy(), 'gene_index': gene_index,
               'tenv_hvg_rank': np.asarray(tenv_rank, dtype=np.int64)}
    for category in CATEGORIES:
        if category in ranks_by_category:
            columns[f'{category}_rank'] = ranks_by_category[category]
        else:
            columns[f'{category}_rank'] = np.full(axis_size, np.nan)
    columns['rank_sum'] = fusion['rank_sum']
    columns['n_sources_measured'] = fusion['n_sources_measured']
    for size in sizes:
        columns[f'p{int(size)}'] = fusion['pn'][int(size)]
    table = pd.DataFrame(columns)
    return table.loc[:, list(candidate_columns(sizes))]


def axis_fill(tenv_rank, p_n, rank_sum, n, gene_index=None,
              n_sources=N_FUSION_SOURCES):
    """One released axis: the declared ``0.5*N + 0.5*N`` fill rule (ADR 0066).

    ``G_N`` takes the top ``0.5*N`` genes of tenv-HVG (``tenv_hvg_rank <= N/2``),
    then ``0.5*N`` more from the **remaining** candidates ordered by ``pN``
    descending, ``rank_sum`` ascending, ``gene_index`` ascending. The result is
    returned as ascending original ``gene_index`` integers: an axis is an
    ordered subset of the frozen space, never renumbered to ``0..N-1``.

    ``p_n`` is the gene's ``pN`` count for this axis size (``p2000`` for N=2000,
    ``p5000`` for N=5000); the differing ``pN`` vectors are the main reason the
    two axes are not nested.
    """
    tenv_rank = np.asarray(tenv_rank, dtype=np.int64)
    axis_size = tenv_rank.size
    tenv_rank = _complete_rank(tenv_rank, axis_size, 'tenv-HVG')
    n = int(n)
    if n < 2 or n % 2:
        raise ValueError('an axis size must be an even integer >= 2')
    if n > axis_size:
        raise ValueError('an axis size cannot exceed the frozen axis')
    half = n // 2
    p_n = np.asarray(p_n, dtype=np.int64)
    rank_sum = np.asarray(rank_sum, dtype=np.float64)
    if p_n.shape != (axis_size,) or rank_sum.shape != (axis_size,):
        raise ValueError('pN and rank_sum must cover the frozen axis')
    n_sources = int(n_sources)
    if n_sources < 1:
        raise ValueError('pN needs at least the tenv-HVG source')
    if p_n.min() < 0 or p_n.max() > n_sources:
        raise ValueError(f'pN must be in 0..{n_sources}')
    if not np.isfinite(rank_sum).all():
        raise ValueError('rank_sum must be finite')
    gene_index = (np.arange(axis_size, dtype=np.int64) if gene_index is None
                  else np.asarray(gene_index, dtype=np.int64))
    if gene_index.shape != (axis_size,):
        raise ValueError('gene_index must cover the frozen axis')
    tenv_half = np.flatnonzero(tenv_rank <= half)
    if len(tenv_half) != half:
        raise ValueError('tenv-HVG rank is not a complete permutation')
    remaining = np.flatnonzero(tenv_rank > half)
    # lexsort's last key is primary: pN descending, rank_sum ascending, index.
    order = np.lexsort((gene_index[remaining], rank_sum[remaining],
                        -p_n[remaining]))
    chosen = remaining[order[:half]]
    axis = np.sort(np.concatenate([tenv_half, chosen]).astype(np.int64))
    if len(np.unique(axis)) != n or len(axis) != n:
        raise ValueError('axis fill did not select N distinct frozen-axis genes')
    return axis


def axis_intersection(left, right):
    """Sizes and Jaccard of two released axes; the manifest reports the intersection."""
    left = set(np.asarray(left, dtype=np.int64).tolist())
    right = set(np.asarray(right, dtype=np.int64).tolist())
    intersection = left & right
    union = left | right
    smaller = min(len(left), len(right))
    return {
        'n_left': len(left), 'n_right': len(right),
        'n_intersection': len(intersection), 'n_union': len(union),
        'jaccard': len(intersection) / len(union) if union else 0.0,
        'share_of_smaller': len(intersection) / smaller if smaller else 0.0,
    }


def axis_composition(axis_genes, categories_measured, categories=None):
    """How many axis genes each category source measures (the composition table).

    ``categories_measured`` is aligned with ``categories`` (default
    :data:`CATEGORIES`). A category that is not in the alignument — or whose
    aligned mask is ``None`` — was excluded from the fusion and contributes
    ``0`` (ADR 0068).
    """
    categories = tuple(CATEGORIES if categories is None else categories)
    if len(categories_measured) != len(categories):
        raise ValueError('a composition needs one mask per active category')
    axis_genes = np.asarray(axis_genes, dtype=np.int64)
    composition = {category: 0 for category in CATEGORIES}
    for category, mask in zip(categories, categories_measured):
        if mask is None:
            continue
        mask = np.asarray(mask, dtype=bool)
        if mask.ndim != 1:
            raise ValueError('a category measured mask must be a vector')
        if axis_genes.size and int(axis_genes.max()) >= mask.size:
            raise ValueError('an axis gene is outside the mask axis')
        if axis_genes.size and int(axis_genes.min()) < 0:
            raise ValueError('axis gene indices must be non-negative')
        composition[category] = int(mask[axis_genes].sum())
    return composition


def _source_coverage_inputs(axes, source_ranks, sizes):
    """Validate the source-coverage inputs shared by the count and percent tables.

    Returns ``(sizes, ranks, genes)``: the requested axis sizes as ints, each
    source's rank (``None`` for a category the fusion excludes) and the axis gene
    indices per size.
    """
    sizes = tuple(int(size) for size in sizes)
    ranks = {}
    for name, rank in source_ranks.items():
        if name != 'tenv_hvg' and name not in CATEGORIES:
            raise ValueError(f'unknown coverage source: {name}')
        if rank is None:
            ranks[name] = None
            continue
        rank = np.asarray(rank, dtype=np.int64)
        if rank.ndim != 1 or rank.size < 2:
            raise ValueError(f'{name}: a source coverage rank must be a vector')
        expected = np.arange(1, rank.size + 1, dtype=np.int64)
        if not np.array_equal(np.sort(rank), expected):
            raise ValueError(f'{name}: a source coverage rank must be a complete '
                             '1..axis_size permutation')
        ranks[name] = rank
    ranked = {name: rank for name, rank in ranks.items() if rank is not None}
    if not ranked:
        raise ValueError('a coverage table needs at least one ranked source')
    axis_size = int(next(iter(ranked.values())).size)
    for name, rank in ranked.items():
        if rank.size != axis_size:
            raise ValueError(f'{name}: source ranks must cover the same frozen axis')
    keyed = {int(key) for key in axes}
    genes = {}
    for size in sizes:
        if size <= 0:
            raise ValueError('a coverage axis size must be positive')
        if size not in keyed:
            raise ValueError(f'no released axis of size {size}')
        genes[size] = np.asarray(axes[size], dtype=np.int64)
        if genes[size].size != size:
            raise ValueError(f'G_{size} does not hold exactly {size} genes')
        if genes[size].min() < 0 or genes[size].max() >= axis_size:
            raise ValueError(f'G_{size} is outside the frozen axis')
    return sizes, ranks, genes


def source_topn_covered_counts(axes, source_ranks, sizes=AXIS_SIZES):
    """Each source's top-``N`` covered gene count on every released axis.

    The integer companion of :func:`source_topn_coverage`: ``source_ranks`` maps
    ``'tenv_hvg'`` and every category to a complete ``1..axis_size`` rank
    permutation indexed by ``gene_index``, and a category the fusion excludes
    (ADR 0068: every dataset non-informative) maps to ``None`` and keeps an empty
    column. For one axis ``G_N`` and source ``s`` the cell is::

        |{g in G_N : rank_s(g) <= N}|

    Because both sets hold ``N`` genes, the count is at once how many of ``s``'s
    top-``N`` genes the axis keeps and how many of the axis's genes come from that
    source's top-``N`` set; the sets overlap across sources, so a row does not sum
    to ``N``.
    """
    sizes, ranks, genes = _source_coverage_inputs(axes, source_ranks, sizes)
    rows = []
    for size in sizes:
        row = {'N': size}
        for name in ('tenv_hvg',) + CATEGORIES:
            rank = ranks.get(name)
            row[name] = (np.nan if rank is None
                         else int((rank[genes[size]] <= size).sum()))
        rows.append(row)
    return pd.DataFrame(rows, columns=list(COUNT_COLUMNS))


def source_topn_coverage(axes, source_ranks, sizes=AXIS_SIZES):
    """Each source's top-``N`` coverage of every released axis, in percent.

    ``source_ranks`` follows the contract of
    :func:`source_topn_covered_counts`. For one axis ``G_N`` and source ``s`` the
    cell is the percentage of ``s``'s top-``N`` set (``rank_s <= N``) that the
    axis contains::

        100 * |{g in G_N : rank_s(g) <= N}| / N

    Because both sets hold ``N`` genes, the count is at once the share of the
    source's top-``N`` the axis keeps and the share of the axis drawn from that
    source's top-``N``; the sets overlap across sources, so a row does not sum to
    100.
    """
    counts = source_topn_covered_counts(axes, source_ranks, sizes=sizes)
    rows = []
    for record in counts.to_dict('records'):
        size = int(record['N'])
        row = {'N': size}
        for name in ('tenv_hvg',) + CATEGORIES:
            value = record[name]
            row[f'%{name}_gN'] = (np.nan if pd.isna(value)
                                  else round(100.0 * int(value) / size, 2))
        rows.append(row)
    return pd.DataFrame(rows, columns=list(COVERAGE_COLUMNS))


def pn_distribution(p_n, n_sources=N_FUSION_SOURCES):
    """Counts of genes at each ``pN`` value ``0..n_sources`` (ordinal, not a score).

    ``n_sources`` is the number of sources actually fused (ADR 0068), not the
    declared maximum; the release passes its published source count.
    """
    p_n = np.asarray(p_n, dtype=np.int64)
    if p_n.ndim != 1:
        raise ValueError('pN must be a vector')
    if p_n.min() < 0 or p_n.max() > n_sources:
        raise ValueError(f'pN must be in 0..{n_sources}')
    return np.bincount(p_n, minlength=n_sources + 1)


def rank_sum_histogram(values, low=0.0, high=None, bins=70):
    """Histogram of ``rank_sum`` over all candidates, with fixed declared edges.

    The evidence figure draws this histogram once per released axis and overlays
    that axis's selected genes on every bar. Fixed equal-width edges (bounded by
    the published source count, ``0..n_sources`` in 70 bins) keep the plotted
    values deterministic and comparable between replays; every candidate falls
    inside the interval and every bin is emitted, including empty ones, so the
    CSV is the complete plotted value set.
    """
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 1 or values.size == 0:
        raise ValueError('rank_sum histogram needs a non-empty vector')
    if not np.isfinite(values).all():
        raise ValueError('rank_sum histogram values must be finite')
    high = float(np.ceil(values.max()) if high is None else high)
    low = float(low)
    bins = int(bins)
    if not low < high:
        raise ValueError('histogram interval must be non-empty')
    if bins <= 0:
        raise ValueError('histogram needs at least one bin')
    if values.min() < low or values.max() > high:
        raise ValueError('a rank_sum value is outside the declared histogram interval')
    edges = np.linspace(low, high, bins + 1)
    counts, edges = np.histogram(values, bins=edges)
    frame = pd.DataFrame({
        'bin_left': edges[:-1],
        'bin_right': edges[1:],
        'bin_center': (edges[:-1] + edges[1:]) / 2.0,
        'count': counts.astype(np.int64),
    })
    if int(frame['count'].sum()) != len(values):
        raise ValueError('rank_sum histogram does not cover every candidate')
    return frame


def kneedle_knee(x, y, sensitivity=1.0):
    """Online Kneedle: first interior local maximum of the normalized arc above the chord.

    Both axes are min-max normalized, ``arc = y_norm - x_norm``; the knee is the
    first interior local maximum reaching ``sensitivity * max(arc)``. A monotone
    or straight curve has no interior maximum and returns ``None`` -- the caller
    must report *no inflection* rather than invent an ``N``.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if len(x) != len(y) or len(x) < 3:
        return None
    span_x, span_y = x[-1] - x[0], y.max() - y.min()
    if span_x <= 0 or span_y <= 0:
        return None
    x_norm = (x - x[0]) / span_x
    y_norm = (y - y.min()) / span_y
    arc = y_norm - x_norm
    peak = float(arc.max())
    if peak <= 1e-12:
        return None
    interior = np.where((arc[1:-1] >= arc[:-2]) & (arc[1:-1] > arc[2:]))[0] + 1
    if not len(interior):
        return None
    accepts = interior[arc[interior] >= sensitivity * peak]
    if not len(accepts):
        return None
    position = int(accepts[0])
    return {'index': position, 'n': int(x[position]), 'arc': float(arc[position]),
            'arc_max': peak, 'knee_strength': float(arc[position])}


def second_difference_peak(y, min_n=100, max_n=None, minimum_curvature=1e-6):
    """Most concave point of the min-max normalized curve inside a declared range.

    The second difference is taken on the uniform rank grid and the search is
    restricted to ``[min_n, max_n]``; the leading positions below ``min_n`` are
    excluded because a rank-sum curve's first genes are strong outliers whose
    curvature would dominate any knee test. Returns ``None`` when the largest
    curvature is below ``minimum_curvature`` (a straight line has zero curvature,
    so it reports no candidate rather than inventing one).
    """
    y = np.asarray(y, dtype=np.float64)
    if len(y) < 3:
        return None
    span = y.max() - y.min()
    if span <= 0:
        return None
    y_norm = (y - y.min()) / span
    curvature = y_norm[2:] - 2.0 * y_norm[1:-1] + y_norm[:-2]
    positions = np.arange(2, len(y))
    upper = len(y) - 1 if max_n is None else min(int(max_n), len(y) - 1)
    lower = max(int(min_n), 2)
    mask = (positions >= lower) & (positions <= upper)
    if not mask.any():
        return None
    masked = np.where(mask, curvature, np.inf)
    packed = int(np.argmin(masked))
    magnitude = float(-masked[packed])
    if magnitude < minimum_curvature:
        return None
    return {'index': int(positions[packed] - 1), 'n': int(positions[packed]),
            'curvature': magnitude}


def detect_inflection(x, y, sensitivity=1.0, agreement_tolerance=0.25,
                      minimum_knee_strength=0.01, min_n=100):
    """Kneedle + second-difference inflection, reported as a diagnostic only.

    ADR 0066: the rank-sum curve gets the declared check and may honestly report
    ``"inflection": false``; it never chooses ``N``. An inflection is accepted
    only when both criteria find a candidate, the two candidate positions agree
    within ``agreement_tolerance`` (relative) and the Kneedle knee clears
    ``minimum_knee_strength`` above the chord.
    """
    kneedle = kneedle_knee(x, y, sensitivity=sensitivity)
    curvature = second_difference_peak(y, min_n=min_n)
    accepted = False
    candidate = None
    if kneedle is not None and curvature is not None:
        left, right = kneedle['n'], curvature['n']
        agreement = abs(left - right) / max(left, right)
        accepted = (agreement <= agreement_tolerance
                    and kneedle['knee_strength'] >= minimum_knee_strength)
        if accepted:
            candidate = int(round((left + right) / 2))
    return {
        'criterion': ('Kneedle online S=1.0 (first interior maximum of the normalized '
                      'arc above the chord) plus the second-difference peak of the same '
                      'normalized curve inside the declared decision range [min_n, n]; '
                      'both must agree within 25% and the knee must exceed 0.01 '
                      'normalized units'),
        'parameters': {'sensitivity': sensitivity,
                       'agreement_tolerance': agreement_tolerance,
                       'minimum_knee_strength': minimum_knee_strength,
                       'min_n': min_n,
                       'minimum_curvature': 1e-6},
        'kneedle': _drop_index(kneedle),
        'second_difference': _drop_index(curvature),
        'inflection': accepted,
        'n': candidate,
    }


def _drop_index(result):
    if result is None:
        return None
    return {key: value for key, value in result.items() if key != 'index'}


def inflection_stability(curves, tolerance=0.10):
    """Stability of the inflection answer under leave-one-source-out replays.

    ``curves`` is the primary ``(x, y)`` rank-sum curve followed by the curves of
    the leave-one-source-out replays. The diagnostic is stable only when every
    replay reproduces the primary's decision and, when an inflection is accepted,
    its ``N`` within ``tolerance`` (relative). Stability is reported, never used.
    """
    if not curves:
        raise ValueError('inflection stability needs at least the primary curve')
    outcomes = [detect_inflection(np.asarray(x, dtype=np.float64),
                                  np.asarray(y, dtype=np.float64))
                for x, y in curves]
    primary = outcomes[0]
    replays = outcomes[1:]
    same = [result['inflection'] == primary['inflection'] for result in replays]
    within = []
    for result in replays:
        if primary['inflection'] and result['inflection']:
            within.append(abs(result['n'] - primary['n'])
                          <= tolerance * max(primary['n'], 1))
        else:
            within.append(result['inflection'] == primary['inflection'])
    return {
        'method': f'leave-one-source-out; N within {tolerance:.0%} of the primary',
        'n_replays': len(replays),
        'n_same_decision': int(sum(same)),
        'n_within_tolerance': int(sum(within)),
        'stable': bool(replays) and all(same) and all(within),
        'primary': primary,
    }


def sha256(path):
    """Streaming sha256 of one file (reader-side digest check)."""
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


class GeneAxis:
    """Public reader for the published ``output/gene-axis/`` release.

    Every read is digest-checked against ``manifest.json``, so a consumer (the
    views of #191, the acceptance of #192) trusts the accepted bytes and never
    re-derives the fill rule. The reader is importable without the external item
    trees and opens only the gene-axis release directory.
    """

    def __init__(self, root):
        self.root = Path(root)

    # ---- structure -----------------------------------------------------
    def manifest(self):
        """The release manifest, after checking its schema and relative paths."""
        path = self.root / 'manifest.json'
        if not path.is_file():
            raise FileNotFoundError(f'no gene-axis manifest at {path}')
        manifest = json.loads(path.read_text())
        if manifest.get('schema') != SCHEMA:
            raise ValueError(f'{path}: unexpected gene-axis schema')
        return manifest

    def axis_sizes(self):
        """The declared released axis sizes, ascending."""
        return tuple(int(size) for size in self.manifest()['axes'])

    def digests(self):
        """Declared sha256 of every output (repo-relative name -> digest)."""
        return {name: record['sha256']
                for name, record in self.manifest()['outputs'].items()}

    # ---- payloads ------------------------------------------------------
    def _declared(self, basename):
        """The single manifest output whose basename is ``basename`` in gene-axis/."""
        manifest = self.manifest()
        matches = [name for name in manifest['outputs']
                   if Path(name).name == basename
                   and Path(name).parent.name == 'gene-axis']
        if len(matches) != 1:
            raise ValueError(f'{basename} is not exactly one declared gene-axis output')
        return matches[0], manifest['outputs'][matches[0]]

    def _verified(self, basename):
        name, declared = self._declared(basename)
        path = self.root / basename
        if not path.is_file():
            raise FileNotFoundError(f'missing declared gene-axis output {name}')
        if sha256(path) != declared['sha256']:
            raise ValueError(f'{name}: digest differs from the release manifest')
        return path

    def candidate_table(self):
        """All frozen-axis candidates with their source ranks and ``pN``."""
        path = self._verified('candidate_rank_sum.csv')
        table = pd.read_csv(path, dtype={'symbol': str}, float_precision='round_trip')
        expected = candidate_columns(self.axis_sizes())
        if list(table.columns) != list(expected):
            raise ValueError(f'{path}: candidate table columns differ from the release contract')
        return table

    def axis_frame(self, size):
        """The released axis as a frame (``gene_index``, ``symbol``), ascending."""
        path = self._verified(f'gene_axis_{int(size)}.csv')
        frame = pd.read_csv(path, dtype={'symbol': str})
        if list(frame.columns) != list(AXIS_COLUMNS):
            raise ValueError(f'{path}: axis file columns differ from the contract')
        return frame

    def axis(self, size):
        """The released axis gene indices (original integers), ascending."""
        return self.axis_frame(size).gene_index.to_numpy(dtype=np.int64)

    def composition(self):
        """The per-category / pN composition table of both released axes."""
        path = self._verified('axis_composition.csv')
        return pd.read_csv(path)

    def coverage(self):
        """Each source's top-``N`` coverage of both released axes (percent)."""
        path = self._verified('source_topn_coverage.csv')
        table = pd.read_csv(path)
        expected = list(COVERAGE_COLUMNS)
        if list(table.columns) != expected:
            raise ValueError(f'{path}: coverage table columns differ from the contract')
        return table
