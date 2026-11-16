"""Cross-disease DEG fusion and the released Stage 2 axes (ADR 0066, ADR 0071).

This module owns every scientific definition of the released gene axis: the
cross-disease fusion of the published per-disease ``DEG1000`` sets, the
``f1000`` disease-replication count, the ``0.5*N + 0.5*N`` axis fill rule with
its tenv-HVG fallback, the ``DEG1000``-based composition and source-coverage
tables, the rank-sum histogram and the diagnostic Kneedle + second-difference
inflection. ``execution/build_gene_axis.py`` reads the pinned inputs and
publishes the products; ``execution/tests/test_gene_axis.py`` exercises these
definitions on hermetic fixtures.

Level 1 (a disease's source fusion and its ``DEG1000`` cut) is the published
product of ADR 0071 decision 4 and lives in :mod:`bases.disease_deg`
(``level_one_rank`` / ``fuse_disease``); this module consumes its
``output/disease-fusion/`` product through
:class:`bases.disease_deg.DiseaseFusion` and never re-reads a within-dataset
measured-gene rank. Level 2 (ADR 0071 decision 5) is::

    f1000(g)    = the number of active diseases whose DEG1000_d contains g
    rank_sum(g) = sum over the active diseases of rank_d(g) / n_d

where ``rank_d(g)`` is the gene's published disease rank when it is in
``DEG1000_d`` and ``n_d + 1`` otherwise, and ``n_d`` is that disease's
axis-resolved DEG-union size -- the same "divide by the number of ranked items"
contract ADR 0071 decision 3 applies to a source. ``rank_sum`` is consequently
in ``[1/n_d, (n_d + 1)/n_d]`` per disease and a gene outside ``DEG1000_d`` is
worse than every gene inside it, never omitted. Ascending is better; ties break
on ``gene_index``. tenv-HVG is **not** a ``rank_sum`` term: it supplies the
axis's first half, and its dispersion rank orders the fallback genes.

It is a pure function module (ADR 0062): no file is opened here except through
:class:`GeneAxis`, the public reader for the published ``output/gene-axis/``
release. Nothing here chooses an axis size: ``N`` is declared (2,000/3,000/4,000/5,000)
and the inflection diagnostic never selects it.

Vocabulary follows ADR 0071 and ``GLOSSARY.md`` -> *Stage 2 gene axis*:
``DEG (gene axis, v2)``, ``DEG1000``, ``f1000``, ``Rank-sum fusion``,
``Axis fill rule``, ``Inflection point``. The retired ``pN`` count has no
definition, column or reader here.

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

from . import disease_deg

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

SCHEMA = 'stage2-gene-axis-release-v3'
#: The axis-construction algorithm marker published in the release manifest.
#: A release whose manifest does not carry it predates the ADR 0071 v2 rule
#: (it was built under the retired ``pN`` fill rule) and its replay is refused
#: with a self-describing message (#207); #208 republishes the products.
ALGORITHM = 'tier-first-deg1000-f1000-four-axes-v3'

# ADR 0072 declares five diseases; measured eligibility selects participants.
# Directory names also prefix candidate rank columns and composition groups.
CATEGORIES = ('nhf', 'icm', 'af', 'genetic_cardiomyopathy', 'tox')
UNPUBLISHED_CATEGORIES = ()
RELEASE_CATEGORIES = CATEGORIES
CATEGORY_LABELS = {'nhf': 'nHF', 'icm': 'ICM', 'af': 'AF',
                   'genetic_cardiomyopathy': 'genetic_cardiomyopathy',
                   'cardiotox_TOP2': 'cardiotox_TOP2',
                   'cardiotox_other': 'cardiotox_other', 'tox': 'TOX'}
CATEGORY_RANK_COLUMNS = tuple(f'{category}_rank' for category in RELEASE_CATEGORIES)

#: Declared released axis sizes (ADR 0066). N is chosen, never read off a curve.
AXIS_SIZES = (2000, 3000, 4000, 5000)
#: Half of an axis comes from tenv-HVG, half from the cross-disease DEG1000
#: order (with the tenv-HVG fallback).
FILL_FRACTION = 0.5
#: The published ADR 0073 ``DEG1000`` cut: a disease keeps at most
#: its first 1,000 DEGs. Re-exported from the disease contract so the two levels
#: can never disagree.
DEG1000_LIMIT = disease_deg.DEG1000_LIMIT


def candidate_columns(categories=None):
    """Contracted candidate-table columns for the requested released roster.

    The tenv-HVG dispersion rank, the published ``DEG1000`` disease rank per
    released disease (empty where the gene is outside that disease's
    ``DEG1000``), the level-2 ``rank_sum`` and ``f1000``. The retired
    ``p2000``/``p5000`` columns have no successor here.
    """
    categories = tuple(RELEASE_CATEGORIES if categories is None else categories)
    if len(set(categories)) != len(categories):
        raise ValueError('a candidate table cannot repeat a category column')
    unknown = set(categories) - set(CATEGORIES)
    if unknown:
        raise ValueError(f'unknown candidate categories: {sorted(unknown)}')
    return (('symbol', 'gene_index', 'tenv_hvg_rank')
            + tuple(f'{category}_rank' for category in categories)
            + ('rank_sum', 'f1000'))


CANDIDATE_COLUMNS = candidate_columns()
AXIS_COLUMNS = ('gene_index', 'symbol')

#: The source-coverage table (long format): one row per released axis and per
#: fusion group, carrying the numerator and the denominator of the plotted rate
#: explicitly so no consumer has to guess which set the denominator is.
#: ``tenv_hvg`` keeps the v1 top-``N`` contract (numerator
#: ``|top_N(tenv) ∩ G_N|``, denominator ``N``); a disease reports
#: ``|DEG1000_d ∩ G_N|`` over ``|DEG1000_d|`` (ADR 0071 decision 8). A group
#: excluded from the fusion keeps an empty row, so the table reads as "not
#: fused", never as "fused but measured nothing".
COVERAGE_COLUMNS = ('N', 'n_axis', 'group', 'group_label', 'numerator',
                    'denominator', 'rate', 'note')

#: ``axis_composition.csv``: one row per released axis and per composition
#: record. ``record`` is ``axis_fill`` (the two declared halves and the tenv-HVG
#: fallback), ``deg1000`` (a disease's ``DEG1000`` coverage of the axis with its
#: own denominator) or ``f1000`` (how many axis genes each disease-replication
#: count has). ``n_axis`` is the axis size; ``share`` is
#: ``n_genes / denominator``.
COMPOSITION_COLUMNS = ('axis', 'n_axis', 'record', 'group', 'n_genes',
                       'denominator', 'share', 'note')

FILL_RULE = (
    'G_N = the top 0.5*N genes of tenv-HVG (tenv_hvg_rank <= N/2), then 0.5*N '
    'more from the remaining candidates that have f1000 >= 1, ordered by f1000 '
    'descending, rank_sum ascending, gene_index ascending; a gene with '
    'f1000 = 0 is never selected by the disease half. When fewer than 0.5*N '
    'remaining candidates have f1000 >= 1, the shortfall is filled from the '
    'unused tenv-HVG genes (the remaining candidates not already selected) in '
    'ascending tenv-HVG rank order -- never from a gene on disease evidence it '
    'does not have -- and the count is recorded in the release manifest. The '
    'union is published ordered by ascending gene_index and is never renumbered.')
F1000_RULE = (
    'f1000(g) is the number of active diseases whose published DEG1000_d '
    'contains g; it is the axis disease-side selection signal and replaces the '
    'retired pN. It is a count of diseases, not a rank or a score.')
RANK_SUM_RULE = (
    'rank_sum(g) = sum over the active diseases of rank_d(g) / n_d, where n_d '
    'is the disease axis-resolved DEG-union size and rank_d(g) is the '
    "disease's published tier-first DEG1000 rank when g is in DEG1000_d and n_d + 1 "
    'otherwise. Ascending is better; ties break on gene_index. tenv-HVG is not '
    'a term: it supplies the first axis half.')
FALLBACK_RULE = (
    'a gene that is not a DEG in any disease enters the axis only as a '
    'tenv-HVG fallback gene, selected in ascending tenv-HVG rank order after '
    'every available f1000 >= 1 candidate; the manifest records how many '
    'second-half slots came from that fallback.')


def _gene_index(gene_index, axis_size):
    """The frozen-axis index vector, validated."""
    gene_index = (np.arange(axis_size, dtype=np.int64) if gene_index is None
                  else np.asarray(gene_index, dtype=np.int64))
    if gene_index.shape != (axis_size,):
        raise ValueError('gene_index must cover the frozen axis')
    return gene_index


def _complete_rank(rank, axis_size, label):
    """Validate one source's frozen-axis rank vector is a complete 1..axis_size permutation."""
    rank = np.asarray(rank, dtype=np.int64)
    if rank.shape != (axis_size,):
        raise ValueError(f'{label} rank must cover the frozen axis')
    if not np.array_equal(np.sort(rank), np.arange(1, axis_size + 1, dtype=np.int64)):
        raise ValueError(f'{label} rank must be a complete 1..{axis_size} permutation')
    return rank


def category_rank(ranks, n_deg, gene_index=None, axis_size=None):
    """Level 1: a category's source rank from its datasets' **DEG-only** ranks.

    ADR 0071 decisions 3 and 4 (this replaced the retired measured-gene contract
    of ADR 0066; #205)::

        level-1 sum = sum(rank_norm)  over the category's informative datasets
        category rank = position in ascending level-1 sum, ties by gene_index

    Each ``ranks`` entry is a dataset's full frozen-axis **DEG rank frame**: its
    DEGs carry ``1..n_DEG`` and every other gene carries ``n_DEG + 1``
    (:func:`bases.deg_rank.complete_deg_rank`), and the matching ``n_deg`` entry is
    that dataset's DEG count. ``rank_norm`` is ``rank / n_DEG``
    (:func:`bases.deg_rank.normalized_deg_rank`), so a non-DEG contributes
    ``(n_DEG + 1) / n_DEG`` -- strictly worse than every DEG of that source and
    never clamped to 1.0. The denominator is the **DEG** count; the measured-gene
    rank and its ``n_measured`` denominator are retired (ADR 0071 decision 3).

    A dataset with ``n_DEG = 0`` has no defined ``rank_norm``; it is
    non-informative and this function rejects it rather than dividing by zero.
    The caller excludes it and records the exclusion (ADR 0068), exactly as the
    disease contract does (:func:`bases.disease_deg.fuse_disease`). The shared
    implementation is :func:`bases.disease_deg.level_one_rank`; the released axis
    consumes the published level-1 product (``output/disease-fusion/``) and never
    re-derives this from a within-dataset frame.

    Returns ``{'rank', 'sum', 'deg', 'measured', 'n_datasets'}``: a 1-based
    ``1..axis_size`` permutation, the level-1 ``sum(rank_norm)``, the DEG-union
    mask (a DEG in at least one dataset) and the dataset count. ``measured`` is
    the retired v1 name for that same union mask, kept because a pre-v2 reader
    still reads it; ``deg`` is the v2 name.
    """
    result = disease_deg.level_one_rank(ranks, n_deg, gene_index=gene_index,
                                        axis_size=axis_size)
    return {'rank': result['rank'], 'sum': result['sum'], 'deg': result['deg'],
            'measured': result['deg'], 'n_datasets': result['n_sources']}


def deg1000_frame(rank, n_deg, axis_size=None, label='disease',
                  limit=DEG1000_LIMIT):
    """Validate one disease's complete frozen-axis ``DEG1000`` rank frame.

    The frame is the level-1 disease rank restricted to the published
    ``DEG1000_d``: a gene inside ``DEG1000_d`` carries its published
    ``disease_rank`` (``1..min(n_deg, limit)``) and every other gene carries
    ``n_deg + 1`` (the ADR 0071 decision-5 "outside ``DEG1000_d``" marker). A
    frame that is not a complete ``1..k`` assignment over
    ``k = min(n_deg, limit)`` genes, or that leaves a non-member at anything
    other than ``n_deg + 1``, is rejected rather than silently re-ordered.
    """
    rank = np.asarray(rank, dtype=np.int64)
    n_deg = int(n_deg)
    limit = int(limit)
    if n_deg <= 0:
        raise ValueError(f'{label}: a DEG1000 frame needs a positive DEG-union size')
    if axis_size is None:
        axis_size = int(rank.size)
    axis_size = int(axis_size)
    if rank.shape != (axis_size,):
        raise ValueError(f'{label}: a DEG1000 rank frame must cover the frozen axis')
    expected = min(n_deg, limit)
    if rank.min() < 1 or rank.max() > n_deg + 1:
        raise ValueError(f'{label}: a DEG1000 rank must be in 1..n_deg+1')
    inside = rank <= expected
    if int(inside.sum()) != expected:
        raise ValueError(
            f'{label}: DEG1000 must hold min(n_deg, {limit}) = {expected} genes, '
            f'found {int(inside.sum())}')
    if not np.array_equal(np.sort(rank[inside]),
                          np.arange(1, expected + 1, dtype=np.int64)):
        raise ValueError(f'{label}: DEG1000 ranks must be exactly 1..{expected}')
    if not inside.all() and not (rank[~inside] == n_deg + 1).all():
        raise ValueError(f'{label}: a gene outside DEG1000 must carry rank n_deg+1')
    return rank


def _check_fusion_inputs(deg1000_ranks, n_deg, categories, gene_index, axis_size,
                         limit):
    """Shared level-2 input validation; returns (frames, counts, categories, size)."""
    categories = tuple(RELEASE_CATEGORIES if categories is None else categories)
    ranks = [np.asarray(rank, dtype=np.int64) for rank in deg1000_ranks]
    counts = [int(value) for value in n_deg]
    if not ranks:
        raise ValueError('a cross-disease fusion needs at least one disease')
    if len(ranks) != len(counts) or len(ranks) != len(categories):
        raise ValueError('a level-2 fusion needs one DEG1000 frame per active disease')
    if len(set(categories)) != len(categories):
        raise ValueError('a fusion source list cannot repeat a category')
    unknown = set(categories) - set(CATEGORIES)
    if unknown:
        raise ValueError(f'unknown fusion categories: {sorted(unknown)}')
    axis_size = int(axis_size if axis_size is not None else ranks[0].size)
    _gene_index(gene_index, axis_size)
    return ranks, counts, categories, axis_size


def f1000_count(deg1000_ranks, n_deg, categories=None, gene_index=None,
                axis_size=None, limit=DEG1000_LIMIT):
    """``f1000``: the number of active diseases whose ``DEG1000_d`` holds a gene.

    ``deg1000_ranks`` is one complete frozen-axis ``DEG1000`` rank frame per
    active disease (:func:`deg1000_frame`) and ``n_deg`` its matching
    DEG-union size. Membership is ``rank <= min(n_deg, limit)``, so a disease
    with a union under the limit keeps its whole union. The count is bounded by
    the number of active diseases, never by a declared constant: a category
    dropped from the fusion (ADR 0068) is passed by omission and cannot be
    counted.
    """
    ranks, counts, categories, axis_size = _check_fusion_inputs(
        deg1000_ranks, n_deg, categories, gene_index, axis_size, limit)
    total = np.zeros(axis_size, dtype=np.int64)
    for position, (frame, count) in enumerate(zip(ranks, counts)):
        frame = deg1000_frame(frame, count, axis_size, categories[position],
                              limit=limit)
        total += frame <= min(int(count), int(limit))
    return total


def cross_disease_rank_sum(deg1000_ranks, n_deg, categories=None, gene_index=None,
                           axis_size=None, limit=DEG1000_LIMIT):
    """Level 2: ``rank_sum`` and ``f1000`` over the active diseases (ADR 0071).

    ADR 0071 decision 5, with the decision-3 denominator contract::

        rank_sum(g) = sum over diseases of rank_d(g) / n_d

    ``rank_d(g)`` is the gene's published ``DEG1000`` disease rank when it is
    inside ``DEG1000_d`` and ``n_d + 1`` when it is not, so a gene outside a
    disease's ``DEG1000`` contributes ``(n_d + 1) / n_d`` -- worse than every
    gene inside that disease and never 0 or an omitted term. ``n_d`` is the
    disease's axis-resolved DEG-**union** size, not the number of genes that
    survived the 1,000 cut, so the same rank position means the same thing in
    every disease whose union has the same size. tenv-HVG is deliberately not a
    term: it is not a disease, and its contribution is the axis's first half.

    ``deg1000_ranks``/``n_deg`` describe the active diseases (a category dropped
    from the fusion is passed by omission, ADR 0068); validation is
    :func:`deg1000_frame`'s. Ascending ``rank_sum`` is better and ties break on
    ``gene_index`` in the fill rule, not here.

    Returns ``{'rank_sum', 'f1000', 'deg1000', 'categories', 'n_deg'}`` where
    ``deg1000`` is the per-disease membership mask list.
    """
    ranks, counts, categories, axis_size = _check_fusion_inputs(
        deg1000_ranks, n_deg, categories, gene_index, axis_size, limit)
    total = np.zeros(axis_size, dtype=np.float64)
    f1000 = np.zeros(axis_size, dtype=np.int64)
    masks = []
    for position, (frame, count) in enumerate(zip(ranks, counts)):
        frame = deg1000_frame(frame, count, axis_size, categories[position],
                              limit=limit)
        total += frame / float(count)
        mask = frame <= min(count, int(limit))
        f1000 += mask
        masks.append(mask)
    return {'rank_sum': total, 'f1000': f1000, 'deg1000': masks,
            'categories': categories, 'n_deg': counts}


def candidate_table(axis, tenv_rank, deg1000_ranks, n_deg, categories=None,
                    gene_index=None, limit=DEG1000_LIMIT):
    """The master candidate table: one row per frozen-axis gene (ADR 0071).

    Columns (:func:`candidate_columns`): ``symbol``, ``gene_index``,
    ``tenv_hvg_rank``, one ``<category>_rank`` per **released** category, the
    level-2 ``rank_sum`` and ``f1000``. A ``<category>_rank`` cell is the
    published ``DEG1000`` disease rank (``1..|DEG1000_d|``) when the gene is in
    that disease's ``DEG1000`` and empty otherwise; a category excluded from the
    fusion (ADR 0068: no ``DEG1000``) keeps its contracted column empty. Rows are
    in frozen-axis order and every vector is indexed by ``gene_index``, so no
    consumer can read a row position as a gene index.
    """
    axis = pd.DataFrame(axis)
    required = {'gene_index', 'symbol'}
    if not required <= set(axis.columns):
        raise ValueError(f'axis table requires {sorted(required)}')
    axis_size = len(axis)
    gene_index = _gene_index(gene_index, axis_size)
    if not np.array_equal(np.asarray(axis.gene_index, dtype=np.int64), gene_index):
        raise ValueError('axis table must be in frozen gene_index order')
    tenv_rank = _complete_rank(tenv_rank, axis_size, 'tenv-HVG')
    active = tuple(RELEASE_CATEGORIES if categories is None else categories)
    frames = [np.asarray(rank, dtype=np.int64) for rank in deg1000_ranks]
    counts = [int(value) for value in n_deg]
    fusion = cross_disease_rank_sum(frames, counts, categories=active,
                                    gene_index=gene_index, axis_size=axis_size,
                                    limit=limit)
    ranks_by_category = {}
    for category, frame, count in zip(active, frames, counts):
        frame = deg1000_frame(frame, count, axis_size, category, limit=limit)
        inside = frame <= min(count, int(limit))
        ranks_by_category[category] = np.where(inside, frame, np.nan)
    columns = {'symbol': axis.symbol.to_numpy(),
               'gene_index': gene_index,
               'tenv_hvg_rank': tenv_rank}
    for category in RELEASE_CATEGORIES:
        if category in ranks_by_category:
            columns[f'{category}_rank'] = ranks_by_category[category].astype(np.float64)
        else:
            columns[f'{category}_rank'] = np.full(axis_size, np.nan)
    columns['rank_sum'] = fusion['rank_sum']
    columns['f1000'] = fusion['f1000']
    table = pd.DataFrame(columns)
    return table.loc[:, list(candidate_columns())]


def axis_fill(tenv_rank, f1000, rank_sum, n, gene_index=None, n_diseases=None):
    """One released axis: the declared ``0.5*N + 0.5*N`` fill rule (ADR 0071).

    ``G_N`` takes the top ``0.5*N`` genes of tenv-HVG (``tenv_hvg_rank <= N/2``),
    then ``0.5*N`` more from the **remaining** candidates that have
    ``f1000 >= 1``, ordered by ``f1000`` descending, ``rank_sum`` ascending,
    ``gene_index`` ascending. When fewer than ``0.5*N`` remaining candidates have
    ``f1000 >= 1``, the shortfall is filled from the unused tenv-HVG genes in
    ascending tenv-HVG rank order -- never from a gene with ``f1000 > 0`` that
    lost the ordering, and never from a non-DEG gene outside that fallback. The
    result is returned as ascending original ``gene_index`` integers: an axis is
    an ordered subset of the frozen space, never renumbered to ``0..N-1``.

    ``f1000``/``rank_sum`` are the level-2 vectors of
    :func:`cross_disease_rank_sum` and ``n_diseases`` the number of active
    diseases (``f1000``'s bound); a ``pN`` vector has no successor here.

    Returns ``{'axis', 'tenv_half', 'from_deg1000', 'fallback',
    'n_deg1000_candidates', 'n_diseases'}``: the axis, the declared tenv half,
    the two second-half contributions and the number of remaining candidates
    that had ``f1000 >= 1``. ``from_deg1000 + fallback == n - tenv_half``.
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
    f1000 = np.asarray(f1000, dtype=np.int64)
    rank_sum = np.asarray(rank_sum, dtype=np.float64)
    if f1000.shape != (axis_size,) or rank_sum.shape != (axis_size,):
        raise ValueError('f1000 and rank_sum must cover the frozen axis')
    if f1000.min() < 0:
        raise ValueError('f1000 cannot be negative')
    n_diseases = int(f1000.max()) if n_diseases is None else int(n_diseases)
    if n_diseases < 0:
        raise ValueError('the active disease count cannot be negative')
    if f1000.max() > n_diseases:
        raise ValueError(f'f1000 must be in 0..{n_diseases} (the active disease count)')
    if not np.isfinite(rank_sum).all():
        raise ValueError('rank_sum must be finite')
    gene_index = _gene_index(gene_index, axis_size)
    tenv_half = np.flatnonzero(tenv_rank <= half)
    if len(tenv_half) != half:
        raise ValueError('tenv-HVG rank is not a complete permutation')
    remaining = np.flatnonzero(tenv_rank > half)
    if len(remaining) < half:
        raise ValueError('the frozen axis cannot fill the declared second half')
    disease_candidates = remaining[f1000[remaining] >= 1]
    # lexsort's last key is primary: f1000 descending, rank_sum ascending, index.
    order = np.lexsort((gene_index[disease_candidates],
                        rank_sum[disease_candidates],
                        -f1000[disease_candidates]))
    chosen = disease_candidates[order[:half]]
    fallback = np.empty(0, dtype=np.int64)
    if len(chosen) < half:
        shortfall = half - len(chosen)
        rest = np.setdiff1d(remaining, chosen)
        rest = rest[np.argsort(tenv_rank[rest], kind='stable')]
        fallback = rest[:shortfall]
        if len(fallback) != shortfall:
            raise ValueError(
                'the tenv-HVG fallback cannot fill the second half: the frozen '
                'axis has no unused candidate left')
        chosen = np.concatenate([chosen, fallback])
    axis = np.sort(np.concatenate([tenv_half, chosen]).astype(np.int64))
    if len(np.unique(axis)) != n or len(axis) != n:
        raise ValueError('axis fill did not select N distinct frozen-axis genes')
    return {'axis': axis, 'tenv_half': int(len(tenv_half)),
            'from_deg1000': int(len(chosen) - len(fallback)),
            'fallback': int(len(fallback)),
            'n_deg1000_candidates': int(len(disease_candidates)),
            'n_diseases': int(n_diseases)}


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


def deg1000_composition(axis_genes, deg1000_masks, deg1000_sizes, categories=None):
    """Per disease ``|DEG1000_d ∩ G_N|`` over ``|DEG1000_d|`` (ADR 0071 decision 8).

    ``deg1000_masks``/``deg1000_sizes`` are aligned with ``categories``; a
    category excluded from the fusion is passed with a ``None`` mask and keeps an
    empty entry (``None``), so a consumer cannot read it as "fused but measured
    nothing". Returns ``{category: {'n_covered', 'n_deg1000', 'share'}}`` where
    ``share = n_covered / n_deg1000``.
    """
    categories = tuple(RELEASE_CATEGORIES if categories is None else categories)
    if len(deg1000_masks) != len(categories) or len(deg1000_sizes) != len(categories):
        raise ValueError('a composition needs one DEG1000 mask and size per active category')
    axis_genes = np.asarray(axis_genes, dtype=np.int64)
    if axis_genes.size and int(axis_genes.min()) < 0:
        raise ValueError('axis gene indices must be non-negative')
    composition = {category: None for category in RELEASE_CATEGORIES}
    for category, mask, size in zip(categories, deg1000_masks, deg1000_sizes):
        if mask is None:
            continue
        mask = np.asarray(mask, dtype=bool)
        size = int(size)
        if mask.ndim != 1:
            raise ValueError('a DEG1000 mask must be a vector')
        if size <= 0 or int(mask.sum()) != size:
            raise ValueError(f'{category}: the DEG1000 mask size and its denominator differ')
        if axis_genes.size and int(axis_genes.max()) >= mask.size:
            raise ValueError('an axis gene is outside the DEG1000 mask axis')
        covered = int(mask[axis_genes].sum()) if axis_genes.size else 0
        composition[category] = {'n_covered': covered, 'n_deg1000': size,
                                 'share': covered / size}
    return composition


def f1000_distribution(f1000, n_diseases=None):
    """Counts of genes at each ``f1000`` value ``0..n_diseases``.

    ``n_diseases`` is the number of active diseases actually fused (a category
    dropped for having no ``DEG1000`` is not counted); the release passes its
    published count so the distribution's bound follows the data, never a
    hard-coded constant.
    """
    f1000 = np.asarray(f1000, dtype=np.int64)
    if f1000.ndim != 1:
        raise ValueError('f1000 must be a vector')
    n_diseases = int(f1000.max()) if n_diseases is None else int(n_diseases)
    if n_diseases < 0:
        raise ValueError('the active disease count cannot be negative')
    if f1000.min() < 0 or f1000.max() > n_diseases:
        raise ValueError(f'f1000 must be in 0..{n_diseases} (the active disease count)')
    return np.bincount(f1000, minlength=n_diseases + 1)


def source_topn_coverage(axes, tenv_rank, deg1000_masks, categories=None,
                         sizes=AXIS_SIZES):
    """Each group's coverage of every released axis, with numerator and denominator.

    Long format (:data:`COVERAGE_COLUMNS`): one row per released axis and per
    fusion group. ``tenv_hvg`` keeps the v1 top-``N`` contract (numerator
    ``|{g in G_N : tenv_hvg_rank(g) <= N}|``, denominator ``N``). A disease
    reports the ADR 0071 decision-8 ``DEG1000`` coverage: numerator
    ``|DEG1000_d ∩ G_N|``, denominator ``|DEG1000_d|``, so the rate is comparable
    between ``G_2000`` and ``G_5000`` and between diseases of different union
    size. A group excluded from the fusion (``None`` mask) keeps a row with empty
    numerator/denominator/rate and a note.
    """
    categories = tuple(RELEASE_CATEGORIES if categories is None else categories)
    masks = list(deg1000_masks)
    if len(masks) != len(categories):
        raise ValueError('a coverage table needs one DEG1000 mask per active category')
    sizes = tuple(int(size) for size in sizes)
    tenv_rank = np.asarray(tenv_rank, dtype=np.int64)
    if tenv_rank.ndim != 1 or tenv_rank.size < 2:
        raise ValueError('a tenv-HVG coverage rank must be a vector')
    if not np.array_equal(np.sort(tenv_rank),
                          np.arange(1, tenv_rank.size + 1, dtype=np.int64)):
        raise ValueError('a tenv-HVG coverage rank must be a complete permutation')
    axis_size = int(tenv_rank.size)
    by_category = {}
    for category, mask in zip(categories, masks):
        if mask is None:
            by_category[category] = None
            continue
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != (axis_size,):
            raise ValueError(f'{category}: a coverage DEG1000 mask must cover the frozen axis')
        if not mask.any():
            raise ValueError(f'{category}: a fused disease needs a non-empty DEG1000')
        by_category[category] = mask
    keyed = {int(key) for key in axes}
    genes_by_size = {}
    for size in sizes:
        if size <= 0:
            raise ValueError('a coverage axis size must be positive')
        if size not in keyed:
            raise ValueError(f'no released axis of size {size}')
        genes = np.asarray(axes[size], dtype=np.int64)
        if genes.size != size:
            raise ValueError(f'G_{size} does not hold exactly {size} genes')
        if genes.min() < 0 or genes.max() >= axis_size:
            raise ValueError(f'G_{size} is outside the frozen axis')
        genes_by_size[size] = genes
    rows = []
    for size in sizes:
        genes = genes_by_size[size]
        hit = int((tenv_rank[genes] <= size).sum())
        rows.append({'N': size, 'n_axis': size, 'group': 'tenv_hvg',
                     'group_label': 'tenv-HVG', 'numerator': hit,
                     'denominator': size,
                     'rate': round(100.0 * hit / size, 2),
                     'note': ('tenv-HVG top-N: |{g in G_N : tenv_hvg_rank(g) <= N}| / N')})
        for category in RELEASE_CATEGORIES:
            mask = by_category.get(category)
            if category not in by_category or mask is None:
                rows.append({'N': size, 'n_axis': size, 'group': category,
                             'group_label': CATEGORY_LABELS[category],
                             'numerator': np.nan, 'denominator': np.nan,
                             'rate': np.nan,
                             'note': ('excluded from cross-disease fusion by eligibility; no participating '
                                      'DEG1000, so no coverage is defined')})
                continue
            covered = int(mask[genes].sum())
            denominator = int(mask.sum())
            rows.append({'N': size, 'n_axis': size, 'group': category,
                         'group_label': CATEGORY_LABELS[category],
                         'numerator': covered, 'denominator': denominator,
                         'rate': round(100.0 * covered / denominator, 2),
                         'note': ('DEG1000 coverage: |DEG1000_d ∩ G_N| / |DEG1000_d|')})
    frame = pd.DataFrame(rows, columns=list(COVERAGE_COLUMNS))
    for column in ('numerator', 'denominator'):
        frame[column] = frame[column].astype('Int64')
    return frame


def rank_sum_histogram(values, low=0.0, high=None, bins=70):
    """Histogram of ``rank_sum`` over all candidates, with fixed equal-width edges.

    The evidence figure draws this histogram once per released axis and overlays
    that axis's selected genes on every bar. Fixed equal-width edges keep the
    plotted values deterministic and comparable between replays; every candidate
    falls inside the interval and every bin is emitted, including empty ones, so
    the CSV is the complete plotted value set. ``high`` defaults to the ceiling
    of the largest value (``rank_sum`` is bounded by the active disease count
    plus each disease's ``1/n_d`` penalty, not by a declared source count).
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
    """Stability of the inflection answer under leave-one-disease-out replays.

    ``curves`` is the primary ``(x, y)`` rank-sum curve followed by the curves of
    the leave-one-disease-out replays. The diagnostic is stable only when every
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
    trees and opens only the gene-axis release directory. The committed release
    predates the ADR 0071 v2 rule until #208 republishes it; a v1 product's
    ``candidate_table``/``coverage`` no longer match the v2 contract and fail
    closed rather than being silently re-read under the new semantics.
    """

    def __init__(self, root):
        self.root = Path(root)

    # ---- structure -----------------------------------------------------
    def manifest(self):
        """The release manifest, after checking its schema."""
        path = self.root / 'manifest.json'
        if not path.is_file():
            raise FileNotFoundError(f'no gene-axis manifest at {path}')
        manifest = json.loads(path.read_text())
        if manifest.get('schema') != SCHEMA:
            raise ValueError(f'{path}: unexpected gene-axis schema')
        if manifest.get('release', {}).get('algorithm') != ALGORITHM:
            raise ValueError(f'{path}: unexpected gene-axis algorithm')
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
        """All frozen-axis candidates with their disease ranks, ``rank_sum`` and ``f1000``."""
        path = self._verified('candidate_rank_sum.csv')
        table = pd.read_csv(path, dtype={'symbol': str}, float_precision='round_trip')
        expected = candidate_columns()
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
        """The per-axis fill / DEG1000 / f1000 composition table of both axes."""
        path = self._verified('axis_composition.csv')
        table = pd.read_csv(path)
        if list(table.columns) != list(COMPOSITION_COLUMNS):
            raise ValueError(f'{path}: composition columns differ from the contract')
        return table

    def coverage(self):
        """Each group's numerator/denominator coverage of both released axes."""
        path = self._verified('source_topn_coverage.csv')
        table = pd.read_csv(path)
        expected = list(COVERAGE_COLUMNS)
        if list(table.columns) != expected:
            raise ValueError(f'{path}: coverage table columns differ from the contract')
        return table
