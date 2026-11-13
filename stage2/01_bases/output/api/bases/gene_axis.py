"""Two-level rank-sum fusion and the released Stage 2 axes (ADR 0066).

This module owns every scientific definition of the released gene axis: the
winning-tier selection inside a disease category, the level-1 category rank, the
level-2 seven-source rank sum, the ``pN`` cross-source replication count, the
``0.5*N + 0.5*N`` axis fill rule, the rank-sum histogram and the diagnostic
Kneedle + second-difference inflection. ``execution/build_gene_axis.py`` reads
the pinned inputs and publishes the products; ``execution/tests/test_gene_axis.py``
exercises these definitions on hermetic fixtures.

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

from .deg_rank import AXIS_SIZE, normalized_rank

SCHEMA = 'stage2-gene-axis-release-v1'

# The six declared disease categories, in ADR 0066 roster order. The product
# directory is the catalog directory under ``output/disease-deg/``; it is also
# the prefix of the candidate table's ``<category>_rank`` column.
CATEGORIES = ('nhf', 'icm', 'af', 'genetic_cardiomyopathy',
              'cardiotox_TOP2', 'cardiotox_other')
CATEGORY_LABELS = {'nhf': 'nHF', 'icm': 'ICM', 'af': 'AF',
                   'genetic_cardiomyopathy': 'genetic_cardiomyopathy',
                   'cardiotox_TOP2': 'cardiotox_TOP2',
                   'cardiotox_other': 'cardiotox_other'}
CATEGORY_RANK_COLUMNS = tuple(f'{category}_rank' for category in CATEGORIES)

#: Declared released axis sizes (ADR 0066). N is chosen, never read off a curve.
AXIS_SIZES = (2000, 5000)
#: Half of an axis comes from tenv-HVG, half from cross-source replication.
FILL_FRACTION = 0.5
#: The seventh fusion source is the published tenv-HVG ranking; the six others
#: are the category class sources. This is the denominator of ``pN``.
N_FUSION_SOURCES = len(CATEGORIES) + 1

def candidate_columns(sizes=AXIS_SIZES):
    """Contracted candidate-table columns for the requested axis sizes."""
    return (('symbol', 'gene_index', 'tenv_hvg_rank') + CATEGORY_RANK_COLUMNS
            + ('rank_sum', 'n_sources_measured')
            + tuple(f'p{int(size)}' for size in sizes))


CANDIDATE_COLUMNS = candidate_columns()
AXIS_COLUMNS = ('gene_index', 'symbol')

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
    by_exception = []
    excluded = sorted({tier for tier in tiers if tier != winning})
    for record, tier in zip(records, tiers):
        isogenic = (record.get('isogenic_exception') or {}).get('applies') is True
        if tier == winning:
            retained.append(record)
        elif isogenic:
            retained.append(record)
            by_exception.append(record['dataset_id'])
    return {
        'winning_tier': winning,
        'datasets': [record['dataset_id'] for record in retained],
        'records': retained,
        'retained_by_exception': by_exception,
        'excluded_tiers': excluded,
        'rule': ('minimum catalog tier, plus records whose isogenic_exception.applies '
                 'is true (ADR 0066 Isogenic exception)'),
    }


def category_rank(ranks, n_measured, gene_index=None, axis_size=None):
    """Level 1: a category's source rank from its datasets' within-dataset ranks.

    ADR 0066 *Candidate table and the two-level rank-sum fusion*::

        level-1 sum = sum(rank_norm)  over the datasets at the winning tier
        category rank = position in ascending level-1 sum, ties by gene_index

    ``rank_norm`` is recomputed from the raw rank and the dataset's
    ``n_measured`` with :func:`bases.deg_rank.normalized_rank` (measured rank
    ``rank/n_measured``; unmeasured ``1.0``) so a category never trusts a
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
                    gene_index=None, axis_size=None, sizes=AXIS_SIZES):
    """Level 2: ``rank_sum`` over the seven fusion sources.

    ADR 0066 / spec D13: ``rank_sum = sum(rank_norm)`` over the six category
    class sources plus the published tenv-HVG ranking, each normalized by the
    frozen-axis size; **ascending is better** and ties break on ``gene_index``.
    ``pN`` counts how many of the seven sources place the gene in their top-``N``
    set; ``n_sources_measured`` counts the sources that measure it (tenv-HVG
    ranks every frozen-axis gene, so it always contributes one).

    ``category_ranks`` are the level-1 category ranks (one complete permutation
    per category, in :data:`CATEGORIES` order), ``categories_measured`` the
    corresponding "source measures the gene" masks, ``tenv_rank`` the published
    tenv-HVG dispersion rank. Returns ``{'rank_sum', 'n_sources_measured', 'pn'}``
    where ``pn`` maps each requested ``N`` to its integer vector.
    """
    if len(category_ranks) != len(CATEGORIES):
        raise ValueError(f'a level-2 fusion needs exactly {len(CATEGORIES)} category ranks')
    if len(categories_measured) != len(CATEGORIES):
        raise ValueError('every category rank needs its measured mask')
    tenv_rank = np.asarray(tenv_rank, dtype=np.int64)
    axis_size = int(axis_size if axis_size is not None else tenv_rank.size)
    gene_index = (np.arange(axis_size, dtype=np.int64) if gene_index is None
                  else np.asarray(gene_index, dtype=np.int64))
    if gene_index.shape != (axis_size,):
        raise ValueError('gene_index must cover the frozen axis')
    tenv_rank = _complete_rank(tenv_rank, axis_size, 'tenv-HVG')
    sources = [_complete_rank(rank, axis_size, category)
               for rank, category in zip(category_ranks, CATEGORIES)]
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
    return {'rank_sum': rank_sum, 'n_sources_measured': n_sources_measured, 'pn': pn}


def candidate_table(axis, tenv_rank, category_ranks, categories_measured,
                    gene_index=None, sizes=AXIS_SIZES):
    """The master candidate table: one row per frozen-axis gene.

    ADR 0066: ``symbol``, ``gene_index``, ``tenv_hvg_rank``, one
    ``<category>_rank`` per category, ``rank_sum``, ``n_sources_measured``,
    ``p2000``, ``p5000``. Rows are in frozen-axis order; every vector is indexed
    by ``gene_index`` so no consumer can read a row position as a gene index.
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
    fusion = fusion_rank_sum(category_ranks, categories_measured, tenv_rank,
                             gene_index=gene_index, axis_size=axis_size, sizes=sizes)
    columns = {'symbol': axis.symbol.to_numpy(), 'gene_index': gene_index,
               'tenv_hvg_rank': np.asarray(tenv_rank, dtype=np.int64)}
    for category, rank in zip(CATEGORIES, category_ranks):
        columns[f'{category}_rank'] = np.asarray(rank, dtype=np.int64)
    columns['rank_sum'] = fusion['rank_sum']
    columns['n_sources_measured'] = fusion['n_sources_measured']
    for size in sizes:
        columns[f'p{int(size)}'] = fusion['pn'][int(size)]
    table = pd.DataFrame(columns)
    return table.loc[:, list(candidate_columns(sizes))]


def axis_fill(tenv_rank, p_n, rank_sum, n, gene_index=None):
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
    if p_n.min() < 0 or p_n.max() > N_FUSION_SOURCES:
        raise ValueError(f'pN must be in 0..{N_FUSION_SOURCES}')
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


def axis_composition(axis_genes, categories_measured):
    """How many axis genes each category source measures (the composition table)."""
    if len(categories_measured) != len(CATEGORIES):
        raise ValueError(f'a composition needs {len(CATEGORIES)} category masks')
    axis_genes = np.asarray(axis_genes, dtype=np.int64)
    composition = {}
    for category, mask in zip(CATEGORIES, categories_measured):
        mask = np.asarray(mask, dtype=bool)
        if mask.ndim != 1:
            raise ValueError('a category measured mask must be a vector')
        if axis_genes.size and int(axis_genes.max()) >= mask.size:
            raise ValueError('an axis gene is outside the mask axis')
        if axis_genes.size and int(axis_genes.min()) < 0:
            raise ValueError('axis gene indices must be non-negative')
        composition[category] = int(mask[axis_genes].sum())
    return composition


def pn_distribution(p_n, n_sources=N_FUSION_SOURCES):
    """Counts of genes at each ``pN`` value ``0..n_sources`` (ordinal, not a score)."""
    p_n = np.asarray(p_n, dtype=np.int64)
    if p_n.ndim != 1:
        raise ValueError('pN must be a vector')
    if p_n.min() < 0 or p_n.max() > n_sources:
        raise ValueError(f'pN must be in 0..{n_sources}')
    return np.bincount(p_n, minlength=n_sources + 1)


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
    def _verified(self, name):
        manifest = self.manifest()
        declared = manifest['outputs'].get(name)
        if not isinstance(declared, dict):
            raise ValueError(f'{name} is not a declared gene-axis output')
        path = self.root.parent.parent / name
        if not path.is_file():
            path = self.root / Path(name).name
        if not path.is_file():
            raise FileNotFoundError(f'missing declared gene-axis output {name}')
        if sha256(path) != declared['sha256']:
            raise ValueError(f'{name}: digest differs from the release manifest')
        return path

    def candidate_table(self):
        """All frozen-axis candidates with their seven source ranks and ``pN``."""
        path = self._verified('output/gene-axis/candidate_rank_sum.csv')
        table = pd.read_csv(path, dtype={'symbol': str}, float_precision='round_trip')
        missing = set(CANDIDATE_COLUMNS) - set(table.columns)
        if missing:
            raise ValueError(f'{path}: candidate table lacks {sorted(missing)}')
        return table

    def axis_frame(self, size):
        """The released axis as a frame (``gene_index``, ``symbol``), ascending."""
        path = self._verified(f'output/gene-axis/gene_axis_{int(size)}.csv')
        frame = pd.read_csv(path, dtype={'symbol': str})
        missing = set(AXIS_COLUMNS) - set(frame.columns)
        if missing:
            raise ValueError(f'{path}: axis file lacks {sorted(missing)}')
        return frame

    def axis(self, size):
        """The released axis gene indices (original integers), ascending."""
        return self.axis_frame(size).gene_index.to_numpy(dtype=np.int64)

    def composition(self):
        """The per-category / pN composition table of both released axes."""
        path = self._verified('output/gene-axis/axis_composition.csv')
        return pd.read_csv(path)
