"""Disease-DEG table/rank contract for the Stage 2 gene axis (ADR 0066, ADR 0071).

Each source dataset publishes one compact DEG table and one within-dataset rank
against this module. The contract is deliberately pure and dataset-agnostic:

* **alignment** — a source symbol resolves to one frozen ``gene_index`` through
  canonical symbol, the pinned identity crosswalk's alias/previous evidence, the
  pinned HGNC alias/previous columns, or Ensembl id. Ambiguous, off-axis and
  unmapped symbols are reported, never guessed.
* **effect score** — ``|log2FC| * -log10(padj)`` per contrast. A missing/``NA``
  ``padj`` (e.g. DESeq2 independent filtering) is non-significant and scores 0;
  tiny ``padj`` values are floored so the score stays finite.
* **multi-contrast collapse** — a gene's score is the **maximum** over its
  contrasts; direction and p-value come from that winning contrast.
* **DEG flag (v2)** — a gene is a DEG of its dataset when one of its contrasts
  passes ``padj < 0.05`` **and** ``|log2FC| >= 2`` (:func:`deg_flag`). The flag
  is pass-any-contrast, ``NA``/non-finite ``padj`` is non-significant, the
  winning contrast is still the maximum-effect one, and direction is carried,
  not filtered on (ADR 0071 *DEG definition*).
* **within-dataset DEG rank (v2)** — a dataset's **DEGs** rank ``1..n_DEG`` by
  descending effect score with ``gene_index`` as the deterministic tie-break; a
  gene that is not a DEG in that dataset gets ``n_DEG + 1``.
* **rank normalization (v2)** — ``rank_norm = rank / n_DEG`` for every gene, so
  a non-DEG's value is ``(n_DEG + 1) / n_DEG``. The denominator is the DEG count,
  never the measured-gene count, which is retired (ADR 0071 decision 3).
* **informativeness** — a dataset with no DEGs, or whose measured effect scores
  are all tied, has no informative ordering (its rank only repeats a constant or
  the ``gene_index`` tie-break). It is published with ``informative: false`` and
  a recorded reason, and excluded from fusion (ADR 0068, kept by ADR 0071).

``rank_scores``/``complete_rank``/``normalized_rank`` stay exported as generic
rank primitives for the pre-v2 fusion reader (``bases.gene_axis``); the DEG
contract itself is the ``*_deg_*`` family below. ``RANK_SEMANTICS`` names the
semantics every published table's ``stats`` carries; ``SCHEMA`` stays the
product-schema label until ticket 8 re-publishes the products. Ranks are
within-dataset only; they are never compared across datasets as if they were
effect sizes.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# The published product schema label stays ``v1`` until ticket 8 re-publishes the
# DEG products on the v2 contract; bumping it here would relabel products that are
# not rebuilt yet and would weaken the gene-axis ``derived_inputs`` fail-closed
# gate. The rank semantics are declared separately so a consumer can tell which
# contract produced a table before that republish.
SCHEMA = 'stage2-disease-deg-v1'
RANK_SEMANTICS = 'deg-only-v2'
AXIS_SIZE = 18164
PADJ_THRESHOLD = 0.05
# ADR 0071 *DEG definition*: padj < 0.05 AND |log2FC| >= 2 (FC >= 4). The padj
# test is strict and the effect-size test is inclusive, so padj == 0.05 is not a
# DEG while |log2FC| == 2 is one.
DEG_LOG2FC_THRESHOLD = 2.0
# -log10(0) is infinite; padj values are floored so the effect score stays finite
# and rank order is deterministic.
PADJ_FLOOR = 1e-300

ALIGNMENT_CATEGORIES = ('canonical', 'crosswalk_alias', 'hgnc_alias', 'ensembl',
                        'ambiguous', 'offaxis', 'unmapped')


class SymbolResolver:
    """Align source gene symbols onto the frozen axis without guessing.

    Precedence, each step only used when the previous one did not resolve:
    1. the axis canonical symbol (``output/gene-identity/genes.csv``);
    2. the identity crosswalk's ``alias``/``previous`` rows
       (``tables/gene_crosswalk.csv``), the repository's per-dataset evidence;
    3. the pinned HGNC complete set's ``alias_symbol``/``prev_symbol`` columns;
    4. Ensembl gene ID, either the axis ``ensembl_id`` column or the HGNC
       ``ensembl_gene_id`` column (version suffix stripped).

    A symbol with more than one on-axis candidate is *ambiguous* and stays
    unresolved; it is never assigned to an arbitrary gene. Off-axis and
    unmapped symbols are likewise reported, never silently dropped.
    """

    def __init__(self, axis, crosswalk=(), hgnc=()):
        axis = pd.DataFrame(axis)
        if not {'gene_index', 'symbol', 'hgnc_id', 'ensembl_id'} <= set(axis.columns):
            raise ValueError('axis table requires gene_index, symbol, hgnc_id, ensembl_id')
        self.axis_size = len(axis)
        self._canonical = {}
        self._by_hgnc = {}
        self._by_ensembl = {}
        for row in axis.itertuples(index=False):
            index = int(row.gene_index)
            self._canonical[str(row.symbol)] = index
            self._by_hgnc[str(row.hgnc_id)] = index
            ensembl = '' if row.ensembl_id is None else str(row.ensembl_id)
            if ensembl:
                self._by_ensembl[ensembl.split('.')[0]] = index
        self._crosswalk_alias = _crosswalk_aliases(crosswalk)
        self._crosswalk_ensembl = _crosswalk_ensembl(crosswalk)
        self._hgnc_alias = _hgnc_aliases(hgnc)
        self._hgnc_ensembl = _hgnc_ensembl(hgnc)

    def resolve(self, symbol):
        """Return ``(gene_index, category)``; ``gene_index`` is -1 when unresolved."""
        symbol = str(symbol)
        if symbol in self._canonical:
            return self._canonical[symbol], 'canonical'
        for mapping, category in ((self._crosswalk_alias, 'crosswalk_alias'),
                                  (self._hgnc_alias, 'hgnc_alias')):
            candidates = mapping.get(symbol)
            if candidates:
                on_axis = {self._by_hgnc[hgnc] for hgnc in candidates if hgnc in self._by_hgnc}
                if len(on_axis) == 1:
                    return next(iter(on_axis)), category
                if len(on_axis) > 1:
                    return -1, 'ambiguous'
        ensembl = symbol.split('.')[0] if symbol.startswith('ENS') else ''
        if ensembl:
            if ensembl in self._by_ensembl:
                return self._by_ensembl[ensembl], 'ensembl'
            for mapping in (self._crosswalk_ensembl, self._hgnc_ensembl):
                candidates = mapping.get(ensembl)
                if candidates:
                    on_axis = {self._by_hgnc[hgnc] for hgnc in candidates if hgnc in self._by_hgnc}
                    if len(on_axis) == 1:
                        return next(iter(on_axis)), 'ensembl'
                    if len(on_axis) > 1:
                        return -1, 'ambiguous'
        known = (symbol in self._crosswalk_alias or symbol in self._hgnc_alias
                 or (ensembl and (ensembl in self._by_ensembl
                                  or ensembl in self._crosswalk_ensembl
                                  or ensembl in self._hgnc_ensembl)))
        return (-1, 'offaxis') if known else (-1, 'unmapped')

    def resolve_many(self, symbols):
        """Resolve a symbol vector; returns ``(indices, categories)``."""
        indices = np.empty(len(symbols), dtype=np.int64)
        categories = []
        for position, symbol in enumerate(symbols):
            index, category = self.resolve(symbol)
            indices[position] = index
            categories.append(category)
        return indices, categories

    def loss_summary(self, indices, categories):
        """Counts per alignment category, including the explicit losses."""
        categories = list(categories)
        counts = {name: sum(1 for category in categories if category == name)
                  for name in ALIGNMENT_CATEGORIES}
        resolved = ('canonical', 'crosswalk_alias', 'hgnc_alias', 'ensembl')
        counts['n_source_symbols'] = len(categories)
        counts['n_resolved'] = sum(counts[name] for name in resolved)
        counts['n_axis_genes_resolved'] = len({int(i) for i in indices if i >= 0})
        return counts


def _iter_rows(rows):
    if isinstance(rows, pd.DataFrame):
        for record in rows.to_dict('records'):
            yield record
    else:
        for row in rows:
            yield row if isinstance(row, dict) else dict(row)


def _crosswalk_aliases(rows):
    """source_gene_id -> set(hgnc_id) from the crosswalk's alias/previous rows."""
    mapping = {}
    for row in _iter_rows(rows):
        if str(row.get('match_type')) not in ('alias', 'previous'):
            continue
        key, hgnc = row.get('source_gene_id'), row.get('hgnc_id')
        if key and hgnc:
            mapping.setdefault(str(key), set()).add(str(hgnc))
    return mapping


def _crosswalk_ensembl(rows):
    """Ensembl source_gene_id -> set(hgnc_id) from the crosswalk's ensembl rows."""
    mapping = {}
    for row in _iter_rows(rows):
        if str(row.get('match_type')) != 'ensembl':
            continue
        key = str(row.get('source_gene_id', '') or '').split('.')[0]
        hgnc = row.get('hgnc_id')
        if key and hgnc:
            mapping.setdefault(key, set()).add(str(hgnc))
    return mapping


def _hgnc_aliases(rows):
    """alias_symbol/prev_symbol -> set(hgnc_id) from the pinned HGNC table."""
    mapping = {}
    for row in _iter_rows(rows):
        hgnc = row.get('hgnc_id')
        if not hgnc:
            continue
        for field in ('alias_symbol', 'prev_symbol'):
            for alias in str(row.get(field, '') or '').split('|'):
                if alias:
                    mapping.setdefault(alias, set()).add(str(hgnc))
    return mapping


def _hgnc_ensembl(rows):
    """ensembl_gene_id (version stripped) -> set(hgnc_id) from the pinned HGNC table."""
    mapping = {}
    for row in _iter_rows(rows):
        ensembl = str(row.get('ensembl_gene_id', '') or '').split('.')[0]
        hgnc = row.get('hgnc_id')
        if ensembl and hgnc:
            mapping.setdefault(ensembl, set()).add(str(hgnc))
    return mapping


def rank_scores(gene_index, scores):
    """Rank 1..n over measured genes by descending score, ties on gene_index."""
    gene_index = np.asarray(gene_index, dtype=np.int64)
    scores = np.asarray(scores, dtype=np.float64)
    if len(gene_index) != len(scores):
        raise ValueError('score vector must match the gene index')
    filled = np.where(np.isfinite(scores), scores, -np.inf)
    order = np.lexsort((gene_index, -filled))
    rank = np.empty(len(gene_index), dtype=np.int64)
    rank[order] = np.arange(1, len(gene_index) + 1)
    return rank


def deg_effect(log2fc, padj):
    """Continuous per-contrast DEG effect score ``|log2FC| * -log10(padj)``.

    ``padj`` NA/NaN means the source's test did not report an adjusted p-value
    (DESeq2 independent filtering) and is treated as non-significant: the
    effect score is exactly 0. Very small padj values are floored at
    ``PADJ_FLOOR`` so the score stays finite. Non-finite ``log2FC`` likewise
    scores 0. The score is continuous, so a non-significant contrast cannot be
    promoted above a significant one by fold change alone.
    """
    log2fc = np.asarray(log2fc, dtype=np.float64)
    padj = np.asarray(padj, dtype=np.float64)
    usable = np.isfinite(log2fc) & np.isfinite(padj)
    clipped = np.clip(np.where(usable, padj, 1.0), PADJ_FLOOR, 1.0)
    score = np.where(usable, np.abs(np.where(usable, log2fc, 0.0)) * -np.log10(clipped), 0.0)
    return np.maximum(score, 0.0)


def deg_flag(log2fc, padj):
    """Boolean DEG flag per ``(gene, contrast)`` row: the v2 definition (ADR 0071).

    A row is a DEG when ``padj < PADJ_THRESHOLD`` **and**
    ``|log2FC| >= DEG_LOG2FC_THRESHOLD`` in the same contrast. Both are required:
    a large fold change with a non-significant (or missing) adjusted p-value is
    not a DEG, and neither is a significant small effect.

    ``NA``/non-finite ``padj`` is non-significant — the source's test did not
    report one (DESeq2 independent filtering) — and a non-finite ``log2FC``
    cannot pass the effect-size test. Direction is not filtered on: an
    up-regulated and a down-regulated gene are equally DEGs.
    """
    log2fc = np.asarray(log2fc, dtype=np.float64)
    padj = np.asarray(padj, dtype=np.float64)
    if log2fc.shape != padj.shape:
        raise ValueError('log2FC and padj must have the same shape')
    significant = np.isfinite(padj) & (padj < PADJ_THRESHOLD)
    large_effect = np.isfinite(log2fc) & (np.abs(log2fc) >= DEG_LOG2FC_THRESHOLD)
    return significant & large_effect


def per_gene_deg(gene_index, log2fc, padj, contrast=None):
    """Collapse per-contrast DEG rows to one row per axis gene.

    The continuous score of a gene is the **maximum** over its contrasts
    (declared collapse rule). The reported direction and p-value come from that
    winning contrast; ``n_contrasts_supporting`` counts contrasts with
    ``padj < PADJ_THRESHOLD``.

    The ``deg`` column is the v2 DEG flag (ADR 0071): true when **any** contrast
    passes :func:`deg_flag`, i.e. ``padj < 0.05`` and ``|log2FC| >= 2`` together;
    ``n_contrasts_deg`` counts the contrasts that pass. The winning contrast is
    unaffected by the DEG flag — collapse always keeps the maximum-effect
    contrast, so a gene can be a DEG through one contrast while its reported
    statistics come from another. Ties on the effect score break lexicographically
    by contrast label, then by input order, so the result is deterministic.

    ``contrast`` is optional; when omitted every row is one contrast.
    """
    gene_index = np.asarray(gene_index, dtype=np.int64)
    log2fc = np.asarray(log2fc, dtype=np.float64)
    padj = np.asarray(padj, dtype=np.float64)
    if contrast is None:
        contrast = np.repeat('contrast', len(gene_index))
    else:
        contrast = np.asarray(contrast, dtype=object)
    effect = deg_effect(log2fc, padj)
    frame = pd.DataFrame({'gene_index': gene_index, 'contrast': contrast,
                          'log2fc': log2fc, 'padj': padj, 'effect': effect,
                          'deg': deg_flag(log2fc, padj)})
    frame = frame.sort_values(['gene_index', 'effect', 'contrast'],
                              ascending=[True, False, True], kind='stable')
    best = frame.drop_duplicates('gene_index', keep='first')
    supporting = frame.assign(significant=(frame.padj < PADJ_THRESHOLD))
    counts = supporting.groupby('gene_index', sort=True).significant.sum()
    grouped = frame.groupby('gene_index', sort=True)
    best = best.assign(n_contrasts_supporting=best.gene_index.map(counts).astype(np.int64))
    best = best.assign(deg=best.gene_index.map(grouped.deg.any()).astype(bool),
                       n_contrasts_deg=best.gene_index.map(grouped.deg.sum()).astype(np.int64))
    sign = np.sign(best.log2fc.to_numpy(dtype=np.float64))
    direction = np.where(np.isfinite(sign) & (best.padj < PADJ_THRESHOLD), sign, 0.0)
    best = best.assign(significant=best.padj < PADJ_THRESHOLD,
                       direction=direction.astype(np.int64))
    return best.reset_index(drop=True)


def complete_rank(measured_gene_index, measured_rank, axis_size=AXIS_SIZE):
    """Raw rank vector over the whole axis; unranked genes get ``n + 1``.

    Generic pre-v2 primitive (kept for ``bases.gene_axis`` until its v2 ticket):
    ``measured`` is any set of ``n`` ranked genes and the fill for everything
    else is ``n + 1``. The v2 DEG contract uses :func:`complete_deg_rank`.
    """
    measured_gene_index = np.asarray(measured_gene_index, dtype=np.int64)
    measured_rank = np.asarray(measured_rank, dtype=np.int64)
    n_measured = len(measured_gene_index)
    if n_measured != len(measured_rank):
        raise ValueError('rank vector must match the measured gene index')
    full = np.full(axis_size, n_measured + 1, dtype=np.int64)
    full[measured_gene_index] = measured_rank
    return full


def normalized_rank(rank, n_measured):
    """``rank / n`` for every gene — generic pre-v2 primitive.

    Kept for the pre-v2 fusion reader (``bases.gene_axis``); the v2 DEG contract
    uses :func:`normalized_deg_rank`, whose denominator is the DEG count and
    which is defined for a dataset with no DEGs. ``n_measured`` must be
    positive here, so a zero-count frame is rejected rather than divided by 0.
    """
    rank = np.asarray(rank, dtype=np.float64)
    n_measured = float(n_measured)
    if n_measured <= 0:
        raise ValueError('a source needs at least one measured gene')
    if not np.isfinite(n_measured):
        raise ValueError('n_measured must be finite')
    return rank / n_measured


def deg_only_rank(gene_index, effect, deg):
    """Rank ``1..n_DEG`` over the DEGs only; every non-DEG gets ``n_DEG + 1``.

    The v2 within-dataset rank (ADR 0071 decision 3): only genes flagged ``deg``
    are ordered, by descending effect score with ``gene_index`` as the
    deterministic tie-break. The returned vector aligns with the input rows: a
    DEG carries its rank in ``1..n_DEG`` and a gene that is not a DEG carries
    ``n_DEG + 1``. With no DEGs every row gets ``1`` (``n_DEG + 1``), which is
    the degenerate ordering the informativeness rule reports on.
    """
    gene_index = np.asarray(gene_index, dtype=np.int64)
    effect = np.asarray(effect, dtype=np.float64)
    deg = np.asarray(deg, dtype=bool)
    if not (len(gene_index) == len(effect) == len(deg)):
        raise ValueError('gene index, effect and DEG flag must have the same length')
    n_deg = int(deg.sum())
    rank = np.full(len(gene_index), n_deg + 1, dtype=np.int64)
    if n_deg:
        selected = np.flatnonzero(deg)
        scores = np.where(np.isfinite(effect[selected]), effect[selected], -np.inf)
        order = np.lexsort((gene_index[selected], -scores))
        rank[selected[order]] = np.arange(1, n_deg + 1, dtype=np.int64)
    return rank


def complete_deg_rank(deg_gene_index, deg_rank, n_deg, axis_size=AXIS_SIZE):
    """DEG-only raw rank over the whole frozen axis; non-DEGs get ``n_DEG + 1``.

    ``deg_gene_index``/``deg_rank`` are the dataset's DEGs and their ``1..n_DEG``
    ranks; every other axis gene — measured but not a DEG, and unmeasured alike —
    carries ``n_DEG + 1``, because the v2 rank has no other denominator or
    ordering to offer them (ADR 0071 decision 3).
    """
    deg_gene_index = np.asarray(deg_gene_index, dtype=np.int64)
    deg_rank = np.asarray(deg_rank, dtype=np.int64)
    n_deg = int(n_deg)
    if n_deg < 0:
        raise ValueError('n_deg must be non-negative')
    if len(deg_gene_index) != len(deg_rank) or len(deg_gene_index) != n_deg:
        raise ValueError('the DEG rank must cover exactly n_DEG genes')
    if n_deg and not np.array_equal(np.sort(deg_rank),
                                    np.arange(1, n_deg + 1, dtype=np.int64)):
        raise ValueError('DEG ranks must be exactly 1..n_DEG')
    full = np.full(int(axis_size), n_deg + 1, dtype=np.int64)
    full[deg_gene_index] = deg_rank
    return full


def normalized_deg_rank(rank, n_deg):
    """``rank / n_DEG`` for every gene, per ADR 0071 decision 3.

    A DEG has ``rank`` in ``1..n_DEG``; a gene that is not a DEG has
    ``rank = n_DEG + 1`` (see :func:`complete_deg_rank`), so its normalized
    value is ``(n_DEG + 1) / n_DEG`` — worse than every DEG, never clamped to
    ``1.0``. The denominator is the **DEG** count, never the measured-gene count.

    ``n_DEG == 0`` has no defined normalized rank (the ADR denominator is
    zero), so the value is ``NaN`` for every gene rather than an invented
    finite number. Such a dataset is non-informative and is excluded from
    fusion, so the ``NaN`` is published but never summed (ADR 0068/0071).
    """
    rank = np.asarray(rank, dtype=np.float64)
    n_deg = int(n_deg)
    if n_deg < 0:
        raise ValueError('n_deg must be non-negative')
    if n_deg == 0:
        return np.full(rank.shape, np.nan, dtype=np.float64)
    return rank / float(n_deg)


INFORMATIVENESS_RULE = (
    'a dataset is informative when it has at least one DEG under the v2 '
    'definition and its measured genes carry more than one distinct continuous '
    'effect score; a dataset with no DEG ranks every gene at n_DEG + 1 and a '
    'dataset whose measured effects are all tied ranks by gene_index alone, so '
    'either carries no ordering evidence (ADR 0068, kept by ADR 0071)')


def dataset_informativeness(effect, deg=None):
    """Decide whether a measured dataset carries an informative ordering.

    The narrow, data-driven rule (ADR 0068, kept by ADR 0071): a dataset that
    measures genes but has no DEG, or whose continuous effect scores are all
    tied, has no ordering information. Its published DEG/rank table is still
    valid — every gene then sits at the v2 ``n_DEG + 1`` rank, or the order is
    the declared ``gene_index`` tie-break — but the dataset must not count as a
    fusion source. A dataset with at least one DEG and two or more distinct
    effect scores is informative.

    ``deg`` is the per-gene v2 DEG flag (:func:`deg_flag`); when omitted only the
    ADR 0068 effect test applies, which is kept for callers that have not
    computed the flag.

    Returns ``{'informative', 'n_measured', 'n_deg', 'n_distinct_effect',
    'reason'}`` where ``reason`` is ``None`` for an informative dataset.
    """
    effect = np.asarray(effect, dtype=np.float64)
    if effect.ndim != 1:
        raise ValueError('effect scores must be a vector')
    if effect.size == 0:
        raise ValueError('informativeness needs at least one measured gene')
    if not np.isfinite(effect).all():
        raise ValueError('effect scores must be finite')
    if deg is None:
        flags = np.ones(effect.size, dtype=bool)
        n_deg = None
    else:
        flags = np.asarray(deg, dtype=bool)
        if flags.shape != effect.shape:
            raise ValueError('the DEG flag must align with the effect scores')
        n_deg = int(flags.sum())
    distinct = int(np.unique(effect).size)
    informative = distinct > 1 and bool(flags.any())
    reason = None
    if not informative:
        if n_deg == 0:
            reason = (f'none of the {effect.size} measured genes is a DEG under the '
                      f'v2 definition (padj < {PADJ_THRESHOLD} and '
                      f'|log2FC| >= {DEG_LOG2FC_THRESHOLD}), so every gene ranks at '
                      f'n_DEG + 1 = 1 and the dataset carries no ordering evidence')
        else:
            reason = (f'all {effect.size} measured effect scores are tied at '
                      f'{float(effect[0])!r}; the within-dataset order is the '
                      f'gene_index tie-break, so the dataset carries no ordering '
                      f'evidence')
    return {'informative': bool(informative), 'n_measured': int(effect.size),
            'n_deg': n_deg, 'n_distinct_effect': distinct, 'reason': reason}




ALIGNMENT_CATEGORIES = ('canonical', 'crosswalk_alias', 'hgnc_alias', 'ensembl',
                        'ambiguous', 'offaxis', 'unmapped')
RESOLVED_CATEGORIES = ('canonical', 'crosswalk_alias', 'hgnc_alias', 'ensembl')


def align_table(resolver, rows, symbol_column='symbol'):
    """Attach ``gene_index`` and ``alignment`` to a per-contrast DEG table.

    Unresolved rows are dropped from the returned table but are reported in the
    loss summary, so a caller can publish how many source rows were lost and why.
    A caller must not silently re-map a dropped row.
    """
    rows = pd.DataFrame(rows).copy()
    if symbol_column not in rows:
        raise ValueError(f'missing source symbol column: {symbol_column}')
    indices, categories = resolver.resolve_many(rows[symbol_column].astype(str).tolist())
    resolved = rows.assign(gene_index=indices, alignment=categories)
    report = resolver.loss_summary(indices, categories)
    table = resolved[resolved.gene_index >= 0].copy()
    if table.gene_index.duplicated().any():
        # Duplicate source symbols can alias to one axis gene; the shared
        # max-effect collapse resolves them, so duplicates are expected here.
        pass
    return table.reset_index(drop=True), report


def dataset_deg_rank(rows, axis_size=AXIS_SIZE):
    """Collapse per-contrast rows and publish the within-dataset DEG rank (v2).

    ``rows`` needs ``gene_index``, ``log2fc`` and ``padj`` columns, plus an
    optional ``contrast`` label. Rows with a non-finite ``log2fc`` are not
    measurements and are dropped. Returns ``(table, rank, stats)``:

    * ``table`` — one row per measured gene, sorted by the DEG rank then
      ``gene_index``, carrying the winning contrast's
      ``log2fc``/``padj``/``direction``/``significant``, the ``deg`` flag and
      ``n_contrasts_deg``, the ``n_contrasts``/``n_contrasts_supporting`` counts,
      ``effect`` and the DEG-only ``rank``/``rank_norm``.
    * ``rank`` — a frozen-axis-length frame (``gene_index``, ``rank``,
      ``rank_norm``) covering all ``axis_size`` genes: a DEG at its ``1..n_DEG``
      rank and every other gene at ``n_DEG + 1``.
    * ``stats`` — measured and DEG counts, the retired-vs-new denominators and
      the informativeness verdict.

    ``rank_norm`` is ``rank / n_DEG``; the measured-gene rank and its
    ``n_measured`` denominator are retired (ADR 0071 decision 3). A dataset with
    zero DEGs ranks every gene at ``n_DEG + 1 = 1``, publishes ``rank_norm`` as
    NaN (the ADR denominator is zero) and records itself non-informative rather
    than failing (ADR 0068, ticket #201).
    """
    rows = pd.DataFrame(rows)
    required = {'gene_index', 'log2fc', 'padj'}
    if not required <= set(rows.columns):
        raise ValueError(f'DEG rows require columns {sorted(required)}')
    if rows.empty:
        raise ValueError('a dataset needs at least one DEG row')
    gene_index = rows.gene_index.to_numpy(dtype=np.int64)
    if gene_index.min() < 0 or gene_index.max() >= axis_size:
        raise ValueError('gene_index outside the frozen axis')
    log2fc = rows.log2fc.to_numpy(dtype=np.float64)
    padj = rows.padj.to_numpy(dtype=np.float64)
    measured = np.isfinite(log2fc)
    contrast = (rows['contrast'].astype(str).to_numpy() if 'contrast' in rows
                else None)
    if not measured.any():
        raise ValueError('a dataset needs at least one measured gene')
    table = per_gene_deg(gene_index[measured], log2fc[measured], padj[measured],
                         None if contrast is None else contrast[measured])
    if 'contrast' not in table:
        table = table.assign(contrast='contrast')
    table['n_contrasts'] = table.gene_index.map(
        pd.Series(gene_index[measured]).value_counts()).astype(np.int64)
    n_measured = len(table)
    n_deg = int(table.deg.sum())
    table['rank'] = deg_only_rank(table.gene_index.to_numpy(),
                                  table.effect.to_numpy(),
                                  table.deg.to_numpy())
    table['rank_norm'] = normalized_deg_rank(table['rank'].to_numpy(), n_deg)
    table = table[['gene_index', 'contrast', 'log2fc', 'padj', 'effect',
                   'direction', 'significant', 'deg', 'n_contrasts',
                   'n_contrasts_supporting', 'n_contrasts_deg', 'rank',
                   'rank_norm']]
    table = table.sort_values(['rank', 'gene_index'], kind='stable').reset_index(drop=True)
    degs = table[table.deg]
    full_rank = complete_deg_rank(degs.gene_index.to_numpy(),
                                  degs['rank'].to_numpy(), n_deg,
                                  axis_size=axis_size)
    full_rank_norm = normalized_deg_rank(full_rank, n_deg)
    rank = pd.DataFrame({'gene_index': np.arange(axis_size, dtype=np.int64),
                         'rank': full_rank,
                         'rank_norm': full_rank_norm})
    informativeness = dataset_informativeness(table.effect.to_numpy(),
                                              table.deg.to_numpy())
    stats = {
        'schema': SCHEMA,
        'rank_semantics': RANK_SEMANTICS,
        'axis_size': int(axis_size),
        'n_contrasts': int(len(rows)),
        'n_measured': int(n_measured),
        'n_unmeasured': int(axis_size - n_measured),
        'n_dropped_nonfinite_log2fc': int((~measured).sum()),
        'n_significant': int(table.significant.sum()),
        'n_deg': n_deg,
        'rank_norm_denominator': 'n_deg',
        'nondeg_rank': int(n_deg + 1),
        'informative': informativeness['informative'],
        'informativeness_reason': informativeness['reason'],
        'n_distinct_effect': informativeness['n_distinct_effect'],
        'informativeness_rule': INFORMATIVENESS_RULE,
    }
    return table, rank, stats
