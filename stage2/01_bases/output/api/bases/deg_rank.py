"""Disease-DEG table/rank contract for the Stage 2 gene axis (ADR 0066).

Tickets #187/#188/#189 produce one compact DEG table and one within-dataset rank
per source dataset against this module; the fusion ticket (#190) consumes those
ranks. The contract is deliberately pure and dataset-agnostic:

* **alignment** — a source symbol resolves to one frozen ``gene_index`` through
  canonical symbol, the pinned identity crosswalk's alias/previous evidence, the
  pinned HGNC alias/previous columns, or Ensembl id. Ambiguous, off-axis and
  unmapped symbols are reported, never guessed.
* **effect score** — ``|log2FC| * -log10(padj)`` per contrast. A missing/``NA``
  ``padj`` (e.g. DESeq2 independent filtering) is non-significant and scores 0;
  tiny ``padj`` values are floored so the score stays finite.
* **multi-contrast collapse** — a gene's score is the **maximum** over its
  contrasts; direction and p-value come from that winning contrast.
* **within-dataset rank** — measured genes rank 1..``n_measured`` by descending
  effect score with ``gene_index`` as the deterministic tie-break; a gene the
  dataset does not measure gets ``n_measured + 1``.
* **rank normalization** — ``rank_norm = rank / n_measured``. The denominator is
  the measured-gene count, never the significant-DEG count.

Ranks are within-dataset only; they are never compared across datasets as if
they were effect sizes.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

SCHEMA = 'stage2-disease-deg-v1'
AXIS_SIZE = 18164
PADJ_THRESHOLD = 0.05
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


def per_gene_deg(gene_index, log2fc, padj, contrast=None):
    """Collapse per-contrast DEG rows to one row per axis gene.

    The continuous score of a gene is the **maximum** over its contrasts
    (declared collapse rule). The reported direction and p-value come from that
    winning contrast; ``n_contrasts_supporting`` counts contrasts with
    ``padj < PADJ_THRESHOLD``. Ties on the effect score break lexicographically
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
                          'log2fc': log2fc, 'padj': padj, 'effect': effect})
    frame = frame.sort_values(['gene_index', 'effect', 'contrast'],
                              ascending=[True, False, True], kind='stable')
    best = frame.drop_duplicates('gene_index', keep='first')
    supporting = frame.assign(significant=(frame.padj < PADJ_THRESHOLD))
    counts = supporting.groupby('gene_index', sort=True).significant.sum()
    best = best.assign(n_contrasts_supporting=best.gene_index.map(counts).astype(np.int64))
    sign = np.sign(best.log2fc.to_numpy(dtype=np.float64))
    direction = np.where(np.isfinite(sign) & (best.padj < PADJ_THRESHOLD), sign, 0.0)
    best = best.assign(significant=best.padj < PADJ_THRESHOLD,
                       direction=direction.astype(np.int64))
    return best.reset_index(drop=True)


def complete_rank(measured_gene_index, measured_rank, axis_size=AXIS_SIZE):
    """Raw rank vector over the whole axis; unmeasured genes get n_measured + 1."""
    measured_gene_index = np.asarray(measured_gene_index, dtype=np.int64)
    measured_rank = np.asarray(measured_rank, dtype=np.int64)
    n_measured = len(measured_gene_index)
    if n_measured != len(measured_rank):
        raise ValueError('rank vector must match the measured gene index')
    full = np.full(axis_size, n_measured + 1, dtype=np.int64)
    full[measured_gene_index] = measured_rank
    return full


def normalized_rank(rank, n_measured):
    """Rank / n_measured for measured genes; unmeasured is exactly 1.0."""
    rank = np.asarray(rank, dtype=np.float64)
    n_measured = float(n_measured)
    if n_measured <= 0:
        raise ValueError('a source needs at least one measured gene')
    return np.where(rank <= n_measured, rank / n_measured, 1.0)




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
    """Collapse per-contrast rows and publish the within-dataset DEG rank.

    ``rows`` needs ``gene_index``, ``log2fc`` and ``padj`` columns, plus an
    optional ``contrast`` label. Rows with a non-finite ``log2fc`` are not
    measurements and are dropped. Returns ``(table, rank, stats)``:

    * ``table`` — one row per measured gene, sorted by ``gene_index``, with the
      winning contrast's ``log2fc``/``padj``/``direction``/``significant``,
      ``n_contrasts`` and ``n_contrasts_supporting`` plus ``effect`` and ``rank``.
    * ``rank`` — a frozen-axis-length frame (``gene_index``, ``rank``,
      ``rank_norm``) covering all ``axis_size`` genes, unmeasured at
      ``n_measured + 1``.
    * ``stats`` — measured/unmeasured counts and the rank_norm denominator.
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
    table['rank'] = rank_scores(table.gene_index.to_numpy(),
                                table.effect.to_numpy())
    n_measured = len(table)
    table['rank_norm'] = table['rank'].to_numpy(dtype=np.float64) / n_measured
    table = table.sort_values(['rank', 'gene_index'], kind='stable').reset_index(drop=True)
    full_rank = complete_rank(table.gene_index.to_numpy(), table['rank'].to_numpy(),
                              axis_size=axis_size)
    full_rank_norm = normalized_rank(full_rank, n_measured)
    rank = pd.DataFrame({'gene_index': np.arange(axis_size, dtype=np.int64),
                         'rank': full_rank,
                         'rank_norm': full_rank_norm})
    stats = {
        'schema': SCHEMA,
        'axis_size': int(axis_size),
        'n_contrasts': int(len(rows)),
        'n_measured': int(n_measured),
        'n_unmeasured': int(axis_size - n_measured),
        'n_dropped_nonfinite_log2fc': int((~measured).sum()),
        'n_significant': int(table.significant.sum()),
        'rank_norm_denominator': 'n_measured',
        'unmeasured_rank': int(n_measured + 1),
    }
    return table, rank, stats
