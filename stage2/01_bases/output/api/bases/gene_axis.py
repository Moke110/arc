"""Cross-source HVG/DEG evidence and rank-sum fusion on the frozen gene axis.

This module owns every scientific definition used by the Stage 2 gene-axis
decision-evidence product (ticket #183, spec §C): symbol/alias alignment to the
frozen 18,164-gene axis, Seurat-style normalized dispersion, the continuous DEG
effect score, within-source ranks, rank-sum fusion, overlap, direction
composition, inflection detection and the tenv non-HVG curve.

It is a pure function module: no file is opened and no truth product is ever
consulted. ``execution/build_gene_axis.py`` reads the inputs and writes the
evidence; ``execution/tests/test_gene_axis.py`` exercises these definitions on
fixtures. Keeping the definitions here (ADR 0062) means the driver, the
verifier and the tests cannot drift apart.

Scope boundary: this is decision *preparation*. Nothing here chooses the axis
size ``N``, a per-source quota, the direction rule or the animal-source weight;
those remain #173's decision. The unweighted rank sum is the primary result and
the animal half-weight column is a reported sensitivity only.

Index convention: every public function that consumes per-gene vectors indexes
them by the frozen ``gene_index`` (0..18163). The tracked ``fused_rank.csv`` is
sorted by ``fused_rank``, so its row order is not the axis order; callers scatter
it back with :func:`align_to_gene_index` before using any ``gene_index``-indexed
array. The fixture tests pin this contract so a consumer cannot silently answer
a row-position question.
"""

import numpy as np
import pandas as pd

SCHEMA = 'stage2-gene-axis-evidence-v1'
AXIS_SIZE = 18164
# Declared tenv HVG set sizes read off the published descriptive ranking. The
# non-HVG curve is reported at each of these; none of them is a chosen axis size.
TENV_HVG_KS = (1000, 2000, 5000)
# Reporting convention only (not an axis decision): a per-source "HVG set" is the
# source's top SOURCE_HVG_K measured genes by its own normalized dispersion, so a
# tenv set and a disease-source set are compared at the same declared size.
SOURCE_HVG_K = 2000
HVG_BINS = 20
PADJ_THRESHOLD = 0.05
# -log10(0) is infinite; padj values are floored so the effect score stays finite
# and rank order is deterministic.
PADJ_FLOOR = 1e-300
# Animal MIRI DEG-only sources: species gap-fill, reported with a half-weight
# sensitivity; the primary fusion is unweighted.
ANIMAL_SOURCES = ('gse160516', 'gse108940')

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


def collapse_expression(values, gene_index, axis_size, rule='mean'):
    """Collapse duplicate source rows onto one axis gene (declared rules).

    ``values`` is genes x samples. Duplicate source symbols can resolve to the
    same axis gene through alias evidence; expression rows are averaged
    (``rule='mean'``) so no single alias row wins arbitrarily.
    """
    values = np.asarray(values, dtype=np.float64)
    if rule != 'mean':
        raise ValueError('unsupported expression collapse rule')
    unique, inverse = np.unique(gene_index, return_inverse=True)
    totals = np.zeros((len(unique), values.shape[1]), dtype=np.float64)
    counts = np.zeros(len(unique), dtype=np.int64)
    np.add.at(totals, inverse, values)
    np.add.at(counts, inverse, 1)
    collapsed = totals / counts[:, None]
    return unique, collapsed


def normalized_dispersion(values, gene_index, bins=HVG_BINS):
    """Seurat-style normalized dispersion on a source's own log expression.

    Rows are genes, columns samples. ``log_mean`` is ``log1p(mean)`` of the
    deposited log values; dispersion is sample variance / mean; log dispersion
    is z-scored inside ``bins`` equal-width ``log_mean`` bins (singleton bins
    score 1, an undefined dispersion sorts last). This is the published tenv
    recipe's family (``bases.hvg``) with one deliberate difference: the tenv
    release recovers linear CP10K with ``expm1`` because it is exactly log1p
    CP10K, while the disease matrices are heterogeneous non-invertible log
    scales (voom log2 abundance, log2(RPKM+1), log2(CPM+1), log2(CP10K+1)), so
    the statistics are computed on the deposited log values directly.

    Ties (including undefined scores) break on the frozen ``gene_index``.
    """
    values = np.asarray(values, dtype=np.float64)
    gene_index = np.asarray(gene_index, dtype=np.int64)
    if values.ndim != 2 or values.shape[0] != len(gene_index):
        raise ValueError('expression matrix rows must match the gene index')
    finite = np.isfinite(values)
    count = finite.sum(axis=1)
    safe = np.where(finite, values, 0.0)
    mean = safe.sum(axis=1) / np.maximum(count, 1)
    deviation = np.where(finite, values - mean[:, None], 0.0)
    m2 = (deviation * deviation).sum(axis=1)
    variance = np.where(count > 1, m2 / np.maximum(count - 1, 1), np.nan)
    with np.errstate(invalid='ignore', divide='ignore'):
        dispersion = np.where(mean > 0, variance / mean, np.nan)
        log_dispersion = np.where(dispersion > 0, np.log(dispersion), np.nan)
    log_mean = np.where(mean > 0, np.log1p(np.maximum(mean, 0.0)), np.nan)
    score = np.full(len(gene_index), np.nan, dtype=np.float64)
    usable = np.isfinite(log_mean) & np.isfinite(log_dispersion)
    if usable.any():
        edges = np.linspace(log_mean[usable].min(), log_mean[usable].max(), bins + 1)
        assignment = np.clip(np.searchsorted(edges, log_mean, side='right') - 1, 0, bins - 1)
        for bin_id in range(bins):
            members = usable & (assignment == bin_id)
            if not members.any():
                continue
            block = log_dispersion[members]
            if len(block) == 1:
                score[members] = 1.0
                continue
            spread = block.std(ddof=1)
            if spread > 0:
                score[members] = (block - block.mean()) / spread
            else:
                score[members] = 0.0
    rank = rank_scores(gene_index, score)
    return pd.DataFrame({
        'gene_index': gene_index,
        'n_obs': count,
        'log_mean': log_mean,
        'log_dispersion': log_dispersion,
        'normalized_dispersion': score,
        'rank': rank,
    })


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


def fuse_ranks(ranks, measured, axis_size=AXIS_SIZE, weights=None, flags=None):
    """Unweighted primary rank-sum plus the animal half-weight sensitivity.

    ``ranks`` maps source id -> raw rank vector over the frozen axis.
    ``measured`` maps source id -> its measured-gene count. ``weights`` maps a
    source id to a weight (default 1.0) and is used only for the reported
    sensitivity column. ``flags`` optionally maps source id -> 0/1 per gene for
    "this source reports the gene as a significant DEG".

    The fused rank is ascending ``rank_sum`` with the frozen ``gene_index`` as
    the deterministic tie-break.
    """
    if not ranks:
        raise ValueError('rank fusion needs at least one source')
    axis_size = int(axis_size)
    names = list(ranks)
    primary = np.array([(weights or {}).get(name, 1.0) for name in names], dtype=np.float64)
    sensitivity = np.array([0.5 if name in ANIMAL_SOURCES else (weights or {}).get(name, 1.0)
                            for name in names], dtype=np.float64)
    if not np.isfinite(primary).all() or not np.isfinite(sensitivity).all():
        raise ValueError('source weights must be finite')
    matrix = np.empty((axis_size, len(names)), dtype=np.float64)
    measured_count = np.zeros(axis_size, dtype=np.int64)
    for column, name in enumerate(names):
        rank = np.asarray(ranks[name], dtype=np.int64)
        if rank.shape != (axis_size,):
            raise ValueError(f'rank vector for {name} must cover the frozen axis')
        n_measured = int(measured[name])
        matrix[:, column] = normalized_rank(rank, n_measured)
        measured_count += (rank <= n_measured).astype(np.int64)
    rank_sum = matrix @ primary
    rank_sum_sensitivity = matrix @ sensitivity
    significant = np.zeros(axis_size, dtype=np.int64)
    for name in (flags or {}):
        significant += np.asarray(flags[name], dtype=np.int64)
    gene_index = np.arange(axis_size, dtype=np.int64)
    fused_order = np.lexsort((gene_index, rank_sum))
    fused_rank = np.empty(axis_size, dtype=np.int64)
    fused_rank[fused_order] = np.arange(1, axis_size + 1)
    sensitivity_order = np.lexsort((gene_index, rank_sum_sensitivity))
    sensitivity_rank = np.empty(axis_size, dtype=np.int64)
    sensitivity_rank[sensitivity_order] = np.arange(1, axis_size + 1)
    result = pd.DataFrame({'gene_index': gene_index})
    for name in names:
        result[name] = np.asarray(ranks[name], dtype=np.int64)
    result['rank_sum'] = rank_sum
    result['rank_sum_animal_half'] = rank_sum_sensitivity
    result['fused_rank'] = fused_rank
    result['fused_rank_animal_half'] = sensitivity_rank
    result['n_sources_measured'] = measured_count
    result['n_sources_deg_significant'] = significant
    return result


def fused_order(fused_rank):
    """Gene indices in ascending fused-rank order."""
    fused_rank = np.asarray(fused_rank, dtype=np.int64)
    if sorted(fused_rank.tolist()) != list(range(1, len(fused_rank) + 1)):
        raise ValueError('fused rank must be a complete 1-based permutation')
    return np.argsort(fused_rank, kind='stable')


def direction_composition(fused_rank, direction, max_abs_effect, sizes, tie_counts=None):
    """Top-N up/down composition under two declared, genuinely different rules.

    All vectors are indexed by the frozen ``gene_index`` (never by the row order
    of the tracked fused table, which is sorted by ``fused_rank``):

    ``direction`` is +1/-1/0 per axis gene (0 = the DEG sources report no
    winning effect for the gene); ``max_abs_effect`` is the gene's maximum
    ``|log2FC|`` over the sources that report it; ``tie_counts`` counts genes
    whose winning effect tied across sources.

    - ``symmetric``: the top-N are the first N genes in ascending
      ``fused_rank`` (symmetric inclusion; every gene counts equally).
    - ``max_abs_log2fc``: the top-N are the N genes with the largest
      ``max_abs_effect``, ties broken by ascending ``fused_rank``. It ranks by
      effect magnitude instead of by the fusion, so the two rules select
      different gene sets.

    Genes whose direction is not significant in any source are counted as
    ``n_unsupported`` in both rules.
    """
    fused_rank = np.asarray(fused_rank, dtype=np.int64)
    direction = np.asarray(direction, dtype=np.int64)
    max_abs_effect = np.asarray(max_abs_effect, dtype=np.float64)
    if not (len(fused_rank) == len(direction) == len(max_abs_effect)):
        raise ValueError('direction vectors must cover the frozen axis')
    if sorted(fused_rank.tolist()) != list(range(1, len(fused_rank) + 1)):
        raise ValueError('fused rank must be a complete 1-based permutation')
    tie_counts = (np.zeros(len(direction), dtype=np.int64) if tie_counts is None
                  else np.asarray(tie_counts, dtype=np.int64))
    fused_order = np.argsort(fused_rank, kind='stable')
    effect_order = np.lexsort((fused_rank, -max_abs_effect))
    rows = []
    for rule, order in (('symmetric', fused_order),
                        ('max_abs_log2fc', effect_order)):
        for size in sizes:
            size = int(size)
            if not 0 < size <= len(fused_rank):
                raise ValueError('direction size outside the fused order')
            chosen = order[:size]
            sign = direction[chosen]
            up, down = sign > 0, sign < 0
            unsupported = sign == 0
            rows.append({
                'n': size, 'rule': rule,
                'n_up': int(up.sum()), 'n_down': int(down.sum()),
                'n_unsupported': int(unsupported.sum()),
                'up_share': up.sum() / size, 'down_share': down.sum() / size,
                'unsupported_share': unsupported.sum() / size,
                'n_direction_ties': int(tie_counts[chosen].sum()),
            })
    return pd.DataFrame(rows)


def align_to_gene_index(values, gene_index, axis_size=AXIS_SIZE):
    """Scatter per-row values (rows in any order) onto a frozen-axis vector.

    The tracked fused table is sorted by ``fused_rank``, so its row order is
    *not* the frozen axis order. Any consumer that indexes an array by
    ``gene_index`` must first put the values back into frozen-axis order with
    this helper; using the row position as if it were a gene index is the
    index-space error this module's tests pin down.

    Every axis gene must appear exactly once.
    """
    values = np.asarray(values)
    gene_index = np.asarray(gene_index, dtype=np.int64)
    if values.ndim != 1 or values.shape[0] != gene_index.shape[0]:
        raise ValueError('values must be a vector matching the gene index')
    if len(np.unique(gene_index)) != len(gene_index):
        raise ValueError('gene index has duplicate rows')
    if gene_index.size and (gene_index.min() < 0 or gene_index.max() >= axis_size):
        raise ValueError('gene index outside the frozen axis')
    filled = np.zeros(axis_size, dtype=values.dtype)
    seen = np.zeros(axis_size, dtype=bool)
    filled[gene_index] = values
    seen[gene_index] = True
    if not seen.all():
        raise ValueError('values do not cover the frozen axis')
    return filled


def fused_rank_stats(genes, fused_rank):
    """Min/median/max fused rank of a gene set; ``fused_rank`` is gene-indexed."""
    genes = np.asarray(sorted(set(map(int, genes))), dtype=np.int64)
    positions = (np.asarray(fused_rank, dtype=np.int64)[genes] if len(genes)
                 else np.empty(0, dtype=np.int64))
    return {
        'fused_rank_min': int(positions.min()) if len(positions) else 0,
        'fused_rank_median': float(np.median(positions)) if len(positions) else 0.0,
        'fused_rank_max': int(positions.max()) if len(positions) else 0,
    }


def overlap_metrics(set_a, set_b, fused_rank):
    """Sizes, Jaccard and fused-rank positions of an intersection.

    ``fused_rank`` is indexed by ``gene_index`` (use
    :func:`align_to_gene_index` on the sorted fused table). For a multi-way
    intersection whose two operands are the same set, Jaccard is undefined; use
    :func:`fused_rank_stats` for that case instead of this function.
    """
    set_a, set_b = set(map(int, set_a)), set(map(int, set_b))
    intersection = set_a & set_b
    union = set_a | set_b
    jaccard = len(intersection) / len(union) if union else 0.0
    return {
        'n_a': len(set_a), 'n_b': len(set_b), 'n_intersection': len(intersection),
        'n_union': len(union), 'jaccard': float(jaccard),
        **fused_rank_stats(intersection, fused_rank),
    }


def transcriptome_hvg_sets(hvg_rank, row_gene_index, ks, axis_size=AXIS_SIZE):
    """Declared per-source HVG sets from a source-transcriptome dispersion rank.

    ``hvg_rank`` is a complete 1..n_source_rows ranking over **all** genes the
    source's expression matrix measures (not only the genes that resolve to the
    frozen axis); ``row_gene_index`` maps each source row to its axis gene
    (``-1`` for unresolved rows). For each declared size ``k`` the source's HVG
    set is its top-``k`` source rows, so ``n_total = min(k, n_source_rows)`` and
    ``n_on_axis`` counts the distinct axis genes in that set. This measures how
    much of the source's own variability evidence lands on the frozen axis, not
    how variable the axis is.
    """
    hvg_rank = np.asarray(hvg_rank, dtype=np.int64)
    row_gene_index = np.asarray(row_gene_index, dtype=np.int64)
    if hvg_rank.shape != row_gene_index.shape:
        raise ValueError('hvg rank and row gene index must cover the source rows')
    n_rows = len(hvg_rank)
    result = {}
    for k in ks:
        k = int(k)
        if k <= 0:
            raise ValueError('declared HVG set size must be positive')
        total = min(k, n_rows)
        members = np.flatnonzero(hvg_rank <= total)
        genes = row_gene_index[members]
        on_axis = np.unique(genes[genes >= 0])
        axis_set = np.zeros(axis_size, dtype=bool)
        axis_set[on_axis] = True
        result[k] = {'set': axis_set, 'n_total': int(total),
                     'n_on_axis': int(len(on_axis))}
    return result


def spearman(left, right):
    """Spearman correlation with average ranks for ties (deterministic)."""
    left, right = pd.Series(left), pd.Series(right)
    return float(left.rank(method='average').corr(right.rank(method='average')))


def empty_rank(order):
    """1-based rank vector (indexed by gene_index) for a gene permutation."""
    order = np.asarray(order, dtype=np.int64)
    rank = np.empty(len(order), dtype=np.int64)
    rank[order] = np.arange(1, len(order) + 1)
    return rank


def deg_count_agreement(fused_rank, significant_count, sizes):
    """Agreement between the fused order and the DEG-count order.

    ``fused_rank`` and ``significant_count`` are indexed by ``gene_index``. The
    DEG-count order is descending ``significant_count`` with the **fused rank**
    as tie-break (declared and deterministic): source counts are small integers,
    so ties are frequent and ordering equal-count genes the way the fusion
    orders them keeps the comparison informative. A frozen-axis ``gene_index``
    tie-break would be a different, much weaker ranking (Spearman 0.679 vs
    0.809 on this product); the product's choice is the fused-rank one.

    Returns ``{'spearman', 'top_n_agreement', 'deg_count_rank'}``.
    """
    fused_rank = np.asarray(fused_rank, dtype=np.int64)
    significant_count = np.asarray(significant_count, dtype=np.int64)
    if fused_rank.shape != significant_count.shape:
        raise ValueError('fused rank and DEG count must cover the frozen axis')
    # lexsort sorts by the last key first: descending count is primary, the
    # ascending fused rank is the tie-break. The result is a gene permutation.
    count_order = np.lexsort((fused_rank, -significant_count))
    count_rank = empty_rank(count_order)
    fused_order = np.argsort(fused_rank, kind='stable')
    agreement = {}
    for size in sizes:
        size = int(size)
        if not 0 < size <= len(fused_rank):
            raise ValueError('agreement size outside the frozen axis')
        overlap = len(set(fused_order[:size].tolist())
                      & set(count_order[:size].tolist()))
        agreement[str(size)] = {'n_overlap': int(overlap), 'share': overlap / size}
    return {'spearman': spearman(fused_rank, count_rank),
            'top_n_agreement': agreement,
            'deg_count_rank': count_rank}


def kneedle_knee(x, y, sensitivity=1.0):
    """Online Kneedle: first local maximum of the normalized arc above the chord.

    Both axes are min-max normalized, ``arc = y_norm - x_norm`` and the knee is
    the first interior local maximum of ``arc`` that reaches
    ``sensitivity * max(arc)``. With ``sensitivity=1.0`` that is the first
    interior maximum attaining the global maximum, i.e. the classic Kneedle
    knee. A monotone or straight curve has no interior maximum and returns
    ``None`` — the caller must then report *no inflection* rather than invent
    an ``N``.
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
        # a straight (or monotone-flat) curve's arc is zero up to float noise
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
    """Most concave point of the min-max normalized curve inside a decision range.

    The second difference is taken on the uniform rank grid and the search is
    restricted to the declared range ``[min_n, max_n]``. The leading positions
    below ``min_n`` are excluded because the first few genes of a rank-sum curve
    are strong outliers whose curvature would dominate any knee test; ``min_n``
    is a declared floor (an axis of fewer than 100 genes is not a candidate).
    Returns ``None`` when the largest curvature is below ``minimum_curvature``
    (a straight line has exactly zero curvature, so it reports no candidate
    rather than inventing one).
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
    """Combined Kneedle + second-difference inflection with an explicit answer.

    An inflection is accepted only when both criteria find a candidate, the two
    candidate positions agree within ``agreement_tolerance`` (relative) and the
    Kneedle knee is at least ``minimum_knee_strength`` above the chord. Anything
    else returns ``inflection: False`` and ``n: None``.
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


def tenv_non_hvg_curve(ranked_genes, tenv_rank, ks=TENV_HVG_KS):
    """Share of the fused top-N lying outside each declared tenv HVG set.

    ``ranked_genes`` is the frozen-axis gene order of the fused ranking (the
    ``gene_index`` values of the fused table sorted by ``fused_rank``), **not**
    a row-position vector: ``tenv_rank`` is indexed by ``gene_index``, so a
    permutation of row positions would silently answer the wrong question.
    ``N`` runs over every position 1..axis_size; membership in a declared set of
    size ``K`` is ``tenv_rank <= K``.
    """
    ranked_genes = np.asarray(ranked_genes, dtype=np.int64)
    tenv_rank = np.asarray(tenv_rank, dtype=np.int64)
    if ranked_genes.shape != tenv_rank.shape:
        raise ValueError('fused order and tenv rank must cover the same axis')
    membership = {k: tenv_rank[ranked_genes] <= k for k in ks}
    n = np.arange(1, len(ranked_genes) + 1, dtype=np.int64)
    frame = pd.DataFrame({'n': n})
    for k in ks:
        cumulative = np.cumsum(~membership[k])
        frame[f'n_non_hvg_k{k}'] = cumulative
        frame[f'non_hvg_share_k{k}'] = cumulative / n
    return frame


def hvg_membership(tenv_rank, k):
    """Boolean membership in the published tenv top-K HVG set."""
    return np.asarray(tenv_rank, dtype=np.int64) <= int(k)
