"""Disease-DEG publication and reader for the Stage 2 gene axis (ADR 0066).

Tickets #187/#188/#189 publish, under ``output/disease-deg/``, one compact
Control-vs-Disease table and one within-dataset rank per source dataset; the
fusion ticket (#190) consumes those level-1 inputs. This module is the single
public boundary for that product family:

* ``align_source_table`` — attach the frozen ``gene_index`` to an item-local DEG
  table (symbol first, then source identifier/Ensembl) through
  :class:`cpst_util.derivation.deg_rank.SymbolResolver`; unresolved rows are reported, never
  guessed.
* ``build_dataset_products`` — collapse an item-local per-contrast DEG table with
  the ticket-01 contract (:func:`cpst_util.derivation.deg_rank.dataset_deg_rank`), carry the
  winning contrast's ``stat``/``p_value``/``base_mean``, the pass-any-contrast
  ``deg`` flag and the frozen ``gene_symbol``, and return the compact table, the
  full-axis DEG-only rank frame and the statistics (ADR 0071 v2, republished by
  #208).
* ``genotype_split`` — the genotype-first donor split the Reichart source
  declares (positive donors leave ``nHF``; the single NCCM donor is not used).
* ``contrast_selection`` — the one-joint-design arm collapse: one fit produces
  several contrast rows, each published dataset selects its declared labels.
* ``DiseaseDeg`` — the reader every later ticket imports. It returns per
  category, and per dataset inside it, the compact table, the rank frame and the
  declared eligibility / shared-control record, all digest-checked against the
  category release manifest.
* the ADR 0071 disease contract (#205): :func:`source_deg_rank` derives one
  source's v2 DEG-only rank from its published compact table,
  :func:`level_one_rank` fuses a disease's sources as ``Σ rank_norm``,
  :func:`fuse_disease` adds the disease rank and the ``DEG1000`` cut,
  :func:`validate_catalog` / :func:`build_catalog` carry the
  ``tier``/``evidence_type``/``n_deg``/``deg1000``/``tiering`` shape on the
  category catalog, and :class:`DiseaseFusion` reads the emitted
  ``output/disease-fusion/`` catalog.

Layout (one product directory per disease category)::

    output/disease-deg/<category>/
        catalog.json                 # category -> dataset -> paths + eligibility
        manifest.json                # candidate-then-replace release manifest
        <dataset>.deg.csv            # compact table, one row per measured gene
        <dataset>.rank.csv           # frozen-axis rank frame (18,164 rows)

and one product directory for the per-disease level-1 fusion + ``DEG1000``
(ADR 0071 decision 4)::

    output/disease-fusion/
        catalog.json                 # per disease: counts, ledger, exclusions
        manifest.json                # candidate-then-replace release manifest
        <directory>.deg1000.csv      # gene_index, disease_rank, disease_rank_sum

``.deg.csv`` columns are ``DEG_COLUMNS``; ``.rank.csv`` columns are
``RANK_COLUMNS``; both are asserted by :func:`validate_table` /
:func:`validate_rank`. Ranks are within-dataset only and must never be read as
effect sizes across datasets.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from . import r_env
from .deg_rank import (
    AXIS_SIZE,
    SCHEMA as RANK_SCHEMA,
    SymbolResolver,
    align_table,
    RANK_SEMANTICS,
    complete_deg_rank,
    dataset_deg_rank,
    dataset_informativeness,
    deg_effect,
    deg_only_rank,
    normalized_deg_rank,
)

CATALOG_SCHEMA = 'stage2-disease-deg-catalog-v2'
RELEASE_SCHEMA = 'stage2-disease-deg-release-v2'
LOCAL_SCHEMA = 'stage2-item-deg-v1'

LOCAL_DEG_COLUMNS = ('source_gene_id', 'gene_symbol', 'contrast', 'log2fc', 'stat',
                     'p_value', 'padj', 'base_mean')
#: The published compact-table columns (ADR 0071 v2, #208). ``deg`` is the
#: authoritative pass-any-contrast DEG flag, ``n_contrasts_deg`` counts the
#: contrasts that pass it, ``rank`` is the DEG-only rank (a DEG at ``1..n_DEG``,
#: every other gene at ``n_DEG + 1``) and ``rank_norm = rank / n_DEG``. The
#: retired measured-gene rank and its ``n_measured`` denominator are gone.
DEG_COLUMNS = ('gene_index', 'gene_symbol', 'contrast', 'log2fc', 'stat', 'p_value',
               'padj', 'base_mean', 'direction', 'significant', 'deg',
               'n_contrasts', 'n_contrasts_supporting', 'n_contrasts_deg',
               'effect', 'rank', 'rank_norm')
RANK_COLUMNS = ('gene_index', 'rank', 'rank_norm')

#: ADR 0073: the per-disease level-1 fusion keeps at most the first
#: 1,000 genes of the disease's DEG union (``DEG1000_d``). The cut is a declared
#: limit, never a property of the data: a disease with fewer DEGs keeps all of
#: them and a disease with more publishes exactly 1,000.
DEG1000_LIMIT = 1000
#: The emitted per-disease fusion product family (ADR 0071 decision 4).
FUSION_SCHEMA = 'stage2-disease-fusion-tier-first-v2'
FUSION_RELEASE_SCHEMA = 'stage2-disease-fusion-release-v2'
DEG1000_COLUMNS = ('gene_index', 'disease_rank', 'disease_rank_sum',
                   'source_tier_highest', 'supporting_sources')
UNION_COLUMNS = DEG1000_COLUMNS
#: The declared ``tiering`` record shape the category catalog carries (#205
#: publishes the shape; #206 fills in the roster decision).
TIERING_KEYS = ('tiers_entered', 'tier_union', 'collection_policy')
COLLECTION_POLICY = 'all_qualifying_declared_sources'
SOURCE_LEDGER_KEYS = ('dataset_id', 'tier', 'informative', 'n_deg',
                      'n_deg_in_deg1000')

#: ADR 0071 decision 2: the three declared source tiers (1 = human
#: atrium/ventricle with cardiomyocyte specificity; 2 = hiPSC/hESC-CM, organoids
#: or engineered tissue; 3 = human atrium/ventricle bulk tissue that is not
#: cardiomyocyte-specific).
TIER_VALUES = (1, 2, 3)
#: The declared clinical (non-genetic) disease categories (ADR 0071 decision 2).
#: Their tier-2 evidence is the clinical, non-genetic form only: a
#: mutation-modelled iPSC-CM record belongs to ``genetic_cardiomyopathy`` and is
#: rejected in one of these categories.
CLINICAL_CATEGORIES = ('nHF', 'ICM', 'AF')
#: The tier-stop record's declared criterion marker (#206): the union compared
#: with the limit is resolved to the frozen 18,164-gene axis.
TIER_CRITERION = 'axis_resolved_deg_union'
#: The declared tier ladder and the clinical tier-2 basis record a catalog (or a
#: registry declaration) may carry (#206).
TIER_LADDER_KEYS = ('tier', 'sources')

#: ADR 0073: record every declared tier regardless of cumulative DEG yield.
TIER_COLLECTION_RULE = (
    'Collect all qualifying declared sources in ascending tier order; measure the '
    'complete cumulative disease DEG union on the frozen axis (axis-resolved, '
    'not genome-wide) after each tier. No DEG-count stopping condition applies. '
    'Missing tiers are absent, never back-filled with invented sources.')

GENOTYPE_GROUPS = ('dcm_genotype_negative', 'dcm_genotype_positive', 'arvc',
                   'normal', 'nccm')


def read_axis(path):
    """Read the frozen 18,164-gene identity table (``gene-identity/genes.csv``)."""
    axis = pd.read_csv(path, dtype={'symbol': str, 'hgnc_id': str,
                                    'ensembl_id': str, 'match_type': str})
    required = {'gene_index', 'symbol', 'hgnc_id', 'ensembl_id'}
    if not required <= set(axis.columns):
        raise ValueError(f'axis table requires {sorted(required)}')
    if len(axis) != AXIS_SIZE:
        raise ValueError(f'axis must hold {AXIS_SIZE} genes, found {len(axis)}')
    if axis.gene_index.tolist() != list(range(AXIS_SIZE)):
        raise ValueError('axis gene_index is not the frozen 0..18163 order')
    return axis


def build_resolver(axis, crosswalk=()):
    """The ticket-01 resolver over the frozen axis and pinned crosswalk evidence.

    The pinned HGNC complete set is deliberately not an input here: it is a
    gitignored external source, so using it would make the published tables
    depend on whether that tree happens to be materialised. Canonical symbol,
    Ensembl id and the tracked crosswalk's alias/previous rows are enough for
    this ticket's sources, and every loss is reported.
    """
    return SymbolResolver(axis, crosswalk=crosswalk)


def read_crosswalk(path):
    if path is None or not Path(path).is_file():
        return pd.DataFrame(columns=['dataset_id', 'source_gene_id', 'hgnc_id',
                                     'match_type'])
    return pd.read_csv(path, dtype=str)


def read_local_deg(path, aliases=None):
    """Read one item-local full DEG table (``local/deg/*_deg.csv``).

    ``aliases`` maps an item's own column names onto the canonical contract; the
    pre-#187 HF items deposit ``log2FC`` and ``mean_expr`` rather than
    ``log2fc``/``base_mean``, and a table with no separate identifier column
    falls back to its gene symbol. A table that is empty or still lacks a
    required column is rejected, never silently repaired.
    """
    frame = pd.read_csv(path)
    if aliases:
        frame = frame.rename(columns=dict(aliases))
    if 'source_gene_id' not in frame and 'gene_symbol' in frame:
        frame = frame.assign(source_gene_id=frame['gene_symbol'])
    missing = set(LOCAL_DEG_COLUMNS) - set(frame.columns)
    if missing:
        raise ValueError(f'{path}: local DEG table lacks {sorted(missing)}')
    if frame.empty:
        raise ValueError(f'{path}: local DEG table is empty')
    return frame


def align_source_table(rows, resolver):
    """Attach ``gene_index`` to item-local DEG rows, symbol first then identifier.

    A row's symbol resolves through the canonical/alias/previous/Ensembl
    precedence; when the symbol is empty or unresolved, the source identifier is
    tried (an Ensembl id resolves through the axis or crosswalk). The returned
    table keeps only resolved rows and carries ``alignment``; the report counts
    every category, so an unresolved row is visible instead of silently dropped.
    """
    rows = pd.DataFrame(rows)
    if rows.empty:
        raise ValueError('an item DEG table needs at least one row')
    symbols = rows['gene_symbol'].astype(str).to_numpy(dtype=object)
    identifiers = rows['source_gene_id'].astype(str).to_numpy(dtype=object)
    indices, categories = resolver.resolve_many(symbols.tolist())
    fallback_index, fallback_category = resolver.resolve_many(identifiers.tolist())
    resolved = np.asarray(indices, dtype=np.int64)
    category = np.asarray(categories, dtype=object)
    fallback_resolved = np.asarray(fallback_index, dtype=np.int64) >= 0
    take = (resolved < 0) & fallback_resolved & (category != 'ambiguous')
    resolved[take] = np.asarray(fallback_index, dtype=np.int64)[take]
    category[take] = np.asarray(fallback_category, dtype=object)[take]
    aligned = rows.assign(gene_index=resolved, alignment=category)
    first_pass = resolver.loss_summary(indices, categories)
    final = resolver.loss_summary(resolved, category)
    report = {
        'n_source_rows': int(len(rows)),
        'n_resolved': int((resolved >= 0).sum()),
        'n_axis_genes_resolved': int(len(set(resolved[resolved >= 0].tolist()))),
        'alignment_categories': {name: int(sum(1 for value in category if value == name))
                                 for name in sorted(set(category.tolist()))},
        # A symbol that is ambiguous under the crosswalk stays ambiguous: the
        # identifier fallback never silently re-labels it.
        'n_recovered_by_identifier': int(take.sum()),
        'n_ambiguous_kept': int(sum(1 for value in category if value == 'ambiguous')),
        'symbol_first_pass': first_pass,
    }
    table = aligned[aligned.gene_index >= 0].copy()
    table['gene_index'] = table['gene_index'].astype(np.int64)
    return table.reset_index(drop=True), report


def winning_contrast_rows(aligned):
    """One row per (gene_index, contrast): the maximum-effect duplicate wins.

    This is the same rule ``deg_rank.per_gene_deg`` applies across contrasts; it
    is applied here within one contrast so two source identifiers that alias to
    one axis gene cannot emit two different statistics for the same cell of the
    published table.
    """
    aligned = pd.DataFrame(aligned)
    if aligned.empty:
        raise ValueError('no aligned rows')
    work = aligned.assign(effect=deg_effect(aligned.log2fc, aligned.padj))
    work = work.sort_values(['gene_index', 'contrast', 'effect'],
                            ascending=[True, True, False], kind='stable')
    return work.drop_duplicates(['gene_index', 'contrast'], keep='first')


def build_dataset_products(local_deg, resolver, axis_symbols, dataset_id,
                           axis_size=AXIS_SIZE):
    """Build one dataset's compact table, full-axis rank and statistics.

    ``local_deg`` is an item-local per-contrast table (``LOCAL_DEG_COLUMNS``);
    ``resolver``/``axis_symbols`` bind the frozen axis. Returns
    ``(table, rank, stats, alignment_report)`` where ``table`` follows
    ``DEG_COLUMNS`` sorted by the DEG-only rank and ``rank`` follows
    ``RANK_COLUMNS`` over the whole frozen axis (a DEG at ``1..n_DEG``, every
    other gene at ``n_DEG + 1``; ADR 0071 decision 3, #208).
    """
    aligned, report = align_source_table(local_deg, resolver)
    stats_rows = winning_contrast_rows(aligned)
    table, rank, stats = dataset_deg_rank(
        aligned[['gene_index', 'log2fc', 'padj', 'contrast']], axis_size=axis_size)
    keep = stats_rows[['gene_index', 'contrast', 'stat', 'p_value', 'base_mean']]
    if keep.duplicated(['gene_index', 'contrast']).any():
        raise ValueError(f'{dataset_id}: duplicate winning contrast rows')
    table = table.merge(keep, on=['gene_index', 'contrast'], how='left',
                        validate='one_to_one')
    symbols = pd.Series(np.asarray(axis_symbols, dtype=object),
                        index=np.arange(axis_size, dtype=np.int64))
    table['gene_symbol'] = table.gene_index.map(symbols)
    table = table[list(DEG_COLUMNS)]
    table['deg'] = table['deg'].astype(bool)
    table['n_contrasts_deg'] = table['n_contrasts_deg'].astype(np.int64)
    table = table.sort_values(['rank', 'gene_index'], kind='stable').reset_index(drop=True)
    stats = dict(stats)
    stats.update({'dataset_id': dataset_id, 'alignment': report,
                  'n_contrasts_published': int(table.contrast.nunique())})
    return table, rank, stats, report


def validate_table(table, axis_size=AXIS_SIZE):
    """Validate a published compact DEG table against the ADR 0071 v2 contract.

    The table must carry the v2 :data:`DEG_COLUMNS`, one row per measured gene in
    the frozen axis, the pass-any-contrast ``deg`` flag (matching the published
    ``n_contrasts_deg``), the DEG-only rank (a DEG at ``1..n_DEG``, every other
    gene at ``n_DEG + 1``) and ``rank_norm = rank / n_DEG``. The effect score and
    the rank are re-derived from the published ``log2fc``/``padj``/``deg`` so a
    tampered table cannot pass; a zero-DEG dataset is representable (every rank is
    ``1`` and ``rank_norm`` is ``NaN``), because it is still published as a
    non-informative source (ADR 0068).
    """
    table = pd.DataFrame(table)
    if tuple(table.columns) != DEG_COLUMNS:
        raise ValueError(f'DEG table columns differ from the contract: {list(table.columns)}')
    if table.empty:
        raise ValueError('DEG table is empty')
    if not table.gene_index.is_unique:
        raise ValueError('DEG table repeats a gene_index')
    if table.gene_index.min() < 0 or table.gene_index.max() >= axis_size:
        raise ValueError('DEG table gene_index outside the frozen axis')
    deg = table['deg'].to_numpy(dtype=bool)
    n_contrasts_deg = table.n_contrasts_deg.to_numpy(dtype=np.int64)
    if (n_contrasts_deg < 0).any():
        raise ValueError('DEG table has a negative DEG-contrast count')
    if not np.array_equal(deg, n_contrasts_deg > 0):
        raise ValueError('DEG table deg flag disagrees with n_contrasts_deg')
    if (table.n_contrasts.to_numpy(dtype=np.int64) < n_contrasts_deg).any():
        raise ValueError('DEG table n_contrasts is below n_contrasts_deg')
    n_deg = int(deg.sum())
    ranks = table['rank'].to_numpy(dtype=np.int64)
    rank_norm = table['rank_norm'].to_numpy(dtype=np.float64)
    if n_deg:
        expected = np.concatenate((np.arange(1, n_deg + 1, dtype=np.int64),
                                   np.full(len(table) - n_deg, n_deg + 1,
                                           dtype=np.int64)))
        if not np.array_equal(np.sort(ranks), expected):
            raise ValueError('DEG table rank frame is not 1..n_DEG plus n_DEG+1')
        if not np.allclose(rank_norm, ranks / float(n_deg)):
            raise ValueError('DEG table rank_norm is not rank / n_DEG')
    else:
        if not (ranks == 1).all():
            raise ValueError('a zero-DEG table must rank every gene at 1')
        if not np.isnan(rank_norm).all():
            raise ValueError('a zero-DEG table has no finite rank_norm')
    effect = deg_effect(table.log2fc.to_numpy(dtype=np.float64),
                        table.padj.to_numpy(dtype=np.float64))
    if not np.allclose(table.effect.to_numpy(dtype=np.float64), effect):
        raise ValueError('DEG table effect is not |log2FC| * -log10(padj)')
    if not np.array_equal(ranks, deg_only_rank(table.gene_index.to_numpy(dtype=np.int64),
                                               effect, deg)):
        raise ValueError('DEG table is not ranked by effect with the DEG-only rank')
    if table.gene_index.to_numpy(dtype=np.int64).tolist() != \
            table.sort_values(['rank', 'gene_index'], kind='stable').gene_index.tolist():
        raise ValueError('DEG table is not sorted by the DEG-only rank')
    return table


def _rank_frame_n_deg(ranks, rank_norm):
    """Recover ``n_DEG`` from a published v2 rank frame and validate it.

    A v2 frame carries a DEG's ``1..n_DEG`` rank and ``n_DEG + 1`` for every other
    gene, so the maximum rank is ``n_DEG + 1`` unless every axis gene is a DEG
    (``max = n_DEG``). Both readings are tested against the published
    ``rank_norm = rank / n_DEG`` and exactly one is consistent.
    """
    ranks = np.asarray(ranks, dtype=np.int64)
    rank_norm = np.asarray(rank_norm, dtype=np.float64)
    if ranks.size == 0 or np.any(ranks < 1):
        raise ValueError('rank frame has a rank below 1')
    maximum = int(ranks.max())
    for candidate in (maximum - 1, maximum):
        if candidate < 0:
            continue
        if candidate == 0:
            if bool((ranks == 1).all()) and bool(np.isnan(rank_norm).all()):
                return 0
            continue
        expected = np.concatenate((np.arange(1, candidate + 1, dtype=np.int64),
                                   np.full(ranks.size - candidate, candidate + 1,
                                           dtype=np.int64)))
        if expected.size != ranks.size:
            continue
        if np.array_equal(np.sort(ranks), expected) and np.allclose(
                rank_norm, ranks / float(candidate)):
            return int(candidate)
    raise ValueError('rank frame is not a v2 DEG-only 1..n_DEG / n_DEG+1 assignment')


def validate_rank(rank, axis_size=AXIS_SIZE):
    """Validate a published full-axis v2 rank frame (ADR 0071 decision 3).

    One row per frozen-axis gene in ``0..axis_size-1`` order, a DEG at ``1..n_DEG``
    and every other gene at ``n_DEG + 1``, with ``rank_norm = rank / n_DEG``. A
    zero-DEG frame (all ranks ``1``, ``rank_norm`` NaN) is valid but
    non-informative.
    """
    rank = pd.DataFrame(rank)
    if tuple(rank.columns) != RANK_COLUMNS:
        raise ValueError(f'rank frame columns differ from the contract: {list(rank.columns)}')
    if len(rank) != axis_size:
        raise ValueError(f'rank frame must hold {axis_size} genes, found {len(rank)}')
    if rank.gene_index.tolist() != list(range(axis_size)):
        raise ValueError('rank frame is not the frozen 0..18163 order')
    _rank_frame_n_deg(rank['rank'].to_numpy(dtype=np.int64),
                      rank.rank_norm.to_numpy(dtype=np.float64))
    return rank


# ---------------------------------------------------------------- v2 disease contract

def source_deg_rank(table, axis_size=AXIS_SIZE):
    """Derive one source's v2 DEG-only rank from its published compact table.

    ADR 0071 decisions 3 and 4: a dataset's DEGs rank ``1..n_DEG`` by descending
    effect with ``gene_index`` tie-break, every other frozen-axis gene takes
    ``n_DEG + 1``, and ``rank_norm = rank / n_DEG``. ``table`` is one published
    ``<dataset>.deg.csv`` read back as a frame; the function accepts both
    published table generations so the same fusion runs before and after the
    ticket-8 republish:

    * a **v2 table** carries the published ``deg`` flag (pass-any-contrast under
      the v2 definition) and, normally, the collapsed ``effect``. The published
      flag is authoritative: the rank is rebuilt over exactly those DEGs.
    * a **v1 table** carries one collapsed winning-contrast row per measured gene
      and no ``deg`` flag. The v2 flag is then derived from that row with
      :func:`cpst_util.derivation.deg_rank.dataset_deg_rank` (the one live implementation), i.e.
      the collapse rule applied to the winning contrast the table already stores.
      That convention can under-count a gene whose non-winning contrast passed
      the threshold; the ticket-8 republish replaces it with the published flag.
      ``flag_source`` reports which convention produced the rank.

    Returns ``{'rank', 'n_deg', 'n_measured', 'informative', 'reason',
    'n_distinct_effect', 'flag_source'}``. ``rank`` is the full frozen-axis
    integer frame (length ``axis_size``); a zero-DEG dataset ranks every gene at
    ``1`` and is flagged non-informative (ADR 0068).
    """
    table = pd.DataFrame(table)
    required = {'gene_index', 'log2fc', 'padj'}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f'published DEG table lacks {sorted(missing)}')
    if table.empty:
        raise ValueError('a published DEG table needs at least one measured gene')
    gene_index = table.gene_index.to_numpy(dtype=np.int64)
    if gene_index.min() < 0 or gene_index.max() >= axis_size:
        raise ValueError('published DEG table gene_index outside the frozen axis')
    if len(np.unique(gene_index)) != len(gene_index):
        raise ValueError('published DEG table repeats a gene_index')
    log2fc = table.log2fc.to_numpy(dtype=np.float64)
    padj = table.padj.to_numpy(dtype=np.float64)
    if 'deg' in table.columns:
        deg = table['deg'].to_numpy(dtype=bool)
        if 'effect' in table.columns:
            effect = table['effect'].to_numpy(dtype=np.float64)
        else:
            effect = deg_effect(log2fc, padj)
        if not np.isfinite(effect).all():
            raise ValueError('a v2 published DEG table has a non-finite effect score')
        rank = deg_only_rank(gene_index, effect, deg)
        n_deg = int(deg.sum())
        full = complete_deg_rank(gene_index[deg], rank[deg], n_deg, axis_size=axis_size)
        verdict = dataset_informativeness(effect, deg)
        flag_source = 'published_deg_flag'
    else:
        columns = ['gene_index', 'log2fc', 'padj']
        if 'contrast' in table.columns:
            columns.append('contrast')
        _compact, frame, stats = dataset_deg_rank(table[columns], axis_size=axis_size)
        full = frame['rank'].to_numpy(dtype=np.int64)
        n_deg = int(stats['n_deg'])
        verdict = {'informative': bool(stats['informative']),
                   'reason': stats['informativeness_reason'],
                   'n_distinct_effect': int(stats['n_distinct_effect'])}
        flag_source = 'collapsed_winning_contrast'
    return {'rank': full, 'n_deg': n_deg, 'n_measured': int(len(table)),
            'informative': bool(verdict['informative']),
            'reason': verdict['reason'],
            'n_distinct_effect': int(verdict['n_distinct_effect']),
            'flag_source': flag_source}


def _deg_rank_frame(rank, n_deg, axis_size, label):
    """Validate one in-memory v2 DEG rank frame; return it as int64."""
    rank = np.asarray(rank, dtype=np.int64)
    n_deg = int(n_deg)
    if rank.shape != (axis_size,):
        raise ValueError(f'{label}: a DEG rank frame must cover the frozen axis')
    if n_deg <= 0:
        raise ValueError(
            f'{label}: a dataset with n_DEG = 0 has no rank_norm (the ADR 0071 '
            'denominator is zero); it is non-informative and must be excluded '
            'from the fusion, never ranked')
    if n_deg > axis_size:
        raise ValueError(f'{label}: n_DEG exceeds the frozen axis size')
    if rank.min() < 1 or rank.max() > n_deg + 1:
        raise ValueError(f'{label}: a DEG rank must be in 1..n_DEG+1')
    if int((rank <= n_deg).sum()) != n_deg:
        raise ValueError(f'{label}: a DEG rank frame is not a complete 1..n_DEG assignment')
    if not np.array_equal(np.sort(rank[rank <= n_deg]),
                          np.arange(1, n_deg + 1, dtype=np.int64)):
        raise ValueError(f'{label}: DEG ranks must be exactly 1..n_DEG')
    return rank


def level_one_rank(rank_frames, n_deg, gene_index=None, axis_size=None):
    """Level 1: a disease's source rank from its datasets' DEG-only ranks.

    ADR 0071 decision 4::

        disease_rank_sum = sum(rank_norm)   over the disease's informative sources
        disease rank     = ascending position of that sum, ties by gene_index

    ``rank_frames`` are the full frozen-axis DEG rank frames
    (:func:`source_deg_rank`), i.e. a DEG at ``1..n_DEG`` and every other gene at
    ``n_DEG + 1``; ``n_deg`` is the matching DEG count per frame.
    :func:`cpst_util.derivation.deg_rank.normalized_deg_rank` supplies ``rank / n_DEG``, so a
    non-DEG contributes ``(n_DEG + 1) / n_DEG`` — worse than every DEG of that
    source and never clamped to 1.0. The denominator is the DEG count, never the
    measured-gene count (ADR 0071 decision 3).

    A zero-DEG source cannot be normalized and is rejected: the caller must
    exclude it as non-informative and record the exclusion (ADR 0068). Every
    returned rank is a complete ``1..axis_size`` permutation over the frozen
    axis; the union mask is true where at least one source ranks the gene as a
    DEG.

    Returns ``{'sum', 'rank', 'deg', 'n_sources'}``.
    """
    frames = [np.asarray(frame, dtype=np.int64) for frame in rank_frames]
    counts = [int(value) for value in n_deg]
    if not frames:
        raise ValueError('a disease needs at least one informative source to fuse')
    if len(frames) != len(counts):
        raise ValueError('every source rank frame needs its n_DEG')
    axis_size = int(axis_size if axis_size is not None else len(frames[0]))
    gene_index = (np.arange(axis_size, dtype=np.int64) if gene_index is None
                  else np.asarray(gene_index, dtype=np.int64))
    if gene_index.shape != (axis_size,):
        raise ValueError('gene_index must cover the frozen axis')
    total = np.zeros(axis_size, dtype=np.float64)
    union = np.zeros(axis_size, dtype=bool)
    for position, (frame, count) in enumerate(zip(frames, counts)):
        frame = _deg_rank_frame(frame, count, axis_size, f'source {position}')
        total += normalized_deg_rank(frame, count)
        union |= frame <= count
    order = np.lexsort((gene_index, total))
    rank = np.empty(axis_size, dtype=np.int64)
    rank[order] = np.arange(1, axis_size + 1, dtype=np.int64)
    return {'sum': total, 'rank': rank, 'deg': union, 'n_sources': len(frames)}


DISEASE_MIN_DEGS = 100


def disease_eligibility(n_deg, published=True):
    """ADR 0072 admission after axis-resolved disease fusion, never per source.

    ``n_deg`` must be the distinct frozen-axis union measured by ``fuse_disease``,
    not a sum of source counts or the truncated DEG1000 count.
    """
    if not published:
        return {'axis_resolved_deg_union': None, 'eligible': False,
                'exclusion_reason': 'declared disease has no published DEG evidence'}
    count = _positive_int(n_deg, 'axis-resolved disease DEG union')
    eligible = count >= DISEASE_MIN_DEGS
    return {'axis_resolved_deg_union': count, 'eligible': eligible,
            'exclusion_reason': (None if eligible else
                f'axis-resolved DEG union {count} is below {DISEASE_MIN_DEGS} '
                '(ADR 0072); upstream evidence and DEG1000 are retained')}


def validate_disease_eligibility(record, required=False):
    """Reject partial or contradictory admission metadata; accept legacy records."""
    keys = ('axis_resolved_deg_union', 'eligible', 'exclusion_reason')
    if not required and not any(key in record for key in keys):
        return
    expected = disease_eligibility(record['n_deg'], record.get('published', True))
    if any(key not in record for key in keys):
        raise ValueError('disease eligibility metadata is incomplete')
    if not isinstance(record['eligible'], bool):
        raise ValueError('disease eligibility must be boolean')
    if record['axis_resolved_deg_union'] != expected['axis_resolved_deg_union']:
        raise ValueError('disease eligibility union differs from measured n_deg')
    if record['eligible'] != expected['eligible']:
        raise ValueError('disease eligibility disagrees with the 100-DEG union gate')
    reason = record['exclusion_reason']
    if record['eligible'] and reason is not None:
        raise ValueError('an eligible disease cannot have an exclusion reason')
    if not record['eligible'] and (not isinstance(reason, str) or not reason):
        raise ValueError('an ineligible disease needs an explicit exclusion reason')


def fuse_disease(sources, axis_size=AXIS_SIZE, limit=DEG1000_LIMIT,
                 gene_index=None):
    """Fuse one disease's sources at level 1 and cut its ``DEG1000`` set.

    ``sources`` is the disease's declared source list; each entry needs
    ``dataset_id``, the source's full-axis DEG rank frame (``rank``) and DEG count
    (``n_deg``), its required ``tier`` and optional ``informative`` plus the recorded
    ``informativeness_reason``.

    ADR 0071 decision 4 / ADR 0068:

    * a source is **excluded** from the sum when it is non-informative
      (``informative is not True``) or has ``n_DEG = 0``; the exclusion and its
      reason are published in ``sources``/``excluded_sources``, never dropped
      silently. A zero-DEG source would make every term ``NaN`` and silently
      remove a gene from every disease that contains it, so it is never summed.
    * ``disease_rank_sum = Σ rank_norm`` over the remaining informative sources;
      the disease rank orders actual highest DEG tier, this sum and ``gene_index``.
    * the **DEG union** is the set of genes that are a DEG in at least one
      informative source; ``DEG1000_d`` is that union intersected with disease
      rank ``<= limit`` (all of it when the union is smaller than the limit).
      The union size is never padded or constrained by the cutoff.
    * a disease with **no informative source** is representable: ``informative``
      is false, ``rank``/``sum`` are ``None``, both masks are empty and the
      record carries the reason. Level 2 must then exclude the disease (ADR 0068)
      rather than divide by a zero DEG count.

    Union genes receive tiers 1..3 and precede every non-DEG (which has no tier).
    The full union therefore occupies contiguous disease ranks 1..n_d, with
    ``deg1000 = min(n_d, limit)``. Numerical rank sums receive no tier weight.

    Returns the fusion record described above; see :func:`deg1000_rows` for the
    published per-gene ``DEG1000`` table.
    """
    limit = int(limit)
    if limit <= 0:
        raise ValueError('the DEG1000 limit must be positive')
    axis_size = int(axis_size)
    prepared, seen = [], set()
    for position, source in enumerate(sources):
        source = dict(source)
        dataset_id = str(source.get('dataset_id') or f'source_{position}')
        if dataset_id in seen:
            raise ValueError(f'{dataset_id}: a disease cannot declare one source twice')
        seen.add(dataset_id)
        tier = source.get('tier')
        tier = _tier(tier, f'{dataset_id}: source tier')
        n_deg = int(source['n_deg'])
        if n_deg < 0:
            raise ValueError(f'{dataset_id}: n_DEG must be non-negative')
        informative = bool(source.get('informative', n_deg > 0))
        if n_deg == 0:
            reason = source.get('informativeness_reason') or (
                'no DEG under the v2 definition: every gene ranks at n_DEG + 1 = 1 '
                'and rank_norm is undefined (n_DEG = 0), so the source is '
                'non-informative and contributes no term')
            prepared.append({'dataset_id': dataset_id, 'tier': tier,
                             'informative': False, 'n_deg': 0, 'reason': reason,
                             'mask': None})
            continue
        if not informative:
            reason = source.get('informativeness_reason') or (
                'declared non-informative: its measured effect scores carry no '
                'ordering evidence (ADR 0068)')
            prepared.append({'dataset_id': dataset_id, 'tier': tier,
                             'informative': False, 'n_deg': n_deg, 'reason': reason,
                             'mask': None})
            continue
        frame = _deg_rank_frame(source['rank'], n_deg, axis_size, dataset_id)
        prepared.append({'dataset_id': dataset_id, 'tier': tier,
                         'informative': True, 'n_deg': n_deg, 'reason': None,
                         'mask': frame <= n_deg, 'rank': frame})
    if not prepared:
        raise ValueError('a disease needs at least one declared source')
    fused = [entry for entry in prepared if entry['informative']]
    if fused:
        level_one = level_one_rank([entry['rank'] for entry in fused],
                                   [entry['n_deg'] for entry in fused],
                                   gene_index=gene_index, axis_size=axis_size)
        total = level_one['sum']
        deg_union = level_one['deg']
        strongest = np.full(axis_size, 4, dtype=np.int64)
        for entry in fused:
            strongest[entry['mask']] = np.minimum(
                strongest[entry['mask']], entry['tier'])
        indices = (np.arange(axis_size, dtype=np.int64) if gene_index is None
                   else np.asarray(gene_index, dtype=np.int64))
        order = np.lexsort((indices, total, strongest))
        disease_rank = np.empty(axis_size, dtype=np.int64)
        disease_rank[order] = np.arange(1, axis_size + 1, dtype=np.int64)
        deg1000_mask = deg_union & (disease_rank <= limit)
    else:
        disease_rank = None
        total = None
        deg_union = np.zeros(axis_size, dtype=bool)
        deg1000_mask = np.zeros(axis_size, dtype=bool)
    source_tier_highest = np.full(axis_size, None, dtype=object)
    if fused:
        source_tier_highest[deg_union] = strongest[deg_union]
    ledger = []
    for entry in prepared:
        hits = (int((entry['mask'] & deg1000_mask).sum())
                if entry['informative'] else 0)
        record = {'dataset_id': entry['dataset_id'], 'tier': entry['tier'],
                  'informative': entry['informative'], 'n_deg': entry['n_deg'],
                  'n_deg_in_deg1000': hits}
        if entry['reason'] is not None:
            record['reason'] = entry['reason']
        ledger.append(record)
    n_deg_union = int(deg_union.sum())
    n_deg1000 = int(deg1000_mask.sum())
    return {
        **disease_eligibility(n_deg_union),
        'informative': bool(fused),
        'n_sources_declared': len(prepared),
        'n_sources_informative': len(fused),
        'n_deg': n_deg_union,
        'deg1000': n_deg1000,
        'limit': limit,
        'rank': disease_rank,
        'sum': total,
        'deg_union': deg_union,
        'source_tier_highest': source_tier_highest,
        'supporting_sources': [sorted(entry['dataset_id'] for entry in fused
                                     if entry['mask'][gene])
                               for gene in range(axis_size)],
        'deg1000_mask': deg1000_mask,
        'sources': ledger,
        'excluded_sources': [{'dataset_id': record['dataset_id'],
                              'tier': record['tier'], 'n_deg': record['n_deg'],
                              'reason': record['reason']}
                             for record in ledger if not record['informative']],
        'reason': (None if fused else
                   'no informative source: the disease has no DEG union and no '
                   'level-1 ordering, so level 2 must exclude it (ADR 0068)'),
    }


def union_rows(fusion):
    """Complete disease DEG evidence, including genes outside the cutoff."""
    genes = np.flatnonzero(fusion['deg_union'])
    if not len(genes):
        return pd.DataFrame(columns=list(UNION_COLUMNS))
    genes = genes[np.argsort(np.asarray(fusion['rank'])[genes], kind='stable')]
    return pd.DataFrame({
        'gene_index': genes.astype(np.int64),
        'disease_rank': np.asarray(fusion['rank'])[genes].astype(np.int64),
        'disease_rank_sum': np.asarray(fusion['sum'])[genes],
        'source_tier_highest': np.asarray(fusion['source_tier_highest'])[genes].astype(np.int64),
        'supporting_sources': [json.dumps(fusion['supporting_sources'][gene],
                                         separators=(',', ':')) for gene in genes],
    }, columns=list(UNION_COLUMNS))


def deg1000_rows(fusion):
    """Tier-first first min(1000,n_d) rows; no padding."""
    return union_rows(fusion).iloc[:fusion['deg1000']].copy()


# ---------------------------------------------------------------- tier rule (#206)

def tier_ladder(datasets):
    """Group a category's declared source records into its ordered tier ladder.

    ADR 0071 decision 2: every source record carries its tier, and the collection
    order is the ascending tier. ``datasets`` is the declared source list (the
    registry records, each carrying ``dataset_id`` and ``tier``); the returned
    ladder is ``[{'tier': 1, 'sources': [...]}, ...]`` with tiers ascending and
    sources in declared order inside a tier. A tier with no declared source is
    absent from the ladder, so AF's ladder is exactly ``[1, 3]`` and the
    TOX's is exactly ``[2]``. Each source appears exactly once.
    """
    if not datasets:
        raise ValueError('a category needs at least one declared source')
    grouped = {}
    seen = set()
    for position, record in enumerate(datasets):
        if not isinstance(record, dict):
            raise ValueError(f'source declaration {position} is not a record')
        dataset_id = record.get('dataset_id')
        if not isinstance(dataset_id, str) or not dataset_id:
            raise ValueError(f'source declaration {position} has no dataset_id')
        if dataset_id in seen:
            raise ValueError(f'{dataset_id}: a category cannot declare one source twice')
        seen.add(dataset_id)
        tier = record.get('tier')
        if not isinstance(tier, int) or isinstance(tier, bool) or tier not in TIER_VALUES:
            raise ValueError(f'{dataset_id}: tier must be an integer in {TIER_VALUES}')
        grouped.setdefault(int(tier), []).append(dataset_id)
    return [{'tier': tier, 'sources': grouped[tier]} for tier in sorted(grouped)]


def _source_deg_genes(entry, label, axis_size=AXIS_SIZE):
    """Validate one measured source DEG set and return it as sorted gene indices."""
    if not isinstance(entry, dict):
        raise ValueError(f'{label}: the measured DEG set must be a record')
    if 'deg' not in entry or 'n_deg' not in entry:
        raise ValueError(f'{label}: the measured DEG set needs deg and n_deg')
    value = np.asarray(entry['deg'])
    if value.dtype == bool:
        if value.shape != (axis_size,):
            raise ValueError(f'{label}: a boolean DEG mask must cover the frozen axis')
        genes = np.flatnonzero(value).astype(np.int64)
    else:
        if value.ndim != 1:
            raise ValueError(f'{label}: deg must be a gene_index vector or an axis mask')
        if not np.issubdtype(value.dtype, np.integer):
            try:
                numeric = value.astype(np.float64)
            except (TypeError, ValueError):
                raise ValueError(f'{label}: deg must be gene_index integers') from None
            if (not np.isfinite(numeric).all()
                    or not np.equal(numeric, np.floor(numeric)).all()):
                raise ValueError(f'{label}: deg must be gene_index integers')
        genes = value.astype(np.int64)
    if genes.size and (genes.min() < 0 or genes.max() >= axis_size):
        raise ValueError(f'{label}: deg holds a gene_index outside the frozen axis')
    if len(np.unique(genes)) != len(genes):
        raise ValueError(f'{label}: deg repeats a gene_index')
    n_deg = entry['n_deg']
    if not isinstance(n_deg, int) or isinstance(n_deg, bool) or n_deg < 0:
        raise ValueError(f'{label}: n_deg must be a non-negative integer')
    if n_deg != len(genes):
        raise ValueError(f'{label}: n_deg ({n_deg}) differs from the measured DEG set '
                         f'({len(genes)})')
    return np.sort(genes)


def validate_clinical_tier_basis(category, record):
    """The ADR 0071 clinical-category tier-2 gate (#206).

    For ``nHF``/``ICM``/``AF`` a tier-2 source must be a model of the clinical,
    non-genetic form: it declares a non-empty ``non_genetic_model_basis`` and it
    is not mutation-modelled (no ``genotype_status == 'pathogenic'`` and no
    ADR 0066 isogenic exception). A mutation-modelled iPSC-CM record belongs to
    ``genetic_cardiomyopathy``, where the same record is accepted.
    """
    if category not in CLINICAL_CATEGORIES:
        return record
    if int(record.get('tier') or 0) != 2:
        return record
    dataset_id = record.get('dataset_id', '?')
    if record.get('genotype_status') == 'pathogenic':
        raise ValueError(
            f'{dataset_id}: a mutation-modelled tier-2 record belongs to '
            'genetic_cardiomyopathy, not to the clinical category '
            f'{category} (ADR 0071 decision 2)')
    if (record.get('isogenic_exception') or {}).get('applies') is True:
        raise ValueError(
            f'{dataset_id}: an isogenic (mutation-modelled) tier-2 record belongs to '
            f'genetic_cardiomyopathy, not to the clinical category {category} '
            '(ADR 0071 decision 2)')
    basis = record.get('non_genetic_model_basis')
    if not isinstance(basis, str) or not basis.strip():
        raise ValueError(
            f'{dataset_id}: a tier-2 record in the clinical category {category} must '
            'declare its non-genetic model basis (non_genetic_model_basis) '
            '(ADR 0071 decision 2)')
    return record


def tier_ledger(category, datasets, source_deg, limit=DEG1000_LIMIT):
    """Collect every declared qualifying source and record cumulative axis unions.

    ``limit`` is only the temporary fusion-cut parameter, independent of source
    collection. ADR 0073 removes all count-based collection stopping conditions.
    Non-informative sources remain in the ledger but contribute no union genes.
    """
    ladder = tier_ladder(datasets)
    declared = [record['dataset_id'] for record in datasets]
    missing = [dataset_id for dataset_id in declared if dataset_id not in source_deg]
    if missing:
        raise ValueError(f'{category}: declared sources without a measured DEG set: '
                         f'{missing}')
    unknown = sorted(set(source_deg) - set(declared))
    if unknown:
        raise ValueError(f'{category}: measured sources that are not declared: {unknown}')
    by_id = {record['dataset_id']: record for record in datasets}
    measured, sources, union = {}, [], set()
    for dataset_id in declared:
        validate_clinical_tier_basis(category, by_id[dataset_id])
        entry = source_deg[dataset_id]
        genes = _source_deg_genes(entry, f'{category}/{dataset_id}')
        informative = bool(entry.get('informative', len(genes) > 0))
        measured[dataset_id] = (genes, informative)
        sources.append({'dataset_id': dataset_id, 'tier': int(by_id[dataset_id]['tier']),
                        'informative': informative, 'n_deg': int(len(genes))})
    entered, tier_union = [], []
    for rung in ladder:
        for dataset_id in rung['sources']:
            genes, informative = measured[dataset_id]
            if informative:
                union.update(int(gene) for gene in genes)
        entered.append(int(rung['tier']))
        tier_union.append({'tier': int(rung['tier']), 'n_deg_union': int(len(union))})
    record = {
        'tiers_entered': entered,
        'tier_union': tier_union,
        'collection_policy': COLLECTION_POLICY,
        'criterion': TIER_CRITERION,
        'rule': TIER_COLLECTION_RULE,
    }
    n_deg = int(len(union))
    _validate_tiering(record, n_deg)
    return {
        'tiering': record,
        'tier_sources': [{'tier': int(rung['tier']),
                          'sources': list(rung['sources'])}
                         for rung in ladder if int(rung['tier']) in entered],
        'n_deg': n_deg,
        'deg1000': min(n_deg, limit),
        'sources': sources,
    }


# ---------------------------------------------------------------- catalog shape

def _walk_strings(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for key, item in value.items():
            yield from _walk_strings(key)
            yield from _walk_strings(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from _walk_strings(item)


def _assert_relative(payload, label):
    """A published catalog may name only repo-relative paths."""
    for text in _walk_strings(payload):
        if text.startswith('/') or (len(text) > 2 and text[1] == ':'
                                    and text[2] in '\\/'):
            raise ValueError(f'{label} contains an absolute path: {text}')


def _positive_int(value, label, allow_zero=True):
    if not isinstance(value, int) or isinstance(value, bool):
        raise ValueError(f'{label} must be an integer')
    if value < 0 or (value == 0 and not allow_zero):
        raise ValueError(f'{label} must be positive' if not allow_zero
                         else f'{label} must be non-negative')
    return value


def _tier(value, label):
    if not isinstance(value, int) or isinstance(value, bool) or not 1 <= value <= 3:
        raise ValueError(f'{label} must be an integer in 1..3')
    return value


def _validate_tiering(tiering, n_deg):
    """Validate the ADR 0073 all-source ledger and full axis-resolved union."""
    if not isinstance(tiering, dict):
        raise ValueError('tiering must be a record')
    missing = [key for key in TIERING_KEYS if key not in tiering]
    if missing:
        raise ValueError(f'tiering lacks {missing}')
    tiers = tiering['tiers_entered']
    if not isinstance(tiers, list) or not tiers:
        raise ValueError('tiering.tiers_entered must be a non-empty list')
    for value in tiers:
        _tier(value, 'tiering.tiers_entered entry')
    if tiers != sorted(set(tiers)):
        raise ValueError('tiering.tiers_entered must be unique and ascending')
    criterion = tiering.get('criterion')
    if criterion is not None and criterion != TIER_CRITERION:
        raise ValueError(f'tiering.criterion must be {TIER_CRITERION!r}')
    rule = tiering.get('rule')
    if rule is not None and (not isinstance(rule, str) or not rule.strip()):
        raise ValueError('tiering.rule must be a non-empty declaration string')
    unions = tiering['tier_union']
    if not isinstance(unions, list) or len(unions) != len(tiers):
        raise ValueError('tiering.tier_union must have one record per entered tier')
    sizes = []
    for expected, record in zip(tiers, unions):
        if not isinstance(record, dict) or set(record) != {'tier', 'n_deg_union'}:
            raise ValueError('a tiering.tier_union record is exactly {tier, n_deg_union}')
        if record['tier'] != expected:
            raise ValueError('tiering.tier_union tiers must match tiers_entered')
        sizes.append(_positive_int(record['n_deg_union'], 'tiering.tier_union size'))
    if sizes != sorted(sizes):
        raise ValueError('the per-tier DEG union must be non-decreasing')
    if sizes[-1] != int(n_deg):
        raise ValueError(
            'the last per-tier DEG union must equal the published DEG union '
            f'({sizes[-1]} != {n_deg})')
    if tiering['collection_policy'] != COLLECTION_POLICY:
        raise ValueError('tiering.collection_policy must collect all qualifying declared sources')
    if {'stopped_at_tier', 'stop_reached', 'limit'} & set(tiering):
        raise ValueError('obsolete count-stopping fields are not a collection contract')
    return tiering


def _validate_tier_sources(tier_sources, tiering, records):
    """Validate the published tier ladder against the ``tiering`` record (#206).

    The ladder is the declaration the ``tiering`` record was computed from, so it
    must name exactly the entered tiers, every published dataset exactly once, and
    each dataset at its own declared tier.
    """
    if not isinstance(tier_sources, list) or not tier_sources:
        raise ValueError('tier_sources must be a non-empty list of tier records')
    declared_tier = {record['dataset_id']: record['tier'] for record in records}
    seen = set()
    tiers = []
    for rung in tier_sources:
        if not isinstance(rung, dict) or set(rung) != set(TIER_LADDER_KEYS):
            raise ValueError('a tier_sources record is exactly {tier, sources}')
        _tier(rung['tier'], 'tier_sources tier')
        tiers.append(rung['tier'])
        sources = rung['sources']
        if not isinstance(sources, list) or not sources:
            raise ValueError(f'tier_sources tier {rung["tier"]} has no sources')
        for dataset_id in sources:
            if not isinstance(dataset_id, str) or not dataset_id:
                raise ValueError('a tier_sources source must be a non-empty string')
            if dataset_id in seen:
                raise ValueError(f'tier_sources repeats the source {dataset_id}')
            if dataset_id not in declared_tier:
                raise ValueError(f'tier_sources names an undeclared source {dataset_id}')
            if int(declared_tier[dataset_id]) != int(rung['tier']):
                raise ValueError(
                    f'{dataset_id}: tier_sources places the source at tier '
                    f'{rung["tier"]} but its record declares tier '
                    f'{declared_tier[dataset_id]}')
            seen.add(dataset_id)
    if tiers != tiering['tiers_entered']:
        raise ValueError('tier_sources tiers must equal the tiering.tiers_entered '
                         'ladder')
    if seen != set(declared_tier):
        raise ValueError('tier_sources must cover every published dataset')


def validate_catalog(catalog, axis_size=AXIS_SIZE):
    """Validate the v2 category-catalog shape (#205; populated by #206/#208).

    On top of ``stage2-disease-deg-catalog-v1`` a category catalog declares
    ``tier`` (1..3), ``evidence_type``, the disease-level ``n_deg`` (DEG-union
    size) and ``deg1000`` (``DEG1000_d`` size), a ``tiering`` record, and per
    dataset record ``n_deg`` / ``n_deg_in_deg1000``.

    The checks are the acceptance contract: the fields are present and typed,
    ``tier`` is in 1..3, ``deg1000 == min(n_deg, 1000)`` (the union occupies the
    first ranks, so the cut is exact), the per-dataset ledger is bounded by both
    its own DEG count and ``deg1000``, the union is not smaller than one source
    nor larger than their sum, every ``DEG1000`` gene is covered by at least one
    source, ``tiering`` agrees with the published union size
    (:func:`_validate_tiering`), a clinical category's tier-2 records pass the
    non-genetic basis gate (:func:`validate_clinical_tier_basis`), an optional
    published ``tier_sources`` ladder agrees with ``tiering`` and the dataset
    records (:func:`_validate_tier_sources`), and no path is absolute.
    """
    if not isinstance(catalog, dict):
        raise ValueError('a category catalog must be a record')
    for key in ('category', 'directory', 'tier', 'evidence_type', 'n_deg',
                'deg1000', 'tiering', 'datasets'):
        if key not in catalog:
            raise ValueError(f'category catalog lacks {key}')
    if not isinstance(catalog['category'], str) or not catalog['category']:
        raise ValueError('category must be a non-empty string')
    if not isinstance(catalog['directory'], str) or not catalog['directory']:
        raise ValueError('directory must be a non-empty string')
    _tier(catalog['tier'], 'category tier')
    if not isinstance(catalog['evidence_type'], str) or not catalog['evidence_type']:
        raise ValueError('evidence_type must be a non-empty string')
    _positive_int(catalog.get('axis_size', axis_size), 'axis_size', allow_zero=False)
    validate_disease_eligibility(catalog)
    n_deg = _positive_int(catalog['n_deg'], 'n_deg')
    deg1000 = _positive_int(catalog['deg1000'], 'deg1000')
    if deg1000 != min(n_deg, DEG1000_LIMIT):
        raise ValueError(
            f'deg1000 must be min(n_deg, {DEG1000_LIMIT}) because the DEG union '
            f'occupies the first ranks ({deg1000} != min({n_deg}, {DEG1000_LIMIT}))')
    records = catalog['datasets']
    if not isinstance(records, list) or not records:
        raise ValueError('datasets must be a non-empty list of records')
    seen, counts, hits = set(), [], 0
    informative_sources = 0
    for record in records:
        if not isinstance(record, dict):
            raise ValueError('a dataset record must be a record')
        dataset_id = record.get('dataset_id')
        if not isinstance(dataset_id, str) or not dataset_id:
            raise ValueError('a dataset record has no dataset_id')
        if dataset_id in seen:
            raise ValueError(f'duplicate dataset id {dataset_id}')
        seen.add(dataset_id)
        _tier(record.get('tier'), f'{dataset_id}: tier')
        if not isinstance(record.get('informative'), bool):
            raise ValueError(f'{dataset_id}: informative must be boolean')
        source_n = _positive_int(record.get('n_deg'), f'{dataset_id}: n_deg')
        source_hits = _positive_int(record.get('n_deg_in_deg1000'),
                                    f'{dataset_id}: n_deg_in_deg1000')
        if source_hits > min(source_n, deg1000):
            raise ValueError(
                f'{dataset_id}: n_deg_in_deg1000 cannot exceed min(n_deg, deg1000)')
        if record['informative']:
            if source_n == 0:
                raise ValueError(f'{dataset_id}: an informative source needs a DEG')
            informative_sources += 1
            counts.append(source_n)
            hits += source_hits
        elif source_hits:
            raise ValueError(
                f'{dataset_id}: a non-informative source is excluded from the '
                'fusion and cannot land in DEG1000')
    if counts and n_deg < max(counts):
        raise ValueError('the DEG union cannot be smaller than one informative '
                         'source DEG set')
    if n_deg > sum(counts):
        raise ValueError('the DEG union cannot exceed the sum of the informative '
                         'source DEG counts')
    if hits < deg1000:
        raise ValueError('every DEG1000 gene must be a DEG in at least one '
                         'informative source')
    if n_deg and not informative_sources:
        raise ValueError('a non-empty DEG union needs an informative source')
    _validate_tiering(catalog['tiering'], n_deg)
    # #206: the dataset records' own tiers are the declared ladder, so they must
    # be exactly the entered tiers (an AF-style gap is representable because a
    # missing tier has no record).
    declared_ladder = [rung['tier'] for rung in tier_ladder(records)]
    if declared_ladder != catalog['tiering']['tiers_entered']:
        raise ValueError(
            'the published dataset tiers do not match tiering.tiers_entered '
            f'({declared_ladder} != {catalog["tiering"]["tiers_entered"]})')
    # ADR 0071 decision 2 (#206): a clinical category's tier-2 evidence must be
    # the clinical, non-genetic form; mutation-modelled iPSC-CM belongs to
    # genetic_cardiomyopathy.
    if catalog['category'] in CLINICAL_CATEGORIES:
        for record in records:
            validate_clinical_tier_basis(catalog['category'], record)
    if 'tier_sources' in catalog:
        _validate_tier_sources(catalog['tier_sources'], catalog['tiering'], records)
    _assert_relative(catalog, 'category catalog')
    return catalog


def validate_fusion_catalog(catalog):
    """Validate the emitted per-disease level-1 fusion catalog.

    Schema ``stage2-disease-fusion-tier-first-v2``: an ``axis`` record, the declared
    ``limit``, rank semantics, and one record per disease with the declared /
    informative source counts, the DEG-union and ``DEG1000`` sizes, the
    per-source ledger and the explicit exclusions. A declared-but-not-yet-
    published category is representable (``published: false``, null sizes); the
    catalog must say so rather than silently omit it.
    """
    if not isinstance(catalog, dict):
        raise ValueError('a fusion catalog must be a record')
    if catalog.get('schema') != FUSION_SCHEMA:
        raise ValueError(f'fusion catalog schema must be {FUSION_SCHEMA}')
    axis = catalog.get('axis')
    if not isinstance(axis, dict):
        raise ValueError('fusion catalog lacks its axis record')
    if not isinstance(axis.get('path'), str) or not axis['path']:
        raise ValueError('fusion catalog axis path must be a repo-relative string')
    if not isinstance(axis.get('sha256'), str) or len(axis['sha256']) != 64:
        raise ValueError('fusion catalog axis sha256 must be 64 hex characters')
    axis_size = _positive_int(axis.get('size'), 'fusion catalog axis size',
                              allow_zero=False)
    if catalog.get('limit') != DEG1000_LIMIT:
        raise ValueError(f'fusion catalog limit must be {DEG1000_LIMIT}')
    if catalog.get('eligibility_min_deg') != DISEASE_MIN_DEGS:
        raise ValueError('disease eligibility threshold must be 100')
    diseases = catalog.get('diseases')
    if not isinstance(diseases, list) or not diseases:
        raise ValueError('fusion catalog diseases must be a non-empty list')
    seen = set()
    for record in diseases:
        if not isinstance(record, dict):
            raise ValueError('a disease record must be a record')
        required = ('category', 'directory', 'published', 'tier', 'evidence_type',
                    'declared_source_count', 'informative_source_count', 'n_deg',
                    'deg1000', 'deg1000_file', 'union_file', 'sources', 'excluded_sources',
                    'not_published')
        missing = [key for key in required if key not in record]
        if missing:
            raise ValueError(f'a disease record lacks {missing}')
        directory = record['directory']
        if not isinstance(directory, str) or not directory:
            raise ValueError('a disease record has no directory')
        if directory in seen:
            raise ValueError(f'duplicate disease directory {directory}')
        seen.add(directory)
        validate_disease_eligibility(
            record, required=True)
        published = record['published']
        if not isinstance(published, bool):
            raise ValueError(f'{directory}: published must be boolean')
        declared = _positive_int(record['declared_source_count'],
                                 f'{directory}: declared_source_count')
        informative_count = _positive_int(record['informative_source_count'],
                                          f'{directory}: informative_source_count')
        sources = record['sources']
        if not isinstance(sources, list):
            raise ValueError(f'{directory}: sources must be a list')
        ledger_ids = set()
        hits = 0
        counted_informative = 0
        for source in sources:
            if not isinstance(source, dict):
                raise ValueError(f'{directory}: a source ledger entry must be a record')
            missing = [key for key in SOURCE_LEDGER_KEYS if key not in source]
            if missing:
                raise ValueError(f'{directory}: a source ledger entry lacks {missing}')
            dataset_id = source['dataset_id']
            if not isinstance(dataset_id, str) or not dataset_id:
                raise ValueError(f'{directory}: a source ledger entry has no dataset_id')
            if dataset_id in ledger_ids:
                raise ValueError(f'{directory}: duplicate source {dataset_id}')
            ledger_ids.add(dataset_id)
            _tier(source['tier'], f'{directory}/{dataset_id}: tier')
            if not isinstance(source['informative'], bool):
                raise ValueError(f'{directory}/{dataset_id}: informative must be boolean')
            source_n = _positive_int(source['n_deg'],
                                     f'{directory}/{dataset_id}: n_deg')
            source_hits = _positive_int(source['n_deg_in_deg1000'],
                                        f'{directory}/{dataset_id}: n_deg_in_deg1000')
            if source['informative']:
                if source_n == 0:
                    raise ValueError(
                        f'{directory}/{dataset_id}: an informative source needs a DEG')
                counted_informative += 1
            elif source_hits != 0:
                raise ValueError(
                    f'{directory}/{dataset_id}: an excluded source cannot land in DEG1000')
            hits += source_hits
        if informative_count != counted_informative:
            raise ValueError(
                f'{directory}: informative_source_count must equal the informative ledger '
                'entries')
        not_published = record['not_published']
        if not isinstance(not_published, list):
            raise ValueError(f'{directory}: not_published must be a list')
        for entry in not_published:
            if (not isinstance(entry, dict) or not isinstance(entry.get('dataset_id'), str)
                    or not entry['dataset_id']):
                raise ValueError(f'{directory}: a not_published entry needs a dataset_id')
        if declared != len(sources) + len(not_published):
            raise ValueError(
                f'{directory}: declared_source_count must be the published plus '
                'not-yet-published source count')
        excluded = record['excluded_sources']
        if not isinstance(excluded, list):
            raise ValueError(f'{directory}: excluded_sources must be a list')
        excluded_ids = set()
        for entry in excluded:
            if (not isinstance(entry, dict) or not isinstance(entry.get('dataset_id'), str)
                    or not entry['dataset_id']):
                raise ValueError(f'{directory}: an excluded entry needs a dataset_id')
            if entry['dataset_id'] in excluded_ids:
                raise ValueError(f'{directory}: duplicate excluded source')
            excluded_ids.add(entry['dataset_id'])
            if not isinstance(entry.get('reason'), str) or not entry['reason']:
                raise ValueError(
                    f'{directory}/{entry["dataset_id"]}: an exclusion needs a reason')
        ledger_excluded = {source['dataset_id'] for source in sources
                           if not source['informative']}
        if excluded_ids != ledger_excluded:
            raise ValueError(
                f'{directory}: excluded_sources must list exactly the non-informative '
                'ledger sources')
        n_deg = record['n_deg']
        deg1000 = record['deg1000']
        deg1000_file = record['deg1000_file']
        if not published:
            if n_deg not in (None, 0) or deg1000 not in (None, 0):
                raise ValueError(f'{directory}: an unpublished disease has no sizes')
            if record['union_file'] is not None:
                raise ValueError(f'{directory}: an unpublished disease has no union file')
            if deg1000_file is not None:
                raise ValueError(f'{directory}: an unpublished disease has no DEG1000 file')
            if informative_count or hits:
                raise ValueError(f'{directory}: an unpublished disease has no fused source')
            if record['tier'] is not None:
                _tier(record['tier'], f'{directory}: tier')
            continue
        n_deg = _positive_int(n_deg, f'{directory}: n_deg')
        deg1000 = _positive_int(deg1000, f'{directory}: deg1000')
        if deg1000 != min(n_deg, DEG1000_LIMIT):
            raise ValueError(
                f'{directory}: deg1000 must be min(n_deg, {DEG1000_LIMIT}) '
                'because the DEG union occupies the first ranks')
        if not isinstance(record['union_file'], str) or not record['union_file']:
            raise ValueError(f'{directory}: a published disease has a full union file')
        if not isinstance(deg1000_file, str) or not deg1000_file:
            raise ValueError(f'{directory}: a published disease has a DEG1000 file')
        if not isinstance(record['category'], str) or not record['category']:
            raise ValueError(f'{directory}: a published disease has a category')
        _tier(record['tier'], f'{directory}: tier')
        if not isinstance(record['evidence_type'], str) or not record['evidence_type']:
            raise ValueError(f'{directory}: a published disease has an evidence_type')
        if n_deg == 0 and informative_count:
            raise ValueError(f'{directory}: a zero-union disease cannot have an '
                             'informative source')
        if n_deg > 0 and not informative_count:
            raise ValueError(f'{directory}: a non-empty union needs an informative source')
        if hits < deg1000:
            raise ValueError(f'{directory}: every DEG1000 gene must be a DEG in at '
                             'least one source')
    _assert_relative(catalog, 'fusion catalog')
    if catalog.get('generated_by') is not None and not isinstance(
            catalog['generated_by'], str):
        raise ValueError('generated_by must be a string when present')
    return catalog


def genotype_split(genotypes):
    """Split Reichart donors by ``Primary.Genetic.Diagnosis`` (ADR 0066).

    ``genotypes`` is a frame with ``donor_id`` and ``genotype``. The split is
    genotype-first, never the study's disease label: a donor whose genotype
    declares a pathogenic variant is *not* part of ``nHF`` regardless of the
    clinical label, the single NCCM donor is not used at all, and normal donors
    supply the ``nHF`` controls. Returns a frame with ``donor_id``,
    ``genotype``, ``genotype_status`` and ``nHF_role`` (``case``, ``control`` or
    ``excluded``/``export_only``), plus a counts summary.
    """
    frame = pd.DataFrame(genotypes)
    if not {'donor_id', 'genotype'} <= set(frame.columns):
        raise ValueError('genotype table requires donor_id and genotype')
    if frame.donor_id.duplicated().any():
        raise ValueError('genotype table repeats a donor')
    roles, statuses = [], []
    for genotype in frame.genotype.astype(str):
        key = genotype.strip().lower()
        if key in ('normal', 'donor', 'healthy', 'non-dcm', 'control'):
            statuses.append('negative')
            roles.append('control')
        elif 'arvc' in key:
            statuses.append('pathogenic')
            roles.append('export_only')
        elif 'nccm' in key:
            statuses.append('pathogenic')
            roles.append('excluded')
        elif key.startswith('dcm') or 'dcm' in key:
            if any(token in key for token in ('pathogenic', 'positive', 'variant',
                                              'mutation', 'carrier')):
                statuses.append('pathogenic')
                roles.append('export_only')
            else:
                statuses.append('negative')
                roles.append('case')
        else:
            raise ValueError(f'unrecognised genotype label: {genotype!r}')
    result = frame.assign(genotype_status=statuses, nHF_role=roles)
    summary = result.groupby(['genotype_status', 'nHF_role'], sort=True).size()
    summary = {f'{status}:{role}': int(count)
               for (status, role), count in summary.items()}
    return result, summary


def contrast_selection(rows, dataset_id, contrasts):
    """Select one dataset's declared contrasts out of one joint-design DEG table.

    A shared-control source is fitted once (``~ arm`` / ``~ drug``); each
    published dataset is a declared selection of contrast labels from that single
    fit. Selecting an undeclared label or producing an empty selection is an
    error, so the joint design is recorded, not re-fitted per dataset.
    """
    rows = pd.DataFrame(rows)
    if 'contrast' not in rows:
        raise ValueError('a joint-design DEG table requires a contrast column')
    wanted = list(contrasts)
    absent = sorted(set(wanted) - set(rows.contrast.astype(str)))
    if absent:
        raise ValueError(f'{dataset_id}: contrasts not produced by the joint design: {absent}')
    selected = rows[rows.contrast.astype(str).isin(wanted)].copy()
    if selected.empty:
        raise ValueError(f'{dataset_id}: contrast selection is empty')
    return selected


# The single live DEG worker (ADR 0066 D9). The `arc` conda env's R wrapper is
# mis-patched on this workstation, so the R entry point is invoked as
# R_HOME=<home> <home>/bin/exec/R. Which R that is, is configuration, not
# science: `cpst_util.derivation.r_env` resolves it (explicit argument, then `ARC_R_HOME`, then
# the running conda env's own `lib/R`, then the `Rscript` on PATH and its
# `R.home()`), and item scripts record only the env name/version, never a host
# path.
DEG_WORKER = Path(__file__).with_name('disease_deg.R')
DEG_SUMMARY_PREFIX = 'DEG_SUMMARY '


def r_entry_point(r_home=None):
    """Return the ``<R_HOME>/bin/exec/R`` entry point and its home.

    ``r_home`` / ``ARC_R_HOME`` select the home explicitly; unset means the
    running conda env's own ``lib/R``, then the ``Rscript`` on ``PATH`` (see
    :mod:`cpst_util.derivation.r_env`). The `arc` env's own ``lib/R`` is the documented setting
    for this workstation.
    """
    return r_env.entry_point(r_home)


def run_deg_worker(matrix, samples, design, out_csv, r_home=None, timeout=7200):
    """Fit one joint design through ``disease_deg.R`` and return its summary.

    ``matrix``/``samples``/``design``/``out_csv`` are paths; the design is the
    JSON-ready dict documented in ``disease_deg.R``. The R worker writes the
    long per-contrast table and prints its fit description as the last
    ``DEG_SUMMARY <json>`` line, which is returned (never parsed for science).
    """
    import json
    import subprocess
    import tempfile

    executable, home = r_entry_point(r_home)
    with tempfile.TemporaryDirectory() as scratch:
        design_path = Path(scratch) / 'design.json'
        design_path.write_text(json.dumps(design, indent=2, sort_keys=True) + '\n')
        command = [str(executable), '--no-save', '--slave', '-f', str(DEG_WORKER),
                   '--args', '--matrix', str(matrix), '--samples', str(samples),
                   '--design', str(design_path), '--out', str(out_csv)]
        env = dict(os.environ)
        env['R_HOME'] = str(home)
        env['OMP_NUM_THREADS'] = env.get('OMP_NUM_THREADS', '2')
        env['MKL_NUM_THREADS'] = env.get('MKL_NUM_THREADS', '2')
        env['OPENBLAS_NUM_THREADS'] = env.get('OPENBLAS_NUM_THREADS', '2')
        completed = subprocess.run(command, capture_output=True, text=True,
                                   env=env, timeout=timeout)
    if completed.returncode != 0:
        raise RuntimeError('disease_deg.R failed:\n' + completed.stderr[-4000:])
    for line in reversed(completed.stdout.splitlines()):
        if line.startswith(DEG_SUMMARY_PREFIX):
            return json.loads(line[len(DEG_SUMMARY_PREFIX):])
    raise RuntimeError('disease_deg.R produced no DEG_SUMMARY line')


def category_root(output_root, directory):
    """The accepted product directory of one disease category (``nhf``, ``icm``, ...)."""
    return Path(output_root) / directory


def catalog_path(output_root, category):
    return category_root(output_root, category) / 'catalog.json'


def manifest_path(output_root, category):
    return category_root(output_root, category) / 'manifest.json'


def catalog_digest_bytes(catalog):
    """Canonical catalog bytes without the self-referential ``catalog_sha256``.

    The rule (shared with the #189 producer) is: serialise the catalog without
    the field using ``json.dumps(..., indent=2, sort_keys=True) + "\\n"`` and
    sha256 those bytes; that digest is stored in ``catalog_sha256`` and repeated
    in the category manifest. Verification pops the field, re-serialises the same
    way and compares, so the digest is a function of the catalog content only.
    """
    body = dict(catalog)
    body.pop('catalog_sha256', None)
    return (json.dumps(body, indent=2, sort_keys=True) + '\n').encode()


def catalog_digest(catalog):
    return hashlib.sha256(catalog_digest_bytes(catalog)).hexdigest()


def write_catalog(path, catalog):
    """Write a catalog with its canonical ``catalog_sha256``; return the digest."""
    digest = catalog_digest(catalog)
    body = dict(catalog)
    body['catalog_sha256'] = digest
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(body, indent=2, sort_keys=True) + '\n')
    return digest


def read_catalog(path):
    """Read a catalog and verify its embedded ``catalog_sha256``."""
    path = Path(path)
    catalog = json.loads(path.read_text())
    declared = catalog.get('catalog_sha256')
    if not isinstance(declared, str) or len(declared) != 64:
        raise ValueError(f'{path}: catalog has no canonical catalog_sha256')
    if declared != catalog_digest(catalog):
        raise ValueError(f'{path}: catalog_sha256 does not match the catalog content')
    return catalog


def build_catalog(category, directory, tier, evidence_type, axis, axis_sha256,
                  source_records, rejected=(), axis_size=AXIS_SIZE,
                  isogenic_exception=None, n_deg=None, deg1000=None, tiering=None,
                  source_deg=None, tier_sources=None):
    """Assemble a category catalog from the declared source/eligibility records.

    Canonical #190 shape: ``datasets`` is a **list of full record objects** (no
    separate id list / record map), and the content digest is carried as
    ``catalog_sha256``.

    The ADR 0071 v2 shape (#205) is opt-in: pass ``n_deg`` (the disease DEG-union
    size), ``deg1000`` (the ``DEG1000_d`` size), ``tiering`` (the tier record) and
    ``source_deg`` (``{dataset_id: {'n_deg': .., 'n_deg_in_deg1000': ..}}``, one
    entry per dataset) together to publish the fusion-derived fields and validate
    the catalog through :func:`validate_catalog`. Each source record must already
    carry its boolean ``informative`` flag (ADR 0068): the union bounds are taken
    over the informative sources only, because an excluded source's DEGs are not
    part of the fused union. ``tier_sources`` is the #206 declared tier ladder the
    ``tiering`` record was computed from (:func:`tier_ledger`); it is published
    alongside the record so the declaration can be audited (when omitted it is
    derived from the dataset records). Omitting the four fusion fields keeps the
    published v1 catalog readable while #208 republishes.
    """
    seen = set()
    records = []
    for record in source_records:
        dataset_id = record['dataset_id']
        if dataset_id in seen:
            raise ValueError(f'{category}: duplicate dataset id {dataset_id}')
        seen.add(dataset_id)
        records.append(record)
    catalog = {
        'schema': CATALOG_SCHEMA,
        'category': category,
        'directory': directory,
        'tier': int(tier),
        'evidence_type': evidence_type,
        'axis_size': int(axis_size),
        'axis': {'path': 'cpst/util/output/gene-identity/genes.csv',
                 'sha256': axis_sha256, 'size': int(axis_size)},
        'rank_contract': RANK_SCHEMA,
        'datasets': records,
        'rejected': list(rejected),
    }
    if isogenic_exception is not None:
        # ADR 0066 Isogenic exception, recorded as data (ticket #188). Absent for
        # every category that has no in-vitro isogenic model.
        catalog['isogenic_exception'] = isogenic_exception
    fusion_fields = (n_deg is not None, deg1000 is not None, tiering is not None,
                     source_deg is not None)
    if any(fusion_fields):
        if not all(fusion_fields):
            raise ValueError('the v2 catalog shape needs n_deg, deg1000, tiering and '
                             'source_deg together')
        ledger = dict(source_deg)
        unknown = sorted(set(ledger) - seen)
        absent = sorted(seen - set(ledger))
        if unknown:
            raise ValueError(f'{category}: source_deg names undeclared datasets {unknown}')
        if absent:
            raise ValueError(f'{category}: source_deg lacks datasets {absent}')
        for record in records:
            entry = ledger[record['dataset_id']]
            record['n_deg'] = int(entry['n_deg'])
            record['n_deg_in_deg1000'] = int(entry['n_deg_in_deg1000'])
        catalog['n_deg'] = int(n_deg)
        catalog.update(disease_eligibility(int(n_deg)))
        catalog['deg1000'] = int(deg1000)
        catalog['tiering'] = tiering
        # The declared ladder is published with the record. Deriving it from the
        # dataset records when the caller omits it makes the two impossible to
        # disagree (and is what a hand-written catalog cannot fake).
        ladder = tier_ladder(records) if tier_sources is None else tier_sources
        catalog['tier_sources'] = [dict(rung) for rung in ladder]
        validate_catalog(catalog, axis_size=axis_size)
    return catalog


def assemble_manifest(category, recipe, inputs, outputs, counts, catalog_sha256):
    """Assemble the per-category release manifest (ADR 0062 shape)."""
    manifest = {
        'schema': RELEASE_SCHEMA,
        'category': category,
        'recipe': recipe,
        'axis': {'path': 'cpst/util/output/gene-identity/genes.csv'},
        'inputs': dict(sorted(inputs.items())),
        'outputs': dict(sorted(outputs.items())),
        'counts': dict(sorted(counts.items())),
        'catalog_sha256': catalog_sha256,
    }
    return manifest


def _sha256(path):
    import hashlib
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


class DiseaseDeg:
    """Reader for the published per-category disease-DEG products.

    Every read is digest-checked against the category ``manifest.json``, so a
    later ticket can trust that the table it fuses is the accepted one. The
    reader is importable without the external item trees: it only reads
    ``output/disease-deg/``.
    """

    def __init__(self, root):
        self.root = Path(root)

    # ---- structure -----------------------------------------------------
    def categories(self):
        """Product directory names (``icm``, ``nhf``, ...), sorted."""
        found = []
        for path in sorted(self.root.iterdir()):
            if path.is_dir() and (path / 'catalog.json').is_file():
                found.append(path.name)
        return found

    def directory(self, key):
        """Resolve a logical category (``nHF``) or a product directory (``nhf``)."""
        if (self.root / key / 'catalog.json').is_file():
            return key
        for name in self.categories():
            catalog = json.loads((self.root / name / 'catalog.json').read_text())
            if catalog.get('category') == key:
                return name
        raise KeyError(f'unknown disease category {key!r}')

    def catalog(self, category):
        path = catalog_path(self.root, self.directory(category))
        if not path.is_file():
            raise FileNotFoundError(f'no disease-DEG catalog for category {category!r}')
        catalog = read_catalog(path)
        if catalog.get('schema') != CATALOG_SCHEMA:
            raise ValueError(f'{path}: unexpected catalog schema')
        if catalog.get('category') != category and catalog.get('directory') != category:
            raise ValueError(f'{path}: catalog category mismatch')
        return catalog

    def disease_eligibility(self, category):
        """Disease-level admission, distinct from a source's eligibility fields."""
        catalog = self.catalog(category)
        if 'n_deg' not in catalog:
            raise ValueError('disease eligibility needs the fused axis-resolved DEG union')
        validate_disease_eligibility(catalog)
        return disease_eligibility(catalog['n_deg'])

    def manifest(self, category):
        path = manifest_path(self.root, self.directory(category))
        if not path.is_file():
            raise FileNotFoundError(f'no disease-DEG manifest for category {category!r}')
        manifest = json.loads(path.read_text())
        if manifest.get('schema') != RELEASE_SCHEMA:
            raise ValueError(f'{path}: unexpected release schema')
        return manifest

    def datasets(self, category):
        """Dataset ids of a category, in catalog order."""
        return [record['dataset_id'] for record in self.catalog(category)['datasets']]

    def record(self, category, dataset_id):
        catalog = self.catalog(category)
        for record in catalog['datasets']:
            if record['dataset_id'] == dataset_id:
                return record
        raise KeyError(f'{category}: unknown dataset {dataset_id!r}')

    def eligibility(self, category, dataset_id):
        """The declared eligibility/shared-control fields of one dataset."""
        record = self.record(category, dataset_id)
        keys = ('category', 'source_item', 'species', 'model_system', 'evidence_type',
                'tier', 'genotype_status', 'time_window_basis', 'exploratory',
                'shared_control_group', 'n_case', 'n_control', 'method')
        return {key: record[key] for key in keys if key in record}

    # ---- payloads ------------------------------------------------------
    def _verified_path(self, category, name):
        manifest = self.manifest(category)
        declared = manifest['outputs'].get(name)
        if not isinstance(declared, dict):
            raise ValueError(f'{category}: {name} is not a declared output')
        path = category_root(self.root, self.directory(category)) / name
        if not path.is_file():
            raise FileNotFoundError(f'{category}: missing declared output {name}')
        if _sha256(path) != declared['sha256']:
            raise ValueError(f'{category}: {name} digest differs from the manifest')
        return path

    def table(self, category, dataset_id):
        record = self.record(category, dataset_id)
        table = pd.read_csv(self._verified_path(category, record['deg']),
                            float_precision='round_trip')
        return validate_table(table)

    def rank(self, category, dataset_id):
        record = self.record(category, dataset_id)
        rank = pd.read_csv(self._verified_path(category, record['rank']),
                           dtype={'gene_index': np.int64},
                           float_precision='round_trip')
        return validate_rank(rank)

    def digest(self, category, name):
        """Declared sha256 of one accepted output."""
        return self.manifest(category)['outputs'][name]['sha256']


class DiseaseFusion:
    """Reader for the emitted per-disease level-1 fusion + ``DEG1000`` product.

    ``output/disease-fusion/`` holds one catalog and one ``<directory>.deg1000.csv``
    per published disease (ADR 0071 decision 4). The catalog is the
    self-describing record of the fusion (declared/informative source counts, the
    DEG-union and ``DEG1000`` sizes, the per-source ledger and every exclusion);
    the CSV is the per-gene ``DEG1000_d`` table in disease-rank order. Both are
    digest-checked against ``manifest.json``. The reader only reads
    ``output/disease-fusion/``, so it imports without the external item trees.
    """

    def __init__(self, root):
        self.root = Path(root)

    def manifest(self):
        path = self.root / 'manifest.json'
        if not path.is_file():
            raise FileNotFoundError(f'no disease-fusion manifest at {path}')
        manifest = json.loads(path.read_text())
        if manifest.get('schema') != FUSION_RELEASE_SCHEMA:
            raise ValueError(f'{path}: unexpected fusion release schema')
        return manifest

    def catalog(self):
        path = self.root / 'catalog.json'
        if not path.is_file():
            raise FileNotFoundError(f'no disease-fusion catalog at {path}')
        catalog = read_catalog(path)
        if catalog.get('schema') != FUSION_SCHEMA:
            raise ValueError(f'{path}: unexpected fusion catalog schema')
        manifest = self.manifest()
        declared = manifest['outputs'].get('catalog.json', {})
        if (declared.get('sha256') != _sha256(path)
                or manifest.get('catalog_sha256') != catalog['catalog_sha256']
                or manifest.get('limit') != catalog.get('limit')):
            raise ValueError('fusion catalog differs from its release manifest')
        return validate_fusion_catalog(catalog)

    def diseases(self):
        """Product directories of the published diseases, in catalog order."""
        return [record['directory'] for record in self.catalog()['diseases']
                if record['published']]

    def eligible_diseases(self):
        """Published fusion participants; ``diseases()`` retains excluded evidence.

        Admission is validated against the complete recorded disease union.
        """
        return [record['directory'] for record in self.catalog()['diseases']
                if disease_eligibility(record['n_deg'], record['published'])['eligible']]

    def record(self, directory):
        for record in self.catalog()['diseases']:
            if record['directory'] == directory:
                return record
        raise KeyError(f'unknown disease-fusion directory {directory!r}')

    def _frame(self, record, field):
        name = record[field]
        manifest = self.manifest()
        declared = manifest['outputs'].get(name)
        if not isinstance(declared, dict):
            raise ValueError(f'{name} is not a declared fusion output')
        path = self.root / name
        if not path.is_file():
            raise FileNotFoundError(f'missing declared fusion output {name}')
        if _sha256(path) != declared['sha256']:
            raise ValueError(f'{name}: digest differs from the fusion manifest')
        frame = pd.read_csv(path, float_precision='round_trip')
        if list(frame.columns) != list(UNION_COLUMNS):
            raise ValueError(f'{name}: unexpected fusion evidence columns')
        return frame

    def union(self, directory):
        """Digest-checked complete union and actual DEG-source tier evidence."""
        record = self.record(directory)
        if not record['published']:
            raise KeyError(f'{directory}: the disease is declared but not published')
        frame = self._frame(record, 'union_file')
        if len(frame) != record['n_deg']:
            raise ValueError('full union size differs from catalog')
        if frame.disease_rank.tolist() != list(range(1, record['n_deg'] + 1)):
            raise ValueError('full union ranks must be contiguous')
        if (not pd.api.types.is_integer_dtype(frame.gene_index) and len(frame)):
            raise ValueError('full union gene_index must be integral')
        if not frame.gene_index.is_unique or not frame.gene_index.between(
                0, self.catalog()['axis']['size'] - 1).all():
            raise ValueError('invalid full union gene_index')
        ledger = {source['dataset_id']: source for source in record['sources']}
        counts = {source: 0 for source in ledger}
        for row in frame.itertuples(index=False):
            support = json.loads(row.supporting_sources)
            if (not isinstance(support, list) or not support
                    or any(not isinstance(source, str) for source in support)
                    or support != sorted(set(support))):
                raise ValueError('union supporting sources must be unique sorted source ids')
            if any(source not in ledger or not ledger[source]['informative']
                   for source in support):
                raise ValueError('union support must name informative declared sources')
            if row.source_tier_highest != min(ledger[source]['tier'] for source in support):
                raise ValueError('highest source tier disagrees with DEG support')
            for source in support:
                counts[source] += 1
        for source, count in counts.items():
            expected = ledger[source]['n_deg'] if ledger[source]['informative'] else 0
            if count != expected:
                raise ValueError('full union source membership differs from source DEG count')
        if not np.isfinite(frame.disease_rank_sum.to_numpy(dtype=float)).all():
            raise ValueError('disease rank sums must be finite')
        ordered = frame.sort_values(['source_tier_highest', 'disease_rank_sum', 'gene_index'])
        if ordered.gene_index.tolist() != frame.gene_index.tolist():
            raise ValueError('union is not ordered by tier, rank sum and gene index')
        return frame

    def deg1000(self, directory):
        """Tier-first cutoff, validated against the retained complete union."""
        record = self.record(directory)
        if not record['published']:
            raise KeyError(f'{directory}: the disease is declared but not published')
        frame = self._frame(record, 'deg1000_file')
        expected = self.union(directory).iloc[:record['deg1000']]
        if not frame.equals(expected.reset_index(drop=True)):
            raise ValueError('DEG1000 differs from the first full-union ranks')
        support = [json.loads(value) for value in frame.supporting_sources]
        for source in record['sources']:
            if sum(source['dataset_id'] in genes for genes in support) != source['n_deg_in_deg1000']:
                raise ValueError('source DEG1000 membership differs from source ledger')
        return frame
