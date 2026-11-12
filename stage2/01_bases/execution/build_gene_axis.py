#!/usr/bin/env python3
"""Build or replay the Stage 2 gene-axis HVG/rank-sum decision evidence.

The product is decision preparation (spec section C, ticket #183): per-source
HVG and DEG ranks over the frozen 18,164-gene axis, an unweighted rank-sum
fusion plus the animal half-weight sensitivity, source overlap, direction
composition, an explicit inflection answer, and the tenv non-HVG curve. It does
not choose the axis size, a source quota, a direction rule or the animal weight
(those are #173), builds no corpus, trains no model, and never opens a truth or
Heliste label file.

Every scientific definition lives in ``output/api/bases/gene_axis.py`` (the
single public implementation, ADR 0062); the ``pert_data/`` cardiotoxicity
adapters live in ``execution/gene_axis_sources.py``. This driver only reads the
inputs, applies those definitions and writes the evidence.

    python stage2/01_bases/execution/build_gene_axis.py
    python stage2/01_bases/execution/build_gene_axis.py --verify

``--verify`` rebuilds the complete evidence in a temporary directory from the
pinned inputs and checks every tracked digest; it reads no truth product and
writes nothing under ``output/``.
"""
import argparse
import hashlib
import json
import sys
from dataclasses import dataclass, replace
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases.gene_axis import (
    ANIMAL_SOURCES, AXIS_SIZE, PADJ_THRESHOLD, SCHEMA, SOURCE_HVG_K, TENV_HVG_KS,
    SymbolResolver, collapse_expression, complete_rank, deg_effect, detect_inflection,
    direction_composition, fuse_ranks, fused_order, hvg_membership, normalized_dispersion,
    overlap_metrics, per_gene_deg, rank_scores, tenv_non_hvg_curve,
)
import gene_axis_sources as sources

STEP = Path(__file__).resolve().parents[1]
DIRECTION_SIZES = (500, 1000, 2000, 3000, 5000)
# Inflection is evaluated on the sorted primary rank-sum curve.
INFLECTION_X = 'fused_rank'
INFLECTION_Y = 'rank_sum'
HISTOGRAM_BINS = 40
CURVE_SAMPLE_STEP = 100
FORBIDDEN = ('tenv-truth', 'heliste')

RECIPE = {
    'hvg': ('Seurat-style normalized dispersion on the item\'s deposited log expression: '
            'per-gene mean over finite samples, sample variance / mean, log dispersion '
            'standardized in 20 equal-width log1p-mean bins, singleton bins score 1, '
            'undefined dispersion sorts last, gene_index tie-break. Same family as '
            'output/api/bases/hvg.py; deliberate difference: no expm1 recovery because '
            'the disease matrices are heterogeneous non-invertible log scales, so the '
            'statistics are computed on the deposited log values directly.'),
    'deg': ('Continuous per-contrast effect score |log2FC| * -log10(padj) with padj floored '
            'at 1e-300; NA padj (DESeq2 independent filtering) is treated as '
            'non-significant with score 0, not dropped from the measured set. A gene '
            'scores the maximum over its source contrasts (declared collapse rule); '
            'direction and padj come from that winning contrast and '
            'n_contrasts_supporting counts padj < 0.05 contrasts.'),
    'alignment': ('Source symbol -> axis gene_index by precedence: canonical axis symbol; '
                  'gene_crosswalk.csv alias/previous rows; pinned HGNC alias_symbol / '
                  'prev_symbol; Ensembl id (version stripped) via the axis ensembl_id or '
                  'the HGNC ensembl_gene_id. Multi-candidate aliases are ambiguous and '
                  'unresolved, never guessed; unresolved symbols are counted and reported.'),
    'expression_collapse': ('Duplicate source symbols resolving to one axis gene are '
                            'averaged on the deposited log expression.'),
    'rank': ('Within-source only. Measured genes are ranked 1..n_measured by descending '
             'score with gene_index as tie-break; unmeasured genes get n_measured + 1 and '
             'rank_norm = 1.0. Ranks are never compared across sources as raw effect sizes.'),
    'fusion': ('Unweighted sum of rank_norm = rank / n_measured (unmeasured = 1.0) over '
               'the fused sources; primary fused rank is ascending rank_sum with gene_index '
               'tie-break. rank_sum_animal_half re-weights the two animal MIRI sources at '
               '0.5 (reported sensitivity, not chosen); the primary fusion is unweighted.'),
    'direction': ('Per gene the significant DEG contrast with the largest effect score '
                  'across all DEG sources gives the direction; exact ties break by the '
                  'declared source order and are counted. The symmetric rule counts genes '
                  'equally, the abs_log2fc rule weights each gene by its winning |log2FC|.'),
    'overlap': ('TenV HVG sets are the published descriptive top-K (K in {1000, 2000, 5000}); '
                'a source HVG set is its top-2000 measured genes; a source DEG set is '
                'padj < 0.05 on the winning contrast.'),
    'inflection': ('Kneedle online S=1.0 (first interior local maximum of the normalized '
                   'arc y_norm - x_norm) plus the second-difference peak; accepted only when '
                   'both candidates agree within 25% and the knee exceeds 0.01 normalized '
                   'units, otherwise an explicit "inflection": false.'),
}

SOURCE_SPECS = [
    {'id': 'tenv', 'column': 'tenv_hvg_rank', 'accession': 'tenv', 'category': 'tenv',
     'species': 'Homo sapiens', 'tissue_model': 'hiPSC-CM TempO-Seq (tenv corpus)',
     'contrast': 'descriptive variability (published tenv HVG rank)',
     'method': 'published descriptive normalized dispersion (output/tenv-hvg, not recomputed)',
     'normalization': 'log1p CP10K', 'ortholog_chain': 'none (human)',
     'fused': True, 'exploratory': False, 'supplementary_only': False,
     'notes': 'variability reference only; no significance or truth input'},
    {'id': 'gse141910', 'column': 'gse141910_rank', 'accession': 'GSE141910',
     'category': 'disease_data', 'species': 'Homo sapiens',
     'tissue_model': 'human left-ventricular free-wall bulk tissue (MAGNet)',
     'contrast': 'HF_vs_NF', 'method': 'limma on deposited voom log2 abundance',
     'ortholog_chain': 'none (human)', 'fused': True, 'exploratory': False,
     'supplementary_only': False, 'notes': ''},
    {'id': 'gse116250', 'column': 'gse116250_rank', 'accession': 'GSE116250',
     'category': 'disease_data', 'species': 'Homo sapiens',
     'tissue_model': 'human left-ventricular bulk tissue (HF / post-MI proxy)',
     'contrast': 'HF_vs_NF (+ DCM_vs_NF, ICM_vs_NF contrasts, max-collapsed)',
     'method': 'limma on log2(RPKM+1)', 'ortholog_chain': 'none (human)',
     'fused': True, 'exploratory': False, 'supplementary_only': False, 'notes': ''},
    {'id': 'gse243674', 'column': 'gse243674_rank', 'accession': 'GSE243674',
     'category': 'disease_data', 'species': 'Homo sapiens',
     'tissue_model': 'hiPSC-CM monolayer, 6 donors, anthracyclines + vehicle',
     'contrast': '10 drug-vs-vehicle contrasts (max-collapsed)',
     'method': 'DESeq2 Wald on hg38 gene counts', 'ortholog_chain': 'none (human)',
     'fused': True, 'exploratory': False, 'supplementary_only': False, 'notes': ''},
    {'id': 'gse128188', 'column': 'gse128188_rank', 'accession': 'GSE128188',
     'category': 'disease_data', 'species': 'Homo sapiens',
     'tissue_model': 'human paired left/right atrial appendage',
     'contrast': 'AF_vs_SR', 'method': 'DESeq2 Wald on donor-pooled counts',
     'ortholog_chain': 'none (human)', 'fused': True, 'exploratory': False,
     'supplementary_only': False, 'notes': ''},
    {'id': 'gse282618', 'column': 'gse282618_rank', 'accession': 'GSE282618',
     'category': 'disease_data', 'species': 'Homo sapiens',
     'tissue_model': 'human post-mortem left/right ventricle (hypertensive heart disease)',
     'contrast': 'HHD_vs_control', 'method': 'DESeq2 Wald on subject-pooled counts',
     'ortholog_chain': 'none (human)', 'fused': True, 'exploratory': True,
     'supplementary_only': False, 'notes': 'exploratory: 3 cases / 2 controls'},
    {'id': 'gse160516', 'column': 'gse160516_rank', 'accession': 'GSE160516',
     'category': 'disease_data', 'species': 'Mus musculus',
     'tissue_model': 'mouse MIRI (pooled IR vs Sham, Affymetrix Clariom S)',
     'contrast': 'IR_vs_Sham (pooled 6h/24h/72h)',
     'method': 'limma on the item\'s full-transcriptome mouse DEG table',
     'ortholog_chain': 'mouse symbol -> homologene mouse2human (inTax 10090) -> human symbol -> frozen axis',
     'fused': True, 'exploratory': False, 'supplementary_only': False,
     'notes': ('animal DEG-only source; half-weight sensitivity column; the local DEG '
               'interface carries no sample count')},
    {'id': 'gse108940', 'column': 'gse108940_rank', 'accession': 'GSE108940',
     'category': 'disease_data', 'species': 'Mus musculus',
     'tissue_model': 'mouse MIRI (WT IR vs WT Sham, GPL7202)',
     'contrast': 'WT_IR_vs_WT_Sham',
     'method': 'limma on the corrected full-transcriptome mouse DEG table',
     'ortholog_chain': 'mouse symbol -> homologene mouse2human (inTax 10090) -> human symbol -> frozen axis',
     'fused': True, 'exploratory': False, 'supplementary_only': False,
     'notes': ('animal DEG-only source; half-weight sensitivity column; n_samples from '
               'the item build script design (6 IR + 6 sham)')},
    {'id': 'dtoxs_ipsc', 'column': 'dtoxs_ipsc_rank', 'accession': 'GSE217421',
     'category': 'pert_data', 'species': 'Homo sapiens',
     'tissue_model': 'DToxS iPSC-CM (GSE217424 SuperSeries)',
     'contrast': '22 cardiotoxic-positive comparison tables (DAU/DOX/EPI/IDA/TRA, max-collapsed)',
     'method': 'upstream GeoMx/RNA-seq Calc DEG tables (logFC, PValue, FDR)',
     'ortholog_chain': 'none (human)',
     'fused': True, 'exploratory': False, 'supplementary_only': False,
     'notes': 'no MTX table in this deposit; never re-downloaded'},
    {'id': 'gse244740', 'column': 'gse244740_rank', 'accession': 'GSE244740',
     'category': 'pert_data', 'species': 'Homo sapiens',
     'tissue_model': 'ScreenSeq hiPSC-CM panel (Rosell-Hidalgo)',
     'contrast': 'DOX / IDA / Mitoxantrone vs same-plate DMSO, all doses pooled (max-collapsed)',
     'method': ('log2(CPM+1) over all source genes, two-sample Welch t-test, BH within '
                'contrast (declared method; no moderation)'),
     'ortholog_chain': 'none (human; Ensembl -> axis symbol)',
     'fused': True, 'exploratory': False, 'supplementary_only': False,
     'notes': 'DEG derived from the already-owned raw counts; 22 same-plate DMSO wells'},
    {'id': 'gse137293', 'column': '', 'accession': 'GSE137293', 'category': 'disease_data',
     'species': 'Sus scrofa', 'tissue_model': 'pig infarct scar vs border zone',
     'contrast': 'scar_vs_border (not disease vs healthy)',
     'method': 'listed only', 'normalization': '', 'ortholog_chain': 'none (supplementary)',
     'fused': False, 'exploratory': False, 'supplementary_only': True,
     'notes': 'supplementary evidence only: an injury-severity gradient, never fused'},
    {'id': 'lincs_lds_1593', 'column': '', 'accession': 'LINCS-LDS-1593',
     'category': 'pert_data', 'species': 'Homo sapiens',
     'tissue_model': 'L1000-978 landmark panel', 'contrast': 'coverage note only',
     'method': 'listed only', 'normalization': 'L1000-978', 'ortholog_chain': 'none',
     'fused': False, 'exploratory': False, 'supplementary_only': True,
     'notes': 'coverage note only (978 landmark genes); never a fusion source'},
]
EXPRESSION_IDS = ('gse141910', 'gse116250', 'gse243674', 'gse128188', 'gse282618')
ANIMAL_IDS = ('gse160516', 'gse108940')
PERT_IDS = ('dtoxs_ipsc', 'gse244740')
FUSED_IDS = ('tenv',) + EXPRESSION_IDS + ANIMAL_IDS + PERT_IDS
DEG_IDS = EXPRESSION_IDS + ANIMAL_IDS + PERT_IDS


@dataclass
class Config:
    """Every path the build reads or writes; all are repo-relative by default."""
    root: Path
    genes: Path
    crosswalk: Path
    hgnc: Path
    tenv_hvg: Path
    orthologs: Path
    disease_root: Path
    pert_root: Path
    output: Path
    figure_data: Path
    intermediates: Path

    def input_roles(self):
        return [
            ('frozen_axis', self.genes),
            ('gene_crosswalk', self.crosswalk),
            ('hgnc_complete_set', self.hgnc),
            ('tenv_hvg_descriptive', self.tenv_hvg),
            ('mouse_human_orthologs', self.orthologs),
            *[(f'{source}_expr', self.disease_root / source.upper() / 'local/expr'
               / f'{source.upper()}_log_expr.tsv.gz') for source in EXPRESSION_IDS],
            *[(f'{source}_deg', self.disease_root / source.upper() / 'local/deg'
               / f'{source.upper()}_deg.csv') for source in EXPRESSION_IDS],
            ('gse160516_deg', self.disease_root / 'GSE160516/source/GSE160516_deg_results.csv'),
            ('gse108940_deg', self.disease_root / 'GSE108940/source/GSE108940_WT_deg_results.csv'),
            ('dtoxs_deg_tar', self.pert_root / 'DToxS-iPSC/source'
             / 'GSE217421_Conv-DEG-Calc-FDR-0.1.tar.gz'),
            ('gse244740_counts', self.pert_root / 'GSE244740/source'
             / 'GSE244740_processed_data_counts.txt.gz'),
            ('gse244740_series', self.pert_root / 'GSE244740/source'
             / 'GSE244740_series_matrix.txt.gz'),
        ]


def repo_root():
    return next(parent for parent in Path(__file__).resolve().parents if (parent / '.git').exists())


def default_config(root=None, output=None, figure_data=None, intermediates=None):
    root = Path(root) if root else repo_root()
    step = root / 'stage2/01_bases'
    return Config(
        root=root,
        genes=step / 'output/gene-identity/genes.csv',
        crosswalk=step / 'output/gene-identity/tables/gene_crosswalk.csv',
        hgnc=root / 'pert_data/Annotation/source/hgnc_complete_set.txt',
        tenv_hvg=step / 'output/tenv-hvg/descriptive.csv',
        orthologs=step / 'output/gene-axis/mouse_human_orthologs.csv',
        disease_root=root / 'disease_data',
        pert_root=root / 'pert_data',
        output=Path(output) if output else step / 'output/gene-axis',
        figure_data=Path(figure_data) if figure_data else step / 'figure/data/gene_axis_source_summary.csv',
        intermediates=Path(intermediates) if intermediates is not None
        else step / 'execution/gene_axis',
    )


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def forbidden_input(path):
    """Raise when a path would open a truth or Heliste label product.

    The evidence product is label-free: the only tenv input is the published
    descriptive HVG ranking. This guard makes the boundary executable and is
    exercised by the leakage test.
    """
    parts = [part.lower() for part in Path(path).parts]
    for token in FORBIDDEN:
        if any(token in part for part in parts):
            raise ValueError(f'forbidden label input: {path}')
    return Path(path)


def read_csv(path, **kwargs):
    return pd.read_csv(forbidden_input(path), **kwargs)


def write_frame(frame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    text = frame.to_csv(index=False, lineterminator='\n', float_format='%.12g')
    path.write_text(text)
    return path


def write_json(payload, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n')
    return path


def _expression_metadata(path):
    import gzip
    with gzip.open(path, 'rt') as handle:
        line = handle.readline().rstrip('\n')
    metadata = {}
    if line.startswith('#'):
        for token in line[1:].split():
            if '=' in token:
                key, value = token.split('=', 1)
                metadata[key] = value
    return metadata


def read_expression(path):
    """Read a disease expression matrix and its declared metadata header."""
    metadata = _expression_metadata(path)
    with Path(path).open('rb') as raw:
        import gzip
        with gzip.open(raw, 'rt') as handle:
            handle.readline()
            frame = pd.read_csv(handle, sep='\t', index_col=0)
    if frame.index.has_duplicates:
        raise ValueError(f'duplicate gene symbols in {path}')
    return metadata, frame


def read_axis(path):
    axis = read_csv(path, dtype={'symbol': str, 'hgnc_id': str,
                                 'ensembl_id': str, 'entrez_id': str})
    if list(axis.columns[:2]) != ['gene_index', 'symbol'] or len(axis) != AXIS_SIZE:
        raise ValueError('unexpected frozen axis layout')
    if sorted(axis.gene_index.astype(int)) != list(range(AXIS_SIZE)):
        raise ValueError('frozen axis gene_index is not 0..18163')
    return axis


def load_resolver(config):
    axis = read_axis(config.genes)
    crosswalk = read_csv(config.crosswalk, dtype=str, keep_default_na=False)
    hgnc = read_csv(config.hgnc, sep='\t', dtype=str, keep_default_na=False,
                    usecols=['hgnc_id', 'symbol', 'alias_symbol', 'prev_symbol',
                             'ensembl_gene_id'])
    return axis, SymbolResolver(axis, crosswalk.to_dict('records'), hgnc.to_dict('records'))


def _alignment_rows(source_id, symbols, indices, categories):
    return pd.DataFrame({'source': source_id, 'source_symbol': symbols,
                         'gene_index': indices, 'resolution': categories})


def _deg_rank(measured_index, effect):
    rank = rank_scores(measured_index, effect)
    return complete_rank(measured_index, rank)


def process_expression_source(config, axis, resolver, spec, intermediates):
    """Per-source HVG rank from the expression matrix and DEG rank from its DEG table."""
    accession = spec['accession']
    expr_path = config.disease_root / accession / 'local/expr' / f'{accession}_log_expr.tsv.gz'
    deg_path = config.disease_root / accession / 'local/deg' / f'{accession}_deg.csv'
    metadata, matrix = read_expression(expr_path)
    symbols = [str(value) for value in matrix.index]
    indices, categories = resolver.resolve_many(symbols)
    alignment = _alignment_rows(spec['id'], symbols, indices, categories)
    loss = resolver.loss_summary(indices, categories)
    keep = indices >= 0
    gene_index, values = collapse_expression(matrix.to_numpy(dtype=np.float64)[keep],
                                             indices[keep], AXIS_SIZE)
    hvg = normalized_dispersion(values, gene_index)
    hvg_top = hvg.sort_values('rank').gene_index.to_numpy(dtype=np.int64)[:SOURCE_HVG_K]
    hvg_set = np.zeros(AXIS_SIZE, dtype=bool)
    hvg_set[hvg_top] = True
    deltas = read_csv(deg_path, dtype={'gene_symbol': str, 'contrast': str})
    deg_indices, deg_categories = resolver.resolve_many([str(v) for v in deltas.gene_symbol])
    deg_rows = pd.DataFrame({'gene_index': deg_indices,
                             'contrast': deltas.contrast.to_numpy(),
                             'log2fc': pd.to_numeric(deltas.log2FC, errors='coerce').to_numpy(),
                             'padj': pd.to_numeric(deltas.padj, errors='coerce').to_numpy()})
    deg_rows = deg_rows[deg_rows.gene_index >= 0]
    per_gene = per_gene_deg(deg_rows.gene_index, deg_rows.log2fc, deg_rows.padj,
                            deg_rows.contrast)
    measured = per_gene.gene_index.to_numpy(dtype=np.int64)
    rank = _deg_rank(measured, per_gene.effect.to_numpy(dtype=np.float64))
    significant = np.zeros(AXIS_SIZE, dtype=np.int64)
    significant[per_gene.loc[per_gene.padj < PADJ_THRESHOLD, 'gene_index'].to_numpy(dtype=np.int64)] = 1
    result = {'rank': rank, 'measured': len(measured), 'significant': significant,
              'alignment': alignment, 'loss': loss, 'per_gene': per_gene,
              'hvg': hvg, 'hvg_set': hvg_set, 'metadata': metadata,
              'n_samples': int(matrix.shape[1]), 'n_contrasts': int(deltas.contrast.nunique()),
              'n_axis_expression': int(len(gene_index))}
    if intermediates is not None:
        write_frame(alignment, intermediates / 'sources' / f'{spec["id"]}_alignment.csv')
        write_frame(hvg, intermediates / 'sources' / f'{spec["id"]}_hvg.csv')
        write_frame(per_gene[['gene_index', 'contrast', 'log2fc', 'padj', 'effect',
                              'significant', 'direction', 'n_contrasts_supporting']],
                    intermediates / 'sources' / f'{spec["id"]}_deg.csv')
    return result


def process_animal_source(config, axis, resolver, spec, intermediates, orthologs):
    """DEG-only mouse source: homologene mouse->human then the frozen-axis resolver."""
    if spec['id'] == 'gse160516':
        deg_path = config.disease_root / 'GSE160516/source/GSE160516_deg_results.csv'
        columns = {'Gene': 'symbol', 'logFC': 'log2fc', 'adj.P.Val': 'padj'}
    else:
        deg_path = config.disease_root / 'GSE108940/source/GSE108940_WT_deg_results.csv'
        columns = {'Gene': 'symbol', 'logFC': 'log2fc', 'adj.P.Val': 'padj'}
    table = read_csv(deg_path, usecols=list(columns), dtype={list(columns)[0]: str})
    table = table.rename(columns=columns)
    table['log2fc'] = pd.to_numeric(table['log2fc'], errors='coerce')
    table['padj'] = pd.to_numeric(table['padj'], errors='coerce')
    mouse = [str(value) for value in table.symbol]
    unique_mouse = sorted(set(mouse))
    mapped = orthologs[orthologs.mouse_symbol.isin(unique_mouse)]
    # One mouse symbol may map to several human symbols; resolve every human
    # candidate and keep the mapping only when exactly one axis gene matches.
    candidates = mapped.groupby('mouse_symbol').human_symbol.apply(list).to_dict()
    mapping = {}
    ambiguous = offaxis = unresolved = 0
    for mouse_symbol, human_symbols in candidates.items():
        hits = set()
        for human in human_symbols:
            index, _ = resolver.resolve(human)
            if index >= 0:
                hits.add(index)
        if len(hits) == 1:
            mapping[mouse_symbol] = next(iter(hits))
        elif len(hits) > 1:
            ambiguous += 1
        else:
            unresolved += 1
    gene_indices = np.asarray([mapping.get(symbol, -1) for symbol in mouse], dtype=np.int64)
    rows = pd.DataFrame({'gene_index': gene_indices, 'log2fc': table.log2fc.to_numpy(),
                         'padj': table.padj.to_numpy(), 'contrast': 'miri'})
    rows = rows[rows.gene_index >= 0]
    per_gene = per_gene_deg(rows.gene_index, rows.log2fc, rows.padj, rows.contrast)
    measured = per_gene.gene_index.to_numpy(dtype=np.int64)
    rank = _deg_rank(measured, per_gene.effect.to_numpy(dtype=np.float64))
    significant = np.zeros(AXIS_SIZE, dtype=np.int64)
    significant[per_gene.loc[per_gene.padj < PADJ_THRESHOLD, 'gene_index'].to_numpy(dtype=np.int64)] = 1
    mapping_frame = pd.DataFrame(sorted(mapping.items()), columns=['mouse_symbol', 'gene_index'])
    result = {'rank': rank, 'measured': len(measured), 'significant': significant,
              'per_gene': per_gene, 'alignment': mapping_frame,
              'loss': {'n_mouse_symbols': len(unique_mouse),
                       'n_mouse_symbols_with_ortholog': int(len(candidates)),
                       'n_mouse_symbols_mapped_to_axis': int(len(mapping)),
                       'n_ambiguous_mouse_symbols': ambiguous,
                       'n_without_axis_hit': unresolved,
                       'n_mouse_rows_dropped': int((gene_indices < 0).sum())},
              'n_samples': 12 if spec['id'] == 'gse108940' else None,
              'n_contrasts': 1, 'metadata': {}}
    if intermediates is not None:
        write_frame(mapping_frame, intermediates / 'sources' / f'{spec["id"]}_ortholog_mapping.csv')
        write_frame(per_gene[['gene_index', 'contrast', 'log2fc', 'padj', 'effect',
                              'significant', 'direction', 'n_contrasts_supporting']],
                    intermediates / 'sources' / f'{spec["id"]}_deg.csv')
    return result


def ensembl_symbol_map(axis, config):
    """Version-stripped Ensembl id -> frozen-axis symbol (axis table first, HGNC fallback)."""
    mapping = {}
    for ensembl, symbol in zip(axis.ensembl_id, axis.symbol):
        if isinstance(ensembl, str) and ensembl:
            mapping[ensembl.split('.')[0]] = symbol
    by_hgnc = dict(zip(axis.hgnc_id, axis.symbol))
    hgnc = read_csv(config.hgnc, sep='\t', dtype=str, keep_default_na=False,
                    usecols=['hgnc_id', 'ensembl_gene_id'])
    for hgnc_id, ensembl in zip(hgnc.hgnc_id, hgnc.ensembl_gene_id):
        key = str(ensembl).split('.')[0]
        if key and key not in mapping and hgnc_id in by_hgnc:
            mapping[key] = by_hgnc[hgnc_id]
    return mapping


def process_pert_sources(config, axis, resolver, intermediates, cache_dir):
    """The two cardiotoxicity sources already owned under ``pert_data/``."""
    results = {}
    dtoxs_frame, dtoxs_provenance = sources.dtoxs_contrasts(
        config.pert_root / 'DToxS-iPSC/source/GSE217421_Conv-DEG-Calc-FDR-0.1.tar.gz',
        cache_dir=cache_dir)
    indices, categories = resolver.resolve_many([str(value) for value in dtoxs_frame.gene_symbol])
    rows = pd.DataFrame({'gene_index': indices, 'contrast': dtoxs_frame.contrast.to_numpy(),
                         'log2fc': dtoxs_frame.log2fc.to_numpy(),
                         'padj': dtoxs_frame.padj.to_numpy()})
    rows = rows[rows.gene_index >= 0]
    per_gene = per_gene_deg(rows.gene_index, rows.log2fc, rows.padj, rows.contrast)
    measured = per_gene.gene_index.to_numpy(dtype=np.int64)
    rank = _deg_rank(measured, per_gene.effect.to_numpy(dtype=np.float64))
    significant = np.zeros(AXIS_SIZE, dtype=np.int64)
    significant[per_gene.loc[per_gene.padj < PADJ_THRESHOLD, 'gene_index'].to_numpy(dtype=np.int64)] = 1
    results['dtoxs_ipsc'] = {
        'rank': rank, 'measured': len(measured), 'significant': significant,
        'per_gene': per_gene, 'n_contrasts': int(dtoxs_frame.contrast.nunique()),
        'n_samples': int(dtoxs_provenance['n_distinct_wells']),
        'loss': {'n_contrast_rows': len(dtoxs_frame),
                 'n_rows_resolved': int((indices >= 0).sum()),
                 'n_symbols_unresolved': int((indices < 0).sum())},
        'provenance': dtoxs_provenance, 'metadata': {}}
    if intermediates is not None:
        write_frame(per_gene[['gene_index', 'contrast', 'log2fc', 'padj', 'effect',
                              'significant', 'direction', 'n_contrasts_supporting']],
                    intermediates / 'sources' / 'dtoxs_ipsc_deg.csv')

    counts = config.pert_root / 'GSE244740/source/GSE244740_processed_data_counts.txt.gz'
    series = config.pert_root / 'GSE244740/source/GSE244740_series_matrix.txt.gz'
    gse_frame, gse_provenance = sources.gse244740_contrasts(
        counts, series, ensembl_symbol_map(axis, config), cache_dir=cache_dir)
    indices, categories = resolver.resolve_many([str(value) for value in gse_frame.gene_symbol])
    rows = pd.DataFrame({'gene_index': indices, 'contrast': gse_frame.contrast.to_numpy(),
                         'log2fc': gse_frame.log2fc.to_numpy(),
                         'padj': gse_frame.padj.to_numpy()})
    rows = rows[rows.gene_index >= 0]
    per_gene = per_gene_deg(rows.gene_index, rows.log2fc, rows.padj, rows.contrast)
    measured = per_gene.gene_index.to_numpy(dtype=np.int64)
    rank = _deg_rank(measured, per_gene.effect.to_numpy(dtype=np.float64))
    significant = np.zeros(AXIS_SIZE, dtype=np.int64)
    significant[per_gene.loc[per_gene.padj < PADJ_THRESHOLD, 'gene_index'].to_numpy(dtype=np.int64)] = 1
    results['gse244740'] = {
        'rank': rank, 'measured': len(measured), 'significant': significant,
        'per_gene': per_gene, 'n_contrasts': int(gse_frame.contrast.nunique()),
        'n_samples': gse_provenance['n_wells'],
        'loss': {'n_contrast_rows': len(gse_frame),
                 'n_rows_resolved': int((indices >= 0).sum()),
                 'n_symbols_unresolved': int((indices < 0).sum())},
        'provenance': gse_provenance, 'metadata': {}}
    if intermediates is not None:
        write_frame(per_gene[['gene_index', 'contrast', 'log2fc', 'padj', 'effect',
                              'significant', 'direction', 'n_contrasts_supporting']],
                    intermediates / 'sources' / 'gse244740_deg.csv')
    return results


SUMMARY_COLUMNS = ('source', 'accession', 'category', 'species', 'tissue_model',
                   'contrast', 'n_samples', 'n_contrasts', 'n_axis_expression',
                   'n_axis_measured', 'axis_coverage_share', 'n_hvg', 'hvg_share',
                   'hvg_k', 'deg_up', 'deg_down', 'deg_up_share', 'deg_down_share',
                   'deg_share', 'method', 'normalization', 'ortholog_chain', 'fused',
                   'exploratory', 'supplementary_only', 'notes')


def _source_summary_row(spec, result, axis_size=AXIS_SIZE):
    row = {'source': spec['id'], 'accession': spec['accession'],
           'category': spec['category'], 'species': spec['species'],
           'tissue_model': spec['tissue_model'], 'contrast': spec['contrast'],
           'method': spec['method'], 'ortholog_chain': spec['ortholog_chain'],
           'fused': spec['fused'], 'exploratory': spec['exploratory'],
           'supplementary_only': spec['supplementary_only'], 'notes': spec['notes']}
    if result is None:
        for key in ('n_samples', 'n_contrasts', 'n_axis_expression', 'n_axis_measured',
                    'axis_coverage_share', 'n_hvg', 'hvg_share', 'hvg_k',
                    'deg_up', 'deg_down', 'deg_up_share', 'deg_down_share', 'deg_share'):
            row[key] = None
        row['normalization'] = spec.get('normalization', '')
        return pd.Series(row, dtype=object)
    per_gene = result.get('per_gene')
    if per_gene is None:
        significant_up = significant_down = None
    else:
        significant_up = int(((per_gene.padj < PADJ_THRESHOLD) & (per_gene.log2fc > 0)).sum())
        significant_down = int(((per_gene.padj < PADJ_THRESHOLD) & (per_gene.log2fc < 0)).sum())
    measured = result['measured']
    hvg = result.get('hvg')
    hvg_set = result.get('hvg_set')
    n_hvg = int(hvg_set.sum()) if hvg_set is not None else None
    n_expression = result.get('n_axis_expression', measured)
    row.update({
        'n_samples': result.get('n_samples'),
        'n_contrasts': result.get('n_contrasts'),
        'n_axis_expression': n_expression,
        'n_axis_measured': measured,
        'axis_coverage_share': measured / axis_size,
        'n_hvg': n_hvg,
        'hvg_share': n_hvg / n_expression if n_hvg is not None and n_expression else None,
        'hvg_k': SOURCE_HVG_K if n_hvg is not None else None,
        'deg_up': significant_up, 'deg_down': significant_down,
        'deg_up_share': significant_up / measured if significant_up is not None and measured else None,
        'deg_down_share': significant_down / measured if significant_down is not None and measured else None,
        'deg_share': ((significant_up + significant_down) / measured
                      if significant_up is not None and measured else None),
        'normalization': spec.get('normalization') or result.get('metadata', {}).get('normalization', ''),
    })
    return pd.Series(row, dtype=object)


def _direction_arrays(results):
    """Per-gene winning direction, |log2FC| and exact-tie count across DEG sources."""
    order = [source for source in DEG_IDS if source in results]
    effects = np.zeros((AXIS_SIZE, len(order)), dtype=np.float64)
    log2fc = np.zeros((AXIS_SIZE, len(order)), dtype=np.float64)
    significant = np.zeros((AXIS_SIZE, len(order)), dtype=bool)
    for column, name in enumerate(order):
        per_gene = results[name]['per_gene']
        index = per_gene.gene_index.to_numpy(dtype=np.int64)
        effects[index, column] = per_gene.effect.to_numpy(dtype=np.float64)
        log2fc[index, column] = per_gene.log2fc.to_numpy(dtype=np.float64)
        significant[index, column] = per_gene.padj.to_numpy(dtype=np.float64) < PADJ_THRESHOLD
    best = np.argmax(effects, axis=1)
    best_effect = effects[np.arange(AXIS_SIZE), best]
    best_log2fc = log2fc[np.arange(AXIS_SIZE), best]
    ties = (np.isclose(effects, best_effect[:, None], rtol=0.0, atol=1e-12)).sum(axis=1) - 1
    direction = np.where(best_effect > 0, np.sign(best_log2fc).astype(np.int64), 0)
    abs_effect = np.where(best_effect > 0, np.abs(best_log2fc), 0.0)
    return direction, abs_effect, ties, significant.any(axis=1)


def _overlap_rows(fused, results, tenv_rank):
    """Long overlap table: tenv sets, source-HVG sets, DEG pairs and the intersection."""
    fused_rank = fused.fused_rank.to_numpy(dtype=np.int64)
    sets = {}
    for k in TENV_HVG_KS:
        sets[f'tenv_hvg_{k}'] = set(np.flatnonzero(hvg_membership(tenv_rank, k)))
    for name in EXPRESSION_IDS:
        if name in results:
            sets[f'{name}_hvg_{SOURCE_HVG_K}'] = set(
                np.flatnonzero(results[name]['hvg_set']))
    for name in DEG_IDS:
        if name in results:
            sets[f'{name}_deg'] = set(np.flatnonzero(results[name]['significant']))
    rows = []

    def add(kind, left, right):
        metrics = overlap_metrics(sets[left], sets[right], fused_rank)
        rows.append({'comparison': kind, 'set_a': left, 'set_b': right,
                     'k_a': _set_k(left), 'k_b': _set_k(right), **metrics})

    def _set_k(name):
        if name.startswith('tenv_hvg_'):
            return int(name.rsplit('_', 1)[1])
        if '_hvg_' in name:
            return int(name.rsplit('_', 1)[1])
        return ''

    for k in TENV_HVG_KS:
        for name in DEG_IDS:
            if name in results:
                add('tenv_hvg_vs_source_deg', f'tenv_hvg_{k}', f'{name}_deg')
    for name in EXPRESSION_IDS:
        if name in results:
            add('tenv_hvg_vs_source_hvg', f'tenv_hvg_{SOURCE_HVG_K}',
                f'{name}_hvg_{SOURCE_HVG_K}')
    names = [name for name in DEG_IDS if name in results]
    for left_index, left in enumerate(names):
        for right in names[left_index + 1:]:
            add('source_pair_deg', f'{left}_deg', f'{right}_deg')
    for label, members in (
            ('all_source_deg', names),
            ('all_source_deg_without_gse108940',
             [name for name in names if name != 'gse108940']),
            ('core_source_deg',
             [name for name in ('gse141910', 'gse116250', 'gse243674', 'dtoxs_ipsc')
              if name in names])):
        intersection = (set.intersection(*[sets[f'{name}_deg'] for name in members])
                        if members else set())
        metrics = overlap_metrics(intersection, intersection, fused_rank)
        rows.append({'comparison': label, 'set_a': 'intersection',
                     'set_b': 'intersection', 'k_a': '', 'k_b': '',
                     'n_a': len(intersection), 'n_b': len(intersection),
                     'n_intersection': len(intersection), 'n_union': len(intersection),
                     'jaccard': 1.0,
                     'fused_rank_min': metrics['fused_rank_min'],
                     'fused_rank_median': metrics['fused_rank_median'],
                     'fused_rank_max': metrics['fused_rank_max']})
    return pd.DataFrame(rows), sets


def _answers(fused, results, sets, direction, significant):
    """Numbers that answer #175's three questions, pinned in the manifest."""
    fused_rank = fused.fused_rank.to_numpy(dtype=np.int64)

    def jaccard(left, right):
        a, b = sets[left], sets[right]
        return len(a & b) / len(a | b) if (a | b) else 0.0

    hf = jaccard('gse141910_deg', 'gse116250_deg')
    remodelling = jaccard('gse128188_deg', 'gse282618_deg')
    cross = [jaccard(f'{left}_deg', f'{right}_deg')
             for left in ('gse141910', 'gse116250') for right in ('gse128188', 'gse282618')]
    tenv = {name: len(sets[f'tenv_hvg_{SOURCE_HVG_K}'] & sets[f'{name}_deg'])
            for name in DEG_IDS if f'{name}_deg' in sets}
    up = significant & (direction > 0)
    down = significant & (direction < 0)
    count = fused.n_sources_deg_significant.to_numpy(dtype=np.int64)
    gene_index = np.arange(AXIS_SIZE)
    count_order = np.lexsort((gene_index, -count))
    count_rank = np.empty(AXIS_SIZE, dtype=np.int64)
    count_rank[count_order] = np.arange(1, AXIS_SIZE + 1)
    order = fused_order(fused_rank)
    agreement = {}
    for size in DIRECTION_SIZES:
        overlap = len(set(order[:size]) & set(count_order[:size]))
        agreement[str(size)] = {'n_overlap': int(overlap), 'share': overlap / size}
    return {
        'q1_chronic_remodelling_vs_hf': {
            'hf_pair_jaccard': hf,
            'remodelling_pair_jaccard': remodelling,
            'cross_mean_jaccard': float(np.mean(cross)),
            'tenv_hvg_2000_intersection': tenv,
            'interpretation': ('HF-source agreement is compared with the AF/HHD '
                               '(chronic-remodelling) pair and the four cross pairs.'),
        },
        'q2_up_down_rank_symmetry': {
            'n_up_significant': int(up.sum()), 'n_down_significant': int(down.sum()),
            'median_fused_rank_up': float(np.median(fused_rank[up])) if up.any() else None,
            'median_fused_rank_down': float(np.median(fused_rank[down])) if down.any() else None,
            'top_n_composition': {},
        },
        'q3_effect_replication_vs_deg_count': {
            'n_genes_supported_by_any_deg': int((count > 0).sum()),
            'n_genes_supported_by_2plus': int((count >= 2).sum()),
            'top_n_agreement': agreement,
            'spearman_fused_vs_count_rank': _spearman(fused_rank, count_rank),
        },
    }


def _spearman(left, right):
    """Spearman correlation with average ranks for ties (deterministic)."""
    left, right = pd.Series(left), pd.Series(right)
    return float(left.rank(method='average').corr(right.rank(method='average')))


def _figure_rows(fused, summary_rows, curve, inflection, overlap_matrix):
    """Complete plotted values for the four-panel summary figure."""
    rows = [{'record': 'figure', 'title': (
        'Gene-axis decision evidence · per-source HVG/DEG support, fused rank-sum '
        'knee, tenv non-HVG share and source overlap')}]
    rows += [{'record': 'panel', 'panel': panel, 'panel_title': title}
             for panel, title in (('a', 'Per-source HVG and DEG shares'),
                                  ('b', 'Fused rank-sum distribution and candidate inflection'),
                                  ('c', 'Tenv non-HVG share vs N'),
                                  ('d', 'DEG / tenv-HVG source overlap (Jaccard)'))]
    for row in summary_rows:
        if not row.get('fused'):
            continue
        rows.append({'record': 'bar', 'panel': 'a', 'source': row['source'],
                     'hvg_share': row.get('hvg_share'),
                     'deg_up_share': row.get('deg_up_share'),
                     'deg_down_share': row.get('deg_down_share')})
    scores = fused.rank_sum.to_numpy(dtype=np.float64)
    counts, edges = np.histogram(scores, bins=HISTOGRAM_BINS)
    for position, count in enumerate(counts):
        rows.append({'record': 'hist', 'panel': 'b',
                     'bin_left': float(edges[position]),
                     'bin_right': float(edges[position + 1]),
                     'count': int(count)})
    if inflection['inflection'] and inflection['n']:
        ordered = fused.sort_values('fused_rank')
        marker_score = float(ordered.rank_sum.iloc[inflection['n'] - 1])
        note = f'Kneedle+second-difference candidate N={inflection["n"]}'
    else:
        marker_score = None
        note = 'no inflection (Kneedle and second-difference criteria do not agree)'
    rows.append({'record': 'marker', 'panel': 'b',
                 'marker_n': inflection['n'], 'marker_score': marker_score,
                 'note': note})
    chosen = sorted({1, *range(CURVE_SAMPLE_STEP, len(fused), CURVE_SAMPLE_STEP),
                     len(fused), *([inflection['n']] if inflection['n'] else [])})
    curve_indexed = curve.set_index('n')
    for k in TENV_HVG_KS:
        for n in chosen:
            rows.append({'record': 'curve', 'panel': 'c', 'k': k, 'n': int(n),
                         'non_hvg_share': float(curve_indexed.loc[n, f'non_hvg_share_k{k}'])})
    for cell in overlap_matrix:
        rows.append({'record': 'matrix', 'panel': 'd',
                     'row_source': cell['row_source'], 'col_source': cell['col_source'],
                     'jaccard': cell['jaccard'],
                     'n_intersection': cell['n_intersection']})
    order = ['record', 'panel', 'panel_title', 'title', 'source', 'row_source',
             'col_source', 'k', 'n', 'hvg_share', 'deg_up_share', 'deg_down_share',
             'bin_left', 'bin_right', 'count', 'marker_n', 'marker_score',
             'non_hvg_share', 'jaccard', 'n_intersection', 'note']
    frame = pd.DataFrame(rows)
    return frame.reindex(columns=[column for column in order if column in frame.columns])


def _overlap_matrix(fused_rank, sets):
    names = [f'tenv_hvg_{SOURCE_HVG_K}'] + [f'{name}_deg' for name in DEG_IDS]
    rows = []
    for left in names:
        for right in names:
            metrics = overlap_metrics(sets[left], sets[right], fused_rank)
            rows.append({'row_source': left, 'col_source': right,
                         'jaccard': metrics['jaccard'],
                         'n_intersection': metrics['n_intersection']})
    return rows


def build(config):
    """Compute and write the complete evidence product; returns its manifest."""
    for role, path in config.input_roles():
        forbidden_input(path)
        if not Path(path).exists():
            raise FileNotFoundError(f'missing input {role}: {path}')
    axis, resolver = load_resolver(config)
    orthologs = read_csv(config.orthologs, dtype={'mouse_symbol': str, 'human_symbol': str})
    tenv = read_csv(config.tenv_hvg)
    tenv = tenv.sort_values('gene_index')
    tenv_rank = tenv['rank'].to_numpy(dtype=np.int64)
    if sorted(tenv_rank.tolist()) != list(range(1, AXIS_SIZE + 1)):
        raise ValueError('published tenv HVG ranking is not a complete 1..18164 ranking')
    intermediates = config.intermediates
    results = {'tenv': {'rank': tenv_rank, 'measured': AXIS_SIZE,
                        'significant': np.zeros(AXIS_SIZE, dtype=np.int64),
                        'n_axis_expression': AXIS_SIZE, 'n_samples': 4206,
                        'n_contrasts': 0, 'metadata': {},
                        'hvg_set': hvg_membership(tenv_rank, SOURCE_HVG_K)}}
    for spec in SOURCE_SPECS:
        if spec['id'] in EXPRESSION_IDS:
            results[spec['id']] = process_expression_source(
                config, axis, resolver, spec, intermediates)
        elif spec['id'] in ANIMAL_IDS:
            results[spec['id']] = process_animal_source(
                config, axis, resolver, spec, intermediates, orthologs)
    results.update(process_pert_sources(config, axis, resolver, intermediates,
                                        intermediates / 'sources' if intermediates else None))

    ranks = {spec['id']: results[spec['id']]['rank'] for spec in SOURCE_SPECS
             if spec['id'] in FUSED_IDS}
    measured = {spec['id']: results[spec['id']]['measured'] for spec in SOURCE_SPECS
                if spec['id'] in FUSED_IDS}
    flags = {spec['id']: results[spec['id']]['significant'] for spec in SOURCE_SPECS
             if spec['id'] in DEG_IDS}
    fused = fuse_ranks(ranks, measured, axis_size=AXIS_SIZE, flags=flags)
    fused = fused.rename(columns={spec['id']: spec['column'] for spec in SOURCE_SPECS
                                  if spec['id'] in FUSED_IDS})
    fused.insert(1, 'symbol', axis.symbol.to_numpy())
    for k in TENV_HVG_KS:
        fused[f'tenv_hvg_member_k{k}'] = hvg_membership(tenv_rank, k).astype(np.int64)
    fused['tenv_hvg_member'] = fused[f'tenv_hvg_member_k{SOURCE_HVG_K}']
    columns = ['gene_index', 'symbol', 'tenv_hvg_rank']
    columns += [spec['column'] for spec in SOURCE_SPECS
                if spec['fused'] and spec['id'] != 'tenv']
    columns += ['rank_sum', 'rank_sum_animal_half', 'fused_rank', 'n_sources_measured',
                'tenv_hvg_member', 'n_sources_deg_significant', 'fused_rank_animal_half']
    columns += [f'tenv_hvg_member_k{k}' for k in TENV_HVG_KS]
    fused = fused.reindex(columns=columns)
    fused = fused.sort_values(['fused_rank', 'gene_index']).reset_index(drop=True)

    overlap, sets = _overlap_rows(fused, results, tenv_rank)
    direction, abs_effect, ties, significant_any = _direction_arrays(results)
    composition = direction_composition(
        fused_order(fused.fused_rank.to_numpy(dtype=np.int64)), direction, abs_effect,
        DIRECTION_SIZES, tie_counts=ties)

    ordered = fused.sort_values('fused_rank')
    x = ordered.fused_rank.to_numpy(dtype=np.int64)
    y = ordered.rank_sum.to_numpy(dtype=np.float64)
    inflection = detect_inflection(x, y)
    leave_one_out = []
    for dropped in FUSED_IDS:
        subset = {name: value for name, value in ranks.items() if name != dropped}
        subset_measured = {name: value for name, value in measured.items() if name != dropped}
        subset_flags = {name: value for name, value in flags.items() if name != dropped}
        alternative = fuse_ranks(subset, subset_measured, axis_size=AXIS_SIZE,
                                 flags=subset_flags)
        alternative_order = alternative.sort_values('fused_rank')
        found = detect_inflection(alternative_order.fused_rank.to_numpy(dtype=np.int64),
                                  alternative_order.rank_sum.to_numpy(dtype=np.float64))
        leave_one_out.append({'dropped_source': dropped, 'inflection': found['inflection'],
                              'n': found['n']})
    inflection['curve'] = {'x': INFLECTION_X, 'y': INFLECTION_Y, 'n': int(len(x))}
    inflection['leave_one_source_out'] = leave_one_out

    curve = tenv_non_hvg_curve(fused_order(fused.fused_rank.to_numpy(dtype=np.int64)),
                               tenv_rank)
    summary_rows = [_source_summary_row(spec, results.get(spec['id']))
                    for spec in SOURCE_SPECS]
    summary = pd.DataFrame(summary_rows, columns=list(SUMMARY_COLUMNS))
    answers = _answers(fused, results, sets, direction, significant_any)
    for size in DIRECTION_SIZES:
        rows = composition[(composition.n == size) & (composition.rule == 'symmetric')]
        answers['q2_up_down_rank_symmetry']['top_n_composition'][str(size)] = {
            'n_up': int(rows.n_up.iloc[0]), 'n_down': int(rows.n_down.iloc[0]),
            'up_share': float(rows.up_share.iloc[0]),
            'down_share': float(rows.down_share.iloc[0])}
    matrix = _overlap_matrix(fused.fused_rank.to_numpy(dtype=np.int64), sets)
    figure = _figure_rows(fused, summary_rows, curve, inflection, matrix)

    outputs = {
        'source_summary.csv': summary,
        'fused_rank.csv': fused,
        'overlap.csv': overlap,
        'tenv_non_hvg_curve.csv': curve,
        'direction_summary.csv': composition,
    }
    config.output.mkdir(parents=True, exist_ok=True)
    output_digests = {}
    for name, frame in outputs.items():
        path = write_frame(frame, config.output / name)
        output_digests[f'output/gene-axis/{name}'] = {
            'sha256': sha256(path), 'rows': int(len(frame))}
    inflection_path = write_json(inflection, config.output / 'inflection.json')
    output_digests['output/gene-axis/inflection.json'] = {
        'sha256': sha256(inflection_path), 'rows': 1}
    figure_path = write_frame(figure, config.figure_data)
    output_digests[f'figure/data/{Path(config.figure_data).name}'] = {
        'sha256': sha256(figure_path), 'rows': int(len(figure))}

    manifest = {
        'schema': SCHEMA,
        'axis': {'size': AXIS_SIZE, 'gene_index': '0-based frozen row order',
                 'genes_sha256': sha256(config.genes)},
        'recipe': RECIPE,
        'constants': {'tenv_hvg_ks': list(TENV_HVG_KS), 'source_hvg_k': SOURCE_HVG_K,
                      'padj_threshold': PADJ_THRESHOLD, 'animal_half_weight': 0.5,
                      'direction_sizes': list(DIRECTION_SIZES),
                      'fused_sources': [spec['id'] for spec in SOURCE_SPECS if spec['fused']],
                      'animal_sources': list(ANIMAL_SOURCES)},
        'inputs': {}, 'outputs': output_digests,
        'sources': {}, 'counts': {}, 'answers': answers,
        'environment': environment_facts(),
        'ortholog_generator': {
            'script': 'stage2/01_bases/execution/build_mouse_orthologs.R',
            'sha256': sha256(step_path(config.root, 'execution/build_mouse_orthologs.R')),
            'command': ('R_HOME=/home/chang/miniconda3/envs/arc/lib/R '
                        '/home/chang/miniconda3/envs/arc/lib/R/bin/exec/R --no-save '
                        '--slave --file=stage2/01_bases/execution/build_mouse_orthologs.R'),
            'r_version': 'R version 4.4.1 (2024-06-14)',
            'homologene_version': '1.4.68.19.3.27',
            'org.Mm.eg.db_version': '3.20.0'},
    }
    for role, path in config.input_roles():
        manifest['inputs'][repo_relative(path, config.root)] = {
            'role': role, 'sha256': sha256(path)}
    for spec in SOURCE_SPECS:
        result = results.get(spec['id'])
        if result is None:
            continue
        entry = {'measured': result['measured'],
                 'n_axis_expression': result.get('n_axis_expression', result['measured']),
                 'n_samples': result.get('n_samples'),
                 'n_contrasts': result.get('n_contrasts')}
        if 'loss' in result:
            entry['alignment'] = result['loss']
        manifest['sources'][spec['id']] = entry
    manifest['counts'] = {
        'n_axis': AXIS_SIZE,
        'n_fused_sources': sum(1 for spec in SOURCE_SPECS if spec['fused']),
        'n_deg_sources': len(DEG_IDS),
        'n_genes_supported_by_any_deg': int((fused.n_sources_deg_significant > 0).sum()),
        'n_genes_supported_by_all_deg': int(
            (fused.n_sources_deg_significant == len(DEG_IDS)).sum()),
        'overlap_rows': int(len(overlap)),
        'direction_rows': int(len(composition)),
        'source_summary_rows': int(len(summary)),
        'tenv_non_hvg_rows': int(len(curve)),
        'fused_rank_rows': int(len(fused)),
        'inflection': inflection['inflection'],
        'inflection_n': inflection['n'],
    }
    write_json(manifest, config.output / 'manifest.json')
    return manifest


def step_path(root, relative):
    return Path(root) / 'stage2/01_bases' / relative


def repo_relative(path, root):
    """Lexical repo-relative path; symlinks are not resolved (worktree inputs).

    Manifest entries must be repo-relative so ``--verify`` works both in this
    worktree (where large inputs are read-only symlinks to the canonical
    checkout) and in the canonical checkout itself.
    """
    import os
    lexical = Path(os.path.abspath(path))
    try:
        return lexical.relative_to(Path(os.path.abspath(root))).as_posix()
    except ValueError:
        return lexical.as_posix()


def environment_facts():
    import platform
    import scipy
    return {
        'python': platform.python_version(),
        'numpy': np.__version__,
        'pandas': pd.__version__,
        'scipy': scipy.__version__,
        'data_access': 'offline; no h5ad read, no truth or Heliste file opened',
        'resource_caps': 'ulimit -v 8G, OMP_NUM_THREADS=2, run under /usr/bin/time -v',
    }


def verify(config):
    """Rebuild the evidence in a temporary directory and check every digest."""
    stored = json.loads((config.output / 'manifest.json').read_text())
    if stored.get('schema') != SCHEMA:
        raise ValueError('unknown gene-axis evidence schema')
    for path, info in stored['inputs'].items():
        current = sha256(Path(config.root) / path)
        if current != info['sha256']:
            raise ValueError(f'input digest differs: {path}')
    output_digests = stored['outputs']
    with tempfile.TemporaryDirectory() as temporary:
        temporary = Path(temporary)
        replay_output = temporary / 'gene-axis'
        replay_figure = temporary / 'figure' / Path(config.figure_data).name
        replay = replace(config, output=replay_output, figure_data=replay_figure,
                         intermediates=None)
        rebuilt = build(replay)
    differences = []
    for relative, info in output_digests.items():
        if relative.startswith('output/gene-axis/'):
            path = config.output / relative.split('output/gene-axis/')[1]
        elif relative.startswith('figure/data/'):
            path = config.figure_data
        else:
            raise ValueError(f'unknown evidence output: {relative}')
        if sha256(path) != info['sha256']:
            differences.append(relative)
    rebuilt.pop('environment')
    stored.pop('environment')
    if rebuilt != stored:
        differences.append('manifest.json')
    if differences:
        raise ValueError(f'gene-axis replay differs: {sorted(set(differences))}')
    return {'result': 'PASS', 'schema': SCHEMA,
            'inputs': len(stored['inputs']), 'outputs': len(output_digests),
            'manifest_sha256': sha256(config.output / 'manifest.json'),
            'inflection': stored['counts']['inflection'],
            'inflection_n': stored['counts']['inflection_n'],
            'n_axis': stored['counts']['n_axis']}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--root', type=Path)
    parser.add_argument('--genes', type=Path)
    parser.add_argument('--crosswalk', type=Path)
    parser.add_argument('--hgnc', type=Path)
    parser.add_argument('--tenv-hvg', type=Path)
    parser.add_argument('--orthologs', type=Path)
    parser.add_argument('--disease-root', type=Path)
    parser.add_argument('--pert-root', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--figure-data', type=Path)
    parser.add_argument('--intermediates', type=Path)
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--record', type=Path)
    args = parser.parse_args()
    config = default_config(args.root, args.output, args.figure_data, args.intermediates)
    overrides = {'genes': args.genes, 'crosswalk': args.crosswalk, 'hgnc': args.hgnc,
                 'tenv_hvg': args.tenv_hvg, 'orthologs': args.orthologs,
                 'disease_root': args.disease_root, 'pert_root': args.pert_root}
    config = replace(config, **{key: value for key, value in overrides.items() if value})
    result = verify(config) if args.verify else build(config)
    text = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.record:
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
