"""Hermetic fixture for the disease-DEG publication tests (ticket #187).

Builds a tiny but contract-shaped repository tree: the frozen 18,164-gene axis,
a crosswalk that carries one alias row and one ambiguous row, three fake external
items with their pinned manifests, and one item-local DEG table per fixture
dataset. Every value is hand-worked in the tests; the fixture only drives the
real alignment/collapse/rank/reader boundary.

The fixture dataset ids deliberately mirror the ticket-187 identifiers so the
tests exercise the same contrast labels (`NICM-Donor`, `ICM-Donor`,
`AMI-Donor`, `AF-control`, `dcm_genotype_negative-normal`) without any network.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

AXIS_SIZE = 18164


def _axis():
    return pd.DataFrame({
        'gene_index': list(range(AXIS_SIZE)),
        'symbol': [f'G{index:05d}' for index in range(AXIS_SIZE)],
        'hgnc_id': [f'HGNC:{index + 1}' for index in range(AXIS_SIZE)],
        'entrez_id': [str(index + 1) for index in range(AXIS_SIZE)],
        'ensembl_id': [f'ENSG{index:011d}' for index in range(AXIS_SIZE)],
        'match_type': ['direct'] * AXIS_SIZE,
    })


def _crosswalk():
    # gene_index 1 gets an alias + a previous symbol; an ambiguous alias maps to
    # two on-axis genes and must stay unresolved.
    return pd.DataFrame([
        {'dataset_id': 'fixture', 'source_gene_id': 'OLD_G1', 'symbol': 'OLD_G1',
         'hgnc_id': 'HGNC:6', 'match_type': 'alias', 'entrez_verified': 'True',
         'in_corpus': 'False'},
        {'dataset_id': 'fixture', 'source_gene_id': 'PREV_G2', 'symbol': 'PREV_G2',
         'hgnc_id': 'HGNC:3', 'match_type': 'previous', 'entrez_verified': 'True',
         'in_corpus': 'False'},
        {'dataset_id': 'fixture', 'source_gene_id': 'AMB', 'symbol': 'AMB',
         'hgnc_id': 'HGNC:4', 'match_type': 'alias', 'entrez_verified': 'True',
         'in_corpus': 'False'},
        {'dataset_id': 'fixture', 'source_gene_id': 'AMB', 'symbol': 'AMB',
         'hgnc_id': 'HGNC:5', 'match_type': 'alias', 'entrez_verified': 'True',
         'in_corpus': 'False'},
    ])


def _row(gene_symbol, source_gene_id, contrast, log2fc, stat, p_value, padj, base_mean):
    return {'source_gene_id': source_gene_id, 'gene_symbol': gene_symbol,
            'contrast': contrast, 'log2fc': log2fc, 'stat': stat,
            'p_value': p_value, 'padj': padj, 'base_mean': base_mean}


def gse270788_deg():
    """One joint-design table: all three arms come out of a single fit."""
    rows = []
    for contrast, scale in (('NICM-Donor', 1.0), ('ICM-Donor', -1.5), ('AMI-Donor', 0.75)):
        rows += [
            _row('G00001', 'ENSG00000000001', contrast, 2.0 * scale, 5.0, 1e-6,
                 1e-4, 8.0),
            _row('G00002', 'ENSG00000000002', contrast, -1.0 * scale, -3.0, 1e-3,
                 0.01, 6.0),
            # A missing padj (independent filtering) is non-significant, not dropped.
            _row('G00003', 'ENSG00000000003', contrast, 0.5, 0.4, None, None, 4.0),
            # An alias symbol resolves through the crosswalk to gene_index 5.
            _row('OLD_G1', 'ENSG00000000005', contrast, 1.5 * scale, 4.0, 1e-4,
                 0.001, 7.0),
            # Off-axis/unmapped rows must be reported, never guessed.
            _row('NOT_ON_AXIS', 'ENSOFFAXIS', contrast, 3.0, 1.0, 0.01, 0.05, 1.0),
            # Ambiguous alias: two on-axis candidates -> unresolved.
            _row('AMB', 'ENSG00000000999', contrast, 1.0, 1.0, 0.01, 0.05, 1.0),
            # A second source id that aliases to gene_index 2: the max-effect row
            # wins within the contrast (ties settle on input order).
            _row('PREV_G2', 'ENSG00000000002', contrast, 0.1, 0.2, 0.5, 0.9, 2.0),
        ]
    return pd.DataFrame(rows)


def reichart_deg():
    return pd.DataFrame([
        _row('G00001', 'ENSG00000000001', 'dcm_genotype_negative-normal', -2.0, -4.0,
             1e-5, 1e-3, 9.0),
        _row('G00002', 'ENSG00000000002', 'dcm_genotype_negative-normal', 0.5, 1.0,
             0.1, 0.3, 5.0),
        _row('G00004', 'ENSG00000000004', 'dcm_genotype_negative-normal', 1.0, 2.0,
             0.01, 0.05, 3.0),
    ])


def gse255612_deg():
    # Two genes share the same effect score: the gene_index tie-break decides.
    return pd.DataFrame([
        _row('G00010', 'ENSG00000000010', 'AF-control', 2.0, 3.0, 1e-4, 1e-3, 6.0),
        _row('G00011', 'ENSG00000000011', 'AF-control', -2.0, -3.0, 1e-4, 1e-3, 6.0),
        _row('G00012', 'ENSG00000000012', 'AF-control', 1.0, 1.0, 0.5, None, 4.0),
    ])


DATASETS = [
    {
        'dataset_id': 'gse270788__nicm', 'category': 'nHF',
        'source_item': 'disease_data/GSE270788',
        'local_deg': 'disease_data/GSE270788/local/deg/gse270788_deg.csv',
        'local_summary': 'disease_data/GSE270788/local/deg/gse270788_deg.json',
        'cells_per_donor': 'disease_data/GSE270788/local/pseudobulk/gse270788_donors.csv',
        'arm_or_class': 'NICM arm vs six Donor controls',
        'species': 'Homo sapiens', 'model_system': 'fixture CM pseudobulk',
        'evidence_type': 'cm_pseudobulk', 'tier': 1,
        'genotype_status': 'not_declared', 'time_window_basis': None,
        'exploratory': False, 'shared_control_group': 'gse270788:donor',
        'n_case': 7, 'n_control': 6, 'contrasts': ['NICM-Donor'],
        'method': {'deposition': 'integer counts', 'model': 'limma-voom',
                   'design': 'one joint fit ~0 + arm', 'test': 'moderated t',
                   'padj': 'BH'},
    },
    {
        'dataset_id': 'gse270788__icm', 'category': 'ICM',
        'source_item': 'disease_data/GSE270788',
        'local_deg': 'disease_data/GSE270788/local/deg/gse270788_deg.csv',
        'local_summary': 'disease_data/GSE270788/local/deg/gse270788_deg.json',
        'cells_per_donor': 'disease_data/GSE270788/local/pseudobulk/gse270788_donors.csv',
        'arm_or_class': 'ICM arm vs six Donor controls',
        'species': 'Homo sapiens', 'model_system': 'fixture CM pseudobulk',
        'evidence_type': 'cm_pseudobulk', 'tier': 1,
        'genotype_status': 'not_declared', 'time_window_basis': None,
        'exploratory': False, 'shared_control_group': 'gse270788:donor',
        'n_case': 5, 'n_control': 6, 'contrasts': ['ICM-Donor'],
        'method': {'deposition': 'integer counts', 'model': 'limma-voom',
                   'design': 'one joint fit ~0 + arm', 'test': 'moderated t',
                   'padj': 'BH'},
    },
    {
        'dataset_id': 'gse270788__ami', 'category': 'ICM',
        'source_item': 'disease_data/GSE270788',
        'local_deg': 'disease_data/GSE270788/local/deg/gse270788_deg.csv',
        'local_summary': 'disease_data/GSE270788/local/deg/gse270788_deg.json',
        'cells_per_donor': 'disease_data/GSE270788/local/pseudobulk/gse270788_donors.csv',
        'arm_or_class': 'AMI arm vs six Donor controls',
        'species': 'Homo sapiens', 'model_system': 'fixture CM pseudobulk',
        'evidence_type': 'cm_pseudobulk', 'tier': 1,
        'genotype_status': 'not_declared', 'time_window_basis': 'author_declared',
        'exploratory': False, 'shared_control_group': 'gse270788:donor',
        'n_case': 5, 'n_control': 6, 'contrasts': ['AMI-Donor'],
        'method': {'deposition': 'integer counts', 'model': 'limma-voom',
                   'design': 'one joint fit ~0 + arm', 'test': 'moderated t',
                   'padj': 'BH'},
    },
    {
        'dataset_id': 'reichart-cm__dcm_genotype_negative', 'category': 'nHF',
        'source_item': 'disease_data/Reichart-CM',
        'local_deg': 'disease_data/Reichart-CM/local/deg/reichart-cm_deg.csv',
        'local_summary': 'disease_data/Reichart-CM/local/deg/reichart-cm_deg.json',
        'cells_per_donor': 'disease_data/Reichart-CM/local/genotype/reichart-cm_donor_genotypes.csv',
        'arm_or_class': 'genotype-negative DCM vs normal donors',
        'species': 'Homo sapiens', 'model_system': 'fixture CM pseudobulk',
        'evidence_type': 'cm_pseudobulk', 'tier': 1, 'genotype_status': 'negative',
        'time_window_basis': None, 'exploratory': False,
        'shared_control_group': 'reichart-cm:normal',
        'n_case': 30, 'n_control': 18,
        'contrasts': ['dcm_genotype_negative-normal'],
        'method': {'deposition': 'integer counts', 'model': 'limma-voom',
                   'design': 'one joint fit ~0 + status', 'test': 'moderated t',
                   'padj': 'BH'},
    },
    {
        'dataset_id': 'gse255612__af', 'category': 'AF',
        'source_item': 'disease_data/GSE255612',
        'local_deg': 'disease_data/GSE255612/local/deg/gse255612_deg.csv',
        'local_summary': 'disease_data/GSE255612/local/deg/gse255612_deg.json',
        'cells_per_donor': 'disease_data/GSE255612/local/pseudobulk/gse255612_donors.csv',
        'arm_or_class': 'AF vs non-AF controls',
        'species': 'Homo sapiens', 'model_system': 'fixture CM pseudobulk',
        'evidence_type': 'cm_pseudobulk', 'tier': 1,
        'genotype_status': 'not_declared', 'time_window_basis': None,
        'exploratory': False, 'shared_control_group': 'gse255612:donor',
        'n_case': 18, 'n_control': 16, 'contrasts': ['AF-control'],
        'method': {'deposition': 'normalised expression', 'model': 'limma',
                   'design': 'one joint fit ~0 + group', 'test': 'moderated t',
                   'padj': 'BH'},
    },
]

REJECTED = [{'candidate': 'GSE141910', 'reason': 'whole-tissue bulk, not CM-derived'}]


def write_fixture(root):
    """Write the fixture repo tree under ``root``; return the sources module path."""
    root = Path(root)
    identity = root / 'stage2/01_bases/output/gene-identity'
    (identity / 'tables').mkdir(parents=True, exist_ok=True)
    _axis().to_csv(identity / 'genes.csv', index=False)
    _crosswalk().to_csv(identity / 'tables/gene_crosswalk.csv', index=False)

    tables = {
        'disease_data/GSE270788/local/deg/gse270788_deg.csv': gse270788_deg(),
        'disease_data/Reichart-CM/local/deg/reichart-cm_deg.csv': reichart_deg(),
        'disease_data/GSE255612/local/deg/gse255612_deg.csv': gse255612_deg(),
    }
    for relative, table in tables.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(path, index=False)
    items = sorted({record['source_item'] for record in DATASETS})
    for item in items:
        base = root / item
        (base / 'manifest.sha256').parent.mkdir(parents=True, exist_ok=True)
        (base / 'manifest.sha256').write_text('00  fixture\n')
    for record in DATASETS:
        for key in ('local_summary', 'cells_per_donor'):
            path = root / record[key]
            path.parent.mkdir(parents=True, exist_ok=True)
            if not path.exists():
                path.write_text('fixture\n')
    sources = root / 'fixture_sources'
    sources.write_text(_sources_module(), encoding='utf-8')
    return root, sources


def _sources_module():
    categories = {
        'nHF': {'directory': 'nhf', 'tier': 1, 'evidence_type': 'cm_pseudobulk',
                'shared_control_group': ['gse270788:donor'], 'rejected': REJECTED},
        'ICM': {'directory': 'icm', 'tier': 1, 'evidence_type': 'cm_pseudobulk',
                'shared_control_group': ['gse270788:donor'], 'rejected': []},
        'AF': {'directory': 'af', 'tier': 1, 'evidence_type': 'cm_pseudobulk',
               'shared_control_group': ['gse255612:donor'], 'rejected': []},
    }
    return (
        '"""Fixture source declaration generated by disease_deg_fixtures.py."""\n'
        f'CATEGORIES = {categories!r}\n'
        f'DATASETS = {DATASETS!r}\n'
        'def datasets(category):\n'
        '    return [record for record in DATASETS if record["category"] == category]\n'
        'def category_record(category):\n'
        '    return CATEGORIES[category]\n'
    )
