#!/usr/bin/env python
"""GSE238242 — CM-pseudobulk AF (persistent) vs sinus rhythm DEG (ADR 0030).

The item's CM pseudobulk (4 sinus-rhythm / 3 persistent-AF donors, 3,049 author-
labelled CM nuclei) is fitted with the one live implementation,
`stage2/01_bases/output/api/bases/disease_deg.R`, in counts mode: limma-voom on
the integer pseudobulk. The declared design is the single AF-vs-SR contrast.

The expected and measured outcome is **zero** v2 DEGs (padj < 0.05 and
|log2FC| >= 2): n = 3 vs 4 donors cannot support FDR-significant differential
expression here, and the item publishes that honest null rather than loosening
the threshold or pooling in non-CM cells.

Usage (arc python, after script/build_pseudobulk.py):
    /home/chang/miniconda3/envs/arc/bin/python script/build_deg.py
"""
from __future__ import annotations

import gzip
import hashlib
import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ITEM = Path(__file__).resolve().parents[1]
REPO = ITEM
while not (REPO / '.git').exists():
    if REPO.parent == REPO:
        raise SystemExit('cannot locate the repository root')
    REPO = REPO.parent
sys.path.insert(0, str(REPO / 'stage2/01_bases/execution'))
import public_api  # noqa: F401  (binds the single public implementation)
from bases import disease_deg  # noqa: E402

LOCAL = ITEM / 'local'
PSEUDOBULK = 'local/pseudobulk/GSE238242_cm_pseudobulk.csv.gz'
PROFILES = 'local/pseudobulk/GSE238242_profiles.csv'
SERIES = 'source/GSE238242_series_matrix.txt.gz'
METADATA = 'source/GSE238242_snAF.metadata.tsv.gz'
PADJ_THRESHOLD = 0.05
DEG_LOG2FC_THRESHOLD = 2.0
FILTER = {'min_count': 10, 'min_total_count': 15, 'large_n': 10, 'min_prop': 0.7}
SPEC = {
    'item': 'GSE238242',
    'category': 'AF',
    'tier': 1,
    'evidence_type': 'cm_pseudobulk',
    'model_system': ('human left-atrial-appendage 10x multiome snRNA, author CM label, '
                     'per-donor CM count pseudobulk'),
    'platform': 'GPL24676 (Illumina NovaSeq 6000)',
    'expected_rna_libraries': 7,
    'expected_total_libraries': 14,
    'contrasts': [{'label': 'AF-SR', 'numerator': 'AF', 'denominator': 'SR'}],
    'series_conclusion': ('the series title/overall design state persistent AF vs sinus '
                          'rhythm, left atrial appendage; the deposit has no paroxysmal arm'),
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write_tsv_gz(frame, path):
    payload = frame.to_csv(sep='\t', index=False, lineterminator='\n').encode()
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode='wb', mtime=0) as handle:
        handle.write(payload)
    Path(path).write_bytes(buffer.getvalue())


def read_series_title(path):
    import gzip as gz
    import re
    with gz.open(path, 'rt', encoding='utf-8-sig') as handle:
        for line in handle:
            if line.startswith('!Series_title'):
                return line.partition('\t')[2].strip().strip('"')
            if line.startswith('!series_matrix_table_begin'):
                break
    raise ValueError(f'{path}: no series title')


def deg_counts(table):
    per_contrast = {}
    for label, rows in table.groupby('contrast'):
        flag = (rows.padj < PADJ_THRESHOLD) & (rows.log2fc.abs() >= DEG_LOG2FC_THRESHOLD)
        per_contrast[str(label)] = {
            'n_deg': int(flag.sum()),
            'n_up': int((flag & (rows.log2fc > 0)).sum()),
            'n_down': int((flag & (rows.log2fc < 0)).sum()),
            'n_padj_only': int((rows.padj < PADJ_THRESHOLD).sum()),
            'n_padj_lt_010': int((rows.padj < 0.10).sum()),
            'n_rows': int(len(rows)),
            'min_padj': float(rows.padj.min()),
        }
    work = table.assign(deg=(table.padj < PADJ_THRESHOLD)
                        & (table.log2fc.abs() >= DEG_LOG2FC_THRESHOLD))
    winner = work.sort_values(['gene_symbol', 'padj', 'p_value'],
                              ascending=[True, True, True], kind='stable'
                              ).drop_duplicates('gene_symbol', keep='first')
    collapsed = winner[winner.deg]
    return {
        'definition': (f'padj < {PADJ_THRESHOLD} AND |log2FC| >= {DEG_LOG2FC_THRESHOLD} '
                       f'in the declared contrast (ADR 0071)'),
        'per_contrast': per_contrast,
        'collapsed': {
            'n_genes': int(len(winner)),
            'n_deg': int(len(collapsed)),
            'n_up': int(collapsed.log2fc.gt(0).sum()),
            'n_down': int(collapsed.log2fc.lt(0).sum()),
            'n_padj_only_genes': int(table.loc[table.padj < PADJ_THRESHOLD,
                                               'gene_symbol'].nunique()),
            'deg_genes': collapsed.gene_symbol.tolist(),
        },
    }


def main():
    for relative in (PSEUDOBULK, PROFILES):
        if not (ITEM / relative).is_file():
            raise SystemExit(f'run script/build_pseudobulk.py first ({relative} missing)')
    item = SPEC['item']
    local = LOCAL / 'deg'
    local.mkdir(parents=True, exist_ok=True)
    pseudobulk = pd.read_csv(ITEM / PSEUDOBULK)
    profiles = pd.read_csv(ITEM / PROFILES)
    if list(pseudobulk.columns[:1]) != ['gene']:
        raise ValueError('the pseudobulk matrix must start with a gene column')
    if len(pseudobulk.columns) - 1 != SPEC['expected_rna_libraries']:
        raise ValueError(f'{len(pseudobulk.columns) - 1} donor profiles, expected '
                         f'{SPEC["expected_rna_libraries"]}')
    if pseudobulk['gene'].duplicated().any():
        raise ValueError('duplicated gene symbol in the pseudobulk matrix')
    profiles = profiles.set_index('profile')
    if sorted(profiles.loc[pseudobulk.columns[1:], 'rhythm']) != ['AF', 'AF', 'AF',
                                                                  'SR', 'SR', 'SR', 'SR']:
        raise ValueError('the profiles do not carry 3 AF / 4 SR donors')
    genes = pseudobulk['gene'].astype(str)
    expression = pd.DataFrame({'source_gene_id': genes, 'gene_symbol': genes})
    sheet = pd.DataFrame({'sample': pseudobulk.columns[1:],
                          'group': [profiles.loc[donor, 'rhythm']
                                    for donor in pseudobulk.columns[1:]],
                          'donor': pseudobulk.columns[1:],
                          'cm_nuclei': [int(profiles.loc[donor, 'cm_nuclei'])
                                        for donor in pseudobulk.columns[1:]],
                          'sex': [str(profiles.loc[donor, 'sex'])
                                  for donor in pseudobulk.columns[1:]]})
    for donor in pseudobulk.columns[1:]:
        expression[donor] = pseudobulk[donor].to_numpy(dtype=np.int64)
    values = expression.iloc[:, 2:].to_numpy()
    if (values < 0).any() or (values != np.round(values)).any():
        raise ValueError('the pseudobulk is not integer counts')
    expr_path = local / f'{item}_expr.tsv.gz'
    samples_path = local / f'{item}_samples.csv'
    write_tsv_gz(expression, expr_path)
    sheet.to_csv(samples_path, index=False)
    design = {'kind': 'counts', 'formula': '~0 + group', 'factor': 'group',
              'contrasts': SPEC['contrasts'], 'filter': FILTER}
    fit_summary = disease_deg.run_deg_worker(expr_path, samples_path, design,
                                             local / f'{item}_deg.csv')
    table = pd.read_csv(local / f'{item}_deg.csv')
    missing = set(disease_deg.LOCAL_DEG_COLUMNS) - set(table.columns)
    if missing:
        raise ValueError(f'DEG worker output lacks {sorted(missing)}')
    counts = deg_counts(table)
    pseudobulk_summary = json.loads((ITEM / 'local/pseudobulk/'
                                     'GSE238242_cm_pseudobulk.json').read_text())
    record = {
        'schema': 'stage2-item-deg-summary-v1',
        'item': item,
        'category': SPEC['category'],
        'tier': SPEC['tier'],
        'evidence_type': SPEC['evidence_type'],
        'model_system': SPEC['model_system'],
        'series': {
            'title': read_series_title(ITEM / SERIES),
            'platform': SPEC['platform'],
            'tissue': 'left atrial appendage',
            'total_libraries': SPEC['expected_total_libraries'],
            'rna_libraries': SPEC['expected_rna_libraries'],
            'atac_libraries_not_fetched': (SPEC['expected_total_libraries']
                                           - SPEC['expected_rna_libraries']),
            'excluded_libraries_reason': (
                'the ATAC arm (~130 MB of the 187 MB RAW.tar) is not part of an RNA '
                'DEG source and was never downloaded'),
            'design': 'single nuclei from left atrial appendages, persistent AF vs sinus rhythm',
            'author_cell_annotation': pseudobulk_summary['author_annotation'],
            'conclusion': SPEC['series_conclusion'],
        },
        'deposition': {
            'series_matrix': {'path': SERIES,
                              'bytes': int((ITEM / SERIES).stat().st_size),
                              'sha256': sha256(ITEM / SERIES)},
            'metadata': {'path': METADATA,
                         'bytes': int((ITEM / METADATA).stat().st_size),
                         'sha256': sha256(ITEM / METADATA)},
            'rna_tables': pseudobulk_summary['deposition']['rna_tables'],
        },
        'pseudobulk': {
            'path': PSEUDOBULK,
            'sha256': pseudobulk_summary['pseudobulk']['matrix_sha256'],
            'rule': pseudobulk_summary['pseudobulk']['rule'],
            'genes': int(len(genes)),
            'profiles': int(len(sheet)),
            'profiles_path': PROFILES,
            'profiles_sha256': pseudobulk_summary['pseudobulk']['profiles_sha256'],
        },
        'samples': {
            'path': str(samples_path.relative_to(ITEM)),
            'sha256': sha256(samples_path),
            'n': int(len(sheet)),
            'group_counts': {str(key): int(value)
                             for key, value in sheet.group.value_counts().items()},
            'cm_nuclei': {str(row.sample): int(row.cm_nuclei) for row in sheet.itertuples()},
            'match_basis': 'author cell table `sample` column joined to the file names',
        },
        'design': {'kind': 'counts', 'formula': '~0 + group', 'factor': 'group',
                   'contrasts': SPEC['contrasts'], 'filter': FILTER},
        'contrasts': {'AF-SR': {'n_case': 3, 'n_control': 4,
                                'numerator': 'AF (persistent)',
                                'denominator': 'SR (sinus rhythm)',
                                'direction': 'positive log2FC = up in persistent AF'}},
        'n_case': 3,
        'n_control': 4,
        'deg_counts': counts,
        'null_outcome_statement': (
            'the measured v2 DEG count is zero: no gene reaches padj < 0.05 in the '
            '3-vs-4-donor contrast, so the item publishes the null explicitly instead of '
            'loosening the threshold, counting padj-only genes, or adding non-CM cells'),
        'sensitivity': {
            'audit_filter_variants': {
                'method': ('the pre-existing AF source audit\'s filter variants on the '
                           'same CM pseudobulk (limma on log2(CPM+1), BH): CPM >= 1 in '
                           '>= 3/7, >= 7/7 and trend=TRUE all give 0 genes at padj<0.05'),
                'padj_lt_005': 0,
                'min_padj_by_filter': {'ge_3_of_7': 0.220, 'ge_7_of_7': 0.150,
                                       'trend': 0.177},
                'provenance': '.scratch/af-source-audit/GSE238242/REPORT.md (untracked '
                              'audit of 2026-11-14)',
                'reading': 'the null is not a filtering artefact',
            },
        },
        'method': {
            'deposition': 'raw integer counts (dense TSV, genes x nuclei)',
            'pseudobulk': pseudobulk_summary['pseudobulk']['rule'],
            'model': "limma-voom (voom's own library-size normalisation, moderated t)",
            'filter': 'edgeR filterByExpr rule written out in the live worker',
            'implementation': 'stage2/01_bases/output/api/bases/disease_deg.R',
            'implementation_sha256': sha256(Path(disease_deg.DEG_WORKER)),
            'derivation': str(Path(__file__).resolve().relative_to(REPO)),
            'derivation_sha256': sha256(Path(__file__).resolve()),
            'test': 'moderated t, BH-adjusted p-value',
            'limma_version': fit_summary['limma_version'],
            'r_version': fit_summary['r_version'],
        },
        'fit_summary': fit_summary,
        'deg': {'path': str((local / f'{item}_deg.csv').relative_to(ITEM)),
                'sha256': sha256(local / f'{item}_deg.csv'),
                'rows': int(len(table)),
                'n_source_genes': int(table.gene_symbol.nunique()),
                'bytes': int((local / f'{item}_deg.csv').stat().st_size)},
        'expr': {'path': str(expr_path.relative_to(ITEM)), 'sha256': sha256(expr_path)},
    }
    summary_path = local / f'{item}_deg.json'
    summary_path.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'result': 'PASS', 'item': item,
                      'n_donors': len(sheet),
                      'n_deg': counts['collapsed']['n_deg'],
                      'min_padj': counts['per_contrast']['AF-SR']['min_padj'],
                      'summary': str(summary_path.relative_to(ITEM))}, indent=2,
                     sort_keys=True))
    return 0


if __name__ == '__main__':
    sys.exit(main())
