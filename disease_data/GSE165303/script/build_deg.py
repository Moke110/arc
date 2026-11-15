#!/usr/bin/env python
"""GSE165303 — human left-ventricle bulk RNA-seq, DCM vs non-failing (tier 3).

The item-local derivation: the deposited RNA-seq master table (Ensembl gene x
101 libraries, on the log2 normalised-count scale, with the authors' own DESeq2
result columns) and the GPL20795 series matrix (the deposit's own sample block)
are bound here and the one live DEG implementation
``bases.disease_deg.R`` (limma on the deposited log2 matrix) is run through
``bases.disease_deg.run_deg_worker``. Nothing below that line re-fits a model.

Verified design (from the deposit, not from the ticket prose):

* the series profiled 177 libraries across HiChIP/Hi-C/ChIP-seq/ATAC-seq/RNA-seq;
  the **RNA-seq arm is 101 bulk left-ventricle libraries: 50 DCM + 51
  non-failing** (the ticket prose "101 bulk DCM samples" is wrong: 101 is the
  total and 50 are DCM).
* samples are biobanked human **left-ventricular** tissue (Duke Human Heart
  Repository; ``!Sample_growth_protocol_ch1``), non-failing donor LV as control.
* the deposited values are **log2(normalised counts)** (the paper's Figure 1D
  legend; 2**value tracks the deposited DESeq2 ``baseMean`` at r = 0.999999),
  not raw counts, and they carry a hard floor at 3.821687527 (= 14.14 on the
  linear scale) applied to ~35% of cells. The exact upstream size-factor /
  pseudocount formula is not stated in the GEO record; the deposit is modelled
  as-is with limma ``lmFit`` (the documented normalised-expression path).
* no clinical covariates (age, sex, cause) are deposited, so the design is the
  single ``~0 + group`` contrast ``DCM-NF``; every library is used.

Usage (arc python):
    /home/chang/miniconda3/envs/arc/bin/python script/build_deg.py
"""
from __future__ import annotations

import gzip
import hashlib
import io
import json
import re
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

MASTER_TABLE = 'source/GSE165303_RNA-Seq_master_table.txt.gz'
SERIES_MATRIX = 'source/GSE165303-GPL20795_series_matrix.txt.gz'
ID_COLUMN = 'ENSEMBLE'
AUTHOR_STAT_COLUMNS = ('baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj')
# The deposit also carries the authors' Entrez gene symbol in a trailing column;
# it is annotation, not a library.
SYMBOL_COLUMN = 'Symbol'
SAMPLE_EXCLUDED_COLUMNS = (ID_COLUMN,) + AUTHOR_STAT_COLUMNS + (SYMBOL_COLUMN,)
PADJ_THRESHOLD = 0.05
DEG_LOG2FC_THRESHOLD = 2.0
EXPECTED_RNA_TITLES = 101
EXPECTED_CASE = 50
EXPECTED_CONTROL = 51

SPEC = {
    'item': 'GSE165303',
    'category': 'nHF',
    'tier': 3,
    'evidence_type': 'bulk',
    'master_table': MASTER_TABLE,
    'series_matrix': SERIES_MATRIX,
    'id_column': ID_COLUMN,
    'expected_tissue': 'heart',
    'expected_rna_samples': EXPECTED_RNA_TITLES,
    'expected_case': EXPECTED_CASE,
    'expected_control': EXPECTED_CONTROL,
    'kind': 'expression',
    'formula': '~0 + group',
    'contrasts': [
        {'label': 'DCM-NF', 'numerator': 'DCM', 'denominator': 'NF'},
    ],
    'contrast_direction': 'positive log2FC = up in DCM left-ventricle tissue vs non-failing',
    'group_map': {'DCM': 'DCM', 'NF': 'NF'},
    'strata': {
        'decision': 'single case/control contrast; the deposit carries no clinical covariate',
        'pooled': False,
        'blocking_covariates': [],
        'why': ('the GPL20795 series matrix declares only tissue (heart), assay (RNA-Seq) '
                'and the DCM/NF sample title; no age, sex or cause field is deposited, so '
                'no covariate can be blocked without importing external metadata'),
    },
    'model_basis': (
        'whole left-ventricle tissue bulk: not cardiomyocyte-specific, so it is a tier-3 '
        'source only (ADR 0071) and enters nHF after tiers 1-2 fall short of 2,000 DEGs'),
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_series_matrix(path):
    """Return ``(series, sample_fields, characteristic_rows)`` from the GPL20795 series matrix."""
    series = {}
    fields = {}
    characteristics = []
    with gzip.open(path, 'rt', encoding='utf-8-sig') as handle:
        for line in handle:
            line = line.rstrip('\r\n')
            if line.startswith('!series_matrix_table_begin'):
                break
            if line.startswith('!Series_'):
                key, _, value = line.partition('\t')
                series.setdefault(key[len('!Series_'):].lower(), []).append(
                    value.strip().strip('"'))
                continue
            if line.startswith('!Sample_characteristics_ch1\t'):
                characteristics.append([value.strip().strip('"')
                                        for value in line.split('\t')[1:]])
                continue
            match = re.match(r'^!Sample_([A-Za-z0-9_]+)\t(.+)$', line)
            if match:
                fields.setdefault(match.group(1).lower(), []).append(
                    [value.strip().strip('"') for value in match.group(2).split('\t')])
    return series, fields, characteristics


def write_tsv_gz(frame, path):
    """Write a reproducible TSV gzip: fixed mtime and the declared column order."""
    payload = frame.to_csv(sep='\t', index=False, lineterminator='\n').encode()
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode='wb', mtime=0) as handle:
        handle.write(payload)
    Path(path).write_bytes(buffer.getvalue())


def deg_counts(table):
    """v2 DEG counts for the single declared contrast."""
    per_contrast = {}
    for label, rows in table.groupby('contrast'):
        flag = (rows.padj < PADJ_THRESHOLD) & (rows.log2fc.abs() >= DEG_LOG2FC_THRESHOLD)
        per_contrast[str(label)] = {
            'n_deg': int(flag.sum()),
            'n_up': int((flag & (rows.log2fc > 0)).sum()),
            'n_down': int((flag & (rows.log2fc < 0)).sum()),
            'n_padj_only': int((rows.padj < PADJ_THRESHOLD).sum()),
            'n_rows': int(len(rows)),
        }
    return {
        'definition': (f'padj < {PADJ_THRESHOLD} AND |log2FC| >= '
                       f'{DEG_LOG2FC_THRESHOLD} (ADR 0071)'),
        'per_contrast': per_contrast,
        'collapsed': {'n_deg': per_contrast[SPEC['contrasts'][0]['label']]['n_deg']},
    }


def characteristic_row(characteristics, key):
    """One ``key: value`` characteristic row, selected by its key (not by position)."""
    for row in characteristics:
        parsed = [value.split(':', 1) for value in row]
        if all(len(part) == 2 and part[0].strip() == key for part in parsed):
            return [part[1].strip() for part in parsed]
    raise ValueError(f'the series matrix has no {key!r} characteristic row')


def rna_titles(fields, characteristics):
    """The deposit's RNA-seq samples: sample title plus its assay characteristic."""
    titles = fields['title'][0]
    accessions = fields['geo_accession'][0]
    method = characteristic_row(characteristics, 'method')
    tissue = characteristic_row(characteristics, 'tissue')
    if len(titles) != len(accessions) or len(method) != len(titles):
        raise ValueError('the series matrix sample block is inconsistent')
    rna = [(accessions[index], titles[index], method[index], tissue[index])
           for index in range(len(titles)) if method[index] == 'RNA-Seq']
    if not rna:
        raise ValueError('no RNA-Seq sample in the GPL20795 series matrix')
    if {value for _gsm, _title, _method, value in rna} != {SPEC['expected_tissue']}:
        raise ValueError('the RNA-seq samples are not the verified tissue')
    return rna


def sample_from_title(title):
    """``DCM-11_RNA`` -> ``DCM-11``: the master table column name."""
    match = re.fullmatch(r'(DCM-\d+|NF-\d+)_RNA', str(title))
    if not match:
        raise ValueError(f'unexpected RNA-seq sample title: {title!r}')
    return match.group(1)


def main():
    for relative in (SPEC['master_table'], SPEC['series_matrix']):
        if not (ITEM / relative).is_file():
            raise SystemExit(f'run script/fetch.sh first ({relative} missing)')
    item = SPEC['item']
    local = ITEM / 'local' / 'deg'
    local.mkdir(parents=True, exist_ok=True)
    master_path = ITEM / SPEC['master_table']
    series_path = ITEM / SPEC['series_matrix']
    series, fields, characteristics = read_series_matrix(series_path)
    series_sample_ids = series.get('sample_id', [''])
    series_total = len(str(series_sample_ids[0]).split())
    if series_total == 0:
        raise ValueError('the series matrix lists no samples')
    rna = rna_titles(fields, characteristics)
    if len(rna) != SPEC['expected_rna_samples']:
        raise ValueError(f'expected {SPEC["expected_rna_samples"]} RNA-seq libraries, '
                         f'found {len(rna)}')
    titles = {sample_from_title(title) for _gsm, title, _method, _tissue in rna}
    if len(titles) != len(rna):
        raise ValueError('the RNA-seq sample titles are not unique')
    sheet = pd.DataFrame({
        'sample': [sample_from_title(title) for _gsm, title, _method, _tissue in rna],
        'group': [SPEC['group_map'][sample_from_title(title).split('-')[0]]
                  for _gsm, title, _method, _tissue in rna],
        'gsm': [gsm for gsm, _title, _method, _tissue in rna],
        'title': [title for _gsm, title, _method, _tissue in rna],
    })
    counts = {str(key): int(value) for key, value in sheet.group.value_counts().items()}
    if counts != {'DCM': SPEC['expected_case'], 'NF': SPEC['expected_control']}:
        raise ValueError(f'arm sizes differ from the verified deposit: {counts}')
    master = pd.read_csv(master_path, sep='\t', compression='gzip')
    missing_author = set(AUTHOR_STAT_COLUMNS) - set(master.columns)
    if missing_author:
        raise ValueError(f'master table lacks author statistic columns: '
                         f'{sorted(missing_author)}')
    if SPEC['id_column'] not in master.columns:
        raise ValueError(f'master table lacks the {SPEC["id_column"]} gene-id column')
    missing_symbol = SYMBOL_COLUMN not in master.columns
    if missing_symbol:
        raise ValueError(f'master table lacks the {SYMBOL_COLUMN} gene-symbol column')
    sample_columns = [column for column in master.columns
                      if column not in SAMPLE_EXCLUDED_COLUMNS]
    if sorted(sample_columns) != sorted(sheet['sample']):
        absent = sorted(set(sheet['sample']) - set(sample_columns))
        extra = sorted(set(sample_columns) - set(sheet['sample']))
        raise ValueError('master-table columns do not match the series-matrix RNA '
                         f'samples (absent={absent}, extra={extra})')
    if len(sample_columns) != len(rna):
        raise ValueError('a master-table column repeats a sample')
    # The fit requires the sample-sheet row order to equal the matrix column
    # order: reorder the declared design to the master table's own column order
    # (the binding above is by name, never by position).
    sheet = sheet.set_index('sample').loc[sample_columns].reset_index()
    if master[SPEC['id_column']].duplicated().any():
        raise ValueError('the master table repeats a gene id')
    values = master[sample_columns].to_numpy(dtype=np.float64)
    if not np.isfinite(values).all():
        raise ValueError('the deposited matrix carries non-finite values')
    if values.max() < 10 or values.min() < 0:
        raise ValueError('the deposited matrix does not look like log2 normalised counts')
    floor = float(values.min())
    floor_cells = int((values == floor).sum())
    scale_check = {
        'mode': ('log2(normalised counts) as deposited, with a fixed floor; the exact '
                 'upstream size-factor/pseudocount formula is not documented in GEO'),
        'evidence': ('2**value tracks the deposited DESeq2 baseMean (Pearson r = '
                     '0.999999 over the 5,408 genes with baseMean > 1,000), so the '
                     'columns are log2 normalised counts, not raw counts'),
        'n_values': int(values.size),
        'min': floor,
        'median': float(np.median(values)),
        'max': float(values.max()),
        'floor_value': floor,
        'floor_linear_scale': float(2.0 ** floor),
        'n_floor_cells': floor_cells,
        'floor_fraction': float(floor_cells / values.size),
        'n_genes': int(master.shape[0]),
    }
    # The scale claim is checked against the deposit's own DESeq2 baseMean
    # (Pearson r on the high-expression genes); a wrong declaration fails here.
    high = master.baseMean.to_numpy(dtype=float) > 1_000
    if high.sum() < 100:
        raise ValueError('too few high-expression genes to verify the deposited scale')
    correlation = float(np.corrcoef(
        np.power(2.0, values[high]).mean(axis=1), master.baseMean.to_numpy()[high])[0, 1])
    if correlation < 0.99:
        raise ValueError('the deposited matrix does not track baseMean on the log2 '
                         f'normalised-count scale (r = {correlation})')
    scale_check['high_expression_correlation'] = correlation
    scale_check['high_expression_genes'] = int(high.sum())
    symbols = master[SYMBOL_COLUMN].astype(object).where(
        master[SYMBOL_COLUMN].notna(), '').astype(str).str.strip()
    identifiers = master[SPEC['id_column']].astype(str)
    identifiers_frame = pd.DataFrame({
        'source_gene_id': identifiers.to_numpy(dtype=object),
        'gene_symbol': symbols.where(symbols != '', identifiers).to_numpy(dtype=object),
    })
    expression = pd.concat(
        [identifiers_frame,
         master[sample_columns].astype(np.float64).reset_index(drop=True)], axis=1)
    expr_path = local / f'{item}_expr.tsv.gz'
    samples_path = local / f'{item}_samples.csv'
    write_tsv_gz(expression, expr_path)
    sheet.to_csv(samples_path, index=False)
    design = {'kind': SPEC['kind'], 'formula': SPEC['formula'], 'factor': 'group',
              'contrasts': SPEC['contrasts']}
    deg_path = local / f'{item}_deg.csv'
    fit_summary = disease_deg.run_deg_worker(expr_path, samples_path, design, deg_path)
    table = pd.read_csv(deg_path)
    missing = set(disease_deg.LOCAL_DEG_COLUMNS) - set(table.columns)
    if missing:
        raise ValueError(f'{item}: DEG worker output lacks {sorted(missing)}')
    produced = sorted(table.contrast.astype(str).unique().tolist())
    declared = sorted(contrast['label'] for contrast in SPEC['contrasts'])
    if produced != declared:
        raise ValueError(f'{item}: worker produced {produced}, declared {declared}')
    deg_counts_record = deg_counts(table)
    # Cross-check against the authors' own DESeq2 statistics carried in the table;
    # the fit is ours (one live implementation), this only records concordance.
    author = master.set_index(SPEC['id_column'])
    our = table.set_index('source_gene_id')
    joined = our.join(author[['padj', 'log2FoldChange', 'baseMean']],
                      how='inner', rsuffix='_author')
    if len(joined) != len(table):
        raise ValueError('the DEG table and the master table do not join one-to-one')
    author_v2 = ((joined.padj_author < PADJ_THRESHOLD)
                 & (joined.log2FoldChange.abs() >= DEG_LOG2FC_THRESHOLD))
    our_v2 = ((joined.padj < PADJ_THRESHOLD)
              & (joined.log2fc.abs() >= DEG_LOG2FC_THRESHOLD))
    author_check = {
        'note': ('the deposit carries the authors\' DESeq2 result; the declared DEGs are '
                 'recomputed by this repository\'s one live implementation (limma on the '
                 'deposited log2 matrix), and this block records agreement only'),
        'n_rows_joined': int(len(joined)),
        'author_n_padj_only': int((joined.padj_author < PADJ_THRESHOLD).sum()),
        'author_n_v2': int(author_v2.sum()),
        'limma_n_v2': int(our_v2.sum()),
        'n_v2_shared': int((author_v2 & our_v2).sum()),
        'log2fc_pearson': float(np.corrcoef(joined.log2fc.to_numpy(dtype=float),
                                            joined.log2FoldChange.to_numpy(dtype=float))[0, 1]),
        'log2fc_spearman': float(joined.log2fc.corr(
            joined.log2FoldChange, method='spearman')),
    }
    contrasts = {}
    for contrast in SPEC['contrasts']:
        contrasts[contrast['label']] = {
            'n_case': int((sheet.group == contrast['numerator']).sum()),
            'n_control': int((sheet.group == contrast['denominator']).sum()),
            'numerator': contrast['numerator'],
            'denominator': contrast['denominator'],
        }
    record = {
        'schema': 'stage2-item-deg-summary-v1',
        'item': item,
        'category': SPEC['category'],
        'tier': SPEC['tier'],
        'evidence_type': SPEC['evidence_type'],
        'model_basis': SPEC['model_basis'],
        'series': {
            'title': series.get('title', [''])[0],
            # GPL20795 is "Illumina HiSeq X Ten" (GEO platform record); the
            # pre-#211 string said "Illumina NovaSeq 6000" and was corrected in
            # the #211 docs sweep.
            'platform': 'GPL20795 (Illumina HiSeq X Ten)',
            'assay': 'bulk RNA-seq',
            'total_libraries_in_series': int(series_total),
            'total_libraries_on_gpl20795': int(len(fields['title'][0])),
            'rna_libraries': int(len(rna)),
            'non_rna_libraries': int(len(fields['title'][0]) - len(rna)),
            'tissue': ('human left-ventricle tissue (Duke Human Heart Repository; '
                       '!Sample_growth_protocol_ch1); non-failing donor LV as control'),
            'case_definition': 'DCM (dilated cardiomyopathy) patient LV',
            'control_definition': 'non-failing (NF) donor LV',
        },
        'deposition': {
            'master_table': {'path': SPEC['master_table'],
                             'sha256': sha256(master_path),
                             'bytes': int(master_path.stat().st_size),
                             'gene_rows': int(master.shape[0]),
                             'columns': int(master.shape[1])},
            'series_matrix': {'path': SPEC['series_matrix'],
                              'sha256': sha256(series_path),
                              'bytes': int(series_path.stat().st_size)},
            'values': scale_check,
            'gene_ids': {
                'id_column': SPEC['id_column'],
                'symbol_column': SYMBOL_COLUMN,
                'n_rows': int(master.shape[0]),
                'n_with_symbol': int((symbols != '').sum()),
                'n_without_symbol_falls_back_to_id': int((symbols == '').sum()),
                'alignment_basis': ('Entrez symbol first (the deposited Symbol column), '
                                    'versioned Ensembl id as the fallback'),
            },
        },
        'samples': {
            'path': str(samples_path.relative_to(ITEM)),
            'sha256': sha256(samples_path),
            'n': int(len(sheet)),
            'match_basis': 'series-matrix RNA-seq sample title -> master-table column',
            'group_source': 'sample title prefix (DCM-/NF-)',
            'group_counts': counts,
        },
        'design': {'kind': SPEC['kind'], 'formula': SPEC['formula'], 'factor': 'group',
                   'contrasts': SPEC['contrasts']},
        'contrasts': contrasts,
        'n_case': EXPECTED_CASE,
        'n_control': EXPECTED_CONTROL,
        'strata': SPEC['strata'],
        'deg_counts': deg_counts_record,
        'author_table_check': author_check,
        'method': {
            'deposition': 'log2(normalised counts) as deposited (floored matrix)',
            'model': 'limma on the deposited log2 matrix (moderated t)',
            'implementation': 'stage2/01_bases/output/api/bases/disease_deg.R',
            'implementation_sha256': sha256(Path(disease_deg.DEG_WORKER)),
            'derivation': str(Path(__file__).resolve().relative_to(REPO)),
            'derivation_sha256': sha256(Path(__file__).resolve()),
            'test': 'moderated t, BH-adjusted p-value',
            'contrast_direction': SPEC['contrast_direction'],
            'limma_version': fit_summary['limma_version'],
            'r_version': fit_summary['r_version'],
        },
        'fit_summary': fit_summary,
        'deg': {'path': str(deg_path.relative_to(ITEM)), 'sha256': sha256(deg_path),
                'rows': int(len(table)),
                'n_source_genes': int(table.source_gene_id.nunique()),
                'bytes': int(deg_path.stat().st_size)},
        'expr': {'path': str(expr_path.relative_to(ITEM)), 'sha256': sha256(expr_path)},
    }
    summary_path = local / f'{item}_deg.json'
    summary_path.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'result': 'PASS', 'item': item, 'n_libraries': int(len(sheet)),
                      'n_deg': deg_counts_record['collapsed']['n_deg'],
                      'deg_sha256': record['deg']['sha256'],
                      'summary': str(summary_path.relative_to(ITEM))}, indent=2,
                     sort_keys=True))
    return 0


if __name__ == '__main__':
    sys.exit(main())
