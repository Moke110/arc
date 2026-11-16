"""Extract the complete author CM pseudobulk DESeq2 table, without refitting.

The immutable workbook is the author test result, not a per-cell annotation.
Iterating its one sheet bounds memory independently of the large GEO R objects.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import openpyxl

ITEM = Path(__file__).resolve().parent.parent

def digest(path):
    return hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()

def build(source, out):
    out.mkdir(parents=True, exist_ok=True)
    book = openpyxl.load_workbook(source, read_only=True, data_only=True)
    rows = book['Cardiomyocytes'].iter_rows(values_only=True)
    if list(next(rows)) != ['Gene', 'log2FoldChange', 'padj']:
        raise ValueError('Unexpected author CM table columns')
    fields = ['gene_symbol', 'source_gene_id', 'contrast', 'log2fc', 'stat',
              'p_value', 'padj', 'base_mean', 'n_case', 'n_control']
    target = out / 'GSE183852_deg.csv'
    seen = set()
    n_up = n_down = n_null = 0
    with target.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        for gene, effect, adjusted in rows:
            if not isinstance(gene, str) or not gene.strip() or gene in seen:
                raise ValueError('Missing/duplicate author gene')
            if not isinstance(effect, (int, float)) or not math.isfinite(effect):
                raise ValueError('Missing/nonfinite author log2FoldChange')
            if adjusted is not None and (not isinstance(adjusted, (int, float)) or
                                         not math.isfinite(adjusted) or not 0 <= adjusted <= 1):
                raise ValueError('Invalid adjusted p value')
            seen.add(gene)
            n_null += adjusted is None
            passed = adjusted is not None and adjusted < 0.05 and abs(effect) >= 2
            n_up += passed and effect > 0
            n_down += passed and effect < 0
            writer.writerow(dict(gene_symbol=gene, source_gene_id=gene,
                contrast='DCM-Donor', log2fc=effect, stat='', p_value='',
                padj='' if adjusted is None else adjusted, base_mean='',
                n_case=13, n_control=25))
    book.close()
    summary = {
        'series': 'GSE183852', 'species': 'Homo sapiens', 'category': 'nHF',
        'platform': '10x human LV snRNA-seq; published cardiomyocyte pseudobulk DESeq2',
        'method': 'Author DESeq2 statistics preserved; no local raw-count refit',
        'contrasts': {'DCM-Donor': {'n_case': 13, 'n_control': 25,
                      'numerator': 'DCM', 'denominator': 'Donor'}},
        'deposition': {'path': 'source/' + source.name, 'sha256': digest(source),
                       'bytes': source.stat().st_size, 'sheet': 'Cardiomyocytes'},
        'deg': {'path': 'local/deg/' + target.name, 'rows': len(seen),
                'n_source_genes': len(seen), 'bytes': target.stat().st_size,
                'sha256': digest(target)},
        'deg_counts': {'definition': 'padj < 0.05 AND |log2FC| >= 2',
                       'collapsed': {'n_genes': len(seen), 'n_deg': n_up + n_down,
                                     'n_up': n_up, 'n_down': n_down}},
        'missing_padj': n_null,
        'limitations': [
            '13 DCM/25 donor are snRNA libraries from GEO; no re-estimation of donor independence.',
            'No genotype stratification or raw-count test refit was performed.',
            'Workbook contains no per-cell labels or unadjusted p values/stat/base mean; absent fields remain empty.',
            'Cell-level R objects remain outside this pinned author-statistic adapter.'],
    }
    (out / 'GSE183852_deg.json').write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=ITEM / 'source/44161_2022_28_MOESM4_ESM.xlsx')
    parser.add_argument('--out', type=Path, default=ITEM / 'local/deg')
    args = parser.parse_args()
    build(args.source, args.out)
