#!/usr/bin/env python
"""GSE255612 — one AF-vs-control limma fit on the donor CM pseudobulk.

The deposit is normalised expression (log1p(CP10K), verified in
build_pseudobulk.py), so the DEG is a documented linear model, not a count
model: the per-donor mean linear-CP10K pseudobulk is put on the log2 scale as
log2(mean_CP10K + 1), one joint design ``~ 0 + group`` is fitted over the 34
donors, and the AF-vs-control contrast is read from that fit (limma moderated t).

Inputs (declared ``local/`` additions):

    local/pseudobulk/gse255612_cm_pseudobulk.csv.gz
    local/pseudobulk/gse255612_donors.csv

Outputs:

    local/deg/gse255612_deg.csv    long table (one contrast; gitignored, regenerable)
    local/deg/gse255612_deg.json   tracked summary (design/method/digests)

Usage (arc python):
    /home/chang/miniconda3/envs/arc/bin/python script/build_deg.py
"""
from __future__ import annotations

import hashlib
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
PSEUDOBULK = LOCAL / 'pseudobulk' / 'gse255612_cm_pseudobulk.csv.gz'
DONORS = LOCAL / 'pseudobulk' / 'gse255612_donors.csv'
DEG = LOCAL / 'deg' / 'gse255612_deg.csv'
MATRIX = LOCAL / 'deg' / 'gse255612_log_matrix.tsv.gz'
SAMPLES = LOCAL / 'deg' / 'gse255612_samples.csv'
SUMMARY = LOCAL / 'deg' / 'gse255612_deg.json'
SCHEMA = 'stage2-item-deg-summary-v1'
DESIGN = {
    'kind': 'expression',
    'formula': '~0 + group',
    'factor': 'group',
    'contrasts': [{'label': 'AF-control', 'numerator': 'AF', 'denominator': 'control'}],
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    if not PSEUDOBULK.is_file() or not DONORS.is_file():
        raise SystemExit('run script/build_pseudobulk.py first')
    donors = pd.read_csv(DONORS)
    groups = donors.group.map(lambda value: 'AF' if str(value).lower() == 'case' else 'control')
    if set(groups) != {'AF', 'control'}:
        raise SystemExit(f'unexpected AF groups: {sorted(set(groups))}')
    samples = pd.DataFrame({'sample': donors.donor_id, 'group': groups,
                            'donor': donors.donor_id})
    LOCAL.joinpath('deg').mkdir(parents=True, exist_ok=True)
    samples.to_csv(SAMPLES, index=False)
    pseudobulk = pd.read_csv(PSEUDOBULK)
    values = pseudobulk[donors.donor_id.tolist()].to_numpy(dtype=np.float64)
    if not np.isfinite(values).all() or (values < 0).any():
        raise SystemExit('pseudobulk has non-finite or negative values')
    logged = pd.DataFrame({
        'source_gene_id': pseudobulk.source_gene_id,
        'gene_symbol': pseudobulk.gene_symbol,
    })
    for position, donor in enumerate(donors.donor_id):
        logged[donor] = np.log2(values[:, position] + 1.0)
    logged.to_csv(MATRIX, index=False, compression={'method': 'gzip', 'mtime': 0}, sep='\t')
    summary = disease_deg.run_deg_worker(MATRIX, SAMPLES, DESIGN, DEG)
    table = pd.read_csv(DEG)
    missing = set(disease_deg.LOCAL_DEG_COLUMNS) - set(table.columns)
    if missing:
        raise SystemExit(f'R output lacks {sorted(missing)}')
    record = {
        'schema': SCHEMA,
        'item': 'GSE255612',
        'category': 'AF',
        'gate': 'cell_type == "Cardiomyocytes"',
        'pseudobulk': {
            'path': 'local/pseudobulk/gse255612_cm_pseudobulk.csv.gz',
            'sha256': sha256(PSEUDOBULK),
            'basis': ('per-donor mean over CM cells of the linear CP10K value; the '
                      'deposit is log1p(CP10K), expm1 recovers the linear value'),
            'cells_per_donor': {row.donor_id: int(row.cells_cm)
                                for row in donors.itertuples()},
        },
        'model_matrix': {
            'path': 'local/deg/gse255612_log_matrix.tsv.gz',
            'sha256': sha256(MATRIX),
            'transform': 'log2(mean_CP10K + 1)',
        },
        'design': DESIGN,
        'method': {
            'model': 'limma on the log2-scale pseudobulk (moderated t)',
            'implementation': 'stage2/01_bases/output/api/bases/disease_deg.R',
            'implementation_sha256': sha256(Path(disease_deg.DEG_WORKER)),
            'test': 'moderated t, BH-adjusted p-value',
            'deposition': 'normalised-expression deposit (log1p(CP10K)), not integer counts',
            'contrast_direction': 'positive log2FC = up in AF vs control',
            'limma_version': summary['limma_version'],
            'r_version': summary['r_version'],
        },
        'contrasts': {'AF-control': {
            'n_case': int((samples.group == 'AF').sum()),
            'n_control': int((samples.group == 'control').sum())}},
        'fit_summary': summary,
        'deg': {'path': 'local/deg/gse255612_deg.csv', 'sha256': sha256(DEG),
                'rows': int(len(table)), 'n_source_genes': int(table.source_gene_id.nunique())},
    }
    SUMMARY.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'result': 'PASS', 'rows': len(table),
                      'n_case': record['contrasts']['AF-control']['n_case'],
                      'n_control': record['contrasts']['AF-control']['n_control'],
                      'deg_sha256': record['deg']['sha256']}, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
