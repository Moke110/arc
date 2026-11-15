#!/usr/bin/env python
"""GSE238242 — cardiomyocyte-nucleus pseudobulk from the snRNA arm (ADR 0030).

GSE238242 (Leblanc et al.) is a paired 10x multiome deposit: 7 snRNA + 7 snATAC
libraries from human left atrial appendages (4 sinus rhythm / 3 persistent AF
donors). ADR 0071's source is the **snRNA arm only**, so the item uses the seven
RNA count tables and never fetches the ~130 MB ATAC arm.

Each RNA table is a dense TSV of raw integer counts, genes x nuclei (38,231 x
2,259 etc.). The author cell table `GSE238242_snAF.metadata.tsv.gz` carries one
row per nucleus with `cell_type`, `sample`, `sex` and `Rhythm`; this script gates
on the author label `cell_type == "CM"` (3,049 of 11,986 nuclei) and sums each
donor's CM counts **before any testing**. The sum is the sufficient statistic for
a count model, so the fit downstream is limma-voom on integer pseudobulk, not a
mean of normalised per-cell values. No cell is re-clustered and no marker-based
gate is invented.

Writes (ADR 0030: the matrix and profile sheet are regenerable `local/`
additions, the JSON summary is tracked):

    local/pseudobulk/GSE238242_cm_pseudobulk.csv.gz  38,231 genes x 7 donors
    local/pseudobulk/GSE238242_profiles.csv          donor -> rhythm / CM nuclei / depth
    local/pseudobulk/GSE238242_cm_pseudobulk.json    tracked summary (digests)

Usage (arc python, after script/fetch.sh && script/verify.sh):
    /home/chang/miniconda3/envs/arc/bin/python script/build_pseudobulk.py
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

SOURCE = ITEM / 'source'
METADATA = SOURCE / 'GSE238242_snAF.metadata.tsv.gz'
LOCAL = ITEM / 'local' / 'pseudobulk'
SCHEMA = 'stage2-item-cm-pseudobulk-v1'
CM_LABEL = 'CM'
EXPECTED_NUCLEI = 11986
EXPECTED_CM_NUCLEI = 3049
EXPECTED_GENES = 38231
# GSM -> donor code, as deposited in the file names and the metadata `sample` column.
LIBRARIES = {
    'GSM7660990': 'CF69', 'GSM7660991': 'CF77', 'GSM7660992': 'CF89',
    'GSM7660993': 'CF91', 'GSM7660994': 'CF93', 'GSM7660995': 'CF97',
    'GSM7660996': 'CF102',
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write_csv_gz(frame, path):
    """Reproducible CSV gzip: fixed mtime so a rebuild is byte-identical."""
    payload = frame.to_csv(index=False, lineterminator='\n').encode()
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode='wb', mtime=0) as handle:
        handle.write(payload)
    Path(path).write_bytes(buffer.getvalue())


def main():
    local = LOCAL
    local.mkdir(parents=True, exist_ok=True)
    for gsm in LIBRARIES:
        candidates = sorted(SOURCE.glob(f'{gsm}_*.counts.tsv.gz'))
        if len(candidates) != 1:
            raise SystemExit(f'run script/fetch.sh first ({gsm} RNA table missing)')
    metadata = pd.read_csv(METADATA, sep='\t', index_col=0)
    if len(metadata) != EXPECTED_NUCLEI:
        raise ValueError(f'metadata has {len(metadata)} nuclei, expected {EXPECTED_NUCLEI}')
    if sorted(metadata.cell_type.unique()) != sorted(
            ['CM', 'FB', 'EC', 'Myeloid', 'Mesothelial', 'Endocardial', 'PC', 'SMC',
             'Lymphoid', 'Neuron', 'Adipo', 'Mast']):
        raise ValueError(f'unexpected author cell types: {sorted(metadata.cell_type.unique())}')
    cm = metadata[metadata.cell_type == CM_LABEL]
    if len(cm) != EXPECTED_CM_NUCLEI:
        raise ValueError(f'the CM gate selects {len(cm)} nuclei, expected {EXPECTED_CM_NUCLEI}')
    if set(cm['sample'].unique()) != set(LIBRARIES.values()):
        raise ValueError('the CM nuclei do not cover exactly the seven RNA donors')
    rhythm = metadata.drop_duplicates('sample').set_index('sample')['Rhythm'].to_dict()

    profiles = []
    sums = {}
    for gsm, donor in LIBRARIES.items():
        path = next(SOURCE.glob(f'{gsm}_*.counts.tsv.gz'))
        # The deposited header row lists the barcodes only (no corner cell): the gene
        # column has no name, so pandas makes it the index and the 2,259 barcodes the
        # columns.
        head = pd.read_csv(path, sep='\t', nrows=1)
        barcodes_in_matrix = head.columns.tolist()
        if len(barcodes_in_matrix) != len(metadata[metadata['sample'] == donor]):
            raise ValueError(f'{gsm}: matrix columns ({len(barcodes_in_matrix)}) disagree '
                             f'with the metadata rows for {donor}')
        barcodes = cm.index[cm['sample'] == donor]
        missing = set(barcodes) - set(barcodes_in_matrix)
        if missing:
            raise ValueError(f'{gsm}: {len(missing)} CM barcodes absent from the matrix')
        if len(set(barcodes_in_matrix)) != len(barcodes_in_matrix):
            raise ValueError(f'{gsm}: duplicated barcode in the matrix header')
        # Read the whole table (38,231 x ~2,259) and select the CM columns: passing
        # every header name as `usecols` trips pandas' ragged-header handling (the
        # corner cell is absent), which silently drops the first barcode column.
        frame = pd.read_csv(path, sep='\t', index_col=0)
        if frame.shape[1] != len(barcodes_in_matrix):
            raise ValueError(f'{gsm}: parsed {frame.shape[1]} columns, expected '
                             f'{len(barcodes_in_matrix)}')
        if frame.shape[0] != EXPECTED_GENES:
            raise ValueError(f'{gsm}: {frame.shape[0]} genes, expected {EXPECTED_GENES}')
        if frame.index.duplicated().any():
            raise ValueError(f'{gsm}: duplicated gene symbol')
        values = frame.to_numpy()
        if (values < 0).any() or (values != np.round(values)).any():
            raise ValueError(f'{gsm}: the deposited table is not raw integer counts')
        summed = frame[barcodes].sum(axis=1).astype(np.int64)
        sums[donor] = summed
        profiles.append({
            'profile': donor, 'gsm': gsm, 'donor': donor,
            'rhythm': 'AF' if rhythm[donor] == 'AF' else 'SR',
            'cm_nuclei': int(len(barcodes)),
            'library_size': int(summed.sum()),
            'genes_detected': int((summed > 0).sum()),
            'nuclei_in_matrix': int(len(barcodes_in_matrix)),
            'sex': str(metadata[metadata['sample'] == donor].sex.iloc[0]),
        })
    genes = sums[LIBRARIES['GSM7660990']].index
    for donor, summed in sums.items():
        if not genes.equals(summed.index):
            raise ValueError(f'{donor}: the RNA tables do not share the gene axis')
    matrix = pd.DataFrame({'gene': genes})
    for gsm, donor in LIBRARIES.items():
        matrix[donor] = sums[donor].to_numpy()
    matrix_path = local / 'GSE238242_cm_pseudobulk.csv.gz'
    write_csv_gz(matrix, matrix_path)
    profile_frame = pd.DataFrame(profiles).sort_values('rhythm', kind='stable')
    profile_path = local / 'GSE238242_profiles.csv'
    profile_frame.to_csv(profile_path, index=False)
    record = {
        'schema': SCHEMA,
        'item': 'GSE238242',
        'arm': 'snRNA only (the 7 snATAC libraries are not fetched or used)',
        'cell_type_column': 'cell_type',
        'cell_type_label': CM_LABEL,
        'author_annotation': {
            'source': 'GSE238242_snAF.metadata.tsv.gz',
            'statement': ('the cell table is an author file shipping beside the count '
                          'archives; the per-cell clustering/integration that produced it '
                          'is not in the deposit (no Seurat/ArchR object, no h5ad), so the '
                          'labels can be trusted but not re-derived or re-checked'),
            'nuclei': int(len(metadata)),
            'cm_nuclei': int(len(cm)),
            'cm_fraction': round(float(len(cm) / len(metadata)), 6),
            'cell_types': sorted(metadata.cell_type.unique().tolist()),
        },
        'pseudobulk': {
            'rule': ('sum of each donor\'s CM raw integer counts (the sufficient statistic '
                     'for a count model); no per-cell normalisation and no mean of '
                     'normalised values, which would add a normalisation choice'),
            'genes': int(matrix.shape[0]),
            'profiles': int(matrix.shape[1] - 1),
            'matrix_path': str(matrix_path.relative_to(ITEM)),
            'matrix_sha256': sha256(matrix_path),
            'matrix_bytes': int(matrix_path.stat().st_size),
            'profiles_path': str(profile_path.relative_to(ITEM)),
            'profiles_sha256': sha256(profile_path),
            'profile_rows': profiles,
        },
        'libraries': {gsm: donor for gsm, donor in LIBRARIES.items()},
        'rhythm_counts': {'SR': int(sum(1 for row in profiles if row['rhythm'] == 'SR')),
                          'AF': int(sum(1 for row in profiles if row['rhythm'] == 'AF'))},
        'deposition': {
            'metadata': {'path': 'source/GSE238242_snAF.metadata.tsv.gz',
                         'bytes': int(METADATA.stat().st_size), 'sha256': sha256(METADATA)},
            'rna_tables': {gsm: {'path': f'source/{path.name}',
                                 'bytes': int(path.stat().st_size),
                                 'sha256': sha256(path)}
                           for gsm, path in ((gsm, next(SOURCE.glob(f'{gsm}_*.counts.tsv.gz')))
                                             for gsm in LIBRARIES)},
        },
    }
    summary_path = local / 'GSE238242_cm_pseudobulk.json'
    summary_path.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'result': 'PASS', 'item': 'GSE238242',
                      'cm_nuclei': int(len(cm)),
                      'profiles': int(matrix.shape[1] - 1),
                      'genes': int(matrix.shape[0]),
                      'summary': str(summary_path.relative_to(ITEM))}, indent=2,
                     sort_keys=True))
    return 0


if __name__ == '__main__':
    sys.exit(main())
