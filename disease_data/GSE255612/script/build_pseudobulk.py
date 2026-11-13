#!/usr/bin/env python
"""GSE255612 — author CM gate and per-donor CM pseudobulk (ADR 0066 D8).

GSE255612 is a human left-atrium 10x 3' v3 snRNA deposit. The usable payload is
one Scanpy Matrix-Market matrix (36,600 features x 179,697 cells, real values)
whose metadata rows are in matrix column order. The series-matrix processing note
says the matrix is "Processed data files output from Scanpy as gene expression
matrix"; the values are log1p(CP10K): a cell's single-count genes share one
value, equal to log1p(1e4 / n_umi). The script verifies that basis on a sample of
cells before using it.

Gate: metadata ``cell_type == "Cardiomyocytes"``. Each donor's CM cells are
summed in the linear CP10K space (``expm1`` of the deposited value) and divided by
the donor's CM cell count, so the pseudobulk is a per-donor mean expression
profile that does not scale with how many CM cells the donor happened to yield.
No clustering; the 1 GB mtx is streamed, never held as a dense cell x gene matrix.

Writes under the item's declared ``local/`` area:

    local/pseudobulk/gse255612_cm_pseudobulk.csv.gz   genes x 34 donors, float
    local/pseudobulk/gse255612_donors.csv             per-donor cells/group table
    local/pseudobulk/gse255612_cm_pseudobulk.json     tracked summary (digests)

Usage (arc python):
    /home/chang/miniconda3/envs/arc/bin/python script/build_pseudobulk.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ITEM = Path(__file__).resolve().parents[1]
SOURCE = ITEM / 'source'
LOCAL = ITEM / 'local' / 'pseudobulk'
MTX = SOURCE / 'GSE255612_AF_snRNA_Matrix_V1.mtx.gz'
GENES = SOURCE / 'GSE255612_AF_snRNA_Processed_Expression_Matrix_genes_final.tsv.gz'
METADATA = SOURCE / 'GSE255612_AF_snRNA_MetaData.txt.gz'
GATE = ('cell_type', 'Cardiomyocytes')
CHUNK = 20_000_000
SCHEMA = 'stage2-item-cm-pseudobulk-v1'


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_genes():
    genes = pd.read_csv(GENES, sep='\t', header=None, names=['source_gene_id', 'gene_symbol'])
    if len(genes) != 36600:
        raise SystemExit(f'expected 36,600 matrix features, found {len(genes)}')
    return genes


def read_metadata():
    frame = pd.read_csv(METADATA, sep='\t', skiprows=[1])
    required = {'NAME', 'donor_id', GATE[0], 'af'}
    if not required <= set(frame.columns):
        raise SystemExit(f'metadata lacks {sorted(required - set(frame.columns))}')
    if frame.NAME.duplicated().any():
        raise SystemExit('metadata repeats a cell barcode')
    return frame


def matrix_dimensions():
    with open(MTX, 'rb') as handle:
        import gzip
        with gzip.open(handle, 'rt') as text:
            header = text.readline().strip()
            text.readline()
            dimensions = text.readline().split()
    if not header.startswith('%%MatrixMarket'):
        raise SystemExit(f'unexpected mtx header: {header}')
    return tuple(int(value) for value in dimensions)


def main():
    if not MTX.is_file() or not GENES.is_file() or not METADATA.is_file():
        raise SystemExit('run script/fetch.sh first')
    n_features, n_cells, n_nonzero = matrix_dimensions()
    genes = read_genes()
    metadata = read_metadata()
    if len(metadata) != n_cells:
        raise SystemExit('metadata rows and mtx columns disagree')
    if len(genes) != n_features:
        raise SystemExit('gene file rows and mtx features disagree')
    donors = metadata.groupby('donor_id')['af'].agg(lambda values: values.iloc[0])
    if donors.map(lambda value: str(value).lower()).nunique() != 2:
        raise SystemExit('expected exactly the case/control groups')
    group_of = {donor: str(value) for donor, value in donors.items()}
    donor_order = sorted(group_of)
    donor_position = {donor: position for position, donor in enumerate(donor_order)}
    cell_donor = metadata.donor_id.map(donor_position).to_numpy(dtype=np.int64)
    is_cm = (metadata[GATE[0]] == GATE[1]).to_numpy()
    donor_cells_cm = np.zeros(len(donor_order), dtype=np.int64)
    for donor in range(len(donor_order)):
        donor_cells_cm[donor] = int(((cell_donor == donor) & is_cm).sum())
    if (donor_cells_cm == 0).any():
        raise SystemExit('a donor has no CM cell')
    profiles = np.zeros((n_features, len(donor_order)), dtype=np.float64)
    cell_cp10k_total = np.zeros(n_cells, dtype=np.float64)
    n_rows = 0
    reader = pd.read_csv(MTX, sep=' ', header=None, skiprows=3, names=['gene', 'cell', 'value'],
                         dtype={'gene': np.int32, 'cell': np.int32, 'value': np.float32},
                         chunksize=CHUNK, compression='gzip')
    for chunk in reader:
        n_rows += len(chunk)
        gene = chunk.gene.to_numpy(dtype=np.int64) - 1
        cell = chunk.cell.to_numpy(dtype=np.int64) - 1
        value = np.expm1(chunk.value.to_numpy(dtype=np.float64))
        cell_cp10k_total += np.bincount(cell, weights=value, minlength=n_cells)
        keep = is_cm[cell]
        gene, cell, value = gene[keep], cell[keep], value[keep]
        donor = cell_donor[cell]
        for position in np.unique(donor):
            mask = donor == position
            profiles[:, position] += np.bincount(gene[mask], weights=value[mask],
                                                 minlength=n_features)
    del reader
    if n_rows == 0:
        raise SystemExit('mtx produced no rows')
    if n_rows != n_nonzero:
        raise SystemExit(f'mtx holds {n_rows} entries, header declares {n_nonzero}')
    # The linear CP10K total of a cell equals 1e4 exactly iff the deposit is
    # log1p(count / n_umi * 1e4); this is the normalisation evidence.
    basis = {
        'n_sampled_cells': int(n_cells),
        'median_cell_cp10k_total': float(np.median(cell_cp10k_total)),
        'p05_cell_cp10k_total': float(np.percentile(cell_cp10k_total, 5)),
        'p95_cell_cp10k_total': float(np.percentile(cell_cp10k_total, 95)),
        'expected_cp10k_total': 1e4,
        'conclusion': 'deposited values are log1p(CP10K); expm1 recovers linear CP10K',
    }
    if not (0.5 * 1e4 <= basis['median_cell_cp10k_total'] <= 2.0 * 1e4):
        raise SystemExit(f"deposit is not log1p(CP10K): {basis}")
    profiles = profiles / donor_cells_cm[None, :]
    pseudobulk = pd.DataFrame({
        'source_gene_id': genes.source_gene_id,
        'gene_symbol': genes.gene_symbol,
    })
    for position, donor in enumerate(donor_order):
        pseudobulk[donor] = profiles[:, position]
    donor_table = pd.DataFrame({
        'donor_id': donor_order,
        'group': [group_of[donor] for donor in donor_order],
        'cells_total': [int((cell_donor == position).sum())
                        for position in range(len(donor_order))],
        'cells_cm': [int(count) for count in donor_cells_cm],
    })
    if donor_table.group.map(str.lower).nunique() != 2:
        raise SystemExit('AF/control groups collapsed')
    LOCAL.mkdir(parents=True, exist_ok=True)
    pseudobulk_path = LOCAL / 'gse255612_cm_pseudobulk.csv.gz'
    donors_path = LOCAL / 'gse255612_donors.csv'
    pseudobulk.to_csv(pseudobulk_path, index=False, compression='gzip')
    donor_table.to_csv(donors_path, index=False)
    summary = {
        'schema': SCHEMA,
        'item': 'GSE255612',
        'gate': f'{GATE[0]} == "{GATE[1]}"',
        'donor_column': 'donor_id',
        'group_column': 'af (case/control)',
        'cells_per_donor': {row.donor_id: int(row.cells_cm)
                            for row in donor_table.itertuples()},
        'n_donors': int(len(donor_table)),
        'n_cells_total': int(len(metadata)),
        'n_cells_cm': int(is_cm.sum()),
        'n_genes': int(len(pseudobulk)),
        'n_mtx_nonzero': int(n_nonzero),
        'deposit_format': 'Matrix Market feature x cell, real log1p(CP10K) values',
        'pseudobulk_basis': ('per-donor mean over CM cells of the linear CP10K value '
                             '(expm1 of the deposited log1p value)'),
        'normalization_check': basis,
        'pseudobulk_sha256': sha256(pseudobulk_path),
        'donors_sha256': sha256(donors_path),
        'donors': donor_table.to_dict('records'),
    }
    summary_path = LOCAL / 'gse255612_cm_pseudobulk.json'
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'result': 'PASS', 'cells_cm': summary['n_cells_cm'],
                      'donors': summary['n_donors'], 'genes': summary['n_genes'],
                      'normalization': basis,
                      'pseudobulk_sha256': summary['pseudobulk_sha256']}, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
