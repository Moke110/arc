#!/usr/bin/env python
"""GSE270788 — author CM gate and per-donor CM pseudobulk (ADR 0066 D8).

The series is a human left-ventricle 10x Multiome deposit whose GEX sub-series
ships one filtered_feature_bc_matrix.h5 per donor (23 GEX donors). GSE218392 is
its SuperSeries, so only GEX payloads are read; every h5 in this deposit also
carries the ATAC peak features, which are dropped here (feature_type ==
"Gene Expression" only).

Gate: ``GSE270788_metadata.csv.gz`` column ``cell.type.l1 == "Cardiomyocyte"``.
The metadata barcode is `<lane-prefix>_<10x barcode>`; the h5 barcode is the raw
10x barcode, so the lane prefix is stripped for the join. Each donor's CM counts
are summed into one pseudobulk profile before any testing. No clustering, no
normalisation, no dense cell x gene matrix: the h5 matrix is read as a CSC over
cells and only the selected CM columns are accumulated.

Writes under the item's declared ``local/`` area (never ``source/``):

    local/pseudobulk/gse270788_cm_pseudobulk.csv.gz   genes x 23 donors, integer
    local/pseudobulk/gse270788_donors.csv             per-donor cells/arm table
    local/pseudobulk/gse270788_cm_pseudobulk.json     tracked summary (digests)

Usage (arc python):
    /home/chang/miniconda3/envs/arc/bin/python script/build_pseudobulk.py
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ITEM = Path(__file__).resolve().parents[1]
SOURCE = ITEM / 'source'
LOCAL = ITEM / 'local' / 'pseudobulk'
GATE = ('cell.type.l1', 'Cardiomyocyte')
ARM_COLUMN = 'HF.etiology'
SAMPLE_COLUMN = 'sample'
SCHEMA = 'stage2-item-cm-pseudobulk-v1'


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def sample_order(name):
    match = re.search(r'MA(\d+)', name)
    return int(match.group(1)) if match else 10 ** 6


def read_metadata():
    path = SOURCE / 'GSE270788_metadata.csv.gz'
    frame = pd.read_csv(path)
    required = {'barcode', SAMPLE_COLUMN, ARM_COLUMN, GATE[0]}
    if not required <= set(frame.columns):
        raise SystemExit(f'metadata lacks {sorted(required - set(frame.columns))}')
    frame['lane_prefix'] = frame.barcode.str.split('_', n=1).str[0]
    frame['raw_barcode'] = frame.barcode.str.split('_', n=1).str[1]
    if frame.raw_barcode.isna().any():
        raise SystemExit('metadata barcode lacks the lane prefix')
    return frame


def feature_table(h5_path):
    with h5py.File(h5_path, 'r') as handle:
        matrix = handle['matrix']
        feature_type = matrix['features/feature_type'][:]
        ids = np.array([value.decode() for value in matrix['features/id'][:]])
        names = np.array([value.decode() for value in matrix['features/name'][:]])
    mask = feature_type == b'Gene Expression'
    return ids[mask], names[mask], int(mask.sum()), int(len(mask))


def count_cm_cells(h5_path, barcodes, wanted):
    """Return the CM column indices of one h5 and the per-cell CM barcode match."""
    barcodes = np.array([value.decode() for value in barcodes])
    # np.isin on a Python set does not compare element-wise; pass a sorted array.
    wanted_array = np.asarray(sorted(wanted), dtype=barcodes.dtype)
    return np.flatnonzero(np.isin(barcodes, wanted_array)), barcodes


def pseudobulk_one_file(h5_path, wanted, n_features):
    with h5py.File(h5_path, 'r') as handle:
        matrix = handle['matrix']
        indptr = matrix['indptr'][:]
        indices = matrix['indices'][:]
        data = matrix['data'][:]
        barcodes = matrix['barcodes'][:]
        feature_type = matrix['features/feature_type'][:]
    # The multiome h5 carries both Gene Expression and Peaks feature rows; map
    # every row to its Gene Expression position (-1 for peaks) before summing.
    mask = feature_type == b'Gene Expression'
    position = np.full(len(mask), -1, dtype=np.int64)
    position[mask] = np.arange(int(mask.sum()), dtype=np.int64)
    columns, decoded = count_cm_cells(h5_path, barcodes, wanted)
    profile = np.zeros(n_features, dtype=np.int64)
    for column in columns:
        start, end = int(indptr[column]), int(indptr[column + 1])
        if start == end:
            continue
        rows = position[indices[start:end]]
        keep = rows >= 0
        if not keep.any():
            continue
        profile += np.bincount(rows[keep], weights=data[start:end][keep],
                               minlength=n_features).astype(np.int64)
    return profile, columns, decoded


def infer_sample(h5_path, barcodes, raw_to_samples):
    """Resolve the donor of an h5 file: filename token first, barcode vote second."""
    match = re.search(r'_(MA\d+)_', h5_path.name)
    if match:
        return match.group(1), 'filename'
    votes = {}
    for value in barcodes:
        token = value.decode() if isinstance(value, bytes) else str(value)
        for sample in raw_to_samples.get(token, ()):
            votes[sample] = votes.get(sample, 0) + 1
    if not votes:
        raise SystemExit(f'cannot resolve the sample of {h5_path.name}')
    ranked = sorted(votes.values(), reverse=True)
    sample = max(votes, key=votes.get)
    if len(ranked) > 1 and ranked[0] <= 1.5 * ranked[1]:
        raise SystemExit(f'ambiguous barcode vote for {h5_path.name}: {votes}')
    if ranked[0] < 0.25 * len(barcodes):
        raise SystemExit(f'weak barcode vote for {h5_path.name}: {votes}')
    return sample, 'barcode_intersection'


def main():
    metadata = read_metadata()
    # The same 10x barcode can recur across lanes, so the join is per sample: a
    # raw barcode maps to a *set* of samples and the h5 is matched to its own
    # sample's CM barcode set. The only file without a sample token in its name
    # is resolved by the barcode vote that scores every candidate sample.
    raw_to_samples = {}
    for sample, raw in zip(metadata[SAMPLE_COLUMN], metadata.raw_barcode):
        raw_to_samples.setdefault(raw, set()).add(sample)
    cm = metadata[metadata[GATE[0]] == GATE[1]]
    if cm.empty:
        raise SystemExit('CM gate selects no cell')
    sample_arms = (metadata.groupby(SAMPLE_COLUMN)[ARM_COLUMN].agg(lambda values: values.iloc[0])
                   if metadata.groupby(SAMPLE_COLUMN)[ARM_COLUMN].nunique().max() == 1
                   else None)
    if sample_arms is None:
        raise SystemExit('a GSE270788 sample carries more than one HF.etiology arm')
    h5_files = sorted(SOURCE.glob('*_filtered_feature_bc_matrix.h5'))
    if len(h5_files) != 23:
        raise SystemExit(f'expected 23 GEX h5 files, found {len(h5_files)}')
    reference_ids = reference_names = None
    profiles, donors = {}, []
    for path in h5_files:
        ids, names, n_gex, n_total = feature_table(path)
        if reference_ids is None:
            reference_ids, reference_names = ids, names
        elif not (np.array_equal(ids, reference_ids) and np.array_equal(names, reference_names)):
            raise SystemExit(f'GEX feature order differs in {path.name}')
        with h5py.File(path, 'r') as handle:
            barcodes = [value.decode() for value in handle['matrix/barcodes'][:]]
        sample, basis = infer_sample(path, barcodes, raw_to_samples)
        wanted = set(cm.loc[cm[SAMPLE_COLUMN] == sample, 'raw_barcode'])
        profile, columns, decoded = pseudobulk_one_file(path, wanted, n_gex)
        matched = {decoded[index] for index in columns}
        donors.append({
            'sample': sample,
            'donor_id': sample,
            'arm': sample_arms.loc[sample],
            'cells_total': int(len(barcodes)),
            'cells_cm': int(len(columns)),
            'cells_cm_metadata': int(len(wanted)),
            'cm_barcodes_matched': int(len(matched & wanted)),
            'source_file': path.name,
            'sample_resolution': basis,
            'gex_features': n_gex,
            'all_features': n_total,
        })
        profiles[sample] = profile
    donors = pd.DataFrame(donors).sort_values(SAMPLE_COLUMN, key=lambda s: s.map(sample_order))
    donors = donors.reset_index(drop=True)
    if donors.cells_cm.duplicated().any() or len(donors) != 23:
        raise SystemExit(f'expected one pseudobulk column per GEX donor; '
                         f'n={len(donors)} dup={donors.cells_cm.duplicated().sum()}\n'
                         + donors.to_string())
    missing = set(cm[SAMPLE_COLUMN].unique()) - set(donors[SAMPLE_COLUMN])
    if missing:
        raise SystemExit(f'CM-bearing samples without a matrix: {sorted(missing)}')
    pseudobulk = pd.DataFrame({
        'source_gene_id': reference_ids,
        'gene_symbol': reference_names,
    })
    for sample in donors[SAMPLE_COLUMN]:
        pseudobulk[sample] = profiles[sample]
    if (pseudobulk[donors[SAMPLE_COLUMN].tolist()] < 0).any().any():
        raise SystemExit('negative pseudobulk count')
    LOCAL.mkdir(parents=True, exist_ok=True)
    pseudobulk_path = LOCAL / 'gse270788_cm_pseudobulk.csv.gz'
    donors_path = LOCAL / 'gse270788_donors.csv'
    pseudobulk.to_csv(pseudobulk_path, index=False, compression='gzip')
    donors.to_csv(donors_path, index=False)
    summary = {
        'schema': SCHEMA,
        'item': 'GSE270788',
        'gate': f'{GATE[0]} == "{GATE[1]}"',
        'arm_column': ARM_COLUMN,
        'cells_per_donor': {row.sample: int(row.cells_cm) for row in donors.itertuples()},
        'n_donors': int(len(donors)),
        'n_genes_gex': int(len(pseudobulk)),
        'n_cells_total': int(donors.cells_total.sum()),
        'n_cells_cm': int(donors.cells_cm.sum()),
        'barcode_join': 'metadata barcode = lane-prefix + "_" + h5 raw barcode',
        'sample_resolution': {row.sample: row.sample_resolution for row in donors.itertuples()},
        'pseudobulk_dtype': 'integer counts summed per donor before testing',
        'pseudobulk_sha256': sha256(pseudobulk_path),
        'donors_sha256': sha256(donors_path),
        'donors': donors.to_dict('records'),
    }
    summary_path = LOCAL / 'gse270788_cm_pseudobulk.json'
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'result': 'PASS', 'cells_cm': summary['n_cells_cm'],
                      'donors': summary['n_donors'],
                      'genes': summary['n_genes_gex'],
                      'pseudobulk_sha256': summary['pseudobulk_sha256']}, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
