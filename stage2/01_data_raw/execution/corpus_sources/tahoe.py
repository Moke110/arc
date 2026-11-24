"""Bounded, resumable Tahoe raw-count projection at the bases release boundary."""
from __future__ import annotations

import ast
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from anndata.io import read_elem


SCHEMA = 'bases-tahoe-shared-counts-v1'


def _digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def _json(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    temporary.replace(path)


def _mapping(mapping, expected_genes):
    mapping = mapping.copy().reset_index(drop=True)
    if len(mapping) != expected_genes or mapping.gene_index.astype(int).tolist() != list(range(expected_genes)):
        raise ValueError('mapping must contain every ordered gene_index exactly once')
    for field in ('symbol', 'hgnc_id', 'native_token_id', 'native_ensembl_id'):
        if mapping[field].astype(str).eq('').any() or mapping[field].astype(str).duplicated().any():
            raise ValueError(f'missing or ambiguous Tahoe mapping: {field}')
    tokens = mapping.native_token_id.astype(int).to_numpy()
    if np.any(tokens < 3):
        raise ValueError('gene mapping includes a marker token')
    return mapping


def _plan(source, mapping, manifest, expected_shards, expected_genes):
    source, manifest = Path(source), Path(manifest)
    mapping = _mapping(mapping, expected_genes)
    declared = {}
    for line in manifest.read_text().splitlines():
        if not line.strip():
            continue
        sha, name = line.split(maxsplit=1)
        name = name.lstrip('* ')
        if name in declared and declared[name] != sha:
            raise ValueError(f'conflicting source seal: {name}')
        declared[name] = sha
    shards = sorted((source / 'data').glob('train-*.parquet'))
    if len(shards) != expected_shards:
        raise ValueError(f'expected {expected_shards} expression shards, found {len(shards)}')
    metadata = sorted((source / 'metadata').glob('*.parquet'))
    records, by_plate = {}, {}
    planned = [*metadata, *shards]
    plan_started = time.perf_counter()
    for file_number, path in enumerate(planned, 1):
        relative = str(path.relative_to(source))
        if relative not in declared:
            raise ValueError(f'source is not declared in immutable manifest: {relative}')
        actual = _digest(path)
        if actual != declared[relative]:
            raise ValueError(f'source SHA256 mismatch: {relative}')
        records[relative] = {'sha256': actual, 'bytes': path.stat().st_size}
        if path in shards:
            parquet = pq.ParquetFile(path)
            records[relative]['cells'] = parquet.metadata.num_rows
            plates = set()
            for batch in parquet.iter_batches(batch_size=4096, columns=['plate'], use_threads=False):
                plates.update(batch.column(0).to_pylist())
            records[relative]['plates'] = sorted(plates)
            for plate in sorted(plates):
                by_plate.setdefault(plate, []).append(path)
        if file_number % 100 == 0:
            print(json.dumps({'phase': 'source_hash_and_plate_scan', 'files_complete': file_number,
                              'files': len(planned), 'seconds': time.perf_counter() - plan_started}), flush=True)
    sample = pq.read_table(source / 'metadata/sample_metadata.parquet').to_pandas()
    if sample['sample'].duplicated().any():
        raise ValueError('ambiguous source well identities')
    sample_meta = {}
    for row in sample.itertuples(index=False):
        treatments = ast.literal_eval(row.drugname_drugconc)
        if len(treatments) != 1:
            raise ValueError(f'{row.sample}: expected one source treatment')
        drug, dose, unit = treatments[0]
        if unit not in ('uM', 'µM', 'μM') or not np.isfinite(float(dose)) or float(dose) < 0:
            raise ValueError(f'{row.sample}: unresolved source dose')
        sample_meta[row.sample] = {'plate': row.plate, 'drug': drug, 'dose_um': float(dose)}
    lines = pq.read_table(source / 'metadata/cell_line_metadata.parquet',
                          columns=['Cell_ID_Cellosaur', 'cell_name']).to_pandas()
    cell_names = dict(zip(lines.Cell_ID_Cellosaur, lines.cell_name))
    line_ids = sorted(cell_names)
    if not line_ids or any(not isinstance(line, str) or not line for line in line_ids):
        raise ValueError('invalid source cell-line identities')
    gene_file = source / 'metadata/gene_metadata.parquet'
    max_token = int(mapping.native_token_id.astype(int).max())
    if gene_file.exists():
        native = pq.read_table(gene_file, columns=['token_id', 'ensembl_id']).to_pandas()
        if native.token_id.duplicated().any() or native.ensembl_id.duplicated().any():
            raise ValueError('ambiguous native Tahoe gene identities')
        native = native.set_index('token_id')
        for row in mapping.itertuples(index=False):
            if int(row.native_token_id) not in native.index or native.at[int(row.native_token_id), 'ensembl_id'] != row.native_ensembl_id:
                raise ValueError('mapping disagrees with immutable native token identity')
        max_token = int(native.index.max())
    lookup = np.full(max_token + 1, -1, dtype=np.int32)
    lookup[mapping.native_token_id.astype(int)] = mapping.gene_index.astype(int)
    evidence = {'schema': SCHEMA, 'source_manifest_sha256': _digest(manifest),
                'source_files': records, 'mapping_sha256': hashlib.sha256(mapping.to_csv(index=False).encode()).hexdigest()}
    input_sha = hashlib.sha256(json.dumps(evidence, sort_keys=True).encode()).hexdigest()
    return source, mapping, by_plate, sample_meta, line_ids, cell_names, lookup, evidence, input_sha


def _plate_layout(plate, sample_meta, line_ids):
    wells = sorted(well for well, value in sample_meta.items() if value['plate'] == plate)
    if not wells:
        raise ValueError(f'unknown source plate: {plate}')
    return wells, {well: index for index, well in enumerate(wells)}, {line: index for index, line in enumerate(line_ids)}


def _entries(batch, plate, well_ix, line_ix, lookup, sample_meta):
    genes, expression = batch.column('genes'), batch.column('expressions')
    offsets = genes.offsets.to_numpy().astype(np.int64)
    other = expression.offsets.to_numpy()
    if not np.array_equal(offsets, other) or np.any(np.diff(offsets) < 1):
        raise ValueError('source gene/count lists disagree or lack marker')
    tokens = genes.values.to_numpy().astype(np.int64)
    values = expression.values.to_numpy()
    if np.any(tokens < 0) or np.any(tokens >= len(lookup)) or np.any(tokens[offsets[:-1]] >= 3):
        raise ValueError('invalid native gene or marker token')
    keep = np.ones(len(tokens), dtype=bool)
    keep[offsets[:-1]] = False
    if np.any(tokens[keep] < 3):
        raise ValueError('marker token appeared inside source gene counts')
    count_values = values[keep]
    if np.any(~np.isfinite(count_values)) or np.any(count_values < 0) or np.any(count_values != np.floor(count_values)):
        raise ValueError('source expressions must be nonnegative integer counts after dropping marker')
    samples, lines, plates = (batch.column(field).to_pylist() for field in ('sample', 'cell_line_id', 'plate'))
    rows = np.full(len(batch), -1, dtype=np.int64)
    for index, (well, line, cell_plate) in enumerate(zip(samples, lines, plates, strict=True)):
        if well not in sample_meta or sample_meta[well]['plate'] != cell_plate or line not in line_ix:
            raise ValueError('source cell has unresolved well/plate/cell-line identity')
        if cell_plate == plate:
            rows[index] = well_ix[well] * len(line_ix) + line_ix[line]
    gene_indices = lookup[tokens[keep]]
    repeated_rows = np.repeat(rows, np.diff(offsets) - 1)
    values = values[keep].astype(np.int64)
    selected = (repeated_rows >= 0) & (gene_indices >= 0)
    return repeated_rows[selected], gene_indices[selected], values[selected], rows[rows >= 0]


def _obs(plate, selected, wells, cells, line_ids, sample_meta, cell_names):
    rows = []
    for index, cell_count in zip(selected, cells, strict=True):
        well, line = wells[index // len(line_ids)], line_ids[index % len(line_ids)]
        treatment = sample_meta[well]
        rows.append({'profile_id': f'{plate}:{well}:{line}', 'plate': plate, 'well': well,
                     'cell_line_id': line, 'cell_name': cell_names[line], 'drug': treatment['drug'],
                     'dose_um': treatment['dose_um'], 'is_control': treatment['drug'] == 'DMSO_TF',
                     'control_group': f'tahoe::{plate}', 'n_cells': int(cell_count)})
    return pd.DataFrame(rows).set_index('profile_id')


def _plate_key(plate):
    return int(plate[5:]) if plate.startswith('plate') and plate[5:].isdigit() else plate
