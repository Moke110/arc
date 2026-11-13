"""Bounded, resumable Tahoe raw-count projection at the bases release boundary."""
from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
import time

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from anndata import AnnData
from anndata.io import read_elem
from scipy.sparse import coo_matrix, csr_matrix


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


def _aggregate(plate, shards, sample_meta, line_ids, lookup, batch_size, independent=False):
    wells, well_ix, line_ix = _plate_layout(plate, sample_meta, line_ids)
    counts = np.zeros((len(wells) * len(line_ids), int(lookup.max()) + 1), dtype=np.int64)
    cells = np.zeros(len(counts), dtype=np.int64)
    for shard_index, path in enumerate(shards):
        for batch in pq.ParquetFile(path).iter_batches(batch_size=batch_size,
                columns=['genes', 'expressions', 'sample', 'cell_line_id', 'plate'], use_threads=False):
            rows, genes, values, cell_rows = _entries(batch, plate, well_ix, line_ix, lookup, sample_meta)
            np.add.at(cells, cell_rows, 1)
            if independent:
                # A separate sparse reduction groups all entries before adding;
                # it does not reuse production's elementwise accumulation.
                contribution = coo_matrix((values, (rows, genes)), shape=counts.shape).tocsr()
                for row in np.flatnonzero(np.diff(contribution.indptr)):
                    lo, hi = contribution.indptr[row:row + 2]
                    counts[row, contribution.indices[lo:hi]] += contribution.data[lo:hi]
            else:
                np.add.at(counts, (rows, genes), values)
        if (shard_index + 1) % 40 == 0:
            print(json.dumps({'plate': plate, 'shards_complete': shard_index + 1, 'shards': len(shards),
                              'independent': independent}), flush=True)
    selected = np.flatnonzero(cells)
    if np.any(counts < 0) or counts.max(initial=0) > np.iinfo(np.int32).max:
        raise ValueError('aggregated counts overflow int32')
    return csr_matrix(counts[selected].astype(np.int32)), cells[selected], selected, wells


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


def _write_h5ad(destination, parts, obs, mapping):
    var = mapping[['symbol', 'hgnc_id']].copy().set_index('symbol')
    var['ensembl_id'] = mapping.native_ensembl_id.to_numpy()
    var['token_id'] = mapping.native_token_id.astype(str).to_numpy()
    var['gene_index'] = mapping.gene_index.astype(int).to_numpy()
    temporary = destination.with_suffix('.h5ad.tmp')
    AnnData(X=csr_matrix((len(obs), len(var)), dtype=np.int32), obs=obs, var=var).write_h5ad(temporary)
    nnz = sum(record['nnz'] for _, record in parts)
    with h5py.File(temporary, 'a') as file:
        del file['X']
        group = file.create_group('X')
        group.attrs['encoding-type'], group.attrs['encoding-version'] = 'csr_matrix', '0.1.0'
        group.attrs['shape'] = np.array([len(obs), len(var)], dtype=np.int64)
        options = {'chunks': (min(max(nnz, 1), 262144),), 'compression': 'gzip', 'compression_opts': 1}
        data = group.create_dataset('data', (nnz,), dtype=np.int32, **options)
        indices = group.create_dataset('indices', (nnz,), dtype=np.int32, **options)
        indptr = group.create_dataset('indptr', (len(obs) + 1,), dtype=np.int64)
        row_offset = nz_offset = 0
        for part, record in parts:
            with np.load(part, allow_pickle=False) as stored:
                size = record['nnz']
                data[nz_offset:nz_offset + size] = stored['data']
                indices[nz_offset:nz_offset + size] = stored['indices']
                indptr[row_offset:row_offset + record['profiles']] = stored['indptr'][:-1] + nz_offset
                row_offset += record['profiles']
                nz_offset += size
        indptr[-1] = nnz
        file.attrs['bases_counts_schema'] = SCHEMA
        file.attrs['preprocessing'] = 'summed raw UMI counts projected before normalization; no normalization stored'
    temporary.replace(destination)


def build_tahoe(source, mapping, source_manifest, output, *, expected_shards=3388, expected_genes=18164, batch_size=1024, manifest='manifest.json', checkpoints=None):
    """Publish shared-axis raw counts; completed, sealed plates resume safely."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    source, mapping, plates, sample_meta, line_ids, cell_names, lookup, evidence, input_sha = _plan(
        source, mapping, source_manifest, expected_shards, expected_genes)
    checkpoints = Path(checkpoints) if checkpoints is not None else output / '.checkpoints'
    checkpoints.mkdir(parents=True, exist_ok=True)
    parts, observations = [], []
    resumed = 0
    for plate in sorted(plates, key=_plate_key):
        part, seal = checkpoints / f'{plate}.npz', checkpoints / f'{plate}.json'
        if seal.exists():
            record = json.loads(seal.read_text())
            if record['input_sha256'] != input_sha or not part.exists() or record['sha256'] != _digest(part):
                raise ValueError(f'{plate}: checkpoint source or payload seal mismatch')
            with np.load(part, allow_pickle=False) as stored:
                cells, selected = stored['cells'], stored['selected']
            wells, _, _ = _plate_layout(plate, sample_meta, line_ids)
            resumed += 1
        else:
            started = time.perf_counter()
            matrix, cells, selected, wells = _aggregate(plate, plates[plate], sample_meta, line_ids, lookup, batch_size)
            temporary = part.with_suffix('.npz.tmp')
            with temporary.open('wb') as handle:
                np.savez_compressed(handle, data=matrix.data, indices=matrix.indices, indptr=matrix.indptr,
                                    cells=cells, selected=selected)
            temporary.replace(part)
            record = {'input_sha256': input_sha, 'sha256': _digest(part), 'profiles': len(cells),
                      'nnz': matrix.nnz, 'cells': int(cells.sum()), 'seconds': time.perf_counter() - started,
                      'shards': [path.name for path in plates[plate]]}
            _json(seal, record)
            print(json.dumps({'plate_complete': plate, **{key: record[key] for key in ('profiles', 'cells', 'seconds')}}), flush=True)
            del matrix
        parts.append((part, record))
        observations.append(_obs(plate, selected, wells, cells, line_ids, sample_meta, cell_names))
    obs = pd.concat(observations)
    if obs.index.has_duplicates:
        raise ValueError('duplicate published profile identities')
    source_cells = sum(record.get('cells', 0) for record in evidence['source_files'].values())
    if int(obs.n_cells.sum()) != source_cells:
        raise ValueError('source cells lost or duplicated during plate aggregation')
    if expected_shards == 3388 and (len(obs) != 67018 or int(obs.is_control.sum()) != 1400):
        raise ValueError(f'unexpected full Tahoe profile/control population: {len(obs)}/{int(obs.is_control.sum())}')
    mapping.to_csv(output / 'token_mapping.csv', index=False)
    _write_h5ad(output / 'tahoe.h5ad', parts, obs, mapping)
    result = {'schema': SCHEMA, 'input_sha256': input_sha, 'inputs': evidence, 'profiles': len(obs),
              'cells': int(obs.n_cells.sum()), 'controls': int(obs.is_control.sum()),
              'treated': int((~obs.is_control).sum()), 'genes': len(mapping), 'plates': len(parts), 'resumed_plates': resumed,
              'source_shards': expected_shards, 'profile_order': 'plate numeric, source well lexical, cell_line_id lexical',
              'h5ad_sha256': _digest(output / 'tahoe.h5ad'), 'mapping_sha256': _digest(output / 'token_mapping.csv')}
    _json(output / manifest, result)
    return result


def verify_tahoe(source, mapping, source_manifest, output, *, expected_shards=3388, expected_genes=18164, batch_size=1024, manifest='manifest.json'):
    """Re-read every source cell and compare every published profile with sparse reduction."""
    output = Path(output)
    release = json.loads((Path(output) / manifest).read_text())
    source, mapping, plates, sample_meta, line_ids, cell_names, lookup, _, input_sha = _plan(
        source, mapping, source_manifest, expected_shards, expected_genes)
    if release['input_sha256'] != input_sha or release['h5ad_sha256'] != _digest(output / 'tahoe.h5ad'):
        raise ValueError('release source or raw-count payload seal mismatch')
    if release['mapping_sha256'] != _digest(output / 'token_mapping.csv'):
        raise ValueError('release token mapping seal mismatch')
    verified_profiles = verified_cells = 0
    with h5py.File(output / 'tahoe.h5ad') as file:
        actual_obs, actual_var = read_elem(file['obs']), read_elem(file['var'])
        if actual_var.hgnc_id.tolist() != mapping.hgnc_id.tolist() or actual_var.token_id.astype(str).tolist() != mapping.native_token_id.astype(str).tolist():
            raise ValueError('published gene axis differs from source mapping')
        matrix = file['X']
        for plate in sorted(plates, key=_plate_key):
            expected, cells, selected, wells = _aggregate(plate, plates[plate], sample_meta, line_ids, lookup, batch_size, independent=True)
            obs = _obs(plate, selected, wells, cells, line_ids, sample_meta, cell_names)
            if not actual_obs.iloc[verified_profiles:verified_profiles + len(obs)].equals(obs):
                # H5AD may encode strings categorically; compare semantic values.
                left = actual_obs.iloc[verified_profiles:verified_profiles + len(obs)].astype(str)
                if not left.equals(obs.astype(str)):
                    raise ValueError(f'{plate}: profile identities or source metadata differ')
            pointers = matrix['indptr'][verified_profiles:verified_profiles + len(obs) + 1]
            actual = csr_matrix((matrix['data'][pointers[0]:pointers[-1]], matrix['indices'][pointers[0]:pointers[-1]], pointers - pointers[0]), shape=expected.shape)
            if (actual != expected).nnz:
                raise ValueError(f'{plate}: raw-count source replay mismatch')
            verified_profiles += len(obs)
            verified_cells += int(cells.sum())
            print(json.dumps({'verified_plate': plate, 'verified_profiles': verified_profiles, 'verified_cells': verified_cells}), flush=True)
        if verified_profiles != len(actual_obs) or verified_cells != release['cells']:
            raise ValueError('source population differs from release')
    return {'schema': SCHEMA, 'verified_profiles': verified_profiles, 'verified_cells': verified_cells,
            'input_sha256': input_sha, 'h5ad_sha256': release['h5ad_sha256']}
