"""Native Tahoe bulk wrappers over the frozen shard-reading and reduction kernels."""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from anndata import AnnData
from anndata.io import read_elem
from scipy.sparse import csr_matrix

from bases import tahoe as projected
from bases.tahoe import _aggregate, _digest, _json, _obs, _plate_key, _plate_layout

SCHEMA = projected.SCHEMA


def _plan(source, mapping, manifest, expected_shards, expected_genes):
    source = Path(source)
    gene_file = source / "metadata/gene_metadata.parquet"
    tokens = pq.read_table(gene_file, columns=["token_id"]).column(0).to_pylist()
    if set(mapping.native_token_id) != set(tokens):
        raise ValueError("native Tahoe map does not cover every deposited gene token")
    # The retained projected planner requires unique annotation columns even
    # though the shard reduction keys exclusively by token. Feed Ensembl IDs
    # for those checks and write the untouched publisher metadata separately.
    planning_map = mapping.assign(symbol=mapping.native_ensembl_id,
                                  hgnc_id=mapping.native_ensembl_id)
    source, _, plates, samples, lines, names, lookup, evidence, _ = projected._plan(
        source, planning_map, manifest, expected_shards, expected_genes
    )
    evidence["mapping_sha256"] = hashlib.sha256(mapping.to_csv(index=False).encode()).hexdigest()
    input_sha = hashlib.sha256(json.dumps(evidence, sort_keys=True).encode()).hexdigest()
    return source, mapping, plates, samples, lines, names, lookup, evidence, input_sha


def _write_h5ad(destination, parts, obs, mapping):
    var = mapping[['symbol']].copy()
    var.index = pd.Index(mapping.native_ensembl_id, name='ensembl_id')
    var['ensembl_id'] = mapping.native_ensembl_id.to_numpy()
    var['token_id'] = mapping.native_token_id.astype(str).to_numpy()
    var['native_column'] = mapping.gene_index.astype(int).to_numpy()
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
        file.attrs['preprocessing'] = 'summed native raw UMI counts; no normalization stored'
    temporary.replace(destination)


def build_tahoe(source, mapping, source_manifest, output, *, expected_shards=3388, expected_genes=18164, batch_size=1024, manifest='manifest.json', checkpoints=None):
    """Publish raw counts; completed, sealed plates resume safely."""
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
        if (actual_var.ensembl_id.tolist() != mapping.native_ensembl_id.tolist()
                or actual_var.token_id.astype(str).tolist() != mapping.native_token_id.astype(str).tolist()):
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
