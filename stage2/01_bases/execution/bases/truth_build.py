"""Checkpointed raw-count DESeq2 release and independent source/statistic replay."""
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import tempfile

import h5py
from anndata.io import read_elem
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from .significance import (FIELDS, SCHEMA, bh_adjust, checked_file, contrast_plan,
                           group_key, membership, publish_statistics, ranked_indices,
                           sha256, write_json)

WORKER = Path(__file__).with_name('deseq2_contrasts.R')


def source_plan(training, native, *, fixture=False):
    """Verify independent training/raw identities and create axis-only raw contrasts."""
    from .corpora import parent_structure, training_header
    training, native = Path(training), Path(native)
    manifest = json.loads((training / 'manifest.json').read_text())
    if (manifest.get('schema') != 'stage2-indexed-training-corpora-v1'
            or bool(manifest.get('fixture')) != fixture):
        raise ValueError('invalid training release for truth building')
    checked_file(training, 'tenv.h5ad', manifest['files']['tenv.h5ad'])
    treated, genes = training_header(training / 'tenv.h5ad')
    if (list(genes.gene_index) != list(range(len(genes)))
            or (not fixture and len(genes) != 18164) or genes.hgnc_id.duplicated().any()):
        raise ValueError('invalid frozen training gene axis')
    with h5py.File(native) as file:
        obs, var = read_elem(file['obs']), read_elem(file['var'])
    if var.hgnc_id.duplicated().any():
        # Native off-axis empty identifiers may repeat, but a measured HGNC may not.
        nonempty = var.hgnc_id.dropna().loc[lambda x: x.astype(str).ne('')]
        if nonempty.duplicated().any():
            raise ValueError('ambiguous native HGNC identity')
    columns = pd.Index(var.hgnc_id).get_indexer(genes.hgnc_id)
    if np.any(columns < 0) or len(set(columns)) != len(genes):
        raise ValueError('missing or ambiguous native axis mapping')
    for row in treated.itertuples(index=False):
        raw = obs.iloc[int(row.source_row)]
        parent, _, _ = parent_structure(raw.smiles)
        if parent != row.drug_parent:
            raise ValueError('drug parent differs from independently normalized native structure')
    groups, samples = contrast_plan(treated, obs)
    return genes.reset_index(drop=True), groups, samples, columns


def raw_axis_rows(matrix, rows, columns):
    """Read one contrast in bounded row chunks; never normalize DESeq2 inputs."""
    rows, columns = np.asarray(rows, dtype=int), np.asarray(columns, dtype=int)
    if (len(rows) == 0 or np.any(np.diff(rows) <= 0) or np.any(rows < 0)):
        raise ValueError('invalid raw contrast rows')
    if isinstance(matrix, h5py.Group):
        if matrix.attrs['encoding-type'] != 'csr_matrix':
            raise ValueError('expected native CSR counts')
        shape = tuple(matrix.attrs['shape'])
    else:
        shape = matrix.shape
    if rows[-1] >= shape[0] or np.any(columns < 0) or np.any(columns >= shape[1]):
        raise ValueError('raw contrast outside native matrix')
    for start in range(0, len(rows), 32):
        selected = rows[start:start + 32]
        if isinstance(matrix, h5py.Group):
            lo, hi = int(selected[0]), int(selected[-1]) + 1
            ptr = matrix['indptr'][lo:hi + 1]
            a, b = int(ptr[0]), int(ptr[-1])
            counts = csr_matrix((matrix['data'][a:b], matrix['indices'][a:b], ptr - a),
                                shape=(hi - lo, shape[1]))[selected - lo][:, columns].toarray()
        else:
            counts = matrix[selected][:, columns]
        if (counts.dtype.kind not in 'iu' or np.any(counts < 0)
                or np.any(counts > np.iinfo(np.int32).max)):
            raise ValueError('invalid DESeq2 raw integer counts')
        yield counts.astype('<i4', copy=False)


def write_part_inputs(directory, native, groups, samples, columns):
    """Seal exact raw count bytes per contrast before the R worker reads them."""
    import hashlib
    directory.mkdir(parents=True, exist_ok=True)
    groups = groups.copy()
    ordered = []
    digests = []
    with h5py.File(native) as file, (directory / 'counts.bin').open('wb') as stream:
        for row in groups.itertuples(index=False):
            part = samples[samples.group_id.eq(row.group_id)].sort_values('source_row')
            digest = hashlib.sha256()
            for counts in raw_axis_rows(file['X'], part.source_row.to_numpy(), columns):
                data = counts.tobytes(order='C')
                digest.update(data)
                stream.write(data)
            digests.append(digest.hexdigest())
            ordered.append(part)
    groups['counts_sha256'] = digests
    groups.to_csv(directory / 'groups.csv', index=False)
    pd.concat(ordered, ignore_index=True).to_csv(directory / 'samples.csv', index=False)
    return groups


def run_worker(directory, rscript, gene_count):
    with (directory / 'worker.log').open('w') as log:
        subprocess.run([str(rscript), '--vanilla', str(WORKER), str(directory), str(gene_count)],
                       stdout=log, stderr=subprocess.STDOUT, check=True)


def build(training, native, output, rscript, *, batch_size=25, workers=2, fixture=False):
    """Build all contrasts; resume only exact input/code-pinned completed parts."""
    if workers not in (1, 2) or batch_size < 1:
        raise ValueError('expected 1 or 2 bounded workers and positive checkpoint batch size')
    training, native, output = Path(training), Path(native), Path(output)
    if output.exists():
        raise FileExistsError(output)
    genes, groups, samples, columns = source_plan(training, native, fixture=fixture)
    provenance = {
        'training_manifest_sha256': sha256(training / 'manifest.json'),
        'tenv_training_sha256': sha256(training / 'tenv.h5ad'),
        'native_counts_sha256': sha256(native), 'axis_hgnc_sha256': membership(genes.hgnc_id.tolist()),
        'worker_sha256': sha256(WORKER), 'builder_sha256': sha256(Path(__file__)),
        'truth_interface_sha256': sha256(Path(__file__).with_name('significance.py')),
        'group_memberships_sha256': membership(samples[samples.condition.eq('treated')].sample_id.tolist()),
        'gene_count': len(genes), 'group_count': len(groups), 'batch_size': batch_size,
        'replication_rule': 'at least 2 treated and 2 matched-control wells; residual design degrees of freedom',
    }
    checkpoint = output.with_name(output.name + '.building')
    checkpoint.mkdir(parents=True, exist_ok=True)
    identity = checkpoint / 'identity.json'
    if identity.exists() and json.loads(identity.read_text()) != provenance:
        raise ValueError('checkpoint input/code identity changed')
    write_json(identity, provenance)
    groups.to_csv(checkpoint / 'planned_groups.csv', index=False)
    samples.to_csv(checkpoint / 'planned_samples.csv', index=False)
    genes.to_csv(checkpoint / 'genes.csv', index=False)

    def process(start):
        part = checkpoint / 'parts' / f'{start // batch_size:05d}'
        seal = part / 'manifest.json'
        expected = groups.iloc[start:start + batch_size]
        if not seal.exists():
            metadata = write_part_inputs(part, native, expected, samples, columns)
            run_worker(part, rscript, len(genes))
            result = pd.read_csv(part / 'results.csv', keep_default_na=False)
            if len(result) != len(metadata):
                raise ValueError('DESeq2 omitted a contrast')
            pd.concat([metadata.reset_index(drop=True), result], axis=1).to_csv(part / 'metadata.csv', index=False)
            (part / 'counts.bin').unlink()
            write_json(seal, {'groups': metadata.group_id.tolist(), 'files': {
                p.name: sha256(p) for p in sorted(part.iterdir()) if p.is_file() and p.name != 'manifest.json'}})
        pin = json.loads(seal.read_text())
        if pin['groups'] != expected.group_id.tolist():
            raise ValueError('checkpoint group order differs')
        for name, digest in pin['files'].items():
            checked_file(part, name, digest)
        print(f'Completed contrasts {start + len(expected)}/{len(groups)}', flush=True)
        return part, pd.read_csv(part / 'metadata.csv', keep_default_na=False)

    with ThreadPoolExecutor(max_workers=workers) as executor:
        parts = list(executor.map(process, range(0, len(groups), batch_size)))
    versions = {(path / 'versions.txt').read_text() for path, _ in parts}
    if len(versions) != 1:
        raise ValueError('DESeq2 runtime versions differ')
    provenance['r_version'], provenance['deseq2_version'] = versions.pop().splitlines()
    statistics = {}
    for field in FIELDS:
        target = checkpoint / f'{field}.npy'
        result = np.lib.format.open_memmap(target, mode='w+', dtype='<f8', shape=(len(groups), len(genes)))
        offset = 0
        for directory, part_groups in parts:
            binary = directory / f'{field}.bin'
            if binary.stat().st_size != len(part_groups) * len(genes) * 8:
                raise ValueError('DESeq2 returned invalid statistic size')
            with binary.open('rb') as stream:
                for _ in range(len(part_groups)):
                    result[offset] = np.fromfile(stream, dtype='<f8', count=len(genes))
                    offset += 1
        result.flush()
        statistics[field] = result
    metadata = pd.concat([part for _, part in parts], ignore_index=True)
    stage = checkpoint / 'release'
    if stage.exists():
        import shutil
        shutil.rmtree(stage)
    publish_statistics(stage, genes, metadata, samples, statistics, provenance, fixture=fixture)
    stage.rename(output)
    return json.loads((output / 'manifest.json').read_text())


def verify(training, native, output, rscript, *, replay_groups=3, fixture=False):
    """Replay every raw membership/BH/ranking plus independently refit representatives."""
    import hashlib
    output, native, training = Path(output), Path(native), Path(training)
    manifest = json.loads((output / 'manifest.json').read_text())
    if (manifest.get('schema') != SCHEMA or manifest.get('dataset') != 'tenv'
            or manifest.get('fixture') != fixture or not manifest.get('complete')):
        raise ValueError('invalid truth release schema')
    for name, digest in manifest['files'].items():
        checked_file(output, name, digest)
    provenance = manifest['provenance']
    for field, path in [('native_counts_sha256', native),
                        ('training_manifest_sha256', training / 'manifest.json'),
                        ('tenv_training_sha256', training / 'tenv.h5ad'),
                        ('worker_sha256', WORKER)]:
        if provenance[field] != sha256(path):
            raise ValueError(f'changed truth input evidence: {field}')
    genes, planned, samples, columns = source_plan(training, native, fixture=fixture)
    stored_genes = pd.read_csv(output / 'genes.csv', keep_default_na=False, dtype=str)
    stored_genes['gene_index'] = stored_genes.gene_index.astype(int)
    pd.testing.assert_frame_equal(genes, stored_genes, check_dtype=False)
    groups = pd.read_csv(output / 'groups.csv', keep_default_na=False)
    pd.testing.assert_frame_equal(planned, groups[planned.columns], check_dtype=False)
    pd.testing.assert_frame_equal(samples, pd.read_csv(output / 'contrast_samples.csv', keep_default_na=False),
                                  check_dtype=False)
    statistics = {field: np.load(output / f'{field}.npy', mmap_mode='r') for field in FIELDS}
    if any(value.shape != (len(groups), len(genes)) for value in statistics.values()):
        raise ValueError('invalid complete truth statistic shape')
    with h5py.File(native) as file:
        for i, row in enumerate(groups.itertuples(index=False)):
            part = samples[samples.group_id.eq(row.group_id)].sort_values('source_row')
            digest = hashlib.sha256()
            for values in raw_axis_rows(file['X'], part.source_row.to_numpy(), columns):
                digest.update(values.tobytes(order='C'))
            if digest.hexdigest() != row.counts_sha256:
                raise ValueError('raw full-axis contrast counts differ')
            expected_padj = bh_adjust(statistics['pvalue'][i])
            if not np.allclose(statistics['padj'][i], expected_padj, rtol=1e-12, atol=1e-15, equal_nan=True):
                raise ValueError('complete-release BH replay differs')
            value = json.loads((output / 'groups' / f'{row.group_id}.json').read_text())
            up, down = ranked_indices(expected_padj, statistics['log2_fold_change'][i])
            if (up != value['up_indices'] or down != value['down_indices']
                    or value['status'] != row.status or value['reason'] != row.reason
                    or bool(np.isfinite(statistics['pvalue'][i]).any()) != (row.status == 'eligible')
                    or int(row.n_significant) != len(up) + len(down)):
                raise ValueError('complete-release ranked truth/status replay differs')
    replayed = []
    if replay_groups:
        candidates = [*groups.index[groups.n_plates.gt(1)],
                      *groups.index[groups.status.eq('eligible') & groups.n_significant.eq(0)],
                      *groups.index[groups.status.eq('eligible')]]
        chosen = list(dict.fromkeys(candidates))[:replay_groups]
        with tempfile.TemporaryDirectory(prefix='.truth-refit-', dir=output.parent) as directory:
            directory = Path(directory)
            selected = groups.loc[chosen, planned.columns].reset_index(drop=True)
            write_part_inputs(directory, native, selected, samples, columns)
            run_worker(directory, rscript, len(genes))
            results = pd.read_csv(directory / 'results.csv', keep_default_na=False)
            for local, original in enumerate(chosen):
                for field in FIELDS:
                    values = np.memmap(directory / f'{field}.bin', dtype='<f8', mode='r',
                                       shape=(len(chosen), len(genes)))[local]
                    if not np.allclose(values, statistics[field][original], rtol=1e-12,
                                       atol=1e-15, equal_nan=True):
                        raise ValueError(f'independent DESeq2 refit differs: {field}')
                if results.iloc[local].status != groups.iloc[original].status:
                    raise ValueError('independent DESeq2 status differs')
                replayed.append(groups.iloc[original].group_id)
    return {'result': 'PASS', 'manifest_sha256': sha256(output / 'manifest.json'),
            'groups_replayed': len(groups), 'treated_samples': int(samples.condition.eq('treated').sum()),
            'gene_count': len(genes), 'independent_refits': replayed,
            'statuses': groups.status.value_counts().to_dict(),
            'tested_empty': int((groups.status.eq('eligible') & groups.n_significant.eq(0)).sum()),
            'multiplate_groups': int(groups.n_plates.gt(1).sum())}
