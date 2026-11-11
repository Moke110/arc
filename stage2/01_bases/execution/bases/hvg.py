"""Seurat normalized dispersion with explicit, isolated training populations.

Only this module's train-only releases may select model genes. Full-tenv
rankings describe a population and are deliberately not selectable.
"""
import hashlib
import json
from pathlib import Path
import tempfile

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

SCHEMA = 'stage2-indexed-hvg-v1'
RECIPE = 'linear CP10K recovered by expm1 of published float32 log1p X; sample variance / mean; log dispersion normalized in 20 equal-width log1p-mean bins; gene-index tie-break'
BATCH = 64
ROLES = {'train', 'validation', 'test'}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8*1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def population_hash(sample_ids):
    ids = sorted(map(str, sample_ids))
    if len(ids) != len(set(ids)):
        raise ValueError('duplicate population sample identity')
    return hashlib.sha256(json.dumps(ids, separators=(',', ':')).encode()).hexdigest()


def make_cohorts(samples, folds=3, seed='bases-cold-parent-v1'):
    """New deterministic cold-parent folds; no historical partitions are reused.

    SHA256(seed|parent) orders parents. Equal contiguous buckets define test
    folds; the first ceil(10%) of remaining parents is validation. Every dose
    and every replicate of a parent has the same role within a fold.
    """
    if folds < 2 or not {'sample_id','drug_parent','dose_um'} <= set(samples):
        raise ValueError('cohort needs sample, parent-dose and at least two folds')
    parents = sorted(set(map(str, samples.drug_parent)),
                     key=lambda p: (hashlib.sha256(f'{seed}|{p}'.encode()).hexdigest(), p))
    if len(parents) < folds + 2:
        raise ValueError('not enough parents for train/validation/test cohorts')
    buckets = np.array_split(np.asarray(parents, dtype=object), folds)
    result = []
    for fold, bucket in enumerate(buckets):
        test = set(bucket)
        remaining = [p for p in parents if p not in test]
        validation = set(remaining[:max(1, int(np.ceil(len(remaining)*.1)))])
        frame = samples.copy()
        frame['fold'] = fold
        frame['role'] = ['test' if str(p) in test else 'validation' if str(p) in validation else 'train'
                         for p in frame.drug_parent]
        result.append(frame)
    return pd.concat(result, ignore_index=True)


def _header(path, expected_genes):
    with h5py.File(path, 'r') as f:
        samples, genes = read_elem(f['obs']), read_elem(f['var'])
        if (not {'sample_id','drug_parent','dose_um'} <= set(samples)
                or samples.sample_id.isna().any() or samples.sample_id.duplicated().any()
                or samples.drug_parent.isna().any() or samples.drug_parent.astype(str).eq('').any()
                or not np.isfinite(samples.dose_um).all() or samples.dose_um.le(0).any()
                or not np.array_equal(genes.gene_index.to_numpy(), np.arange(expected_genes))):
            raise ValueError('training samples or frozen indexed gene axis differ')
        shape = f['X'].shape if isinstance(f['X'], h5py.Dataset) else tuple(f['X'].attrs['shape'])
        if shape != (len(samples), expected_genes):
            raise ValueError('training expression shape differs')
    columns = ['sample_id','drug_parent','dose_um'] + (['group_id'] if 'group_id' in samples else [])
    return samples[columns].reset_index(drop=True), genes


def _validate_cohorts(samples, cohorts):
    required = {'sample_id', 'fold', 'role'}
    if not required <= set(cohorts) or cohorts[list(required)].isna().any().any():
        raise ValueError('cohort requires sample_id, fold, role')
    fold_values = pd.to_numeric(cohorts.fold, errors='raise').to_numpy()
    if not np.isfinite(fold_values).all() or np.any(fold_values != np.floor(fold_values)):
        raise ValueError('cohort fold must be an integer')
    if not set(cohorts.role) <= ROLES:
        raise ValueError('cohort has an unknown role')
    samples = samples.assign(sample_id=samples.sample_id.astype(str))
    result = []
    for fold, part in cohorts.groupby('fold', sort=True):
        if (part.sample_id.astype(str).duplicated().any()
                or set(part.sample_id.astype(str)) != set(samples.sample_id)):
            raise ValueError('cohort must partition the complete training population exactly once')
        merged = samples.merge(part[['sample_id','role']].assign(sample_id=part.sample_id.astype(str)),
                               on='sample_id', validate='one_to_one')
        dose_keys = merged.dose_um.map(lambda dose:format(float(dose),'.12g'))
        if merged.groupby([merged.drug_parent,dose_keys]).role.nunique().gt(1).any():
            raise ValueError('parent-dose group crosses cohort roles')
        if set(merged.role) != ROLES or merged.role.eq('train').sum() < 2:
            raise ValueError('cohort requires populated roles and at least two training wells')
        merged['fold'] = int(fold)
        result.append(merged)
    if not result:
        raise ValueError('empty cohort release')
    return pd.concat(result, ignore_index=True)


def _linear_batches(matrix, rows):
    """Read only bounded rows, including CSR-backed X without layer loading."""
    from scipy.sparse import csr_matrix
    for start in range(0, len(rows), BATCH):
        take = rows[start:start+BATCH]
        if isinstance(matrix, h5py.Dataset):
            value = matrix[take].astype(np.float32)
        else:
            n_genes = int(matrix.attrs['shape'][1])
            slabs = []
            for row in take:
                lo, hi = matrix['indptr'][row:row+2]
                slabs.append(csr_matrix((matrix['data'][lo:hi], matrix['indices'][lo:hi],
                                        [0, hi-lo]), shape=(1,n_genes)))
            from scipy.sparse import vstack
            value = vstack(slabs).toarray().astype(np.float32)
        if not np.isfinite(value).all() or np.any(value < 0):
            raise ValueError('HVG requires finite nonnegative log1p CP10K expression')
        linear = np.expm1(value)
        if (not np.isfinite(linear).all()
                or not np.allclose(linear.sum(axis=1,dtype=np.float64),10_000,rtol=1e-5,atol=.01)):
            raise ValueError('HVG input is not full-axis log1p CP10K')
        yield linear


def _ranking(matrix, rows, gene_count):
    if len(rows) < 2:
        raise ValueError('HVG requires at least two population wells')
    mean, m2 = np.zeros(gene_count), np.zeros(gene_count)
    total = 0
    for linear in _linear_batches(matrix, np.sort(rows)):
        batch_mean = linear.mean(axis=0, dtype=np.float64)
        centered = linear - batch_mean
        n, delta = len(linear), batch_mean - mean
        merged = total + n
        m2 += np.einsum('ij,ij->j', centered, centered) + delta**2 * (total*n/merged)
        mean += delta * (n/merged)
        total = merged
    safe_mean = np.where(mean == 0, 1e-12, mean)
    log_mean = np.log1p(safe_mean)
    dispersion = m2 / ((total-1)*safe_mean)
    log_dispersion = np.full(gene_count, np.nan)
    np.log(dispersion, out=log_dispersion, where=dispersion > 0)
    groups = pd.Series(log_dispersion).groupby(pd.cut(log_mean, 20), observed=True)
    avg, sd = groups.transform('mean').to_numpy(), groups.transform('std').to_numpy()
    score = np.full(gene_count, np.nan)
    np.divide(log_dispersion-avg, sd, out=score, where=sd > 0)
    score[(groups.transform('count').to_numpy() == 1) & np.isfinite(log_dispersion)] = 1
    order = np.argsort(-np.nan_to_num(score, nan=-np.inf), kind='stable')
    rank = np.empty_like(order)
    rank[order] = np.arange(1,gene_count+1)
    return pd.DataFrame({'gene_index':np.arange(gene_count), 'log_mean':log_mean,
                         'log_dispersion':log_dispersion, 'normalized_dispersion':score, 'rank':rank})


def build_hvg(training, output, *, cohorts=None, folds=3, seed='bases-cold-parent-v1', expected_genes=18164):
    """Publish descriptive and explicit train-only rankings, without truth access."""
    training, output = Path(training), Path(output)
    samples, _ = _header(training, expected_genes)
    generated = cohorts is None
    cohorts = _validate_cohorts(samples, make_cohorts(samples, folds, seed) if generated else cohorts)
    output.mkdir(parents=True, exist_ok=False)
    cohorts.to_csv(output/'cohorts.csv', index=False)
    populations = {'descriptive': {'kind':'descriptive', 'sample_ids':sorted(samples.sample_id.astype(str))}}
    for fold, part in cohorts.groupby('fold', sort=True):
        populations[f'fold_{int(fold)}'] = {'kind':'train-only', 'fold':int(fold),
            'sample_ids':sorted(part.loc[part.role.eq('train'),'sample_id'].astype(str))}
    lookup = pd.Series(np.arange(len(samples)), index=samples.sample_id.astype(str))
    with h5py.File(training, 'r') as f:
        for name, population in populations.items():
            ids = population.pop('sample_ids')
            rows = lookup.loc[ids].to_numpy(dtype=int)
            _ranking(f['X'], rows, expected_genes).to_csv(output/f'{name}.csv', index=False)
            population.update(n_wells=len(ids), population_sha256=population_hash(ids),
                              ranking=f'{name}.csv', sha256=sha256(output/f'{name}.csv'))
    manifest = {'schema':SCHEMA, 'training_sha256':sha256(training), 'gene_count':expected_genes,
                'recipe':RECIPE, 'cohorts':{'file':'cohorts.csv','sha256':sha256(output/'cohorts.csv'),
                    'algorithm':'SHA256(seed|parent) sorted buckets; 10% remaining parents validation' if generated else 'explicit supplied complete group-isolated partition',
                    'seed':seed if generated else None},
                'populations':populations}
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True)+'\n')
    return manifest


class HVGRelease:
    def __init__(self, root):
        self.root = Path(root)
        self.manifest = json.loads((self.root/'manifest.json').read_text())
        if self.manifest.get('schema') != SCHEMA:
            raise ValueError('unknown HVG schema')
        cohort = self.manifest['cohorts']
        if sha256(self.root/cohort['file']) != cohort['sha256']:
            raise ValueError('HVG cohort seal differs')
        self.cohorts = pd.read_csv(self.root/cohort['file'],keep_default_na=False,dtype={'sample_id':str,'drug_parent':str})
        base_samples = self.cohorts[['sample_id','drug_parent','dose_um']].drop_duplicates()
        if base_samples.sample_id.duplicated().any():
            raise ValueError('HVG sample metadata differs across folds')
        _validate_cohorts(base_samples, self.cohorts)
        self._rankings = {}
        for name, population in self.manifest['populations'].items():
            if population['kind'] == 'train-only':
                ids = self.cohorts.loc[self.cohorts.fold.eq(population['fold']) & self.cohorts.role.eq('train'),'sample_id']
            elif name == 'descriptive' and population['kind'] == 'descriptive':
                ids = base_samples.sample_id
            else:
                raise ValueError('unknown HVG population kind')
            if len(ids) != population['n_wells'] or population_hash(ids) != population['population_sha256']:
                raise ValueError('HVG population differs from its cohort')
            path = self.root/population['ranking']
            if sha256(path) != population['sha256']:
                raise ValueError('HVG ranking seal differs')
            frame = pd.read_csv(path)
            n = self.manifest['gene_count']
            if (not np.array_equal(frame.gene_index,np.arange(n))
                    or sorted(frame['rank']) != list(range(1,n+1))):
                raise ValueError('invalid indexed HVG ranking')
            self._rankings[name] = frame

    def ranking(self, population):
        return self._rankings[population].copy()

    def select(self, population, size, *, training_sample_ids):
        info = self.manifest['populations'][population]
        if info['kind'] != 'train-only':
            raise ValueError('descriptive HVG population cannot select model genes')
        if population_hash(training_sample_ids) != info['population_sha256']:
            raise ValueError('model training population differs from HVG population')
        if isinstance(size, bool) or not isinstance(size,(int,np.integer)) or not 0 < size <= self.manifest['gene_count']:
            raise ValueError('invalid HVG selection size')
        return self._rankings[population].sort_values('rank').gene_index.to_numpy(dtype=int)[:size]


def verify_hvg(training, release):
    """Recompute published ranks from the exact expression and cohort population."""
    current = HVGRelease(release)
    if sha256(training) != current.manifest['training_sha256']:
        raise ValueError('HVG training expression seal differs')
    with tempfile.TemporaryDirectory() as tmp:
        rebuilt = build_hvg(training,Path(tmp)/'hvg',cohorts=current.cohorts,
                            expected_genes=current.manifest['gene_count'])
        for name, info in current.manifest['populations'].items():
            other = rebuilt['populations'][name]
            if info['sha256'] != other['sha256'] or info['population_sha256'] != other['population_sha256']:
                raise ValueError('HVG source replay differs')
    return {'result':'PASS', 'gene_count':current.manifest['gene_count'],
            'populations':{name:info['n_wells'] for name,info in current.manifest['populations'].items()},
            'manifest_sha256':sha256(Path(release)/'manifest.json')}
