"""Bounded, step-owned predictor-input release for the PRnet entrant (#163).

The corpus is read only through the accepted `bases` public readers
(`bases.corpora.training_header`, `bases.fold_seed.FoldSeedIndex`) plus bounded
`h5py` row slabs of the corpus layers.  No dense wells x genes matrix is built:
only the declared bounded prefixes of the train/validation/test roles are read,
in 64-row slabs, as float32.  The Morgan-2048 representation is the accepted
chemical representation of the retired adapter input product: RDKit Morgan,
radius 2, 2048 bits, chirality on, binary, computed from the corpus's normalized
parent SMILES, keyed by `drug_key` and ordered by `drug_key`.

This module imports cleanly in the pinned PRnet environment because the
`bases` readers (which need the shared environment's newer `anndata`) are
imported lazily inside `build_release`.

    /home/chang/miniconda3/envs/arc/bin/python \
        execution/script/prepare_inputs.py --fold-row 0 --seed-row 0 \
        --train-wells 512 --validation-wells 128 --test-wells 64
"""

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import sys

import h5py
import numpy as np
import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from artifacts import sha256, write_json, write_rows  # noqa: E402

SCHEMA = 'stage2-02-models-tenv-inputs-v1'
TARGET_KIND = 'delta'
ROLES = ('train', 'validation', 'test')
CASE_FIELDS = ['case_id', 'drug_key', 'dose_um', 'group_id', 'control_group',
               'source_row', 'morgan_index']
MORGAN_SPEC = {
    'radius': 2,
    'bits': 2048,
    'chirality': True,
    'binary': True,
    'invariants': 'atom',
}
SLAB = 64
EXPECTED_GENES = 18164


class InputError(ValueError):
    """A declared input-contract precondition was violated."""


def repo_root():
    for candidate in Path(__file__).resolve().parents:
        if (candidate / '.git').exists():
            return candidate
    raise RuntimeError('repository root not found')


def public_readers():
    """The accepted Stage 2 foundation readers, imported lazily."""
    root = repo_root()
    sys.path.insert(0, str(root / 'stage2/01_bases/output/api'))
    from bases import corpora, fold_seed
    return corpora, fold_seed


def drug_fingerprints(obs):
    """Binary Morgan-2048 per `drug_key`, ordered by `drug_key`."""
    if not {'drug_key', 'parent_smiles'} <= set(obs.columns):
        raise InputError('drug table needs drug_key and parent_smiles')
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=MORGAN_SPEC['radius'], fpSize=MORGAN_SPEC['bits'],
        includeChirality=MORGAN_SPEC['chirality'])
    rows, fingerprints = [], []
    frame = obs.loc[:, ['drug_key', 'parent_smiles']].copy()
    if frame.drug_key.isna().any() or frame.parent_smiles.isna().any():
        raise InputError('missing drug identity or parent structure')
    for drug_key, group in frame.groupby('drug_key', sort=True, observed=True):
        smiles = sorted(set(str(value) for value in group.parent_smiles))
        if len(smiles) != 1:
            raise InputError(
                f'drug {drug_key} has inconsistent parent structures: {smiles}')
        molecule = Chem.MolFromSmiles(smiles[0])
        if molecule is None:
            raise InputError(f'drug {drug_key} has an unparsable parent structure')
        fingerprint = np.asarray(
            generator.GetFingerprintAsNumPy(molecule), dtype='uint8')
        if fingerprint.shape != (MORGAN_SPEC['bits'],) or not np.isin(
                fingerprint, (0, 1)).all():
            raise InputError(f'drug {drug_key} did not produce a binary fingerprint')
        rows.append({'drug_key': str(drug_key), 'morgan_index': len(rows)})
        fingerprints.append(fingerprint)
    if not rows:
        raise InputError('no drugs to fingerprint')
    return {'drugs': rows,
            'morgan': np.stack(fingerprints).astype('uint8'),
            'rdkit_version': rdBase.rdkitVersion}


def role_selection(cohorts, fold, limits):
    """Deterministic bounded prefixes of one fold's role cohorts.

    A role's wells are ordered by ascending `sample_id` and its prefix is the
    first ``limits[role]`` identities.  Selection never reads a name or a fold
    order beyond the requested fold number.
    """
    if not {'sample_id', 'fold', 'role'} <= set(cohorts.columns):
        raise InputError('cohort table needs sample_id, fold and role')
    if isinstance(fold, bool) or not isinstance(fold, int) or fold < 0:
        raise InputError(f'invalid fold identity: {fold!r}')
    unknown = set(cohorts.role.astype(str)) - set(ROLES)
    if unknown:
        raise InputError(f'unknown cohort roles: {sorted(unknown)}')
    if not {'train', 'validation', 'test'} <= set(cohorts.role.astype(str)):
        raise InputError('a fold must declare all three roles')
    if not isinstance(limits, dict) or set(limits) - set(ROLES):
        raise InputError(f'limits must be a mapping over {set(ROLES)}')
    scoped = cohorts.loc[cohorts.fold.astype(int) == fold]
    if scoped.empty:
        raise InputError(f'fold {fold} has no cohort rows')
    if scoped.sample_id.duplicated().any():
        raise InputError('duplicate sample identity inside a fold')
    selected = {}
    for role in ROLES:
        available = sorted(str(value) for value in
                           scoped.loc[scoped.role.astype(str) == role, 'sample_id'])
        if not available:
            raise InputError(f'fold {fold} has no {role} wells')
        if role in limits:
            limit = limits[role]
            if isinstance(limit, bool) or not isinstance(limit, int) or limit < 0:
                raise InputError(f'invalid {role} limit: {limit!r}')
            if limit > len(available):
                raise InputError(
                    f'{role} limit {limit} exceeds the {len(available)} available wells')
        else:
            limit = len(available)
        selected[role] = available[:limit]
    return selected


def _read_rows(dataset, positions):
    """Read selected rows of a dense corpus layer in bounded sorted slabs."""
    positions = np.asarray(positions, dtype='int64')
    if positions.ndim != 1 or len(positions) == 0:
        raise InputError('empty row selection')
    if positions.min() < 0 or positions.max() >= dataset.shape[0]:
        raise InputError('row selection out of range')
    values = np.empty((len(positions), dataset.shape[1]), dtype='float32')
    for start in range(0, len(positions), SLAB):
        part = positions[start:start + SLAB]
        order = np.argsort(part, kind='stable')
        block = np.asarray(dataset[part[order], :], dtype='float32')
        restored = np.empty_like(block)
        restored[order] = block
        values[start:start + len(part)] = restored
    if not np.isfinite(values).all():
        raise InputError('non-finite corpus layer value in the bounded selection')
    return values


def build_release(corpus, fold_product, output, *, fold_row, seed_row,
                  declared_folds=2, declared_seeds=2, limits):
    """Build the bounded predictor-input release and its integrity manifest."""
    corpora, fold_seed = public_readers()
    index = fold_seed.FoldSeedIndex.load(fold_product)
    folds = index.select_folds(declared_folds)
    seeds = index.select_seeds(declared_seeds)
    if not 0 <= fold_row < len(folds):
        raise InputError(f'fold row {fold_row} outside the declared prefix')
    if not 0 <= seed_row < len(seeds):
        raise InputError(f'seed row {seed_row} outside the declared prefix')
    fold = int(folds.iloc[fold_row]['fold'])
    seed = int(seeds.iloc[seed_row]['seed'])

    cohorts_path = Path(fold_product) / 'cohorts.csv'
    declared = (index.manifest.get('cohorts') or {}).get('sha256')
    if declared != sha256(cohorts_path):
        raise InputError('cohorts.csv digest differs from the fold product manifest')
    cohorts = pd.read_csv(cohorts_path, keep_default_na=False)
    selection = role_selection(cohorts, fold, limits)

    obs, var = corpora.training_header(corpus)
    if len(var) != EXPECTED_GENES:
        raise InputError(f'expected {EXPECTED_GENES} genes, found {len(var)}')
    obs = obs.reset_index(drop=True)
    positions = {}
    for role, sample_ids in selection.items():
        lookup = {str(value): row for row, value in enumerate(obs['sample_id'])}
        missing = [value for value in sample_ids if value not in lookup]
        if missing:
            raise InputError(f'{role} wells missing from the corpus: {missing[:3]}')
        positions[role] = np.asarray([lookup[value] for value in sample_ids],
                                     dtype='int64')

    combined = np.concatenate([positions[role] for role in ROLES])
    drugs = drug_fingerprints(obs.loc[combined, ['drug_key', 'parent_smiles']])
    drug_index = {row['drug_key']: row['morgan_index'] for row in drugs['drugs']}

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    files = {'drugs': 'drugs.csv', 'morgan': 'morgan.npy'}
    write_rows(output / 'drugs.csv', drugs['drugs'])
    np.save(output / 'morgan.npy', drugs['morgan'])
    with h5py.File(corpus, 'r') as handle:
        if set(handle['layers']) != {'control', 'delta'}:
            raise InputError('unexpected corpus layers')
        control = handle['layers/control']
        delta = handle['layers/delta']
        role_counts = {}
        for role in ROLES:
            rows = positions[role]
            role_counts[role] = int(len(rows))
            baseline = _read_rows(control, rows)
            np.save(output / f'{role}_baseline.npy', baseline)
            files[f'{role}_baseline'] = f'{role}_baseline.npy'
            if role != 'test':
                target = _read_rows(delta, rows)
                np.save(output / f'{role}_delta.npy', target)
                files[f'{role}_delta'] = f'{role}_delta.npy'
            records = []
            for row_position, sample_id in zip(rows, selection[role]):
                source = obs.iloc[int(row_position)]
                dose = float(source['dose_um'])
                if not np.isfinite(dose) or dose <= 0:
                    raise InputError(f'{role} well has an invalid treated dose')
                drug_key = str(source['drug_key'])
                if drug_key not in drug_index:
                    raise InputError(f'{role} well has no Morgan identity')
                records.append({
                    'case_id': str(sample_id),
                    'drug_key': drug_key,
                    'dose_um': dose,
                    'group_id': str(source['group_id']),
                    'control_group': str(source['control_group']),
                    'source_row': int(source['source_row']),
                    'morgan_index': drug_index[drug_key],
                })
            write_rows(output / f'{role}_cases.csv', records)
            files[f'{role}_cases'] = f'{role}_cases.csv'

    axis = '\n'.join(str(value) for value in var['hgnc_id']).encode()
    manifest = {
        'schema': SCHEMA,
        'target_kind': TARGET_KIND,
        'n_genes': int(len(var)),
        'gene_axis_sha256': hashlib.sha256(axis).hexdigest(),
        'identity': {
            'fold_row': fold_row, 'fold': fold, 'seed_row': seed_row, 'seed': seed,
        },
        'declared': {'folds': declared_folds, 'seeds': declared_seeds},
        'selection': {
            'order': 'ascending case_id prefix per role inside the fold',
            'limits': {role: len(selection[role]) for role in ROLES},
        },
        'roles': role_counts,
        'morgan': dict(MORGAN_SPEC),
        'files': files,
        'sha256': {name: sha256(output / relative)
                   for name, relative in files.items()},
        'source': {
            'corpus_sha256': sha256(corpus),
            'folds_sha256': sha256(Path(fold_product) / 'folds.csv'),
            'seeds_sha256': sha256(Path(fold_product) / 'seeds.csv'),
            'cohorts_sha256': declared,
            'fold_product_manifest_sha256': sha256(Path(fold_product) / 'manifest.json'),
        },
        'environment': {
            'python': sys.version,
            'executable': sys.executable,
            'numpy': np.__version__,
            'pandas': pd.__version__,
            'rdkit': drugs['rdkit_version'],
        },
        'created_utc': datetime.now(timezone.utc).isoformat(),
    }
    write_json(output / 'manifest.json', manifest)
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    root = repo_root()
    parser.add_argument('--corpus', default=str(
        root / 'stage2/01_bases/output/training-corpora/tenv.h5ad'))
    parser.add_argument('--fold-product', default=str(
        root / 'stage2/01_bases/output/fold_n_seed'))
    parser.add_argument('--output', default=str(
        root / 'stage2/02_models_tenv/output/bounded/inputs'))
    parser.add_argument('--fold-row', type=int, required=True)
    parser.add_argument('--seed-row', type=int, required=True)
    parser.add_argument('--folds', type=int, default=2)
    parser.add_argument('--seeds', type=int, default=2)
    parser.add_argument('--train-wells', type=int, default=512)
    parser.add_argument('--validation-wells', type=int, default=128)
    parser.add_argument('--test-wells', type=int, default=64)
    args = parser.parse_args(argv)
    manifest = build_release(
        args.corpus, args.fold_product, args.output,
        fold_row=args.fold_row, seed_row=args.seed_row,
        declared_folds=args.folds, declared_seeds=args.seeds,
        limits={'train': args.train_wells,
                'validation': args.validation_wells,
                'test': args.test_wells})
    print(json.dumps({key: manifest[key] for key in
                      ('n_genes', 'identity', 'roles')}, sort_keys=True))
    return 0


if __name__ == '__main__':
    sys.exit(main())
