"""Freeze the paired ADR0077 views into source-specific round1 inputs."""
import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP / 'output/api'))
sys.path.insert(0, str(STEP.parent / '01_bases/output/api'))
from bases import release, view
from tsg_round1.inputs import INPUT_SCHEMA, REQUIRED, SOURCE_POPULATIONS, validate_population
from tsg_round1.protocol import SPLIT_SEED, TRAINING_SEED, make_splits


def _freeze_source(source, destination, generation):
    from rdkit import Chem, DataStructs, rdBase
    from rdkit.Chem import rdFingerprintGenerator

    obs = source.obs.copy()
    obs['profile_id'] = obs.index
    obs['view_row'] = np.arange(len(obs))
    treated = np.flatnonzero(obs.role.eq('treated').to_numpy())
    population = obs.iloc[treated].copy()
    validate_population(source.source, population, source.axis_size)
    groups = source.control_groups()
    controls = [groups[str(group)] for group in population.matched_control_group]
    vehicles = np.unique(np.concatenate(controls))
    all_rows = np.r_[treated, vehicles]
    target_rows = np.full(len(obs), -1, dtype=int)
    target_rows[all_rows] = np.arange(len(all_rows))

    destination.mkdir()
    for name, reader in (('expression', source.X), ('control', source.control)):
        matrix = np.lib.format.open_memmap(destination / f'{name}.npy', mode='w+',
                                           dtype=np.float32, shape=(len(all_rows), source.axis_size))
        for start in range(0, len(obs), view.BATCH):
            stop = min(start + view.BATCH, len(obs))
            selected = target_rows[start:stop] >= 0
            if selected.any():
                matrix[target_rows[start:stop][selected]] = reader(start, stop)[selected]
        matrix.flush()

    population.to_csv(destination / 'population.csv', index=False)
    obs.iloc[all_rows].to_csv(destination / 'all_obs.csv', index=False)
    pd.DataFrame({'gene_index': source.gene_index, 'symbol': source.var.symbol}).to_csv(
        destination / 'genes.csv', index=False)
    make_splits(population).to_csv(destination / 'splits.csv', index=False)
    release.write_json(destination / 'matched_controls.json',
                       {str(i): target_rows[rows].tolist() for i, rows in enumerate(controls)})

    parents = population[['drug_parent_inchikey', 'drug_parent_smiles']].drop_duplicates().sort_values(
        'drug_parent_inchikey')
    if parents.drug_parent_inchikey.duplicated().any():
        raise ValueError(f'{source.source}: parent has conflicting structures')
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
    fingerprints = np.zeros((len(parents), 2048), dtype=np.uint8)
    for bits, smiles in zip(fingerprints, parents.drug_parent_smiles, strict=True):
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(f'{source.source}: invalid parent SMILES')
        DataStructs.ConvertToNumpyArray(generator.GetFingerprint(molecule), bits)
    np.save(destination / 'morgan.npy', fingerprints)
    parents.to_csv(destination / 'parents.csv', index=False)

    config = dict(source=source.source, view_generation=generation, split_seed=SPLIT_SEED,
                  training_seed=TRAINING_SEED, n_treated=len(population), n_parents=len(parents),
                  n_matched_vehicles=len(vehicles), preparation_sha256=release.sha256(Path(__file__)),
                  morgan=dict(rdkit=rdBase.rdkitVersion, radius=2, bits=2048, include_chirality=True))
    release.write_json(destination / 'config.json', config)
    release.write_json(destination / 'manifest.json', dict(
        schema=INPUT_SCHEMA, view_generation=generation,
        outputs={name: release.sha256(destination / name) for name in sorted(REQUIRED)}))


def prepare(view_dir, destination):
    """Verify the paired release before writing a fresh input candidate."""
    view_dir, destination = Path(view_dir), Path(destination)
    if destination.exists() and any(destination.iterdir()):
        raise ValueError('round1 inputs require a new empty candidate')
    if view_dir.resolve() in (destination.resolve(), *destination.resolve().parents):
        raise ValueError('round1 inputs cannot reside inside the source view')
    with view.open_view(view_dir) as paired:
        manifest = paired.manifest()
        if manifest.get('fixture') or paired.size != 5000:
            raise ValueError('formal round1 requires the paired 5k release')
        generation = {key: manifest[key] for key in ('schema', 'recipe', 'normalization', 'control', 'axis')}
        generation['manifest_sha256'] = release.sha256(view_dir / 'manifest.json')
        generation['payload_sha256'] = {source: paired.digests()[f'{source}.h5ad']
                                        for source in SOURCE_POPULATIONS}
        readers = [paired.source(source) for source in SOURCE_POPULATIONS]
        for reader in readers:
            validate_population(reader.source, reader.obs[reader.obs.role == 'treated'], reader.axis_size)
        destination.mkdir(parents=True, exist_ok=True)
        for reader in readers:
            _freeze_source(reader, destination / reader.source, generation)
    release.write_json(destination / 'manifest.json', dict(
        schema='tsg-dose-round1-axis-local-pair-v1', view_generation=generation,
        sources={name: release.sha256(destination / name / 'manifest.json')
                 for name in SOURCE_POPULATIONS}))
    return destination


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--view', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    args = parser.parse_args()
    prepare(args.view, args.candidate)
