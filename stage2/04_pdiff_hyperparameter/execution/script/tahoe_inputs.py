"""Build and verify the Step 04 Tahoe view of the frozen pseudobulk corpus."""
from pathlib import Path
import shutil
import tempfile

from anndata.io import read_elem
import h5py
import numpy as np
import pandas as pd
from rdkit import rdBase
from rdkit.Chem import rdFingerprintGenerator
from scipy.sparse import csr_matrix

from artifacts import read_json, sha256, write_json
from inputs import CASE_COLUMNS, SEEDS, normalized_structure, split_drugs
from metrics import truth_labels
from training_inputs import STEP03, json_hash, predictor_paths, request_for, verify_files

STEP02 = Path('stage2/02_dataset_corpus')
CORPUS = STEP02 / 'output/corpus/v2/tahoe.h5ad'
GENES = CORPUS.parent / 'tables/genes.csv'
DRUGS = Path('pert_data/Tahoe-100M/source/metadata/drug_metadata.parquet')
ROLES = ('train', 'validation', 'test')
N_GENES = 18164
ZERO_GENE = 'ENSG00000182230'
ROW_BATCH = 128
EXPRESSION = 'project native log1p_cp10k onto the CM axis; subtract same-plate, same-cell-line DMSO_TF mean'
CASE_METADATA = [*CASE_COLUMNS, 'plate', 'well', 'cell_line_id']


def source_files(runtime, *, verify_corpus=True):
    reference = runtime / STEP03 / 'output/input/manifest.json'
    reference_pin = read_json(runtime / STEP03 / 'execution/input_record/verification.json')['manifest_sha256']
    if sha256(reference) != reference_pin:
        raise ValueError('Step 03 reference manifest differs from the frozen identity')
    genes_pin = read_json(reference)['files']['genes.csv']
    if sha256(reference.with_name('genes.csv')) != genes_pin or sha256(runtime / GENES) != genes_pin:
        raise ValueError('Tahoe gene axis differs from the frozen ordered CM axis')
    corpus_pin = read_json(runtime / STEP02 / 'execution/build_record/tahoe_corpus_move.json')['file_sha256']
    vendor_manifest = runtime / DRUGS.parents[2] / 'manifest.sha256'
    drug_pin = dict(line.split(maxsplit=1)[::-1] for line in vendor_manifest.read_text().splitlines())[
        'metadata/drug_metadata.parquet']
    pins = {str(CORPUS): corpus_pin, str(GENES): genes_pin, str(DRUGS): drug_pin}
    for relative in ((CORPUS, DRUGS) if verify_corpus else (DRUGS,)):
        if sha256(runtime / relative) != pins[str(relative)]:
            raise ValueError(f'Tahoe source hash mismatch: {relative}')
    return {'corpus_version': 'v2', 'sha256': pins, 'step03_manifest_sha256': reference_pin}


def gene_projection(genes, var):
    target = list(genes.ensembl_id)
    if len(target) != N_GENES or len(set(target)) != N_GENES or '' in target:
        raise ValueError('invalid ordered CM gene axis')
    selected = set(target)
    columns = {}
    for index, gene in enumerate(var.ensembl_id.astype(str)):
        if gene in selected:
            if gene in columns:
                raise ValueError(f'ambiguous Tahoe gene identity: {gene}')
            columns[gene] = index
    if selected - columns.keys() != {ZERO_GENE}:
        raise ValueError('Tahoe must retain exactly the historical FAM153B zero slot')
    present = np.array([i for i, gene in enumerate(target) if gene in columns])
    return present, np.array([columns[target[i]] for i in present]), target.index(ZERO_GENE)


def source_metadata(obs, drug_metadata):
    required = ['plate', 'well', 'cell_line_id', 'drug', 'is_control', 'dose_um']
    if obs.index.name != 'profile_id' or obs.index.has_duplicates or obs[required].isna().any().any():
        raise ValueError('invalid Tahoe pseudobulk metadata')
    obs = obs.reset_index().copy()
    for field in ('profile_id', 'plate', 'well', 'cell_line_id', 'drug'):
        obs[field] = obs[field].astype(str)
    if not obs.profile_id.eq(obs.plate + ':' + obs.well + ':' + obs.cell_line_id).all():
        raise ValueError('invalid Tahoe profile identity')
    obs['drug'] = obs.drug.str.strip()
    blank = obs.get('is_blank', pd.Series(False, index=obs.index)).astype(bool)
    control = obs.is_control.astype(bool)
    dmso = ~blank & control & obs.drug.eq('DMSO_TF')
    if (~blank & ~control & obs.drug.eq('DMSO_TF')).any():
        raise ValueError('DMSO_TF profile declared as treatment')
    treated = obs.loc[~blank & ~control]
    controls, control_rows, control_indices = [], [], {}
    for (plate, line), part in obs.loc[dmso].groupby(['plate', 'cell_line_id'], sort=True):
        index = len(controls)
        controls.append(dict(control_index=index, plate=plate, cell_line_id=line,
                             n_wells=len(part), control_case_ids=json_hash(part.profile_id.tolist())))
        control_rows.append(part.index.to_numpy())
        control_indices[plate, line] = index

    metadata = drug_metadata.dropna(subset=['drug']).copy()
    metadata['drug'] = metadata.drug.astype(str).str.strip()
    metadata = metadata.groupby('drug').canonical_smiles.agg(lambda values: sorted(set(values.dropna().str.strip())))
    structures, structure_errors = {}, {}
    for drug in sorted(treated.drug.unique()):
        candidates = [value for value in metadata.get(drug, []) if value]
        if not candidates:
            structure_errors[drug] = 'missing_structure'
            continue
        try:
            parents = [normalized_structure(value) for value in candidates]
        except ValueError:
            structure_errors[drug] = 'invalid_structure'
            continue
        if len({parent for _, parent, _ in parents}) != 1:
            raise ValueError(f'inconsistent normalized Tahoe structure: {drug}')
        structures[drug] = parents[0]

    eligible, exclusions = [], []
    for row in treated.itertuples():
        match = control_indices.get((row.plate, row.cell_line_id))
        reason = structure_errors.get(row.drug)
        if not np.isfinite(row.dose_um) or row.dose_um <= 0:
            reason = 'invalid_treated_dose'
        if match is None and reason is None:
            reason = 'missing_matched_control'
        if reason:
            exclusions.append(dict(case_id=row.profile_id, drug_key=row.drug, plate=row.plate,
                                   cell_line_id=row.cell_line_id, reason=reason))
            continue
        eligible.append(dict(case_id=row.profile_id, drug_key=row.drug, dose_um=float(row.dose_um),
                             control_index=match,
                             group_id=json_hash([row.drug, float(row.dose_um), row.cell_line_id])[:24],
                             plate=row.plate, well=row.well, cell_line_id=row.cell_line_id,
                             source_row=row.Index))
    cases = pd.DataFrame(eligible)
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
    drugs, fingerprints = [], []
    for drug in sorted(cases.drug_key.unique()):
        mol, parent, connectivity = structures[drug]
        drugs.append(dict(drug_key=drug, parent_smiles=parent, parent_connectivity=connectivity,
                          morgan_index=len(drugs)))
        fingerprints.append(generator.GetFingerprintAsNumPy(mol))
    counts = {'raw': len(obs), 'declared_treated': len(treated), 'included_treated': len(cases),
              'declared_dmso_tf': int(dmso.sum()), 'other_controls': int((~blank & control & ~dmso).sum()),
              'blank': int(blank.sum()), 'control_strata': len(controls),
              'declared_drugs': int(treated.drug.nunique()), 'structure_eligible_drugs': len(structures),
              'matched_drugs': len(drugs)}
    excluded = pd.DataFrame(exclusions, columns=['case_id', 'drug_key', 'plate', 'cell_line_id', 'reason'])
    counts['excluded_treated_by_reason'] = excluded.reason.value_counts().sort_index().to_dict()
    return (cases, pd.DataFrame(controls), control_rows, pd.DataFrame(drugs),
            np.asarray(fingerprints, dtype=np.uint8), excluded, counts)


def batches(length):
    for start in range(0, length, ROW_BATCH):
        yield slice(start, min(start + ROW_BATCH, length))


def projected_rows(layer, native_shape, source_rows, present, native_cols):
    if layer.attrs['encoding-type'] != 'csr_matrix' or tuple(layer.attrs['shape']) != native_shape:
        raise ValueError('Tahoe log1p_cp10k must be the native CSR layer')
    ptr = layer['indptr'][:]
    for part in batches(len(source_rows)):
        requested = source_rows[part]
        values = np.zeros((len(requested), N_GENES), dtype=np.float32)
        # Never read CSR data for a native row span larger than the chunk bound.
        boundaries = np.flatnonzero(np.diff(requested // ROW_BATCH)) + 1
        for indices in np.split(np.arange(len(requested)), boundaries):
            rows = requested[indices]
            start, stop = rows[0], rows[-1] + 1
            lo, hi = ptr[start], ptr[stop]
            slab = csr_matrix((layer['data'][lo:hi], layer['indices'][lo:hi], ptr[start:stop + 1] - lo),
                              shape=(stop - start, native_shape[1]))
            values[np.ix_(indices, present)] = slab[rows - start][:, native_cols].toarray()
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError('invalid native log1p_cp10k values')
        yield part, values


def fold_cases(cases, split, fold):
    membership = split[split.fold.eq(fold)].set_index('drug_key').role
    roles = cases.drug_key.map(membership)
    for role in ROLES:
        yield role, cases[roles.eq(role)]


def save_folds(root, cases, drugs, delta):
    folds = split_drugs(drugs)
    folds.to_csv(root / 'folds.csv', index=False)
    coverage = []
    for fold in range(3):
        base = root / f'fold_{fold}'
        (base / 'evaluator').mkdir(parents=True)
        paths = predictor_paths(fold)
        fold_counts = {'fold': fold}
        for role, selected in fold_cases(cases, folds, fold):
            if selected.empty:
                raise ValueError(f'empty Tahoe fold {fold} {role}')
            selected[CASE_COLUMNS].to_csv(root / paths[f'{role}_cases'], index=False)
            fold_counts[role] = {'drugs': int(selected.drug_key.nunique()),
                                 'wells': len(selected), 'groups': int(selected.group_id.nunique())}
            if role == 'test':
                save_truth(base / 'evaluator', selected, delta)
            else:
                target = np.lib.format.open_memmap(root / paths[f'{role}_expression'], mode='w+',
                                                  dtype='float32', shape=(len(selected), N_GENES))
                indices = selected.index.to_numpy()
                for part in batches(len(selected)):
                    target[part] = delta[indices[part]]
                target.flush()
                del target
        coverage.append(fold_counts)
    return coverage


def evaluator_groups(selected):
    groups = selected.drop_duplicates('group_id')[['group_id', 'drug_key', 'dose_um', 'cell_line_id']].copy()
    groups['n_wells'] = groups.group_id.map(selected.group_id.value_counts())
    return groups


def save_truth(directory, selected, delta):
    groups = evaluator_groups(selected)
    members = selected.groupby('group_id', sort=False).groups
    groups.to_csv(directory / 'groups.csv', index=False)
    shape = (len(groups), N_GENES)
    truth, positive, negative = [
        np.lib.format.open_memmap(directory / f'{name}.npy', mode='w+', dtype=dtype, shape=shape)
        for name, dtype in [('truth', 'float32'), ('positive', 'bool'), ('negative', 'bool')]
    ]
    for index, group in enumerate(groups.group_id):
        total = np.zeros(N_GENES, dtype=np.float64)
        for part in batches(len(members[group])):
            total += delta[members[group][part]].sum(axis=0, dtype=np.float64)
        value = (total / len(members[group])).astype(np.float32)
        truth[index] = value
        positive[index], negative[index] = truth_labels(value)
    for array in (truth, positive, negative):
        array.flush()


def build_tahoe(runtime, output):
    output = output.resolve()
    if output.exists():
        raise ValueError('Tahoe release requires a fresh output directory')
    sources = source_files(runtime)
    genes = pd.read_csv(runtime / GENES, keep_default_na=False)
    drug_metadata = pd.read_parquet(runtime / DRUGS, columns=['drug', 'canonical_smiles'])
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.tahoe-build-', dir=output.parent) as temporary:
        root = Path(temporary) / 'release'
        root.mkdir()
        with h5py.File(runtime / CORPUS, 'r') as source:
            obs, var = read_elem(source['obs']), read_elem(source['var'])
            present, native_cols, zero_index = gene_projection(genes, var)
            cases, controls, control_rows, drugs, fingerprints, exclusions, counts = source_metadata(obs, drug_metadata)
            layer = source['layers/log1p_cp10k']
            native_shape = len(obs), len(var)
            control_values = np.zeros((len(controls), N_GENES), dtype=np.float32)
            control_sums = np.zeros_like(control_values, dtype=np.float64)
            raw_controls = np.concatenate(control_rows)
            indices = np.repeat(np.arange(len(controls)), controls.n_wells)
            order = np.argsort(raw_controls)
            for part, values in projected_rows(layer, native_shape, raw_controls[order], present, native_cols):
                np.add.at(control_sums, indices[order[part]], values)
            control_values[:] = control_sums / controls.n_wells.to_numpy()[:, None]
            np.save(root / 'controls.npy', control_values)
            delta = np.lib.format.open_memmap(Path(temporary) / 'delta.npy', mode='w+', dtype='float32',
                                             shape=(len(cases), N_GENES))
            for part, values in projected_rows(layer, native_shape, cases.source_row.to_numpy(), present, native_cols):
                delta[part] = values - control_values[cases.control_index.to_numpy()[part]]
            shutil.copyfile(runtime / GENES, root / 'genes.csv')
            cases[CASE_METADATA].to_csv(root / 'wells.csv', index=False)
            controls.to_csv(root / 'controls.csv', index=False)
            exclusions.to_csv(root / 'exclusions.csv', index=False)
            drugs.to_csv(root / 'drugs.csv', index=False)
            np.save(root / 'morgan.npy', fingerprints)
            coverage = save_folds(root, cases, drugs, delta)
            del delta
        files = {str(path.relative_to(root)): sha256(path) for path in sorted(root.rglob('*')) if path.is_file()}
        manifest = {'schema': 'tahoe-input-v1', 'target_kind': 'delta', 'gene_count': N_GENES,
                    'source_sha256': sources['sha256'][str(CORPUS)], 'sources': sources,
                    'gene_projection': {'axis': 'ordered CM Ensembl IDs', 'mapped': len(present),
                                        'zero_gene': ZERO_GENE, 'zero_index': zero_index},
                    'expression': EXPRESSION, 'controls': 'nonblank is_control rows with drug DMSO_TF',
                    'grouping': 'drug × dose_um × cell_line_id', 'row_counts': counts,
                    'drug_count': len(drugs), 'well_count': len(cases),
                    'group_count': int(cases.group_id.nunique()), 'folds': coverage,
                    'split_seed': 20261105, 'training_seeds': SEEDS,
                    'split_rule': 'Step 03 greedy outer block balance; seeded ~15% validation from remaining drugs',
                    'rdkit_version': rdBase.rdkitVersion,
                    'morgan': {'radius': 2, 'bits': 2048, 'chirality': True, 'binary': True},
                    'normalization': 'Cleanup; FragmentParent; Uncharger; canonical isomeric SMILES; parent InChIKey first block',
                    'builder_sha256': sha256(Path(__file__)), 'files': files}
        write_json(root / 'manifest.json', manifest)
        root.rename(output)
    return manifest


def release_manifest(root):
    manifest = read_json(root / 'manifest.json')
    if (manifest['schema'], manifest['target_kind'], manifest['gene_count']) != ('tahoe-input-v1', 'delta', N_GENES):
        raise ValueError('unexpected Tahoe release contract')
    return manifest


def prepare_tahoe(runtime, input_root, output, fold):
    root = input_root.resolve()
    manifest = release_manifest(root)
    paths = predictor_paths(fold)
    hashes = {key: manifest['files'][relative] for key, relative in paths.items()}
    consumer = {'schema': 'step04-input-v1', 'dataset': 'tahoe', 'fold': fold, 'n_genes': N_GENES,
                'target_kind': 'delta', 'source_input_sha256': sha256(root / 'manifest.json'),
                'corpus_sha256': manifest['source_sha256'], 'preprocessing': manifest['expression'],
                'split_sha256': manifest['files']['folds.csv'], 'gene_sha256': hashes['genes'],
                'predictor_files': paths, 'file_sha256': hashes}
    request = request_for(root, consumer)
    from training_inputs import verify_request
    verify_request(request, runtime)
    write_json(output.with_suffix('.input.json'), consumer)
    write_json(output, request)


def verify_table(path, expected, description):
    if path.read_text() != expected.to_csv(index=False):
        raise ValueError(f'Tahoe {description} differs from the frozen rule')


def verify_predictor_metadata(runtime, root, consumer):
    sources = source_files(runtime, verify_corpus=False)
    if consumer['gene_sha256'] != sources['sha256'][str(GENES)]:
        raise ValueError('Tahoe gene axis differs from the frozen ordered CM axis')
    paths, hashes = consumer['predictor_files'], consumer['file_sha256']
    verify_files(root, {'manifest': 'manifest.json', 'split': 'folds.csv'},
                 {'manifest': consumer['source_input_sha256'], 'split': consumer['split_sha256']})
    manifest = release_manifest(root)
    if (manifest['sources'] != sources
            or manifest['source_sha256'] != sources['sha256'][str(CORPUS)]
            or consumer['corpus_sha256'] != manifest['source_sha256']
            or consumer['preprocessing'] != EXPRESSION or manifest['expression'] != EXPRESSION
            or consumer['split_sha256'] != manifest['files']['folds.csv']
            or any(hashes[key] != manifest['files'][name] for key, name in paths.items())):
        raise ValueError('Tahoe request differs from frozen source and release identities')
    # Only the pseudobulk observation metadata is needed to reconstruct the
    # pairing and fold. Full corpus hashing belongs to build and verify-tahoe.
    with h5py.File(runtime / CORPUS, 'r') as source:
        obs = read_elem(source['obs'])
    drug_metadata = pd.read_parquet(runtime / DRUGS, columns=['drug', 'canonical_smiles'])
    cases, controls, _, drugs, _, _, _ = source_metadata(obs, drug_metadata)
    split = split_drugs(drugs)
    verify_table(root / 'controls.csv', controls, 'matched control metadata')
    verify_table(root / paths['drugs'], drugs, 'drug identities')
    verify_table(root / 'folds.csv', split, 'parent-connectivity block membership')
    for role, selected in fold_cases(cases, split, consumer['fold']):
        verify_table(root / paths[f'{role}_cases'], selected[CASE_COLUMNS], 'case membership or control pairing')


def verify_matrix(path, shape, dtype, zero_index):
    matrix = np.load(path, mmap_mode='r', allow_pickle=False)
    if matrix.shape != shape or matrix.dtype != np.dtype(dtype):
        raise ValueError(f'Tahoe shape/dtype mismatch: {path.name}')
    if zero_index is not None and np.any(matrix[:, zero_index]):
        raise ValueError(f'Tahoe historical FAM153B slot is nonzero: {path.name}')
    return matrix


def verify_tahoe(runtime, release):
    root = release.resolve()
    manifest = release_manifest(root)
    relative_files = {str(path.relative_to(root)) for path in root.rglob('*') if path.is_file()} - {'manifest.json'}
    if relative_files != manifest['files'].keys():
        raise ValueError('Tahoe manifest inventory mismatch')
    verify_files(root, {path: path for path in relative_files}, manifest['files'])
    sources = source_files(runtime)
    if manifest['sources'] != sources or manifest['source_sha256'] != sources['sha256'][str(CORPUS)]:
        raise ValueError('Tahoe source identities differ from the frozen inputs')
    genes = pd.read_csv(runtime / GENES, keep_default_na=False)
    if sha256(root / 'genes.csv') != sources['sha256'][str(GENES)]:
        raise ValueError('Tahoe gene axis differs from the frozen ordered CM axis')
    with h5py.File(runtime / CORPUS, 'r') as source:
        obs, var = read_elem(source['obs']), read_elem(source['var'])
    present, _, zero_index = gene_projection(genes, var)
    expected_projection = {'axis': 'ordered CM Ensembl IDs', 'mapped': len(present),
                           'zero_gene': ZERO_GENE, 'zero_index': zero_index}
    if manifest['gene_projection'] != expected_projection:
        raise ValueError('Tahoe gene projection identity mismatch')
    drug_metadata = pd.read_parquet(runtime / DRUGS, columns=['drug', 'canonical_smiles'])
    cases, controls, _, drugs, fingerprints, exclusions, counts = source_metadata(obs, drug_metadata)
    for name, frame, description in [('wells', cases[CASE_METADATA], 'matched treatment metadata'),
                                     ('controls', controls, 'matched control metadata'),
                                     ('drugs', drugs, 'drug identities'),
                                     ('exclusions', exclusions, 'exclusion metadata')]:
        verify_table(root / f'{name}.csv', frame, description)
    if not np.array_equal(np.load(root / 'morgan.npy', mmap_mode='r', allow_pickle=False), fingerprints):
        raise ValueError('Tahoe fingerprints differ from the frozen rule')
    verify_matrix(root / 'controls.npy', (len(controls), N_GENES), 'float32', zero_index)
    split = split_drugs(drugs)
    verify_table(root / 'folds.csv', split, 'parent-connectivity block membership')
    expected_counts = []
    for fold in range(3):
        paths = predictor_paths(fold)
        fold_counts = {'fold': fold}
        for role, selected in fold_cases(cases, split, fold):
            verify_table(root / paths[f'{role}_cases'], selected[CASE_COLUMNS], 'case membership or control pairing')
            fold_counts[role] = {'drugs': int(selected.drug_key.nunique()),
                                 'wells': len(selected), 'groups': int(selected.group_id.nunique())}
            if role != 'test':
                verify_matrix(root / paths[f'{role}_expression'], (len(selected), N_GENES), 'float32', zero_index)
                continue
            directory = root / f'fold_{fold}/evaluator'
            groups = evaluator_groups(selected)
            verify_table(directory / 'groups.csv', groups, 'evaluator group membership')
            shape = (len(groups), N_GENES)
            for name, dtype in [('truth', 'float32'), ('positive', 'bool'), ('negative', 'bool')]:
                verify_matrix(directory / f'{name}.npy', shape, dtype, zero_index)
        expected_counts.append(fold_counts)
    if manifest['row_counts'] != counts or manifest['folds'] != expected_counts:
        raise ValueError('Tahoe coverage differs from the frozen inputs')
    return manifest
