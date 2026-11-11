"""Exact, versioned reuse of accepted measured truth with rebuilt training inputs.

No fit is performed. Path changes, raw hashes and approximate expression are
never treated as mathematical equivalence. The original truth stays immutable.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

import h5py
from anndata.io import read_elem
import numpy as np
import pandas as pd

from .corpora import NORMALIZATION, training_header
from .significance import (SCHEMA as TRUTH_SCHEMA, PartitionTruth, checked_file, group_key, membership,
                           sha256, write_json)
from .truth_build import WORKER, VERIFICATION_SCHEMA, raw_axis_rows, source_plan, verify as verify_truth

SCHEMA = 'stage2-truth-successor-v1'
FULL_TENV_SAMPLES = 4206
FULL_TENV_GROUPS = 1377
DEFAULT_RSCRIPT = Path('/home/chang/miniconda3/envs/biomni_e1/bin/Rscript')
REPLICATION = 'at least 2 treated and 2 matched-control wells; residual design degrees of freedom'


def recipe():
    return {'gene_count': 18164, 'fdr': .05, 'test': 'two-sided DESeq2 Wald',
            'method_description': 'raw-count DESeq2 two-sided Wald; within-contrast BH over exactly indexed axis',
            'lfc_threshold': 0, 'fit_type': 'parametric', 'independent_filtering': False,
            'bh_family': 'entire fixed indexed axis, including untestable genes',
            'plate_design': '~ plate + condition for multi-plate; otherwise ~ condition',
            'size_factors': 'ratio or poscounts fallback', 'replication_rule': REPLICATION,
            'prediction_significance': False, 'normalization': NORMALIZATION}


def runtime_seal(rscript):
    rscript = Path(rscript).resolve(strict=True)
    command = [str(rscript), '--vanilla', '-e',
               'cat(as.character(getRversion()), as.character(packageVersion("DESeq2")), sep="\\n")']
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    versions = result.stdout.strip().splitlines()
    if len(versions) != 2:
        raise ValueError('unexpected R/DESeq2 runtime version report')
    return {'r_version': versions[0], 'deseq2_version': versions[1],
            'rscript_sha256': sha256(rscript), 'rscript_path': str(rscript)}


def binding(path):
    path = Path(path).resolve(strict=True)
    return {'path': str(path), 'sha256': sha256(path)}


def producer_binding(native):
    path = Path(native).parent / 'tenv_manifest.json'
    manifest = json.loads(path.read_text())
    if manifest.get('files', {}).get('tenv.h5ad') != sha256(native):
        raise ValueError('native producer manifest does not seal these count bytes')
    return binding(path)


def read_ids(file, ids, columns, positions):
    """Read at most one contrast or row batch, restoring requested stable-ID order."""
    rows = np.array([positions[value] for value in ids], dtype=int)
    if len(set(rows)) != len(rows):
        raise ValueError('duplicate native profile in exact count comparison')
    order = np.argsort(rows, kind='stable')
    values = np.concatenate(list(raw_axis_rows(file['X'], rows[order], columns)))
    return values[np.argsort(order)]


def normalized(counts):
    values = counts.astype('float64')
    totals = values.sum(axis=1, dtype='float64')
    if np.any(totals <= 0):
        raise ValueError('invalid zero library in exact training replay')
    return np.log1p(values * (10000 / totals[:, None])).astype('float32')


def array_rows(dataset, positions):
    positions = np.asarray(positions, dtype=int)
    order = np.argsort(positions, kind='stable')
    return dataset[positions[order]][np.argsort(order)]


def sample_identity(old, new):
    if (old.sample_id.duplicated().any() or new.sample_id.duplicated().any()
            or set(old.sample_id) != set(new.sample_id)):
        raise ValueError('successor sample identity differs')
    fields = ['sample_id', 'source_profile_id', 'drug_parent', 'dose_um']
    fields += [field for field in ('source', 'group_id', 'cell_identity', 'time_h', 'vehicle',
                                  'source_inchikey', 'source_parent_key') if field in old]
    if any(field not in new for field in fields):
        raise ValueError('successor sample metadata field is missing')
    left = old.set_index('sample_id', drop=False).sort_index()[fields].copy()
    right = new.set_index('sample_id', drop=False).sort_index()[fields].copy()
    for field in ('dose_um', 'time_h'):
        if field in left:
            left[field], right[field] = pd.to_numeric(left[field]), pd.to_numeric(right[field])
    try:
        pd.testing.assert_frame_equal(left, right, check_dtype=False, check_categorical=False)
    except AssertionError as error:
        raise ValueError('successor sample identity/covariates differ') from error
    for frame in (old, new):
        if 'group_id' in frame and any(group_key(r.drug_parent, r.dose_um) != r.group_id
                                       for r in frame.itertuples(index=False)):
            raise ValueError('training parent-dose key differs from the statistical contrast')
    controls = {}
    for row in old.itertuples(index=False):
        controls[row.sample_id] = sorted(json.loads(row.matched_control_ids))
    for row in new.itertuples(index=False):
        if sorted(json.loads(row.matched_control_ids)) != controls[row.sample_id]:
            raise ValueError('successor matched-control membership differs')
    text = left.to_csv(index=False, lineterminator='\n')
    return hashlib.sha256(text.encode()).hexdigest()


def compare_layers(original_training, training, old_native, columns, original_obs, new_obs):
    """Exact finite values in every layer, also independently replayed from integers."""
    digests = {name: hashlib.sha256() for name in ('X', 'control', 'delta')}
    with h5py.File(old_native) as raw, h5py.File(original_training / 'tenv.h5ad') as old, \
            h5py.File(training / 'tenv.h5ad') as new:
        raw_obs = read_elem(raw['obs'])
        source_positions = {value: index for index, value in enumerate(raw_obs.index)}
        new_positions = {value: index for index, value in enumerate(new_obs.sample_id)}
        means = {}
        for start in range(0, len(original_obs), 32):
            part = original_obs.iloc[start:start + 32]
            expected_x = normalized(read_ids(raw, part.source_profile_id.tolist(), columns, source_positions))
            controls = []
            for row in part.itertuples(index=False):
                key = tuple(json.loads(row.matched_control_ids))
                if key not in means:
                    mean = normalized(read_ids(raw, list(key), columns, source_positions)).mean(axis=0, dtype='float64')
                    means[key] = mean.astype('float32')
                controls.append(means[key])
            expected_control = np.stack(controls)
            expected = {'X': expected_x, 'control': expected_control,
                        'delta': expected_x - expected_control}
            positions = [new_positions[value] for value in part.sample_id]
            for name, reference in expected.items():
                path = 'X' if name == 'X' else f'layers/{name}'
                left = old[path][start:start + len(part)]
                right = array_rows(new[path], positions)
                if (not np.isfinite(left).all() or not np.isfinite(right).all()
                        or not np.array_equal(left, right) or not np.array_equal(left, reference)):
                    raise ValueError(f'exact mathematical training layer differs: {name}, batch {start}')
                digests[name].update(reference.astype('<f4', copy=False).tobytes(order='C'))
    return {name: value.hexdigest() for name, value in digests.items()}


def accepted_input_seals(report, manifest):
    fields = ('native_counts_sha256','training_manifest_sha256',
              'tenv_training_sha256','worker_sha256')
    if report.get('input_seals') != {field: manifest['provenance'][field] for field in fields}:
        raise ValueError('accepted verification report input seals differ from fitted truth')


def accepted_refits(report, groups, fixture):
    if fixture:
        return
    ids = report.get('independent_refits', [])
    if len(set(ids)) < 3 or not set(ids).issubset(set(groups.group_id)):
        raise ValueError('accepted truth lacks three valid independent refits')
    selected = groups.set_index('group_id').loc[ids]
    masks = [lambda x: x.status.eq('eligible') & x.n_plates.eq(1) & x.n_significant.gt(0),
             lambda x: x.status.eq('eligible') & x.n_plates.gt(1),
             lambda x: x.status.eq('eligible') & x.n_significant.eq(0)]
    if any(mask(groups).any() and not mask(selected).any() for mask in masks):
        raise ValueError('accepted refits do not cover ordinary, multi-plate and tested-empty contrasts')


def equivalence_proof(truth, original_training, training, original_native, native,
                      accepted_verification, *, rscript=DEFAULT_RSCRIPT, fixture=False,
                      runtime=None, config=None, worker=WORKER):
    """Compute complete accepted-statistic, integer-count and exact-expression proof."""
    truth, original_training, training, original_native, native = map(
        Path, (truth, original_training, training, original_native, native))
    manifest = json.loads((truth / 'manifest.json').read_text())
    accepted = json.loads(Path(accepted_verification).read_text())
    if (accepted.get('schema') != VERIFICATION_SCHEMA or accepted.get('dataset') != 'tenv'
            or accepted.get('complete') is not True or manifest.get('schema') != TRUTH_SCHEMA
            or manifest.get('dataset') != 'tenv' or not manifest.get('complete')
            or bool(manifest.get('fixture')) != fixture
            or (not fixture and (manifest.get('group_count') != FULL_TENV_GROUPS
                or accepted.get('treated_samples') != FULL_TENV_SAMPLES))
            or accepted.get('result') != 'PASS' or accepted.get('manifest_sha256') != sha256(truth / 'manifest.json')
            or accepted.get('groups_replayed') != manifest.get('group_count')
            or accepted.get('gene_count') != manifest.get('gene_count')
            or (not fixture and len(set(accepted.get('independent_refits', []))) < 3)):
        raise ValueError('truth successor requires complete accepted truth and independent raw refits')
    accepted_input_seals(accepted, manifest)
    expected = recipe()
    if fixture:
        expected['gene_count'] = manifest['gene_count']
    requested = {**expected, **(config or {})}
    if (requested != expected or manifest.get('gene_count') != expected['gene_count']
            or manifest.get('fdr') != expected['fdr']
            or manifest.get('method') != expected['method_description']
            or manifest.get('independent_filtering') is not False
            or manifest.get('prediction_significance') is not False
            or manifest['provenance'].get('replication_rule') != REPLICATION
            or sha256(worker) != manifest['provenance'].get('worker_sha256')):
        raise ValueError('statistical recipe or worker seal changed; refit is required')
    if runtime is not None and not fixture:
        raise ValueError('production runtime seals must be read from the actual R executable')
    runtime = runtime if runtime is not None else runtime_seal(rscript)
    if any(runtime.get(field) != manifest['provenance'].get(field)
           for field in ('r_version', 'deseq2_version')):
        raise ValueError('R/DESeq2 runtime recipe seal changed')
    # This replays every original raw count digest, membership, BH value, status
    # and ranked list. The acceptance report additionally pins independent refits.
    replay = verify_truth(original_training, original_native, truth, rscript,
                          replay_groups=0, fixture=fixture)
    old_genes, old_groups, old_samples, old_columns = source_plan(
        original_training, original_native, fixture=fixture)
    new_genes, new_groups, new_samples, new_columns = source_plan(training, native, fixture=fixture)
    try:
        pd.testing.assert_frame_equal(old_genes, new_genes, check_dtype=False, check_categorical=False)
        pd.testing.assert_frame_equal(old_groups, new_groups, check_dtype=False, check_categorical=False)
        fields = ['group_id', 'sample_id', 'source_profile_id', 'condition', 'plate']
        order = ['group_id', 'condition', 'source_profile_id']
        pd.testing.assert_frame_equal(old_samples[fields].sort_values(order).reset_index(drop=True),
            new_samples[fields].sort_values(order).reset_index(drop=True), check_dtype=False, check_categorical=False)
    except AssertionError as error:
        raise ValueError('exact axis/contrast/treated/control/design membership differs') from error
    old_obs, _ = training_header(original_training / 'tenv.h5ad')
    new_obs, _ = training_header(training / 'tenv.h5ad')
    identity_digest = sample_identity(old_obs, new_obs)
    if not fixture and (len(old_obs) != FULL_TENV_SAMPLES or len(old_groups) != FULL_TENV_GROUPS):
        raise ValueError('truth successor requires the complete frozen tenv population')
    for directory in (original_training, training):
        training_manifest = json.loads((directory / 'manifest.json').read_text())
        if training_manifest.get('normalization') != NORMALIZATION:
            raise ValueError('normalization recipe differs; exact proof is unavailable')
    reference = pd.read_csv(truth / 'groups.csv', keep_default_na=False).set_index('group_id')
    accepted_refits(accepted, reference.reset_index(), fixture)
    contrasts = []
    with h5py.File(original_native) as old_raw, h5py.File(native) as new_raw:
        new_positions = {value: index for index, value in enumerate(read_elem(new_raw['obs']).index)}
        for key, part in old_samples.groupby('group_id', sort=True):
            part = part.sort_values('source_row')
            left = np.concatenate(list(raw_axis_rows(old_raw['X'], part.source_row.to_numpy(), old_columns)))
            right = read_ids(new_raw, part.source_profile_id.tolist(), new_columns, new_positions)
            digest = hashlib.sha256(left.tobytes(order='C')).hexdigest()
            if not np.array_equal(left, right) or digest != reference.loc[key].counts_sha256:
                raise ValueError(f'exact integer contrast counts differ: {key}')
            row = reference.loc[key]
            contrasts.append({'group_id': key, 'counts_sha256': digest,
                'sample_ids_sha256': row.sample_ids_sha256, 'source_ids_sha256': row.source_ids_sha256,
                'control_ids_sha256': row.control_ids_sha256, 'design': row.design,
                'n_treated': int(row.n_treated), 'n_controls': int(row.n_controls)})
    layers = compare_layers(original_training, training, original_native, old_columns, old_obs, new_obs)
    axis_digest = hashlib.sha256(old_genes.to_csv(index=False, lineterminator='\n').encode()).hexdigest()
    bindings = {'truth_manifest': binding(truth / 'manifest.json'),
                'accepted_verification': binding(accepted_verification),
                'original_training_manifest': binding(original_training / 'manifest.json'),
                'original_training_tenv': binding(original_training / 'tenv.h5ad'),
                'training_manifest': binding(training / 'manifest.json'),
                'training_tenv': binding(training / 'tenv.h5ad'),
                'original_native': binding(original_native), 'native': binding(native),
                'original_native_manifest': producer_binding(original_native),
                'native_manifest': producer_binding(native)}
    if 'rscript_path' in runtime:
        bindings['rscript'] = binding(runtime['rscript_path'])
        if bindings['rscript']['sha256'] != runtime['rscript_sha256']:
            raise ValueError('R executable changed while sealing successor')
    return {'schema': SCHEMA, 'fixture': fixture, 'complete': True,
            'method': 'reuse immutable accepted statistic bytes only after exact mathematical input equivalence',
            'bindings': bindings, 'runtime': runtime, 'worker_sha256': sha256(worker),
            'statistical_recipe': expected,
            'source_truth_files': manifest['files'],
            'proof': {'gene_count': len(old_genes), 'axis_sha256': axis_digest,
                      'sample_count': len(old_obs), 'sample_identity_sha256': identity_digest,
                      'raw_contrasts_compared': len(contrasts), 'contrasts': contrasts,
                      'exact_layer_digests': layers, 'original_complete_replay': replay}}


def create_successor(truth, original_training, training, original_native, native,
                     accepted_verification, output, **options):
    output = Path(output)
    if output.exists():
        raise FileExistsError(output)
    proof = equivalence_proof(truth, original_training, training, original_native, native,
                              accepted_verification, **options)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.truth-successor-', dir=output.parent) as temporary:
        stage = Path(temporary)
        worker = Path(options.get('worker', WORKER))
        (stage / 'deseq2_contrasts.R').write_bytes(worker.read_bytes())
        if sha256(stage / 'deseq2_contrasts.R') != proof['worker_sha256']:
            raise ValueError('worker bytes changed while sealing successor')
        write_json(stage / 'manifest.json', proof)
        stage.rename(output)
    return {'result': 'PASS', 'manifest_sha256': sha256(output / 'manifest.json'),
            'groups_compared': proof['proof']['raw_contrasts_compared'],
            'samples_compared': proof['proof']['sample_count'], 'gene_count': proof['proof']['gene_count']}


def check_successor(envelope, training, *, allow_fixture=False, full_payload=False):
    """Validate immutable proof/input seals before exposing partition truth."""
    root = Path(envelope)
    value = json.loads((root / 'manifest.json').read_text())
    if (value.get('schema') != SCHEMA or not value.get('complete')
            or (value.get('fixture') and not allow_fixture)
            or set(value.get('proof', {}).get('exact_layer_digests', {})) != {'X','control','delta'}):
        raise ValueError('invalid complete truth successor proof')
    for item in value['bindings'].values():
        path = Path(item['path'])
        if path.resolve() != path or sha256(path) != item['sha256']:
            raise ValueError('truth successor input seal changed')
    selected = (Path(training) / 'manifest.json').resolve()
    if str(selected) != value['bindings']['training_manifest']['path']:
        raise ValueError('truth successor is bound to another training release')
    truth = Path(value['bindings']['truth_manifest']['path']).parent
    original = json.loads((truth / 'manifest.json').read_text())
    if (original.get('fixture') != value.get('fixture') or not original.get('complete')
            or (not value.get('fixture') and original.get('gene_count') != 18164)):
        raise ValueError('invalid successor gene family or source fixture status')
    expected = recipe()
    if value.get('fixture'):
        expected['gene_count'] = original['gene_count']
    accepted = json.loads(Path(value['bindings']['accepted_verification']['path']).read_text())
    if (value.get('statistical_recipe') != expected
            or value.get('worker_sha256') != original['provenance']['worker_sha256']
            or sha256(root / 'deseq2_contrasts.R') != value['worker_sha256']
            or any(value['runtime'].get(field) != original['provenance'][field]
                   for field in ('r_version','deseq2_version'))
            or accepted.get('manifest_sha256') != value['bindings']['truth_manifest']['sha256']
            or accepted.get('schema') != VERIFICATION_SCHEMA or accepted.get('dataset') != 'tenv'
            or accepted.get('complete') is not True
            or (not value.get('fixture') and (original['group_count'] != FULL_TENV_GROUPS
                or accepted.get('treated_samples') != FULL_TENV_SAMPLES))
            or accepted.get('result') != 'PASS'
            or accepted.get('groups_replayed') != original['group_count']
            or accepted.get('gene_count') != original['gene_count']
            or (not value.get('fixture') and len(set(accepted.get('independent_refits', []))) < 3)
            or value['proof'].get('gene_count') != original['gene_count']
            or value['proof'].get('raw_contrasts_compared') != original['group_count']
            or len(value['proof'].get('contrasts', [])) != original['group_count']):
        raise ValueError('truth successor recipe/runtime/acceptance proof differs')
    accepted_input_seals(accepted, original)
    checked_file(truth, 'groups.csv', original['files']['groups.csv'])
    groups = pd.read_csv(truth / 'groups.csv', keep_default_na=False).set_index('group_id')
    accepted_refits(accepted, groups.reset_index(), value['fixture'])
    if value['proof'].get('axis_sha256') != sha256(truth / 'genes.csv'):
        raise ValueError('truth successor fixed axis proof differs')
    observed = set()
    for contrast in value['proof']['contrasts']:
        key = contrast['group_id']
        if key in observed or key not in groups.index:
            raise ValueError('truth successor has duplicate/missing contrasts')
        observed.add(key)
        row = groups.loc[key]
        fields = ('counts_sha256','sample_ids_sha256','source_ids_sha256',
                  'control_ids_sha256','design','n_treated','n_controls')
        if any(contrast[field] != row[field] for field in fields):
            raise ValueError('truth successor contrast proof differs')
    if value['source_truth_files'] != original['files']:
        raise ValueError('truth successor statistic seals differ')
    if full_payload:
        for name, digest in original['files'].items():
            checked_file(truth, name, digest)
    return value, truth


class SuccessorTruth(PartitionTruth):
    """Partition reader authorized for an exactly proven successor training release."""
    def __init__(self, envelope, training, assignments, role, *, allow_fixture=False):
        value, original = check_successor(envelope, training, allow_fixture=allow_fixture)
        super().__init__(original, assignments, role, allow_fixture=allow_fixture)
        self.successor_manifest = value
        self.successor_manifest_sha256 = sha256(Path(envelope) / 'manifest.json')


def verify_successor(envelope, training, *, rscript=DEFAULT_RSCRIPT, allow_fixture=False,
                     runtime=None):
    value, truth = check_successor(envelope, training, allow_fixture=allow_fixture, full_payload=True)
    inputs = value['bindings']
    recomputed = equivalence_proof(truth,
        Path(inputs['original_training_manifest']['path']).parent,
        Path(inputs['training_manifest']['path']).parent,
        Path(inputs['original_native']['path']), Path(inputs['native']['path']),
        Path(inputs['accepted_verification']['path']), rscript=rscript,
        fixture=value['fixture'], runtime=runtime, worker=Path(envelope) / 'deseq2_contrasts.R')
    if recomputed != value:
        raise ValueError('truth successor complete equivalence proof replay differs')
    return {'result': 'PASS', 'manifest_sha256': sha256(Path(envelope) / 'manifest.json'),
            'groups_compared': value['proof']['raw_contrasts_compared'],
            'gene_count': value['proof']['gene_count']}
