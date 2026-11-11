"""Verify fold-local ADR0060 inputs and freeze the independent 24-run matrix."""
import argparse
import hashlib
import json
import sys
import h5py
from pathlib import Path
import numpy as np
import pandas as pd

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / '.git').exists())
STEP = REPO / 'stage2/09_xpert_optimize'
sys.path.insert(0, str(STEP / 'execution/01_corpus'))
from build import DEFAULT_NATIVE, normalize_counts, raw_batches

MODELS = ('XPert_G512', 'XPert_G1024', 'XPert_G2048',
          'XPert_MLP512', 'XPert_MLP1024', 'XPert_MLP2048')


def sha256(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def checked(root, manifest, relative):
    path = root / relative
    if path.resolve() != path.absolute() or sha256(path) != manifest['files'][relative]:
        raise ValueError(f'input payload differs: {path}')
    return path


def freeze(inputs, landmark, tenv_release, output, fixture=False, native=DEFAULT_NATIVE):
    landmark_manifest = json.loads((landmark / 'predictor/manifest.json').read_text())
    landmark_sha = sha256(landmark / 'predictor/manifest.json')
    if landmark_manifest['fixture'] != fixture:
        raise ValueError('fixture flag must match landmark release')
    if not fixture:
        accepted = json.loads((STEP / 'execution/landmark_record/corpus.json').read_text())
        if accepted['manifest_sha256']['predictor'] != landmark_sha:
            raise ValueError('landmark release differs from accepted source')
    folds_path = checked(landmark / 'predictor', landmark_manifest, 'tenv_folds.csv')
    folds = pd.read_csv(folds_path, keep_default_na=False)
    pd.testing.assert_frame_equal(folds, pd.read_csv(tenv_release / 'folds.csv', keep_default_na=False))
    if folds.groupby(['fold', 'parent_connectivity']).role.nunique().max() != 1:
        raise ValueError('parent drug crosses roles')
    runs, release_hashes, gates = [], {}, {}
    for fold in (0, 1):
        root = inputs / f'fold_{fold}'
        release = json.loads((root / 'manifest.json').read_text())
        release_hashes[str(fold)] = sha256(root / 'manifest.json')
        if (release['fixture'] != fixture or release['ready_folds'] != [fold]
                or release['landmark_predictor_sha256'] != landmark_sha
                or any(landmark_manifest['source_hashes'].get(name) != digest
                       for name, digest in release['source_hashes'].items())
                or release['rank_recipe']['source_weights'] != {'tahoe': .5, 'tenv': .5}):
            raise ValueError('rerun release provenance differs')
        for relative, digest in release['panels_sha256'].items():
            if sha256(root / relative) != digest:
                raise ValueError('rank/panel payload differs')
        previous = set()
        view_records = {}
        for axis, count in [('G512', 512), ('G1024', 1024), ('G2048', 2048), ('G18164', 18164)]:
            view = root / 'views' / axis
            predictor, evaluator = view / 'predictor', view / 'evaluator'
            manifest = json.loads((predictor / 'manifest.json').read_text())
            evaluation = json.loads((evaluator / 'manifest.json').read_text())
            if (manifest['axis'] != axis or manifest['ready_folds'] != [fold]
                    or evaluation['predictor_manifest_sha256'] != sha256(predictor / 'manifest.json')):
                raise ValueError('view identity differs')
            for side, record in [(predictor, manifest), (evaluator, evaluation)]:
                for relative in record['files']:
                    checked(side, record, relative)
            genes = pd.read_csv(predictor / 'genes.csv')
            panel_path = root / f'fold_{fold}/{axis}.csv'
            if (sha256(panel_path) != manifest['panel_sha256']
                    or sha256(panel_path) != sha256(predictor / 'genes.csv')):
                raise ValueError('gene panel differs')
            if (not fixture and len(genes) != count) or not previous <= set(genes.hgnc_id):
                raise ValueError('axes must be correctly sized and nested')
            previous = set(genes.hgnc_id)
            ranks = pd.read_csv(root / f'fold_{fold}/ranks.csv')
            expected = ranks.iloc[np.sort(ranks.rank_sum.sort_values(kind='stable').index[:len(genes)])]
            if axis != 'G18164' and expected.hgnc_id.tolist() != genes.hgnc_id.tolist():
                raise ValueError('gene panel does not follow equal-source training ranks')
            if axis == 'G18164':
                cm = pd.read_csv(checked(landmark / 'predictor', landmark_manifest, 'cm_genes.csv'))
                if genes.hgnc_id.tolist() != cm.hgnc_id.tolist():
                    raise ValueError('MLP axis must retain the full ordered CM axis')
            for source, roles in [('tahoe', ('train', 'validation')), ('tenv', ('train', 'validation', 'test'))]:
                prefix = source if source == 'tahoe' else f'tenv/fold_{fold}'
                for role in roles:
                    relative = f'{prefix}/{role}.csv'
                    original = checked(landmark / 'predictor', landmark_manifest, relative)
                    if sha256(predictor / relative) != sha256(original):
                        raise ValueError('landmark cohort changed')
                    cases = pd.read_csv(predictor / relative, keep_default_na=False)
                    if source == 'tenv':
                        membership = folds.loc[folds.fold.eq(fold)].set_index('drug_key').role
                        if not cases.drug_key.map(membership).eq(role).all():
                            raise ValueError('tenv drug partition differs')
                    if source == 'tenv':
                        reference_name = 'test_cases' if role == 'test' else role
                        reference = pd.read_csv(tenv_release / f'fold_{fold}/{reference_name}.csv', keep_default_na=False)
                        pd.testing.assert_frame_equal(cases[reference.columns], reference)
                    expression = (evaluator if role == 'test' else predictor) / f'{prefix}/{role}_expression.npy'
                    values = np.load(expression, mmap_mode='r', allow_pickle=False)
                    if values.shape != (len(cases), len(genes)) or values.dtype != np.float32:
                        raise ValueError('normalized expression shape/dtype differs')
                    observed = cases
                    if role == 'test':
                        observed = cases.merge(pd.read_csv(evaluator / f'{prefix}/test_rows.csv'),
                                               on='case_id', validate='one_to_one')
                    positions = np.unique([0, len(cases) - 1])
                    mapping = pd.read_csv(predictor / f'{source}/mapping.csv')
                    rows = observed.source_row.iloc[positions].to_numpy(dtype=int)
                    if np.any(np.diff(rows) <= 0):
                        raise ValueError('native row order differs')
                    with h5py.File(native / f'{source}.h5ad', 'r') as raw:
                        normalized = np.concatenate([normalize_counts(counts)[0]
                            for _, counts in raw_batches(raw['X'], rows, mapping.native_column.to_numpy())])
                    np.testing.assert_allclose(values[positions], normalized, rtol=1e-6, atol=1e-7)
            if list(predictor.rglob('test_expression.npy')) or list(predictor.rglob('test_rows.csv')):
                raise ValueError('test truth leaked into predictor')
            view_records[axis] = {'predictor': str(predictor.absolute()),
                                  'predictor_sha256': sha256(predictor / 'manifest.json'),
                                  'evaluator': str(evaluator.absolute()),
                                  'evaluator_sha256': sha256(evaluator / 'manifest.json'),
                                  'gene_count': len(genes)}
        gates[str(fold)] = {'membership': 'verified', 'parent_drug_exclusion': 'verified',
                            'nested_training_axes': 'verified', 'full_axis_mlp': 'verified',
                            'predictor_evaluator_separation': 'verified', 'normalized_axis_samples': 'verified'}
        for model in MODELS:
            axis = model.removeprefix('XPert_') if '_G' in model else 'G18164'
            for seed in (17, 42):
                runs.append({'id': f'{model}-fold-{fold}-seed-{seed}', 'model': model,
                             'fold': fold, 'seed': seed, **view_records[axis]})
    payload = {'schema': 'step09-round1-rerun-plan-v1', 'fixture': fixture,
               'landmark_predictor_sha256': landmark_sha, 'split_source_sha256': sha256(folds_path),
               'input_manifest_sha256': release_hashes, 'gates': gates, 'runs': runs}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n')
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', type=Path, default=DEFAULT_NATIVE)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--landmark', type=Path, default=STEP / 'output/corpus/tahoe_tenv_bulk_l1000')
    parser.add_argument('--tenv-release', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--fixture', action='store_true')
    args = parser.parse_args()
    freeze(args.inputs, args.landmark, args.tenv_release, args.output, args.fixture, args.native)


if __name__ == '__main__':
    main()
