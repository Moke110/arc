"""Tahoe release and predictor CLI through a compact CSR fixture."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

CLI = Path(__file__).resolve().parents[1] / 'script/train.py'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class TahoeArtifactsTest(unittest.TestCase):
    def fixture(self, root):
        step02 = root / 'stage2/02_dataset_corpus'
        corpus = step02 / 'output/corpus/v2'
        (corpus / 'tables').mkdir(parents=True)
        genes = pd.DataFrame({'symbol': ['a', 'b', 'FAM153B', *[f'g{i}' for i in range(18161)]],
                              'ensembl_id': ['a', 'b', 'ENSG00000182230', *[f'id{i}' for i in range(18161)]]})
        genes.to_csv(corpus / 'tables/genes.csv', index=False)
        reference = root / 'stage2/03_modes_tenv/output/input'
        reference.mkdir(parents=True)
        shutil.copyfile(corpus / 'tables/genes.csv', reference / 'genes.csv')
        manifest = reference / 'manifest.json'
        manifest.write_text(json.dumps({'files': {'genes.csv': sha(reference / 'genes.csv')}}))
        pin = root / 'stage2/03_modes_tenv/execution/input_record'
        pin.mkdir(parents=True)
        (pin / 'verification.json').write_text(json.dumps({'manifest_sha256': sha(manifest)}))
        names = [f'd{i}' for i in range(9)] + ['unknown']
        records, values = [], []
        for plate, line, well, drug, dose, value in [
            ('p1', 'A', 'c1', 'DMSO_TF', 0, 2),
            ('p1', 'A', 'c2', 'DMSO_TF', 0, 4),
            ('p1', 'B', 'c1', 'DMSO_TF', 0, 1),
            ('p2', 'A', 'c1', 'DMSO_TF', 0, 0.5),
            ('p1', 'A', 'other', 'Other', 0, 100),
            ('p2', 'C', 'unpaired', 'd0', 5, 10),
            ('p1', 'B', 'repeat', 'd0', 5, 7),
            ('p2', 'A', 'repeat', 'd0', 0.5, 8),
            *[('p1', 'A', str(i), name, 5, 5) for i, name in enumerate(names)],
        ]:
            records.append({'profile_id': f'{plate}:{well}:{line}', 'plate': plate, 'well': well,
                            'cell_line_id': line, 'drug': drug, 'dose_um': float(dose),
                            'is_control': dose == 0})
            values.append((0.5, value))
        obs = pd.DataFrame(records).set_index('profile_id')
        matrix = csr_matrix(np.asarray(values, dtype='float32'))
        native = csr_matrix((matrix.data, matrix.indices, matrix.indptr), shape=(len(obs), 18163))
        var = genes.drop(index=2).iloc[[1, 0, *range(2, 18163)]].set_index('symbol')
        ad.AnnData(X=native, obs=obs, var=var, layers={'log1p_cp10k': native}).write_h5ad(corpus / 'tahoe.h5ad')
        build = step02 / 'execution/build_record'
        build.mkdir(parents=True)
        (build / 'tahoe_corpus_move.json').write_text(json.dumps({'file_sha256': sha(corpus / 'tahoe.h5ad')}))
        source = root / 'pert_data/Tahoe-100M'
        drugs = source / 'source/metadata/drug_metadata.parquet'
        drugs.parent.mkdir(parents=True)
        smiles = ['C[C@H](O)F', 'C[C@@H](O)F', 'CCC', 'CCCC', 'CCCCC', 'O', 'N', 'CCO', 'CCN']
        pd.DataFrame({'drug': names[:-1], 'canonical_smiles': smiles}).to_parquet(drugs)
        (source / 'manifest.sha256').write_text(f'{sha(drugs)}  metadata/drug_metadata.parquet\n')

    def command(self, root, *args, expected_error=None):
        completed = subprocess.run([sys.executable, str(CLI), '--runtime-root', str(root), *map(str, args)],
                                   capture_output=True, text=True,
                                   env={**os.environ, 'OPENBLAS_NUM_THREADS': '2'})
        if expected_error:
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn(expected_error, completed.stderr)
        else:
            self.assertEqual(completed.returncode, 0, completed.stderr)
        return completed

    def test_release_matches_controls_and_freezes_case_roles(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.fixture(root)
            release, duplicate = root / 'release', root / 'duplicate'
            self.command(root, 'build-tahoe', release)
            self.command(root, 'verify-tahoe', release)
            self.command(root, 'build-tahoe', duplicate)
            manifest_path = release / 'manifest.json'
            self.assertEqual(sha(manifest_path), sha(duplicate / 'manifest.json'))
            manifest = json.loads(manifest_path.read_text())
            self.assertEqual(manifest['gene_projection']['zero_gene'], 'ENSG00000182230')
            exclusions = pd.read_csv(release / 'exclusions.csv').set_index('case_id')
            self.assertEqual(exclusions.loc['p1:9:A', 'reason'], 'missing_structure')
            self.assertEqual(exclusions.loc['p2:unpaired:C', 'reason'], 'missing_matched_control')
            memberships = pd.read_csv(release / 'folds.csv')
            for fold in range(3):
                base = release / f'fold_{fold}'
                roles = memberships[memberships.fold.eq(fold)].set_index('drug_key').role
                self.assertEqual(roles['d0'], roles['d1'])
                role = roles['d0']
                test = role == 'test'
                cases = pd.read_csv(base / ('test_cases.csv' if test else f'{role}.csv'))
                expression = np.load(base / ('evaluator/truth.npy' if test else f'{role}_expression.npy'))
                for line_case, expected in [('p1:0:A', 2), ('p1:repeat:B', 6), ('p2:repeat:A', 7.5)]:
                    index = cases.index[cases.case_id.eq(line_case)][0]
                    if test:
                        groups = pd.read_csv(base / 'evaluator/groups.csv')
                        index = groups.index[groups.group_id.eq(cases.loc[index, 'group_id'])][0]
                    np.testing.assert_array_equal(expression[index, :3], [expected, 0, 0])
            shutil.rmtree(release / 'fold_0/evaluator')
            request = root / 'request.json'
            self.command(root, 'prepare-tahoe', release, request, '--fold', 0)
            self.assertNotIn('evaluator', request.read_text())
            self.command(root, 'verify-request', request)


    def test_verifier_rejects_rehashed_parent_leak_bad_pairing_and_reordered_genes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.fixture(root)
            release = root / 'release'
            self.command(root, 'build-tahoe', release)
            manifest_path = release / 'manifest.json'
            manifest = json.loads(manifest_path.read_text())
            split = pd.read_csv(release / 'folds.csv')
            first = split.fold.eq(0) & split.drug_key.eq('d1')
            split.loc[first, 'role'] = 'test' if split.loc[first, 'role'].iat[0] == 'train' else 'train'
            wells = pd.read_csv(release / 'wells.csv')
            wells.loc[wells.case_id.eq('p1:0:A'), 'control_index'] = 1
            genes = pd.read_csv(release / 'genes.csv')
            for name, frame, error in [
                ('folds.csv', split, 'parent-connectivity block membership'),
                ('wells.csv', wells, 'matched treatment metadata'),
                ('genes.csv', genes.iloc[::-1], 'frozen ordered CM axis'),
            ]:
                with self.subTest(name=name):
                    path = release / name
                    before = path.read_bytes()
                    old_hash = manifest['files'][name]
                    frame.to_csv(path, index=False)
                    manifest['files'][name] = sha(path)
                    manifest_path.write_text(json.dumps(manifest))
                    self.command(root, 'verify-tahoe', release, expected_error=error)
                    path.write_bytes(before)
                    manifest['files'][name] = old_hash
                    manifest_path.write_text(json.dumps(manifest))

    def test_predictor_rejects_rehashed_gene_order_and_in_range_control(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.fixture(root)
            release = root / 'release'
            self.command(root, 'build-tahoe', release)
            request_path = root / 'request.json'
            self.command(root, 'prepare-tahoe', release, request_path, '--fold', 0)
            request = json.loads(request_path.read_text())
            manifest_path = release / 'manifest.json'
            manifest = json.loads(manifest_path.read_text())
            cases_key = next(
                f'{role}_cases' for role in ('train', 'validation', 'test')
                if 'p1:0:A' in pd.read_csv(request['predictor_files'][f'{role}_cases']).case_id.values
            )
            for key, error in [('genes', 'frozen ordered CM axis'),
                               (cases_key, 'case membership or control pairing')]:
                with self.subTest(key=key):
                    consumer = request['input']
                    relative = consumer['predictor_files'][key]
                    path = release / relative
                    before = path.read_bytes()
                    frame = pd.read_csv(path)
                    if key == 'genes':
                        frame = frame.iloc[::-1]
                    else:
                        frame.loc[frame.case_id.eq('p1:0:A'), 'control_index'] = 1
                    frame.to_csv(path, index=False)
                    manifest['files'][relative] = sha(path)
                    manifest_path.write_text(json.dumps(manifest))
                    consumer['source_input_sha256'] = request['identity']['source_input_sha256'] = sha(manifest_path)
                    consumer['file_sha256'][key] = request['file_sha256'][key] = sha(path)
                    if key == 'genes':
                        consumer['gene_sha256'] = request['identity']['gene_sha256'] = sha(path)
                    request['identity']['input_sha256'] = hashlib.sha256(
                        json.dumps(consumer, sort_keys=True, separators=(',', ':')).encode()
                    ).hexdigest()
                    forged = root / 'forged.json'
                    forged.write_text(json.dumps(request))
                    self.command(root, 'prepare-tahoe', release, root / 'bad.json', '--fold', 0,
                                 expected_error=error)
                    self.command(root, 'verify-request', forged, expected_error=error)
                    self.command(root, 'fit', forged, root / 'unused-recipe', root / 'unused-output',
                                 expected_error=error)
                    path.write_bytes(before)
                    manifest['files'][relative] = consumer['file_sha256'][key] = request['file_sha256'][key] = sha(path)
                    consumer['gene_sha256'] = request['identity']['gene_sha256'] = manifest['files']['genes.csv']

    def test_predictor_rejects_parent_leak_from_distinct_drugs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.fixture(root)
            release = root / 'release'
            self.command(root, 'build-tahoe', release)
            manifest_path = release / 'manifest.json'
            manifest = json.loads(manifest_path.read_text())
            split = pd.read_csv(release / 'folds.csv')
            role = split[split.fold.eq(0) & split.drug_key.eq('d1')].role.iat[0]
            target_role = 'test' if role == 'train' else 'train'
            original = 'fold_0/test_cases.csv' if role == 'test' else f'fold_0/{role}.csv'
            original_cases = pd.read_csv(release / original)
            original_cases[original_cases.drug_key.ne('d0')].to_csv(release / original, index=False)
            manifest['files'][original] = sha(release / original)
            relative = 'fold_0/test_cases.csv' if target_role == 'test' else 'fold_0/train.csv'
            path = release / relative
            cases = pd.read_csv(path)
            cases.loc[0, 'drug_key'] = 'd0'
            cases.to_csv(path, index=False)
            manifest['files'][relative] = sha(path)
            manifest_path.write_text(json.dumps(manifest))
            self.command(root, 'prepare-tahoe', release, root / 'request.json', '--fold', 0,
                         expected_error='case membership or control pairing')


if __name__ == '__main__':
    unittest.main()
