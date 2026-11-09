import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))

from baseline import inference_request, score_prediction
from artifacts import sha256, write_csv
import numpy as np


class ValidationBaseline(unittest.TestCase):
    def test_directional_scores_use_equal_group_weight_and_separate_truth(self):
        cases = [
            {'case_id': 'a1', 'group_id': 'a', 'drug_key': 'a', 'dose_um': '0.2', 'control_index': '0'},
            {'case_id': 'a2', 'group_id': 'a', 'drug_key': 'a', 'dose_um': '0.2', 'control_index': '0'},
            {'case_id': 'b1', 'group_id': 'b', 'drug_key': 'b', 'dose_um': '1', 'control_index': '0'},
        ]
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            files = {}
            for key in ('train_cases', 'validation_cases', 'test_cases', 'genes', 'drugs', 'morgan', 'controls'):
                path = root / key
                if key == 'validation_cases':
                    write_csv(path, cases)
                else:
                    path.touch()
                files[key] = str(path)
            files['validation_expression'] = str(root / 'unavailable')
            request = inference_request({
                'identity': {'gene_sha256': sha256(files['genes'])},
                'target_kind': 'delta', 'n_genes': 2, 'predictor_files': files,
                'file_sha256': {key: sha256(path) for key, path in files.items()
                                if key != 'validation_expression'},
            })
            self.assertNotIn('validation_expression', request['predictor_files'])
            observed = np.array([[2, -2], [2, -2], [2, 0]], dtype=np.float32)
            predicted = np.array([[9, 0], [-5, 0], [0, 1]], dtype=np.float32)
            group, summary = score_prediction(request, observed, predicted)
            self.assertEqual([row['group_id'] for row in group], ['a', 'b'])
            self.assertEqual(summary[0]['positive_ap100'], 0.75)
            self.assertEqual(summary[0]['negative_ap100'], 1.0)
            self.assertEqual(summary[0]['mean_ap100'], 0.875)
            self.assertEqual(summary[0]['negative_valid_groups'], 1)
            self.assertEqual(group[1]['negative_reason'], 'empty_direction')
            self.assertEqual([row['group_count'] for row in summary], [2, 1, 1])

    def test_request_rejects_reordered_gene_case_or_group_evidence(self):
        from baseline import verify_inference

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            contents = {
                'genes': 'symbol\nA\nB\n',
                'validation_cases': 'case_id,group_id\na,g1\nb,g2\n',
                'train_cases': '', 'test_cases': '', 'controls': '', 'morgan': '', 'drugs': '',
            }
            files, hashes = {}, {}
            for key, content in contents.items():
                path = root / key
                path.write_text(content)
                files[key], hashes[key] = str(path), sha256(path)
            request = inference_request({
                'identity': {'gene_sha256': hashes['genes']}, 'target_kind': 'delta',
                'n_genes': 2, 'predictor_files': files, 'file_sha256': hashes,
            })
            self.assertEqual([row['case_id'] for row in verify_inference(request)], ['a', 'b'])
            for key, changed in (
                ('genes', 'symbol\nB\nA\n'),
                ('validation_cases', 'case_id,group_id\nb,g2\na,g1\n'),
                ('validation_cases', 'case_id,group_id\na,g2\nb,g1\n'),
            ):
                with self.subTest(key=key, changed=changed):
                    (root / key).write_text(changed)
                    with self.assertRaisesRegex(ValueError, 'hash mismatch'):
                        verify_inference(request)
                    (root / key).write_text(contents[key])
