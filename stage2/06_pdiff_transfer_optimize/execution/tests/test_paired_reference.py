import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))

from paired_reference import score_prediction
from artifacts import write_csv


class PairedReference(unittest.TestCase):
    def test_group_delta_and_equal_direction_and_group_weighting(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'validation.csv'
            write_csv(path, [
                {'case_id': 'a1', 'group_id': 'a', 'drug_key': 'a', 'dose_um': '0.2'},
                {'case_id': 'a2', 'group_id': 'a', 'drug_key': 'a', 'dose_um': '0.2'},
                {'case_id': 'b1', 'group_id': 'b', 'drug_key': 'b', 'dose_um': '0.2'},
                {'case_id': 'c1', 'group_id': 'c', 'drug_key': 'c', 'dose_um': '1'},
                {'case_id': 'd1', 'group_id': 'd', 'drug_key': 'd', 'dose_um': '10'},
            ])
            request = {'role': 'validation', 'n_genes': 3, 'predictor_files': {'validation_cases': str(path)}}
            observed = np.array([[2, -2, 0], [2, -2, 0], [1, 0, 0], [2, 2, 2], [2, -2, 0]],
                                dtype=np.float32)
            predicted = np.array([[5, -5, 0], [-1, 1, 0], [-1, 0, 0], [1, 0, 0], [0, 0, 0]],
                                 dtype=np.float32)
            groups, doses = score_prediction(request, observed, predicted)
            self.assertEqual([row['group_id'] for row in groups], ['a', 'b', 'c', 'd'])
            self.assertEqual(groups[0]['n_wells'], 2)
            self.assertEqual(groups[1]['negative_reason'], 'empty_direction')
            self.assertEqual(groups[2]['pcc_reason'], 'constant_truth')
            self.assertEqual(groups[3]['pcc_reason'], 'constant_prediction')
            self.assertEqual([row['dose_um'] for row in doses], ['0.2', '1.0', '10.0'])
            self.assertAlmostEqual(doses[0]['mean_ap100'], 5 / 6)
            self.assertAlmostEqual(doses[0]['pcc'], 0)
            self.assertEqual(doses[0]['pcc_valid_groups'], 2)
            self.assertEqual(doses[0]['negative_valid_groups'], 1)
            self.assertIsNone(doses[1]['pcc'])
            self.assertEqual(doses[1]['pcc_undefined_reasons'], {'constant_truth': 1})
            self.assertEqual(doses[2]['pcc_undefined_reasons'], {'constant_prediction': 1})
