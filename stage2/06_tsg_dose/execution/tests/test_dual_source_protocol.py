from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'output/api'))
from tsg_round1.protocol import make_splits, score_prediction, summary


class DualSourceProtocolTests(unittest.TestCase):
    def test_fold_audit_covers_rows_and_isolates_parents_and_doses(self):
        obs = pd.DataFrame(
            [(f'p{p}', d, line) for p in range(12)
             for d in (.05, .5, 5.) for line in ('a', 'b')],
            columns=['drug_parent_inchikey', 'dose_um', 'cell_line'],
        )
        splits = make_splits(obs)
        cold = splits[splits.scenario == 'cold-drug']
        self.assertEqual(sorted(cold[cold.partition == 'test'].row), list(range(72)))
        auxiliary = splits[splits.scenario == 'cold-drug-dose']
        for fold, dose in enumerate((.05, .5, 5.)):
            group = auxiliary[auxiliary.fold == fold]
            self.assertEqual(sorted(group.row), list(range(72)))
            parts = {partition: group[group.partition == partition]
                     for partition in ('train', 'validation', 'test')}
            tests = parts['test']
            self.assertEqual(set(tests.dose_um), {dose})
            self.assertEqual(set(group[group.partition.str.startswith('excluded')].partition),
                             {'excluded-held-parent', 'excluded-held-dose'})
            for left, right in [('train', 'validation'), ('train', 'test'),
                                ('validation', 'test')]:
                self.assertTrue(set(parts[left].drug_parent).isdisjoint(parts[right].drug_parent))
            self.assertNotIn(dose, parts['train'].dose_um.values)
            self.assertNotIn(dose, parts['validation'].dose_um.values)
        self.assertEqual(len(auxiliary[auxiliary.partition == 'test']), 24)

    def test_sample_means_weight_repeats_and_record_undefined_counts(self):
        obs = pd.DataFrame({'drug_parent_inchikey': ['a', 'a', 'b', 'c', 'd'],
                            'dose_um': [.2] * 5})
        control = np.tile([.2, .5, 1.], (5, 1))
        truth = control + [[0, 1, 2], [0, 1, 2], [0, 1, 2], [0, 0, 0], [0, 1, 2]]
        predicted = control + [[0, 1, 2], [0, 1, 2], [2, 1, 0], [0, 0, 0], [0, 1, np.nan]]
        result = summary(score_prediction(obs, truth, predicted, control))
        self.assertAlmostEqual(result['pcc'], 1 / 3)
        self.assertEqual(result['undefined_pcc_count'], 2)
        self.assertEqual(result['constant_truth_count'], 1)
        self.assertEqual(result['nonfinite_count'], 1)
