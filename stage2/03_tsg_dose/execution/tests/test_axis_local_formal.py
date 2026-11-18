"""Formal launcher boundaries stop on failed origins and never retrain completed phases."""
import importlib.util
from pathlib import Path
import unittest

path = Path(__file__).resolve().parents[1] / 'axis_local/formal.py'
spec = importlib.util.spec_from_file_location('axis_local_formal', path)
formal = importlib.util.module_from_spec(spec)
spec.loader.exec_module(formal)


class PhaseTests(unittest.TestCase):
    def test_resume_runs_only_missing_phase(self):
        self.assertEqual(formal.next_action(None, 'validation'), 'train')
        self.assertEqual(formal.next_action('validation-trained-awaiting-replay', 'validation'), 'replay')
        self.assertIsNone(formal.next_action('validation-complete', 'validation'))
        self.assertEqual(formal.next_action('validation-complete', 'test'), 'evaluate')
        self.assertEqual(formal.next_action('test-awaiting-replay', 'test'), 'replay')
        self.assertIsNone(formal.next_action('complete', 'test'))
        for status in ('failed', 'interrupted'):
            with self.assertRaisesRegex(ValueError, status):
                formal.next_action(status, 'validation')


class EstimandTests(unittest.TestCase):
    def test_sample_pool_and_equal_fold_spread_remain_distinct(self):
        import sys
        import pandas as pd
        sys.path.insert(0, str(path.parent))
        from formal_report import summarize_folds

        def scores(values):
            return pd.DataFrame(dict(pcc=values, spearman=values, rmse=values, r2=values,
                                     constant_truth=False, constant_prediction=False,
                                     nonfinite=False))
        folds = {0: scores([0.0, 0.0, 0.0]), 1: scores([0.9]), 2: scores([0.6])}
        pooled, dots, spread = summarize_folds(folds)
        self.assertAlmostEqual(pooled['pcc'], 0.3)
        self.assertEqual([row['pcc'] for row in dots], [0.0, 0.9, 0.6])
        pcc = next(row for row in spread if row['metric'] == 'pcc')
        self.assertAlmostEqual(pcc['mean'], 0.5)
        self.assertAlmostEqual(pcc['sd'], (0.42 / 2) ** 0.5)
        self.assertEqual(pcc['n_defined_folds'], 3)
