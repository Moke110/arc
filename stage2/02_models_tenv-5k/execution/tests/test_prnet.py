"""Approved native prediction and stopping seams (#226)."""
import sys
from pathlib import Path
import unittest
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '.git').exists())
sys.path.insert(0, str(ROOT / 'stage2/02_models_tenv-5k/output/api'))
from tenv5k.learning import EarlyStopping

class StopTests(unittest.TestCase):
    def test_equal_is_not_improvement_and_fiftieth_stale_stops(self):
        stop = EarlyStopping()
        self.assertTrue(stop.observe(3.0, 1))
        for epoch in range(2, 51):
            self.assertFalse(stop.observe(3.0, epoch))
            self.assertFalse(stop.stopped)
        self.assertFalse(stop.observe(3.0, 51))
        self.assertTrue(stop.stopped)
        self.assertEqual(stop.best_epoch, 1)

    def test_strict_improvement_resets_patience(self):
        stop = EarlyStopping()
        stop.observe(4, 1)
        for epoch in range(2, 50): stop.observe(5, epoch)
        self.assertTrue(stop.observe(3, 50))
        self.assertEqual(stop.stale, 0)
        self.assertFalse(stop.stopped)


class BudgetTests(unittest.TestCase):
    def test_cap_and_nonfinite_fail_closed(self):
        stop = EarlyStopping(max_epochs=3)
        stop.observe(3, 1)
        stop.observe(2, 2)
        stop.observe(1, 3)
        self.assertTrue(stop.stopped)
        with self.assertRaises(ValueError): EarlyStopping().observe(float('nan'), 1)

# Native imports are opt-in so these tests remain separately resource-capped.
import os
@unittest.skipUnless(os.environ.get('PRNET_SOURCE'), 'set PRNET_SOURCE under a resource cap')
class NativePredictionTests(unittest.TestCase):
    def test_sampled_treated_expression_replays_without_target_labels(self):
        import numpy as np
        import torch
        from tenv5k.model_prnet import build_model, sampled_prediction
        class Background:
            control = np.ones((3, 5000), dtype=np.float32)
            def fingerprints(self, rows): return np.zeros((len(rows), 2048), np.float32)
        model = build_model(os.environ['PRNET_SOURCE'], 17)
        rows = np.arange(3)
        before = torch.random.get_rng_state().clone()
        first = sampled_prediction(model, Background(), rows, seed=123)
        np.testing.assert_array_equal(first, sampled_prediction(model, Background(), rows, seed=123))
        self.assertTrue(torch.equal(before, torch.random.get_rng_state()))
        other = sampled_prediction(model, Background(), rows, seed=124)
        self.assertFalse(np.array_equal(first, other))
        self.assertLess(first.min(), 0)  # Native Normal sample, unlike nonnegative decoder mean.
        self.assertEqual(first.shape, (3, 5000))


class FailedCandidateTests(unittest.TestCase):
    def test_invalid_formal_budget_fails_without_prediction_or_manifest(self):
        import json
        import tempfile
        from tenv5k.learner import train
        with tempfile.TemporaryDirectory() as directory:
            candidate = Path(directory) / 'candidate'
            with self.assertRaises(ValueError):
                train('/absent-input', candidate, '/absent-source', adapter='prnet',
                      scenario='cold-drug', fold=0, seed=17, max_epochs=2)
            self.assertEqual(json.loads((candidate / 'terminal.json').read_text())['status'], 'failed')
            self.assertFalse((candidate / 'prediction.npy').exists())
            self.assertFalse((candidate / 'manifest.json').exists())


@unittest.skipUnless(os.environ.get('PRNET_SMOKE'), 'set PRNET_SMOKE after real-input smoke')
class SmokeAcceptanceTests(unittest.TestCase):
    def test_best_checkpoint_is_recorded_curve_minimum_and_config_tampering_rejected(self):
        import json
        import shutil
        import tempfile
        import pandas as pd
        import torch
        from tenv5k.learner import accept_replay
        source = Path(os.environ['PRNET_SMOKE'])
        curve = pd.read_csv(source / 'learning_curve.csv')
        checkpoint = torch.load(source / 'checkpoint.pt', map_location='cpu', weights_only=False)
        self.assertEqual(checkpoint['best_epoch'], int(curve.iloc[curve.validation_mse.argmin()].epoch))
        self.assertAlmostEqual(checkpoint['validation_mse'], curve.validation_mse.min())
        with tempfile.TemporaryDirectory() as directory:
            candidate = Path(directory) / 'corrupt'
            shutil.copytree(source, candidate)
            config = json.loads((candidate / 'config.json').read_text())
            config['seed'] = 42
            (candidate / 'config.json').write_text(json.dumps(config))
            with self.assertRaisesRegex(ValueError, 'configuration mismatch'):
                accept_replay('/tmp/arc-tenv5k-inputs', candidate, os.environ['PRNET_SOURCE'], adapter='prnet')

if __name__ == '__main__': unittest.main()
