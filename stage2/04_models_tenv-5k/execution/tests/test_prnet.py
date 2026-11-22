"""Approved native prediction and stopping seams (#226)."""
import sys
from pathlib import Path
import unittest
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '.git').exists())
sys.path.insert(0, str(ROOT / 'stage2/04_models_tenv-5k/output/api'))
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
    def test_registered_lower_lr_training_keeps_float32_state_and_predictions(self):
        import numpy as np
        import torch
        from tenv5k.model_prnet import CONFIG, build_model, make_optimizer, make_scheduler, train_epoch, sampled_prediction
        class Data:
            control = np.ones((2, 5000), np.float64)
            expression = np.full((2, 5000), 1.1, np.float64)
            def fingerprints(self, rows): return np.ones((len(rows), 2048), np.float64)
        model = build_model(os.environ['PRNET_SOURCE'], 17)
        optimizer = make_optimizer(model)
        self.assertEqual(optimizer.param_groups[0]['lr'], 1e-4)
        self.assertEqual(CONFIG['method_identity'], 'prnet-native-f32-lr1e4-v1')
        self.assertIsInstance(make_scheduler(optimizer), torch.optim.lr_scheduler.ReduceLROnPlateau)
        self.assertTrue(all(t.dtype == torch.float32 for t in model.state_dict().values() if t.is_floating_point()))
        loss = train_epoch(model, Data(), np.arange(2), optimizer, np.random.default_rng(17), 'cpu')
        self.assertTrue(np.isfinite(loss))
        self.assertTrue(all(p.grad.dtype == torch.float32 for p in model.parameters() if p.grad is not None))
        self.assertTrue(all(s[k].dtype == torch.float32 for s in optimizer.state.values() for k in ('exp_avg', 'exp_avg_sq')))
        prediction = sampled_prediction(model, Data(), np.arange(2), seed=42)
        self.assertEqual(prediction.dtype, np.float32)
        self.assertTrue(np.isfinite(prediction).all())

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
    def test_historical_fourth_fold_is_not_a_registered_training_identity(self):
        import tempfile
        from tenv5k.learner import train
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, 'unregistered run identity'):
                train('/absent-input', Path(directory) / 'candidate', '/absent-source',
                      adapter='prnet', scenario='cold-drug', fold=3, seed=17,
                      smoke_rows=2, max_epochs=1)

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
    def test_best_validation_is_float32_and_historical_identity_is_rejected(self):
        import tempfile
        import numpy as np
        import torch
        from tenv5k.inputs import Inputs
        from tenv5k.learner import replay
        from tenv5k.model_prnet import build_model, sampled_prediction
        source = Path(os.environ['PRNET_SMOKE'])
        saved = torch.load(source / 'checkpoint.pt', map_location='cpu', weights_only=False)
        config = saved['config']
        data = Inputs(os.environ['TENV5K_INPUTS'])
        model = build_model(os.environ['PRNET_SOURCE'], config['seed'], config['device'])
        model.load_state_dict(saved['state_dict'])
        rows = np.asarray(config['validation_order'])
        prediction = sampled_prediction(model, data, rows, seed=config['validation_seed'], device=config['device'])
        self.assertEqual(prediction.dtype, np.float32)
        self.assertTrue(all(t.dtype == torch.float32 for t in saved['state_dict'].values() if t.is_floating_point()))
        error = prediction - np.asarray(data.expression[rows], dtype=np.float32)
        self.assertEqual(saved['validation_mse'], float(np.mean(error * error, dtype=np.float32)))
        with tempfile.TemporaryDirectory() as directory:
            checkpoint = Path(directory) / 'historical.pt'
            saved['config'].pop('method_identity')
            torch.save(saved, checkpoint)
            with self.assertRaisesRegex(ValueError, 'numerical method identity'):
                replay(os.environ['TENV5K_INPUTS'], checkpoint, os.environ['PRNET_SOURCE'], adapter='prnet')

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
                accept_replay(os.environ.get('TENV5K_INPUTS', '/tmp/arc-tenv5k-inputs'), candidate, os.environ['PRNET_SOURCE'], adapter='prnet')

if __name__ == '__main__': unittest.main()
