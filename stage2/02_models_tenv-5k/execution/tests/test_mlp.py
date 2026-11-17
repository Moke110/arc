"""Approved prediction seam: chemical-only treated expression and exact replay."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'output/api'))


class MLPTests(unittest.TestCase):
    def test_prediction_uses_only_morgan_and_preserves_row_order(self):
        from tenv5k.model_mlp import build_model, sampled_prediction
        class ChemicalInputs:
            def fingerprints(self, rows):
                result = np.zeros((len(rows), 2048), np.float32)
                result[:, 0] = np.asarray(rows) % 2
                return result
            @property
            def expression(self):
                raise AssertionError('prediction accessed treated labels')
            @property
            def control(self):
                raise AssertionError('chemical baseline accessed control')
        model = build_model(None, 17, 'cpu')
        data = ChemicalInputs()
        result = sampled_prediction(model, data, np.array([1, 0, 1]), seed=42)
        self.assertEqual(result.shape, (3, 5000))
        self.assertTrue(np.isfinite(result).all())
        np.testing.assert_array_equal(result[0], result[2])
        self.assertFalse(np.array_equal(result[0], result[1]))
        np.testing.assert_array_equal(result, sampled_prediction(model, data, np.array([1, 0, 1]), seed=2026))
        self.assertTrue((result < 0).any(), 'treated prediction must retain linear negative outputs')

    def test_training_learns_expression_and_retains_singleton_batch(self):
        from tenv5k.model_mlp import CONFIG, build_model, make_optimizer, train_epoch, sampled_prediction
        import torch
        class Inputs:
            expression = np.ones((65, 5000), np.float32)
            def fingerprints(self, rows):
                return np.ones((len(rows), 2048), np.float32)
        model = build_model(None, 17, 'cpu')
        data, rows = Inputs(), np.arange(65)
        optimizer = make_optimizer(model)
        self.assertEqual(optimizer.param_groups[0]['lr'], 1e-4)
        self.assertEqual(CONFIG['method_identity'], 'mlp-morgan-project-f32-lr1e4-v1')
        self.assertTrue(all(t.dtype == torch.float32 for t in model.state_dict().values() if t.is_floating_point()))
        before = float(np.mean((sampled_prediction(model, data, rows) - 1) ** 2))
        loss = train_epoch(model, data, rows, optimizer, np.random.default_rng(17), 'cpu')
        after = float(np.mean((sampled_prediction(model, data, rows) - 1) ** 2))
        self.assertTrue(np.isfinite(loss))
        self.assertLess(after, before)
        self.assertEqual(sampled_prediction(model, data, rows).dtype, np.float32)
        self.assertTrue(all(p.grad.dtype == torch.float32 for p in model.parameters() if p.grad is not None))
        self.assertTrue(all(s[k].dtype == torch.float32 for s in optimizer.state.values() for k in ('exp_avg', 'exp_avg_sq')))
