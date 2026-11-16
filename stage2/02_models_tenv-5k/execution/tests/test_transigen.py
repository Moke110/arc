"""Approved model prediction seam: native sampling without treated-label access."""
import sys
import os
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'output/api'))

class TranSiGenTests(unittest.TestCase):
    def test_prediction_ignores_treated_labels_and_replays_native_sampling(self):
        from tenv5k.model_transigen import Adapter, sampled_prediction
        import torch
        from types import SimpleNamespace
        root = next(p for p in Path(__file__).resolve().parents if (p / '.git').exists())
        source = Path(os.environ.get('TRANSIGEN_SOURCE', str(root / 'pert_model/TranSiGen/source')))
        if not (source / 'src/model.py').is_file():
            self.skipTest('materialize native source or set TRANSIGEN_SOURCE')
        adapter = Adapter(source)
        model = adapter.build(17, 'cpu')
        class Inputs(SimpleNamespace):
            def fingerprints(self, rows):
                return np.zeros((len(rows), 2048), np.float32)
            @property
            def expression(self):
                raise AssertionError('prediction accessed treated labels')
        inputs = Inputs(control=np.ones((3,5000), np.float32))
        model._tenv_adapter = adapter
        before = torch.random.get_rng_state().clone()
        first = sampled_prediction(model, inputs, np.arange(3), seed=42, batch_size=2)
        self.assertTrue(torch.equal(before, torch.random.get_rng_state()))
        torch.manual_seed(42)
        replay = adapter.predict(model, inputs, np.arange(3), 2)
        np.testing.assert_array_equal(first, replay)
        torch.manual_seed(43)
        other = adapter.predict(model, inputs, np.arange(3), 2)
        self.assertEqual(first.shape, (3,5000))
        self.assertTrue(np.isfinite(first).all())
        self.assertFalse(np.array_equal(first, other))
