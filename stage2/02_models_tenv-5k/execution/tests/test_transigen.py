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

    def test_full_native_training_and_sampling_preserve_declared_precision(self):
        import torch
        from types import SimpleNamespace
        from tenv5k.model_transigen import Adapter, CONFIG
        root = next(p for p in Path(__file__).resolve().parents if (p / '.git').exists())
        source = Path(os.environ.get('TRANSIGEN_SOURCE', str(root / 'pert_model/TranSiGen/source')))
        if not (source / 'src/model.py').is_file():
            self.skipTest('materialize native source or set TRANSIGEN_SOURCE')
        adapter = Adapter(source)
        model = adapter.build(17, 'cpu')
        self.assertEqual(CONFIG['method_identity'], 'transigen-native-f64-v1')
        # Native initialization is the independent reference (float32); casting
        # must add no random draws and must retain every initialized value.
        after_build = torch.get_rng_state().clone()
        torch.manual_seed(17)
        native = adapter.native(n_genes=5000, n_latent=100, n_en_hidden=[1200],
            n_de_hidden=[800], features_dim=2048, features_embed_dim=[400],
            init_w=True, beta=0.1, dropout=0.1, random_seed=17, device=torch.device('cpu'))
        self.assertTrue(torch.equal(after_build, torch.get_rng_state()))
        for name, value in native.state_dict().items():
            torch.testing.assert_close(model.state_dict()[name], value.double() if value.is_floating_point() else value, rtol=0, atol=0)
        del native

        self.assertTrue(all(v.dtype == torch.float64 for v in model.state_dict().values() if v.is_floating_point()))
        # Compare the sampling boundary to native CPU float32 draws, including
        # the subsequent RNG state; direct float64 draws consume another stream.
        mu = torch.linspace(-1, 1, 300, dtype=torch.float64).reshape(3, 100)
        logvar = torch.zeros_like(mu)
        torch.manual_seed(42)
        noise = torch.FloatTensor(mu.size()).normal_()
        expected_state = torch.get_rng_state().clone()
        expected = noise.double() + mu
        torch.manual_seed(42)
        actual = model.sample_latent(mu, logvar)
        self.assertEqual(actual.dtype, torch.float64)
        torch.testing.assert_close(actual, expected, rtol=0, atol=0)
        self.assertTrue(torch.equal(torch.get_rng_state(), expected_state))
        class Data(SimpleNamespace):
            def fingerprints(self, rows): return np.ones((len(rows), 2048), np.float32)
        data = Data(control=np.ones((2,5000), np.float32), expression=np.full((2,5000), 1.1, np.float32))
        optimizer = adapter.optimizer(model)
        loss = adapter.train_epoch(model, data, np.arange(2), optimizer, np.random.default_rng(17))
        self.assertTrue(np.isfinite(loss))
        self.assertTrue(all(p.grad.dtype == torch.float64 and torch.isfinite(p.grad).all() for p in model.parameters() if p.grad is not None))
        self.assertTrue(all(s[k].dtype == torch.float64 for s in optimizer.state.values() for k in ('exp_avg', 'exp_avg_sq')))
        prediction = adapter.predict(model, data, np.arange(2))
        self.assertEqual(prediction.dtype, np.float64)
        self.assertTrue(np.isfinite(prediction).all())

class PrecisionReplayTests(unittest.TestCase):
    def test_serialized_values_are_checked_after_declared_conversion(self):
        from tenv5k.learner import serialize_prediction
        value = np.array([[1e100]], np.float64)
        self.assertEqual(serialize_prediction(value, {'prediction_dtype':'float64'}).dtype, np.float64)
        with self.assertRaisesRegex(ValueError, 'after serialization'):
            serialize_prediction(value, {'prediction_dtype':'float32'})

    @unittest.skipUnless(os.environ.get('TRANSIGEN_SMOKE'), 'set TRANSIGEN_SMOKE after native smoke')
    def test_replay_rejects_old_precision_even_with_same_source(self):
        import torch
        import tempfile
        from tenv5k.learner import replay
        saved = torch.load(Path(os.environ['TRANSIGEN_SMOKE'])/'checkpoint.pt', map_location='cpu', weights_only=False)
        with tempfile.TemporaryDirectory() as directory:
            for key in ('method_identity', 'compute_dtype', 'prediction_dtype', 'initialization_dtype', 'noise_dtype'):
                original = saved['config'].pop(key)
                checkpoint = Path(directory)/'old.pt'
                torch.save(saved, checkpoint)
                with self.assertRaisesRegex(ValueError, 'numerical method identity'):
                    replay('/tmp/arc-tenv5k-inputs', checkpoint, os.environ['TRANSIGEN_SOURCE'], adapter='transigen')
                saved['config'][key] = original
            saved['state_dict'] = {k: v.float() if v.is_floating_point() else v for k,v in saved['state_dict'].items()}
            torch.save(saved, checkpoint)
            with self.assertRaisesRegex(ValueError, 'state precision'):
                replay('/tmp/arc-tenv5k-inputs', checkpoint, os.environ['TRANSIGEN_SOURCE'], adapter='transigen')
