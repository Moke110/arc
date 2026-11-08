"""Exercise the shared Tahoe and Transfer initialization seam."""
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
from tahoe_bridge import load_initialization, seal_initialization, state_digest


class FixtureModel:
    def __init__(self):
        self.keys = [('recipient-a', 10.0), ('recipient-b', 0.2), ('recipient-b', 1.0), ('recipient-c', 1.0)]
        self.weights = {
            'native.cov_encoder.drug_encoder.0.weight': torch.tensor([[21.], [32.], [32.], [46.], [0.]]),
            'native.cov_encoder.dose_indices': torch.tensor([3, 1, 2, 2, 0], dtype=torch.int32),
            'native.cov_encoder.dose_encoder.weight': torch.zeros((4, 1)),
            'adapter.weight': torch.tensor([[999.]]),
        }

    def state_dict(self):
        return self.weights

    def load_state_dict(self, weights, strict):
        assert strict and weights.keys() == self.weights.keys()
        self.weights = weights


class BridgeTest(unittest.TestCase):
    def test_tahoe_and_transfer_resolve_same_mapped_state_with_recipient_lookup(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            checkpoint = root / 'best.pt'
            weights = {
                'native.cov_encoder.drug_encoder.0.weight': torch.tensor([[10.], [10.], [10.], [0.]]),
                'native.cov_encoder.dose_indices': torch.tensor([1, 2, 3, 0], dtype=torch.int32),
                'native.cov_encoder.dose_encoder.weight': torch.tensor([[6.], [0.], [10.], [20.]]),
                'adapter.weight': torch.tensor([[45.]]),
            }
            identity = {'dataset': 'tahoe', 'fold': 0, 'seed': 17, 'stage': 'formal'}
            selected = {'identity': identity, 'best_step': 12, 'best_validation_loss': 2.5,
                        'recipe_sha256': 'a' * 64, 'model_state_sha256': state_digest(weights)}
            torch.save({'format': 'step04-ema-prediction-v1', **selected, 'model_state': weights}, checkpoint)
            source = {**selected, 'checkpoint_sha256': hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
                      'path': str(checkpoint), 'doses_um': [0.05, 0.5, 5.0],
                      'dose_indices': [1, 2, 3, 0]}
            request = {'identity': {'dataset': 'tenv', 'fold': 0}, 'n_genes': 18164}
            directory = root / 'bridge'
            directory.mkdir()
            digest = seal_initialization(FixtureModel(), directory, source, request)
            tahoe, transfer = FixtureModel(), FixtureModel()
            first = load_initialization(tahoe, directory, digest, request)
            second = load_initialization(transfer, directory, digest, request)
            self.assertEqual(first, second)
            for model in (tahoe, transfer):
                self.assertEqual(state_digest(model.weights), first['mapped_state_sha256'])
                self.assertEqual(model.weights['native.cov_encoder.drug_encoder.0.weight'].flatten().tolist(),
                                 [21., 32., 32., 46., 0.])
                self.assertEqual(model.weights['native.cov_encoder.dose_indices'].tolist(), [3, 1, 2, 2, 0])
                self.assertEqual(model.weights['adapter.weight'].item(), 45.)
                dose = model.weights['native.cov_encoder.dose_encoder.weight'].flatten()
                self.assertAlmostEqual(dose[1].item(), 6.020599913, places=5)
                self.assertAlmostEqual(dose[2].item(), 13.01029996, places=5)
                self.assertEqual(dose[3].item(), 20.)
                self.assertEqual(dose[0].item(), 6.)
            wrong = {**request, 'identity': {'dataset': 'tenv', 'fold': 1}}
            with self.assertRaisesRegex(ValueError, 'recipient'):
                load_initialization(FixtureModel(), directory, digest, wrong)
            with checkpoint.open('ab') as stream:
                stream.write(b'altered checkpoint')
            with self.assertRaisesRegex(ValueError, 'checkpoint checksum'):
                load_initialization(FixtureModel(), directory, digest, request)


if __name__ == '__main__':
    unittest.main()
