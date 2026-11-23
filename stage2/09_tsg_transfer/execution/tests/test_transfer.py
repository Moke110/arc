from pathlib import Path
import sys
import unittest

import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'output/api'))

from tsg_transfer.transfer import initialize_from_tahoe, verify_bridge


class Recipient(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.backbone = torch.nn.Linear(8, 4)
        self.dose_embedding = torch.nn.Sequential(
            torch.nn.Linear(1, 32), torch.nn.ReLU(), torch.nn.Linear(32, 32),
        )


class HandoffTest(unittest.TestCase):
    def test_recipient_preserves_source_embedding_at_seen_and_held_doses(self):
        source, recipient = Recipient(), Recipient()
        transform_kind = 'training-log10-dose-population-standardization-v1'
        source_transform = {'kind': transform_kind, 'mean': 0.2, 'std': 1.2}
        target_transform = {'kind': transform_kind, 'mean': -0.31, 'std': 0.4}
        source_before = {key: value.clone() for key, value in source.state_dict().items()}
        doses = (0.05, 0.2, 0.5, 1.0, 5.0, 10.0)
        proof = initialize_from_tahoe(source, recipient, source_transform, target_transform,
                                      doses, held_dose=10.0)
        for key, value in source_before.items():
            torch.testing.assert_close(source.state_dict()[key], value, rtol=0, atol=0)
            if key not in ('dose_embedding.0.weight', 'dose_embedding.0.bias'):
                torch.testing.assert_close(recipient.state_dict()[key], value, rtol=0, atol=0)

        log_doses = torch.log10(torch.tensor(doses)).reshape(-1, 1)
        expected = source.dose_embedding((log_doses - 0.2) / 1.2)
        observed = recipient.dose_embedding((log_doses + 0.31) / 0.4)
        torch.testing.assert_close(observed, expected, rtol=1e-6, atol=1e-6)
        self.assertEqual(proof['held_dose_um'], 10.0)
        self.assertEqual([row['dose_um'] for row in proof['checks']], list(doses))
        verify_bridge(proof, source_transform, target_transform, doses, held_dose=10.0)
        proof['checks'][-1]['max_fraction_of_tolerance'] = 1.1
        with self.assertRaisesRegex(ValueError, 'numerical evidence'):
            verify_bridge(proof, source_transform, target_transform, doses, held_dose=10.0)

    def test_handoff_rejects_missing_held_dose_before_writing_weights(self):
        source, recipient = Recipient(), Recipient()
        transform = {'kind': 'training-log10-dose-population-standardization-v1',
                     'mean': 0.0, 'std': 1.0}
        before = {key: value.clone() for key, value in recipient.state_dict().items()}
        with self.assertRaisesRegex(ValueError, 'held dose'):
            initialize_from_tahoe(source, recipient, transform, transform, (0.2, 1.0),
                                  held_dose=10.0)
        for key, value in before.items():
            torch.testing.assert_close(recipient.state_dict()[key], value, rtol=0, atol=0)


if __name__ == '__main__':
    unittest.main()
