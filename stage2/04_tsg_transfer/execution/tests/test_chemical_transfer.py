"""Exercise the chemical-only handoff at its explicit module boundary."""
from pathlib import Path
import sys
import unittest

import torch
from torch import nn

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'output/api'))
from tsg_transfer.chemical import initialize_chemical


class ChemicalTransferTest(unittest.TestCase):
    def test_preserves_recipient_dose_and_response_parameters_and_buffers(self):
        def model(value):
            network = nn.Module()
            network.feat_embeddings = nn.Sequential(nn.Linear(2, 3), nn.BatchNorm1d(3))
            network.dose_embedding = nn.Linear(1, 2)
            network.response = nn.Linear(5, 2)
            with torch.no_grad():
                for tensor in network.state_dict().values():
                    tensor.fill_(value)
            return network

        source = model(7)
        recipient = model(3)
        copied = initialize_chemical(source, recipient)
        self.assertIn('feat_embeddings.1.running_mean', copied)
        for name, tensor in recipient.state_dict().items():
            self.assertTrue(torch.all(tensor == (7 if name in copied else 3)), name)
            self.assertEqual(name in copied, name.startswith('feat_embeddings.'))
