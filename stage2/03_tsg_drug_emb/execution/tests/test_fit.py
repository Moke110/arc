"""Check native FiLM equivalence and width-independent initialization."""
import os
from pathlib import Path
import sys
import unittest

import torch

EXECUTION = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXECUTION))
from fit import ARMS, WIDTHS, RepresentationAdapter, initialization
from tsg_dose.model import DoseAdapter, fit_transform


SOURCE = Path(os.environ.get('TSG_DRUG_EMB_NATIVE_SOURCE',
                             EXECUTION.parents[2] / 'pert_model/TranSiGen/source'))


@unittest.skipUnless((SOURCE / 'src/model.py').is_file(), 'requires verified native TranSiGen source')
class ModelContractTests(unittest.TestCase):
    def test_all_arms_have_matching_downstream_tensors_and_rng(self):
        torch.set_num_threads(2)
        transform = fit_transform('film', [.2, 1, 10])
        reference = DoseAdapter(SOURCE, 'film', transform).build(42, 'cpu')
        expected = initialization(reference)
        rng = torch.get_rng_state()
        x1 = torch.zeros((2, 5000))
        for arm in ARMS:
            with self.subTest(arm=arm):
                width = WIDTHS[arm]
                model = RepresentationAdapter(SOURCE, transform, width).build(42, 'cpu')
                observed = initialization(model)
                self.assertEqual(observed['downstream_sha256'], expected['downstream_sha256'])
                self.assertEqual(observed['tensors']['feat_embeddings.0.weight']['shape'], [400, width])
                self.assertTrue(torch.equal(torch.get_rng_state(), rng))
                model.eval()
                with torch.no_grad():
                    result = model(x1, torch.zeros((2, width + 1)))[3]
                self.assertEqual(result.shape, (2, 5000))
                self.assertTrue(torch.isfinite(result).all())

    def test_2048_arm_preserves_reference_forward_exactly(self):
        torch.set_num_threads(2)
        transform = fit_transform('film', [.2, 1, 10])
        reference = DoseAdapter(SOURCE, 'film', transform).build(42, 'cpu').eval()
        model = RepresentationAdapter(SOURCE, transform, 2048).build(42, 'cpu').eval()
        rng = torch.get_rng_state()
        controls = torch.ones((2, 5000))
        packed = torch.ones((2, 2049))
        with torch.no_grad():
            expected = reference(controls, packed)
            torch.set_rng_state(rng)
            observed = model(controls, packed)
        self.assertTrue(all(torch.equal(a, b) for a, b in zip(expected, observed, strict=True)))


if __name__ == '__main__':
    unittest.main()
