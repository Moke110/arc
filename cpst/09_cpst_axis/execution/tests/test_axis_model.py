"""Paired initialization preserves invariant parameters and the training stream."""

import os
from pathlib import Path
import sys
import unittest

import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from axis_model import AxisAdapter, RESIZED_PARAMETERS
from chem_runs import AlphaAdapter, WIDTHS
from run_measured import DEFAULT_SOURCE

SOURCE = Path(os.environ.get("CARPS_TRANSIGEN_SOURCE", DEFAULT_SOURCE))


class PairedConstruction(unittest.TestCase):
    @unittest.skipUnless((SOURCE / "src/model.py").is_file(),
                         "native TranSiGen source is not materialized")
    def test_gene_width_changes_only_four_projections(self):
        def original():
            return AlphaAdapter(SOURCE, {"kind": "test"}, chemical_width=WIDTHS["UniMol"])
        old = original().build(17, "cpu")
        old_rng = torch.get_rng_state()
        adapted = AxisAdapter(SOURCE, {"kind": "test"}).build(17, "cpu", n_genes=2000)
        before, after = old.state_dict(), adapted.state_dict()
        observed = {name for name in before if before[name].shape != after[name].shape}
        self.assertEqual(observed, RESIZED_PARAMETERS)
        self.assertEqual(after["encoder_x1.0.weight"].shape, (1200, 2000))
        self.assertEqual(after["decoder_delta.4.bias"].shape, (2000,))
        for name in before.keys() - observed:
            self.assertTrue(torch.equal(before[name], after[name]), name)
        self.assertTrue(torch.equal(old_rng, torch.get_rng_state()))
