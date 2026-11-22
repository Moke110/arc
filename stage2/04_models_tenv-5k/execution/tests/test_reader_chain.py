"""Real public Bases reader to train-only Mean and matched-control scoring."""
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / '01_bases/execution/tests'))
sys.path.insert(0, str(ROOT / '01_bases/output/api'))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'output/api'))
import view_fixtures
from bases import view
from tenv5k.protocol import mean_prediction, score_prediction

class ReaderChainTest(unittest.TestCase):
    def test_real_reader_mean_uses_original_control_rows_and_frozen_delta(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture = view_fixtures.normalization_fixture(root)
            release = root / 'view'
            view.build(release, corpus_dir=fixture['corpus'], axis_dir=fixture['axis'], size=2, fixture=True)
            with view.open_view(release) as opened:
                source = opened.source('tenv')
                obs = source.obs.copy()
                # A treated training well and an independently scored row. The
                # two original DMSO controls enter exactly once each.
                x = source.X()
                controls = {0: source.matched_control_rows(0).tolist()}
                prediction = mean_prediction(x, obs, controls, [0], [3], 'Mean', 'cold-drug')
                np.testing.assert_allclose(source.delta(), x - source.control(), atol=0, rtol=0)
                wells, parents = score_prediction(obs.iloc[[3]], x[[3]], prediction, source.control(3,4), prediction_ids=obs.index[[3]])
                self.assertEqual(len(wells), 1)
                self.assertEqual(len(parents), 1)
                self.assertTrue(np.isfinite(wells.rmse).all())
