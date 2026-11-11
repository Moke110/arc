import sys
from pathlib import Path
import tempfile
import unittest
import numpy as np
from PIL import Image
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from bases.metrics import retrieval, pooled_pr, METRICS
from bases.plots import render_retrieval, replay_retrieval

class PlotContractTests(unittest.TestCase):
    def test_tenv_panel_has_replayable_source_and_rejects_tahoe(self):
        result=retrieval(np.array([4.,-3.,-2.,1.,0.]),np.array([1,-1,1,-1,0,1]),np.arange(5))
        runs=[{'label':'Worked example','summary':{name:result[name].ap for name in METRICS},
               'curves':{name:pooled_pr([result[name]]) for name in METRICS}}]
        with tempfile.TemporaryDirectory() as tmp:
            target=Path(tmp)/'panel.png'
            report=render_retrieval(runs,target,source='tenv',diagnostic=True)
            self.assertTrue(target.exists())
            self.assertEqual(report['source'],'tenv')
            self.assertTrue(report['diagnostic'])
            self.assertIn('Diagnostic',report['title'])
            self.assertIn('not model performance',report['title'])
            self.assertTrue(target.with_suffix('.json').exists())
            self.assertTrue(target.with_suffix('.csv').exists())
            replay=Path(tmp)/'replay.png'
            replay_retrieval(target.with_suffix('.json'),replay)
            self.assertEqual(target.read_bytes(),replay.read_bytes())
            with Image.open(target) as image:
                self.assertGreater(image.width,image.height*2)
            with self.assertRaisesRegex(ValueError,'tenv'):
                render_retrieval(runs,Path(tmp)/'tahoe.png',source='tahoe')
