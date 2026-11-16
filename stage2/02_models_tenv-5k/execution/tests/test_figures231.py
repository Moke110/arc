"""CSV-only rendering survives absence of every model/input payload."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import pandas as pd
RENDER=Path(__file__).resolve().parents[2]/'figure/script/render_figures.py'
class FigureTests(unittest.TestCase):
    def test_csv_only_exact_redraw(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);data=root/'csv';data.mkdir()
            pd.DataFrame([dict(scenario=s,method=m,metric=k,value=.2,n_parents=3,n_undefined=0,zero=0) for s in ('cold-drug','warm-start') for m in ('PRnet','TranSiGen','Mean') for k in ('pcc','spearman','rmse','r2')]).to_csv(data/'metrics_values.csv',index=False)
            pd.DataFrame([dict(scenario='cold-drug',candidate='PRnet',reference='Mean',metric='pcc',estimate=.1,low=-.1,high=.3,zero=0)]).to_csv(data/'paired_differences_values.csv',index=False)
            pd.DataFrame([dict(scenario='cold-drug',method='PRnet',fold=0,seed=17,epoch=e,train_loss=v,validation_mse=v,best_epoch=2) for e,v in [(1,2.),(2,1.)]]).to_csv(data/'learning_curve_values.csv',index=False)
            for folder in ('a','b'):
                subprocess.run([sys.executable,str(RENDER),'--figure-data',str(data),'--out-dir',str(root/folder)],check=True)
            self.assertEqual({p.name:p.read_bytes() for p in (root/'a').iterdir()},{p.name:p.read_bytes() for p in (root/'b').iterdir()})
            self.assertEqual(len(list((root/'a').glob('*.png'))),3)
if __name__=='__main__':unittest.main()
