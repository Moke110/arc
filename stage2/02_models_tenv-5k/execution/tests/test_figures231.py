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
            pd.DataFrame([dict(scenario=s,method=m,metric=k,value=.2,n_parents=3,n_undefined=0,zero=0,applicable=True,status='complete',n_constant_truth=1,n_constant_prediction=2,constant_truth_ids='drug-x',reason='truth constant excludes correlation, errors retained') for s in ('cold-drug','warm-start') for m in ('PRnet','TranSiGen','Mean') for k in ('pcc','spearman','rmse','r2')]).to_csv(data/'metrics_values.csv',index=False)
            metric_values=pd.read_csv(data/'metrics_values.csv')
            metric_values=pd.concat([metric_values,pd.DataFrame([dict(scenario='cold-drug',method='Mean_drug',metric=k,value=float('nan'),n_parents=0,n_undefined=3,zero=0,applicable=False,status='not-applicable',n_constant_truth=1,n_constant_prediction=0,constant_truth_ids='drug-x',reason='unseen test drug') for k in ('pcc','spearman','rmse','r2')])],ignore_index=True)
            metric_values.to_csv(data/'metrics_values.csv',index=False)
            pd.DataFrame([dict(scenario='cold-drug',candidate='PRnet',reference='Mean',metric='pcc',estimate=.1,low=-.1,high=.3,zero=0,n_excluded_constant=1,candidate_constant_prediction=2,reference_constant_prediction=0,unavailable_reference='Mean_drug: unseen test drug')]).to_csv(data/'paired_differences_values.csv',index=False)
            pd.DataFrame([dict(scenario='cold-drug',method='PRnet',fold=0,seed=17,epoch=e,train_loss=v,validation_mse=v,best_epoch=2) for e,v in [(1,2.),(2,1.)]]).to_csv(data/'learning_curve_values.csv',index=False)
            for folder in ('a','b'):
                subprocess.run([sys.executable,str(RENDER),'--figure-data',str(data),'--out-dir',str(root/folder)],check=True)
            self.assertEqual({p.name:p.read_bytes() for p in (root/'a').iterdir()},{p.name:p.read_bytes() for p in (root/'b').iterdir()})
            self.assertEqual(len(list((root/'a').glob('*.png'))),3)
            import importlib.util
            spec=importlib.util.spec_from_file_location('renderer231',RENDER);renderer=importlib.util.module_from_spec(spec);spec.loader.exec_module(renderer)
            text=renderer.metric_annotation(metric_values[metric_values.method=='Mean_drug'].iloc[0])
            self.assertIn('N/A',text);self.assertIn('unseen test drug',text)
            text=renderer.metric_annotation(metric_values[metric_values.method=='PRnet'].iloc[0])
            self.assertIn('truth constant=1',text);self.assertIn('pred constant=2',text)
if __name__=='__main__':unittest.main()
