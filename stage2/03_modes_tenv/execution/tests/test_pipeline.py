"""Artifact CLI seam with small frozen inputs and real adapter subprocesses."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1] / 'script'))
from artifacts import write_csv, write_json, sha256, read_json
from benchmark import run_one, freeze_predictions, verify_barrier, validate_run, evaluate
SCRIPT = Path(__file__).resolve().parents[1] / 'script' / 'benchmark.py'


def fixture(root):
    root.mkdir(parents=True)
    write_csv(root / 'genes.csv',[{'gene':'a'},{'gene':'b'},{'gene':'c'},{'gene':'d'}])
    write_csv(root / 'folds.csv',[{'fold':i,'drug_key':f'd{i}','role':'test'} for i in range(3)])
    write_csv(root / 'drugs.csv',[{'drug_key':d,'morgan_index':i,'parent_connectivity':d} for i,d in enumerate(['d0','d1','d2','train','train2','validation'])])
    np.save(root / 'controls.npy',np.zeros((1,4),dtype='float32'))
    np.save(root / 'morgan.npy',np.zeros((6,2048),dtype='uint8'))
    for fold in range(3):
        directory=root / f'fold_{fold}'
        train=[{'case_id':'train1','drug_key':'train','dose_um':1,'control_index':0,'group_id':'g1'},
               {'case_id':'train2','drug_key':'train','dose_um':1,'control_index':0,'group_id':'g1'},
               {'case_id':'train3','drug_key':'train2','dose_um':1,'control_index':0,'group_id':'g2'}]
        test=[{'case_id':'test1','drug_key':f'd{fold}','dose_um':1,'control_index':0,'group_id':'test'},
              {'case_id':'test2','drug_key':f'd{fold}','dose_um':1,'control_index':0,'group_id':'test'}]
        write_csv(directory / 'train.csv',train)
        write_csv(directory / 'validation.csv',[{'case_id':'validation','drug_key':'validation','dose_um':1,'control_index':0,'group_id':'validation'}])
        write_csv(directory / 'test_cases.csv',test)
        np.save(directory / 'train_expression.npy',np.array([[4,2,-2,-4],[4,2,-2,-4],[0,0,0,0]],dtype='float32'))
        np.save(directory / 'validation_expression.npy',np.array([[4,2,-2,-4]],dtype='float32'))
        write_json(directory / 'predictor.json',{'target_kind':'delta'})
        write_csv(directory / 'evaluator' / 'groups.csv',[{'group_id':'test','drug_key':f'd{fold}','dose_um':1,'n_wells':2}])
        np.save(directory / 'evaluator' / 'truth.npy',np.array([[2,1,-1,-2]],dtype='float32'))
        np.save(directory / 'evaluator' / 'positive.npy',np.array([[1,1,0,0]],bool))
        np.save(directory / 'evaluator' / 'negative.npy',np.array([[0,0,1,1]],bool))
    write_json(root / 'manifest.json',{'training_seeds':[17,42,2026],'synthetic':True})

class Pipeline(unittest.TestCase):
    def test_two_baselines_persist_and_replay_three_figures(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); fixture(root / 'input')
            result=subprocess.run([sys.executable,str(SCRIPT),'synthetic','--input',str(root / 'input'),'--output',str(root / 'output')],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            manifest=json.loads((root / 'output' / 'predictions.json').read_text())
            self.assertEqual(len(manifest['runs']),18)
            self.assertTrue(all(r['status']=='completed' for r in manifest['runs']))
            expected=np.array([2,1,-1,-2])
            np.testing.assert_array_equal(np.load(root / 'output/runs/delta_mean/fold_0/seed_17/prediction.npy')[0],expected)
            from artifacts import rows
            no_change=[r for r in rows(root/'output/per_run.csv') if r['model']=='no_change']
            self.assertEqual(float(no_change[0]['AUPRC_positive@100']),.5)
            self.assertEqual(int(no_change[0]['PCC_valid_groups']),0)
            for figure in ('metrics.png','pr_panels.png','best_run.png'):
                self.assertGreater((root / 'output' / figure).stat().st_size,1000)
            hashes={p.name:sha256(p) for p in (root / 'output').glob('*.csv')}
            result=subprocess.run([sys.executable,str(SCRIPT),'replay','--output',str(root / 'output')],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual(hashes,{p.name:sha256(p) for p in (root / 'output').glob('*.csv')})
            # Completion and metric definition are different: one constant-truth
            # fold leaves six defined correlation runs out of nine completions.
            np.save(root/'input/fold_0/evaluator/truth.npy',np.zeros((1,4),dtype='float32'))
            summary=evaluate(root/'input',root/'output',workers=1)
            self.assertEqual(summary[1]['PCC_valid_runs'],6)
            import xml.etree.ElementTree as ET
            svg=ET.parse(root/'output/metrics.svg')
            visible=' '.join(''.join(node.itertext()) for node in svg.iter() if node.tag.endswith('}text'))
            self.assertIn('delta_mean (completed 9/9; defined 6/9)',visible)
            self.assertIn('no_change (completed 9/9; defined 0/9)',visible)
            failed=root/'output/runs/delta_mean/fold_2/seed_2026/run.json'
            run=read_json(failed);run['status']='timeout';write_json(failed,run)
            freeze_predictions(root/'input',root/'output',['no_change','delta_mean'])
            summary=evaluate(root/'input',root/'output',workers=1)
            self.assertFalse(summary[1]['complete_ranking'])
            self.assertEqual(summary[1]['completed_runs'],8)

    def test_barrier_and_stale_predictions_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); fixture(root / 'input')
            run_one(root/'input',root/'output','no_change',0,17)
            with self.assertRaisesRegex(ValueError,'barrier'):
                freeze_predictions(root/'input',root/'output',['no_change'])
            path=root/'output/runs/no_change/fold_0/seed_17/case_ids.json'
            write_json(path,['test2','test1'])
            with self.assertRaisesRegex(ValueError,'checksum'):
                validate_run(root/'input',root/'output','no_change',0,17)

    def test_test_truth_and_global_indices_cannot_enter_predictor_cases(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); fixture(root/'input')
            write_csv(root/'input/fold_0/test_cases.csv',[{'case_id':'leak','drug_key':'d0','dose_um':1,'control_index':0,'group_id':'test','global_row_index':12}])
            with self.assertRaisesRegex(ValueError,'leakage'):
                run_one(root/'input',root/'output','no_change',0,17)

    def test_adapter_wrong_shape_nonfinite_and_case_order_fail(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); fixture(root/'input')
            adapter=root/'adapter.py'
            write_json(root/'adapter.identity.json',{'source_files':['adapter.py'],'environment_files':[sys.executable],'config_files':['adapter.py']})
            adapter.write_text("import argparse,json,numpy as np\nfrom pathlib import Path\np=argparse.ArgumentParser();p.add_argument('--request');r=json.loads(Path(p.parse_args().request).read_text());o=Path(r['output_directory']);np.save(o/'prediction.npy',np.zeros((1,4)));(o/'case_ids.json').write_text(json.dumps(['test1','test2']));(o/'checkpoint.json').write_text(json.dumps({'identity':r['identity']}))\n")
            with self.assertRaisesRegex(ValueError,'shape'):
                run_one(root/'input',root/'output','fixture',0,17,[sys.executable,str(adapter)])
            adapter.write_text(adapter.read_text().replace('np.zeros((1,4))','np.full((2,4),np.nan)'))
            with self.assertRaisesRegex(ValueError,'nonfinite'):
                run_one(root/'input',root/'output','fixture',0,17,[sys.executable,str(adapter)])
            adapter.write_text(adapter.read_text().replace('np.full((2,4),np.nan)','np.zeros((2,4))').replace("['test1','test2']","['test2','test1']"))
            with self.assertRaisesRegex(ValueError,'misordered'):
                run_one(root/'input',root/'output','fixture',0,17,[sys.executable,str(adapter)])

if __name__=='__main__': unittest.main()
