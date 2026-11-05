"""Native PerturbDiff through the approved artifact subprocess boundary."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np
from test_prnet import write_csv

STEP = Path(__file__).resolve().parents[1]
PYTHON = str(STEP.parents[1]/'_envs/perturbdiff/bin/python')

class PerturbDiffArtifact(unittest.TestCase):
    def test_joint_training_reload_and_prediction(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rng = np.random.default_rng(91)
            n_genes = 64
            np.save(root/'controls.npy',rng.uniform(.2,1,(2,n_genes)).astype('float32'))
            np.save(root/'morgan.npy',rng.integers(0,2,(3,2048),dtype='uint8'))
            write_csv(root/'drugs.csv',[{'drug_key':f'd{i}','morgan_index':i} for i in range(3)])
            files={k:str(root/f'{k}.npy') for k in ('controls','morgan')}
            files['drugs']=str(root/'drugs.csv')
            for j,role in enumerate(('train','validation','test')):
                cases=[{'case_id':f'{role}{i}','drug_key':f'd{j}','dose_um':1+i//2,
                        'control_index':i%2,'group_id':f'{role}_{i//2}'} for i in range(4)]
                write_csv(root/f'{role}.csv',cases)
                files[f'{role}_cases']=str(root/f'{role}.csv')
                if role!='test':
                    np.save(root/f'{role}.npy',rng.uniform(-.1,.1,(4,n_genes)).astype('float32'))
                    files[f'{role}_expression']=str(root/f'{role}.npy')
            request={'identity':{'model':'perturbdiff','fold':0,'seed':17},'n_genes':n_genes,
                     'target_kind':'delta','predictor_files':files,'output_directory':str(root/'output')}
            (root/'request.json').write_text(json.dumps(request))
            config=json.loads((STEP/'models/perturbdiff/config.json').read_text())
            config.update(max_steps=12,validation_interval=12,batch_groups=2,population_size=2,device='cpu',sampling_start_time=2)
            (root/'config.json').write_text(json.dumps(config))
            result=subprocess.run([PYTHON,str(STEP/'script/perturbdiff_adapter.py'),'--request',str(root/'request.json'),
                                   '--config',str(root/'config.json')],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            output=root/'output'
            pred=np.load(output/'prediction.npy')
            self.assertEqual(pred.shape,(4,n_genes))
            self.assertTrue(np.isfinite(pred).all())
            checkpoint=json.loads((output/'checkpoint.json').read_text())
            self.assertEqual(checkpoint['identity'],request['identity'])
            self.assertEqual(checkpoint['author_objective'],'native GaussianDiffusion MSE + energy MMD')
            self.assertEqual(checkpoint['architecture'],{'depth':8,'hidden_size':512,'diffusion_steps':1000})
            for name in ('adapter_parameter_change_l2','baseline_sensitivity_l2','dose_sensitivity_l2'):
                self.assertGreater(checkpoint['behavior'][name],0,name)
            self.assertTrue(checkpoint['behavior']['predictor_parameters_changed'])
            self.assertEqual(checkpoint['behavior']['reload_max_abs_error'],0)
            self.assertEqual(json.loads((output/'case_ids.json').read_text()),[f'test{i}' for i in range(4)])

if __name__=='__main__':
    unittest.main()
