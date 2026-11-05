"""CPA through the approved isolated artifact seam."""
import csv
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np

EXECUTION=Path(__file__).resolve().parent.parent
PYTHON='/home/chang/miniconda3/envs/arc-stage2-bench-cpa/bin/python'


class CPAArtifacts(unittest.TestCase):
    def test_native_training_changes_adapter_and_baseline_controls_prediction(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            rng=np.random.default_rng(9)
            n_genes=64
            files={}
            for name,array in [('controls',rng.uniform(.5,1.5,(2,n_genes)).astype('float32')),
                               ('morgan',rng.integers(0,2,(3,2048),dtype='uint8'))]:
                np.save(root/f'{name}.npy',array);files[name]=str(root/f'{name}.npy')
            def table(name,records):
                path=root/f'{name}.csv'
                with path.open('w',newline='') as stream:
                    writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
                return str(path)
            files['drugs']=table('drugs',[{'drug_key':f'd{i}','morgan_index':i} for i in range(3)])
            for role,drug,n in [('train',0,8),('validation',1,4),('test',2,4)]:
                files[role+'_cases']=table(role,[{'case_id':f'{role}{i}','drug_key':f'd{drug}',
                    'dose_um':float(i+1),'control_index':i%2,'group_id':f'g{role}{i}'} for i in range(n)])
                if role!='test':
                    np.save(root/f'{role}_expression.npy',rng.uniform(-.1,.1,(n,n_genes)).astype('float32'))
                    files[role+'_expression']=str(root/f'{role}_expression.npy')
            output=root/'run'
            request={'identity':{'model':'cpa','fold':0,'seed':17},'predictor_files':files,
                     'target_kind':'delta','n_genes':n_genes,'output_directory':str(output)}
            (root/'request.json').write_text(json.dumps(request))
            config=json.loads((EXECUTION/'models/cpa/config.json').read_text())
            config.update(max_epochs=2,batch_size=4,device='cpu',n_hidden_encoder=32,n_hidden_decoder=32,n_latent=8,n_hidden_adv=16)
            (root/'config.json').write_text(json.dumps(config))
            result=subprocess.run([PYTHON,str(EXECUTION/'script/cpa_adapter.py'),'--request',str(root/'request.json'),'--config',str(root/'config.json')],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            checkpoint=json.loads((output/'checkpoint.json').read_text())
            self.assertEqual(checkpoint['author_objective'],'CPA native Gaussian reconstruction NLL + native adversarial disentanglement')
            self.assertEqual(checkpoint['identity'],request['identity'])
            self.assertEqual(checkpoint['stop_reason'],'max_epochs')
            self.assertTrue(np.isfinite(checkpoint['best_validation_loss']))
            for name in ['adapter_parameter_change_l2','predictor_parameter_change_l2','baseline_sensitivity_l2','adversary_parameter_change_l2']:
                self.assertGreater(checkpoint['behavior'][name],0,name)
            self.assertEqual(checkpoint['behavior']['reload_max_abs_error'],0.)
            self.assertEqual(np.load(output/'prediction.npy').shape,(4,n_genes))
            self.assertTrue(np.isfinite(np.load(output/'prediction.npy')).all())
            self.assertEqual(json.loads((output/'case_ids.json').read_text()),[f'test{i}' for i in range(4)])

if __name__=='__main__': unittest.main()
