"""Native chemCPA qualification through the artifact subprocess seam."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np

STEP = Path(__file__).resolve().parents[1]
PYTHON = '/home/chang/miniconda3/envs/arc-stage2-bench-chemcpa/bin/python'
ADAPTER = STEP/'script/chemcpa_adapter.py'


class ChemCPASubprocess(unittest.TestCase):
    def test_full_axis_preserves_adversarial_training_and_cold_prediction(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rng = np.random.default_rng(13)
            genes = 18164
            np.save(root/'controls.npy',rng.uniform(.1,1.,(2,genes)).astype('float32'))
            np.save(root/'morgan.npy',rng.integers(0,2,(4,2048),dtype='uint8'))
            def csv_file(name,records):
                path = root/name
                with path.open('w',newline='') as stream:
                    writer = csv.DictWriter(stream,fieldnames=list(records[0]))
                    writer.writeheader()
                    writer.writerows(records)
                return str(path)
            files = {'controls':str(root/'controls.npy'),'morgan':str(root/'morgan.npy'),
                     'drugs':csv_file('drugs.csv',[{'drug_key':f'd{i}','morgan_index':i} for i in range(4)])}
            for role,count in [('train',8),('validation',2),('test',2)]:
                records = [{'case_id':f'{role}{i}','drug_key':f'd{i%2 if role=="train" else (2 if role=="validation" else 3)}',
                            'dose_um':.1+i,'control_index':i%2,'group_id':f'{role}{i}'} for i in range(count)]
                files[role+'_cases'] = csv_file(role+'.csv',records)
                if role != 'test':
                    path = root/(role+'_expression.npy')
                    np.save(path,rng.normal(.02,.01,(count,genes)).astype('float32'))
                    files[role+'_expression'] = str(path)
            output = root/'output'
            output.mkdir()
            request = {'identity':{'model':'chemcpa','fold':0,'seed':42,'input_sha256':'synthetic'},
                       'n_genes':genes,'target_kind':'delta','predictor_files':files,'output_directory':str(output)}
            (root/'request.json').write_text(json.dumps(request))
            config = json.loads((STEP/'models/chemcpa/config.json').read_text())
            config.update(max_epochs=2,batch_size=2,device='cpu')
            (root/'config.json').write_text(json.dumps(config))
            result = subprocess.run([PYTHON,str(ADAPTER),'--request',str(root/'request.json'),
                                     '--config',str(root/'config.json')],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            prediction = np.load(output/'prediction.npy',mmap_mode='r')
            self.assertEqual(prediction.shape,(2,genes))
            self.assertTrue(np.isfinite(prediction).all())
            self.assertEqual(json.loads((output/'case_ids.json').read_text()),['test0','test1'])
            checkpoint = json.loads((output/'checkpoint.json').read_text())
            self.assertEqual(checkpoint['identity'],request['identity'])
            self.assertEqual(checkpoint['adapter_parameters'],2048*194+194)
            self.assertEqual(checkpoint['author_objective'],'GaussianNLL - reg_adversary * drug_BCE')
            self.assertGreater(checkpoint['adversary_steps'],0)
            self.assertGreater(checkpoint['autoencoder_steps'],0)
            self.assertGreater(checkpoint['history'][0]['adversary_gradient_penalty'],0.)
            for name in ['adapter_parameter_change_l2','predictor_parameter_change_l2','adversary_parameter_change_l2',
                         'doser_parameter_change_l2','baseline_sensitivity_l2','dose_sensitivity_l2']:
                self.assertGreater(checkpoint['behavior'][name],0.,name)
            self.assertEqual(checkpoint['behavior']['reload_max_abs_error'],0.)
            self.assertTrue(np.isfinite(checkpoint['best_validation_loss']))
            for blob,digest in checkpoint['blobs'].items():
                self.assertEqual(hashlib.sha256((output/blob).read_bytes()).hexdigest(),digest)


if __name__ == '__main__':
    unittest.main()
