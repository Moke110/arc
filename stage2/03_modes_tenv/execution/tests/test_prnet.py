"""PRnet at the agreed isolated fit/save/reload/predict artifact seam."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np

STEP = Path(__file__).resolve().parents[1]
ADAPTER = STEP / 'script/prnet_adapter.py'
PYTHON = '/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python'


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class PRnetSubprocess(unittest.TestCase):
    def test_full_axis_training_reload_and_allowed_inputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            n_genes = 18164
            rng = np.random.default_rng(44)
            controls = rng.uniform(.1, 1., (2, n_genes)).astype('float32')
            np.save(root/'controls.npy', controls)
            np.save(root/'morgan.npy', rng.integers(0, 2, (3, 2048), dtype='uint8'))
            write_csv(root/'drugs.csv', [{'drug_key':f'd{i}', 'morgan_index':i} for i in range(3)])
            files = {k:str(root/f'{k}.npy') for k in ('controls','morgan')}
            files['drugs'] = str(root/'drugs.csv')
            for role, drug, count in [('train',0,6), ('validation',1,2), ('test',2,2)]:
                cases = [{'case_id':f'{role}{i}', 'drug_key':f'd{drug}', 'dose_um':1+i,
                          'control_index':i%2, 'group_id':f'{role}{i}'} for i in range(count)]
                path = root/f'{role}.csv'
                write_csv(path, cases)
                files[f'{role}_cases'] = str(path)
                if role != 'test':
                    delta = rng.uniform(-.05,.05,(count,n_genes)).astype('float32')
                    np.save(root/f'{role}_expression.npy', delta)
                    files[f'{role}_expression'] = str(root/f'{role}_expression.npy')
            output = root/'output'
            output.mkdir()
            request = {'identity':{'model':'prnet','fold':0,'seed':17,'input_sha256':'fixture'},
                       'n_genes':n_genes,'target_kind':'delta','predictor_files':files,
                       'output_directory':str(output)}
            (root/'request.json').write_text(json.dumps(request))
            config = json.loads((STEP/'models/prnet/config.json').read_text())
            config.update(max_epochs=2,batch_size=2,device='cpu')
            (root/'pilot_config.json').write_text(json.dumps(config))
            process = subprocess.run([PYTHON,str(ADAPTER),'--request',str(root/'request.json'),
                                      '--config',str(root/'pilot_config.json')],capture_output=True,text=True)
            self.assertEqual(process.returncode,0,process.stdout+process.stderr)
            prediction = np.load(output/'prediction.npy',mmap_mode='r')
            self.assertEqual(prediction.shape,(2,n_genes))
            self.assertTrue(np.isfinite(prediction).all())
            self.assertEqual(json.loads((output/'case_ids.json').read_text()),['test0','test1'])
            checkpoint = json.loads((output/'checkpoint.json').read_text())
            self.assertEqual(checkpoint['identity'],request['identity'])
            self.assertEqual(checkpoint['author_objective'],'GaussianNLLLoss')
            self.assertEqual(checkpoint['adapter_parameters'],2048*1024+1024)
            for key in ['adapter_parameter_change_l2','predictor_parameter_change_l2','baseline_sensitivity_l2']:
                self.assertGreater(checkpoint['behavior'][key],0,key)
            self.assertEqual(checkpoint['behavior']['reload_max_abs_error'],0.)
            self.assertTrue(np.isfinite(checkpoint['best_validation_loss']))
            self.assertEqual(len(checkpoint['history']),2)
            self.assertEqual(checkpoint['stop_reason'],'max_epochs')
            for name,digest in checkpoint['blobs'].items():
                self.assertEqual(hashlib.sha256((output/name).read_bytes()).hexdigest(),digest)
            progress = json.loads((output/'progress.json').read_text())
            self.assertEqual(progress['phase'],'complete')
            self.assertGreater(progress['training_seconds'],0.)


if __name__ == '__main__':
    unittest.main()
