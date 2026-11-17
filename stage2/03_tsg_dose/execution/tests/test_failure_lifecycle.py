"""Opt-in bounded native integration checks at evaluate/replay/resume seams."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
import numpy as np
import torch
EXECUTION=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(EXECUTION))
from run_matrix import train_or_resume,execution_registration,registration_path,run_path,write
from tsg_dose.learner import manifest,verify_candidate


@unittest.skipUnless(os.environ.get('TSG_DOSE_LIFECYCLE_INPUTS'),'requires bounded real input freeze')
class FailureLifecycle(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory();cls.root=Path(cls.temp.name)
        cls.inputs=Path(os.environ['TSG_DOSE_LIFECYCLE_INPUTS'])
        cls.source=Path(os.environ['TSG_DOSE_LIFECYCLE_SOURCE'])
        cls.base=cls.root/'base'
        cls.command=[sys.executable,str(EXECUTION/'transigen.py'),'--inputs',str(cls.inputs),'--source',str(cls.source)]
        subprocess.run(cls.command+['--candidate',str(cls.base),'--device','cuda','--smoke-rows','8','--max-epochs','2'],check=True,capture_output=True,text=True)

    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()

    def candidate(self,name):
        run=dict(scenario='cold-drug',fold=0,method='no-dose',seed=17)
        root=self.root/name;path=run_path(root,run);path.parent.mkdir(parents=True)
        shutil.copytree(self.base,path)
        frozen=dict(device='cuda',input_manifest_sha256='fixture')
        registration=execution_registration(run,frozen,'a'*64)
        write(registration_path(root,run),registration)
        return root,path,run,frozen

    def assert_failed_and_not_resumable(self,root,path,run,frozen,phase,error):
        terminal=json.loads((path/'terminal.json').read_text())
        self.assertEqual(terminal['status'],'failed')
        self.assertEqual(terminal['failure_phase'],phase)
        self.assertEqual(terminal['previous_terminal']['status'],'validation-complete')
        self.assertIn(error,terminal['error'])
        evidence=json.loads((path/'failure.json').read_text())
        self.assertEqual(evidence['previous_terminal'],terminal['previous_terminal'])
        self.assertIn('failure.json',verify_candidate(path)['outputs'])
        with self.assertRaisesRegex(ValueError,'failed scientific attempt'):
            train_or_resume(SimpleNamespace(candidate=root,resume=True),run,frozen)

    def test_nonfinite_evaluation_becomes_terminal_failure_and_blocks_resume(self):
        root,path,run,frozen=self.candidate('nonfinite')
        checkpoint=torch.load(path/'checkpoint.pt',map_location='cpu',weights_only=False)
        for value in checkpoint['state_dict'].values():
            if value.is_floating_point():value.fill_(float('nan'))
        torch.save(checkpoint,path/'checkpoint.pt');manifest(path,self.inputs)
        result=subprocess.run(self.command+['--candidate',str(path),'--evaluate-only'],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assert_failed_and_not_resumable(root,path,run,frozen,'evaluation','non-finite')

    def test_exact_replay_mismatch_becomes_terminal_failure_and_blocks_resume(self):
        root,path,run,frozen=self.candidate('replay-mismatch')
        prediction=np.load(path/'validation_prediction.npy');prediction[0,0]+=np.float32(1.)
        np.save(path/'validation_prediction.npy',prediction);manifest(path,self.inputs)
        result=subprocess.run(self.command+['--candidate',str(path),'--replay-only'],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assert_failed_and_not_resumable(root,path,run,frozen,'replay','replay differs')


if __name__=='__main__':unittest.main()
