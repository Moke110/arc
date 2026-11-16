"""Real lifecycle entry points, observed before native construction; no GPU/model fit."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
import torch
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'output/api'))
from tenv5k import learner,model_transigen as native

EXPECTED={'policy':'author train_TranSiGen setup_seed CUDA',
          'cudnn_enabled':False,'cudnn_benchmark':False,'cudnn_deterministic':True}


def flags():
    return (torch.backends.cudnn.enabled,torch.backends.cudnn.benchmark,torch.backends.cudnn.deterministic)


class BackendContractTests(unittest.TestCase):
    def setUp(self):
        self.previous=flags();self.deterministic=torch.are_deterministic_algorithms_enabled()
        torch.backends.cudnn.enabled=True;torch.backends.cudnn.benchmark=True;torch.backends.cudnn.deterministic=False
    def tearDown(self):
        torch.backends.cudnn.enabled,torch.backends.cudnn.benchmark,torch.backends.cudnn.deterministic=self.previous
        torch.use_deterministic_algorithms(self.deterministic)

    def test_train_entry_applies_and_records_author_cuda_flags_before_native_build(self):
        data=SimpleNamespace(config={'early_stopping':{'validation_seed':22617}},rows=lambda *args:np.array([0,1]))
        observed=[]
        def observe(*args):observed.append(flags());raise RuntimeError('native construction boundary reached')
        with tempfile.TemporaryDirectory() as temporary,patch.object(learner,'Inputs',return_value=data),patch.object(native,'build_model',side_effect=observe):
            with self.assertRaisesRegex(RuntimeError,'construction boundary'):
                learner.train('unused',temporary,'unused',adapter='transigen',scenario='cold-drug',fold=0,seed=17,device='cuda',smoke_rows=2,max_epochs=1)
            self.assertEqual(observed,[(False,False,True)])
            self.assertEqual(json.loads((Path(temporary)/'config.json').read_text())['backend'],EXPECTED)

    def test_cpu_train_entry_preserves_backend_flags(self):
        data=SimpleNamespace(config={'early_stopping':{'validation_seed':22617}},rows=lambda *args:np.array([0,1]))
        before=flags();observed=[]
        def observe(*args):observed.append(flags());raise RuntimeError('native construction boundary reached')
        with tempfile.TemporaryDirectory() as temporary,patch.object(learner,'Inputs',return_value=data),patch.object(native,'build_model',side_effect=observe):
            with self.assertRaisesRegex(RuntimeError,'construction boundary'):
                learner.train('unused',temporary,'unused',adapter='transigen',scenario='cold-drug',fold=0,seed=17,device='cpu',smoke_rows=2,max_epochs=1)
            self.assertEqual(observed,[before]);self.assertEqual(flags(),before)

    def test_fresh_process_replay_entry_applies_flags_and_rejects_unrecorded_cuda_policy(self):
        script='''import json,sys
from unittest.mock import patch
from types import SimpleNamespace
sys.path.insert(0,sys.argv[1])
import torch
from tenv5k import learner,model_transigen as native
expected=json.loads(sys.argv[2])
torch.backends.cudnn.enabled=True
torch.backends.cudnn.benchmark=True
torch.backends.cudnn.deterministic=False
saved={'input_manifest_sha256':'frozen','config':{'adapter':'transigen','native_source_sha256':native.SOURCE_DIGEST,'device':'cuda','seed':17,'backend':expected}}
observed=[]
def build(*args):
    observed.append([torch.backends.cudnn.enabled,torch.backends.cudnn.benchmark,torch.backends.cudnn.deterministic])
    raise RuntimeError('native construction boundary reached')
with patch.object(learner,'Inputs',return_value=SimpleNamespace()),patch.object(learner,'sha256',return_value='frozen'),patch.object(torch,'load',return_value=saved),patch.object(native,'build_model',side_effect=build):
    try:learner.replay('unused','unused','unused',adapter='transigen')
    except RuntimeError as error:assert 'construction boundary' in str(error)
    assert observed==[[False,False,True]],observed
    saved['config'].pop('backend')
    try:learner.replay('unused','unused','unused',adapter='transigen')
    except ValueError as error:assert 'backend' in str(error)
    else:raise AssertionError('unrecorded old CUDA checkpoint was relabeled')
    saved['config']['device']='cpu'
    torch.backends.cudnn.enabled=True
    torch.backends.cudnn.benchmark=True
    torch.backends.cudnn.deterministic=False
    try:learner.replay('unused','unused','unused',adapter='transigen')
    except RuntimeError as error:assert 'construction boundary' in str(error)
    assert observed[-1]==[True,True,False],observed
print(json.dumps({'observed':observed,'legacy_cuda':'rejected'}))
'''
        result=subprocess.run([sys.executable,'-c',script,str(Path(__file__).resolve().parents[2]/'output/api'),json.dumps(EXPECTED)],capture_output=True,text=True,env={**os.environ,'OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'2'})
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(json.loads(result.stdout)['observed'],[[False,False,True],[True,True,False]])

if __name__=='__main__':unittest.main()
