"""CMonge adapter checked at the shared subprocess/artifact seam."""
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'script'))
sys.path.insert(0,str(Path(__file__).resolve().parent))
from test_pipeline import fixture
from artifacts import read_json, write_json, rows
from benchmark import run_one, validate_run
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'.git').exists())
ENV=ROOT/'stage2/_envs/cmonge/bin/python'
ADAPTER=ROOT/'stage2/03_modes_tenv/execution/script/cmonge_adapter.py'

class CMongeRuntime(unittest.TestCase):
    @unittest.skipUnless(ENV.exists(),'isolated CMonge environment required')
    def test_joint_adapter_author_loss_reload_and_input_participation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); fixture(root/'input')
            fingerprints=np.zeros((6,2048),dtype='uint8')
            fingerprints[:,0]=1
            for i in range(6):fingerprints[i,i+1]=1
            np.save(root/'input/morgan.npy',fingerprints)
            config=root/'config.json'
            write_json(config,{'iterations':3,'batch_size':4,'eval_every':1,'max_groups':2,'qualification':True})
            command=[str(ENV),str(ADAPTER),'--config',str(config)]
            try:
                run=run_one(root/'input',root/'output','cmonge',0,17,command)
            except Exception:
                log=root/'output/runs/cmonge/fold_0/seed_17/adapter.log'
                if log.exists():print(log.read_text(),file=sys.stderr)
                raise
            checkpoint=read_json(root/'output/runs/cmonge/fold_0/seed_17/checkpoint.json')
            self.assertGreater(checkpoint['adapter_parameter_change'],0)
            self.assertGreater(checkpoint['baseline_sensitivity'],0)
            self.assertGreater(checkpoint['dose_sensitivity'],0)
            self.assertTrue(np.isfinite(checkpoint['best_validation_loss']))
            self.assertEqual(checkpoint['reload_max_absolute_error'],0)
            self.assertEqual(checkpoint['author_objective'],'Sinkhorn divergence + Monge gap (Euclidean)')
            self.assertEqual(validate_run(root/'input',root/'output','cmonge',0,17,command)['status'],'completed')

if __name__=='__main__':unittest.main()
