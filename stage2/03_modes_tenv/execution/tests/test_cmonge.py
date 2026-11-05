"""CMonge adapter checked at the shared subprocess/artifact seam."""
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'script'))
sys.path.insert(0,str(Path(__file__).resolve().parent))
from test_pipeline import fixture
from artifacts import read_json, write_json, write_csv, rows
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
            # More than one held-out condition makes stochastic selection visible.
            validation=rows(root/'input/fold_0/validation.csv')
            validation.append(dict(validation[0],case_id='validation2',group_id='validation2'))
            write_csv(root/'input/fold_0/validation.csv',validation)
            np.save(root/'input/fold_0/validation_expression.npy',np.array([[4,2,-2,-4],[2,1,-1,-2]],dtype='float32'))
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
            self.assertEqual(checkpoint['validation_batch_policy'], 'one seeded random held-out condition per validation event')
            self.assertEqual([event['validation_conditions_evaluated'] for event in checkpoint['history']], [1, 1, 1])
            valid_groups={row['group_id'] for row in rows(root/'input/fold_0/validation.csv')}
            self.assertTrue(all(event['validation_group_id'] in valid_groups for event in checkpoint['history']))
            run_one(root/'input',root/'replay','cmonge',0,17,command)
            replay=read_json(root/'replay/runs/cmonge/fold_0/seed_17/checkpoint.json')
            self.assertEqual(checkpoint['history'],replay['history'])
            np.testing.assert_array_equal(np.load(root/'output/runs/cmonge/fold_0/seed_17/prediction.npy'),
                                          np.load(root/'replay/runs/cmonge/fold_0/seed_17/prediction.npy'))
            self.assertEqual(validate_run(root/'input',root/'output','cmonge',0,17,command)['status'],'completed')

if __name__=='__main__':unittest.main()
