"""Fail-closed matrix planning/resume seams, without model imports."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
SCRIPT=Path(__file__).resolve().parents[1]/'run_matrix.py'
spec=importlib.util.spec_from_file_location('matrix231',SCRIPT)
matrix=importlib.util.module_from_spec(spec);spec.loader.exec_module(matrix)

class MatrixTests(unittest.TestCase):
    def test_plan_is_complete_and_unique(self):
        runs=matrix.plan()
        self.assertEqual(len(runs),54)
        self.assertEqual(len({tuple(r.values()) for r in runs}),54)
        self.assertEqual({r['seed'] for r in runs},{17,42,2026})
        self.assertEqual({r['fold'] for r in runs},{0,1,2})
    def test_import_is_unavailable_for_every_current_method(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as temporary:
            args=SimpleNamespace(candidate=Path(temporary))
            for run in matrix.plan():
                with self.assertRaisesRegex(ValueError, 'fresh training'):
                    matrix.import_run(args,run)
            target=Path(temporary)/'imported';target.mkdir()
            (target/'import-receipt.json').write_text('{}')
            with self.assertRaisesRegex(ValueError,'fresh training'):
                matrix.verify_execution(target,matrix.plan()[0],{'schema':'tenv5k-matrix-v3-threefold-f32-lr1e4'})

    def test_historical_acceptance_cannot_publish(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for acceptance in ({}, {'schema':'tenv5k-matrix-acceptance-v2','runs':90},
                               {'schema':matrix.ACCEPTANCE_SCHEMA,'runs':53}):
                (root/'acceptance.json').write_text(json.dumps(acceptance))
                with self.assertRaisesRegex(ValueError,'current three-fold acceptance'):
                    matrix.publish(SimpleNamespace(candidate=root))

    def test_failed_smoke_missing_digest_resume_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            for status,smoke in [('failed',False),('complete',True),('complete',False)]:
                (root/'terminal.json').write_text(json.dumps(dict(status=status,smoke=smoke)))
                with self.assertRaises(ValueError): matrix.verify_run(root,dict(adapter='mlp',scenario='cold-drug',fold=0,seed=17),'input')
    def test_matrix_requires_every_registered_identity(self):
        with self.assertRaises(ValueError): matrix.require_complete(plan=matrix.plan(),found=matrix.plan()[:-1])
        with self.assertRaises(ValueError): matrix.require_complete(plan=matrix.plan(),found=matrix.plan()[:-1]+matrix.plan()[:1])

    def test_resume_preserves_interrupted_evidence_and_replays_trained_identity(self):
        from types import SimpleNamespace
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);run=matrix.plan()[0];target=matrix.run_path(root,run);target.mkdir(parents=True)
            (target/'partial.bin').write_bytes(b'precious interrupted evidence')
            args=SimpleNamespace(candidate=root,device='cpu')
            action=matrix.recover_existing(args,run,'input',lambda:None)
            self.assertEqual(action,'train')
            self.assertFalse(target.exists())
            self.assertEqual(next((root/'attempts').rglob('partial.bin')).read_bytes(),b'precious interrupted evidence')
            target.mkdir(parents=True)
            config={**run,'device':'cpu','smoke_rows':None,'max_epochs':1000,'patience':50,'strict_lower':True}
            (target/'config.json').write_text(json.dumps(config))
            (target/'terminal.json').write_text(json.dumps(dict(status='trained-awaiting-fresh-process-replay',smoke=False)))
            calls=[]
            with patch.object(matrix,'verify_run',return_value=config):
                self.assertEqual(matrix.recover_existing(args,run,'input',lambda:calls.append('replay-only')),'complete')
            self.assertEqual(calls,['replay-only'])
            self.assertTrue(list((root/'attempts').rglob('terminal.before-replay.json')))

if __name__=='__main__':unittest.main()

class ScientificAcceptanceTests(unittest.TestCase):
    def test_each_method_requires_new_identity_float32_prediction_and_checkpoint(self):
        import numpy as np
        import torch
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for adapter in ('prnet','transigen','mlp'):
                config=dict(method_identity=adapter+'-f32-lr1e4',compute_dtype='float32',
                            initialization_dtype='float32',prediction_dtype='float32',learning_rate=1e-4)
                method=dict(scientific_identity=config['method_identity'],scientific_configuration=config.copy())
                def save(dtype=torch.float32,prediction_dtype=np.float32):
                    torch.save(dict(config=config.copy(),state_dict={'weight':torch.ones(1,dtype=dtype),'counter':torch.ones(1,dtype=torch.int64)}),root/'checkpoint.pt')
                    np.save(root/'prediction.npy',np.ones((2,3),dtype=prediction_dtype))
                save()
                matrix.verify_scientific_run(root,config,method)
                save(prediction_dtype=np.float64)
                with self.assertRaisesRegex(ValueError,'prediction.*float32'):
                    matrix.verify_scientific_run(root,config,method)
                save(dtype=torch.float64)
                with self.assertRaisesRegex(ValueError,'state.*float32'):
                    matrix.verify_scientific_run(root,config,method)
                save();config['learning_rate']=1e-3
                with self.assertRaisesRegex(ValueError,'scientific configuration'):
                    matrix.verify_scientific_run(root,config,method)
                config['learning_rate']=1e-4;config['method_identity']='historical'
                with self.assertRaisesRegex(ValueError,'method identity'):
                    matrix.verify_scientific_run(root,config,method)
