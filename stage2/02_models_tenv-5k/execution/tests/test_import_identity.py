"""Public matrix identity/import acceptance boundaries; no scientific training."""
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'run_matrix.py'
spec = importlib.util.spec_from_file_location('matrix_import', SCRIPT)
matrix = importlib.util.module_from_spec(spec)
spec.loader.exec_module(matrix)


class ImportIdentityTests(unittest.TestCase):
    def test_plan_declares_shared_and_three_method_identities(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ('inputs', 'means'):
                (root / name).mkdir()
                (root / name / 'manifest.json').write_text('{}')
            args = SimpleNamespace(inputs=root/'inputs', means=root/'means', device='cpu',
                                   prnet_source=root/'prnet', transigen_source=root/'transigen')
            frozen = matrix.frozen_plan(args)
            self.assertEqual(frozen['schema'], 'tenv5k-matrix-v3-threefold-f32-lr1e4')
            self.assertEqual(set(frozen['methods']), {'prnet', 'transigen', 'mlp'})
            self.assertEqual(frozen['methods']['prnet']['scientific_identity'], 'prnet-native-f32-lr1e4-v1')
            self.assertEqual(frozen['methods']['transigen']['scientific_identity'], 'transigen-native-f32-lr1e4-v1')
            self.assertIn('output/api/tenv5k/protocol.py', frozen['shared']['evaluation_code_sha256'])
            self.assertNotIn('output/api/tenv5k/model_transigen.py', frozen['methods']['prnet']['code_sha256'])

    def test_historical_candidate_cannot_resume_as_current_matrix(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for name in ('inputs','means'):
                (root/name).mkdir();(root/name/'manifest.json').write_text('{}')
            candidate=root/'matrix';candidate.mkdir()
            args=SimpleNamespace(inputs=root/'inputs',means=root/'means',candidate=candidate,
                                 device='cpu',prnet_source=root/'prnet',transigen_source=root/'transigen',resume=True)
            current=matrix.frozen_plan(args)
            for schema in ('tenv5k-matrix-v1','tenv5k-matrix-v2'):
                (candidate/'plan.json').write_text(json.dumps({**current,'schema':schema}))
                with self.assertRaisesRegex(ValueError,'resume plan/config/code identity changed'):
                    matrix.execute(args)
                self.assertFalse((candidate/'registrations').exists())
                self.assertFalse((candidate/'runs').exists())

    def test_import_cannot_relabel_transigen_or_change_shared_scoring(self):
        with tempfile.TemporaryDirectory() as directory:
            args = SimpleNamespace(candidate=Path(directory), import_original=Path(directory)/'old',
                                   import_runtime=Path(directory)/'runtime', compatibility_audit=None)
            with self.assertRaisesRegex(ValueError, 'fresh training'):
                matrix.import_run(args, dict(adapter='transigen', scenario='cold-drug', fold=0, seed=17))

class ImportedEvidenceTests(unittest.TestCase):
    """Verifier fixtures represent replay receipts, not substitute model training."""
    def setUp(self):
        import copy
        import sys
        sys.path.insert(0,str(SCRIPT.parent))
        import import_matrix
        self.imports=import_matrix
        self.temporary=tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root=Path(self.temporary.name)
        self.original=json.loads((SCRIPT.parent/'tests/fixtures/original-plan-551a3c9.json').read_text())
        self.target=copy.deepcopy(self.original)
        self.target['schema']='tenv5k-matrix-v2'
        code=self.original['code_sha256']
        names=['execution/prnet.py']+['output/api/tenv5k/'+n+'.py' for n in ('model_prnet','learner','learning','inputs','protocol','native_import')]
        self.method=dict(scientific_identity='prnet-native-f32-v1',code_sha256={n:code[n] for n in names})
        self.target['methods']={'prnet':self.method}
        self.target['shared']=dict(input_manifest_sha256=self.original['input_manifest_sha256'],
            means_manifest_sha256=self.original['means_manifest_sha256'],
            evaluation_code_sha256={n:code[n] for n in ('output/api/tenv5k/protocol.py','output/api/tenv5k/comparison.py')})
        self.audit=dict(source_plan_sha256=import_matrix.ORIGINAL_PLAN_SHA256,target_method=self.method,
                        source_runtime='551a3c9',reviewed_by='test fixture unchanged PRnet runtime',changed_code={})
        self.run=matrix.plan()[0]

    def fixture(self):
        import shutil
        import subprocess
        root=self.root/'run';root.mkdir()
        config={**self.run,'device':'cuda','smoke_rows':None,'max_epochs':1000,'patience':50,'strict_lower':True}
        for name in ('checkpoint.pt','prediction.npy','identity.csv','per_well.csv','per_drug.csv','learning_curve.csv'):
            (root/name).write_bytes(b'fixed immutable historical evidence\n')
        (root/'config.json').write_text(json.dumps(config))
        (root/'environment.json').write_text(json.dumps(dict(python='3.fixture',torch='fixture',numpy='fixture',device='cuda',cuda='fixture')))
        (root/'terminal.json').write_text(json.dumps(dict(status='complete',smoke=False,replay='exact native sampled prediction')))
        outputs=matrix._tree_identity(root)
        (root/'manifest.json').write_text(json.dumps(dict(input_manifest_sha256=self.original['input_manifest_sha256'],outputs=outputs)))
        provenance=root/'import-provenance';provenance.mkdir()
        shutil.copy2(SCRIPT.parent/'tests/fixtures/original-plan-551a3c9.json',provenance/'original-plan.json')
        (provenance/'compatibility-audit.json').write_text(json.dumps(self.audit))
        repo=SCRIPT.parents[3]
        for name in self.original['code_sha256']:
            path=provenance/'runtime'/name;path.parent.mkdir(parents=True,exist_ok=True)
            path.write_bytes(subprocess.check_output(['git','-C',str(repo),'show','551a3c9:stage2/02_models_tenv-5k/'+name]))
        replay=dict(status='exact',prediction_sha256=outputs['prediction.npy'],checkpoint_sha256=outputs['checkpoint.pt'],
                    input_manifest_sha256=self.original['input_manifest_sha256'],config_sha256=outputs['config.json'],environment_sha256=outputs['environment.json'])
        (provenance/'fresh-process-replay.json').write_text(json.dumps(replay))
        receipt=dict(schema='tenv5k-import-v1',run=self.run,origin='imported; not retrained',
            source_plan_sha256=self.imports.ORIGINAL_PLAN_SHA256,source_runtime_commit=self.imports.FROZEN_COMMIT,
            source_manifest_sha256=matrix.digest(root/'manifest.json'),source_config_sha256=outputs['config.json'],
            source_environment=json.loads((root/'environment.json').read_text()),target_method=self.method,
            target_shared=self.target['shared'],artifacts_sha256=matrix._tree_identity(root))
        (root/'import-receipt.json').write_text(json.dumps(receipt))
        return root

    def test_import_metadata_preserved_and_portable_verification_accepts(self):
        root=self.fixture()
        original={n:(root/n).read_bytes() for n in ('config.json','terminal.json','manifest.json')}
        result=self.imports.verify_import(root,self.run,self.target,matrix)
        self.assertEqual(result['origin'],'imported; not retrained')
        self.assertEqual({n:(root/n).read_bytes() for n in original},original)
        import shutil
        portable=self.root/'published-copy';shutil.copytree(root,portable)
        self.assertEqual(self.imports.verify_import(portable,self.run,self.target,matrix)['origin'],'imported; not retrained')
        self.assertEqual({n:(portable/n).read_bytes() for n in original},original)
        with self.assertRaisesRegex(ValueError,'relabelled'):
            matrix.record_execution(root,self.run,self.target)

    def test_compatibility_rejects_forged_plan_shared_scoring_adapter_and_unaudited_learner(self):
        import copy
        self.imports.validate_compatibility(self.original,self.target,self.audit,self.imports.ORIGINAL_PLAN_SHA256)
        for case in ('plan','input','scoring','adapter','learner'):
            target=copy.deepcopy(self.target)
            expected=self.imports.ORIGINAL_PLAN_SHA256
            if case=='plan':expected='f'*64
            elif case=='input':target['input_manifest_sha256']='other input'
            elif case=='scoring':target['shared']['evaluation_code_sha256']['output/api/tenv5k/protocol.py']='other scoring'
            elif case=='adapter':target['methods']['prnet']['code_sha256']['output/api/tenv5k/model_prnet.py']='other model'
            else:target['methods']['prnet']['code_sha256']['output/api/tenv5k/learner.py']='unreviewed learner'
            with self.subTest(case=case), self.assertRaises(ValueError):
                self.imports.validate_compatibility(self.original,target,self.audit,expected)

    def test_wrong_identity_and_artifact_provenance_receipt_tampering_rejected(self):
        root=self.fixture()
        for name in ('checkpoint.pt','prediction.npy','config.json','terminal.json','manifest.json',
                     'import-provenance/runtime/output/api/tenv5k/learner.py','import-provenance/fresh-process-replay.json'):
            path=root/name;before=path.read_bytes();path.write_bytes(before+b'changed')
            with self.subTest(name=name),self.assertRaises((ValueError,json.JSONDecodeError)):
                self.imports.verify_import(root,self.run,self.target,matrix)
            path.write_bytes(before)
        receipt=json.loads((root/'import-receipt.json').read_text())
        receipt['target_method']['scientific_identity']='renamed'
        (root/'import-receipt.json').write_text(json.dumps(receipt))
        with self.assertRaises(ValueError):self.imports.verify_import(root,self.run,self.target,matrix)
        for run in (dict(adapter='prnet',scenario='warm-start',fold=0,seed=17),dict(adapter='prnet',scenario='cold-drug',fold=9,seed=17)):
            with self.subTest(run=run),self.assertRaises(ValueError):self.imports.verify_import(root,run,self.target,matrix)

    def test_replay_prediction_and_environment_mismatch_rejected_even_with_rehashed_receipt(self):
        root=self.fixture()
        replay_path=root/'import-provenance/fresh-process-replay.json'
        for field in ('prediction_sha256','environment_sha256','config_sha256'):
            replay=json.loads(replay_path.read_text());before=replay.copy();replay[field]='mismatch'
            replay_path.write_text(json.dumps(replay))
            receipt=json.loads((root/'import-receipt.json').read_text())
            receipt['artifacts_sha256']['import-provenance/fresh-process-replay.json']=matrix.digest(replay_path)
            (root/'import-receipt.json').write_text(json.dumps(receipt))
            with self.subTest(field=field),self.assertRaises(ValueError):self.imports.verify_import(root,self.run,self.target,matrix)
            replay_path.write_text(json.dumps(before))

    def test_original_runtime_git_identity_and_disk_bytes_fail_closed(self):
        from unittest.mock import patch
        runtime=self.root/'runtime';runtime.mkdir()
        with patch('subprocess.check_output',side_effect=['wrong-commit\n','551a3c9-original\n']):
            with self.assertRaisesRegex(ValueError,'frozen runtime'):
                self.imports.verify_runtime(runtime,self.original,matrix.digest)
        path=runtime/'stage2/02_models_tenv-5k/execution/prnet.py'
        path.parent.mkdir(parents=True);path.write_text('changed frozen runtime')
        with patch('subprocess.check_output',return_value=self.imports.FROZEN_COMMIT+'\n'):
            with self.assertRaisesRegex(ValueError,'runtime bytes changed'):
                self.imports.verify_runtime(runtime,{'code_sha256':{'execution/prnet.py':self.original['code_sha256']['execution/prnet.py']}},matrix.digest)

    def test_new_training_requires_preregistered_identity_and_keeps_origin_distinct(self):
        import shutil
        root=self.fixture()
        (root/'import-receipt.json').unlink();shutil.rmtree(root/'import-provenance')
        with self.assertRaisesRegex(ValueError,'preregistered'):
            matrix.record_execution(root,self.run,self.target)
        registration=matrix.execution_registration(self.run,self.target,'1'*64)
        from tenv5k.attempt_identity import attempt_binding
        import torch
        binding=attempt_binding(registration)
        config=json.loads((root/'config.json').read_text());config['execution_binding']=binding
        current=dict(method_identity='prnet-native-f32-lr1e4-v1',compute_dtype='float32',
                     initialization_dtype='float32',prediction_dtype='float32',learning_rate=1e-4)
        config.update(current)
        self.target['methods']['prnet']['scientific_identity']=current['method_identity']
        self.target['methods']['prnet']['scientific_configuration']=current
        # Updating the method requires a newly registered attempt, before binding.
        registration=matrix.execution_registration(self.run,self.target,'1'*64)
        binding=attempt_binding(registration);config['execution_binding']=binding
        (root/'config.json').write_text(json.dumps(config))
        torch.save(dict(config=config,execution_binding=binding,state_dict={'weight':torch.ones(1)}),root/'checkpoint.pt')
        import numpy as np
        np.save(root/'prediction.npy',np.ones((1,3),dtype=np.float32))
        manifest=json.loads((root/'manifest.json').read_text());manifest['execution_binding']=binding
        for name in ('config.json','checkpoint.pt','prediction.npy'):manifest['outputs'][name]=matrix.digest(root/name)
        (root/'manifest.json').write_text(json.dumps(manifest))
        matrix.record_execution(root,self.run,self.target,registration)
        self.assertEqual(matrix.verify_execution(root,self.run,self.target)['origin'],'trained')
        receipt=json.loads((root/'execution-receipt.json').read_text())
        receipt['registration']['runtime_code_sha256']['execution/prnet.py']='foreign runtime'
        (root/'execution-receipt.json').write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError,'execution identity'):
            matrix.verify_execution(root,self.run,self.target)

    def test_registered_slot_cannot_recover_historical_complete_run_as_fresh(self):
        import shutil
        root=self.fixture()
        (root/'import-receipt.json').unlink();shutil.rmtree(root/'import-provenance')
        candidate=self.root/'matrix';target=matrix.run_path(candidate,self.run)
        target.parent.mkdir(parents=True);shutil.move(str(root),target)
        registration=matrix.execution_registration(self.run,self.target,'2'*64)
        slot=matrix.registration_path(candidate,self.run);slot.parent.mkdir(parents=True)
        slot.write_text(json.dumps(registration,indent=2)+'\n')
        args=SimpleNamespace(candidate=candidate,device='cuda')
        with self.assertRaisesRegex(ValueError,'attempt|registration'):
            action=matrix.recover_existing(args,self.run,self.original['input_manifest_sha256'],lambda:None)
            if action=='complete':matrix.record_execution(target,self.run,self.target,registration)

    def test_new_checkpoint_binding_cannot_be_replaced_by_legacy_checkpoint(self):
        import shutil
        import torch
        from tenv5k.attempt_identity import attempt_binding,verify_artifact_binding
        root=self.fixture();(root/'import-receipt.json').unlink();shutil.rmtree(root/'import-provenance')
        registration=matrix.execution_registration(self.run,self.target,'3'*64)
        binding=attempt_binding(registration)
        config=json.loads((root/'config.json').read_text());config['execution_binding']=binding
        (root/'config.json').write_text(json.dumps(config))
        # Config and manifest bind the new attempt, but checkpoint is historical.
        torch.save(dict(config={k:v for k,v in config.items() if k!='execution_binding'},state_dict={}),root/'checkpoint.pt')
        manifest=json.loads((root/'manifest.json').read_text());manifest['execution_binding']=binding
        for name in ('config.json','checkpoint.pt'):manifest['outputs'][name]=matrix.digest(root/name)
        (root/'manifest.json').write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError,'checkpoint attempt'):
            matrix.record_execution(root,self.run,self.target,registration)
        torch.save(dict(config=config,execution_binding=binding,state_dict={}),root/'checkpoint.pt')
        manifest['outputs']['checkpoint.pt']=matrix.digest(root/'checkpoint.pt')
        (root/'manifest.json').write_text(json.dumps(manifest))
        self.assertEqual(verify_artifact_binding(root,registration),binding)
        manifest['execution_binding']={**binding,'attempt_nonce':'4'*64}
        (root/'manifest.json').write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError,'manifest attempt'):
            verify_artifact_binding(root,registration)

    def test_current_method_cannot_use_historical_metadata_import_authorization(self):
        import copy
        root=self.root/'inputs';root.mkdir();(root/'manifest.json').write_text('{}')
        args=SimpleNamespace(inputs=root,means=root,device='cuda',prnet_source=root,transigen_source=root)
        target=copy.deepcopy(self.target);target['methods']['prnet']=matrix.frozen_plan(args)['methods']['prnet']
        audit=copy.deepcopy(self.audit);audit['target_method']=target['methods']['prnet']
        paths={'output/api/tenv5k/learner.py','execution/prnet.py','output/api/tenv5k/attempt_identity.py'}
        audit['changed_code']={p:'Reviewed metadata-only binding, no model RNG or mathematics change' for p in paths}
        with self.assertRaisesRegex(ValueError,'scientific identity'):
            self.imports.validate_compatibility(self.original,target,audit,self.imports.ORIGINAL_PLAN_SHA256)
        target['methods']['prnet']['code_sha256']['execution/prnet.py']='arbitrary CLI change'
        with self.assertRaisesRegex(ValueError,'scientific identity'):
            self.imports.validate_compatibility(self.original,target,audit,self.imports.ORIGINAL_PLAN_SHA256)

    def test_failed_attempt_keeps_its_registration_and_original_bytes(self):
        candidate=self.root/'failed-matrix';target=matrix.run_path(candidate,self.run);target.mkdir(parents=True)
        (target/'terminal.json').write_text(json.dumps(dict(status='failed')))
        (target/'partial.bin').write_bytes(b'original interrupted attempt')
        registration=matrix.execution_registration(self.run,self.target,'5'*64)
        args=SimpleNamespace(candidate=candidate,device='cuda')
        self.assertEqual(matrix.recover_existing(args,self.run,self.original['input_manifest_sha256'],lambda:None,registration),'train')
        archived=next((candidate/'attempts').rglob('execution-registration.json'))
        self.assertEqual(json.loads(archived.read_text()),registration)
        self.assertEqual((archived.parent/'partial.bin').read_bytes(),b'original interrupted attempt')


if __name__ == '__main__':
    unittest.main()
