"""Public process/publication boundaries required by the pretraining review."""
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
STEP=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(STEP/'output/api'))
spec=importlib.util.spec_from_file_location('matrix231_review',STEP/'execution/run_matrix.py')
matrix=importlib.util.module_from_spec(spec);spec.loader.exec_module(matrix)


def contents(root):
    return {str(p.relative_to(root)):p.read_bytes() for p in root.rglob('*') if p.is_file()}


class ReviewFixTests(unittest.TestCase):
    def test_native_source_and_sibling_import_are_byte_immutable_in_fresh_process(self):
        original=Path('/home/chang/arc/pert_model')
        if not (original/'TranSiGen/source/src/model.py').is_file(): self.skipTest('native source unavailable')
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            for model in ('PRnet','TranSiGen'):
                shutil.copytree(original/model/'source',root/model,ignore=shutil.ignore_patterns('__pycache__','*.pyc','.git'))
            before=contents(root)
            script="import sys;sys.path.insert(0,sys.argv[1]);from tenv5k.model_prnet import build_model;from tenv5k.model_transigen import Adapter;build_model(sys.argv[2],17,'cpu');Adapter(sys.argv[3])"
            environment={**os.environ,'OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'2','NUMBA_CACHE_DIR':'/tmp/arc-tenv5k-test-numba','MPLCONFIGDIR':'/tmp/arc-tenv5k-test-matplotlib'}
            environment.pop('PYTHONDONTWRITEBYTECODE',None)
            # Same native loading boundary used by fit and the fresh replay CLI.
            for _ in range(2):
                subprocess.run([sys.executable,'-c',script,str(STEP/'output/api'),str(root/'PRnet'),str(root/'TranSiGen')],check=True,env=environment)
            self.assertEqual(set(contents(root)),set(before))
            self.assertEqual(contents(root),before)

    def test_publication_preflights_every_destination_and_rolls_back_io_failure(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);step=root/'step';candidate=root/'candidate'
            for relative,value in [('runs/new','newrun'),('comparison/new','newcomparison'),('figures/metrics.png','newpng'),('comparison/metrics_values.csv','newcsv'),('comparison/paired_differences_values.csv','newpair'),('comparison/learning_curve_values.csv','newcurve'),('resources/run.time','newtime'),('acceptance.json','{}'),('plan.json','{}')]:
                path=candidate/relative;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(value)
            for relative,value in [('output/comparison/old','oldcomparison'),('figure/data/metrics_values.csv','oldcsv'),('figure/metrics.png','oldpng'),('log/231-acceptance.json','oldaccept')]:
                path=step/relative;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(value)
            before=contents(step);args=SimpleNamespace(candidate=candidate)
            with patch.object(matrix,'STEP',step):
                with self.assertRaises(ValueError):matrix.publish(args)
                self.assertEqual(contents(step),before)
                shutil.rmtree(step/'output/comparison')
                before=contents(step)
                copytree=shutil.copytree
                def failing_copy(source,target,*a,**kw):
                    if Path(source)==candidate/'comparison':raise OSError('injected staging copy failure')
                    return copytree(source,target,*a,**kw)
                with patch.object(matrix.shutil,'copytree',side_effect=failing_copy):
                    with self.assertRaises(OSError):matrix.publish(args)
                self.assertEqual(contents(step),before)
                replace=os.replace
                calls=[]
                def fail_second_rename(source,target):
                    calls.append(str(target))
                    if len(calls)==2:raise OSError('injected second rename failure')
                    return replace(source,target)
                with patch.object(matrix.os,'replace',side_effect=fail_second_rename):
                    with self.assertRaises(OSError):matrix.publish(args)
                self.assertEqual(contents(step),before)
                def failing_replace(source,target):
                    if Path(target)==step/'figure/metrics.png':raise OSError('injected publication interruption')
                    return replace(source,target)
                with patch.object(matrix.os,'replace',side_effect=failing_replace):
                    with self.assertRaises(OSError):matrix.publish(args)
                self.assertEqual(contents(step),before)
                rename=os.rename
                def failing_rollback_rename(source,target):
                    if Path(target)==step/'figure/metrics.png':raise OSError('injected rollback I/O failure')
                    return rename(source,target)
                with patch.object(matrix.os,'replace',side_effect=failing_replace),patch.object(matrix.os,'rename',side_effect=failing_rollback_rename):
                    with self.assertRaises((OSError,RuntimeError)):matrix.publish(args)
                journal=root/'.tenv5k-publication-journal.json'
                self.assertTrue(journal.exists())
                state=json.loads(journal.read_text())
                prior_png=next(e for e in state['entries'] if Path(e['target'])==step/'figure/metrics.png')
                self.assertTrue(Path(state['transaction']).is_dir())
                self.assertEqual(Path(prior_png['backup']).read_bytes(),b'oldpng')
                # Next public attempt restores old bytes before staging anything.
                def verify_restored_then_fail_copy(source,target,*a,**kw):
                    self.assertEqual(contents(step),before)
                    raise OSError('stop after proving recovery restored old release')
                with patch.object(matrix.shutil,'copytree',side_effect=verify_restored_then_fail_copy):
                    with self.assertRaises(OSError):matrix.publish(args)
                self.assertEqual(contents(step),before)
                self.assertFalse(journal.exists())
                matrix.publish(args)
                self.assertEqual((step/'output/runs/new').read_text(),'newrun')
                self.assertEqual((step/'figure/metrics.png').read_text(),'newpng')

    def test_interrupted_publication_restores_only_journaled_destinations(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);step=root/'step';step.mkdir();candidate=root/'candidate';candidate.mkdir()
            (candidate/'acceptance.json').write_text('{}')
            transaction=root/'transaction';transaction.mkdir()
            target=step/'figure/metrics.png';target.parent.mkdir();target.write_bytes(b'new interrupted png')
            backup=transaction/'backup';backup.write_bytes(b'old accepted png')
            untouched=step/'unrelated';untouched.write_bytes(b'do not touch')
            journal=root/'.tenv5k-publication-journal.json'
            journal.write_text(json.dumps(dict(transaction=str(transaction),entries=[dict(target=str(target),backup=str(backup),staged=str(transaction/'already-renamed'),existed=True)])))
            with patch.object(matrix,'STEP',step):
                with self.assertRaises((ValueError,FileNotFoundError)):matrix.publish(SimpleNamespace(candidate=candidate))
            self.assertEqual(target.read_bytes(),b'old accepted png')
            self.assertEqual(untouched.read_bytes(),b'do not touch')
            self.assertFalse(journal.exists())

    def test_memory_preflight_refuses_before_heavy_command(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            free=root/'free';free.write_text('#!/bin/sh\necho "Mem: 47 40 0 0 7 11"\n');free.chmod(0o755)
            marker=root/'launched'
            result=subprocess.run(['bash',str(STEP/'execution/capped.sh'),'touch',str(marker)],env={**os.environ,'PATH':str(root)+':'+os.environ['PATH']},capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertFalse(marker.exists());self.assertIn('12 GiB',result.stderr)

if __name__=='__main__':unittest.main()
