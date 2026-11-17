"""Runtime identities include owned source, never candidate provenance artifacts."""
import importlib.util
from pathlib import Path
import shutil
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[1]/'run_matrix.py'
spec = importlib.util.spec_from_file_location('runtime_identity', SCRIPT)
matrix = importlib.util.module_from_spec(spec)
spec.loader.exec_module(matrix)


class RuntimeIdentityTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.step = Path(temporary.name)/'step'
        for directory in ('execution', 'output/api'):
            shutil.copytree(SCRIPT.parent.parent/directory, self.step/directory,
                            ignore=shutil.ignore_patterns('tests', '__pycache__', 'checkpoints'))
        inputs=Path(temporary.name)/'inputs';inputs.mkdir();(inputs/'manifest.json').write_text('{}')
        self.args=SimpleNamespace(inputs=inputs,means=inputs,device='cpu',prnet_source=inputs,
                                  transigen_source=inputs,candidate=self.step/'execution/checkpoints/matrix-v2')
        self.patch=patch.object(matrix,'STEP',self.step);self.patch.start();self.addCleanup(self.patch.stop)

    def test_archived_runtime_does_not_change_the_frozen_plan(self):
        before=matrix.frozen_plan(self.args)
        archive=self.args.candidate/'runs/cold-drug/prnet/0/17/import-provenance/runtime/execution/prnet.py'
        archive.parent.mkdir(parents=True);archive.write_text('archived original runtime\n')
        self.assertEqual(matrix.frozen_plan(self.args),before)

    def test_other_checkpoint_archives_and_arbitrary_execution_candidates_are_not_source(self):
        self.args.candidate=self.step/'execution/custom-candidate'
        before=matrix.frozen_plan(self.args)
        for relative in ('custom-candidate/import-provenance/runtime/replay.py',
                         'checkpoints/previous-matrix/import-provenance/runtime/learner.py'):
            path=self.step/'execution'/relative;path.parent.mkdir(parents=True);path.write_text('archived source bytes\n')
        self.assertEqual(matrix.frozen_plan(self.args),before)

    def test_actual_execution_and_public_api_source_changes_still_change_identity(self):
        for relative in ('execution/run_matrix.py','output/api/tenv5k/protocol.py'):
            with self.subTest(source=relative):
                before=matrix.frozen_plan(self.args)
                path=self.step/relative;original=path.read_text();path.write_text(original+'\n# real source changed\n')
                after=matrix.frozen_plan(self.args)
                self.assertNotEqual(after['code_sha256'][relative],before['code_sha256'][relative])
                path.write_text(original)
        nested=self.step/'output/api/tenv5k/helper/new_reader.py'
        nested.parent.mkdir();nested.write_text('reader_version = 1\n')
        self.assertIn('output/api/tenv5k/helper/new_reader.py',matrix.frozen_plan(self.args)['code_sha256'])

    def test_candidate_cannot_hide_real_execution_or_public_api_source(self):
        for path in (self.step,self.step.parent,self.step/'execution',self.step/'execution/run_matrix.py',
                     self.step/'output/api',self.step/'output/api/tenv5k/candidate'):
            with self.subTest(candidate=str(path)):
                self.args.candidate=path
                with self.assertRaisesRegex(ValueError,'overlaps runtime source'):
                    matrix.frozen_plan(self.args)


if __name__=='__main__':unittest.main()
