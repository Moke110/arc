"""Approved candidate publication seam: failure preserves accepted products."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest import mock

PATH = Path(__file__).resolve().parents[1] / 'publish_axis_local.py'

class PublicationTests(unittest.TestCase):
    def test_invalid_candidate_leaves_all_accepted_views_untouched(self):
        spec = importlib.util.spec_from_file_location('axis_local_publish', PATH)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / 'stage2/01_bases/output'
            for size in (2, 3, 4, 5):
                path = output / f'view-{size}k'
                path.mkdir(parents=True)
                (path / 'accepted').write_text(str(size))
            candidate = root / 'invalid'
            candidate.mkdir()
            with self.assertRaises((ValueError, FileNotFoundError)):
                module.install(root, candidate, root / 'recovery')
            for size in (2, 3, 4, 5):
                self.assertEqual((output / f'view-{size}k/accepted').read_text(), str(size))


class TransactionFailureTests(unittest.TestCase):
    def exercise(self, failure=None):
        spec = importlib.util.spec_from_file_location('axis_local_transaction', PATH)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / module.STEP / 'output'
            for size in (2, 3, 4, 5):
                old = output / f'view-{size}k'
                old.mkdir(parents=True)
                (old / 'accepted').write_bytes(f'accepted-{size}'.encode())
                (old / 'nested').mkdir()
                (old / 'nested/science').write_bytes(bytes([size]) * 41)
            candidate = root / 'candidate'
            candidate.mkdir()
            (candidate / 'manifest.json').write_text('{}')
            (candidate / 'new-payload').write_bytes(b'complete paired candidate')
            recovery = root / 'recovery'
            staged = recovery / 'new-view-5k'
            before = {p.relative_to(output).as_posix(): p.read_bytes()
                      for p in output.rglob('*') if p.is_file()}
            original_rename = Path.rename
            original_copy = module.shutil.copytree

            def rename(source, target):
                if ((failure in ('move-2', 'move-3', 'move-4', 'move-5')
                     and source == output / ('view-' + failure[-1] + 'k'))
                        or (failure == 'promote' and source == staged)):
                    raise OSError('injected transaction failure: ' + failure)
                return original_rename(source, target)

            def validate(root_, path):
                # Scientific validation is separately tested at the real public
                # view boundary; here inject faults at its filesystem boundary.
                if ((failure == 'staged-validation' and path == staged)
                        or (failure == 'installed-validation' and path == output / 'view-5k')):
                    raise ValueError('injected transaction validation failure')
                return {'result': 'PASS'}

            def copy(source, target):
                if failure == 'copy':
                    target.mkdir()
                    (target / 'partial').write_bytes(b'partial outside output')
                    raise OSError('injected transaction copy failure')
                return original_copy(source, target)

            with mock.patch.object(Path, 'rename', rename), mock.patch.object(module, 'validate', validate), mock.patch.object(module.shutil, 'copytree', copy):
                if failure:
                    with self.assertRaises((OSError, ValueError)):
                        module.install(root, candidate, recovery)
                else:
                    module.install(root, candidate, recovery)
            after = {p.relative_to(output).as_posix(): p.read_bytes()
                     for p in output.rglob('*') if p.is_file()}
            if failure:
                self.assertEqual(after, before, 'every original byte must survive rollback')
                self.assertEqual({p.name for p in output.iterdir()},
                                 {f'view-{size}k' for size in (2, 3, 4, 5)})
            else:
                self.assertEqual({p.name for p in output.iterdir()}, {'view-5k'})
                self.assertEqual(after, {'view-5k/' + p.name: p.read_bytes()
                                        for p in candidate.iterdir()})
                for size in (2, 3, 4, 5):
                    self.assertEqual((recovery / f'view-{size}k/accepted').read_bytes(),
                                     f'accepted-{size}'.encode())

    def test_success_installs_only_complete_candidate_and_retains_prior_bytes(self):
        self.exercise()


for _failure in ('copy', 'staged-validation', 'move-2', 'move-3', 'move-4', 'move-5',
                 'promote', 'installed-validation'):
    def _test(self, failure=_failure):
        self.exercise(failure)
    setattr(TransactionFailureTests, 'test_failure_at_' + _failure.replace('-', '_'), _test)


class FuturePreparationTests(unittest.TestCase):
    def test_both_future_entrypoints_reject_historical_generation_before_writing(self):
        import json
        repo = Path(__file__).resolve().parents[4]
        for step in ('02_models_tenv-5k', '03_tsg_dose'):
            path = repo / 'stage2' / step / 'execution/axis_local/prepare.py'
            spec = importlib.util.spec_from_file_location('future_' + step, path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            with tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                old = root / 'historical'
                old.mkdir()
                (old / 'manifest.json').write_text(json.dumps({'schema': 'stage2-view-v2'}))
                with self.assertRaisesRegex(ValueError, 'schema'):
                    module.prepare(old, root / 'new-inputs')
                self.assertFalse((root / 'new-inputs').exists())

if __name__ == '__main__':
    unittest.main()
