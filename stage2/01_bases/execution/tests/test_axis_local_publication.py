"""Approved candidate publication seam: failure preserves accepted products."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

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
