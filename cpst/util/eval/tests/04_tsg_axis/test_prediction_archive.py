"""Exercise the Step 04 prediction-archive output boundary on disposable runs."""
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

ROOT = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(ROOT / 'cpst/04_tsg_axis/execution'))
from prediction_archive import _retained_output_check
from cpst_util.derivation.release import ReleaseError, read_json, sha256, write_json


class PredictionArchiveTests(unittest.TestCase):
    def test_only_the_pinned_run_can_omit_its_checkpoint(self):
        with TemporaryDirectory() as scratch:
            matrix = Path(scratch)
            for name in ('allowed', 'other'):
                run = matrix / name
                run.mkdir()
                (run / 'prediction.npy').write_bytes(b'prediction')
                write_json(run / 'manifest.json', {'outputs': {
                    'checkpoint.pt': 'a' * 64,
                    'prediction.npy': sha256(run / 'prediction.npy'),
                }})
            manifest = matrix / 'allowed/manifest.json'
            check = _retained_output_check(matrix, {'allowed/manifest.json': sha256(manifest)})
            check(manifest.parent, read_json(manifest))
            with self.assertRaisesRegex(ReleaseError, 'checkpoint.pt'):
                check(matrix / 'other', read_json(matrix / 'other/manifest.json'))
            (manifest.parent / 'prediction.npy').unlink()
            with self.assertRaisesRegex(ReleaseError, 'prediction.npy'):
                check(manifest.parent, read_json(manifest))
            write_json(manifest, {'outputs': {'checkpoint.pt': 'b' * 64}})
            with self.assertRaisesRegex(ValueError, 'manifest'):
                check(manifest.parent, read_json(manifest))

    def test_archive_adapter_checks_full_then_recomputes_retained_report(self):
        from argparse import Namespace
        from unittest.mock import patch
        import prediction_archive as archive
        from cpst_util.derivation.release import check_outputs

        with TemporaryDirectory() as scratch:
            work = Path(scratch)
            matrix = work / 'matrix'
            run = matrix / 'runs/axis/fold-0/seed-17'
            run.mkdir(parents=True)
            manifest = run / 'manifest.json'
            weight = run / 'checkpoint.pt'
            weight.write_bytes(b'weights')
            (run / 'prediction.npy').write_bytes(b'predictions')
            write_json(manifest, {'outputs': {name: sha256(run / name) for name in ('checkpoint.pt', 'prediction.npy')}})
            name = manifest.relative_to(matrix).as_posix()
            score = work / 'score'
            result = {'matrix_manifest_sha256': sha256(write_json(matrix / 'manifest.json', {})),
                      'report_manifest_sha256': sha256(write_json(score / 'manifest.json', {})),
                      'csv_only_figure_sha256': {'metric.png': 'a' * 64}}
            acceptance = write_json(work / 'acceptance.json', dict(
                result, result='PASS', fit_source_commit=archive.FIT_SOURCE,
                matrix_fits=[{'axis': 'axis', 'fold': 0, 'seed': 17, 'run_manifest_sha256': sha256(manifest)}],
            ))
            args = Namespace(matrix=matrix, report=score, run_manifest=[name])
            calls = []

            def report(_args, *, archive_check):
                calls.append(weight.is_file())
                (archive_check or check_outputs)(run, read_json(manifest))
                return result

            with (patch.object(archive, 'ACCEPTED', acceptance),
                  patch.object(archive, 'ACCEPTED_SHA256', sha256(acceptance)),
                  patch.object(archive, 'run_paths', return_value={('axis', 0, 17): run}),
                  patch.object(archive, 'report', side_effect=report)):
                disposed = archive.dispose(args)
                self.assertEqual(calls, [True, False])
                args.disposition_sha256 = disposed['disposition_sha256']
                verified = archive.verify(args)
                self.assertEqual(calls, [True, False, False])
                self.assertEqual(verified['verification'], 'predictions-and-figures')
                self.assertEqual(verified['disposed_weights'], 1)
                args.disposition_sha256 = '0' * 64
                with self.assertRaisesRegex(ReleaseError, 'receipt digest'):
                    archive.verify(args)
