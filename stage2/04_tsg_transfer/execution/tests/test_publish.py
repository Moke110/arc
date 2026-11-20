"""Check the published comparison against a worked three-fold example."""
import json
import sys
from tempfile import TemporaryDirectory
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from publish import comparison, verify_report_files
from report_resources import verify_phase_output


class ComparisonTest(unittest.TestCase):
    def test_fold_means_weight_folds_equally_and_reject_unpaired_identities(self):
        values = ((0.1, -0.2, 0.4), (0.3, 0.2, 0.1), (0.5, 0.4, 0.6))
        receipts = [
            dict(fold=fold, path=path, files={'identity.csv': f'identities-{fold}'},
                 metrics=dict(pcc=score, n_samples=3 ** fold, undefined_pcc_count=0))
            for fold, scores in enumerate(values)
            for path, score in zip(('tenv-only', 'tahoe-only', 'transfer'), scores, strict=True)
        ]
        scores, differences = comparison(receipts)
        means = {row['path']: row['pcc'] for row in scores if row['fold'] == 'mean'}
        self.assertAlmostEqual(means['tenv-only'], 0.3)
        self.assertAlmostEqual(means['tahoe-only'], 0.4 / 3)
        self.assertAlmostEqual(means['transfer'], 1.1 / 3)
        self.assertAlmostEqual(differences[-1]['transfer_minus_tenv_only'], 0.2 / 3)
        self.assertAlmostEqual(differences[-1]['transfer_minus_tahoe_only'], 0.7 / 3)
        receipts[-1]['files']['identity.csv'] = 'different'
        with self.assertRaisesRegex(ValueError, 'ordered test identities'):
            comparison(receipts)


    def test_report_rejects_unexpected_archive_entries(self):
        with TemporaryDirectory() as temporary:
            stage, archive = (Path(temporary) / name for name in ('stage', 'archive'))
            for folder in (stage, archive):
                folder.mkdir()
                (folder / 'receipt').write_text('equal')
            verify_report_files(stage, archive)
            (archive / 'unverified').mkdir()
            with self.assertRaisesRegex(ValueError, 'inventory differs'):
                verify_report_files(stage, archive)

    def test_export_and_replay_are_bound_after_candidate_relocation(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            attempt = root / 'attempt'
            attempt.mkdir()
            candidate = root / 'relocated'
            relative = Path('tahoe-only/fold-0/test/manifest.json')
            (candidate / relative).parent.mkdir(parents=True)
            receipt = {'path': 'tahoe-only', 'metrics': {'pcc': 0.12}}
            (candidate / relative).write_text(json.dumps(receipt))
            (attempt / 'stdout.txt').write_text(json.dumps(str(Path('/original') / relative)))
            entry = 'python study/execution/transfer.py'
            flags = '--fold 0 --path tahoe-only --candidate /original'
            (attempt / 'command.txt').write_text(f'{entry} evaluate {flags}')
            verify_phase_output(attempt, candidate, 'tahoe-only', 'evaluate', 0)

            (attempt / 'command.txt').write_text(f'{entry} replay-test {flags}')
            (attempt / 'stdout.txt').write_text(json.dumps(receipt))
            verify_phase_output(attempt, candidate, 'tahoe-only', 'replay-test', 0)
            (attempt / 'stdout.txt').write_text(json.dumps({'path': 'transfer'}))
            with self.assertRaisesRegex(ValueError, 'stdout differs'):
                verify_phase_output(attempt, candidate, 'tahoe-only', 'replay-test', 0)
            with self.assertRaisesRegex(ValueError, 'wrong fold'):
                verify_phase_output(attempt, candidate, 'tahoe-only', 'replay-test', 1)


if __name__ == '__main__':
    unittest.main()
