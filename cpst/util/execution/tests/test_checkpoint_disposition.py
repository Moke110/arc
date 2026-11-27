"""Exercise the accepted release and prediction-archive boundaries."""
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from cpst_util.checkpoint_disposition import dispose_checkpoints, verify_prediction_archive
from cpst_util.derivation.release import ReleaseError, check_outputs, read_json, sha256, write_json


class CheckpointDispositionTests(unittest.TestCase):
    def setUp(self):
        work = TemporaryDirectory()
        self.addCleanup(work.cleanup)
        self.root = Path(work.name) / "bundle"
        self.root.mkdir()
        self.acceptance = write_json(self.root / 'acceptance.json', {'result': 'PASS'})
        self.receipt = self.root / 'disposition.json'
        self.runs = ('runs/new/manifest.json',)
        for run in ('new', 'deferred', 'historical'):
            folder = self.root / 'runs' / run
            folder.mkdir(parents=True)
            for name in ('checkpoint.pt', 'prediction.npy'):
                (folder / name).write_bytes(f'{run} {name}'.encode())
            write_json(folder / 'manifest.json', {'outputs': {
                name: sha256(folder / name) for name in ('checkpoint.pt', 'prediction.npy')
            }})

    def seal(self):
        return {name: sha256(self.root / name) for name in self.runs}

    def dispose(self, *, verify=lambda: None, retained=lambda entries: None):
        return dispose_checkpoints(self.root, receipt=self.receipt.name,
                                   acceptance=self.acceptance,
                                   acceptance_sha256=sha256(self.acceptance),
                                   runs=self.seal(), verify_full=verify, verify_retained=retained)

    def verify(self, pin):
        return verify_prediction_archive(self.root, receipt=self.receipt.name,
                                         receipt_sha256=pin, acceptance=self.acceptance,
                                         acceptance_sha256=sha256(self.acceptance), runs=self.seal())

    def test_only_authorized_weight_is_removed_and_default_verify_fails(self):
        manifest = self.root / self.runs[0]
        before = manifest.read_bytes()
        pin = self.dispose()
        self.assertEqual(manifest.read_bytes(), before)
        self.assertFalse((manifest.parent / 'checkpoint.pt').exists())
        self.assertTrue((self.root / 'runs/historical/checkpoint.pt').exists())
        with self.assertRaisesRegex(ReleaseError, 'missing release payload'):
            check_outputs(manifest.parent, read_json(manifest))
        self.assertEqual(self.verify(pin)['status'], 'complete')
        (manifest.parent / 'prediction.npy').unlink()
        with self.assertRaisesRegex(ReleaseError, 'missing release payload'):
            self.verify(pin)

    def test_missing_weight_fails_before_acceptance(self):
        (self.root / 'runs/new/checkpoint.pt').unlink()
        with self.assertRaisesRegex(ReleaseError, 'missing release payload'):
            self.dispose()
        self.assertFalse(self.receipt.exists())

    def test_failed_full_gate_preserves_files(self):
        def reject():
            raise ReleaseError('full gate failed')
        with self.assertRaisesRegex(ReleaseError, 'full gate failed'):
            self.dispose(verify=reject)
        self.assertTrue((self.root / 'runs/new/checkpoint.pt').exists())

    def test_pending_receipt_and_altered_pin_fail(self):
        pin = self.dispose()
        data = read_json(self.receipt)
        data['status'] = 'pending'
        write_json(self.receipt, data)
        with self.assertRaisesRegex(ReleaseError, 'receipt digest'):
            self.verify(pin)
        with self.assertRaisesRegex(ReleaseError, 'complete'):
            self.verify(sha256(self.receipt))

    def test_symbolic_checkpoint_cannot_be_disposed(self):
        weight = self.root / 'runs/new/checkpoint.pt'
        alias = weight.with_name('other.pt')
        weight.rename(alias)
        weight.symlink_to(alias.name)
        with self.assertRaisesRegex(ReleaseError, 'symbolic link'):
            self.dispose()
        self.assertTrue(weight.is_symlink())

    def test_traversal_manifest_cannot_be_authorized(self):
        self.runs = ('runs/new/../historical/manifest.json',)
        with self.assertRaisesRegex(ReleaseError, 'invalid disposition path'):
            self.dispose()
        self.assertTrue((self.root / 'runs/historical/checkpoint.pt').exists())

    def test_invalid_retained_digest_prevents_removal(self):
        (self.root / 'runs/new/prediction.npy').write_bytes(b'corrupted')
        with self.assertRaisesRegex(ReleaseError, 'digest differs'):
            self.dispose()
        self.assertTrue((self.root / 'runs/new/checkpoint.pt').exists())

    def test_failed_retained_verification_restores_checkpoint(self):
        weight = self.root / 'runs/new/checkpoint.pt'
        before = weight.read_bytes()

        def reject(entries):
            self.assertFalse(weight.exists())
            self.assertEqual(entries[0]['size'], len(before))
            raise ReleaseError('prediction gate failed')

        with self.assertRaisesRegex(ReleaseError, 'prediction gate failed'):
            self.dispose(retained=reject)
        self.assertEqual(weight.read_bytes(), before)
        self.assertEqual(read_json(self.receipt)['status'], 'pending')
        with self.assertRaisesRegex(ReleaseError, 'complete'):
            self.verify(sha256(self.receipt))

    def test_second_batch_requires_pinned_first_batch_and_verifies_both(self):
        first_pin = self.dispose()
        previously_disposed = self.seal()
        self.receipt = self.root / 'second.json'
        self.runs = ('runs/deferred/manifest.json',)
        deferred = self.root / 'runs/deferred/checkpoint.pt'
        kwargs = dict(root=self.root, receipt=self.receipt.name,
                      acceptance=self.acceptance, acceptance_sha256=sha256(self.acceptance),
                      runs=self.seal(), verify_full=lambda: None,
                      verify_retained=lambda entries: self.assertEqual(len(entries), 2))
        with self.assertRaisesRegex(ReleaseError, 'receipt digest'):
            dispose_checkpoints(**kwargs, previous=('disposition.json', '0' * 64, previously_disposed))
        self.assertTrue(deferred.is_file())
        second_pin = dispose_checkpoints(**kwargs,
                                         previous=('disposition.json', first_pin, previously_disposed))
        self.runs += tuple(previously_disposed)
        self.assertEqual(len(self.verify(second_pin)['weights']), 2)
        self.assertEqual(sha256(self.root / 'disposition.json'), first_pin)

    def test_retry_restores_checkpoint_after_interrupted_prediction_gate(self):
        import subprocess
        import sys
        code = """
import os
from pathlib import Path
import sys
from cpst_util.checkpoint_disposition import dispose_checkpoints
from cpst_util.derivation.release import sha256

root = Path(sys.argv[1])
acceptance = root / 'acceptance.json'
name = 'runs/new/manifest.json'
dispose_checkpoints(root, receipt='disposition.json', acceptance=acceptance,
                    acceptance_sha256=sha256(acceptance), runs={name: sha256(root / name)},
                    verify_full=lambda: None, verify_retained=lambda entries: os._exit(15))
"""
        run = subprocess.run([sys.executable, '-c', code, str(self.root)], check=False)
        self.assertEqual(run.returncode, 15)
        self.assertFalse((self.root / 'runs/new/checkpoint.pt').exists())
        with self.assertRaisesRegex(ReleaseError, 'complete'):
            self.verify(sha256(self.receipt))
        verified = []
        pin = self.dispose(verify=lambda: verified.append(
            (self.root / 'runs/new/checkpoint.pt').read_bytes()))
        self.assertEqual(verified, [b'new checkpoint.pt'])
        self.assertEqual(self.verify(pin)['status'], 'complete')
        self.assertFalse((self.root / 'runs/new/checkpoint.pt').exists())

    def test_repeat_completed_request_reverifies_without_requiring_removed_weights(self):
        first_pin = self.dispose()
        def reject_full():
            self.fail('full gate should not run after completed disposal')
        verified = []
        pin = self.dispose(verify=reject_full, retained=lambda entries: verified.extend(entries))
        self.assertEqual(pin, first_pin)
        self.assertEqual([entry['path'] for entry in verified], ['runs/new/checkpoint.pt'])
