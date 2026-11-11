"""Sealed-envelope swap safety for the truth successor re-issue (#157 review S3).

A failed candidate must leave the previously accepted envelope usable: the
re-issue decides only after both the create and the verify commands return 0, so
a failure must not publish, must not touch the stable bytes and must not leave
partial candidate bytes at the stable path.  These tests pin that seam with a
fake runner and tiny tmpdir fixtures -- no payload, no R, no h5ad.

The failure is injected through ``main(argv, runner=...)``; production passes the
real subprocess runner and ``sys.argv``.
"""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
import reissue_truth_successor as reissue


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def snapshot(root):
    """Every file under ``root`` -> digest; detects changed, added or lost bytes."""
    return {str(path.relative_to(root)): digest(path)
            for path in sorted(root.rglob('*')) if path.is_file()}


class CandidateFailure(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.step = self.root / 'step'
        output = self.step / 'output'
        self.stable = output / 'tenv-truth-successor'
        self.stable.mkdir(parents=True)
        (output / 'training-corpora').mkdir()
        (self.step / 'log/records').mkdir(parents=True)
        # A stable envelope whose schema ``check_successor()`` rejects, so the
        # re-issue proceeds to build a candidate.  Its bytes are the thing the
        # test protects: they stand in for the accepted products.
        (self.stable / 'manifest.json').write_text(json.dumps(
            {'schema': 'accepted-stable-envelope', 'bindings': {},
             'worker_sha256': 'a' * 64}))
        (self.stable / 'deseq2_contrasts.R').write_text('print("stable")\n')
        (self.stable / 'proof.json').write_text('{"accepted": true}\n')
        (self.step / 'log/records/truth_complete_verify.json').write_text(
            '{"result": "PASS"}\n')
        self.candidate = self.root / 'candidate'
        self.record = self.root / 'record.json'
        self.seals = self.step / 'log/historical/seals'

    def run_with(self, runner):
        argv = ['--step-root', str(self.step), '--candidate', str(self.candidate),
                '--record', str(self.record)]
        return reissue.main(argv, runner=runner)

    def assert_failed_without_swap(self, before, manifest_before):
        record = json.loads(self.record.read_text())
        self.assertEqual(record['result'], 'CANDIDATE_FAILED')
        self.assertEqual(record['stable_untouched'], manifest_before)
        # No byte changed, no partial candidate appeared at the stable path.
        self.assertEqual(snapshot(self.stable), before)
        self.assertEqual(digest(self.stable / 'manifest.json'), manifest_before)
        self.assertEqual((self.stable / 'deseq2_contrasts.R').read_text(),
                         'print("stable")\n')
        # The swap (preserve old, rename candidate) happens only after PASS.
        self.assertFalse(self.seals.exists())

    def test_verify_failure_leaves_accepted_envelope_untouched(self):
        before = snapshot(self.stable)
        manifest_before = digest(self.stable / 'manifest.json')

        def runner(command, record, label, output_path=None):
            record.setdefault('commands', []).append({'label': label})
            if label == 'create candidate':
                candidate = Path(command[command.index('--output') + 1])
                candidate.mkdir(parents=True)
                (candidate / 'manifest.json').write_text('{"candidate": true}\n')
                (candidate / 'deseq2_contrasts.R').write_text('print("candidate")\n')
                return ''
            raise RuntimeError('verify candidate failed: fixture mismatch')

        self.assertEqual(self.run_with(runner), 1)
        self.assert_failed_without_swap(before, manifest_before)
        # The rejected candidate stays outside the stable path for inspection.
        self.assertTrue(self.candidate.is_dir())

    def test_create_failure_leaves_accepted_envelope_untouched(self):
        before = snapshot(self.stable)
        manifest_before = digest(self.stable / 'manifest.json')

        def runner(command, record, label, output_path=None):
            record.setdefault('commands', []).append({'label': label})
            raise RuntimeError('create candidate failed: R worker missing')

        self.assertEqual(self.run_with(runner), 1)
        self.assert_failed_without_swap(before, manifest_before)
        self.assertFalse(self.candidate.exists())


if __name__ == '__main__':
    unittest.main()
