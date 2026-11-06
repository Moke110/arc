"""Transaction crash seams without models, live services or formal mutation."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
import joint_v2_apply as apply


class Transaction(unittest.TestCase):
    def fixture(self, root):
        formal, replacement = root / 'formal', root / 'replacement'
        for directory, value in ((formal, 'old'), (replacement, 'new')):
            directory.mkdir()
            for name in ('queue.json', 'commands.json', 'expected_identities.json'):
                (directory / name).write_text(value)
            (directory / 'unchanged.json').write_text('same')
        source = root / 'runs/perturbdiff'
        source.mkdir(parents=True)
        (source / 'native.pt').write_bytes(b'legacy\x00binary')
        destination = root / 'legacy/perturbdiff-v1/runs'
        journal = {'formal': str(formal), 'replacement': str(replacement),
            'old_formal': apply.inventory(formal), 'new_formal': apply.inventory(replacement),
            'moves': [{'source': str(source), 'destination': str(destination),
                       'inventory': apply.inventory(source)}], 'status': 'prepared'}
        path = root / 'journal.json'
        apply.publish(path, journal)
        return path, journal, source, destination

    def test_crashes_after_move_and_exchange_recover_without_overwrite(self):
        for crash in ('archive_move', 'exchange'):
            with self.subTest(crash=crash), tempfile.TemporaryDirectory() as directory:
                path, journal, source, destination = self.fixture(Path(directory))
                def fault(phase):
                    if phase == crash:
                        raise RuntimeError('simulated power loss')
                with patch.object(apply, 'verify_identities'):
                    with self.assertRaisesRegex(RuntimeError, 'power loss'):
                        apply.finish(path, journal, lambda: None, fault)
                    self.assertFalse(source.exists())
                    self.assertEqual((destination / 'native.pt').read_bytes(), b'legacy\x00binary')
                    recovered = apply.finish(path, json.loads(path.read_text()), lambda: None)
                self.assertEqual(recovered['status'], 'committed')
                self.assertEqual(apply.inventory(Path(journal['formal'])), journal['new_formal'])
                self.assertEqual(apply.inventory(Path(journal['replacement'])), journal['old_formal'])

    def test_existing_destination_or_changed_legacy_refused(self):
        for scenario in ('destination', 'modified'):
            with self.subTest(scenario=scenario), tempfile.TemporaryDirectory() as directory:
                path, journal, source, destination = self.fixture(Path(directory))
                if scenario == 'destination':
                    destination.mkdir(parents=True)
                else:
                    (source / 'native.pt').write_bytes(b'changed')
                with self.assertRaises(ValueError):
                    apply.finish(path, journal, lambda: None)
                self.assertTrue(source.exists())
                self.assertEqual((Path(journal['formal']) / 'queue.json').read_text(), 'old')

    def test_review_pin_rejects_changed_document(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            hashes = {}
            for name in apply.NAMES:
                (root / name).write_text('{}')
                hashes[name] = apply.digest(root / name)
            apply.publish(root / 'staging-sha256.json', hashes)
            sha = apply.digest(root / 'staging-sha256.json')
            self.assertEqual(apply.reviewed_stage(root, sha), hashes)
            (root / 'budget-ledger.json').write_text('{"reset":true}')
            with self.assertRaisesRegex(ValueError, 'document hash'):
                apply.reviewed_stage(root, sha)
            with self.assertRaisesRegex(ValueError, 'inventory SHA'):
                apply.reviewed_stage(root, '0' * 64)

    def test_host_proof_requires_expected_namespace_before_service_access(self):
        with self.assertRaisesRegex(ValueError, 'host PID namespace'):
            apply.host_stopped(Path('/nonexistent'), 'not-the-host', 'anything')


if __name__ == '__main__':
    unittest.main()
