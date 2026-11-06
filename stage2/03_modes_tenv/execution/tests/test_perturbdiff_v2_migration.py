"""Compact stopped-boundary evidence contract; no datasets/models imported."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'script'
sys.path.insert(0, str(SCRIPT))
import perturbdiff_v2_migration as migration


class Reconciliation(unittest.TestCase):
    def fixture(self, root):
        output = root / 'output'
        (output / 'scheduling').mkdir(parents=True)
        stopped = root / 'stopped'
        (stopped / 'after').mkdir(parents=True)
        jobs = [{'model': 'perturbdiff', 'fold': f, 'seed': s}
                for f in range(3) for s in (17, 42, 2026)]
        rows = {migration.key(j): {'consumed_seconds': None} for j in jobs}
        (stopped / 'stop.json').write_text(json.dumps({'process_exit_verified': True, 'ledger': rows}))
        proof = {'process_exit_verified': True, 'stop_sha256': migration.digest(stopped / 'stop.json'),
                 'nonblocking_scheduler_lock_acquired': True, 'old_scheduler_exists': False,
                 'old_worker_exists': False, 'service': 'Restart=no\nActiveState=inactive\n'}
        (stopped / 'post-stop-proof.json').write_text(json.dumps(proof))
        (stopped / 'after/boundary.json').write_text(json.dumps({'source_state_hashes': {}}))
        return {'jobs': jobs}, stopped, output

    def test_absent_attempt_and_artifact_required_for_zero(self):
        with tempfile.TemporaryDirectory() as d:
            queue, stopped, output = self.fixture(Path(d))
            rows, _ = migration.reconcile(queue, stopped, output)
            self.assertEqual(len(rows), 9)
            self.assertTrue(all(r['consumed_seconds'] == 0 and r['evidence']['zero_debit_proof'] for r in rows.values()))
            attempt = output / 'scheduling/perturbdiff_fold1_seed17.attempt0.log'
            attempt.touch()
            with self.assertRaisesRegex(ValueError, 'contradicted'):
                migration.reconcile(queue, stopped, output)
            attempt.unlink()
            orphan = output / 'scheduling/renamed-orphan.json'
            orphan.write_text(json.dumps({'identity': {'model': 'perturbdiff', 'fold': 1, 'seed': 17}}))
            with self.assertRaisesRegex(ValueError, 'contradicted'):
                migration.reconcile(queue, stopped, output)
            orphan.unlink()
            run = output / 'runs/perturbdiff/fold_1/seed_17'
            run.mkdir(parents=True)
            with self.assertRaisesRegex(ValueError, 'contradicted'):
                migration.reconcile(queue, stopped, output)

    def test_stop_debit_preserved_and_changed_state_refused(self):
        with tempfile.TemporaryDirectory() as d:
            queue, stopped, output = self.fixture(Path(d))
            state = output / 'scheduling/perturbdiff_fold0_seed17.json'
            state.write_text(json.dumps({'status': 'completed', 'training_seconds': 12475.943569269002}))
            stop = migration.read(stopped / 'stop.json')
            stop['ledger']['perturbdiff_fold0_seed17']['consumed_seconds'] = 12591.858767883998
            (stopped / 'stop.json').write_text(json.dumps(stop))
            proof = migration.read(stopped / 'post-stop-proof.json')
            proof['stop_sha256'] = migration.digest(stopped / 'stop.json')
            (stopped / 'post-stop-proof.json').write_text(json.dumps(proof))
            (stopped / 'after/boundary.json').write_text(json.dumps({'source_state_hashes': {
                'scheduling/' + state.name: migration.digest(state)}}))
            rows, _ = migration.reconcile(queue, stopped, output)
            self.assertEqual(rows['perturbdiff_fold0_seed17']['remaining_seconds'], 1808.1412321160024)
            state.write_text('{}')
            with self.assertRaisesRegex(ValueError, 'changed'):
                migration.reconcile(queue, stopped, output)

    def test_unverified_stop_and_duplicate_slots_refused(self):
        with tempfile.TemporaryDirectory() as d:
            queue, stopped, output = self.fixture(Path(d))
            queue['jobs'][-1] = queue['jobs'][0]
            with self.assertRaisesRegex(ValueError, 'nine unique'):
                migration.reconcile(queue, stopped, output)
            proof = migration.read(stopped / 'post-stop-proof.json')
            proof['old_worker_exists'] = True
            (stopped / 'post-stop-proof.json').write_text(json.dumps(proof))
            with self.assertRaisesRegex(ValueError, 'verified stopped'):
                migration.reconcile(queue, stopped, output)


if __name__ == '__main__':
    unittest.main()
