"""Joint reconciliation public seam, without models or datasets."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))
from test_perturbdiff_v2_migration import Reconciliation
import joint_v2_migration as migration

class JointReconciliation(unittest.TestCase):
    def test_both_models_all_attempt_walls_and_orphan_rejection(self):
        with tempfile.TemporaryDirectory() as d:
            queue, stopped, output = Reconciliation().fixture(Path(d))
            queue['jobs'] += [{'model': 'cmonge', 'fold': f, 'seed': s}
                              for f in range(3) for s in (17,42,2026)]
            path = output / 'scheduling/cmonge_fold0_seed17.json'
            path.write_text(json.dumps(dict(status='completed',training_seconds=9,
                attempts=[dict(wall_seconds=7),dict(wall_seconds=5)])))
            run = output / 'runs/cmonge/fold_0/seed_17'
            run.mkdir(parents=True)
            (run / 'native.msgpack').write_bytes(b'preserve me')
            ledger, _ = migration.reconcile(queue, stopped, output)
            self.assertEqual(len(ledger),18)
            self.assertEqual(ledger['cmonge_fold0_seed17']['consumed_seconds'],12)
            self.assertEqual(ledger['cmonge_fold0_seed17']['remaining_seconds'],14388)
            self.assertEqual(ledger['cmonge_fold2_seed2026']['consumed_seconds'],0)
            (output/'scheduling/renamed.json').write_text(json.dumps({'run':'runs/cmonge/fold_2/seed_2026'}))
            with self.assertRaisesRegex(ValueError,'contradicted'):
                migration.reconcile(queue,stopped,output)

if __name__ == '__main__': unittest.main()
