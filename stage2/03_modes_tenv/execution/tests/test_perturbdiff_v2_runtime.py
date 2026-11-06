"""Artifact seam coverage for isolated runtime gates and terminal watchdog policy."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[1] / 'script'
sys.path.insert(0, str(SCRIPT))
import perturbdiff_v2_runtime as runtime
from perturbdiff_v2_budget import initialize, write


class Runtime(unittest.TestCase):
    def test_prepared_migration_refuses_before_identity_or_dispatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); formal = root / 'formal'; formal.mkdir()
            (formal / 'queue.json').write_text(json.dumps({'output_root': str(root)}))
            (root / 'perturbdiff-v2-apply.json').write_text(json.dumps({'status': 'prepared'}))
            with patch.object(runtime, 'run_identity', side_effect=AssertionError('identity reached')):
                with self.assertRaisesRegex(ValueError, 'committed reviewed'):
                    runtime.preflight(formal, 'reviewed')

    def test_changed_committed_generation_refuses(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); formal = root / 'formal'; formal.mkdir()
            (formal / 'queue.json').write_text(json.dumps({'output_root': str(root)}))
            journal = dict(status='committed', reviewed_sha256='reviewed', formal=str(formal),
                           new_formal=runtime.inventory(formal))
            (root / 'perturbdiff-v2-apply.json').write_text(json.dumps(journal))
            (formal / 'budget-ledger.json').write_text('{}')
            with self.assertRaisesRegex(ValueError, 'generation changed'):
                runtime.preflight(formal, 'reviewed')

    def test_protected_account_timeout_overrides_artifacts_and_reservation_not_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); account = root / 'account'; source = root / 'source.json'
            key = 'perturbdiff_fold0_seed17'
            source.write_text(json.dumps(dict(schema='perturbdiff-v2-migration-budget-v1',
                boundary_verified=True, ledger={key: dict(consumed_seconds=12475.943569269002,
                cap_seconds=14400, account_path=str(account), identity=dict(fold=0, seed=17),
                reconciliation_required=False, evidence=['reviewed debit'])})))
            initialize(source, hashlib.sha256(source.read_bytes()).hexdigest(), key, account)
            status, state = runtime.outcome(account)
            self.assertIsNone(status)
            state['attempts'] = [dict(prior_training_seconds=12475.943569269002, status='reserved')]
            state['consumed_seconds'] = 14400
            write(account / 'account.json', state)
            self.assertIsNone(runtime.outcome(account)[0])
            state['attempts'][0].update(status='timeout', returncode=0)
            write(account / 'account.json', state)
            self.assertEqual(runtime.outcome(account)[0], 'timeout')
            state['attempts'][0].update(status='exited', wall_seconds=1, returncode=5)
            state['consumed_seconds'] = 12476.943569269002
            write(account / 'account.json', state)
            self.assertEqual(runtime.outcome(account)[0], 'failed')


if __name__ == '__main__':
    unittest.main()
