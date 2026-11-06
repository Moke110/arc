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
    def test_foreign_account_refuses_before_reuse_or_dispatch(self):
        for mismatch in ('ledger', 'key'):
            for previous_status in (None, 'completed', 'timeout'):
                with self.subTest(mismatch=mismatch, previous_status=previous_status), tempfile.TemporaryDirectory() as tmp:
                    root = Path(tmp); formal = root / 'formal'; formal.mkdir()
                    output = root / 'output'; (output / 'scheduling').mkdir(parents=True)
                    account = root / 'account'
                    job = dict(model='perturbdiff', fold=0, seed=17, command=['fake_adapter'])
                    job_key = runtime.key(job)
                    foreign_key = 'perturbdiff_fold1_seed17'
                    row = dict(consumed_seconds=0, cap_seconds=14400, account_path=str(account),
                               identity=job, reconciliation_required=False, evidence=['reviewed fixture'])
                    ledger = dict(schema='perturbdiff-v2-migration-budget-v1', boundary_verified=True,
                                  ledger={job_key: row, foreign_key: dict(row, identity=dict(job, fold=1))})
                    source = root / 'source.json'; source.write_text(json.dumps(ledger))
                    initialize(source, hashlib.sha256(source.read_bytes()).hexdigest(),
                               foreign_key if mismatch == 'key' else job_key, account)
                    if mismatch == 'ledger':
                        ledger['ledger'][job_key] = dict(row, consumed_seconds=14400)
                    (formal / 'budget-ledger.json').write_text(json.dumps(ledger))
                    queue = dict(input_root=str(root / 'input'), output_root=str(output))
                    directory = runtime.run_directory(output, job['model'], job['fold'], job['seed'])
                    directory.mkdir(parents=True)
                    if previous_status:
                        (directory / 'run.json').write_text(json.dumps(dict(identity=job, status=previous_status)))
                    with patch.object(runtime, 'validate_run') as validate, patch.object(runtime.subprocess, 'run') as launch, patch('schedule.process_snapshot', side_effect=AssertionError('admission reached')):
                        with self.assertRaisesRegex(ValueError, 'account .*binding'):
                            runtime.run_v2(formal, 'reviewed', queue, job, job)
                        validate.assert_not_called()
                        launch.assert_not_called()
                    self.assertFalse((directory / 'request.json').exists())

    def test_bound_exhausted_account_publishes_timeout_without_dispatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); formal = root / 'formal'; formal.mkdir()
            output = root / 'output'; (output / 'scheduling').mkdir(parents=True)
            account = root / 'account'
            job = dict(model='perturbdiff', fold=0, seed=17, command=['fake_adapter'])
            ledger = dict(schema='perturbdiff-v2-migration-budget-v1', boundary_verified=True,
                          ledger={runtime.key(job): dict(consumed_seconds=14400, cap_seconds=14400,
                          account_path=str(account), identity=job, reconciliation_required=False,
                          evidence=['reviewed exhausted fixture'])})
            (formal / 'budget-ledger.json').write_text(json.dumps(ledger))
            queue = dict(input_root=str(root / 'input'), output_root=str(output))
            with patch.object(runtime.subprocess, 'run') as launch:
                runtime.run_v2(formal, 'reviewed', queue, job, job)
                runtime.run_v2(formal, 'reviewed', queue, job, job)
                launch.assert_not_called()
            directory = runtime.run_directory(output, job['model'], job['fold'], job['seed'])
            state = json.loads((directory / 'run.json').read_text())
            self.assertEqual(state['status'], 'timeout')
            self.assertEqual(state['training_seconds'], 14400)

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
