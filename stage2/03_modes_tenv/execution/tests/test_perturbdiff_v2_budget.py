"""Artifact CLI tests; no native imports or datasets."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest

CLI = Path(__file__).resolve().parents[1] / 'script/perturbdiff_v2_budget.py'
KEY = 'perturbdiff_fold0_seed17'


class BudgetCLI(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.account = self.root / 'account'
        self.source = self.root / 'migration.json'
        self.request = self.root / 'request.json'
        self.request.write_text(json.dumps({'identity': {'fold': 0, 'seed': 17}}))

    def tearDown(self):
        self.temp.cleanup()

    def cli(self, *args, ok=True):
        result = subprocess.run([sys.executable, str(CLI), *map(str, args)],
                                text=True, capture_output=True)
        self.assertEqual(result.returncode == 0, ok, result.stderr)
        return result

    def initialize(self, consumed=100, cap=14400):
        self.source.write_text(json.dumps({'schema': 'perturbdiff-v2-migration-budget-v1',
            'boundary_verified': True, 'ledger': {KEY: {'consumed_seconds': consumed,
            'cap_seconds': cap, 'account_path': str(self.account.resolve()), 'identity': {'fold': 0, 'seed': 17}, 'reconciliation_required': False, 'evidence': ['reviewed source hashes']}}}))
        digest = hashlib.sha256(self.source.read_bytes()).hexdigest()
        return self.cli('initialize', '--source', self.source, '--sha256', digest,
            '--key', KEY, '--account', self.account, ok=consumed is not None and cap == 14400)

    def test_detached_inference_deadline_is_terminal_timeout(self):
        self.initialize(100)
        (self.root / 'progress.json').write_text(json.dumps({'phase': 'inference'}))
        self.request.write_text(json.dumps({'identity': {'fold': 0, 'seed': 17},
            'output_directory': str(self.root), 'inference_timeout_seconds': .05}))
        self.cli('launch', '--account', self.account, '--request', self.request,
            '--', sys.executable, '-c', 'import time;time.sleep(30)', '{request}')
        deadline = time.monotonic() + 5
        while time.monotonic() < deadline:
            state = json.loads((self.account / 'account.json').read_text())
            if state['attempts'] and state['attempts'][-1]['status'] == 'timeout':
                break
            time.sleep(.02)
        else:
            self.fail((self.account / 'watchdog.log').read_text())
        self.assertEqual(state['consumed_seconds'], 14400)
        self.assertEqual(state['attempts'][-1]['stop_reason'], 'inference_timeout')

    def test_missing_debit_and_cap_reset(self):
        self.initialize(None)
        self.initialize(100, 28800)
        self.assertFalse(self.account.exists())

    def test_detached_retry_charges_and_injects(self):
        self.initialize(12475.943569269002)
        for index in range(2):
            code = 'import json,sys,time; r=json.load(open(sys.argv[1])); assert r["prior_training_seconds"]>=12475.943569269002; time.sleep(.12)'
            launched = json.loads(self.cli('launch', '--account', self.account, '--request', self.request,
                '--', sys.executable, '-c', code, '{request}').stdout)
            self.assertFalse(launched['dispatch_confirmed'])
            deadline = time.monotonic() + 5
            while time.monotonic() < deadline:
                state = json.loads((self.account / 'account.json').read_text())
                if len(state['attempts']) == index + 1 and state['attempts'][-1]['status'] == 'exited':
                    break
                time.sleep(.02)
            else:
                self.fail((self.account / 'watchdog.log').read_text())
        self.assertGreater(state['consumed_seconds'], 12476.18)
        self.assertGreater(state['attempts'][1]['prior_training_seconds'], state['attempts'][0]['prior_training_seconds'])

    def test_ambiguous_reservation_and_tampered_source_refused(self):
        self.initialize()
        path = self.account / 'account.json'
        state = json.loads(path.read_text()); state['consumed_seconds'] = 14400
        path.write_text(json.dumps(state))
        self.cli('supervise', '--account', self.account, '--request', self.request,
                 '--', sys.executable, '-c', 'pass', '{request}', ok=False)
        self.source.write_text('{}')
        (self.account / 'migration.json').write_text('{}')
        self.cli('inspect', '--account', self.account, ok=False)

    def test_timeout_is_terminal_and_no_credit_reset(self):
        self.initialize(14399.8)
        self.cli('launch', '--account', self.account, '--request', self.request,
                 '--', sys.executable, '-c', 'import time; time.sleep(30)', '{request}')
        deadline = time.monotonic() + 5
        while time.monotonic() < deadline:
            state = json.loads((self.account / 'account.json').read_text())
            if state['attempts'] and state['attempts'][-1]['status'] == 'timeout':
                break
            time.sleep(.02)
        else:
            self.fail('watchdog failed to enforce cumulative cap')
        self.assertEqual(state['consumed_seconds'], 14400)
        # Corruption that tries to retain history but reset account credit fails.
        state['consumed_seconds'] = 14399.8
        (self.account / 'account.json').write_text(json.dumps(state))
        self.cli('inspect', '--account', self.account, ok=False)

    def test_pid_reuse_identity_rejected(self):
        # Exercise the same exact /proc identity comparison used for inspection.
        import importlib.util
        spec = importlib.util.spec_from_file_location('budget', CLI)
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        import os
        identity = module.process_identity(os.getpid())
        self.assertTrue(module.matches(identity))
        identity['start_ticks'] += 1
        self.assertFalse(module.matches(identity))


if __name__ == '__main__':
    unittest.main()
