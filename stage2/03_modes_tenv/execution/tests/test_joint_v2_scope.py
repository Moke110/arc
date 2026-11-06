"""Opt-in host scope contract: tiny Python/sleep workers, no training or data."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import unittest
from unittest.mock import patch
import uuid

import test_joint_v2_budget as budget_cli

CLI = budget_cli.CLI


@unittest.skipUnless(os.environ.get('ARC_TEST_HOST_SCOPES') == '1', 'requires host user systemd')
class OwnedScope(unittest.TestCase):
    setUp = budget_cli.BudgetCLI.setUp
    tearDown = budget_cli.BudgetCLI.tearDown
    cli = budget_cli.BudgetCLI.cli
    initialize = budget_cli.BudgetCLI.initialize

    def command(self, stay=False):
        self.scope = 'arc-tenv-v2-test-' + uuid.uuid4().hex + '.scope'
        self.addCleanup(subprocess.run, ['systemctl', '--user', 'kill', '--kill-whom=all',
            '--signal=SIGKILL', self.scope], capture_output=True, timeout=5)
        code = ('import subprocess,sys,pathlib,time; '
            'p=subprocess.Popen([sys.executable,"-c","import time;time.sleep(30)"],start_new_session=True); '
            'pathlib.Path(sys.argv[1]).write_text(str(p.pid)); ' + ('time.sleep(30)' if stay else ''))
        return ['systemd-run', '--user', '--scope', '--unit=' + self.scope,
            '-p', 'MemoryMax=64M', '-p', 'MemorySwapMax=0', '--', sys.executable,
            '-c', code, str(self.root / 'child.pid'), '{request}']

    def assert_empty(self):
        report = subprocess.run(['systemctl', '--user', 'show', self.scope,
            '--property=ActiveState,LoadState'], capture_output=True, text=True, check=False)
        values = dict(line.split('=', 1) for line in report.stdout.splitlines() if '=' in line)
        self.assertTrue(values.get('ActiveState') in ('inactive', 'failed') or
                        values.get('LoadState') == 'not-found', report.stdout)
        pid = int((self.root / 'child.pid').read_text())
        stat = Path(f'/proc/{pid}/stat')
        # Reparented zombies cannot execute; systemd need not reap them itself.
        self.assertTrue(not stat.exists() or stat.read_text().rsplit(')', 1)[1].split()[0] == 'Z')

    def wait_terminal(self):
        end = time.monotonic() + 9
        while time.monotonic() < end:
            state = json.loads((self.account / 'account.json').read_text())
            if state['attempts'] and state['attempts'][-1]['status'] in ('exited', 'timeout'):
                return state
            time.sleep(.02)
        self.fail((self.account / 'watchdog.log').read_text())

    def test_launcher_exit_cleans_daemon_before_refund(self):
        self.initialize(14398)
        self.cli('launch', '--account', self.account, '--request', self.request, '--', *self.command())
        state = self.wait_terminal()
        self.assert_empty()
        self.assertEqual(state['attempts'][-1]['status'], 'exited')
        self.assertLess(state['consumed_seconds'], 14400)

    def test_deadline_cleans_daemon_without_refund(self):
        self.initialize(14399)
        self.cli('launch', '--account', self.account, '--request', self.request, '--', *self.command(stay=True))
        state = self.wait_terminal()
        self.assert_empty()
        self.assertEqual(state['attempts'][-1]['status'], 'timeout')
        self.assertEqual(state['consumed_seconds'], 14400)

    def test_exception_after_launcher_exit_cleans_daemon_without_refund(self):
        self.initialize(14398)
        spec = importlib.util.spec_from_file_location('scope_budget', CLI)
        budget = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(budget)
        command = self.command()
        original_sleep = time.sleep
        injected = False
        def fail_after_launcher_exit(_):
            nonlocal injected
            if injected:
                return original_sleep(_)
            end = time.monotonic() + 5
            while time.monotonic() < end:
                state = json.loads((self.account / 'account.json').read_text())
                stat = Path('/proc') / str(state['attempts'][-1]['process']['pid']) / 'stat'
                if (self.root / 'child.pid').exists() and (not stat.exists() or
                        stat.read_text().rsplit(')', 1)[1].split()[0] == 'Z'):
                    injected = True
                    raise OSError('injected watchdog polling failure after launcher exit')
                original_sleep(.01)
            raise AssertionError('test launcher did not exit')
        with patch.object(budget.time, 'sleep', side_effect=fail_after_launcher_exit):
            with self.assertRaisesRegex(OSError, 'injected watchdog'):
                budget.supervise(self.account, self.request, command)
        self.assert_empty()
        state = json.loads((self.account / 'account.json').read_text())
        self.assertEqual(state['consumed_seconds'], 14400)
        self.assertNotIn(state['attempts'][-1]['status'], ('exited', 'timeout'))
        with self.assertRaisesRegex(ValueError, 'outstanding reservation'):
            budget.supervise(self.account, self.request, command)
        self.assertEqual(len(json.loads((self.account / 'account.json').read_text())['attempts']), 1)

    def test_unverifiable_cleanup_retains_reservation_and_blocks_retry(self):
        self.initialize(14398)
        spec = importlib.util.spec_from_file_location('scope_budget', CLI)
        budget = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(budget)
        command = self.command()
        original_run = subprocess.run
        def unavailable_scope_report(argv, **kwargs):
            if argv[:3] == ['systemctl', '--user', 'show']:
                raise subprocess.TimeoutExpired(argv, 5)
            return original_run(argv, **kwargs)
        with patch.object(budget.subprocess, 'run', side_effect=unavailable_scope_report):
            with self.assertRaises(subprocess.TimeoutExpired):
                budget.supervise(self.account, self.request, command)
        self.assert_empty()
        state = json.loads((self.account / 'account.json').read_text())
        self.assertEqual(state['consumed_seconds'], 14400)
        self.assertNotIn(state['attempts'][-1]['status'], ('exited', 'timeout'))
        with self.assertRaisesRegex(ValueError, 'outstanding reservation'):
            budget.supervise(self.account, self.request, command)


class CMongeOwnedScope(OwnedScope):
    model = 'cmonge'
