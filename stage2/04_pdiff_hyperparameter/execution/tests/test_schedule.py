"""Exercise scheduler admission through the CLI and emitted records."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'script' / 'schedule.py'


class SchedulerCLI(unittest.TestCase):
    def setUp(self):
        root = tempfile.TemporaryDirectory()
        self.addCleanup(root.cleanup)
        self.root = Path(root.name)
        self.queue = self.root / 'queue.json'
        self.fixture = self.root / 'fixture.json'
        self.record = self.root / 'record'
        self.job = dict(id='profile', command=[sys.executable, '-c', 'pass'], kind='gpu',
                        profile=dict(host_peak_gib=3, gpu_peak_gib=2), disk_gib=2)
        self.queue.write_text(json.dumps(dict(record_dir='record', jobs=[self.job])))

    def invoke(self, *args, env=None):
        return subprocess.run([sys.executable, str(SCRIPT), args[0], str(self.queue), *args[1:]],
                              capture_output=True, text=True, timeout=10, env=env)

    def sample(self, **changes):
        sample = dict(observable=True, age_seconds=0, available_gib=32, gpu_free_gib=12,
                      gpu_used_gib=4, gpu_utilization_percent=25, disk_free_gib=30,
                      external_heavy=0, external_peak_gib=0, external_pids=[])
        self.fixture.write_text(json.dumps(sample | changes))

    def wait_status(self, predicate):
        deadline = time.monotonic() + 18
        while time.monotonic() < deadline:
            result = self.invoke('status')
            if result.returncode == 0 and predicate(evidence := json.loads(result.stdout)):
                return evidence
            time.sleep(0.05)
        self.fail('scheduler status did not reach expected state: ' + result.stdout + result.stderr)

    def fake_scopes(self):
        # Replace systemd launch in the simulated workstation only. The public
        # scheduler still starts and supervises a real, bounded subprocess.
        for name, body in {
            'systemd-run': 'while [ "$1" != -- ]; do shift; done; shift; exec "$@"',
            'systemctl': 'printf "inactive\\n"',
        }.items():
            executable = self.root / name
            executable.write_text('#!/bin/sh\n' + body + '\n')
            executable.chmod(0o755)
        return dict(os.environ, PATH=str(self.root) + ':' + os.environ['PATH'])

    def start(self):
        self.process = subprocess.Popen(
            [sys.executable, str(SCRIPT), 'run', str(self.queue),
             '--fixture', str(self.fixture)], env=self.fake_scopes(),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.addCleanup(self.finish, self.process)
        return self.process

    def finish(self, process):
        if process.poll() is None:
            process.terminate()
        process.communicate(timeout=20)

    def test_pause_drains_and_resume_runs_only_the_remaining_queue(self):
        self.sample()
        self.job['command'] = [sys.executable, '-c', 'import time; time.sleep(5)']
        self.queue.write_text(json.dumps(dict(record_dir='record', jobs=[
            self.job, self.job | dict(id='second', kind='training', profile=None)])))
        process = self.start()
        active = self.wait_status(lambda status: status['active'] == ['profile'])
        self.assertEqual(active['queue'], ['second'])
        self.assertEqual(self.invoke('pause').returncode, 0)
        drained = self.wait_status(lambda status: status['paused'] and not status['active'])
        self.assertEqual(drained['jobs']['profile']['status'], 'completed')
        self.assertEqual(drained['queue'], ['second'])
        attempt, = drained['jobs']['profile']['attempts']
        self.assertGreater(attempt['peak_host_rss_gib'], 0)
        self.assertEqual(attempt['memory_swap_max_gib'], 8)
        self.assertEqual(process.wait(timeout=10), 0)
        self.assertEqual(self.invoke('resume').returncode, 0)
        self.start()
        completed = self.wait_status(lambda status: status['jobs']['second']['status'] == 'completed')
        self.assertEqual(len(completed['jobs']['profile']['attempts']), 1)

    def test_persistent_pressure_only_stops_owned_worker(self):
        self.sample()
        unrelated = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(45)'])
        self.addCleanup(unrelated.wait, 5)
        self.addCleanup(unrelated.terminate)
        self.job['command'] = [sys.executable, '-c', 'import time; time.sleep(40)']
        self.queue.write_text(json.dumps(dict(record_dir='record', jobs=[self.job])))
        worker = self.start()
        self.wait_status(lambda status: status['active'] == ['profile'])
        self.sample(available_gib=10, external_heavy=1, external_pids=[unrelated.pid])
        stopped = self.wait_status(lambda status: status['jobs']['profile']['status'] == 'interrupted')
        attempt, = stopped['jobs']['profile']['attempts']
        self.assertEqual(attempt['reason'], 'less than 12 GiB available RAM')
        self.assertIsNone(unrelated.poll())
        self.assertEqual(worker.wait(timeout=10), 0)
        retry = self.invoke('retry', 'profile', env=self.fake_scopes())
        self.assertEqual(retry.returncode, 0, retry.stderr)
        self.assertEqual(json.loads(retry.stdout)['queue'], ['profile'])

    def test_measured_compatible_workers_share_two_slots(self):
        self.sample()
        self.job['command'] = [sys.executable, '-c', 'import time; time.sleep(6)']
        second = self.job | dict(id='cpu', kind='cpu', profile=dict(host_peak_gib=2, gpu_peak_gib=0))
        third = second | dict(id='third')
        self.queue.write_text(json.dumps(dict(record_dir='record', jobs=[self.job, second, third])))
        self.start()
        active = self.wait_status(lambda status: len(status['active']) == 2)
        self.assertEqual(active['active'], ['profile', 'cpu'])
        self.assertEqual(active['queue'], ['third'])
        self.assertEqual(self.invoke('pause').returncode, 0)

    def test_exclusive_worker_yields_to_external_pressure(self):
        self.sample()
        self.job['kind'] = 'training'
        self.job['command'] = [sys.executable, '-c', 'import time; time.sleep(40)']
        self.queue.write_text(json.dumps(dict(record_dir='record', jobs=[self.job])))
        self.start()
        self.wait_status(lambda status: status['active'] == ['profile'])
        self.sample(external_heavy=1, external_pids=[2468])
        stopped = self.wait_status(lambda status: status['jobs']['profile']['status'] == 'interrupted')
        self.assertEqual(stopped['jobs']['profile']['attempts'][0]['reason'],
                         'external heavy worker conflicts with exclusive training')

    def test_scheduler_interrupt_retains_attempt(self):
        self.sample()
        self.job['command'] = [sys.executable, '-c', 'import time; time.sleep(40)']
        self.queue.write_text(json.dumps(dict(record_dir='record', jobs=[self.job])))
        worker = self.start()
        self.wait_status(lambda status: status['active'] == ['profile'])
        worker.terminate()
        stopped = self.wait_status(lambda status: status['jobs']['profile']['status'] == 'interrupted')
        attempt, = stopped['jobs']['profile']['attempts']
        self.assertEqual(attempt['reason'], 'scheduler interrupted')
        self.assertEqual(worker.wait(timeout=10), 130)

    def test_preflight_uses_elapsed_age_when_wall_clock_jumps(self):
        self.sample(age_seconds=0.5, timestamp=time.time() - 60)
        result = self.invoke('preflight', '--fixture', str(self.fixture))
        self.assertEqual(result.returncode, 0, result.stderr)
        evidence = json.loads(result.stdout)['sample']
        self.assertGreater(evidence['observed_at'] - evidence['timestamp'], 59)
        self.assertLess(evidence['observed_monotonic'] - evidence['sample_started_monotonic'], 5)

    def test_preflight_gates_on_live_resources(self):
        for changes, reason in [
            ({}, 'resource budgets satisfied'),
            ({'age_seconds': 6}, 'stale telemetry'),
            ({'observable': False}, 'missing telemetry'),
            ({'available_gib': 11}, 'less than 12 GiB available RAM'),
            ({'gpu_utilization_percent': 95}, 'GPU utilization at or above 90 percent'),
            ({'external_heavy': 1, 'external_pids': [123]}, 'unprofiled external heavy job'),
            ({'disk_free_gib': 2}, 'insufficient free disk'),
        ]:
            with self.subTest(reason=reason):
                self.sample(**changes)
                result = self.invoke('preflight', '--fixture', str(self.fixture))
                evidence = json.loads(result.stdout)
                self.assertEqual(evidence['decisions']['profile']['reason'], reason)
                self.assertEqual(result.returncode, 0 if not changes else 2)


if __name__ == '__main__':
    unittest.main()
