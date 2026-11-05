"""Resource policy is tested with telemetry, never by exhausting the host."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'script' / 'schedule.py'


class SchedulerContract(unittest.TestCase):
    def policy(self, snapshot, profile=None, active=None):
        with tempfile.TemporaryDirectory() as td:
            source = Path(td) / 'request.json'
            source.write_text(json.dumps(dict(snapshot=snapshot, profile=profile,
                                             active=active or [])))
            p = subprocess.run([sys.executable, str(SCRIPT), 'admit', str(source)],
                               capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            return json.loads(p.stdout)

    def test_unknown_profile_starts_alone_and_respects_external_work(self):
        snapshot = dict(available_gib=30, gpu_free_gib=12, external_heavy=0,
                        external_peak_gib=0, observable=True)
        self.assertTrue(self.policy(snapshot)['admit'])
        snapshot['external_heavy'] = 1
        snapshot['external_peak_gib'] = 3
        self.assertFalse(self.policy(snapshot)['admit'])

    def test_known_profiles_need_headroom_and_host_visibility(self):
        s = dict(available_gib=24, gpu_free_gib=10, external_heavy=1,
                 external_peak_gib=8, observable=True)
        self.assertTrue(self.policy(s, dict(host_peak_gib=8, gpu_peak_gib=4))['admit'])
        self.assertFalse(self.policy(s, dict(host_peak_gib=12, gpu_peak_gib=4))['admit'])
        s['observable'] = False
        self.assertFalse(self.policy(s, dict(host_peak_gib=1))['admit'])

    def test_cumulative_budget_and_independent_retry_limits(self):
        with tempfile.TemporaryDirectory() as td:
            state = Path(td) / 'state.json'
            def transition(reason, seconds):
                p = subprocess.run([sys.executable, str(SCRIPT), 'transition',
                                    str(state), '--reason', reason,
                                    '--seconds', str(seconds)], capture_output=True, text=True)
                self.assertEqual(p.returncode, 0, p.stderr)
                return json.loads(p.stdout)
            first = transition('environment', 20)
            self.assertEqual(first['status'], 'pending')
            self.assertEqual(first['environment_retries'], 1)
            second = transition('oom', 30)
            self.assertEqual(second['batch_scale'], 0.5)
            self.assertEqual(second['training_seconds'], 50)
            self.assertEqual(transition('oom', 10)['status'], 'failed')
            # A terminal state cannot be reopened by another launch.
            self.assertEqual(transition('environment', 0)['status'], 'failed')
            state.unlink()
            self.assertEqual(transition('oom', 14400)['status'], 'timeout')

    def test_pressure_stops_only_newest_owned_task(self):
        with tempfile.TemporaryDirectory() as td:
            source = Path(td) / 'pressure.json'
            source.write_text(json.dumps(dict(
                history=[dict(available_gib=11, gpu_free_gib=2),
                         dict(available_gib=10, gpu_free_gib=1)],
                tasks=[dict(id='older', started=1, owned=True),
                       dict(id='unrelated', started=3, owned=False),
                       dict(id='newer', started=2, owned=True)])))
            p = subprocess.run([sys.executable, str(SCRIPT), 'pressure', str(source)],
                               capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertEqual(json.loads(p.stdout)['stop'], 'newer')


if __name__ == '__main__':
    unittest.main()
