"""Exercise the Step 10 resource audit for an evidenced idle interruption."""

import json
from hashlib import sha256
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dose_interruption import Scope, active_samples, scope_windows, verify_attempt


class ScopeWindowsTest(unittest.TestCase):
    def test_bounds_are_read_from_two_pinned_same_boot_scope_journals(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            journals = []
            for index, (start, end) in enumerate(((10, 40), (160, 200))):
                unit = f'run-example{index}.scope'
                path = root / (unit + '.jsonl')
                records = (
                    ('39f53479d3a045ac8e11786248231fbf', start,
                     f'Started {unit} - python dose_matrix.py --output /fit/matrix --logs /fit/matrix-audit.'),
                    ('ae8f7b866b0347b9af31fe1c80b127c0', end, f'{unit}: Consumed 12s CPU time.'),
                )
                path.write_text(''.join(json.dumps({
                    'USER_UNIT': unit, '_BOOT_ID': '1' * 32, 'MESSAGE_ID': event,
                    'USER_INVOCATION_ID': str(index) * 32, 'CPU_USAGE_NSEC': '1',
                    '__MONOTONIC_TIMESTAMP': str(second * 10**6),
                    '__REALTIME_TIMESTAMP': str((1000 + second) * 10**6), 'MESSAGE': message,
                }) + '\n' for event, second, message in records))
                journals.append({'unit': unit, 'path': path.name, 'sha256': sha256(path.read_bytes()).hexdigest(),
                                 'monitor_session': str(index) * 32})
            windows = scope_windows(root, journals)
            self.assertEqual([(window.start, window.stop) for window in windows], [(10, 40), (160, 200)])
            first = root / journals[0]['path']
            rows = [json.loads(line) for line in first.read_text().splitlines()]
            rows[-1]['USER_INVOCATION_ID'] = 'f' * 32
            first.write_text(''.join(json.dumps(row) + '\n' for row in rows))
            journals[0]['sha256'] = sha256(first.read_bytes()).hexdigest()
            with self.assertRaisesRegex(ValueError, 'invocation'):
                scope_windows(root, journals)
            journals[0]['sha256'] = '0' * 64
            with self.assertRaisesRegex(ValueError, 'journal digest'):
                scope_windows(root, journals)


class ActiveCoverageTest(unittest.TestCase):
    def test_idle_gap_is_audit_evidence_and_active_gap_remains_a_failure(self):
        scopes = [Scope('first', 10, 39, 1010, 1039), Scope('last', 160, 191, 1160, 1191)]
        supervisor = [
            {'monotonic_seconds': tick, 'wall_time': f'1970-01-01T00:{(1000+tick)//60:02d}:{(1000+tick)%60:02d}+00:00',
             'gpu_total_used_mib_and_utilization_percent': '16000,1000,5',
             'windows_host_total_and_available_kib': '25000000,18000000',
             'wsl': {'MemAvailable_kib': 14 * 1024**2}}
            for tick in (8, 18, 28, 38, 155, 165, 175, 185, 195)
        ]
        matrix = [
            {'monitor_session': session, 'monotonic_seconds': tick, 'scope_memory_bytes': 1024}
            for session, ticks in (('first', (11, 21, 31)), ('last', (162, 172, 182))) for tick in ticks
        ]
        selected, quality = active_samples(supervisor, matrix, scopes, 1189)
        self.assertEqual(selected, supervisor)
        self.assertEqual([scope['coordinator_count'] for scope in quality], [4, 5])
        with self.assertRaisesRegex(ValueError, 'cadence'):
            active_samples(supervisor[:2] + supervisor[3:], matrix, scopes, 1189)

    def test_phase_measurement_cannot_cross_into_the_idle_window(self):
        scopes = [Scope('first', 10, 39, 1010, 1039), Scope('last', 160, 191, 1160, 1191)]
        timing = ['Elapsed (wall clock) time (h:mm:ss or m:ss): 0:02.00']
        verify_attempt(Path('.'), 'train', Path('train-1036000000000.resources'), timing, scopes)
        with self.assertRaisesRegex(ValueError, 'inside one evidenced scope'):
            verify_attempt(Path('.'), 'train', Path('train-1038000000000.resources'), timing, scopes)


if __name__ == '__main__':
    unittest.main()
