"""Resource policy tests use synthetic telemetry and tiny owned children only."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import resources


def snapshot(**changes):
    value = dict(memory_available_bytes=20 * resources.GIB, host_processes_observable=True,
                 gpu_observable=True, peers=[], processes=[])
    value.update(changes)
    return value


class ResourceTest(unittest.TestCase):
    def test_external_training_entrypoints_and_inline_code_redaction(self):
        commands = [
            ('python\0/other/project/train.py\0--token\0SECRET', True),
            ('python3.11\0-m\0package.trainer\0--token\0SECRET', True),
            ('python\0-X\0dev\0/other/train.py', True),
            ('python\0/other/worker.py', True),
            ('python\0/other/cmonge_v2_gpu_adapter.py', True),
            ('python\0/other/workflow.py\0fit', True),
            ('python\0/other/workflow.py\0train', True),
            ('python\0/other/workflow.py\0benchmark', True),
            ('python\0/other/pretrain.py', True),
            ('torchrun\0--nnodes\0' + '1\0--token\0SECRET', True),
            ('accelerate\0launch\0/other/job.py', True),
            ('deepspeed\0/other/job.py', True),
            ('python\0/other/monitor.py\0--label\0train\0--token\0SECRET', False),
            ('python\0-c\0print("SECRET cpsp/train.py")', False),
            ('bash\0-c\0python /other/train.py SECRET', False),
            ('python\0/other/monitor.py\0--input\0/other/train.py', False),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'meminfo').write_text('MemAvailable: 20971520 kB\n')
            (root/'1').mkdir()
            (root/'1/comm').write_text('systemd')
            (root/'1/cmdline').write_bytes(b'init')
            (root/'1/status').write_text('PPid:\t0\n')
            process = root/'12'
            process.mkdir()
            (process/'status').write_text('PPid:\t1\nVmRSS:\t100 kB\n')
            result = subprocess.CompletedProcess([], 0, '0, 100, 16000, 4\n', '')
            for command, expected in commands:
                with self.subTest(command=command), patch.object(resources.subprocess, 'run', return_value=result):
                    (process/'cmdline').write_bytes(command.encode())
                    value = resources.collect_snapshot(proc_root=root)
                    self.assertEqual(bool(value['peers']), expected)
                    self.assertNotIn('SECRET', json.dumps(value))
                    if expected:
                        with self.assertRaisesRegex(resources.ResourceError, 'external_long_training'):
                            resources.preflight(value)
                    else:
                        resources.preflight(value)

    def test_heavy_python_rss_threshold_is_independent_of_repository_location(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'meminfo').write_text('MemAvailable: 20971520 kB\n')
            for pid, parent in [(1, 0), (12, 1)]:
                process = root/str(pid)
                process.mkdir()
                (process/'comm').write_text('systemd' if pid == 1 else 'python')
                (process/'cmdline').write_bytes(b'init')
                (process/'status').write_text(f'PPid:\t{parent}\n')
            process = root/'12'
            result = subprocess.CompletedProcess([], 0, '0, 100, 16000, 4\n', '')
            cases = [
                ('python\0/home/arc/stage2/joint_v2_migration.py', 40 * 1024, False),
                ('python\0/home/arc/cpsp/joint_v2_apply.py', 40 * 1024, False),
                ('python\0/other/job.py', 1024 * 1024 - 1, False),
                ('python\0/other/job.py', 1024 * 1024, True),
                ('python\0-c\0print("SECRET")', 1024 * 1024 + 1, True),
                ('python\0-m\0other.job', 1024 * 1024, True),
            ]
            for command, rss_kib, expected in cases:
                with self.subTest(command=command, rss_kib=rss_kib), patch.object(
                        resources.subprocess, 'run', return_value=result):
                    (process/'cmdline').write_bytes(command.encode())
                    (process/'status').write_text(f'PPid:\t1\nVmRSS:\t{rss_kib} kB\n')
                    value = resources.collect_snapshot(proc_root=root)
                    self.assertEqual(bool(value['peers']), expected)
                    self.assertNotIn('SECRET', json.dumps(value))
                    if expected:
                        self.assertFalse(value['peers'][0]['long_training'])
                        with self.assertRaisesRegex(resources.ResourceError, 'serial_execution'):
                            resources.preflight(value)
                        resources.preflight(value, formal=False)
                        self.assertEqual(resources.collect_snapshot(
                            exclude_pids=(12,), proc_root=root)['peers'], [])
                    else:
                        resources.preflight(value)

    def test_policy_and_unknown_visibility(self):
        resources.preflight(snapshot())
        for changes, message in [
            ({'memory_available_bytes': 11 * resources.GIB}, 'available_ram'),
            ({'host_processes_observable': False}, 'inventory_unobservable'),
            ({'gpu_observable': False}, 'gpu_telemetry'),
            ({'peers': [dict(long_training=True)]}, 'external_long_training'),
            ({'peers': [dict(long_training=False)]}, 'serial_execution')]:
            with self.assertRaisesRegex(resources.ResourceError, message):
                resources.preflight(snapshot(**changes))
        resources.preflight(snapshot(peers=[dict(long_training=False)]), formal=False)
        with self.assertRaisesRegex(resources.ResourceError, 'two_heavy'):
            resources.preflight(snapshot(peers=[dict(long_training=False)] * 2), formal=False)

    def test_proc_fixture_redaction_gpu_and_owned_exclusion(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'meminfo').write_text('MemAvailable: 20971520 kB\n')
            for pid, parent, command in [(1, 0, 'init'), (12, 1, 'python\0/home/arc/stage2/fit.py\0--token\0SECRET'),
                                         (13, 12, 'python\0/home/arc/stage2/worker.py')]:
                path = root/str(pid)
                path.mkdir()
                (path/'comm').write_text('systemd' if pid == 1 else 'python')
                (path/'cmdline').write_bytes(command.encode())
                (path/'status').write_text(f'PPid:\t{parent}\nVmRSS:\t100 kB\n')
            result = subprocess.CompletedProcess([], 0, '0, 100, 16000, 4\n', '')
            with patch.object(resources.subprocess, 'run', return_value=result) as query:
                value = resources.collect_snapshot(proc_root=root)
                self.assertTrue(value['host_processes_observable'])
                self.assertEqual(len(value['peers']), 2)
                self.assertTrue(value['peers'][0]['long_training'])
                self.assertNotIn('SECRET', json.dumps(value))
                self.assertEqual(query.call_args.kwargs['timeout'], 3)
                owned = resources.collect_snapshot(exclude_pids=(12,), proc_root=root)
                self.assertEqual(owned['peers'], [])
                (root/'1/comm').write_text('python')
                self.assertFalse(resources.collect_snapshot(proc_root=root)['host_processes_observable'])

    def test_guard_stops_only_owned_child_and_records_peak(self):
        calls = []
        def telemetry(exclude_pids):
            calls.append(exclude_pids)
            if not exclude_pids:
                return snapshot()
            pid = exclude_pids[0]
            return snapshot(memory_available_bytes=resources.GIB,
                            processes=[dict(pid=pid, ppid=1, rss_bytes=1234),
                                       dict(pid=pid + 100000, ppid=pid, rss_bytes=456)])
        with tempfile.TemporaryDirectory() as tmp:
            result = resources.run_guarded([sys.executable, '-c', 'import time; time.sleep(30)'],
                                           tmp, snapshot_fn=telemetry, interval=.01)
            self.assertEqual(result['stop_reason'], 'available_ram_below_12_gib')
            self.assertEqual(result['peak_tree_rss_bytes'], 1690)
            self.assertNotEqual(result['exit_code'], 0)
            self.assertEqual(len((Path(tmp)/'resources.jsonl').read_text().splitlines()), 2)

    def test_rejected_preflight_does_not_launch(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(resources.subprocess, 'Popen') as launch:
            with self.assertRaises(resources.ResourceError):
                resources.run_guarded(['unused'], tmp, snapshot_fn=lambda **kw: snapshot(gpu_observable=False))
            launch.assert_not_called()


if __name__ == '__main__':
    unittest.main()
