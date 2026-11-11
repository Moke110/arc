"""Check that callers can immediately consume both capped worker logs."""

import os
from pathlib import Path
import subprocess
import time
from tempfile import TemporaryDirectory
import unittest


LAUNCHER = Path(__file__).resolve().parents[1] / "script" / "run_capped.sh"


class CappedLauncherTest(unittest.TestCase):
    def test_cancel_pending_blocks_new_worker_without_creating_run_logs(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            work = root / 'formal'
            work.mkdir()
            (work / 'cancel-pending.json').write_text('{"status":"finish_current_then_stop"}')
            logs = work / 'resources' / 'pending-model'
            called = root / 'worker-called'
            result = subprocess.run([LAUNCHER, logs, '/usr/bin/touch', called],
                                    capture_output=True, text=True, timeout=5)
            self.assertEqual(result.returncode, 75)
            self.assertFalse(called.exists())
            self.assertFalse(logs.exists())
            self.assertIn('cancelled', result.stderr)

    def test_live_script_rewrite_does_not_launch_worker_twice(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            bin_path = root / "bin"
            bin_path.mkdir()
            programs = {
                "systemd-run": 'while [[ "$1" != "--" ]]; do shift; done; shift; exec "$@"',
                "systemctl": "exit 1",
                "nvidia-smi": "echo '0, 1, 0'",
            }
            for name, body in programs.items():
                tool = bin_path / name
                tool.write_text("#!/bin/bash\n" + body + "\n")
                tool.chmod(0o755)
            launcher = root / "launcher.sh"
            original = LAUNCHER.read_text()
            launcher.write_text(original)
            launcher.chmod(0o755)
            worker = root / "worker.sh"
            worker.write_text('#!/bin/bash\necho run >> "$1"\ntouch "$2"\n'
                              'while [[ ! -f "$3" ]]; do sleep .01; done\n')
            worker.chmod(0o755)
            process = subprocess.Popen(
                [launcher, root / "logs", worker, root / "calls", root / "started", root / "release"],
                env={**os.environ, "PATH": f"{bin_path}:{os.environ['PATH']}",
                     "XDG_RUNTIME_DIR": str(root)},
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            )
            try:
                deadline = time.monotonic() + 5
                while not (root / "started").exists() and time.monotonic() < deadline:
                    time.sleep(.01)
                self.assertTrue((root / "started").exists())
                addition = ('gpu_monitor=$(command -v nvidia-smi || true)\n'
                            'if [[ -z "$gpu_monitor" && -x /usr/lib/wsl/lib/nvidia-smi ]]; then\n'
                            '    gpu_monitor=/usr/lib/wsl/lib/nvidia-smi\nfi\n')
                launcher.write_text(original.replace('scope="arc-step09-$$.scope"',
                                                     addition + 'scope="arc-step09-$$.scope"'))
                (root / "release").touch()
                self.assertEqual(process.wait(timeout=5), 0)
                self.assertEqual((root / "calls").read_text().splitlines(), ["run"])
            finally:
                if process.poll() is None:
                    process.kill()
                    process.wait()

    def test_waits_for_both_log_writers_and_preserves_worker_exit(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            bin_path = root / "bin"
            bin_path.mkdir()
            programs = {
                "systemd-run": 'while [[ "$1" != "--" ]]; do shift; done; shift; exec "$@"',
                "systemctl": "exit 1",
                "nvidia-smi": "echo '0, 1, 0'",
                "tee": 'case "$1" in */stdout.log|*/stderr.log) sleep .4;; esac; exec /usr/bin/tee "$@"',
            }
            for name, body in programs.items():
                tool = bin_path / name
                tool.write_text("#!/usr/bin/env bash\n" + body + "\n")
                tool.chmod(0o755)
            logs = root / "logs"
            result = subprocess.run(
                [LAUNCHER, logs, "/bin/bash", "-c", "echo final-out; echo final-err >&2; exit 7"],
                env={**os.environ, "PATH": f"{bin_path}:{os.environ['PATH']}", "XDG_RUNTIME_DIR": str(root)},
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                timeout=10,
            )
            self.assertEqual(result.returncode, 7)
            self.assertEqual((logs / "stdout.log").read_text(), "final-out\n")
            self.assertEqual((logs / "stderr.log").read_text(), "final-err\n")
            start, finish = map(float, (logs / "wall-clock.tsv").read_text().split())
            self.assertLessEqual(start, finish)
