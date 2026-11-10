"""Check that callers can immediately consume both capped worker logs."""

import os
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
import unittest


LAUNCHER = Path(__file__).resolve().parents[1] / "script" / "run_capped.sh"


class CappedLauncherTest(unittest.TestCase):
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
