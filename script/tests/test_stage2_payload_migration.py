import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "migrate_stage2_payloads.py"


class PayloadMigrationTest(unittest.TestCase):
    def test_preflight_conflict_then_move_link_and_restore(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "active"
            root.mkdir()
            subprocess.run(["git", "init", "-q", str(root)], check=True)
            (root / ".gitignore").write_text("*.npy\n")
            old = root / "stage2/04_tsg_dose/execution/checkpoints/inputs.npy"
            old.parent.mkdir(parents=True)
            old.write_bytes(b"original checkpoint")
            inode = old.stat().st_ino
            new = root / "stage2/06_tsg_dose/execution/checkpoints/inputs.npy"
            plan = root.parent / "plan.json"

            def run(action, *args, success=True):
                result = subprocess.run(
                    [sys.executable, str(SCRIPT), action, "--root", str(root),
                     "--plan", str(plan), *map(str, args)], text=True, capture_output=True,
                )
                self.assertEqual(result.returncode == 0, success, result.stderr)
                return result

            run("plan")
            new.parent.mkdir(parents=True)
            new.write_bytes(b"collision")
            self.assertIn("destination exists", run("apply", success=False).stderr)
            self.assertEqual(old.read_bytes(), b"original checkpoint")
            new.unlink()
            run("apply")
            self.assertFalse((root / "stage2/04_tsg_dose").exists())
            self.assertEqual(new.stat().st_ino, inode)
            snapshot = root.parent / "snapshot"
            subprocess.run(["git", "init", "-q", str(snapshot)], check=True)
            run("link-snapshot", "--target", snapshot)
            replay = snapshot / old.relative_to(root)
            self.assertTrue(os.path.samefile(replay, new))
            self.assertEqual(replay.resolve(), replay)
            run("reverse")
            self.assertEqual(old.stat().st_ino, inode)
            self.assertFalse(new.exists())


if __name__ == "__main__":
    unittest.main()
