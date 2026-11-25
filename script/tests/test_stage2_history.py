"""Archive installation on a small Git checkout with materialized payloads."""
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "stage2_history.py"
OLD_RELEASE = Path("stage2/03_data_train/03.1_tenv_hvg-deg-5k")
NEW_RELEASE = Path("stage2/03_data_train/03.1_tenv_hvg-5k")
RETIRED = Path("stage2/05_tsg_loss")


class ArchiveInstallTest(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(TemporaryDirectory()))
        self.git("init", "-q", "-b", "main")
        self.git("config", "user.name", "Fixture")
        self.git("config", "user.email", "fixture@example.invalid")
        self.write(".gitignore", b"*.h5\n/tmp/\n")
        self.write(RETIRED / "README.md", b"original evidence\n")
        self.write(OLD_RELEASE / "output/manifest.json", b"original manifest\n")
        self.git("add", ".")
        self.git("commit", "-qm", "Predecessor")
        self.old = self.git("rev-parse", "HEAD").strip()
        self.retired = self.write(RETIRED / "output/frozen.h5", b"sealed data")
        self.active = self.write(OLD_RELEASE / "output/active.h5", b"input data")
        self.state = self.root / "tmp/cutover"
        self.target = self.root / "tmp/target"
        self.git("worktree", "add", "-qb", "dev", str(self.target))
        self.git("rm", "-r", str(RETIRED), root=self.target)
        (self.target / OLD_RELEASE).rename(self.target / NEW_RELEASE)
        self.write(".gitignore", b"*.h5\n/tmp/\n/stage2/history/\n", root=self.target)
        self.git("add", ".", root=self.target)
        self.git("commit", "-qm", "Successor", root=self.target)

    def git(self, *arguments, root=None):
        return subprocess.check_output(("git", *arguments), cwd=root or self.root, text=True)

    def write(self, path, body, root=None):
        path = (root or self.root) / path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(body)
        return path

    def execute(self, action, *args, check=True):
        return subprocess.run(
            (sys.executable, str(SCRIPT), action, "--root", str(self.root),
             "--state", str(self.state), *args),
            text=True, capture_output=True, check=check,
        )

    def test_preserves_tracked_evidence_and_moves_ignored_inodes(self):
        retired_inode, active_inode = self.retired.stat().st_ino, self.active.stat().st_ino
        self.execute("prepare", "--target", "dev")
        self.assertEqual(self.retired.stat().st_ino, retired_inode)
        self.assertEqual(self.active.stat().st_ino, active_inode)
        self.git("merge", "--no-overwrite-ignore", "--no-ff", "--no-edit", "dev")
        # Model an interruption after a payload was linked but before its old name was removed.
        (self.state / "archive" / RETIRED / "output/frozen.h5").hardlink_to(self.retired)
        self.execute("install")
        self.execute("install")
        history = self.root / "stage2/history"
        self.assertEqual((history / RETIRED / "README.md").read_bytes(), b"original evidence\n")
        self.assertEqual((history / RETIRED / "output/frozen.h5").stat().st_ino, retired_inode)
        self.assertEqual((self.root / NEW_RELEASE / "output/active.h5").stat().st_ino, active_inode)
        self.assertFalse((self.root / RETIRED).exists())
        self.assertFalse((self.root / OLD_RELEASE).exists())
        self.execute("verify")
        self.execute("undo")
        self.assertFalse(history.exists())
        self.assertEqual(self.retired.stat().st_ino, retired_inode)
        self.assertEqual(self.active.stat().st_ino, active_inode)
        self.assertFalse((self.root / NEW_RELEASE / "output/active.h5").exists())
        self.execute("install")
        self.execute("verify")

    def test_prepare_refuses_effectively_unignored_archive(self):
        self.write(".gitignore", b"*.h5\n/tmp/\n/stage2/history/\n!/stage2/history/\n", root=self.target)
        self.git("commit", "-qam", "Negate archive ignore", root=self.target)
        failed = self.execute("prepare", "--target", "dev", check=False)
        self.assertEqual(failed.returncode, 1)
        self.assertIn("target Git ignore does not cover", failed.stderr)
        self.assertFalse(self.state.exists())
        self.assertTrue(self.retired.exists())

    def test_prepare_refuses_ignored_payload_collision_before_merge(self):
        collision = self.write(NEW_RELEASE / "output/active.h5", b"other data")
        failed = self.execute("prepare", "--target", "dev", check=False)
        self.assertEqual(failed.returncode, 1)
        self.assertIn("destination is occupied", failed.stderr)
        self.assertEqual(self.git("rev-parse", "HEAD").strip(), self.old)
        self.assertEqual(collision.read_bytes(), b"other data")
        self.assertFalse(self.state.exists())
        self.assertTrue(self.active.exists())

    def test_prepare_refuses_ignored_tracked_destination_before_merge(self):
        self.write(".gitignore", b"*.h5\n/tmp/\n/stage2/03_data_train/03.1_tenv_hvg-5k/\n")
        self.git("commit", "-qam", "Old ignored local namespace")
        collision = self.write(NEW_RELEASE / "output/manifest.json", b"other manifest")
        failed = self.execute("prepare", "--target", "dev", check=False)
        self.assertEqual(failed.returncode, 1)
        self.assertIn("destination is occupied", failed.stderr)
        self.assertEqual(collision.read_bytes(), b"other manifest")
        self.assertFalse(self.state.exists())

    def test_refuses_an_occupied_canonical_destination_before_moving_payloads(self):
        self.execute("prepare", "--target", "dev")
        self.git("merge", "--no-overwrite-ignore", "--no-ff", "--no-edit", "dev")
        collision = self.write(NEW_RELEASE / "output/active.h5", b"other data")
        failed = self.execute("install", check=False)
        self.assertEqual(failed.returncode, 1)
        self.assertIn("payload identity or bytes differ", failed.stderr)
        self.assertEqual(collision.read_bytes(), b"other data")
        self.assertTrue(self.retired.exists())
        self.assertTrue(self.active.exists())
        self.assertFalse((self.root / "stage2/history").exists())

    def test_refuses_a_file_in_place_of_the_archive_before_moving_payloads(self):
        self.execute("prepare", "--target", "dev")
        self.git("merge", "--no-overwrite-ignore", "--no-ff", "--no-edit", "dev")
        collision = self.write("stage2/history", b"other data")
        failed = self.execute("install", check=False)
        self.assertEqual(failed.returncode, 1)
        self.assertIn("archive directory", failed.stderr)
        self.assertEqual(collision.read_bytes(), b"other data")
        self.assertTrue(self.retired.exists())
        self.assertTrue(self.active.exists())


if __name__ == "__main__":
    unittest.main()
