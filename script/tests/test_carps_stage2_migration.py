"""Exercise the same-filesystem CarPS cutover through its CLI."""
import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "migrate_carps_stage2.py"
OLD = "stage2"
NEW = "stage2_CarPS-VAE"


class Stage2CutoverTest(unittest.TestCase):
    def test_cutover_preserves_archive_and_links_and_updates_tracked_files(self):
        scratch = SCRIPT.parent.parent / "tmp"
        scratch.mkdir(exist_ok=True)
        with TemporaryDirectory(dir=scratch) as repo:
            self.check_cutover(Path(repo))

    def check_cutover(self, root):
        def run(*args):
            result = subprocess.run(args, cwd=root, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            return result.stdout

        def git(*args):
            return run("git", "-c", "user.name=Fixture", "-c", "user.email=fixture@example.com", *args)

        git("init", "-q", "-b", "main")
        stage = root / OLD
        history = stage / "history" / OLD
        stage.mkdir()
        (root / ".gitignore").write_text(f"/tmp/\n/{OLD}/history/\n/*.pt\n*.bin\n")
        tracked = stage / "entry.py"
        tracked.write_text("old route")
        pinned = stage / "receipt.json"
        pinned.write_text('{"path":"stage2/old"}\n')
        unrelated = root / "unrelated.txt"
        unrelated.write_text("keep the local deletion")
        git("add", ".")
        git("commit", "-qm", "Predecessor")
        git("switch", "-qc", "dev")
        git("mv", OLD, NEW)
        (root / NEW / tracked.name).write_text("new route")
        (root / NEW / "new.py").write_text("new entry")
        (root / ".gitignore").write_text(f"/tmp/\n/{NEW}/history/\n/*.pt\n*.bin\n")
        git("add", ".")
        git("commit", "-qm", "Canonical Stage 2")
        git("switch", "-q", "main")

        history.mkdir(parents=True, exist_ok=True)
        archive = history / "archive.txt"
        archive.write_bytes(b"preserved archive\0")
        weight = root / "external.pt"
        weight.write_bytes(b"external checkpoint")
        link = stage / "checkpoint.bin"
        link.symlink_to(weight)
        unrelated.unlink()
        old_dir, old_receipt, old_archive, old_link = (
            path.lstat().st_ino for path in (stage, pinned, archive, link)
        )
        state = root / "tmp/cutover"
        result = json.loads(run(sys.executable, str(SCRIPT), "--root", str(root), "--work", str(state),
                                "deploy", "--target", "dev"))
        self.assertEqual(result["status"], "PASS")
        moved = root / NEW
        self.assertEqual(moved.lstat().st_ino, old_dir)
        self.assertEqual((moved / "history" / OLD / archive.name).lstat().st_ino, old_archive)
        self.assertEqual((moved / pinned.name).lstat().st_ino, old_receipt)
        self.assertEqual((moved / link.name).lstat().st_ino, old_link)
        self.assertEqual(os.readlink(moved / link.name), str(weight))
        self.assertEqual((moved / link.name).read_bytes(), b"external checkpoint")
        self.assertEqual((moved / tracked.name).read_text(), "new route")
        self.assertEqual((moved / "new.py").read_text(), "new entry")
        self.assertFalse(stage.exists())
        self.assertEqual(git("status", "--porcelain"), " D unrelated.txt\n")
        command = [sys.executable, str(SCRIPT), "--root", str(root), "--work", str(state), "verify"]
        self.assertEqual(json.loads(run(*command))["status"], "PASS")
        archived_file = moved / "history" / OLD / archive.name
        archived_file.write_bytes(b"modified archive")
        result = subprocess.run(command, cwd=root, capture_output=True, text=True)
        self.assertEqual(result.returncode, 1)
        self.assertIn(str(archived_file.relative_to(moved)), result.stderr)


if __name__ == "__main__":
    unittest.main()
