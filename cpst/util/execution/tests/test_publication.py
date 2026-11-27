import errno
import tempfile
import threading
import unittest
from pathlib import Path
from unittest.mock import patch

from cpst_util.publication import publication_lock, publish_files, replace_release


class PublicationTest(unittest.TestCase):
    def test_same_destination_excludes_other_publishers_and_recovers(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "release"
            previous = output.with_name(".release.previous")
            previous.mkdir()
            (previous / "manifest.json").write_text("accepted")
            entered = threading.Event()
            finished = threading.Event()

            def second_publisher():
                entered.set()
                with publication_lock(output):
                    finished.set()

            with publication_lock(output):
                worker = threading.Thread(target=second_publisher)
                worker.start()
                self.assertTrue(entered.wait(2))
                self.assertFalse(finished.wait(0.05))
                self.assertEqual((output / "manifest.json").read_text(), "accepted")

                def fail_replacement(candidate, destination):
                    destination.rename(previous)
                    raise OSError(errno.EIO, "simulated swap failure")

                with self.assertRaisesRegex(OSError, "simulated swap failure"):
                    replace_release(Path(temp) / "candidate", output, fail_replacement)
                self.assertEqual((output / "manifest.json").read_text(), "accepted")
            worker.join(2)
            self.assertTrue(finished.is_set())

    def test_metadata_io_failure_restores_the_accepted_files_and_allows_retry(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            release = root / "release"
            release.mkdir()
            manifest = release / "manifest.json"
            receipt = release / "receipt.json"
            archive = release / "originals.tar.gz"
            manifest.write_bytes(b"accepted manifest")
            receipt.write_bytes(b"accepted receipt")
            inodes = {path: path.stat().st_ino for path in (manifest, receipt)}
            updates = {manifest: b"candidate manifest", archive: b"candidate archive", receipt: b"candidate receipt"}
            replace = Path.replace

            def fail_after_replace(candidate, target):
                replacing_receipt = target == receipt and candidate.read_bytes() == updates[receipt]
                replace(candidate, target)
                if replacing_receipt:
                    raise OSError(errno.EIO, "simulated write failure")

            with publication_lock(release):
                with patch.object(Path, "replace", fail_after_replace):
                    with self.assertRaisesRegex(OSError, "simulated write failure"):
                        publish_files(updates, staging=root)
                self.assertEqual(manifest.read_bytes(), b"accepted manifest")
                self.assertEqual(receipt.read_bytes(), b"accepted receipt")
                self.assertEqual({path: path.stat().st_ino for path in inodes}, inodes)
                self.assertFalse(archive.exists())
                publish_files(updates, staging=root)
                self.assertEqual({path: path.read_bytes() for path in updates}, updates)
