import errno
import tempfile
import threading
import unittest
from pathlib import Path

from cpst_util.publication import publication_lock, replace_release


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
