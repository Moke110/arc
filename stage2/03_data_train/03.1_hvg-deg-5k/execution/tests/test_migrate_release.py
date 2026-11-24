"""Regression tests for preserving the release and predecessor during staging."""
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from migrate_release import stage_migration


class MigrationPathsTest(unittest.TestCase):
    def setUp(self):
        scratch = self.enterContext(tempfile.TemporaryDirectory())
        self.root = Path(scratch)
        self.predecessor = self.root / "predecessor"
        self.predecessor.mkdir()
        self.manifest = self.predecessor / "manifest.json"
        self.manifest.write_bytes(b"sealed predecessor")
        self.candidate = self.root / "candidate"

    def test_rejects_receipt_inside_candidate(self):
        with self.assertRaisesRegex(ValueError, "receipt must be outside"):
            stage_migration(self.predecessor, self.candidate, self.candidate / "manifest.json")
        self.assertFalse(self.candidate.exists())

    def test_rejects_receipt_alias_inside_predecessor(self):
        alias = self.root / "alias"
        alias.symlink_to(self.predecessor, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, "receipt must be outside"):
            stage_migration(self.predecessor, self.candidate, alias / "manifest.json")
        self.assertEqual(self.manifest.read_bytes(), b"sealed predecessor")
        self.assertFalse(self.candidate.exists())

    def test_rejects_nested_candidate(self):
        with self.assertRaisesRegex(ValueError, "predecessor and candidate must be separate"):
            stage_migration(self.predecessor, self.predecessor / "candidate",
                            self.root / "receipt.json")
        self.assertEqual(list(self.predecessor.iterdir()), [self.manifest])

    def test_rejects_existing_external_receipt(self):
        receipt = self.root / "receipt.json"
        receipt.write_bytes(b"previous receipt")
        with self.assertRaises(FileExistsError):
            stage_migration(self.predecessor, self.candidate, receipt)
        self.assertEqual(receipt.read_bytes(), b"previous receipt")
        self.assertFalse(self.candidate.exists())
