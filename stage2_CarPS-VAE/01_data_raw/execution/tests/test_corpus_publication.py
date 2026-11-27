import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import build_corpus
from native_corpus import CorpusError, sha256


class CorpusPublicationTest(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.output = self.root / "output"
        self.legacy = self.output / "native"
        self.legacy.mkdir(parents=True)
        self.predecessor = self.root / "different-legacy"
        self.predecessor.mkdir()
        for directory, fingerprint in ((self.legacy, "accepted"), (self.predecessor, "different")):
            (directory / "manifest.json").write_text(json.dumps({
                "schema": "stage2-native-corpus-v1", "fingerprint": fingerprint,
            }))

    def run_corpus(self, action, *options):
        argv = ["build_corpus.py", action, "--output", str(self.output),
                "--staging", str(self.root / "staging"), *options]
        with (patch.object(sys, "argv", argv),
              patch.object(build_corpus, "DEFAULT", self.output),
              patch.object(build_corpus, "REPO", self.root)):
            build_corpus.main()

    def migrate(self, predecessor):
        self.run_corpus("migrate", "--predecessor", str(predecessor))

    def archived_legacy(self):
        payload = self.legacy / "payload"
        payload.write_text("accepted")
        manifest = self.legacy / "manifest.json"
        manifest.write_text(json.dumps({
            "schema": "stage2-native-corpus-v1", "outputs": {"payload": sha256(payload)},
        }))
        history = self.root / "history"
        original = self.legacy.relative_to(self.root)
        destination = Path("history") / original
        archive = self.root / destination
        archive.mkdir(parents=True)
        archived = {}
        for path in (payload, manifest):
            target = archive / path.name
            target.write_bytes(path.read_bytes())
            archived[(original / path.name).as_posix()] = {
                "destination": (destination / path.name).as_posix()
            }
        (history / "archive-manifest.json").write_text(json.dumps({
            "schema": "stage2-local-history-v1",
            "path_map": {original.as_posix(): destination.as_posix()},
            "entries": archived,
        }))
        return archive

    def test_verify_migrated_release_needs_no_archived_predecessor(self):
        self.legacy.rename(self.root / "unavailable-predecessor")
        (self.output / "manifest.json").write_text(json.dumps({
            "inputs": {}, "predecessor": {"manifest_sha256": "past"},
        }))
        argv = ["build_corpus.py", "verify", "--output", str(self.output),
                "--staging", str(self.root / "staging")]
        with (patch.object(sys, "argv", argv),
              patch.object(build_corpus, "DEFAULT", self.output),
              patch.object(build_corpus, "verify_product", return_value="accepted") as verify):
            build_corpus.main()
        verify.assert_called_once_with(self.output)

    def test_archive_cannot_authorize_replacement_of_an_installed_legacy_corpus(self):
        archive = self.archived_legacy()
        manifest = self.legacy / "manifest.json"
        saved = sha256(manifest)
        with (patch.object(build_corpus, "publish_candidate") as build,
              self.assertRaisesRegex(CorpusError, "installed legacy corpus.*scratch output")):
            self.migrate(self.legacy)
        build.assert_not_called()
        self.assertEqual(sha256(manifest), saved)
        self.assertEqual(sha256(archive / "manifest.json"), saved)

    def test_installed_publish_verifies_pinned_generation_without_build(self):
        self.legacy.rename(self.root / "unavailable-legacy")
        manifest = self.output / "manifest.json"
        manifest.write_text(json.dumps({
            "schema": "stage2-native-corpus-v2", "predecessor": {"manifest_sha256": "past"},
        }))
        expected = sha256(manifest)
        with (patch.object(build_corpus, "verify_product", return_value=expected) as verify,
              patch.object(build_corpus, "verify_downstream") as downstream,
              patch.object(build_corpus, "publish_candidate") as candidate):
            self.run_corpus("publish")
        verify.assert_called_once_with(self.output)
        downstream.assert_called_once_with(self.output, self.output)
        candidate.assert_not_called()
        self.assertEqual(sha256(manifest), expected)

    def test_installed_replacement_requires_scratch_output_before_build(self):
        self.legacy.rename(self.root / "unavailable-legacy")
        (self.output / "manifest.json").write_text('{"schema":"stage2-native-corpus-v2"}')
        with (patch.object(build_corpus, "publish_candidate") as candidate,
              self.assertRaisesRegex(CorpusError, "scratch output")):
            self.run_corpus("publish", "--tahoe-built", str(self.root / "other"))
        candidate.assert_not_called()
