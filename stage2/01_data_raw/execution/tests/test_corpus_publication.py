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
              patch.object(build_corpus, "REPO", self.root),
              patch.object(build_corpus, "HISTORY", self.root / "history")):
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

    def test_different_predecessor_cannot_replace_the_accepted_legacy_release(self):
        with (patch.object(build_corpus, "migrate") as build,
              self.assertRaisesRegex(CorpusError, "installed legacy corpus")):
            self.migrate(self.predecessor)
        build.assert_not_called()
        self.assertTrue((self.legacy / "manifest.json").is_file())

    def test_unarchived_legacy_refuses_migration_before_build(self):
        with (patch.object(build_corpus, "migrate") as build,
              self.assertRaisesRegex(CorpusError, "verified archive handoff")):
            self.migrate(self.legacy)
        build.assert_not_called()
        self.assertTrue(self.legacy.exists())
        self.assertFalse((self.root / "history").exists())

    def test_failed_downstream_read_leaves_accepted_legacy_release_and_api(self):
        self.archived_legacy()
        api = self.output / "api"
        api.mkdir()
        reader = api / "reader.py"
        reader.write_text("# accepted reader\n")
        with (patch.object(build_corpus, "migrate"),
              patch.object(build_corpus, "verify_product"),
              patch.object(build_corpus, "verify_equivalence"),
              patch.object(build_corpus, "verify_downstream",
                           side_effect=CorpusError("failed downstream read")),
              self.assertRaisesRegex(CorpusError, "failed downstream read")):
            self.migrate(self.legacy)
        self.assertTrue((self.legacy / "manifest.json").is_file())
        self.assertEqual(reader.read_text(), "# accepted reader\n")

    def test_verified_archive_allows_migration_and_keeps_archive_read_only(self):
        archive = self.archived_legacy()
        self.predecessor.joinpath("manifest.json").write_bytes((archive / "manifest.json").read_bytes())
        inventory = self.root / "history/archive-manifest.json"
        before = {path: sha256(path) for path in (inventory, *archive.iterdir())}
        api = self.output / "api"
        api.mkdir()
        reader = api / "reader.py"
        reader.write_text("# accepted reader\n")

        def stage(predecessor, candidate):
            self.assertEqual(predecessor.root, self.legacy)
            (candidate / "manifest.json").write_text("successor")

        with (patch.object(build_corpus, "migrate", side_effect=stage),
              patch.object(build_corpus, "verify_product", return_value="successor"),
              patch.object(build_corpus, "verify_equivalence"),
              patch.object(build_corpus, "verify_downstream") as downstream):
            self.migrate(self.predecessor)
        downstream.assert_called_once()
        self.assertFalse(self.legacy.exists())
        self.assertEqual({path: sha256(path) for path in before}, before)
        self.assertEqual(reader.read_text(), "# accepted reader\n")

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
