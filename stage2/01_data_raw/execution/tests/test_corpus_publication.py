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

    def migrate(self, predecessor):
        argv = ["build_corpus.py", "migrate", "--predecessor", str(predecessor),
                "--output", str(self.output), "--staging", str(self.root / "staging")]
        with (patch.object(sys, "argv", argv),
              patch.object(build_corpus, "DEFAULT", self.output)):
            build_corpus.main()

    def test_different_predecessor_cannot_replace_the_accepted_legacy_release(self):
        with (patch.object(build_corpus, "migrate") as build,
              self.assertRaisesRegex(CorpusError, "installed legacy corpus")):
            self.migrate(self.predecessor)
        build.assert_not_called()
        self.assertTrue((self.legacy / "manifest.json").is_file())

    def test_failed_downstream_read_leaves_accepted_legacy_release_and_api(self):
        api = self.output / "api"
        api.mkdir()
        reader = api / "reader.py"
        reader.write_text("# accepted reader\n")
        with (patch.object(build_corpus, "migrate"),
              patch.object(build_corpus, "verify_product"),
              patch.object(build_corpus, "verify_equivalence"),
              patch.object(build_corpus, "verify_downstream",
                           side_effect=CorpusError("failed downstream read"), create=True),
              self.assertRaisesRegex(CorpusError, "failed downstream read")):
            self.migrate(self.legacy)
        self.assertTrue((self.legacy / "manifest.json").is_file())
        self.assertEqual(reader.read_text(), "# accepted reader\n")

    def test_matching_external_predecessor_preserves_actual_release_and_api(self):
        payload = self.legacy / "payload"
        payload.write_text("accepted")
        manifest = json.dumps({
            "schema": "stage2-native-corpus-v1", "outputs": {"payload": sha256(payload)},
        })
        for root in (self.legacy, self.predecessor):
            (root / "manifest.json").write_text(manifest)
        api = self.output / "api"
        api.mkdir()
        (api / "reader.py").write_text("# accepted reader\n")

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
        archived = self.root / "historical/native-v1"
        self.assertEqual((archived / "manifest.json").read_text(), manifest)
        self.assertEqual((archived / "payload").read_text(), "accepted")
        self.assertEqual((api / "reader.py").read_text(), "# accepted reader\n")
