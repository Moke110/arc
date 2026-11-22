"""Build/verify seam: the `full` view is a digest-equal alias of the accepted release."""
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXECUTION = HERE.parent
for path in (str(EXECUTION), str(HERE)):
    if path not in sys.path:
        sys.path.insert(0, path)

import build  # noqa: E402
import tsg_view_fixtures as fixtures  # noqa: E402
import views  # noqa: E402

REPO = Path(__file__).resolve().parents[4]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class FullViewAliasTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        views.bind(REPO)

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        fixture = fixtures.grid_fixture(self.root)
        self.accepted, self.accepted_manifest, self.fixture = fixtures.build_view(
            fixture, size=2)

    def build(self, output=None, staging=None, source="tenv"):
        return build.build_product(
            view="full", source=source,
            output=output or self.root / "views",
            staging=staging or self.root / "staging",
            repo_root_path=REPO, size=2, fixture=True,
            corpus=self.fixture["corpus"], axis=self.fixture["axis"],
            accepted=self.accepted, provenance_root=None)

    def test_full_view_is_published_as_a_digest_equal_alias(self):
        receipt = self.build()
        accepted_manifest = json.loads((self.accepted / "manifest.json").read_text())
        product = json.loads((self.root / "views/full/tenv.json").read_text())
        self.assertTrue(receipt["evidence"]["digest_equal"])
        self.assertEqual(receipt["evidence"]["accepted_manifest_sha256"],
                         sha256(self.accepted / "manifest.json"))
        self.assertEqual(receipt["evidence"]["rebuilt_manifest_sha256"],
                         receipt["evidence"]["accepted_manifest_sha256"])
        self.assertEqual(product["schema"], "tsg-view-product-v1")
        self.assertEqual(product["product"], "tenv-full-5k-view")
        self.assertEqual(product["view"], "full")
        self.assertIsNone(product["payload"])
        self.assertEqual(product["alias"]["manifest_sha256"],
                         receipt["evidence"]["accepted_manifest_sha256"])
        for name in ("tenv", "tahoe_bulk"):
            self.assertEqual(product["alias"]["payload_sha256"][name],
                             accepted_manifest["outputs"][f"{name}.h5ad"]["sha256"])
        self.assertEqual(product["rows"], 20)
        self.assertEqual(product["treated"], 16)
        self.assertEqual(product["conditions"], 9)
        self.assertEqual(sorted(path.name for path in (self.root / "views/full").iterdir()),
                         ["tenv.json", "tenv.profile_ids.csv"])
        self.assertFalse(any((self.root / "views/full").glob("*.h5ad")))
        views.modules()["release"].check_relative_paths(product)
        # No rebuild staging or duplicate payload survives a successful build.
        self.assertEqual([path.name for path in (self.root / "staging").iterdir()], [])

    def test_full_view_fails_closed_when_the_accepted_payload_drifts(self):
        payload = self.accepted / "tenv.h5ad"
        with payload.open("ab") as handle:
            handle.write(b"drift")
        with self.assertRaisesRegex(views.HarnessError,
                                    "release payload digest differs: tenv.h5ad"):
            self.build()
        self.assertFalse((self.root / "views/full/tenv.json").exists())

    def test_full_view_fails_closed_when_the_accepted_manifest_changes(self):
        manifest_path = self.accepted / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["recipe"] = manifest["recipe"] + " drift"
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        with self.assertRaisesRegex(views.HarnessError, "full-view rebuild differs"):
            self.build()

    def test_both_sources_publish_into_the_same_alias_directory(self):
        self.build(source="tenv")
        self.build(source="tahoe_bulk")
        names = sorted(path.name for path in (self.root / "views/full").iterdir())
        self.assertEqual(names, ["tahoe_bulk.json", "tahoe_bulk.profile_ids.csv",
                                 "tenv.json", "tenv.profile_ids.csv"])


if __name__ == "__main__":
    unittest.main()
