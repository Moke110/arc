"""Worked examples through published view builder and reader."""
import json
import sys
import tempfile
import unittest
from hashlib import sha256
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData

STEP = Path(__file__).resolve().parents[2]
STAGE = STEP.parents[1]
sys.path[:0] = [str(item) for item in (
    STEP / "execution", STEP / "output/api",
    STAGE / "01_data_raw/output/api",
    STAGE / "02_gene_axis/02.1_hvg-deg/output/api",
)]

from native_corpus import open_corpus, seal_corpus
from prediction_axis import COMMON_SHA256, PredictionAxis
from build_views import publish_views, replay_views
from expression_view import PAYLOAD_SCHEMA, open_view
from migrate_release import stage_migration


class ExpressionViewTest(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        corpus = self.root / "counts"
        corpus.mkdir()
        roles = ["treated", "vehicle", "vehicle", "treated", "treated", "vehicle", "treated"]
        obs = pd.DataFrame({
            "role": roles, "matched_control_group": ["A", "A", "A", "A", "B", "B", ""],
            "cell_line": "CM", "drug_parent_inchikey": ["P" if r == "treated" else "" for r in roles],
            "dose_um": [1.0 if r == "treated" else 0.0 for r in roles],
            "view_exclusion_reason": "",
        }, index=pd.Index(["t1", "v1", "v2", "t2", "t3", "v3", "unmatched"], name="profile_id"))
        vectors = np.zeros((7, 18164), dtype=np.int32)
        vectors[:, 1] = 100000
        vectors[:, [5, 2]] = [[3, 0], [1, 0], [0, 9], [0, 2], [1, 4], [10, 0], [1, 1]]
        for source in ("tenv", "tahoe_bulk"):
            AnnData(X=vectors, obs=obs, var=pd.DataFrame(index=[str(i) for i in range(18164)])).write_h5ad(corpus / f"{source}.h5ad")
        seal_corpus(corpus, dict.fromkeys(("tenv", "tahoe_bulk"), dict(enumerate(range(18164)))),
                    18164, COMMON_SHA256, {}, {})
        self.corpus = self.enterContext(open_corpus(corpus))
        self.axis = PredictionAxis((5, 2), ("five", "two"), "a" * 64)
        self.output = self.root / "views"

    def test_raw_first_control_and_individual_vehicle_rows(self):
        publish_views(self.output, self.corpus, self.axis)
        with open_view(self.output, source="tenv") as view:
            self.assertEqual(list(view.obs.index), ["t1", "v1", "v2", "t2", "t3", "v3"])
            self.assertEqual(list(view.var.gene_index), [5, 2])
            np.testing.assert_allclose(view.read(0, 1, "control")[0], np.log1p([1000, 9000]))
            np.testing.assert_allclose(view.read(1, 3), np.log1p([[10000, 0], [0, 10000]]))
            self.assertEqual(view.exclusions.set_index("profile_id").loc["unmatched", "reason"],
                             "no_matched_vehicle")
            np.testing.assert_allclose(view.read(0, 1, "delta"), view.read(0, 1) - view.read(0, 1, "control"))

    def test_flat_release_is_replayable(self):
        publish_views(self.output, self.corpus, self.axis)
        self.assertEqual({path.name for path in self.output.iterdir()}, {
            "manifest.json", "tenv.h5ad", "tahoe_bulk.h5ad", "tenv.members.csv",
            "tahoe_bulk.members.csv", "tenv.exclusions.csv", "tahoe_bulk.exclusions.csv",
        })
        for source in ("tenv", "tahoe_bulk"):
            with open_view(self.output, source=source, family="full") as view:
                self.assertEqual(view.shape, (6, 2))
        with self.assertRaisesRegex(ValueError, "unknown expression view arm"):
            open_view(self.output, source="tenv", family="merge")
        self.assertEqual(replay_views(self.output, self.corpus, self.axis)["result"], "PASS")

    def test_migration_preserves_payloads_and_cleans_up_when_receipt_cannot_be_written(self):
        genes = tuple(range(5000))
        axis = PredictionAxis(genes, tuple(map(str, genes)), self.axis.sha256)
        publish_views(self.output, self.corpus, axis)
        manifest_path = self.output / "manifest.json"
        predecessor = json.loads(manifest_path.read_bytes())
        predecessor["schema"] = PAYLOAD_SCHEMA
        (self.output / "full").mkdir()
        for name in list(predecessor["outputs"]):
            if name.endswith((".h5ad", ".members.csv")):
                old_name = f"full/{name}"
                (self.output / name).rename(self.output / old_name)
                predecessor["outputs"][old_name] = predecessor["outputs"].pop(name)
        raw_manifest = json.dumps(predecessor).encode()
        manifest_path.write_bytes(raw_manifest)

        candidate = self.root / "candidate"
        occupied = self.root / "file"
        occupied.touch()
        with self.assertRaises(OSError):
            stage_migration(self.output, candidate, occupied / "receipt.json")
        self.assertFalse(candidate.exists())
        self.assertEqual(manifest_path.read_bytes(), raw_manifest)

        receipt_path = self.root / "receipt.json"
        stage_migration(self.output, candidate, receipt_path)
        receipt = json.loads(receipt_path.read_bytes())
        self.assertEqual(receipt["previous_manifest"]["sha256"],
                         sha256(raw_manifest).hexdigest())
        self.assertEqual(receipt["release_manifest"]["sha256"],
                         sha256((candidate / "manifest.json").read_bytes()).hexdigest())
        self.assertEqual((candidate / "tenv.h5ad").read_bytes(),
                         (self.output / "full/tenv.h5ad").read_bytes())

    def test_reader_rejects_stale_digest_and_modified_payload(self):
        publish_views(self.output, self.corpus, self.axis)
        with self.assertRaisesRegex(ValueError, "generation"):
            open_view(self.output, source="tenv", corpus_sha256="0" * 64)
        with open_view(self.output, source="tenv") as view:
            with self.assertRaisesRegex(ValueError, "expression window"):
                view.read(0, 257)
        path = self.output / "tenv.h5ad"
        with path.open("ab") as output:
            output.write(b"unexpected")
        with self.assertRaisesRegex(ValueError, "payload digest"):
            open_view(self.output, source="tenv")

    def test_canonical_reader_rejects_stale_corpus_and_coherent_wrong_recipe(self):
        from unittest.mock import patch
        import expression_view

        publish_views(self.output, self.corpus, self.axis)
        corpus = self.root / "accepted-corpus.json"
        corpus.write_text('{"schema":"stage2-native-corpus-v1"}')
        path = self.output / "manifest.json"
        manifest = json.loads(path.read_bytes())
        with (patch.object(expression_view, "DEFAULT", self.output),
              patch.object(expression_view, "CORPUS_MANIFEST", corpus)):
            with self.assertRaisesRegex(ValueError, "requested corpus or axis"):
                open_view(source="tenv")
        manifest["normalization"] = "average normalized vehicles"
        path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError, "expression view recipe"):
            open_view(self.output, source="tenv", corpus_sha256=self.corpus.manifest_sha256)

    def test_canonical_reader_keeps_the_sealed_view_pin_after_migration(self):
        from unittest.mock import patch
        import expression_view

        from shutil import copyfile

        archived = self.root / "historical/native-v1"
        archived.mkdir(parents=True)
        manifest = {
            **self.corpus.manifest,
            "schema": "stage2-native-corpus-v1", "outputs": {},
            "sources": {name: {key: value for key, value in source.items()
                               if key != "manifest_sha256"}
                        for name, source in self.corpus.manifest["sources"].items()},
        }
        for name, source in manifest["sources"].items():
            file = f"{name}.h5ad"
            copyfile(self.corpus.root / name / file, archived / file)
            manifest["outputs"][file] = source["sha256"]
        archived_manifest = archived / "manifest.json"
        archived_manifest.write_text(json.dumps(manifest))
        axis = self.root / "gene_axis_5000.csv"
        axis.write_text("fixture")
        digest = sha256(axis.read_bytes()).hexdigest()
        axis.with_name("manifest.json").write_text(json.dumps({"axes": {"5000": digest}}))
        accepted_axis = PredictionAxis(self.axis.gene_indices, self.axis.symbols, digest)
        with open_corpus(archived) as previous:
            old_pin = previous.manifest_sha256
            publish_views(self.output, previous, accepted_axis)
        successor_manifest = self.corpus.root / "manifest.json"
        manifest = self.corpus.manifest
        manifest["predecessor"] = {
            "manifest_sha256": old_pin,
            "unchanged_sources": {"tenv": "tenv", "tahoe_bulk": "tahoe_bulk"},
        }
        successor_manifest.write_text(json.dumps(manifest))
        successor = self.enterContext(open_corpus(self.corpus.root))
        self.assertEqual(replay_views(self.output, successor, accepted_axis)["result"], "PASS")
        moved_archive = archived.rename(self.root / "portable-v1")
        with self.assertRaisesRegex(ValueError, "requested corpus or axis"):
            replay_views(self.output, successor, accepted_axis)
        self.assertEqual(replay_views(self.output, successor, accepted_axis,
                                     historical_manifest=moved_archive / "manifest.json")["result"], "PASS")
        moved_archive.rename(archived)
        with (patch.object(expression_view, "DEFAULT", self.output),
              patch.object(expression_view, "CORPUS_MANIFEST", successor_manifest),
              patch.object(expression_view, "HISTORICAL_CORPUS_MANIFEST", archived_manifest),
              patch.object(expression_view, "CANONICAL_AXIS", axis)):
            with open_view(source="tenv") as view:
                self.assertEqual(view.read(0, 1).shape, (1, 2))
            archived_manifest.write_text("{}")
            with self.assertRaisesRegex(ValueError, "requested corpus or axis"):
                open_view(source="tenv")
