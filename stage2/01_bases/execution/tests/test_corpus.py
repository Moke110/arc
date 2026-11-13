"""Corpus contract: roles, vehicle rows, matched controls, frozen indices, replay."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases import corpus, release

import release_fixtures


class CorpusReleaseTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.axis, self.crosswalk = release_fixtures.write_axis(self.root)
        self.fixture = release_fixtures.write_corpus(self.root, self.axis, self.crosswalk)
        self.release = self.fixture["corpus"]

    def manifest(self):
        return json.loads((self.release / "manifest.json").read_text())

    def corpus_kwargs(self):
        return dict(raw_counts_dir=self.fixture["raw_counts"],
                    axis_path=self.axis,
                    drug_metadata=self.fixture["drug_metadata"], fixture=True)

    def test_row_sets_roles_and_vehicle_rows(self):
        manifest = self.manifest()
        self.assertEqual(manifest["schema"], corpus.SCHEMA)
        self.assertEqual(manifest["sources"]["tenv"]["counts"],
                         {"profiles": 3, "treated": 1, "vehicle": 2,
                          "raw_profiles": 4})
        self.assertEqual(manifest["sources"]["tahoe_bulk"]["counts"],
                         {"profiles": 2, "treated": 1, "vehicle": 1,
                          "raw_profiles": 3})
        with corpus.open_corpus(self.release) as opened:
            tenv = opened.source("tenv")
            self.addCleanup(tenv.close)
            self.assertEqual(tenv.obs.role.tolist(),
                             ["treated", "vehicle", "vehicle"])
            self.assertEqual(tenv.obs.vehicle.tolist(), ["DMSO", "DMSO", "Media"])
            self.assertEqual(tenv.obs.plate.tolist(), ["Plate1", "Plate1", "Plate1"])
            tahoe = opened.source("tahoe_bulk")
            self.addCleanup(tahoe.close)
            self.assertEqual(tahoe.obs.role.tolist(), ["treated", "vehicle"])
            self.assertEqual(tahoe.obs.cell_line.tolist(), ["line1", "line1"])

    def test_matched_control_group_follows_the_source_rule(self):
        with corpus.open_corpus(self.release) as opened:
            tenv = opened.source("tenv")
            self.addCleanup(tenv.close)
            self.assertEqual(tenv.obs.matched_control_group.tolist(),
                             ["tenv::Plate1::DMSO", "tenv::Plate1::DMSO",
                              "tenv::Plate1::Media"])
            # The treated row resolves only to the same-plate DMSO vehicle,
            # never to the Media well.
            self.assertEqual(tenv.matched_control_rows(0).tolist(), [1])
            self.assertEqual(sorted(tenv.control_groups()), [
                "tenv::Plate1::DMSO", "tenv::Plate1::Media"])
            tahoe = opened.source("tahoe_bulk")
            self.addCleanup(tahoe.close)
            self.assertEqual(tahoe.obs.matched_control_group.tolist(),
                             ["tahoe::plate1::L1", "tahoe::plate1::L1"])
            self.assertEqual(tahoe.matched_control_rows(0).tolist(), [1])
            self.assertEqual(tahoe.control_rows("tahoe::plate1::L1").tolist(), [1])

    def test_treated_row_without_control_is_excluded_with_reason(self):
        manifest = self.manifest()
        self.assertEqual(manifest["sources"]["tahoe_bulk"]["excluded_by_reason"],
                         {"missing_matched_control": 1})
        self.assertEqual(manifest["sources"]["tenv"]["excluded_by_reason"],
                         {"technical_blank": 1})
        self.assertEqual(manifest["controls"]["tahoe_bulk"]["treated_without_control"], 1)
        exclusions = pd.read_csv(self.release / "exclusions.csv", keep_default_na=False)
        reasons = dict(zip(exclusions.source + ":" + exclusions.reason,
                           exclusions.profile_id))
        self.assertEqual(reasons["tahoe_bulk:missing_matched_control"],
                         "plate2:B:L1")
        self.assertEqual(reasons["tenv:technical_blank"], "temp_env::Plate1-A4")

    def test_corpus_has_no_normalization_or_delta(self):
        with h5py.File(self.release / "tenv.h5ad") as handle:
            self.assertEqual(handle["X"].dtype, np.dtype("int32"))
            self.assertEqual(set(handle["layers"]), set())
            obs = read_elem(handle["obs"])
            self.assertFalse({"control", "delta", "normalized"} & set(obs.columns))
            self.assertEqual(read_elem(handle["var"]).gene_index.tolist(), [0, 1])

    def test_counts_are_raw_integers_on_the_frozen_axis(self):
        with corpus.open_corpus(self.release) as opened:
            tenv = opened.source("tenv")
            self.addCleanup(tenv.close)
            self.assertEqual(tenv.counts().tolist(), [[1, 2], [3, 4], [5, 6]])
            self.assertEqual(tenv.gene_index.tolist(), [0, 1])
            tahoe = opened.source("tahoe_bulk")
            self.addCleanup(tahoe.close)
            self.assertEqual(tahoe.counts().tolist(), [[2, 7], [1, 0]])
            self.assertEqual(tahoe.gene_index.tolist(), [0, 1])
            self.assertEqual(tahoe.var.symbol.tolist(), ["A", "B"])

    def test_reader_never_returns_a_renumbered_axis(self):
        with corpus.open_corpus(self.release) as opened:
            for name in opened.sources:
                source = opened.source(name)
                self.addCleanup(source.close)
                self.assertEqual(list(source.var.gene_index),
                                 list(range(opened.axis_size)))

    def test_build_twice_produces_identical_digests(self):
        second = self.root / "second"
        corpus.build(second, **self.corpus_kwargs())
        corpus.verify(second, fixture=True)
        self.assertEqual(
            json.loads((second / "manifest.json").read_text()), self.manifest())

    def test_replay_release_reproduces_every_digest(self):
        result = corpus.replay_release(
            self.release, root=self.root / "staging", **self.corpus_kwargs())
        self.assertEqual(result["result"], "PASS")
        self.assertEqual(result["manifest_sha256"],
                         release.sha256(self.release / "manifest.json"))

    def test_manifest_has_no_absolute_paths(self):
        text = (self.release / "manifest.json").read_text()
        self.assertNotIn(str(self.root), text)
        self.assertNotIn("/tmp/", text)
        release.check_relative_paths(self.manifest())

    def test_tampered_payload_fails_verification(self):
        h5ad = self.release / "tenv.h5ad"
        with h5py.File(h5ad, "a") as handle:
            handle["X"][0, 0] = 999
        with self.assertRaises(release.ReleaseError):
            corpus.verify(self.release, fixture=True)


if __name__ == "__main__":
    unittest.main()
