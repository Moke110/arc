"""Raw-count release contract: integer payloads, bridge evidence, seals, replay.

Every assertion is against the public build/replay/reader boundary on a tiny
hermetic source tree; the multi-GB accepted products are the acceptance step,
not this test.
"""
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
from bases import raw_counts, release

import release_fixtures


class RawCountsReleaseTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.axis, self.crosswalk = release_fixtures.write_axis(self.root)
        self.fixture = release_fixtures.write_raw_release(
            self.root, self.axis, self.crosswalk)
        self.release = self.fixture["raw_counts"]

    def rebuild_kwargs(self):
        return dict(
            mode="source-accepted-replay",
            data_root=self.fixture["data_root"],
            axis_path=self.axis,
            crosswalk_path=self.crosswalk,
            tenv_payload=self.fixture["tenv_payload"],
            tahoe_payload=self.fixture["tahoe_payload"],
            fixture=True,
        )

    def manifest(self):
        return json.loads((self.release / "manifest.json").read_text())

    def test_release_is_integer_counts_without_normalization(self):
        manifest = self.manifest()
        self.assertEqual(manifest["schema"], raw_counts.SCHEMA)
        self.assertEqual(manifest["counts"]["axis_size"], 2)
        self.assertEqual(manifest["sources"]["tenv"]["counts"],
                         {"profiles": 4, "genes": 2, "treated": 1,
                          "vehicle": 2, "blank": 1})
        self.assertEqual(manifest["sources"]["tahoe"]["counts"],
                         {"profiles": 3, "genes": 2, "treated": 2,
                          "vehicle": 1, "cells": 3})
        with h5py.File(self.release / "tenv.h5ad") as handle:
            self.assertEqual(handle["X"].dtype, np.dtype("int32"))
            self.assertEqual(tuple(handle["X"].shape), (4, 2))
            self.assertEqual(set(handle["layers"]), set())
            self.assertEqual(read_elem(handle["obs"]).is_blank.tolist(),
                             [False, False, False, True])
        with h5py.File(self.release / "tahoe.h5ad") as handle:
            self.assertEqual(handle["X"]["data"].dtype, np.dtype("int32"))
            self.assertEqual(tuple(int(v) for v in handle["X"].attrs["shape"]), (3, 2))
            self.assertEqual(set(handle["layers"]), set())
        obs = pd.read_csv(self.release / "tenv_axis_mapping.csv")
        self.assertNotIn("control", obs.columns)
        self.assertNotIn("delta", obs.columns)

    def test_reader_returns_counts_role_and_bridge_evidence(self):
        with raw_counts.open_tenv(self.release) as reader:
            self.assertEqual(reader.shape, (4, 2))
            self.assertEqual(reader.counts().tolist(), [[1, 2], [3, 4], [5, 6], [0, 1]])
            bridge = reader.bridge()
            self.assertEqual(bridge.gene_index.tolist(), [0, 1])
            self.assertEqual(bridge.native_column.tolist(), [0, 1])
            self.assertTrue(bridge.resolution.eq("exact_hgnc").all())
        with raw_counts.open_tahoe(self.release) as reader:
            self.assertEqual(reader.shape, (3, 2))
            self.assertEqual(reader.counts().tolist(), [[2, 7], [1, 0], [8, 0]])
            mapping = reader.bridge()
            self.assertEqual(mapping.gene_index.tolist(), [0, 1])
            self.assertEqual(mapping.native_token_id.tolist(), ["3", "5"])
            self.assertEqual(mapping.native_ensembl_id.tolist(), ["E1", "E2"])
            self.assertEqual(reader.obs.is_control.tolist(), [False, True, False])

    def test_build_twice_produces_identical_digests(self):
        accepted = {name: record["sha256"]
                    for name, record in self.manifest()["outputs"].items()}
        second = self.root / "second"
        raw_counts.build(second, **self.rebuild_kwargs())
        raw_counts.verify(second)
        self.assertEqual(
            {name: record["sha256"] for name, record in
             json.loads((second / "manifest.json").read_text())["outputs"].items()},
            accepted)
        self.assertEqual(
            json.loads((second / "manifest.json").read_text()), self.manifest())

    def test_replay_release_reproduces_every_digest(self):
        result = raw_counts.replay_release(
            self.release, root=self.root / "staging", **self.rebuild_kwargs())
        self.assertEqual(result["result"], "PASS")
        self.assertEqual(result["manifest_sha256"],
                         release.sha256(self.release / "manifest.json"))
        self.assertEqual(sorted(result["outputs"]), sorted(
            self.manifest()["outputs"]))

    def test_manifest_has_no_absolute_paths(self):
        text = (self.release / "manifest.json").read_text()
        self.assertNotIn(str(self.root), text)
        self.assertNotIn("/tmp/", text)
        release.check_relative_paths(self.manifest())
        with self.assertRaises(release.ReleaseError):
            release.check_relative_paths({"inputs": {"/etc/passwd": {}}})

    def test_tampered_payload_fails_verification(self):
        path = self.release / "token_mapping.csv"
        original = path.read_text()
        path.write_text(original.replace("E1", "E9", 1))
        with self.assertRaises(release.ReleaseError):
            raw_counts.verify(self.release)

    def test_payload_digest_guard_rejects_unaccepted_bytes(self):
        candidate = self.root / "candidate"
        with self.assertRaises(raw_counts.RawCountsError):
            raw_counts._copy_payload(
                candidate, self.fixture["tenv_payload"], ("tenv.h5ad",),
                {"tenv.h5ad": "0" * 64}, "tenv")

    def test_full_rebuild_mode_rebuilds_the_same_payloads_from_source(self):
        candidate = self.root / "full"
        raw_counts.build(candidate, mode="source-rebuild",
                         data_root=self.fixture["data_root"],
                         axis_path=self.axis, crosswalk_path=self.crosswalk,
                         tenv_payload=self.fixture["tenv_payload"],
                         tahoe_payload=self.fixture["tahoe_payload"], fixture=True,
                         gene_metadata=self.fixture["data_root"] /
                         "Tahoe-100M/source/metadata/gene_metadata.parquet",
                         tahoe_expected_shards=2, tahoe_expected_genes=2)
        raw_counts.verify(candidate)
        for name in raw_counts.PRODUCT_FILES:
            self.assertEqual(release.sha256(candidate / name),
                             release.sha256(self.release / name), name)


if __name__ == "__main__":
    unittest.main()
