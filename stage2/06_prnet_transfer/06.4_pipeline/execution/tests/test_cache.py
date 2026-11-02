"""Derived Tahoe log1p memmap cache: full-array checksum, sample check, fail-closed use."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

import _support
from ptlib import cache as cache_mod
from ptlib import data as data_mod
from ptlib import h5io


class CacheTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = _support.small_fixture_dir()
        cls.source = cls.dir / "tahoe_transfer.h5ad"
        cls.source_sha = _support.artifact_hashes(cls.dir)["tahoe_transfer.h5ad"]
        cls.tmp = tempfile.TemporaryDirectory()
        cls.cache_path = Path(cls.tmp.name) / "tahoe_log1p.npy"
        cls.record = cache_mod.build_cache(cls.source, cls.cache_path, cls.source_sha)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_build_record(self):
        self.assertEqual(self.record["shape"], [128, 512])
        self.assertEqual(self.record["dtype"], "float32")
        self.assertEqual(self.record["sample_max_abs_diff"], 0.0)
        self.assertEqual(self.record["source_sha256"], self.source_sha)
        self.assertEqual(self.record["bytes"], self.cache_path.stat().st_size)
        self.assertTrue(self.cache_path.exists())
        sidecar = json.loads(Path(str(self.cache_path) + ".json").read_text())
        self.assertEqual(sidecar["layer_sha256"], self.record["layer_sha256"])
        self.assertNotEqual(sidecar["file_sha256"], None)

    def test_cache_values_equal_the_frozen_layer(self):
        opened = cache_mod.open_cache(
            self.cache_path, self.source_sha, self.source, verify="deep"
        )
        rows = np.array([0, 5, 17, 100])
        frozen, _ = h5io.read_layer_rows(self.source, rows, 512)
        np.testing.assert_array_equal(opened.rows(rows), frozen)
        self.assertEqual(opened.sample_max_abs_diff, 0.0)

    def test_wrong_source_hash_fails_closed(self):
        with self.assertRaisesRegex(cache_mod.CacheError, "different frozen artifact"):
            cache_mod.open_cache(self.cache_path, "0" * 64, self.source)

    def test_corrupted_cache_fails_deep_verification(self):
        copy = Path(self.tmp.name) / "corrupt.npy"
        copy.write_bytes(self.cache_path.read_bytes())
        Path(str(copy) + ".json").write_text(
            Path(str(self.cache_path) + ".json").read_text()
        )
        data = bytearray(copy.read_bytes())
        data[-1] ^= 0xFF
        copy.write_bytes(bytes(data))
        with self.assertRaisesRegex(cache_mod.CacheError, "checksum"):
            cache_mod.open_cache(copy, self.source_sha, self.source, verify="deep")

    def test_hash_log1p_layer_matches_build_checksum(self):
        self.assertEqual(
            cache_mod.hash_log1p_layer(self.source), self.record["layer_sha256"]
        )

    def test_rebind_to_a_byte_identical_copy(self):
        copy_dir = Path(self.tmp.name) / "copy"
        copy_dir.mkdir()
        import shutil

        shutil.copy(self.source, copy_dir / "tahoe_transfer.h5ad")
        new_sha = "b" * 64  # e.g. the frozen input has its own file hash
        record = cache_mod.rebind_cache(self.cache_path, copy_dir / "tahoe_transfer.h5ad", new_sha)
        self.assertEqual(record["source_sha256"], new_sha)
        self.assertEqual(record["layer_sha256"], self.record["layer_sha256"])
        self.assertEqual(record["rebound_from"]["source_sha256"], self.source_sha)
        opened = cache_mod.open_cache(self.cache_path, new_sha, verify="metadata")
        self.assertEqual(opened.record["rebound_from"]["source_sha256"], self.source_sha)

    def test_rebind_refuses_a_different_matrix(self):
        import h5py
        import shutil

        copy_dir = Path(self.tmp.name) / "tampered"
        copy_dir.mkdir()
        shutil.copy(self.source, copy_dir / "tahoe_transfer.h5ad")
        with h5py.File(copy_dir / "tahoe_transfer.h5ad", "r+") as handle:
            handle["layers/log1p_cp10k/data"][0] += 1.0
        with self.assertRaisesRegex(cache_mod.CacheError, "does not match the frozen layer"):
            cache_mod.rebind_cache(
                self.cache_path, copy_dir / "tahoe_transfer.h5ad", "c" * 64
            )

    def test_dataset_uses_cache_transparently(self):
        frozen = data_mod.open_frozen_input(self.dir)
        opened = cache_mod.open_cache(self.cache_path, self.source_sha, self.source)
        with_cache = data_mod.TahoeDataset(
            frozen=frozen, structures=frozen.structures(), max_train_rows=8,
            max_valid_rows=4, seed=3, cache=opened,
        ).prepare()
        without_cache = data_mod.TahoeDataset(
            frozen=frozen, structures=frozen.structures(), max_train_rows=8,
            max_valid_rows=4, seed=3,
        ).prepare()
        np.testing.assert_allclose(
            with_cache.target(with_cache.rows_train), without_cache.target(without_cache.rows_train)
        )
        self.assertEqual(with_cache.cache_record["mode"], "cache")


if __name__ == "__main__":
    unittest.main()
