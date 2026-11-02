"""Parquet/CSV table cache: the pinned env has no parquet engine, so fail closed."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

import _support
from ptlib import tables


class TableCacheTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = _support.small_fixture_dir()
        cls.sha = _support.artifact_hashes(cls.dir)
        cls.tmp = tempfile.TemporaryDirectory()
        cls.cache = Path(cls.tmp.name) / "input_cache"
        for name in ("folds.parquet", "consensus_membership.parquet",
                     "tahoe_control_groups.parquet"):
            tables.convert_parquet(cls.dir / name, cls.cache, cls.sha[name])

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_cache_manifest_records_source_hashes(self):
        payload = tables.read_table.__globals__  # noqa: F841  (documentation)
        import json

        manifest = json.loads((self.cache / "manifest.json").read_text())
        self.assertEqual(manifest["schema"], tables.CACHE_SCHEMA)
        self.assertEqual(manifest["tables"]["folds.parquet"]["source_sha256"],
                         self.sha["folds.parquet"])

    def test_fallback_returns_the_same_values(self):
        with mock.patch.object(tables, "parquet_engine_available", return_value=(False, "none")):
            via_cache = tables.read_table(
                self.dir / "folds.parquet", expected_sha256=self.sha["folds.parquet"],
                cache_dir=self.cache,
            )
        direct = pd.read_parquet(self.dir / "folds.parquet")
        pd.testing.assert_frame_equal(
            via_cache.reset_index(drop=True), direct.reset_index(drop=True)
        )

    def test_stale_cache_is_rejected(self):
        with mock.patch.object(tables, "parquet_engine_available", return_value=(False, "none")):
            with self.assertRaisesRegex(tables.TableCacheError, "different source"):
                tables.read_table(
                    self.dir / "folds.parquet", expected_sha256="0" * 64,
                    cache_dir=self.cache,
                )

    def test_missing_cache_fails_with_an_actionable_message(self):
        with mock.patch.object(tables, "parquet_engine_available", return_value=(False, "none")):
            with self.assertRaisesRegex(tables.TableCacheError, "convert_frozen_tables"):
                tables.read_table(
                    self.dir / "folds.parquet", expected_sha256=self.sha["folds.parquet"],
                    cache_dir=Path(self.tmp.name) / "absent-cache",
                )


class SmokeGateTests(unittest.TestCase):
    """`smoke_real.py` must refuse to run without ticket #3's certified frozen input."""

    def test_refuses_missing_manifest(self):
        import smoke_real

        with tempfile.TemporaryDirectory() as tmp:
            empty = Path(tmp) / "not-frozen"
            empty.mkdir()
            code = smoke_real.main(
                ["--data-dir", str(empty), "--out-root", str(Path(tmp) / "out")]
            )
        self.assertEqual(code, 2)

    def test_refuses_tampered_artifacts(self):
        import smoke_real

        data_dir = _support.copy_fixture("smoke-gate-tampered")
        (data_dir / "folds.parquet").write_bytes((data_dir / "folds.parquet").read_bytes() + b"x")
        code = smoke_real.main(
            ["--data-dir", str(data_dir), "--out-root", str(data_dir.parent / "gate-out")]
        )
        self.assertEqual(code, 2)


if __name__ == "__main__":
    unittest.main()
