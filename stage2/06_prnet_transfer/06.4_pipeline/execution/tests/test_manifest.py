"""Run-manifest schema, frozen-artifact verification (fail closed) and resume checks."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import _support
from ptlib import manifest


class HashVerificationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        (self.dir / "a.bin").write_bytes(b"alpha")
        (self.dir / "b.bin").write_bytes(b"beta")
        self.hashes = {
            "a.bin": manifest.sha256_file(self.dir / "a.bin"),
            "b.bin": manifest.sha256_file(self.dir / "b.bin"),
        }

    def tearDown(self):
        self.tmp.cleanup()

    def test_verify_passes_and_reports_bytes(self):
        record = manifest.verify_artifacts(
            self.dir, self.hashes, required=("a.bin", "b.bin")
        )
        self.assertEqual(set(record["verified"]), {"a.bin", "b.bin"})
        self.assertEqual(record["bytes"], 9)
        self.assertEqual(record["cached"], [])

    def test_tampered_file_fails_closed(self):
        (self.dir / "b.bin").write_bytes(b"BETA")
        with self.assertRaisesRegex(manifest.ManifestError, "sha256 mismatch"):
            manifest.verify_artifacts(self.dir, self.hashes, required=("a.bin", "b.bin"))

    def test_missing_entry_fails_closed(self):
        with self.assertRaisesRegex(manifest.ManifestError, "not listed"):
            manifest.verify_artifacts(
                self.dir, {"a.bin": self.hashes["a.bin"]}, required=("a.bin", "b.bin")
            )

    def test_hash_cache_is_used_and_recorded(self):
        cache = self.dir / "cache.json"
        first = manifest.verify_artifacts(
            self.dir, self.hashes, required=("a.bin", "b.bin"), hash_cache=cache
        )
        self.assertEqual(first["cached"], [])
        second = manifest.verify_artifacts(
            self.dir, self.hashes, required=("a.bin", "b.bin"), hash_cache=cache
        )
        self.assertEqual(sorted(second["cached"]), ["a.bin", "b.bin"])

    def test_hash_cache_is_invalidated_by_a_changed_file(self):
        cache = self.dir / "cache.json"
        manifest.verify_artifacts(
            self.dir, self.hashes, required=("a.bin",), hash_cache=cache
        )
        (self.dir / "a.bin").write_bytes(b"alpha2")
        new_hashes = {"a.bin": manifest.sha256_file(self.dir / "a.bin")}
        record = manifest.verify_artifacts(
            self.dir, new_hashes, required=("a.bin",), hash_cache=cache
        )
        self.assertEqual(record["cached"], [])


class NormalizeManifestTests(unittest.TestCase):
    def test_flat_mapping(self):
        digest = "a" * 64
        self.assertEqual(manifest.normalize_artifact_manifest({"x.h5ad": digest}),
                         {"x.h5ad": digest})

    def test_artifacts_mapping_with_records(self):
        digest = "b" * 64
        out = manifest.normalize_artifact_manifest(
            {"artifacts": {"x.h5ad": {"sha256": digest, "bytes": 3}}}
        )
        self.assertEqual(out, {"x.h5ad": digest})

    def test_files_list(self):
        digest = "c" * 64
        out = manifest.normalize_artifact_manifest(
            {"files": [{"path": "dir/x.h5ad", "sha256": digest}]}
        )
        self.assertEqual(out, {"x.h5ad": digest})

    def test_missing_hashes_rejected(self):
        with self.assertRaises(manifest.ManifestError):
            manifest.normalize_artifact_manifest({"artifacts": {"x.h5ad": {"md5": "x"}}})


class SchemaTests(unittest.TestCase):
    def _base(self, kind: str = "train") -> dict:
        payload = {
            "schema": manifest.SCHEMA,
            "created_utc": manifest.now_utc(),
            "kind": kind,
            "stage": "pretrain",
            "status": "completed",
            "pipeline_version": "06.4.0",
            "command": ["train.py"],
            "git": {"sha": "deadbeef", "branch": "b", "dirty": False},
            "data": {"dir": "/data", "artifacts": {}},
            "locked_config": {"path": "/data/locked_config.json", "sha256": "x"},
            "config": {"structure": {}},
            "config_sha256": "x",
            "structure": {"n_genes": 18164},
            "seed": 1,
            "device": {"resolved": "cpu"},
            "init": {"mode": "random"},
            "checkpoint": {"path": "checkpoint/model.pt", "sha256": "y"},
            "noise_policy": {"train": "seeded"},
            "timings": {"wall_s": 1.0},
            "peak_rss_mb": 10.0,
            "epochs_run": 1,
            "stop_reason": "max_epochs",
        }
        return payload

    def test_valid_train_manifest(self):
        manifest.validate_manifest(self._base())

    def test_missing_required_field_rejected(self):
        payload = self._base()
        del payload["stop_reason"]
        with self.assertRaisesRegex(manifest.ManifestError, "stop_reason"):
            manifest.validate_manifest(payload)

    def test_unknown_kind_rejected(self):
        with self.assertRaises(manifest.ManifestError):
            manifest.build_manifest("deploy", status="completed")

    def test_bad_status_rejected(self):
        payload = self._base()
        payload["status"] = "running"
        with self.assertRaises(manifest.ManifestError):
            manifest.validate_manifest(payload)

    def test_predict_manifest_requires_rows(self):
        payload = self._base("predict")
        for key in ("structure", "seed", "init", "epochs_run", "stop_reason", "checkpoint"):
            payload.pop(key, None)
        payload.update(data={"dir": "/data", "artifacts": {}}, config_sha256="x",
                       checkpoint={"path": "m.pt", "sha256": "y"}, rows=3)
        manifest.validate_manifest(payload)
        del payload["rows"]
        with self.assertRaises(manifest.ManifestError):
            manifest.validate_manifest(payload)

    def test_write_read_round_trip_with_numpy(self):
        import numpy as np

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "manifest.json"
            manifest.write_json(path, {"value": np.float32(0.5), "n": np.int64(3)})
            loaded = manifest.read_json(path)
        self.assertEqual(loaded["n"], 3)
        self.assertAlmostEqual(loaded["value"], 0.5)


class ResumeTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def _complete_run(self) -> None:
        checkpoint_dir = self.dir / "checkpoint"
        checkpoint_dir.mkdir()
        (checkpoint_dir / "model.pt").write_bytes(b"weights")
        payload = SchemaTests()._base()
        payload.update(
            created_utc=manifest.now_utc(),
            schema=manifest.SCHEMA,
            checkpoint={
                "path": "checkpoint/model.pt",
                "sha256": manifest.sha256_file(checkpoint_dir / "model.pt"),
            },
        )
        manifest.write_json(self.dir / "manifest.json", payload)

    def test_complete_unit_detected(self):
        self._complete_run()
        complete, reason = manifest.unit_complete(self.dir)
        self.assertTrue(complete, reason)

    def test_config_fingerprint_mismatch_forces_rerun(self):
        self._complete_run()
        complete, reason = manifest.unit_complete(self.dir, expected_config_sha256="other")
        self.assertFalse(complete)
        self.assertIn("config", reason)

    def test_tampered_checkpoint_forces_rerun(self):
        self._complete_run()
        (self.dir / "checkpoint" / "model.pt").write_bytes(b"tampered")
        complete, reason = manifest.unit_complete(self.dir)
        self.assertFalse(complete)
        self.assertIn("sha256", reason)

    def test_missing_manifest_is_incomplete(self):
        complete, _ = manifest.unit_complete(self.dir)
        self.assertFalse(complete)


class ConcurrentHashCacheTests(unittest.TestCase):
    """The hash memo is shared by every worker of a parallel training queue.

    Ticket #7 lost one candidate to this: two simultaneous units both wrote a
    fixed ``hash_cache.json.tmp``, the first ``Path.replace`` moved it away and the
    second raised ``FileNotFoundError``, failing that unit through no fault of its
    data.  Saves must therefore use unique temp files.
    """

    def test_concurrent_saves_all_succeed_without_leftovers(self):
        import concurrent.futures
        import json as _json

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "hash_cache.json"
            payloads = [
                {"entry": index, "sha256": f"{index:064d}"} for index in range(8)
            ]

            with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
                futures = [
                    pool.submit(manifest._save_hash_cache, path, payload)
                    for _ in range(20)
                    for payload in payloads
                ]
                for future in futures:
                    future.result()  # the old fixed-tmp code raised FileNotFoundError

            stored = _json.loads(path.read_text(encoding="utf-8"))
            self.assertIn(stored, payloads)
            leftovers = sorted(
                entry.name for entry in Path(tmp).iterdir() if entry.name != path.name
            )
            self.assertEqual(leftovers, [])


if __name__ == "__main__":
    unittest.main()
