"""Resume-safe queue: completed units skipped, tampered config re-run, status CSV."""

from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

import _support
from ptlib import queue as queue_mod


def _unit(unit_id: str, seed: int, stage: str = "pretrain", cache=None) -> dict:
    args = [
        "--seed", str(seed),
        "--device", "cpu",
        "--epochs", "1",
        "--min-epochs", "1",
        "--batch-size", "8",
        "--max-train-rows", "16",
        "--max-valid-rows", "8",
        "--no-model-card",
    ]
    if cache is not None:
        # Hermetic: keep the unit off the developer's real derived cache under
        # output/cache/ (bound to the frozen input hash, not to this fixture).
        args += ["--tahoe-cache", str(cache)]
    return {"unit_id": unit_id, "stage": stage, "args": args}


class QueueValidationTests(unittest.TestCase):
    def _write(self, payload) -> Path:
        tmp = tempfile.NamedTemporaryFile("w", suffix=".json", delete=False)
        json.dump(payload, tmp)
        tmp.close()
        self.addCleanup(lambda: Path(tmp.name).unlink(missing_ok=True))
        return Path(tmp.name)

    def test_duplicate_unit_id_rejected(self):
        path = self._write({"units": [_unit("a", 1), _unit("a", 2)]})
        with self.assertRaises(queue_mod.QueueError):
            queue_mod.load_queue(path)

    def test_unsafe_unit_id_rejected(self):
        path = self._write({"units": [_unit("../escape", 1)]})
        with self.assertRaises(queue_mod.QueueError):
            queue_mod.load_queue(path)

    def test_queue_owns_out_and_data_dir(self):
        path = self._write({"units": [{"unit_id": "a", "stage": "pretrain",
                                       "args": ["--out", "/tmp/x"]}]})
        with self.assertRaises(queue_mod.QueueError):
            queue_mod.load_queue(path)

    def test_finetune_requires_fold(self):
        path = self._write({"units": [{"unit_id": "a", "stage": "finetune", "args": []}]})
        with self.assertRaises(queue_mod.QueueError):
            queue_mod.load_queue(path)

    def test_dict_args_are_converted(self):
        path = self._write({"units": [{"unit_id": "a", "stage": "pretrain",
                                       "args": {"seed": 5, "no_model_card": True,
                                                "max_train_rows": 4}}]})
        units = queue_mod.load_queue(path)
        self.assertIn("--seed", units[0].args)
        self.assertIn("--no-model-card", units[0].args)
        self.assertIn("--max-train-rows", units[0].args)


class QueueResumeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = _support.small_fixture_dir()
        cls.tmp = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def _options(self, queue_path: Path, out_root: Path, **kwargs) -> queue_mod.QueueOptions:
        return queue_mod.QueueOptions(
            queue=queue_path,
            out_root=out_root,
            data_dir=self.dir,
            workers=kwargs.pop("workers", 2),
            force=kwargs.pop("force", False),
            dry_run=kwargs.pop("dry_run", False),
            omp_threads=1,
            python=_support.sys.executable,
            train_script=_support.EXECUTION / "train.py",
            **kwargs,
        )

    def test_resume_skips_completed_and_reruns_changed_config(self):
        queue_path = Path(self.tmp.name) / "queue.json"
        out_root = Path(self.tmp.name) / "runs"
        cache = Path(self.tmp.name) / "cache" / "tahoe_log1p_cp10k_v1.npy"
        queue_path.write_text(json.dumps({"units": [_unit("u1", 11, cache=cache),
                                                    _unit("u2", 12, cache=cache)]}))

        first = queue_mod.run_queue(self._options(queue_path, out_root))
        self.assertEqual(first["completed"], 2, first)
        self.assertEqual(first["failed"], 0, first)
        for unit_id in ("u1", "u2"):
            manifest = json.loads((out_root / unit_id / "manifest.json").read_text())
            self.assertEqual(manifest["status"], "completed")
            self.assertTrue((out_root / unit_id / "checkpoint" / "model.pt").exists())
            # both stages save/load the identical payload keys
            import torch

            payload = torch.load(out_root / unit_id / "checkpoint" / "model.pt",
                                 map_location="cpu", weights_only=False)
            self.assertEqual(set(payload), {
                "model", "state_dict", "config", "seed", "n_genes",
                "best_valid_auprc", "epoch",
            })

        second = queue_mod.run_queue(self._options(queue_path, out_root))
        self.assertEqual(second["skipped_complete"], 2, second)
        self.assertEqual(second["completed"], 0, second)
        self.assertEqual(second["failed"], 0, second)

        with (out_root / "queue_status.csv").open() as handle:
            status_rows = list(csv.DictReader(handle))
        self.assertEqual({row["unit_id"] for row in status_rows}, {"u1", "u2"})
        self.assertEqual({row["status"] for row in status_rows}, {"skipped_complete"})
        for row in status_rows:
            self.assertEqual(row["exit_code"], "0")
            self.assertNotEqual(row["wall_s"], "")

        # change one unit's config → only that unit re-runs
        queue_path.write_text(json.dumps({"units": [_unit("u1", 99, cache=cache),
                                                    _unit("u2", 12, cache=cache)]}))
        third = queue_mod.run_queue(self._options(queue_path, out_root))
        self.assertEqual(third["completed"], 1, third)
        self.assertEqual(third["skipped_complete"], 1, third)
        manifest = json.loads((out_root / "u1" / "manifest.json").read_text())
        self.assertEqual(manifest["seed"], 99)

    def test_dry_run_does_not_launch(self):
        queue_path = Path(self.tmp.name) / "dry.json"
        out_root = Path(self.tmp.name) / "dry-runs"
        queue_path.write_text(json.dumps({"units": [_unit("d1", 5)]}))
        summary = queue_mod.run_queue(
            self._options(queue_path, out_root, dry_run=True)
        )
        self.assertEqual(summary["planned"], 1)
        self.assertFalse((out_root / "d1").exists())


if __name__ == "__main__":
    unittest.main()
