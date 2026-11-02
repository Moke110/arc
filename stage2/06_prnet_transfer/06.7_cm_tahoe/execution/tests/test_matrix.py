"""Unit tests for the ticket-#7 matrix definition and the per-fold LR selection."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matrix as matrix_mod  # noqa: E402

LOCKED = {
    "n_folds": 5,
    "seeds": [20261021, 20261022, 20261023],
    "lrs": {"cm": [1e-4, 3e-4, 1e-3], "pretrain": 1e-3},
}


class LrLabelTests(unittest.TestCase):
    def test_flat_labels(self):
        self.assertEqual(matrix_mod.lr_label(1e-4), "1e-4")
        self.assertEqual(matrix_mod.lr_label(3e-4), "3e-4")
        self.assertEqual(matrix_mod.lr_label(1e-3), "1e-3")

    def test_rejects_non_positive(self):
        for bad in (0, -1e-4, float("nan"), float("inf")):
            with self.assertRaises(matrix_mod.MatrixError):
                matrix_mod.lr_label(bad)


class UnitIdTests(unittest.TestCase):
    def test_roundtrip(self):
        unit_id = matrix_mod.unit_id_for(3, 20261022, 3e-4)
        self.assertEqual(unit_id, "cm_f3_s20261022_lr3e-4_tahoe")
        parsed = matrix_mod.parse_unit_id(unit_id)
        self.assertEqual(parsed["fold"], 3)
        self.assertEqual(parsed["seed"], 20261022)
        self.assertEqual(parsed["lr_label"], "3e-4")

    def test_rejects_foreign_ids(self):
        for bad in ("fold0_seed1", "cm_f0_s1_lr1e-4_random", "cm_f0_s1_tahoe"):
            with self.assertRaises(matrix_mod.MatrixError):
                matrix_mod.parse_unit_id(bad)


class MatrixSpecTests(unittest.TestCase):
    def test_matrix_is_45_and_ordered(self):
        spec = matrix_mod.matrix_spec(LOCKED)
        self.assertEqual(spec["n_candidates"], 45)
        self.assertEqual(spec["init"], "tahoe")
        self.assertEqual(spec["pipeline_init"], "pretrained")
        cands = matrix_mod.candidates(spec)
        self.assertEqual(len(cands), 45)
        self.assertEqual(len({c.unit_id for c in cands}), 45)
        # fold-major, then ascending LR, then seed
        self.assertEqual([c.fold for c in cands[:9]], [0] * 9)
        self.assertEqual([c.lr for c in cands[:3]], [1e-4] * 3)
        self.assertEqual([c.seed for c in cands[:3]], LOCKED["seeds"])

    def test_train_args_carry_the_warm_start(self):
        spec = matrix_mod.matrix_spec(LOCKED)
        candidate = matrix_mod.candidates(spec)[4]
        args = matrix_mod.train_args(
            candidate, "/checkpoints/tahoe_seed20261022/checkpoint/model.pt",
            hash_cache="/tmp/hash.json",
        )
        self.assertEqual(args, [
            "--fold", str(candidate.fold),
            "--seed", str(candidate.seed),
            "--lr", repr(float(candidate.lr)),
            "--init", "pretrained",
            "--init-checkpoint", "/checkpoints/tahoe_seed20261022/checkpoint/model.pt",
            "--device", "cuda",
            "--hash-cache", "/tmp/hash.json",
        ])
        # no --out/--data-dir/--stage: those belong to the queue
        for owned in ("--out", "--data-dir", "--stage"):
            self.assertNotIn(owned, args)

    def test_seed_checkpoint_mapping(self):
        spec = matrix_mod.matrix_spec(LOCKED)
        with tempfile.TemporaryDirectory() as tmp:
            for seed in LOCKED["seeds"]:
                path = Path(tmp) / f"tahoe_seed{seed}" / "checkpoint" / "model.pt"
                path.parent.mkdir(parents=True)
                path.write_bytes(b"fixture")
            mapping = matrix_mod.seed_checkpoint_map(spec["seeds"], tmp)
            self.assertEqual(sorted(mapping), LOCKED["seeds"])
            for seed in LOCKED["seeds"]:
                self.assertEqual(
                    mapping[seed],
                    Path(tmp) / f"tahoe_seed{seed}" / "checkpoint" / "model.pt",
                )

    def test_seed_checkpoint_map_fails_closed(self):
        spec = matrix_mod.matrix_spec(LOCKED)
        with self.assertRaisesRegex(matrix_mod.MatrixError, "no ticket-#5 formal checkpoint"):
            matrix_mod.seed_checkpoint_map(spec["seeds"], "/nonexistent/runs")

    def test_rejects_drifting_configs(self):
        bad = []
        bad.append({**LOCKED, "seeds": [1, 2]})
        bad.append({**LOCKED, "lrs": {"cm": [1e-4, 3e-4]}})
        bad.append({**LOCKED, "lrs": {"cm": [3e-4, 1e-4, 1e-3]}})
        bad.append({**LOCKED, "seeds": [1, 1, 2]})
        bad.append({**LOCKED, "n_folds": 0})
        for payload in bad:
            with self.assertRaises(matrix_mod.MatrixError):
                matrix_mod.matrix_spec(payload)


class DeriveStatusTests(unittest.TestCase):
    def _manifest(self, **overrides):
        manifest = {
            "status": "completed",
            "checkpoint": {"reload_identical": True},
            "best_valid_auprc": 0.1,
            "epochs_run": 3,
            "stop_reason": "early_stop_no_auprc_improvement_20",
        }
        manifest.update(overrides)
        return manifest

    def test_ok(self):
        status, _ = matrix_mod.derive_status(self._manifest(), {"history": [{}, {}, {}]}, True, True)
        self.assertEqual(status, "ok")

    def test_capped(self):
        status, _ = matrix_mod.derive_status(
            self._manifest(stop_reason="max_epochs"), {"history": [{}, {}, {}]}, True, True
        )
        self.assertEqual(status, "capped")

    def test_failed_variants(self):
        cases = [
            (None, None, True, True),
            (self._manifest(status="failed"), None, True, True),
            (self._manifest(), None, False, True),
            (self._manifest(), None, True, False),
            (self._manifest(checkpoint={"reload_identical": False}), None, True, True),
            (self._manifest(best_valid_auprc=None), None, True, True),
            (self._manifest(stop_reason="no_finite_auprc"), None, True, True),
            (self._manifest(), {"history": [{}, {}]}, True, True),
        ]
        for manifest, history, exists, sha_ok in cases:
            status, detail = matrix_mod.derive_status(manifest, history, exists, sha_ok)
            self.assertEqual(status, "failed", detail)


class SelectLrTests(unittest.TestCase):
    def setUp(self):
        self.spec = matrix_mod.matrix_spec(LOCKED)
        self.cands = matrix_mod.candidates(self.spec)

    def _records(self, value):
        records = {}
        for candidate in self.cands:
            auprc = value(candidate) if callable(value) else value
            records[candidate.unit_id] = {"status": "ok", "best_valid_auprc": auprc}
        return records

    def test_picks_highest_three_seed_mean(self):
        def auprc(candidate):
            return {1e-4: 0.10, 3e-4: 0.12, 1e-3: 0.11}[candidate.lr]

        folds = matrix_mod.select_lr(self.cands, self._records(auprc))
        self.assertEqual(len(folds), 5)
        for fold in folds:
            self.assertTrue(fold["complete"])
            self.assertEqual(fold["selected_lr"], 3e-4)
            self.assertFalse(fold["tie_break_applied"])
            self.assertAlmostEqual(fold["selected_mean_auprc"], 0.12, places=12)

    def test_tie_break_takes_lower_lr(self):
        def auprc(candidate):
            return {1e-4: 0.12, 3e-4: 0.12, 1e-3: 0.05}[candidate.lr]

        folds = matrix_mod.select_lr(self.cands, self._records(auprc))
        for fold in folds:
            self.assertEqual(fold["selected_lr"], 1e-4)
            self.assertTrue(fold["tie_break_applied"])

    def test_mean_uses_all_three_seeds(self):
        def auprc(candidate):
            base = 0.10 if candidate.lr == 1e-3 else 0.09
            return base + 0.001 * candidate.seed_index

        folds = matrix_mod.select_lr(self.cands, self._records(auprc))
        fold = folds[0]
        expected = (0.10 + 0.0 + 0.10 + 0.001 + 0.10 + 0.002) / 3
        self.assertAlmostEqual(fold["selected_mean_auprc"], expected, places=15)

    def test_incomplete_lr_is_not_eligible(self):
        records = self._records(0.1)
        victim = next(c for c in self.cands if c.fold == 0 and c.lr == 1e-4 and c.seed_index == 1)
        records[victim.unit_id] = {"status": "failed", "best_valid_auprc": None}
        folds = {row["fold"]: row for row in matrix_mod.select_lr(self.cands, records)}
        self.assertEqual(folds[0]["selected_lr"], 3e-4)
        scores = {score["lr"]: score for score in folds[0]["lr_scores"]}
        self.assertEqual(scores[1e-4]["n_seeds_ok"], 2)
        self.assertFalse(scores[1e-4]["complete"])
        self.assertIsNone(scores[1e-4]["mean_auprc"])

    def test_no_complete_lr_marks_fold_incomplete(self):
        records = self._records(0.1)
        for candidate in [c for c in self.cands if c.fold == 2 and c.seed_index == 0]:
            records[candidate.unit_id] = {"status": "failed", "best_valid_auprc": None}
        folds = {row["fold"]: row for row in matrix_mod.select_lr(self.cands, records)}
        self.assertFalse(folds[2]["complete"])
        self.assertIsNone(folds[2]["selected_lr"])
        self.assertTrue(folds[0]["complete"])

    def test_capped_counts_as_complete(self):
        records = self._records(0.1)
        victim = next(c for c in self.cands if c.fold == 0)
        records[victim.unit_id] = {"status": "capped", "best_valid_auprc": 0.1}
        folds = {row["fold"]: row for row in matrix_mod.select_lr(self.cands, records)}
        self.assertTrue(folds[0]["complete"])


if __name__ == "__main__":
    unittest.main()
