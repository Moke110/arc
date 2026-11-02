"""Memory-bounded fit/predict paths (ticket #8's 9.7 GiB allocation fix).

The baselines must never materialise an ``n_observations × 18,164`` array.  These
tests pin the two structural mechanisms on the synthetic fixture:

* :func:`baselines.iter_layer_blocks` streams the dense CM log1p layer in gene
  blocks and reconstructs exactly the full row read;
* :func:`baselines.plate_control_matrix` reproduces
  ``CMControlPool.mean_for_rows`` from ``n_plates × n_genes`` plate means, so no
  ``n_rows × n_genes`` control matrix is built;
* :func:`baselines.build_fit_table` aggregates the row Δ condition by condition
  from those blocks and never calls the full-row layer reader;
* :func:`baselines.condition_predictions` + the blocked
  :func:`baselines.write_prediction_dir` produce byte-identical prediction files
  to writing a dense row-level Δ (the evaluator-sanity path), so #9 scores the
  same numbers either way.

A rising peak here is a defect (``docs/agents/execution-environment.md``); the
matching end-to-end number is measured by ``make mem-check`` (row-replicated
fixture), not by this offline suite.
"""

from __future__ import annotations

import unittest
from pathlib import Path
from unittest import mock

import numpy as np

import _support
import baselines
from ptlib import data as data_mod
from ptlib import h5io
from ptlib import manifest as manifest_mod
from ptlib import metrics as metrics_mod


def _fold_zero():
    frozen = _support.open_frozen()
    folds = _support.folds_table(frozen)
    structures = frozen.structures()
    source = data_mod.CMPredictSource(frozen, folds, 0, structures)
    pool = data_mod.CMControlPool(frozen)
    return frozen, folds, structures, source, pool


class BlockedLayerReadTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        (
            cls.frozen,
            cls.folds,
            cls.structures,
            cls.source,
            cls.pool,
        ) = _fold_zero()
        cls.rows = np.asarray(cls.source.rows_for("train"), dtype=np.int64)
        cls.n_genes = int(cls.frozen.n_genes)

    def test_blocks_reconstruct_the_full_row_read(self):
        reference, sorted_rows = h5io.read_layer_rows(
            self.frozen.cm_h5ad, self.rows, self.n_genes
        )
        np.testing.assert_array_equal(sorted_rows, self.rows)
        rebuilt = np.empty((self.rows.size, self.n_genes), dtype=np.float32)
        block_sizes = []
        for g0, g1, block in baselines.iter_layer_blocks(
            self.frozen.cm_h5ad, self.rows, self.n_genes
        ):
            self.assertEqual(block.dtype, np.float32)
            self.assertEqual(block.shape, (self.rows.size, g1 - g0))
            self.assertLessEqual(block.shape[1], baselines.GENE_BLOCK)
            block_sizes.append(block.shape[1])
            rebuilt[:, g0:g1] = block
        self.assertEqual(sum(block_sizes), self.n_genes)
        np.testing.assert_array_equal(rebuilt, reference)

    def test_blocks_are_bounded_in_row_count_times_gene_block(self):
        # The bounded temporary is n_rows × GENE_BLOCK (float32/float64), not
        # n_rows × n_genes: assert the block contract itself.
        for _, g1, block in baselines.iter_layer_blocks(
            self.frozen.cm_h5ad, self.rows, self.n_genes
        ):
            self.assertLessEqual(block.size, self.rows.size * baselines.GENE_BLOCK)
            self.assertLess(block.size, self.rows.size * self.n_genes)

    def test_unsorted_rows_are_refused(self):
        bad = self.rows[::-1]
        with self.assertRaises(baselines.BaselineError):
            next(iter(baselines.iter_layer_blocks(self.frozen.cm_h5ad, bad, self.n_genes)))


class CompactControlTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.frozen, cls.folds, cls.structures, cls.source, cls.pool = _fold_zero()
        cls.rows = np.asarray(cls.source.rows_for("train"), dtype=np.int64)

    def test_plate_matrix_matches_mean_for_rows(self):
        matrix, per_row = baselines.plate_control_matrix(self.pool, self.rows)
        self.assertEqual(matrix.shape, (len(self.pool.plate_means), self.frozen.n_genes))
        self.assertLess(matrix.size, self.rows.size * self.frozen.n_genes)
        reference = self.pool.mean_for_rows(self.rows)
        np.testing.assert_array_equal(matrix[per_row], reference)


class StreamingFitTableTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.frozen, cls.folds, cls.structures, cls.source, cls.pool = _fold_zero()
        cls.rows = np.asarray(cls.source.rows_for("train"), dtype=np.int64)
        cls.table = baselines.build_fit_table(
            cls.frozen,
            cls.folds,
            0,
            structures=cls.structures,
            control_pool=cls.pool,
            source=cls.source,
        )

    def _reference_condition_delta(self) -> np.ndarray:
        layer, sorted_rows = h5io.read_layer_rows(
            self.frozen.cm_h5ad, self.rows, self.frozen.n_genes
        )
        np.testing.assert_array_equal(sorted_rows, self.rows)
        row_delta = layer.astype(np.float64) - self.pool.mean_for_rows(
            self.rows
        ).astype(np.float64)
        return metrics_mod.aggregate_rows(
            row_delta, self.table.row_conditions, self.table.n_conditions
        )

    def test_streamed_condition_means_match_the_reference(self):
        reference = self._reference_condition_delta()
        np.testing.assert_allclose(
            self.table.condition_delta.astype(np.float64), reference, rtol=1e-5, atol=1e-6
        )
        np.testing.assert_allclose(
            self.table.pooled_mean, reference.mean(axis=0), rtol=1e-5, atol=1e-6
        )

    def test_fit_table_does_not_call_the_full_row_layer_reader(self):
        def _refuse(*args, **kwargs):  # pragma: no cover - fired only on regression
            raise AssertionError(
                "build_fit_table must read the layer in gene blocks, not full rows"
            )

        with mock.patch.object(baselines.h5io, "read_layer_rows", _refuse):
            table = baselines.build_fit_table(
                self.frozen,
                self.folds,
                0,
                structures=self.structures,
                control_pool=self.pool,
                source=self.source,
            )
        np.testing.assert_array_equal(table.condition_delta, self.table.condition_delta)
        self.assertEqual(table.train_delta_sha256, self.table.train_delta_sha256)

    def test_condition_delta_sha256_is_block_size_independent(self):
        # The record's train_delta_sha256 hashes the condition-level Δ, so it
        # does not change when GENE_BLOCK is tuned.
        original = baselines.GENE_BLOCK
        try:
            baselines.GENE_BLOCK = 2048
            table = baselines.build_fit_table(
                self.frozen,
                self.folds,
                0,
                structures=self.structures,
                control_pool=self.pool,
                source=self.source,
            )
        finally:
            baselines.GENE_BLOCK = original
        self.assertEqual(table.train_delta_sha256, self.table.train_delta_sha256)


class BlockedPredictionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.frozen, cls.folds, cls.structures, cls.source, cls.pool = _fold_zero()
        cls.table = baselines.build_fit_table(
            cls.frozen,
            cls.folds,
            0,
            structures=cls.structures,
            control_pool=cls.pool,
            source=cls.source,
        )
        cls.rows = np.asarray(cls.source.rows_for("test"), dtype=np.int64)
        cls.meta = cls.source.meta(cls.rows)
        cls.c_star = np.asarray(cls.source.c_star, dtype=np.float32)
        cls.context = baselines.condition_context(
            cls.frozen, cls.source, cls.structures, cls.rows
        )
        cls.tmp = _support.tmp_root() / "memory"
        cls.tmp.mkdir(parents=True, exist_ok=True)

    def test_condition_predictions_are_the_row_wrapper_source(self):
        for name in baselines.BASELINES:
            params, _ = baselines.fit_baseline(
                name, self.table, self.structures, self.frozen
            )
            condition_delta, info = baselines.condition_predictions(
                name, params, self.table, self.context, self.frozen
            )
            self.assertEqual(
                condition_delta.shape,
                (self.context["condition_keys"].size, self.frozen.n_genes),
            )
            self.assertEqual(condition_delta.dtype, np.float32)
            self.assertLess(
                condition_delta.size,
                self.rows.size * self.frozen.n_genes,
                f"{name}: condition matrix should be smaller than the row matrix",
            )
            row_level = baselines.predict_baseline(
                name, params, self.table, self.context, self.frozen
            )
            np.testing.assert_array_equal(
                row_level.delta,
                condition_delta[np.asarray(self.context["row_conditions"], dtype=np.int64)],
            )
            self.assertIsInstance(info, dict)

    def _write(self, name: str, **kwargs) -> Path:
        out_dir = self.tmp / name / "pred"
        baselines.write_prediction_dir(
            out_dir,
            baseline=name,
            role="primary",
            fold=0,
            meta=self.meta,
            c_star=self.c_star,
            fit_record={},
            frozen=self.frozen,
            **kwargs,
        )
        return out_dir

    def test_blocked_writer_equals_a_dense_row_level_write(self):
        params, _ = baselines.fit_baseline(
            "delta_mean", self.table, self.structures, self.frozen
        )
        condition_delta, _ = baselines.condition_predictions(
            "delta_mean", params, self.table, self.context, self.frozen
        )
        inverse = np.asarray(self.context["row_conditions"], dtype=np.int64)
        dense = condition_delta[inverse]

        blocked_dir = self._write(
            "blocked", condition_delta=condition_delta, row_conditions=inverse
        )
        dense_dir = self._write("dense", delta=dense)
        for relative in ("pred_expr.npy", "delta.npy", "meta.csv"):
            self.assertEqual(
                manifest_mod.sha256_file(blocked_dir / relative),
                manifest_mod.sha256_file(dense_dir / relative),
                f"{relative} differs between the blocked and dense writers",
            )
        # and the two files are the pipeline contract: pred_expr − delta == c*
        pred_expr = np.load(blocked_dir / "pred_expr.npy")
        stored = np.load(blocked_dir / "delta.npy")
        self.assertEqual(pred_expr.shape, (self.rows.size, self.frozen.n_genes))
        self.assertEqual(stored.shape, pred_expr.shape)
        self.assertLess(
            float(np.max(np.abs(pred_expr.astype(np.float64) - stored.astype(np.float64)
                               - self.c_star.astype(np.float64)))),
            1e-5,
        )

    def test_writer_records_file_hashes_and_block_config(self):
        params, _ = baselines.fit_baseline(
            "no_change", self.table, self.structures, self.frozen
        )
        condition_delta, _ = baselines.condition_predictions(
            "no_change", params, self.table, self.context, self.frozen
        )
        pred_dir = self._write(
            "manifest",
            condition_delta=condition_delta,
            row_conditions=np.asarray(self.context["row_conditions"], dtype=np.int64),
        )
        payload = manifest_mod.read_json(pred_dir / "manifest.json")
        self.assertEqual(
            payload["predictions"]["pred_expr_sha256"],
            manifest_mod.sha256_file(pred_dir / "pred_expr.npy"),
        )
        self.assertEqual(
            payload["predictions"]["delta_sha256"],
            manifest_mod.sha256_file(pred_dir / "delta.npy"),
        )
        memory = payload["config"]["memory"]
        self.assertEqual(memory["row_block_rows"], baselines.ROW_BLOCK)
        self.assertEqual(memory["gene_block_genes"], baselines.GENE_BLOCK)
        self.assertIn("memmap", memory["prediction_arrays"])

    def test_writer_refuses_ambiguous_input(self):
        dense = np.zeros((self.rows.size, self.frozen.n_genes), dtype=np.float32)
        with self.assertRaises(baselines.BaselineError):
            baselines.write_prediction_dir(
                self.tmp / "ambiguous" / "pred",
                baseline="no_change",
                role="primary",
                fold=0,
                meta=self.meta,
                c_star=self.c_star,
                fit_record={},
                frozen=self.frozen,
                delta=dense,
                condition_delta=dense,
                row_conditions=np.zeros(self.rows.size, dtype=np.int64),
            )


class BlockConstantsTests(unittest.TestCase):
    def test_block_constants_are_small_and_positive(self):
        self.assertIsInstance(baselines.GENE_BLOCK, int)
        self.assertIsInstance(baselines.ROW_BLOCK, int)
        self.assertGreater(baselines.GENE_BLOCK, 0)
        self.assertGreater(baselines.ROW_BLOCK, 0)
        self.assertLess(baselines.GENE_BLOCK, 18164)
        self.assertLess(baselines.ROW_BLOCK, 4096)


if __name__ == "__main__":
    unittest.main()
