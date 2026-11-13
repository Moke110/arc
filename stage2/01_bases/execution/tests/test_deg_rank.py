"""Disease-DEG table/rank contract: alignment, effect score, collapse, ranks.

These are the pure functions tickets #187/#188/#189 import; every case here is a
hand-worked fixture, not a re-derivation from the implementation's internals.
"""
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases import deg_rank


def axis_fixture():
    return pd.DataFrame({
        "gene_index": np.arange(8, dtype=np.int64),
        "symbol": [f"G{index}" for index in range(8)],
        "hgnc_id": [f"HGNC:{index}" for index in range(8)],
        "ensembl_id": [f"ENSG{index:011d}" for index in range(8)],
    })


def crosswalk_fixture():
    return pd.DataFrame([
        {"dataset_id": "x", "source_gene_id": "X1", "hgnc_id": "HGNC:1",
         "match_type": "alias"},
        {"dataset_id": "x", "source_gene_id": "AMB", "hgnc_id": "HGNC:2",
         "match_type": "alias"},
        {"dataset_id": "x", "source_gene_id": "AMB", "hgnc_id": "HGNC:3",
         "match_type": "alias"},
    ])


def hgnc_fixture():
    return pd.DataFrame([
        {"hgnc_id": "HGNC:4", "alias_symbol": "HALIAS", "prev_symbol": "",
         "ensembl_gene_id": "ENSG00000000004"},
        {"hgnc_id": "HGNC:5", "alias_symbol": "", "prev_symbol": "",
         "ensembl_gene_id": "ENSG00000000005"},
        {"hgnc_id": "HGNC:99", "alias_symbol": "OFFAXIS", "prev_symbol": "",
         "ensembl_gene_id": "ENSG00000000099"},
    ])


class DegEffectTest(unittest.TestCase):
    def test_effect_score_is_continuous_and_floor_bounded(self):
        effect = deg_rank.deg_effect(
            np.array([2.0, -2.0, 1.0, 3.0, np.nan]),
            np.array([1e-4, 1e-4, np.nan, 0.0, 1e-4]))
        np.testing.assert_allclose(effect[:2], [8.0, 8.0])
        # A missing padj is non-significant: exactly zero, not dropped.
        self.assertEqual(effect[2], 0.0)
        # padj=0 is floored at 1e-300, so the score stays finite and large.
        self.assertTrue(np.isfinite(effect[3]))
        np.testing.assert_allclose(effect[3], 900.0)
        # A non-finite log2FC scores zero.
        self.assertEqual(effect[4], 0.0)

    def test_na_padj_is_non_significant_not_a_missing_row(self):
        collapsed = deg_rank.per_gene_deg(
            np.array([0]), np.array([1.0]), np.array([np.nan]))
        self.assertEqual(len(collapsed), 1)
        self.assertEqual(float(collapsed.effect.iloc[0]), 0.0)
        self.assertFalse(bool(collapsed.significant.iloc[0]))
        self.assertEqual(int(collapsed.direction.iloc[0]), 0)

    def test_multi_contrast_collapses_by_maximum_effect(self):
        collapsed, _, _ = deg_rank.dataset_deg_rank(pd.DataFrame({
            "gene_index": [5, 5, 5, 6],
            "log2fc": [0.5, -3.0, 2.0, 1.0],
            "padj": [0.01, 1e-3, np.nan, 0.2],
            "contrast": ["a", "b", "c", "a"],
        }), axis_size=8)
        gene5 = collapsed[collapsed.gene_index.eq(5)].iloc[0]
        # |log2FC| * -log10(padj): 0.5*2=1, 3*3=9, 2*0=0 -> contrast b wins.
        self.assertEqual(gene5.contrast, "b")
        np.testing.assert_allclose(gene5.effect, 9.0)
        self.assertEqual(gene5.log2fc, -3.0)
        self.assertEqual(int(gene5.direction), -1)
        self.assertEqual(int(gene5.n_contrasts_supporting), 2)
        self.assertEqual(int(gene5.n_contrasts), 3)

    def test_rank_ties_break_on_gene_index(self):
        rank = deg_rank.rank_scores(
            np.array([4, 2, 9, 1]), np.array([1.0, 1.0, 1.0, 1.0]))
        self.assertEqual(rank.tolist(), [3, 2, 4, 1])

    def test_complete_rank_unmeasured_is_n_measured_plus_one(self):
        full = deg_rank.complete_rank(np.array([1, 3, 7]),
                                      np.array([1, 2, 3]), axis_size=8)
        self.assertEqual(full.tolist(), [4, 1, 4, 2, 4, 4, 4, 3])

    def test_rank_norm_denominator_is_n_measured_not_significant_count(self):
        # Three measured genes, one significant: rank_norm divides by 3, not 1.
        table, full, stats = deg_rank.dataset_deg_rank(
            pd.DataFrame({"gene_index": [0, 1, 2], "log2fc": [2.0, 1.0, 0.5],
                          "padj": [1e-4, 0.5, np.nan]}), axis_size=8)
        self.assertEqual(stats["n_measured"], 3)
        self.assertEqual(stats["n_significant"], 1)
        self.assertEqual(stats["rank_norm_denominator"], "n_measured")
        self.assertEqual(table["rank"].tolist(), [1, 2, 3])
        np.testing.assert_allclose(table["rank_norm"], [1 / 3, 2 / 3, 1.0])
        self.assertEqual(full["rank"].tolist()[:3], [1, 2, 3])
        self.assertEqual(full["rank"].iloc[3:].tolist(), [4, 4, 4, 4, 4])
        np.testing.assert_allclose(full["rank_norm"].iloc[3:], 1.0)

    def test_dataset_deg_rank_full_axis_and_determinism(self):
        rows = pd.DataFrame({
            "gene_index": [7, 7, 2, 0],
            "log2fc": [1.0, -4.0, 0.5, np.nan],
            "padj": [0.5, 1e-5, 0.2, 0.1],
        })
        first = deg_rank.dataset_deg_rank(rows, axis_size=8)
        second = deg_rank.dataset_deg_rank(rows, axis_size=8)
        pd.testing.assert_frame_equal(first[0], second[0])
        pd.testing.assert_frame_equal(first[1], second[1])
        self.assertEqual(first[2], second[2])
        table, full, stats = first
        self.assertEqual(table.gene_index.tolist(), [7, 2])
        self.assertEqual(full.gene_index.tolist(), list(range(8)))
        self.assertEqual(stats["n_dropped_nonfinite_log2fc"], 1)
        self.assertEqual(stats["n_unmeasured"], 6)
        self.assertEqual(int(stats["unmeasured_rank"]), 3)

    def test_dataset_deg_rank_rejects_out_of_axis_indices(self):
        with self.assertRaises(ValueError):
            deg_rank.dataset_deg_rank(
                pd.DataFrame({"gene_index": [8], "log2fc": [1.0], "padj": [0.1]}),
                axis_size=8)


class SymbolResolutionTest(unittest.TestCase):
    def setUp(self):
        self.resolver = deg_rank.SymbolResolver(
            axis_fixture(), crosswalk_fixture().to_dict("records"),
            hgnc_fixture().to_dict("records"))

    def test_exact_alias_ensemble_and_gene_index_space(self):
        self.assertEqual(self.resolver.resolve("G0"), (0, "canonical"))
        self.assertEqual(self.resolver.resolve("X1"), (1, "crosswalk_alias"))
        self.assertEqual(self.resolver.resolve("HALIAS"), (4, "hgnc_alias"))
        self.assertEqual(self.resolver.resolve("ENSG00000000005"), (5, "ensembl"))

    def test_ambiguous_offaxis_and_unmapped_are_recorded_not_guessed(self):
        self.assertEqual(self.resolver.resolve("AMB"), (-1, "ambiguous"))
        self.assertEqual(self.resolver.resolve("OFFAXIS"), (-1, "offaxis"))
        self.assertEqual(self.resolver.resolve("NOPE"), (-1, "unmapped"))

    def test_align_table_drops_unresolved_rows_and_reports_losses(self):
        rows = pd.DataFrame({"symbol": ["G0", "X1", "AMB", "OFFAXIS", "NOPE", "HALIAS"]})
        table, report = deg_rank.align_table(self.resolver, rows)
        self.assertEqual(table.gene_index.tolist(), [0, 1, 4])
        self.assertEqual(report["n_source_symbols"], 6)
        self.assertEqual(report["n_resolved"], 3)
        self.assertEqual(report["ambiguous"], 1)
        self.assertEqual(report["offaxis"], 1)
        self.assertEqual(report["unmapped"], 1)

    def test_alignment_accepts_an_alias_that_collapses_two_source_rows(self):
        rows = pd.DataFrame({"symbol": ["X1", "G1"],
                             "log2fc": [1.0, 3.0], "padj": [1e-4, 1e-3]})
        table, _ = deg_rank.align_table(self.resolver, rows)
        collapsed, _, stats = deg_rank.dataset_deg_rank(table, axis_size=8)
        self.assertEqual(len(collapsed), 1)
        self.assertEqual(collapsed.gene_index.tolist(), [1])
        self.assertEqual(int(collapsed.n_contrasts.iloc[0]), 2)
        self.assertEqual(stats["n_measured"], 1)


if __name__ == "__main__":
    unittest.main()
