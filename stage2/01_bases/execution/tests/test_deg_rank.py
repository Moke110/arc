"""Disease-DEG table/rank contract: alignment, effect score, collapse, DEG ranks.

These are the pure functions the disease-DEG builders import; every case here is
a hand-worked fixture, not a re-derivation from the implementation's internals.
ADR 0071 v2: the DEG flag is ``padj < 0.05`` and ``|log2FC| >= 2`` in at least one
contrast, and the within-dataset rank orders the dataset's DEGs only (a non-DEG
gets ``n_DEG + 1``, ``rank_norm = rank / n_DEG``). The measured-gene rank and its
``n_measured`` denominator are retired, so no case here asserts them.
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


class DegFlagTest(unittest.TestCase):
    """ADR 0071 *DEG definition*: padj < 0.05 AND |log2FC| >= 2, per contrast."""

    def test_threshold_boundaries_are_strict_padj_and_inclusive_fold_change(self):
        log2fc = np.array([3.0, 3.0, 2.0, 1.999, -2.0, -1.999])
        padj = np.array([0.05, 0.049999, 0.01, 1e-6, 0.01, 1e-6])
        flag = deg_rank.deg_flag(log2fc, padj)
        # padj == 0.05 is not a DEG (strict); padj just below is one.
        self.assertFalse(bool(flag[0]))
        self.assertTrue(bool(flag[1]))
        # |log2FC| == 2 is a DEG (inclusive); just under it is not.
        self.assertTrue(bool(flag[2]))
        self.assertFalse(bool(flag[3]))
        # Sign is carried, not filtered: a down-regulated DEG is a DEG.
        self.assertTrue(bool(flag[4]))
        self.assertFalse(bool(flag[5]))

    def test_na_and_non_finite_padj_are_not_degs(self):
        log2fc = np.array([5.0, 5.0, 5.0, 5.0])
        padj = np.array([np.nan, np.inf, -np.inf, np.nan])
        self.assertFalse(deg_rank.deg_flag(log2fc, padj).any())

    def test_non_finite_log2fc_is_not_a_deg(self):
        log2fc = np.array([np.nan, np.inf, -np.inf])
        padj = np.array([1e-6, 1e-6, 1e-6])
        self.assertFalse(deg_rank.deg_flag(log2fc, padj).any())

    def test_both_conditions_are_required_in_the_same_contrast(self):
        # Large fold change but non-significant, and significant but small:
        # neither alone is a DEG.
        flag = deg_rank.deg_flag(np.array([4.0, 0.5]), np.array([0.4, 1e-4]))
        self.assertFalse(bool(flag[0]))
        self.assertFalse(bool(flag[1]))

    def test_na_padj_is_non_significant_not_a_missing_row(self):
        collapsed = deg_rank.per_gene_deg(
            np.array([0]), np.array([1.0]), np.array([np.nan]))
        self.assertEqual(len(collapsed), 1)
        self.assertEqual(float(collapsed.effect.iloc[0]), 0.0)
        self.assertFalse(bool(collapsed.significant.iloc[0]))
        self.assertFalse(bool(collapsed.deg.iloc[0]))
        self.assertEqual(int(collapsed.n_contrasts_deg.iloc[0]), 0)
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
        # Contrast b is also a DEG (padj 1e-3, |log2FC| 3): one DEG contrast.
        self.assertTrue(bool(gene5.deg))
        self.assertEqual(int(gene5.n_contrasts_deg), 1)

    def test_deg_is_pass_any_contrast_and_collapse_keeps_max_effect(self):
        # Contrast a: significant small effect. Contrast b: large but not
        # significant. Contrast c: significant large effect -> the gene is a DEG
        # and the collapse still reports contrast c (the maximum effect).
        table, _, stats = deg_rank.dataset_deg_rank(pd.DataFrame({
            "gene_index": [2, 2, 2],
            "log2fc": [0.5, 4.0, -2.1],
            "padj": [0.02, 0.5, 0.03],
            "contrast": ["a", "b", "c"],
        }), axis_size=4)
        row = table.iloc[0]
        self.assertEqual(row.contrast, "c")
        self.assertTrue(bool(row.deg))
        self.assertEqual(int(row.n_contrasts_deg), 1)
        self.assertEqual(int(row.n_contrasts_supporting), 2)
        self.assertEqual(int(stats["n_deg"]), 1)
        np.testing.assert_allclose(row.effect, 2.1 * -np.log10(0.03))


class DegRankTest(unittest.TestCase):
    """ADR 0071 decision 3: DEGs rank 1..n_DEG, non-DEGs n_DEG + 1."""

    def test_deg_only_rank_orders_degs_and_fills_non_degs(self):
        rank = deg_rank.deg_only_rank(
            np.array([4, 2, 9, 1]),
            np.array([1.0, 9.0, 4.0, 6.0]),
            np.array([True, True, False, False]))
        # DEGs are genes 4 and 2; non-DEGs (9, 1) get n_DEG + 1 = 3.
        self.assertEqual(rank.tolist(), [2, 1, 3, 3])

    def test_deg_rank_ties_break_on_gene_index(self):
        # Two DEGs with the same effect score: the lower gene_index wins.
        rank = deg_rank.deg_only_rank(
            np.array([7, 3]), np.array([4.0, 4.0]), np.array([True, True]))
        self.assertEqual(rank.tolist(), [2, 1])

    def test_dataset_deg_rank_is_deg_only_and_uses_the_n_deg_denominator(self):
        # Two DEGs (genes 0 and 4) and one measured non-DEG (gene 2).
        # Effects: gene0 = 3*4 = 12, gene4 = 2*3 = 6, gene2 = 0.5*2 = 1.
        table, full, stats = deg_rank.dataset_deg_rank(pd.DataFrame({
            "gene_index": [0, 4, 2],
            "log2fc": [3.0, 2.0, 0.5],
            "padj": [1e-4, 1e-3, 0.01],
        }), axis_size=8)
        self.assertEqual(stats["n_measured"], 3)
        self.assertEqual(stats["n_deg"], 2)
        self.assertEqual(stats["rank_norm_denominator"], "n_deg")
        self.assertEqual(stats["rank_semantics"], deg_rank.RANK_SEMANTICS)
        self.assertEqual(stats["rank_semantics"], "deg-only-v2")
        self.assertEqual(stats["nondeg_rank"], 3)
        self.assertEqual(table.gene_index.tolist(), [0, 4, 2])
        self.assertEqual(table["rank"].tolist(), [1, 2, 3])
        np.testing.assert_allclose(table["rank_norm"], [0.5, 1.0, 1.5])
        self.assertEqual(table["deg"].tolist(), [True, True, False])
        # The full-axis frame covers every gene: the two DEGs at 1..2 and
        # every other gene at n_DEG + 1 = 3.
        self.assertEqual(full["rank"].tolist(), [1, 3, 3, 3, 2, 3, 3, 3])
        np.testing.assert_allclose(full["rank_norm"].tolist(),
                                   [0.5, 1.5, 1.5, 1.5, 1.0, 1.5, 1.5, 1.5])

    def test_deg_rank_does_not_depend_on_the_measured_gene_count(self):
        degs = pd.DataFrame({"gene_index": [0, 4], "log2fc": [3.0, 2.0],
                             "padj": [1e-4, 1e-3]})
        extra = pd.DataFrame({"gene_index": [1, 2, 3, 5, 6, 7],
                              "log2fc": [0.2, 0.1, -0.3, 0.0, 0.4, -0.2],
                              "padj": [0.9, 0.8, 0.7, 1.0, 0.95, 0.6]})
        only_degs, _, _ = deg_rank.dataset_deg_rank(degs, axis_size=8)
        with_extra, _, stats = deg_rank.dataset_deg_rank(
            pd.concat([degs, extra], ignore_index=True), axis_size=8)
        self.assertEqual(stats["n_measured"], 8)
        self.assertEqual(stats["n_deg"], 2)
        # The same two DEGs keep ranks 1 and 2 and rank_norm 0.5 / 1.0: the
        # measured-gene count never enters the denominator.
        self.assertEqual(with_extra[with_extra.deg].gene_index.tolist(), [0, 4])
        self.assertEqual(with_extra[with_extra.deg]["rank"].tolist(),
                         only_degs["rank"].tolist())
        np.testing.assert_allclose(with_extra[with_extra.deg]["rank_norm"].to_numpy(),
                                   only_degs["rank_norm"].to_numpy())
        # Every non-DEG, measured or not, gets n_DEG + 1 = 3.
        self.assertTrue((with_extra[~with_extra.deg]["rank"] == 3).all())

    def test_ties_break_on_gene_index_in_the_dataset_rank(self):
        # Genes 1 and 2 are DEGs with identical effect scores (|log2FC|=2,
        # padj=1e-4 -> 8): the lower gene_index ranks first.
        table, _, _ = deg_rank.dataset_deg_rank(pd.DataFrame({
            "gene_index": [2, 1],
            "log2fc": [2.0, -2.0],
            "padj": [1e-4, 1e-4],
        }), axis_size=4)
        self.assertEqual(table.gene_index.tolist(), [1, 2])
        self.assertEqual(table["rank"].tolist(), [1, 2])

    def test_zero_deg_dataset_is_non_informative_and_does_not_crash(self):
        # Every padj is 1, so no gene passes the v2 definition: every rank is
        # n_DEG + 1 = 1, rank_norm has no defined denominator (NaN), and the
        # dataset records itself non-informative (ADR 0068, ticket #201).
        table, full, stats = deg_rank.dataset_deg_rank(pd.DataFrame({
            "gene_index": [2, 0, 1],
            "log2fc": [1.0, -1.0, 0.5],
            "padj": [1.0, 1.0, 1.0],
        }), axis_size=8)
        self.assertEqual(stats["n_deg"], 0)
        self.assertEqual(stats["nondeg_rank"], 1)
        self.assertEqual(table["rank"].tolist(), [1, 1, 1])
        self.assertTrue(np.isnan(table["rank_norm"]).all())
        self.assertEqual(full["rank"].tolist(), [1] * 8)
        self.assertTrue(np.isnan(full["rank_norm"]).all())
        self.assertFalse(stats["informative"])
        self.assertEqual(stats["n_distinct_effect"], 1)
        self.assertIn("is a DEG", stats["informativeness_reason"])

    def test_zero_deg_with_varied_effects_is_still_non_informative(self):
        # Two distinct effect scores but neither gene is a DEG (small effects):
        # under v2 the rank carries no ordering, so the dataset is not a source.
        _, _, stats = deg_rank.dataset_deg_rank(pd.DataFrame({
            "gene_index": [0, 1],
            "log2fc": [1.5, -1.0],
            "padj": [0.01, 1e-4],
        }), axis_size=4)
        self.assertEqual(stats["n_deg"], 0)
        self.assertEqual(stats["n_distinct_effect"], 2)
        self.assertFalse(stats["informative"])
        self.assertIn("is a DEG", stats["informativeness_reason"])

    def test_exactly_one_deg_dataset_behaves(self):
        _, full, stats = deg_rank.dataset_deg_rank(pd.DataFrame({
            "gene_index": [3, 5],
            "log2fc": [3.0, 0.2],
            "padj": [1e-4, 0.9],
        }), axis_size=8)
        self.assertEqual(stats["n_deg"], 1)
        self.assertEqual(stats["nondeg_rank"], 2)
        self.assertEqual(full["rank"].tolist(), [2, 2, 2, 1, 2, 2, 2, 2])
        np.testing.assert_allclose(full["rank_norm"].tolist(),
                                   [2.0, 2.0, 2.0, 1.0, 2.0, 2.0, 2.0, 2.0])
        self.assertTrue(stats["informative"])

    def test_determinism_identical_input_gives_identical_output(self):
        rows = pd.DataFrame({
            "gene_index": [7, 7, 2, 0, 5],
            "log2fc": [1.0, -4.0, 0.5, np.nan, 2.5],
            "padj": [0.5, 1e-5, 0.2, 0.1, 0.02],
        })
        first = deg_rank.dataset_deg_rank(rows, axis_size=8)
        second = deg_rank.dataset_deg_rank(rows, axis_size=8)
        pd.testing.assert_frame_equal(first[0], second[0])
        pd.testing.assert_frame_equal(first[1], second[1])
        self.assertEqual(first[2], second[2])
        table, full, stats = first
        self.assertEqual(table.gene_index.tolist(), [7, 5, 2])
        self.assertEqual(table["rank"].tolist(), [1, 2, 3])
        self.assertEqual(full.gene_index.tolist(), list(range(8)))
        self.assertEqual(stats["n_dropped_nonfinite_log2fc"], 1)
        self.assertEqual(stats["n_unmeasured"], 5)
        self.assertEqual(stats["n_deg"], 2)
        self.assertEqual(stats["nondeg_rank"], 3)

    def test_dataset_deg_rank_rejects_out_of_axis_indices(self):
        with self.assertRaises(ValueError):
            deg_rank.dataset_deg_rank(
                pd.DataFrame({"gene_index": [8], "log2fc": [1.0], "padj": [0.1]}),
                axis_size=8)

    def test_all_tied_effect_scores_are_non_informative(self):
        result = deg_rank.dataset_informativeness(np.zeros(4))
        self.assertFalse(result["informative"])
        self.assertEqual(result["n_distinct_effect"], 1)
        self.assertIn("tied", result["reason"])
        informative = deg_rank.dataset_informativeness(np.array([0.0, 0.0, 1.0]))
        self.assertTrue(informative["informative"])
        self.assertIsNone(informative["reason"])


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
        # The collapsed gene's winning row has |log2FC| 3 at padj 1e-3: a DEG.
        self.assertTrue(bool(collapsed.deg.iloc[0]))
        self.assertEqual(stats["n_deg"], 1)


if __name__ == "__main__":
    unittest.main()
