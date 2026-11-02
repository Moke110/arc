"""Rule-engine tests: thresholds are literals taken from the locked rule set.

Every expected value here is written by hand from the rule text, never recomputed
with the production helper (no tautological assertions).
"""

import math
import unittest

import numpy as np
import pandas as pd

import rules


def _tempenv_frame(**over):
    """One treated / one DMSO / one media / one blank well, all QC-clean."""
    base = {
        "well": ["temp_env::Plate1-B02", "temp_env::Plate1-B16", "temp_env::Plate1-B17",
                 "temp_env::Plate1-B12", "temp_env::Plate1-B03"],
        "role": ["vehicle_control_dmso", "vehicle_control_dmso", "control_media",
                 "technical_blank", "treated_compound"],
        "vehicle": ["DMSO", "DMSO", "Media", "Blank", "DMSO"],
        "lib_size_panel": [900_000, 500_000, 700_000, 300_000, 800_000],
        "n_genes_panel": [12_000, 11_000, 11_500, 7_000, 11_800],
        "meanr_vs_vehicle_pool": [0.95, 0.90, 0.91, 0.05, 0.93],
    }
    base.update(over)
    return pd.DataFrame(base)


REFS = {
    "tempenv_gene_filter_min_mean": 10.0,
    "tempenv_n_genes_retained": 11_497.0,
    "tempenv_sample_expression_floor": 15.0,
    "tempenv_vehicle_corr_floor": 0.85,
    "tempenv_blank_ngenes_threshold": 9242.0,
    "tempenv_blank_ngenes_median": 7067.5,
    "tempenv_blank_ngenes_mad": 489.0,
}


class TestModifiedZ(unittest.TestCase):
    def test_known_values(self):
        # [1,2,3,4,100]: median 3, |dev| = [2,1,0,1,97], MAD = 1 -> scale 1.4826
        z = rules.modified_z(np.array([1.0, 2.0, 3.0, 4.0, 100.0]))
        self.assertAlmostEqual(z[0], (1 - 3) / (1.4826 * 1.0), places=9)
        self.assertAlmostEqual(z[2], 0.0, places=12)
        self.assertAlmostEqual(z[4], (100 - 3) / (1.4826 * 1.0), places=9)

    def test_zero_mad_is_nan_not_inf(self):
        z = rules.modified_z(np.array([5.0, 5.0, 5.0]))
        self.assertTrue(np.all(np.isnan(z)))


class TestAuthorGeneFilter(unittest.TestCase):
    def test_mean_threshold_is_per_gene_over_all_wells(self):
        # gene 0: total 100 over 10 wells -> mean 10 (kept, boundary is >= 10)
        # gene 1: total 99  over 10 wells -> mean 9.9 (dropped)
        totals = np.array([100.0, 99.0, 10_000.0])
        mask = rules.author_gene_filter(totals, n_wells=10, min_mean=10.0)
        self.assertEqual(mask.tolist(), [True, False, True])


class TestBlankSeparation(unittest.TestCase):
    def test_threshold_is_median_plus_three_mad_scales(self):
        blanks = np.array([7000.0, 7100.0, 6900.0, 7050.0])
        med, mad, thr = rules.blank_separation(blanks, k=3.0)
        # median of [6900,7000,7050,7100] is 7025; deviations [125,25,25,75] give MAD 50
        self.assertAlmostEqual(med, 7025.0)
        self.assertAlmostEqual(mad, 50.0)
        self.assertAlmostEqual(thr, 7025.0 + 3.0 * 1.4826 * 50.0, places=9)


class TestTempEnvRules(unittest.TestCase):
    def test_clean_wells_pass(self):
        out = rules.evaluate_tempenv(_tempenv_frame(), REFS)
        self.assertEqual(out["qc_status"].tolist(), ["pass"] * 5)
        self.assertTrue((out["qc_reason"] == "").all())

    def test_low_mean_expression_excludes_any_role(self):
        df = _tempenv_frame(lib_size_panel=[20_000, 500_000, 700_000, 300_000, 800_000])
        out = rules.evaluate_tempenv(df, REFS)
        # 20_000 / 11_497 = 1.74 < 15
        self.assertEqual(out["qc_status"].iloc[0], "fail")
        self.assertIn("src_sample_expression_floor", out["qc_reason"].iloc[0])
        self.assertEqual(out["qc_status"].iloc[1], "pass")

    def test_expression_floor_is_strict_at_the_boundary(self):
        # exactly 15.0 mean -> kept (rule is strictly lower than 15)
        df = _tempenv_frame(lib_size_panel=[15.0 * 11_497, 500_000, 700_000, 300_000, 800_000])
        out = rules.evaluate_tempenv(df, REFS)
        self.assertEqual(out["qc_status"].iloc[0], "pass")

    def test_vehicle_correlation_floor_applies_to_dmso_controls_only(self):
        df = _tempenv_frame(meanr_vs_vehicle_pool=[0.80, 0.90, 0.80, 0.05, 0.80])
        out = rules.evaluate_tempenv(df, REFS)
        # DMSO well r=0.80 -> excluded as an outlier vehicle control
        self.assertEqual(out["qc_status"].iloc[0], "fail")
        self.assertIn("src_vehicle_reproducibility", out["qc_reason"].iloc[0])
        # media / blank / treated below the same floor are NOT excluded by T2
        self.assertNotIn("src_vehicle_reproducibility", out["qc_reason"].iloc[2])
        self.assertNotIn("src_vehicle_reproducibility", out["qc_reason"].iloc[4])

    def test_vehicle_correlation_floor_is_flagged_not_enforced_off_controls(self):
        df = _tempenv_frame(meanr_vs_vehicle_pool=[0.95, 0.90, 0.91, 0.05, 0.80])
        out = rules.evaluate_tempenv(df, REFS)
        self.assertEqual(out["qc_status"].iloc[4], "pass")
        self.assertIn("vehicle_corr_below_source_floor", out["qc_flags"].iloc[4])

    def test_blank_separation_floor_excludes_near_blank_library(self):
        df = _tempenv_frame(n_genes_panel=[12_000, 11_000, 11_500, 8_000, 9_000])
        out = rules.evaluate_tempenv(df, REFS)
        # 9,000 <= 9,242 threshold (strictly greater keeps)
        self.assertEqual(out["qc_status"].iloc[4], "fail")
        self.assertIn("derived_blank_separation_floor", out["qc_reason"].iloc[4])
        self.assertEqual(out["qc_status"].iloc[0], "pass")

    def test_boundary_of_blank_floor_is_kept(self):
        df = _tempenv_frame(n_genes_panel=[12_000, 11_000, 11_500, 8_000, 9_243])
        out = rules.evaluate_tempenv(df, REFS)
        self.assertEqual(out["qc_status"].iloc[4], "pass")

    def test_missing_metric_is_unresolved_not_pass(self):
        df = _tempenv_frame(lib_size_panel=[np.nan, 500_000, 700_000, 300_000, 800_000])
        out = rules.evaluate_tempenv(df, REFS)
        self.assertEqual(out["qc_status"].iloc[0], "unresolved")
        self.assertIn("unresolved_missing_metric", out["qc_reason"].iloc[0])

    def test_failure_lists_every_failing_reason(self):
        # well 4 (treated) fails T1 (20,000 / 11,497 = 1.74) and T3 (6,000 <= 9,242)
        df = _tempenv_frame(lib_size_panel=[900_000, 500_000, 700_000, 300_000, 20_000],
                            n_genes_panel=[12_000, 11_000, 11_500, 7_000, 6_000])
        out = rules.evaluate_tempenv(df, REFS)
        reasons = out["qc_reason"].iloc[4]
        self.assertIn("src_sample_expression_floor", reasons)
        self.assertIn("derived_blank_separation_floor", reasons)
        self.assertEqual(out["qc_status"].iloc[3], "pass")  # blank exempt from T3


def _tahoe_frame(**over):
    base = {
        "observation_id": ["tahoe::plate1:smp_1495:CVCL_0399",
                           "tahoe::plate1:smp_1496:CVCL_0399",
                           "tahoe::plate1:smp_1495:CVCL_0058"],
        "role": ["treated_compound", "vehicle_control_dmso_tf", "treated_compound"],
        "cell_name": ["A549", "A549", "A-427"],
        "n_cells": [1500, 1200, 40],
        "mean_gene_count": [1400.0, 1350.0, 1200.0],
        "mean_tscp_count": [2300.0, 2250.0, 2100.0],
        "mean_pcnt_mito": [0.06, 0.07, 0.05],
        "well_z_tscp": [0.4, -0.2, 0.1],
        "has_structure": [True, False, True],
    }
    base.update(over)
    return pd.DataFrame(base)


TREFS = {
    "tahoe_n_cells_min": 100,
    "tahoe_source_n_cells_floor": 50,
    "tahoe_excluded_cell_lines": ["NCI-H2122", "NCI-H596", "NCI-H661"],
    "tahoe_mito_max_pcnt": 20.0,
    "tahoe_min_mean_gene_count": 250.0,
    "tahoe_well_z_abs_max": 3.0,
}


class TestTahoeRules(unittest.TestCase):
    def test_source_rules_pass_and_low_cell_count_fails(self):
        out = rules.evaluate_tahoe(_tahoe_frame(), TREFS)
        self.assertEqual(out["qc_status"].tolist(), ["pass", "pass", "fail"])
        self.assertIn("project_min_cells_per_profile", out["qc_reason"].iloc[2])

    def test_source_floor_50_is_flagged_for_sensitivity(self):
        out = rules.evaluate_tahoe(_tahoe_frame(n_cells=[60, 1200, 1500]), TREFS)
        # 60 < 100 so H1 fails, but >= the source floor of 50 and therefore flagged
        self.assertEqual(out["qc_status"][0], "fail")
        self.assertIn("project_min_cells_per_profile", out["qc_reason"].iloc[0])
        self.assertIn("n_cells_between_source_and_project_floor", out["qc_flags"].iloc[0])

    def test_excluded_cell_line_fails_with_source_reason(self):
        out = rules.evaluate_tahoe(
            _tahoe_frame(cell_name=["NCI-H661", "A549", "A549"]), TREFS)
        self.assertEqual(out["qc_status"].iloc[0], "fail")
        self.assertIn("src_cell_line_excluded", out["qc_reason"].iloc[0])

    def test_mito_over_twenty_percent_fails(self):
        out = rules.evaluate_tahoe(_tahoe_frame(mean_pcnt_mito=[25.0, 7.0, 5.0]), TREFS)
        self.assertEqual(out["qc_status"].iloc[0], "fail")
        self.assertIn("src_mito_fraction_ceiling", out["qc_reason"].iloc[0])

    def test_low_mean_gene_count_fails(self):
        out = rules.evaluate_tahoe(_tahoe_frame(mean_gene_count=[100.0, 1350.0, 1200.0]), TREFS)
        self.assertEqual(out["qc_status"].iloc[0], "fail")
        self.assertIn("src_min_genes_per_cell", out["qc_reason"].iloc[0])

    def test_per_plate_robust_z_outlier_fails(self):
        out = rules.evaluate_tahoe(_tahoe_frame(well_z_tscp=[4.2, 0.0, 0.0]), TREFS)
        self.assertEqual(out["qc_status"].iloc[0], "fail")
        self.assertIn("derived_plate_depth_outlier", out["qc_reason"].iloc[0])

    def test_no_structure_is_flagged_not_failed(self):
        out = rules.evaluate_tahoe(_tahoe_frame(), TREFS)
        self.assertEqual(out["qc_status"].iloc[1], "pass")
        self.assertIn("no_structure", out["qc_flags"].iloc[1])

    def test_tahoe_missing_metric_is_unresolved(self):
        out = rules.evaluate_tahoe(_tahoe_frame(n_cells=[None, 1200, 1500]), TREFS)
        self.assertEqual(out["qc_status"].iloc[0], "unresolved")


class TestRulesHash(unittest.TestCase):
    def test_hash_is_stable_and_threshold_sensitive(self):
        h1 = rules.canonical_rules_sha256()
        h2 = rules.canonical_rules_sha256()
        self.assertEqual(h1, h2)
        self.assertEqual(len(h1), 64)
        import copy
        tweaked = copy.deepcopy(rules.RULES)
        tweaked["tempenv"]["thresholds"]["sample_expression_floor"] = 15.0001
        h3 = rules.canonical_rules_sha256(tweaked)
        self.assertNotEqual(h1, h3)

    def test_every_rule_carries_evidence_and_no_forbidden_basis(self):
        forbidden = ("deg count", "deg_count", "response strength", "model performance",
                     "auc", "auprc")
        for src, block in rules.RULES.items():
            for name, rule in block["rules"].items():
                with self.subTest(rule=name):
                    self.assertIn("evidence", rule)
                    self.assertTrue(rule["evidence"]["kind"])
                    self.assertTrue(rule["evidence"]["citation"])
                    text = (rule["evidence"].get("quote", "") or "").lower()
                    for word in forbidden:
                        self.assertNotIn(word, text)


if __name__ == "__main__":
    unittest.main()
