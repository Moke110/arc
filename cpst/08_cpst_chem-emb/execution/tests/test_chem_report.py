"""Exercise pairing and standard deviation at the report boundary."""

from pathlib import Path
import sys
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from chem_report import METRICS, paired_scores


def observations(roster):
    values = {("KPGT", 17): .4, ("KPGT", 42): .6, ("UniMol", 17): .5, ("UniMol", 42): .9}
    return pd.DataFrame([{"arm": arm, "fold": fold, "seed": seed, "metric": metric,
                          "value": values[arm, seed], "n_conditions": 4, "n_eligible": 4, "n_undefined": 0}
                         for arm, fold, seed in reversed(roster) for metric in METRICS])


class ChemicalReportTest(unittest.TestCase):
    def test_difference_uses_matching_fold_and_seed_and_sample_sd(self):
        roster = [("KPGT", 0, 17), ("KPGT", 0, 42), ("UniMol", 0, 17), ("UniMol", 0, 42)]
        paired, by_fold, summary = paired_scores(observations(roster), roster)
        self.assertEqual(paired.loc[paired.arm.eq("UniMol")].seed.unique().tolist(), [17, 42])
        for table in (by_fold, summary):
            record = table.loc[table.arm.eq("UniMol") & table.metric.eq("Δ PCC")].iloc[0]
            self.assertAlmostEqual(record.score_difference, .2)
            self.assertAlmostEqual(record.score_difference_sd, .2 / 2**.5)
            self.assertEqual(record.n_paired, 2)

    def test_report_rejects_coordinate_without_kpgt_reference(self):
        roster = [("KPGT", 0, 17), ("UniMol", 0, 17), ("UniMol", 0, 42)]
        with self.assertRaisesRegex(ValueError, "KPGT"):
            paired_scores(observations(roster), roster)


if __name__ == "__main__":
    unittest.main()
