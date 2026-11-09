"""The dose intervention compares matching strata under loss-selected replay."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))

from dose_reset import ARM, paired_summaries
from screen import screen_decision


def summary(dose, positive, negative, mean):
    return {'dose_um': dose, 'group_count': 1, 'positive_ap100': positive,
            'negative_ap100': negative, 'mean_ap100': mean}


class DoseComparison(unittest.TestCase):
    def test_large_ten_micromolar_gain_does_not_override_overall_gate(self):
        arms = {
            ARM: [summary('pooled', 0.31, 0.19, 0.25), summary('10.0', 0.4, 0.4, 0.4)],
            'Transfer': [summary('pooled', 0.29, 0.2, 0.245), summary('10.0', 0.2, 0.2, 0.2)],
            'CM': [summary('pooled', 0.35, 0.25, 0.3), summary('10.0', 0.5, 0.5, 0.5)],
        }
        pooled, high = paired_summaries(arms)
        self.assertEqual(high['gain']['mean_ap100'], 0.2)
        self.assertEqual(pooled['gain'], {'positive_ap100': 0.02, 'negative_ap100': -0.01,
                                          'mean_ap100': 0.005})
        self.assertEqual(pooled['cm_minus_transfer']['mean_ap100'], 0.055)
        self.assertEqual(screen_decision({17: pooled['gain']['mean_ap100']})['status'],
                         'below_threshold')


if __name__ == '__main__':
    unittest.main()
