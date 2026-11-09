"""The Step 06 screen uses overall, paired validation effects."""
from decimal import Decimal
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))

from screen import paired_gain, screen_decision


class ScreenDecision(unittest.TestCase):
    def test_exact_gate_only_promotes_positive_replication(self):
        effect = paired_gain(0.29, 0.28)
        self.assertEqual(effect, Decimal('0.010'))
        self.assertEqual(screen_decision({17: effect})['next_seeds'], [42, 2026])
        self.assertEqual(screen_decision({17: effect, 42: 0.001, 2026: 0})['status'], 'inconsistent')
        self.assertEqual(screen_decision({17: effect, 42: 0.001, 2026: None})['status'], 'incomplete')
        self.assertEqual(screen_decision({17: effect, 42: 0.001, 2026: 0.002})['status'], 'promoted')
        self.assertEqual(screen_decision({17: 0.009999})['status'], 'below_threshold')


if __name__ == '__main__':
    unittest.main()
