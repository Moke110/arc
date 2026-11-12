"""Epoch-based early stopping at the ADR 0063 contract seam.

The tests state the two independent rules, the 100-epoch ceiling, the
lowest-validation-loss selection and the recorded dominance of the 5-epoch
rule over the 20-epoch rule.  Everything is driven by the public
``EarlyStopper.update`` contract with plain numbers; no model is needed.
"""

import sys
from pathlib import Path
import unittest

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP / 'execution/script'))

from early_stopping import RULES, DOMINANCE, EarlyStopper  # noqa: E402


def run(stopper, values):
    decisions = []
    for epoch, value in enumerate(values, start=1):
        decisions.append(stopper.update(epoch, value))
        if decisions[-1]['stop']:
            break
    return decisions


class DefaultRules(unittest.TestCase):
    def test_rules_match_adr_0063(self):
        self.assertEqual(RULES['no_strict_decrease_epochs'], 5)
        self.assertEqual(RULES['above_running_min_epochs'], 20)
        self.assertEqual(RULES['max_epochs'], 100)
        self.assertIn('5-epoch', DOMINANCE)
        self.assertIn('dominates', DOMINANCE)

    def test_strictly_improving_trajectory_runs_to_the_epoch_ceiling(self):
        stopper = EarlyStopper(max_epochs=4)
        decisions = run(stopper, [1.0, .9, .8, .7])
        last = decisions[-1]
        self.assertTrue(last['stop'])
        self.assertEqual(last['reason'], 'max_epochs')
        self.assertEqual(last['best_epoch'], 4)
        self.assertEqual(last['best_validation_loss'], .7)
        self.assertEqual(stopper.record()['observed']['no_strict_decrease'], 0)

    def test_five_consecutive_non_decreases_stop_at_the_fifth(self):
        # Epoch 1 sets the minimum; epochs 2..6 do not strictly decrease.
        stopper = EarlyStopper()
        decisions = run(stopper, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5])
        last = decisions[-1]
        self.assertTrue(last['stop'])
        self.assertEqual(last['reason'], 'no_strict_decrease')
        self.assertEqual(last['epoch'], 6)
        self.assertEqual(last['no_strict_decrease'], 5)
        self.assertEqual(last['best_epoch'], 1)
        self.assertAlmostEqual(last['best_validation_loss'], 1.0)

    def test_a_strict_improvement_resets_the_counter(self):
        stopper = EarlyStopper()
        decisions = run(stopper, [1.0, 1.1, 1.2, 1.3, 1.4, .9])
        self.assertFalse(any(d['stop'] for d in decisions))
        self.assertEqual(stopper.record()['observed']['no_strict_decrease'], 0)
        self.assertEqual(stopper.record()['observed']['best_epoch'], 6)

    def test_a_tie_counts_as_no_strict_decrease_but_not_as_above_the_minimum(self):
        stopper = EarlyStopper()
        run(stopper, [1.0] + [1.0] * 5)
        record = stopper.record()['observed']
        self.assertEqual(record['no_strict_decrease'], 5)
        self.assertEqual(record['above_running_min'], 0)

    def test_the_twenty_epoch_rule_is_independent(self):
        # Disable the 5-epoch rule to observe the 20-epoch rule alone.
        stopper = EarlyStopper(no_strict_decrease_epochs=1000,
                               above_running_min_epochs=20)
        decisions = run(stopper, [1.0] + [2.0] * 20)
        last = decisions[-1]
        self.assertTrue(last['stop'])
        self.assertEqual(last['reason'], 'above_running_min')
        self.assertEqual(last['epoch'], 21)
        self.assertEqual(last['above_running_min'], 20)
        self.assertEqual(last['best_epoch'], 1)

    def test_default_configuration_shows_the_five_epoch_rule_dominates(self):
        stopper = EarlyStopper()
        decisions = run(stopper, [1.0] + [2.0] * 25)
        last = decisions[-1]
        self.assertEqual(last['reason'], 'no_strict_decrease')
        self.assertEqual(last['epoch'], 6)
        self.assertLess(last['above_running_min'], RULES['above_running_min_epochs'])
        record = stopper.record()
        self.assertEqual(record['dominance'], DOMINANCE)
        self.assertTrue(record['dominance_observed'])

    def test_dominance_accounting_records_the_first_firing_epochs(self):
        stopper = EarlyStopper()
        run(stopper, [1.0] + [2.0] * 8)
        record = stopper.record()
        self.assertEqual(record['observed']['first_no_strict_decrease_epoch'], 6)
        self.assertIsNone(record['observed']['first_above_running_min_epoch'])
        self.assertTrue(record['dominance_observed'])

    def test_dominance_accounting_can_be_false_for_a_custom_recipe(self):
        # A 5-epoch threshold that never fires while a 3-epoch above-minimum
        # rule does: the recorded accounting must expose the inversion instead
        # of the old tautology `above <= no_strict_decrease`.
        stopper = EarlyStopper(no_strict_decrease_epochs=1000,
                               above_running_min_epochs=3)
        run(stopper, [1.0, 2.0, 2.0, 2.0])
        record = stopper.record()
        self.assertEqual(record['observed']['first_above_running_min_epoch'], 4)
        self.assertIsNone(record['observed']['first_no_strict_decrease_epoch'])
        self.assertFalse(record['dominance_observed'])
        self.assertEqual(record['stop_reason'], 'above_running_min')

    def test_invalid_recipe_is_rejected(self):
        for kwargs in ({'max_epochs': 0}, {'no_strict_decrease_epochs': 0},
                       {'above_running_min_epochs': 0}):
            with self.assertRaises(ValueError):
                EarlyStopper(**kwargs)
        with self.assertRaises(ValueError):
            EarlyStopper().update(1, float('nan'))


if __name__ == '__main__':
    unittest.main()
