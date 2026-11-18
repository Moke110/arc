"""Scientific acceptance of the paired held-parent/held-dose intersections."""
from pathlib import Path
import sys
import unittest

import pandas as pd

STEP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(STEP / 'output/api'))
sys.path.insert(0, str(STEP / 'execution/axis_local'))

from auxiliary import audit_coverage
from tsg_round1.protocol import make_splits


class AuxiliaryAcceptance(unittest.TestCase):
    def setUp(self):
        population = pd.DataFrame([
            {'profile_id': f'{parent}-{dose}-{line}', 'drug_parent_inchikey': parent,
             'dose_um': dose, 'cell_line': line}
            for parent in 'abcdef' for dose in (.05, .5, 5.) for line in ('a', 'b')
        ]).set_index('profile_id')
        self.ledger = make_splits(population)
        self.scored = {
            fold: self.ledger.query(
                'scenario == "cold-drug-dose" and `fold` == @fold and partition == "test"'
            )[['profile_id', 'drug_parent', 'dose_um']].reset_index(drop=True)
            for fold in range(3)
        }

    def test_coverage_records_only_intersections_and_the_two_excluded_strips(self):
        folds = audit_coverage(self.ledger, self.scored)
        self.assertEqual([fold['held_dose_um'] for fold in folds], [.05, .5, 5.])
        self.assertEqual([fold['n_test_samples'] for fold in folds], [4, 4, 4])
        self.assertEqual([fold['n_excluded_held_parent'] for fold in folds], [8, 8, 8])
        self.assertEqual([fold['n_excluded_held_dose'] for fold in folds], [8, 8, 8])
        self.assertEqual([fold['test_population_fraction'] for fold in folds], [1 / 9] * 3)

    def test_a_seen_dose_from_a_held_parent_cannot_enter_scored_rows(self):
        foreign = self.ledger.loc[
            self.ledger.scenario.eq('cold-drug-dose') & self.ledger.fold.eq(0)
            & self.ledger.partition.eq('excluded-held-parent'),
            ['profile_id', 'drug_parent', 'dose_um'],
        ].head(1)
        self.scored[0] = pd.concat([foreign, self.scored[0].iloc[1:]], ignore_index=True)
        with self.assertRaisesRegex(ValueError, 'fold 0 scored identities'):
            audit_coverage(self.ledger, self.scored)
