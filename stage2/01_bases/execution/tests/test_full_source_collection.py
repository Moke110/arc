"""ADR0073 collection through the public ledger/catalog boundary."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import public_api
from bases import disease_deg

class FullCollectionTests(unittest.TestCase):
    def test_later_tiers_contribute_after_two_thousand_degs(self):
        declared = [dict(dataset_id='a', tier=1), dict(dataset_id='b', tier=2),
                    dict(dataset_id='c', tier=3)]
        measured = {name: dict(deg=genes, n_deg=len(genes), informative=True)
                    for name, genes in [('a', list(range(2100))),
                                        ('b', [2099, 2100]), ('c', [2101])]}
        result = disease_deg.tier_ledger('genetic_cardiomyopathy', declared, measured)
        self.assertEqual(result['n_deg'], 2102)
        self.assertEqual(result['tiering']['tiers_entered'], [1, 2, 3])
        self.assertEqual(result['tiering']['tier_union'], [
            dict(tier=1,n_deg_union=2100),dict(tier=2,n_deg_union=2101),
            dict(tier=3,n_deg_union=2102)])
        self.assertEqual(result['tiering']['collection_policy'], 'all_qualifying_declared_sources')
        self.assertTrue({'stop_reached','stopped_at_tier','limit'}.isdisjoint(result['tiering']))

if __name__ == '__main__': unittest.main()
