"""Check the published comparison against a worked three-fold example."""
import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from publish import comparison


class ComparisonTest(unittest.TestCase):
    def test_fold_means_weight_folds_equally_and_reject_unpaired_identities(self):
        values = ((0.1, -0.2, 0.4), (0.3, 0.2, 0.1), (0.5, 0.4, 0.6))
        receipts = [
            dict(fold=fold, path=path, files={'identity.csv': f'identities-{fold}'},
                 metrics=dict(pcc=score, n_samples=3 ** fold, undefined_pcc_count=0))
            for fold, scores in enumerate(values)
            for path, score in zip(('tenv-only', 'tahoe-only', 'transfer'), scores, strict=True)
        ]
        scores, differences = comparison(receipts)
        means = {row['path']: row['pcc'] for row in scores if row['fold'] == 'mean'}
        self.assertAlmostEqual(means['tenv-only'], 0.3)
        self.assertAlmostEqual(means['tahoe-only'], 0.4 / 3)
        self.assertAlmostEqual(means['transfer'], 1.1 / 3)
        self.assertAlmostEqual(differences[-1]['transfer_minus_tenv_only'], 0.2 / 3)
        self.assertAlmostEqual(differences[-1]['transfer_minus_tahoe_only'], 0.7 / 3)
        receipts[-1]['files']['identity.csv'] = 'different'
        with self.assertRaisesRegex(ValueError, 'ordered test identities'):
            comparison(receipts)


if __name__ == '__main__':
    unittest.main()
