from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[2]
for path in (STEP / 'output/api', STEP.parent / '04_tsg_dose/output/api'):
    sys.path.insert(0, str(path))

from tsg_transfer.tenv import export_test, freeze_split, load_split, verify_test


def fixture(root):
    obs = pd.DataFrame(
        [(f'row-{parent}-{dose}', f'parent-{parent}', dose, 'CM')
         for parent in range(12) for dose in (0.2, 1., 10.)],
        columns=['profile_id', 'drug_parent_inchikey', 'dose_um', 'cell_line'],
    ).set_index('profile_id', drop=False)
    (root / 'manifest.json').write_text('{}')
    return SimpleNamespace(obs=obs, source='tenv', directory=root)


class TenvSplitTest(unittest.TestCase):
    def test_three_folds_preserve_the_two_excluded_strips(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = fixture(root)
            split = root / 'split'
            freeze_split(data, split)
            held_parents = set()
            for fold, dose in enumerate((0.2, 1., 10.)):
                rows = load_split(data, split, fold)
                self.assertEqual(set(rows), {
                    'train', 'validation', 'test', 'excluded-held-parent', 'excluded-held-dose',
                })
                self.assertEqual(len(set(np.concatenate(list(rows.values())))), len(data.obs))
                observed = {role: data.obs.iloc[positions] for role, positions in rows.items()}
                held = set(observed['test'].drug_parent_inchikey)
                self.assertFalse(held_parents.intersection(held))
                held_parents.update(held)
                self.assertEqual(set(observed['test'].dose_um), {dose})
                self.assertEqual(set(observed['excluded-held-parent'].drug_parent_inchikey), held)
                self.assertEqual(set(observed['excluded-held-dose'].dose_um), {dose})
                self.assertFalse(set(observed['excluded-held-parent'].dose_um).intersection({dose}))
                self.assertFalse(set(observed['excluded-held-dose'].drug_parent_inchikey).intersection(held))
                for role in ('train', 'validation'):
                    self.assertFalse(set(observed[role].drug_parent_inchikey).intersection(held))
                    self.assertFalse(set(observed[role].dose_um).intersection({dose}))
            self.assertEqual(held_parents, set(data.obs.drug_parent_inchikey))

    def test_test_receipt_binds_predictions_to_ordered_test_profiles(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = fixture(root)
            data.expression = np.tile(np.array([1, 2, 5], dtype=np.float32), (len(data.obs), 1))
            data.control = np.zeros_like(data.expression)
            split = root / 'split'
            freeze_split(data, split)
            rows = load_split(data, split, 0)['test']
            output = root / 'test'
            lineage = {'fit_receipt_sha256': 'a' * 64}
            prediction = data.expression[rows].copy()
            prediction[1] = prediction[1][::-1]
            path = export_test(data, split, 0, output, prediction, path='tenv-only',
                               lineage=lineage, prediction_ids=data.obs.index[rows])
            receipt = verify_test(data, split, 0, output, path='tenv-only', lineage=lineage,
                                  require_replay=False)
            self.assertEqual(receipt['metrics']['n_samples'], 4)
            self.assertAlmostEqual(receipt['metrics']['pcc'], 7 / 13)
            self.assertEqual(path, output / 'manifest.json')
            with self.assertRaisesRegex(ValueError, 'replay'):
                verify_test(data, split, 0, output, path='tenv-only', lineage=lineage)
            with self.assertRaises(ValueError):
                export_test(data, split, 0, root / 'wrong', prediction, path='tenv-only',
                            lineage=lineage, prediction_ids=data.obs.index[rows][::-1])
            (output / 'prediction.npy').write_bytes(b'corrupt')
            with self.assertRaises(ValueError):
                verify_test(data, split, 0, output, path='tenv-only', lineage=lineage,
                            require_replay=False)

    def test_changed_split_bytes_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = fixture(root)
            split = root / 'split'
            freeze_split(data, split)
            with (split / 'splits.csv').open('a') as output:
                output.write('foreign identity\n')
            with self.assertRaises(ValueError):
                load_split(data, split, 0)


if __name__ == '__main__':
    unittest.main()
