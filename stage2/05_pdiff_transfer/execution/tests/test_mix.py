"""Exercise the pooled predictor input through its artifact and batch interface."""
import subprocess
import tempfile
from pathlib import Path
import sys
import unittest

sys.dont_write_bytecode = True
STEP = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(STEP / 'execution/script'),
                str(STEP.parent / '03_modes_tenv/execution/script')]

import numpy as np

from artifacts import sha256, write_csv, write_json
from mix_inputs import MixInputs


class MixTest(unittest.TestCase):
    def test_pool_preserves_source_control_pairs_and_exposes_only_tenv_holdout(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = self.source(root / 'tenv', ('train', 'validation', 'test'), 10, '0.2')
            source = self.source(root / 'tahoe', ('train', 'validation', 'test'), 50, '0.5')
            for name in ('validation_cases', 'validation_expression', 'test_cases'):
                source[name].unlink()
            pooled = MixInputs({'tenv': target, 'tahoe': source, 'n_genes': 2})

            self.assertEqual([case['group_id'] for case in pooled.cases['train']],
                             ['tenv:group-0', 'tahoe:group-0'])
            self.assertEqual([case['case_id'] for case in pooled.cases['validation']], ['case-1'])
            self.assertEqual([case['case_id'] for case in pooled.cases['test']], ['case-2'])
            control, fingerprint, dose, treated = pooled.batch('train', [1, 0, 1], 'cpu')
            np.testing.assert_array_equal(control, [[51, 52], [11, 12], [51, 52]])
            np.testing.assert_array_equal(fingerprint[:, 0], [50, 10, 50])
            np.testing.assert_allclose(dose, [0.5, 0.2, 0.5])
            np.testing.assert_array_equal(treated, [[56, 58], [16, 18], [56, 58]])
            self.assertEqual(len(pooled.batch('test', [0], 'cpu')), 3)

    def test_fit_rejects_a_changed_plan_before_opening_sources(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            plan = root / 'plan.json'
            write_json(plan, {'stage': 'formal'})
            digest = sha256(plan)
            write_json(plan, {'stage': 'smoke'})
            result = subprocess.run([sys.executable, str(STEP / 'execution/script/mix.py'),
                                     'fit', str(root), 'cold-drug', '17', '--plan-sha256', digest],
                                    capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 1)
            self.assertIn('frozen artifact path/hash mismatch: plan.json', result.stderr)
            self.assertFalse((root / 'slots').exists())

    @staticmethod
    def source(root, roles, offset, dose):
        root.mkdir()
        files = {key: root / f'{key}.npy' for key in ('controls', 'morgan')}
        files['drugs'] = root / 'drugs.csv'
        write_csv(files['drugs'], [{'drug_key': f'drug-{index}', 'morgan_index': index}
                                  for index in range(len(roles))])
        np.save(files['controls'], np.array([[1, 2], [offset + 1, offset + 2]], dtype='float32'))
        np.save(files['morgan'], np.full((3, 2048), offset, dtype='float32'))
        for index, role in enumerate(roles):
            files[f'{role}_cases'] = root / f'{role}.csv'
            write_csv(files[f'{role}_cases'], [{'case_id': f'case-{index}', 'drug_key': f'drug-{index}',
                                              'group_id': f'group-{index}', 'dose_um': dose,
                                              'control_index': 1}])
            files[f'{role}_expression'] = root / f'{role}.npy'
            if role != 'test':
                np.save(files[f'{role}_expression'], np.array([[5, 6]], dtype='float32'))
        return files


if __name__ == '__main__':
    unittest.main()
