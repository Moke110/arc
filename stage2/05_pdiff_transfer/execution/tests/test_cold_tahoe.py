"""Check the filtered Tahoe release boundary without the real data."""
import csv
import tempfile
import unittest
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))

from cold_tahoe import verify_roles
from training_inputs import predictor_paths


class ColdTahoeTests(unittest.TestCase):
    def test_rejects_forbidden_parent_in_either_source_role(self):
        fields = ('case_id', 'drug_key', 'dose_um', 'control_index', 'group_id')
        paths = predictor_paths(0)
        original = {
            'train': [dict(zip(fields, ('1', 'alias-a', '0.5', '2', 'group1'))),
                      dict(zip(fields, ('2', 'b', '5', '3', 'group2')))],
            'validation': [dict(zip(fields, ('3', 'alias-c', '0.5', '0', 'group3'))),
                           dict(zip(fields, ('4', 'd', '5', '3', 'group4')))],
            'test': [dict(zip(fields, ('5', 'e', '5', '1', 'group5')))],
        }
        parents = {'alias-a': 'A', 'b': 'B', 'alias-c': 'C', 'd': 'D', 'e': 'E'}
        target = {'train': {'B'}, 'validation': {'C'}, 'test': {'A', 'E'}}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'fold_0').mkdir()

            def write_role(role, cases):
                with (root / paths[f'{role}_cases']).open('w') as stream:
                    writer = csv.DictWriter(stream, fieldnames=fields)
                    writer.writeheader()
                    writer.writerows(cases)

            for role in ('train', 'validation'):
                write_role(role, original[role][1:])
            write_role('test', original['test'])
            counts = verify_roles(root, original, parents, target)
            self.assertEqual(counts['train']['retained']['wells'], 1)
            self.assertEqual(counts['train']['excluded']['groups'], 1)
            self.assertEqual(counts['validation']['retained_target_overlap']['validation'], 0)
            self.assertEqual(counts['test']['retained_target_overlap']['test'], 1)
            for role in ('train', 'validation'):
                with self.subTest(role=role):
                    write_role(role, original[role])
                    with self.assertRaisesRegex(ValueError, 'forbidden parent'):
                        verify_roles(root, original, parents, target)
                    write_role(role, original[role][1:])


if __name__ == '__main__':
    unittest.main()
