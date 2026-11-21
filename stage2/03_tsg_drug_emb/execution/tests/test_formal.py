"""Formal receipt coordinates must preserve the paired experiment."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from formal import ARMS, validate_fit_grid


def fixture():
    return [
        {'arm': arm, 'fold': fold, 'order': {'train': [fold], 'validation': [fold + 3],
                                             'test': [fold + 6]},
         'transform': {'fixed': fold}, 'downstream_sha256': 'matched', 'rng_sha256': {'cpu': 'matched'}}
        for fold in range(3) for arm in ARMS
    ]


class GridTests(unittest.TestCase):
    def test_accepts_exactly_paired_complete_grid(self):
        validate_fit_grid(fixture())

    def test_rejects_missing_duplicate_or_individually_reassigned_fit(self):
        for replaced in ([], [fixture()[1]], [dict(fixture()[0], order={'train': [7]})],
                         [dict(fixture()[0], downstream_sha256='different')]):
            with self.subTest(replaced=replaced), self.assertRaises(ValueError):
                validate_fit_grid(replaced + fixture()[1:])


if __name__ == '__main__':
    unittest.main()
