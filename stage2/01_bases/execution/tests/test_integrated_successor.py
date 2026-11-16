"""ADR 0072 successor contracts at the public product boundary."""
import json
import sys
import subprocess
import tempfile
import unittest
from pathlib import Path

EXECUTION = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXECUTION))
import public_api
from bases import gene_axis


class SuccessorContractTests(unittest.TestCase):
    def test_declared_roster_retains_af_and_replaces_class_terms_with_one_tox(self):
        self.assertEqual(gene_axis.CATEGORIES,
                         ('nhf', 'icm', 'af', 'genetic_cardiomyopathy', 'tox'))
        self.assertEqual(gene_axis.CATEGORY_LABELS['tox'], 'TOX')

    def test_publish_refuses_to_replace_an_existing_candidate(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sentinel = root / 'retained.txt'
            sentinel.write_text('authoritative bytes')
            result = subprocess.run(
                [sys.executable, str(EXECUTION / 'build_integrated_successor.py'),
                 'publish', '--candidate', str(root), '--inputs', str(root),
                 '--tox-candidate', str(root)], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('candidate exists', result.stderr)
            self.assertEqual(sentinel.read_text(), 'authoritative bytes')

    def test_reader_refuses_retired_schema_instead_of_adapting_legacy_roster(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / 'manifest.json').write_text(json.dumps({'schema': 'stage2-gene-axis-release-v1'}))
            with self.assertRaisesRegex(ValueError, 'schema'):
                gene_axis.GeneAxis(root).candidate_table()


if __name__ == '__main__':
    unittest.main()
