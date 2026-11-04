"""Regression at the real frozen-source adapter and public release boundary."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


EXECUTION = Path(__file__).resolve().parent.parent


class SirenkoSource(unittest.TestCase):
    def test_complete_source_preserves_blanks_and_zero_anomalies(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / 'sirenko.json'
            process = subprocess.run([sys.executable, str(EXECUTION/'sirenko_adapter.py'),
                                      '--output', str(output)], capture_output=True, text=True)
            self.assertEqual(process.returncode, 0, process.stderr)
            bundle = json.loads(output.read_text())
            audit = bundle['source_audits']['sirenko_2013']
            self.assertEqual(audit['compound_labels'], 131)
            self.assertEqual(audit['comparisons'], 524)
            self.assertEqual(audit['raw_at_least_one'], 524)
            self.assertEqual(audit['raw_complete_pairs'], 514)
            self.assertEqual(audit['raw_finite_values'], 1038)
            self.assertEqual(audit['raw_blank_cells'], 10)
            self.assertEqual(audit['raw_double_zero_comparisons'], 71)
            self.assertEqual([row['raw_selector'] for row in audit['raw_missing_normalized_zero']],
                             ['24hr_beat_count!F120', '24hr_beat_count!F121'])
            self.assertEqual({record['condition']['condition_vehicle_percent']
                              for record in bundle['records']}, {.1, .3})
            release = Path(temp)/'release'
            built = subprocess.run([sys.executable, str(EXECUTION/'corpus.py'), 'build',
                                    '--input', str(output), '--output', str(release)],
                                   capture_output=True, text=True)
            self.assertEqual(built.returncode, 0, built.stderr)
            coverage = json.loads((release/'coverage.json').read_text())
            self.assertEqual(coverage['sample_rows'], 0)
            self.assertEqual(coverage['excluded_comparisons'], 524)
            self.assertEqual(len(json.loads((release/'manifest.json').read_text())['partitions']), 2)
            comparisons = json.loads((release/'comparisons.json').read_text())
            self.assertTrue(all('unevaluable_fold_change' in row['exclusion_reasons'] for row in comparisons))


if __name__ == '__main__':
    unittest.main()
