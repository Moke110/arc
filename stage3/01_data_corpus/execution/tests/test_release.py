"""All configured real sources through the highest public frozen-input seam."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

EXECUTION = Path(__file__).resolve().parent.parent


class CompleteRelease(unittest.TestCase):
    def test_all_five_sources_are_accounted_for_without_invented_labels(self):
        with tempfile.TemporaryDirectory() as temp:
            bundle = Path(temp)/'input.json'
            run = subprocess.run([sys.executable,str(EXECUTION/'release.py'),'bundle','--output',str(bundle)],
                                 capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stderr)
            data = json.loads(bundle.read_text())
            self.assertEqual(len(data['sources']),5)
            self.assertEqual(len(data['records']),6096)
            release = Path(temp)/'release'
            for args in [('build','--input',str(bundle),'--output',str(release)),
                         ('verify','--release',str(release))]:
                result = subprocess.run([sys.executable,str(EXECUTION/'corpus.py'),*args],
                                        capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stderr)
            coverage = json.loads((release/'coverage.json').read_text())
            self.assertEqual(coverage['sample_rows'],0)
            self.assertEqual(coverage['class_rows'],{'-1':0,'0':0,'1':0})
            self.assertEqual(coverage['excluded_comparisons'],6096)
            self.assertEqual(len(json.loads((release/'manifest.json').read_text())['partitions']),9)
            self.assertEqual(set(json.loads((release/'source_audits.json').read_text())),
                             {'sirenko_2013','heliste_2021','linn_2024','liu_2021','duan_2013'})


if __name__ == '__main__':
    unittest.main()
