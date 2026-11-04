"""Real-source regression through the public curation/build/verify workflow."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'.git').exists())
EXEC = ROOT/'stage3/01_data_corpus/execution'

class ViabilitySources(unittest.TestCase):
    def test_linn_full_six_curves_retain_uncertain_doses(self):
        with tempfile.TemporaryDirectory() as temp:
            bundle = Path(temp)/'bundle.json'
            result = subprocess.run([sys.executable,str(EXEC/'viability_adapters.py'),'--source','linn_2024','--output',str(bundle)],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            data=json.loads(bundle.read_text())
            self.assertEqual(len(data['records']),48)
            self.assertEqual(len({r['condition']['condition_cell_type'] for r in data['records']}),2)
            self.assertTrue(all(r['dose'] is None for r in data['records']))
            self.assertTrue(all(r['quality_reasons']==['exact_source_dose_unresolved'] for r in data['records']))
            self.assertEqual(sum(r['drug']=='remdesivir' for r in data['records']),16)
            release=Path(temp)/'release'
            for args in [('build','--input',str(bundle),'--output',str(release)),('verify','--release',str(release))]:
                run=subprocess.run([sys.executable,str(EXEC/'corpus.py'),*args],capture_output=True,text=True)
                self.assertEqual(run.returncode,0,run.stderr)
            coverage=json.loads((release/'coverage.json').read_text())
            self.assertEqual(coverage['excluded_comparisons'],48)
            self.assertEqual(coverage['sample_rows'],0)

if __name__=='__main__': unittest.main()
