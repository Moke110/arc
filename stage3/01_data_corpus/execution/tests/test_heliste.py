"""Real-source regression through the public adapter and corpus workflow."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'.git').exists())
EXECUTION = ROOT/'stage3/01_data_corpus/execution'

class HelisteWorkflow(unittest.TestCase):
    def test_full_source_values_are_retained_and_unresolved_comparisons_excluded(self):
        with tempfile.TemporaryDirectory() as tmp:
            bundle, release = Path(tmp)/'heliste.json', Path(tmp)/'release'
            subprocess.run([sys.executable,str(EXECUTION/'heliste_adapter.py'),'--output',str(bundle)],check=True)
            data = json.loads(bundle.read_text())
            audit = data['source_audits']['heliste_2021']
            self.assertEqual(audit['distinct_exact_labels'],689)
            self.assertEqual(audit['numeric_values_by_oxygen'],{'CTRL':2756,'HYP':2756})
            self.assertEqual(audit['dataset_s3_hit_labels'],48)
            self.assertEqual(len(data['records']),5512)
            self.assertEqual({r['dose'] for r in data['records']},{4e-8,2e-7,1e-6,5e-6})
            self.assertEqual({r['phenotype'] for r in data['records']},{'ctg_atp_viability'})
            self.assertEqual({m['representation'] for r in data['records'] for m in r['measurements']},{'normalized_only'})
            subprocess.run([sys.executable,str(EXECUTION/'corpus.py'),'build','--input',str(bundle),'--output',str(release)],check=True)
            subprocess.run([sys.executable,str(EXECUTION/'corpus.py'),'verify','--release',str(release)],check=True)
            manifest = json.loads((release/'manifest.json').read_text())
            self.assertEqual(len(manifest['partitions']),2)
            self.assertTrue(all(p['sample_rows']==0 for p in manifest['partitions']))

if __name__ == '__main__':
    unittest.main()
