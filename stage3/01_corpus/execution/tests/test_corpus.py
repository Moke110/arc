"""Public offline Heliste build/verify acceptance seam (spec #52)."""
import csv
import json
import hashlib
import shutil
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import pyarrow.parquet as pq

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'.git').exists())
COMMAND = ROOT/'stage3/01_corpus/execution/corpus.py'

class CorpusWorkflow(unittest.TestCase):
    def command(self, *args, succeeds=True, script=COMMAND):
        result = subprocess.run([sys.executable, str(script), *map(str,args)], capture_output=True, text=True)
        if succeeds:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0)
        return result

    def test_source_fidelity_and_complete_accounting(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp)/'release'
            self.command('build','--output',output)
            self.command('verify','--output',output)
            for condition in ['baseline','H-R']:
                directory = output/condition
                rows = pq.read_table(directory/'corpus.parquet').to_pylist()
                with (directory/'corpus.csv').open() as stream:
                    csv_rows = list(csv.DictReader(stream))
                report = json.loads((directory/'chemical_resolution_report.json').read_text())
                self.assertEqual(report['source_label_rows'],689)
                self.assertEqual(report['source_values'],2756)
                self.assertEqual(len(rows)+report['pending_values'],2756)
                self.assertEqual(len(rows),len(csv_rows))
                self.assertGreater(len(rows),0)
                literal_cell = 'Sheet1!B5' if condition=='baseline' else 'Sheet1!F5'
                literal_value = 102.64832002399746 if condition=='baseline' else 103.74662895104134
                first = next(r for r in rows if r['source_cell']==literal_cell)
                self.assertEqual(first['source_label'],'Sulindac (Clinoril)')
                self.assertEqual(first['readout']['source_reported'],literal_value)
                for parquet_row,csv_row in zip(rows,csv_rows):
                    self.assertEqual(parquet_row['readout'],json.loads(csv_row.pop('readout')))
                    self.assertEqual(parquet_row['dose_um'],float(csv_row.pop('dose_um')))
                    self.assertEqual(csv_row.pop('rep'),'')
                    self.assertEqual({k:v for k,v in parquet_row.items() if k not in {'readout','dose_um','rep'}},csv_row)
                self.assertEqual(len({r['source_cell'] for r in rows}),len(rows))
                # Every resolved label has all four cells, even shared structures.
                self.assertEqual(len(rows),report['resolved_labels']*4)
                self.assertEqual({r['dose_um'] for r in rows},{0.04,0.2,1.0,5.0})
                self.assertTrue(all(r['rep'] is None and set(r['readout'])=={'source_reported'} for r in rows))
            self.command('build','--output',output)
            path = output/'baseline/corpus.csv'
            path.write_text(path.read_text().replace('source_reported','invented',1))
            self.command('verify','--output',output,succeeds=False)
            self.command('build','--output',output,succeeds=False)
            manifest_path = output/'manifest.json'
            manifest = json.loads(manifest_path.read_text())
            manifest['files']['baseline/corpus.csv'] = hashlib.sha256(path.read_bytes()).hexdigest()
            manifest_path.write_text(json.dumps(manifest))
            self.command('verify','--output',output,succeeds=False)

    def test_rehashed_false_identity_evidence_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            clone = Path(temp)/'repo'
            clone.mkdir()
            (clone/'.git').mkdir()
            shutil.copytree(ROOT/'pheno_data/Heliste-2021',clone/'pheno_data/Heliste-2021')
            shutil.copy(ROOT/'pheno_data/registry.json',clone/'pheno_data/registry.json')
            script = clone/'stage3/01_corpus/execution/corpus.py'
            script.parent.mkdir(parents=True)
            shutil.copy(COMMAND,script)
            evidence_path = clone/'pheno_data/Heliste-2021/local/chemical_resolution_v1.json'
            evidence = json.loads(evidence_path.read_text())
            record = next(r for r in evidence['records'] if r['status']=='resolved')
            attempt = next(a for a in record['attempts'] if a['query']==record['evidence']['query'])
            attempt['response']['PropertyTable']['Properties'][0]['CID'] = -1
            evidence_path.write_text(json.dumps(evidence))
            log_path = clone/'pheno_data/Heliste-2021/author_log.json'
            log = json.loads(log_path.read_text())
            addition = next(a for a in log['local_additions'] if a['path']=='local/chemical_resolution_v1.json')
            addition.update(bytes=evidence_path.stat().st_size,sha256=hashlib.sha256(evidence_path.read_bytes()).hexdigest())
            log_path.write_text(json.dumps(log))
            result = self.command('build','--output',clone/'release',script=script,succeeds=False)
            self.assertIn('Identity property evidence is not unique or differs',result.stderr)

if __name__ == '__main__':
    unittest.main()
