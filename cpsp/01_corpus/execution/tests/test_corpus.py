"""Public offline Heliste build/verify acceptance seam (spec #52)."""
import csv
from collections import Counter
import importlib.util
import math
import json
import hashlib
import shutil
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'.git').exists())
COMMAND = ROOT/'cpsp/01_corpus/execution/corpus.py'

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
                    self.assertEqual(parquet_row['direction'],int(csv_row.pop('direction')))
                    self.assertEqual(csv_row.pop('rep'),'')
                    self.assertEqual({k:v for k,v in parquet_row.items() if k not in {'readout','dose_um','rep','direction'}},csv_row)
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

    def test_magnitude_direction_boundaries(self):
        spec = importlib.util.spec_from_file_location('heliste_corpus', COMMAND)
        corpus = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(corpus)
        for value,expected in [(math.nextafter(80,-math.inf),-1),(80,0),
                               (100,0),(120,0),(math.nextafter(120,math.inf),1)]:
            with self.subTest(value=value):
                self.assertEqual(corpus.direction_from_readout(value),expected)

    def test_direction_exports_and_rehashed_false_label_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp)/'release'
            self.command('build','--output',output)
            self.command('verify','--output',output)
            manifest_path = output/'manifest.json'
            manifest = json.loads(manifest_path.read_text())
            self.assertEqual(manifest['schema_version'],2)
            for condition,expected in [('baseline',{-1:194,0:2490,1:60}),
                                       ('H-R',{-1:178,0:2539,1:27})]:
                table = pq.read_table(output/condition/'corpus.parquet')
                self.assertEqual(table.schema.field('direction').type,pa.int8())
                self.assertEqual(Counter(table.column('direction').to_pylist()),expected)
                annotation = json.loads((output/condition/'source_annotation.json').read_text())
                policy = annotation['direction_policy']
                self.assertEqual(policy['version'],1)
                self.assertEqual(policy['reference_value'],100)
                self.assertEqual(policy['threshold_percent'],20)
                self.assertEqual(policy['zero_interval'],[80,120])
                self.assertTrue(policy['zero_interval_inclusive'])
                self.assertEqual(policy['normalization_status'],'working_assumption')
                self.assertFalse(policy['statistical_significance_claim'])
            # Change both exports consistently and recompute their checksums:
            # verification must still reconstruct the correct source-derived label.
            directory = output/'baseline'
            table = pq.read_table(directory/'corpus.parquet')
            values = table.to_pylist()
            values[0]['direction'] = 1
            pq.write_table(pa.Table.from_pylist(values,schema=table.schema),directory/'corpus.parquet')
            csv_path = directory/'corpus.csv'
            with csv_path.open() as stream:
                reader = csv.DictReader(stream)
                fields = reader.fieldnames
                rows = list(reader)
            rows[0]['direction'] = '1'
            with csv_path.open('w',newline='') as stream:
                writer = csv.DictWriter(stream,fieldnames=fields,lineterminator='\n')
                writer.writeheader()
                writer.writerows(rows)
            for name in ['baseline/corpus.csv','baseline/corpus.parquet']:
                manifest['files'][name] = hashlib.sha256((output/name).read_bytes()).hexdigest()
            manifest_path.write_text(json.dumps(manifest))
            result = self.command('verify','--output',output,succeeds=False)
            self.assertIn('release differs from offline rebuild',result.stderr)

    def test_desalted_parents_preserve_source_rows_and_permanent_charge(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp)/'release'
            self.command('build','--output',output)
            self.command('verify','--output',output)
            rows = pq.read_table(output/'baseline/corpus.parquet').to_pylist()
            coordination = {'Gadodiamide (Omniscan)','Sodium orthovanadate','Bismuth Subsalicylate','Vitamin B12','Cisplatin','Carboplatin','Oxaliplatin (Eloxatin)','Sodium Nitroprusside'}
            self.assertTrue(all('.' not in r['canonical_smiles'] for r in rows if r['source_label'] not in coordination))
            for label,metal in [('Vitamin B12','Co'),('Oxaliplatin (Eloxatin)','Pt'),('Sodium Nitroprusside','Fe')]:
                drug = next(r for r in rows if r['source_label']==label)
                self.assertIn(metal,drug['canonical_smiles'])
                self.assertNotIn('Na',drug['canonical_smiles'])
            amiloride = [r for r in rows if r['source_label'] in {'Amiloride hydrochloride (Midamor)','Amiloride hydrochloride dihydrate'}]
            self.assertEqual(len(amiloride),8)
            self.assertEqual(len({r['canonical_smiles'] for r in amiloride}),1)
            self.assertEqual(amiloride[0]['canonical_smiles'],'NC(N)=NC(=O)c1nc(Cl)c(N)nc1N')
            meth = [r for r in rows if r['source_label']=='Methscopolamine (Pamine)']
            self.assertEqual(len(meth),4)
            self.assertIn('[N+]',meth[0]['canonical_smiles'])
            self.assertNotIn('Br',meth[0]['canonical_smiles'])
            self.assertEqual(len([r for r in rows if r['source_label']=='Quinine hydrochloride dihydrate']),4)
            self.assertEqual(len([r for r in rows if r['source_label']=='Formoterol hemifumarate']),0)
            vitamin = next(r for r in rows if r['source_label']=='Vitamin D2')
            self.assertNotIn('[2H]',vitamin['canonical_smiles'])

    def test_rehashed_false_identity_evidence_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            clone = Path(temp)/'repo'
            clone.mkdir()
            (clone/'.git').mkdir()
            shutil.copytree(ROOT/'pheno_data/Heliste-2021',clone/'pheno_data/Heliste-2021')
            shutil.copy(ROOT/'pheno_data/registry.json',clone/'pheno_data/registry.json')
            script = clone/'cpsp/01_corpus/execution/corpus.py'
            script.parent.mkdir(parents=True)
            shutil.copy(COMMAND,script)
            evidence_path = clone/'pheno_data/Heliste-2021/local/chemical_resolution_v2.json'
            evidence = json.loads(evidence_path.read_text())
            record = next(r for r in evidence['records'] if r['status']=='resolved')
            attempt = next(a for a in record['original_reference']['attempts'] if a['query']==record['original_reference']['evidence']['query'])
            attempt['response']['PropertyTable']['Properties'][0]['CID'] = -1
            evidence_path.write_text(json.dumps(evidence))
            log_path = clone/'pheno_data/Heliste-2021/author_log.json'
            log = json.loads(log_path.read_text())
            addition = next(a for a in log['local_additions'] if a['path']=='local/chemical_resolution_v2.json')
            addition.update(bytes=evidence_path.stat().st_size,sha256=hashlib.sha256(evidence_path.read_bytes()).hexdigest())
            log_path.write_text(json.dumps(log))
            result = self.command('build','--output',clone/'release',script=script,succeeds=False)
            self.assertIn('Identity property evidence is not unique or differs',result.stderr)
            evidence = json.loads((ROOT/'pheno_data/Heliste-2021/local/chemical_resolution_v2.json').read_text())
            quinine = next(r for r in evidence['records'] if r['source_label']=='Quinine hydrochloride dihydrate')
            quinine['parent_attempts'].pop(0)
            evidence_path.write_text(json.dumps(evidence))
            addition.update(bytes=evidence_path.stat().st_size,sha256=hashlib.sha256(evidence_path.read_bytes()).hexdigest())
            log_path.write_text(json.dumps(log))
            result = self.command('build','--output',clone/'release',script=script,succeeds=False)
            self.assertIn('Missing candidate synonym evidence',result.stderr)


if __name__ == '__main__':
    unittest.main()
