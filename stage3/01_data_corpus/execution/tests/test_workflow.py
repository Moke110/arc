"""Acceptance at the public frozen-input build and offline verify boundary."""
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pyarrow.parquet as pq

SCRIPT = Path(__file__).resolve().parent.parent / 'corpus.py'

class Workflow(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.directory = Path(self.temp.name)
        self.asset = self.directory/'source.txt'
        self.asset.write_text('test-only known evidence\n')
        self.bundle = {'schema_version':1,'sources':[{'source':'fixture','test_only':True,
            'assets':[{'path':str(self.asset),'sha256':hashlib.sha256(self.asset.read_bytes()).hexdigest()}]}],
            'records':[self.record('a',1.3,.01),self.record('b',1.3,.01)]}
    def tearDown(self):
        self.temp.cleanup()
    def record(self, identity, fc, p):
        return {'comparison_id':identity,'source':'fixture','condition':{
            'condition_cell_type':'hiPSC-CM','condition_drug_time':24,'condition_hypoxia_time':0,
            'condition_reoxygen_time':0,'condition_drug_schedule':'direct'},
            'phenotype':'cell_viability','assay_method':'MTT','batch':'campaign1','drug':'evidence-parent-A','dose':1e-6,
            'identity_evidence':'fixture parent mapping','dose_evidence':'1 uM = 1e-6 mol/L',
            'control_evidence':'same experimental condition and vehicle; normalized control known',
            'measurements':[{'measurement_id':'m-'+identity,'value':13 if identity=='a' else 14,
                'assay_method':'MTT','unit':'OD','representation':'raw','asset':str(self.asset),
                'selector':'fixture '+identity,'experimental_unit':identity,'normalized_value':130}],
            'statistics':{'treated_mean':fc*10,'control_mean':10,'scale_compatible':True,'control_matched':True,
                'method':'Welch','sidedness':'two-sided','p_adjusted':p,'correction':'BH','family_matches':True,
                'family_complete':True,'evidence':'test-only complete family'}}
    def cli(self, *args, success=True):
        process=subprocess.run([sys.executable,str(SCRIPT),*args],capture_output=True,text=True)
        if success:
            self.assertEqual(process.returncode,0,process.stderr)
        else:
            self.assertNotEqual(process.returncode,0,process.stdout)
        return process
    def build(self,bundle=None,release='release',success=True,test_only=True):
        input_path=self.directory/(release+'-input.json')
        input_path.write_text(json.dumps(self.bundle if bundle is None else bundle))
        self.cli('build','--input',str(input_path),'--output',str(self.directory/release),
                 *(['--test-only'] if test_only else []),success=success)
        return self.directory/release
    def rows(self,release):
        return [row for path in release.glob('*.parquet') for row in pq.read_table(path).to_pylist()]
    def test_frozen_evidence_becomes_one_native_list_and_cannot_overwrite_release(self):
        release=self.build()
        rows=self.rows(release)
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0]['readout'],[13.,14.])
        self.assertEqual(rows[0]['endpoint'],1)
        self.assertEqual([m['measurement_id'] for m in rows[0]['annotations']],['m-a','m-b'])
        self.cli('verify','--release',str(release),'--test-only')
        self.build()
        changed=copy.deepcopy(self.bundle)
        changed['records'][0]['measurements'][0]['value']=15
        self.build(changed,success=False)
        self.cli('verify','--release',str(release),'--test-only')

    def test_test_sources_are_rejected_in_production_and_missing_evidence_fails(self):
        self.build(success=False,test_only=False)
        malformed=copy.deepcopy(self.bundle)
        malformed['records'][0]['measurements'][0]['asset']='absent.xlsx'
        self.build(malformed,release='broken',success=False)

    def test_empty_after_qc_preserves_blanks_and_unknown_timing(self):
        self.bundle['records']=self.bundle['records'][:1]
        record=self.bundle['records'][0]
        record['statistics']['control_matched']=False
        record['measurements'][0]['value']=None
        record['measurements'][0]['normalized_value']=0
        record['condition']['condition_hypoxia_time']=None
        release=self.build()
        self.assertEqual(self.rows(release),[])
        manifest=json.loads((release/'manifest.json').read_text())
        self.assertEqual(manifest['data_status'],'empty_after_QC')
        self.assertIsNone(manifest['partitions'][0]['condition']['condition_hypoxia_time'])
        comparisons=json.loads((release/'comparisons.json').read_text())
        self.assertIsNone(comparisons[0]['measurements'][0]['value'])
        self.assertEqual(comparisons[0]['measurements'][0]['normalized_value'],0)
        self.assertIn('unevaluable_fold_change',comparisons[0]['exclusion_reasons'])

    def test_corrupted_release_is_rejected(self):
        release=self.build()
        (release/'coverage.json').write_text('{}')
        self.cli('verify','--release',str(release),'--test-only',success=False)

if __name__=='__main__':
    unittest.main()
