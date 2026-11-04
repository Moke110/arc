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
        selected=copy.deepcopy(self.bundle if bundle is None else bundle)
        for record in selected['records']:
            key=lambda other:(other['source'],other['condition'],other['phenotype'],other['assay_method'],other.get('batch'))
            members=sorted(other['comparison_id'] for other in selected['records'] if key(other)==key(record))
            record['statistics'].setdefault('family_members',members)
        input_path=self.directory/(release+'-input.json')
        input_path.write_text(json.dumps(selected))
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

    def test_full_bh_family_includes_small_effect_before_label_selection(self):
        self.bundle['records']=[self.record('a',1.3,None),self.record('b',.7,None),self.record('c',1.1,None)]
        for record,p in zip(self.bundle['records'],(.01,.04,.20)):
            record['drug']='parent-'+record['comparison_id']
            record['statistics'].update(p_raw=p,correction='none')
        release=self.build()
        comparisons=json.loads((release/'comparisons.json').read_text())
        self.assertEqual([r['p_adjusted'] for r in comparisons],[.03,.06,.20])
        self.assertEqual([r['endpoint'] for r in comparisons],[1,0,0])
        families=json.loads((release/'families.json').read_text())
        self.assertEqual(families[0]['eligible_members'],['a','b','c'])

    def test_independent_paired_and_summary_routes_use_actual_experimental_units(self):
        self.bundle['records']=[self.record('a',2,None),self.record('b',2,None),self.record('c',2,None)]
        for record in self.bundle['records']:
            record['batch']=record['comparison_id']
            record['statistics'].update(correction='none',independent_units=True)
        self.bundle['records'][0]['statistics'].update(design='independent',
            treated_values=[1,2,3],control_values=[0,1,2],treated_unit_ids=['t1','t2','t3'],control_unit_ids=['c1','c2','c3'])
        self.bundle['records'][1]['statistics'].update(design='paired',
            treated_values=[2,4,7],control_values=[1,2,3],treated_unit_ids=['u1','u2','u3'],control_unit_ids=['u1','u2','u3'])
        self.bundle['records'][2]['statistics'].update(design='independent_summary',
            treated_summary={'mean':2,'sd':1,'n':3},control_summary={'mean':1,'sd':1,'n':3})
        release=self.build()
        comparisons=json.loads((release/'comparisons.json').read_text())
        self.assertAlmostEqual(comparisons[0]['p_raw'],.2878641347266908)
        self.assertAlmostEqual(comparisons[1]['p_raw'],.11808289631180317)
        self.assertAlmostEqual(comparisons[2]['p_raw'],.2878641347266908)
        self.assertEqual([r['statistical_origin'] for r in comparisons],['recomputed_Welch','recomputed_paired','recomputed_Welch_summary'])

    def test_endpoint_boundaries_missing_tests_and_invalid_denominators(self):
        cases=[('small',1.1,None,0),('large',1.3,None,None),('low',.8,.05,0),
               ('high',1.2,.05,0),('lowmissing',.8,None,None),('highmissing',1.2,None,None),
               ('zero',0.,.01,-1),('badcontrol',1.1,.01,None)]
        self.bundle['records']=[]
        for identity,fc,p,endpoint in cases:
            record=self.record(identity,fc,p)
            record['drug']='parent-'+identity
            if identity=='badcontrol':
                record['statistics']['control_mean']=0
            self.bundle['records'].append(record)
        release=self.build()
        comparisons={r['comparison_id']:r for r in json.loads((release/'comparisons.json').read_text())}
        for identity,_,_,endpoint in cases:
            self.assertEqual(comparisons[identity]['endpoint'],endpoint,identity)
        self.assertIsNone(comparisons['small']['p_adjusted'])
        self.assertIn('missing_required_adjusted_p',comparisons['large']['exclusion_reasons'])
        self.assertIn('unevaluable_fold_change',comparisons['badcontrol']['exclusion_reasons'])

    def test_technical_units_unknown_pairing_and_selective_p_cannot_create_significance(self):
        self.bundle['records']=[self.record('technical',1.3,None),self.record('paired',1.3,None),self.record('selective',1.3,None)]
        for record in self.bundle['records']:
            record['statistics'].update(correction='none',independent_units=True,design='independent',
                treated_values=[1,2,3],control_values=[0,1,2],treated_unit_ids=['t','t','t'],control_unit_ids=['c','c','c'])
        self.bundle['records'][1]['statistics'].update(design='paired',
            treated_unit_ids=['u1','u2','u3'],control_unit_ids=['u1','u2','other'])
        self.bundle['records'][2]['statistics'].update(p_raw=.001,family_complete=False)
        release=self.build()
        comparisons=json.loads((release/'comparisons.json').read_text())
        self.assertTrue(all(record['endpoint'] is None for record in comparisons))
        self.assertTrue(all(record['p_adjusted'] is None for record in comparisons))

    def test_cross_source_and_endpoint_stay_separate_and_canonical_times_share_partition(self):
        other=copy.deepcopy(self.bundle['sources'][0]);other['source']='fixture2'
        self.bundle['sources'].append(other)
        self.bundle['records'].append(self.record('c',.7,.01))
        self.bundle['records'].append(self.record('d',1.3,.01))
        self.bundle['records'][-1]['source']='fixture2'
        self.bundle['records'][1]['condition']['condition_drug_time']=24.0
        release=self.build()
        self.assertEqual(len(self.rows(release)),3)
        coverage=json.loads((release/'coverage.json').read_text())
        self.assertEqual(coverage['drug_dose_union'],1)
        self.assertEqual(len(list(release.glob('*.parquet'))),1)

    def test_declared_family_cannot_omit_small_effect_record(self):
        self.bundle['records']=[self.record('a',1.3,None),self.record('b',.7,None)]
        for record,p in zip(self.bundle['records'],(.01,.04)):
            record['statistics'].update(p_raw=p,correction='none',family_members=['a','b','c'])
        self.build(success=False)

    def test_coverage_separates_sources_partitions_and_unknown_independence(self):
        self.bundle['sources'].append({'source':'not-acquired','test_only':True,
                                      'status':'unavailable','assets':[],'reason':'bounded download failure'})
        self.bundle['records'][0]['measurements'][0]['experimental_unit']=None
        excluded=self.bundle['records'][1]
        excluded['drug']=None
        excluded['statistics']['control_matched']=False
        excluded['quality_reasons']=['pharmacology_unresolved']
        self.bundle['source_audits']={'fixture':{'all_doses':2,'control_gap':'not published'}}
        release=self.build()
        self.assertEqual(json.loads((release/'source_audits.json').read_text()),self.bundle['source_audits'])
        coverage=json.loads((release/'coverage.json').read_text())
        self.assertEqual(coverage['source_status_counts'],{'archived':1,'unavailable':1})
        self.assertEqual(coverage['archived_source_count'],1)
        self.assertEqual(coverage['unavailable_source_count'],1)
        source=coverage['per_source']['fixture']
        self.assertEqual(source['measurement_entries'],2)
        self.assertEqual(source['finite_source_values'],2)
        self.assertEqual(source['known_unit_count'],1)
        self.assertEqual(source['unknown_unit_measurements'],1)
        self.assertEqual(source['experimental_units_status'],'partially_known')
        self.assertIsNone(source['independent_experimental_units'])
        self.assertEqual(source['independence_unknown_entries'],2)
        self.assertEqual(source['identity_eligible_comparisons'],1)
        self.assertEqual(source['dose_eligible_comparisons'],2)
        self.assertEqual(source['fc_computable_comparisons'],1)
        self.assertEqual(source['admitted_comparisons'],1)
        self.assertEqual(source['excluded_comparisons'],1)
        self.assertEqual(source['class_rows'],{'-1':0,'0':0,'1':1})
        self.assertEqual(source['drug_dose_union'],1)
        self.assertEqual(coverage['partition_coverage'][0]['measurement_entries'],2)
        self.assertEqual(coverage['partition_coverage'][0]['excluded_comparisons'],1)
        self.assertEqual(coverage['per_source']['not-acquired']['assessment_status'],'not_run')
        self.assertEqual(source['exclusion_reasons_overlapping']['unresolved_drug_identity'],1)
        self.assertEqual(source['exclusion_reasons_overlapping']['unevaluable_fold_change'],1)

if __name__=='__main__':
    unittest.main()
