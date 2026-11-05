"""Behavior tests at the frozen-input CLI seam (ADR 0042)."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import h5py
import numpy as np
import pandas as pd
from anndata.io import write_elem

SCRIPT = Path(__file__).resolve().parent.parent / 'script' / 'inputs.py'


class InputReleaseTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.corpus = self.root / 'corpus.h5ad'
        obs = []
        counts = []
        for plate in ['p0', 'p1']:
            obs.append(dict(well=f'{plate}-control', plate=plate, drug_key='DMSO',
                            smiles='', dose_um=0., is_control=True, is_blank=False,
                            compound_name='DMSO', vehicle='DMSO'))
            counts.append([1, 3, 996])
            for i, smiles in enumerate(['CCO', 'CCN', 'CCC', 'CCCl', 'CCBr', 'CCF']):
                obs.append(dict(well=f'{plate}-{i}', plate=plate, drug_key=f'd{i}',
                                smiles=smiles, dose_um=10., is_control=False,
                                is_blank=False, compound_name=f'd{i}', vehicle='DMSO'))
                counts.append([3, 1, 996])
        with h5py.File(self.corpus, 'w') as f:
            f.create_dataset('X', data=np.array(counts, dtype=np.int32))
            write_elem(f, 'obs', pd.DataFrame(obs).set_index('well', drop=False))
            write_elem(f, 'var', pd.DataFrame({'hgnc_id': ['HGNC:1','HGNC:2','HGNC:3']}, index=['A','B','C']))
        self.genes = self.root / 'genes.csv'
        pd.DataFrame({'symbol':['B','A'], 'hgnc_id':['HGNC:2','HGNC:1']}).to_csv(self.genes,index=False)
        self.output = self.root / 'release'

    def run_cli(self, *args, ok=True):
        p = subprocess.run([sys.executable, str(SCRIPT), *map(str,args)], capture_output=True,text=True)
        if ok:
            self.assertEqual(p.returncode, 0, p.stdout+p.stderr)
        else:
            self.assertNotEqual(p.returncode,0)
        return p

    def build(self):
        self.run_cli('build','--corpus',self.corpus,'--genes',self.genes,'--output',self.output,
                     '--expected-genes','2','--expected-drugs','6')

    def test_cli_freezes_projected_normalization_and_safe_test_cases(self):
        self.build()
        self.run_cli('verify','--input',self.output)
        controls=np.load(self.output/'controls.npy',mmap_mode='r')
        np.testing.assert_allclose(controls[0],np.log1p([7500,2500]),rtol=1e-6)
        tested=[]
        for fold in range(3):
            base=self.output/f'fold_{fold}'
            cases=pd.read_csv(base/'test_cases.csv')
            self.assertEqual(set(cases.columns),{'case_id','drug_key','dose_um','control_index','group_id'})
            self.assertTrue(cases.case_id.is_unique)
            tested.extend(cases.drug_key.unique())
            truth=np.load(base/'evaluator/truth.npy')
            np.testing.assert_allclose(truth, np.tile(np.log1p([2500,7500])-controls[0],(len(truth),1)),rtol=1e-5)
            roles=[set(pd.read_csv(base/f'{role}.csv').drug_key) for role in ['train','validation','test_cases']]
            self.assertFalse(roles[0]&roles[1] or roles[0]&roles[2] or roles[1]&roles[2])
        self.assertEqual(sorted(tested),[f'd{i}' for i in range(6)])
        original=json.loads((self.output/'manifest.json').read_text())
        self.build()
        self.assertEqual(original,json.loads((self.output/'manifest.json').read_text()))

    def test_cli_rejects_malformed_control_axis_even_with_updated_hash(self):
        import hashlib
        self.build()
        np.save(self.output/'controls.npy', np.array([[float('nan')]],dtype=np.float32))
        path=self.output/'manifest.json'
        manifest=json.loads(path.read_text())
        manifest['files']['controls.npy']=hashlib.sha256((self.output/'controls.npy').read_bytes()).hexdigest()
        path.write_text(json.dumps(manifest))
        self.run_cli('verify','--input',self.output,ok=False)

    def test_cli_rejects_missing_dmso_and_tampered_cases(self):
        self.build()
        cases=self.output/'fold_0/test_cases.csv'
        frame=pd.read_csv(cases)
        frame['global_row_index']=range(len(frame))
        frame.to_csv(cases,index=False)
        self.run_cli('verify','--input',self.output,ok=False)
        with h5py.File(self.corpus,'r+') as f:
            del f['obs/compound_name']
            write_elem(f['obs'],'compound_name',np.array(['Media']*len(f['X']),dtype=object))
        self.run_cli('build','--corpus',self.corpus,'--genes',self.genes,'--output',self.output,
                     '--expected-genes','2','--expected-drugs','6',ok=False)

if __name__=='__main__':
    unittest.main()
