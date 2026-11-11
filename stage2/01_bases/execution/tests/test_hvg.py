"""Public release seam: held-out expression never selects model genes."""
import sys
from pathlib import Path
import tempfile
import unittest
import h5py
import numpy as np
import pandas as pd
from anndata.io import write_elem
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases.hvg import build_hvg, HVGRelease, verify_hvg, make_cohorts
from bases.hvg_qa import verify_population_isolation

class HVGReleaseTests(unittest.TestCase):
    def fixture(self, root, heldout_scale=1):
        path = root/'training.h5ad'
        linear = np.array([[100,100,0,9800],[200,200,0,9600],[400,400,0,9200],[800,800,0,8400]] * 3, dtype=np.float32)
        if heldout_scale != 1:
            linear[4:] = [0,0,10000,0]
        samples = pd.DataFrame({'sample_id':[f's{i}' for i in range(12)],'drug_parent':list('AAAABBBBCCCC'),'dose_um':[1.]*12})
        with h5py.File(path,'w') as f:
            f.create_dataset('X',data=np.log1p(linear))
            write_elem(f,'obs',samples.set_index('sample_id',drop=False))
            write_elem(f,'var',pd.DataFrame({'gene_index':range(4)},index=[str(i) for i in range(4)]))
        cohorts = samples.assign(fold=0,role=['train']*4+['validation']*4+['test']*4)
        return path, cohorts

    def test_release_train_only_selection_ignores_heldout_and_rejects_leakage(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); first=root/'first'; second=root/'second'; first.mkdir(); second.mkdir()
            training, cohorts=self.fixture(first)
            build_hvg(training,first/'hvg',cohorts=cohorts,expected_genes=4)
            release=HVGRelease(first/'hvg')
            rank=release.select('fold_0',4,training_sample_ids=['s0','s1','s2','s3'])
            self.assertEqual(rank.tolist(),[3,0,1,2])
            self.assertLess(rank.tolist().index(0),rank.tolist().index(1))
            with self.assertRaisesRegex(ValueError,'descriptive'):
                release.select('descriptive',2,training_sample_ids=['s0','s1'])
            with self.assertRaisesRegex(ValueError,'population'):
                release.select('fold_0',2,training_sample_ids=['s0','s4'])
            changed, _=self.fixture(second,heldout_scale=99)
            build_hvg(changed,second/'hvg',cohorts=cohorts,expected_genes=4)
            np.testing.assert_array_equal(rank,HVGRelease(second/'hvg').select('fold_0',4,training_sample_ids=['s0','s1','s2','s3']))
            self.assertEqual(verify_hvg(training,first/'hvg')['result'],'PASS')
            isolation=verify_population_isolation(training,first/'hvg')
            self.assertEqual(isolation['result'],'PASS')
            self.assertEqual(isolation['folds']['fold_0']['heldout_wells'],8)
            with self.assertRaisesRegex(ValueError,'seal'):
                verify_hvg(changed,first/'hvg')
            broken=cohorts.copy();broken.loc[1,'role']='test'
            with self.assertRaisesRegex(ValueError,'group'):
                build_hvg(training,first/'broken',cohorts=broken,expected_genes=4)

    def test_new_cold_parent_cohorts_are_deterministic_and_keep_all_doses(self):
        samples = pd.DataFrame([{'sample_id':f'{p}-{d}-{r}','drug_parent':p,'dose_um':float(d)}
                                for p in 'ABCDEF' for d in (1,10) for r in range(2)])
        cohorts=make_cohorts(samples,folds=3,seed='fixture-v1')
        self.assertEqual(cohorts.groupby(['fold','drug_parent']).role.nunique().max(),1)
        self.assertEqual(cohorts.groupby(['fold','role']).size().unstack().to_dict('records'),
                         [{'test':8,'train':12,'validation':4}]*3)
        changed=make_cohorts(samples.iloc[::-1],folds=3,seed='fixture-v1')
        pd.testing.assert_frame_equal(cohorts.sort_values(['fold','sample_id']).reset_index(drop=True),
                                      changed.sort_values(['fold','sample_id']).reset_index(drop=True))
