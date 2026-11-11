"""Exact successor binding at the public envelope/reader seam; no R fitting."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

import h5py
from anndata.io import write_elem
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases.corpora import NORMALIZATION, parent_structure
from bases.significance import contrast_plan, publish_statistics, sha256, write_json
from bases.truth_build import WORKER, verify, write_part_inputs
from bases.truth_successor import create_successor, SuccessorTruth, verify_successor


class SuccessorContract(unittest.TestCase):
    def fixture(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        parent = parent_structure('CCO')[0]
        genes = pd.DataFrame({'gene_index':[0,1,2], 'hgnc_id':['HGNC:1','HGNC:2','HGNC:3']}, index=['0','1','2'])
        raw_obs = pd.DataFrame({'plate':['P']*4, 'is_control':[True,True,False,False],
            'is_blank':[False]*4, 'compound_name':['DMSO','DMSO','ethanol','ethanol'],
            'drug_key':['','',parent,parent], 'dose_um':[0.,0.,1.,1.],
            'vehicle':['DMSO']*4, 'time_h':[24.]*4, 'smiles':['','','CCO','CCO']},
            index=['c1','c2','t1','t2'])
        counts = np.array([[1,3,0],[3,1,0],[2,2,0],[1,3,0]], dtype='int32')
        normalized = np.log1p(np.array([[2500.,7500.,0],[7500.,2500.,0],[5000.,5000.,0],[2500.,7500.,0]])).astype('float32')
        control = normalized[:2].mean(axis=0, dtype='float64').astype('float32')
        versions = []
        for number, order in enumerate(([0,1,2,3],[3,2,1,0])):
            directory = root / f'v{number+1}'
            directory.mkdir()
            native = directory / 'tenv.h5ad'
            with h5py.File(native, 'w') as file:
                write_elem(file, 'obs', raw_obs.iloc[order])
                write_elem(file, 'var', genes.drop(columns='gene_index'))
                file.create_dataset('X', data=counts[order])
            write_json(directory/'tenv_manifest.json', {'schema':'fixture-native', 'files':{'tenv.h5ad':sha256(native)}})
            training = directory / 'training'
            training.mkdir()
            source = {name:index for index,name in enumerate(raw_obs.index[order])}
            names = ['t1','t2'] if number == 0 else ['t2','t1']
            obs = pd.DataFrame({'sample_id':['s-'+n for n in names], 'source_profile_id':names,
                'source_row':[source[n] for n in names], 'drug_parent':[parent]*2,
                'drug_key':[parent]*2, 'dose_um':[1.]*2,
                'matched_control_rows':[json.dumps([source['c1'],source['c2']])]*2,
                'matched_control_ids':['["c1","c2"]']*2}, index=['s-'+n for n in names])
            treated_positions = [list(raw_obs.index).index(n) for n in names]
            with h5py.File(training/'tenv.h5ad', 'w') as file:
                write_elem(file, 'obs', obs)
                write_elem(file, 'var', genes)
                write_elem(file, 'uns', {'schema':'stage2-indexed-training-corpora-v1'})
                file.create_dataset('X', data=normalized[treated_positions])
                file.create_group('layers').create_dataset('control', data=np.tile(control,(2,1)))
                file['layers'].create_dataset('delta', data=normalized[treated_positions]-control)
            self.seal_training(training)
            versions.append((native,training,obs))
        old_native, old_training, old_obs = versions[0]
        groups,samples = contrast_plan(old_obs,raw_obs)
        workspace = root/'fit-input'
        groups = write_part_inputs(workspace,old_native,groups,samples,np.arange(3))
        groups['status']='eligible';groups['reason']='';groups['n_significant']=1
        stats={'pvalue':np.array([[.001,.5,np.nan]]), 'padj':np.array([[.003,.75,np.nan]]),
               'log2_fold_change':np.array([[2.,-1.,np.nan]])}
        provenance={'native_counts_sha256':sha256(old_native),
            'training_manifest_sha256':sha256(old_training/'manifest.json'),
            'tenv_training_sha256':sha256(old_training/'tenv.h5ad'),
            'worker_sha256':sha256(WORKER), 'r_version':'fixture-R', 'deseq2_version':'fixture-DESeq2',
            'replication_rule':'at least 2 treated and 2 matched-control wells; residual design degrees of freedom'}
        truth=root/'truth'
        publish_statistics(truth,genes.reset_index(drop=True),groups,samples,stats,provenance,fixture=True)
        accepted=root/'accepted.json'
        write_json(accepted,verify(old_training,old_native,truth,'unused',replay_groups=0,fixture=True))
        kwargs={'truth':truth,'original_training':old_training,'training':versions[1][1],
                'original_native':old_native,'native':versions[1][0],'accepted_verification':accepted,
                'output':root/'successor','fixture':True,
                'runtime':{'r_version':'fixture-R','deseq2_version':'fixture-DESeq2','rscript_sha256':'fixture'}}
        return kwargs,parent,versions

    def seal_training(self, training):
        write_json(training/'manifest.json', {'schema':'stage2-indexed-training-corpora-v1',
            'normalization':NORMALIZATION,'fixture':True,'files':{'tenv.h5ad':sha256(training/'tenv.h5ad')}})

    def test_repacking_rows_preserves_exact_counts_layers_and_partition_truth(self):
        args,parent,_=self.fixture()
        result=create_successor(**args)
        self.assertEqual(result['result'],'PASS')
        replay=verify_successor(args['output'],args['training'],allow_fixture=True,runtime=args['runtime'])
        self.assertEqual(replay['result'],'PASS')
        assignment=pd.DataFrame({'sample_id':['s-t1','s-t2'],'role':['test','test']})
        bound=SuccessorTruth(args['output'],args['training'],assignment,'test',allow_fixture=True)
        self.assertEqual(bound.truth(parent,1,['s-t2','s-t1']).tolist(),[1,0,0])
        envelope=json.loads((args['output']/'manifest.json').read_text())
        self.assertEqual(envelope['proof']['raw_contrasts_compared'],1)
        self.assertEqual(set(envelope['proof']['exact_layer_digests']),{'X','control','delta'})

    def test_one_ulp_expression_difference_is_not_equivalence(self):
        args,_,_=self.fixture()
        with h5py.File(args['training']/'tenv.h5ad','r+') as file:
            value=file['layers/delta'][0,0]
            file['layers/delta'][0,0]=np.nextafter(value,np.float32(np.inf))
        self.seal_training(args['training'])
        with self.assertRaisesRegex(ValueError,'exact|layer'):
            create_successor(**args)

    def test_identical_normalized_profiles_do_not_hide_different_integer_counts(self):
        args,_,_=self.fixture()
        with h5py.File(args['native'],'r+') as file:
            file['X'][:]=file['X'][:]*2
        write_json(args['native'].parent/'tenv_manifest.json', {'schema':'fixture-native','files':{'tenv.h5ad':sha256(args['native'])}})
        with self.assertRaisesRegex(ValueError,'count'):
            create_successor(**args)

    def test_changed_recipe_is_not_reuse(self):
        args,_,_=self.fixture()
        args['config']={'fdr':.01}
        with self.assertRaisesRegex(ValueError,'recipe'):
            create_successor(**args)

    def test_runtime_and_worker_seals_are_not_caller_relabeling(self):
        args,_,_=self.fixture()
        args['runtime']={**args['runtime'],'deseq2_version':'changed'}
        with self.assertRaisesRegex(ValueError,'runtime'):
            create_successor(**args)
        args,_,_=self.fixture()
        altered=args['output'].parent/'altered-worker.R'
        altered.write_bytes(WORKER.read_bytes()+b'\n# a different worker seal\n')
        args['worker']=altered
        with self.assertRaisesRegex(ValueError,'worker'):
            create_successor(**args)

    def test_incomplete_control_membership_is_not_an_envelope(self):
        args,_,_=self.fixture()
        from anndata.io import read_elem
        with h5py.File(args['training']/'tenv.h5ad','r+') as file:
            obs=read_elem(file['obs'])
            obs['matched_control_rows']='[3]'
            obs['matched_control_ids']='["c1"]'
            write_elem(file,'obs',obs)
        self.seal_training(args['training'])
        with self.assertRaisesRegex(ValueError,'control|DMSO'):
            create_successor(**args)

    def test_partial_and_fixture_reports_cannot_publish_production_truth(self):
        args,_,_=self.fixture()
        report=json.loads(args['accepted_verification'].read_text())
        report['groups_replayed']=0
        write_json(args['accepted_verification'],report)
        with self.assertRaisesRegex(ValueError,'complete accepted'):
            create_successor(**args)
        args,_,_=self.fixture()
        args['fixture']=False
        args['runtime']=None
        with self.assertRaisesRegex(ValueError,'complete accepted'):
            create_successor(**args)


class RepresentativeReplayContract(unittest.TestCase):
    def test_three_refits_cover_ordinary_multiplate_and_tested_empty(self):
        from bases.truth_build import representative_indices
        groups=pd.DataFrame({'status':['eligible']*4,'n_plates':[1,2,2,1],'n_significant':[3,2,1,0]})
        self.assertEqual(representative_indices(groups,3),[0,1,3])


if __name__=='__main__': unittest.main()
