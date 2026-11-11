"""Opt-in real-R release boundary for insufficient replication and axis-only inputs."""
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

import h5py
from anndata.io import write_elem
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


@unittest.skipUnless(os.environ.get('ARC_DESEQ2_TEST') == '1', 'explicit capped R integration job')
class RealWorkerRelease(unittest.TestCase):
    def test_raw_count_release_preserves_untestable_status(self):
        from bases.corpora import parent_structure
        from bases.significance import PartitionTruth
        from bases.truth_build import build, verify
        parent = parent_structure('CCO')[0]
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            native = directory / 'native.h5ad'
            training = directory / 'training'
            training.mkdir()
            obs = pd.DataFrame({
                'plate': ['plate'] * 3, 'is_control': [False, True, True],
                'is_blank': [False] * 3, 'compound_name': ['ethanol', 'DMSO', 'DMSO'],
                'vehicle': ['DMSO'] * 3, 'time_h': [24.] * 3,
                'drug_key': [parent, '', ''], 'dose_um': [1., 0., 0.],
                'smiles': ['CCO', 'CS(C)=O', 'CS(C)=O'],
            }, index=['well1', 'dmso1', 'dmso2'])
            genes = pd.DataFrame({'gene_index': [0, 1, 2], 'hgnc_id': ['HGNC:1','HGNC:2','HGNC:3']},
                                 index=['0','1','2'])
            with h5py.File(native, 'w') as file:
                write_elem(file, 'obs', obs)
                write_elem(file, 'var', genes.drop(columns='gene_index'))
                file.create_dataset('X', data=np.array([[10,20,30],[20,30,40],[30,40,50]], dtype='int32'))
            treated = pd.DataFrame({'sample_id': ['s1'], 'source_profile_id': ['well1'],
                'source_row': [0], 'drug_parent': [parent], 'drug_key': [parent], 'dose_um': [1.],
                'matched_control_rows': ['[1,2]'], 'matched_control_ids': ['["dmso1","dmso2"]']}, index=['s1'])
            path = training / 'tenv.h5ad'
            with h5py.File(path, 'w') as file:
                write_elem(file, 'obs', treated)
                write_elem(file, 'var', genes)
                write_elem(file, 'uns', {'schema': 'stage2-indexed-training-corpora-v1'})
                file.create_group('layers').create_dataset('control', data=np.zeros((1,3)))
                file['layers'].create_dataset('delta', data=np.zeros((1,3)))
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            (training / 'manifest.json').write_text(json.dumps({'schema':'stage2-indexed-training-corpora-v1',
                'fixture':True, 'files':{'tenv.h5ad':digest}}))
            output = directory / 'truth'
            rscript = '/home/chang/miniconda3/envs/biomni_e1/bin/Rscript'
            build(training, native, output, rscript, workers=1, batch_size=1, fixture=True)
            result = verify(training, native, output, rscript, replay_groups=0, fixture=True)
            self.assertEqual(result['statuses'], {'untestable': 1})
            assignment = pd.DataFrame({'sample_id':['s1'], 'role':['test']})
            labels = PartitionTruth(output, assignment, 'test', allow_fixture=True)
            self.assertIsNone(labels.truth(parent, 1, ['s1']))
            metadata = pd.read_csv(output/'groups.csv')
            self.assertEqual(metadata.iloc[0].reason, 'insufficient_replication')
            self.assertEqual(result['gene_count'], 3)


if __name__ == '__main__':
    unittest.main()
