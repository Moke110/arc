import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd
from stage2.evaluation.evaluate import Labels, membership
from stage2.evaluation.inventory import REPO, sha


class LabelIntegrityTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(dir=REPO / 'stage2/evaluation/tests')
        self.root = Path(self.temporary.name)
        self.release = self.root / 'tenv'
        self.release.mkdir()
        pd.DataFrame({'hgnc_id': [f'HGNC:{i}' for i in range(18165)],
                      'symbol': [f'G{i}' for i in range(18164)] + ['HDGFL3']}).to_csv(self.release / 'genes.csv', index=False)
        self.groups = pd.DataFrame({'group_id': ['g'], 'status': ['eligible'],
                                    'case_ids_sha256': [membership(['a', 'b'])]})
        self.groups.to_csv(self.release / 'groups.csv', index=False)
        self.padj = np.ones((1, 18165)); self.padj[0, 0] = .01
        np.save(self.release / 'padj.npy', self.padj)
        self.direction = np.zeros_like(self.padj); self.direction[0, 0] = 1
        np.save(self.release / 'log2_fold_change.npy', self.direction)
        self.seal()

    def seal(self):
        manifest = {'schema': 'stage2-significance-v1', 'dataset': 'tenv', 'fdr': .05,
                    'files': {p.name: sha(p) for p in self.release.iterdir() if p.name != 'manifest.json'}}
        (self.release / 'manifest.json').write_text(json.dumps(manifest))

    def tearDown(self):
        self.temporary.cleanup()

    def test_payload_and_manifest_mutation_are_rejected_after_loading(self):
        labels = Labels('tenv', self.root)
        self.assertEqual(labels.truth('g', ['a', 'b'])[0], 1)
        self.padj[0, 0] = .5
        np.save(self.release / 'padj.npy', self.padj)
        with self.assertRaisesRegex(ValueError, 'changed'):
            labels.verify()
        self.seal()
        labels = Labels('tenv', self.root)
        path = self.release / 'manifest.json'
        path.write_text(path.read_text() + '\n')
        with self.assertRaisesRegex(ValueError, 'changed'):
            labels.verify()

    def test_contrast_membership_and_rehashed_eligibility_are_checked(self):
        labels = Labels('tenv', self.root)
        with self.assertRaisesRegex(ValueError, 'membership'):
            labels.truth('g', ['a'])
        self.groups['status'] = 'untestable'
        self.groups.to_csv(self.release / 'groups.csv', index=False)
        self.seal()
        with self.assertRaisesRegex(ValueError, 'eligibility'):
            Labels('tenv', self.root).truth('g', ['a', 'b'])

    def test_boundary_fdr_and_unmeasured_genes(self):
        self.padj[0, 0] = .05
        self.padj[0, 1] = np.nan
        np.save(self.release / 'padj.npy', self.padj)
        self.seal()
        self.assertEqual(np.count_nonzero(Labels('tenv', self.root).truth('g', ['a', 'b'])), 0)
