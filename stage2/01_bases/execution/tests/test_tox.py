"""Unified TOX candidate publication through the public product/reader seam."""
import sys
import tempfile
import unittest
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import public_api  # noqa: F401
from bases import disease_deg, tox


def axis_fixture():
    return pd.DataFrame({'gene_index': np.arange(18164),
                         'symbol': [f'G{i}' for i in range(18164)],
                         'hgnc_id': [f'HGNC:{i}' for i in range(18164)],
                         'ensembl_id': [f'ENSG{i}' for i in range(18164)]})


def source(source_id, genes, tier=2):
    return {'source': source_id, 'tier': tier, 'item': 'fixture/' + source_id,
            'contrasts': pd.DataFrame({
                'gene_symbol': genes, 'contrast': 'TOX|matched',
                'log2fc': [3. - i * .2 for i in range(len(genes))], 'padj': [0.001] * len(genes),
                'p_value': [0.0001] * len(genes), 'stat': [8.] * len(genes),
                'base_mean': [10.] * len(genes)}),
            'contrast_records': [{'contrast': 'TOX|matched', 'stratum': 'matched',
                                  'n_case': 2, 'n_control': 2}],
            'method': {'name': 'fixture known contrasts'},
            'audit': [{'sample': 'drug1', 'role': 'treated', 'source_id': 'D1',
                       'stratum': 'matched', 'matched_controls': ['v1', 'v2']},
                      {'sample': 'drug2', 'role': 'treated', 'source_id': 'D2',
                       'stratum': 'matched', 'matched_controls': ['v1', 'v2']},
                      {'sample': 'v1', 'role': 'vehicle', 'source_id': 'vehicle',
                       'stratum': 'matched', 'matched_controls': ['v1', 'v2']},
                      {'sample': 'v2', 'role': 'vehicle', 'source_id': 'vehicle',
                       'stratum': 'matched', 'matched_controls': ['v1', 'v2']}],
            'member_drugs': ['D1', 'D2'], 'member_source_ids': ['D1', 'D2']}


class CandidateTest(unittest.TestCase):
    def test_candidate_reader_rejects_membership_changed_after_open(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            tox.assemble_candidate(root, [source('one', ['G1', 'G2'])],
                                   axis_fixture(), {}, 'f' * 64)
            reader = tox.ToxCandidate(root)
            (root / 'tox/one.audit.json').write_text('[]\n')
            with self.assertRaisesRegex(tox.ReleaseError, 'digest'):
                reader.audit('one')

    def test_a_reference_drug_cannot_substitute_for_matched_vehicle(self):
        with tempfile.TemporaryDirectory() as temporary:
            record = source('one', ['G1', 'G2'])
            record['audit'][0]['matched_controls'] = ['drug2']
            with self.assertRaisesRegex(tox.ReleaseError, 'vehicle'):
                tox.assemble_candidate(Path(temporary), [record], axis_fixture(), {}, 'f' * 64)

    def test_original_source_membership_pools_all_drugs_at_each_matched_time(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / tox.cardiotox.GSE243674_COUNTS
            path.parent.mkdir(parents=True)
            columns = ['VEH.1.3h', 'DOX.1.3h', 'EPI.1.3h', 'TRZ.1.3h',
                       'VEH.1.24h', 'DOX.1.24h', 'EPI.1.24h']
            pd.DataFrame([[1] * len(columns)], index=['101'], columns=columns).to_csv(
                path, sep=' ', compression='gzip')
            members = pd.DataFrame([dict(source_id=n, parent_inchikey=n, drug_name=n)
                                    for n in ['DOX', 'EPI']])
            audit = tox.legacy_membership(root, 'gse243674', members)
            selected = [r for r in audit if r['role'] == 'treated']
            self.assertEqual({r['sample'] for r in selected},
                             {'DOX.1.3h', 'EPI.1.3h', 'DOX.1.24h', 'EPI.1.24h'})
            self.assertEqual(next(r for r in selected if r['sample'] == 'EPI.1.24h')
                             ['matched_controls'], ['VEH.1.24h'])
            self.assertEqual(next(r for r in audit if r['sample'] == 'TRZ.1.3h')['role'],
                             'excluded')

    def test_plan_uses_only_positive_admitted_drugs_and_one_independent_dataset(self):
        ledger = {'sources': [dict(source='one', tier=2, disposition='admitted', accession='GSE1'),
                              dict(source='duplicate', tier=2, disposition='excluded', accession='GSE1')],
                  'drugs': [dict(source='one', source_id='yes', disposition='admitted',
                                 parent_inchikey='KEY', small_molecule=True, evidence_label='positive'),
                            dict(source='one', source_id='no', disposition='excluded',
                                 parent_inchikey='NEG', small_molecule=True, evidence_label='negative'),
                            dict(source='one', source_id='unknown', disposition='pending',
                                 parent_inchikey='UNK', small_molecule=True, evidence_label='pending')]}
        plan = tox.source_plan(ledger)
        self.assertEqual([p['source'] for p in plan], ['one'])
        self.assertEqual(plan[0]['drugs'].source_id.tolist(), ['yes'])
        ledger['sources'][1]['disposition'] = 'admitted'
        with self.assertRaisesRegex(tox.ReleaseError, 'independent'):
            tox.source_plan(ledger)

    def test_dataset_pooling_publishes_one_source_vote_and_replays(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            axis = axis_fixture()
            records = [source('studyA', ['G1', 'G2']),
                       source('studyB', ['G2', 'G3'])]
            tox.assemble_candidate(root, records, axis, {}, 'f' * 64)
            reader = disease_deg.DiseaseDeg(root)
            self.assertEqual(reader.categories(), ['tox'])
            self.assertEqual(reader.datasets('TOX'), ['studyA', 'studyB'])
            catalog = reader.catalog('TOX')
            self.assertEqual(catalog['n_deg'], 3)
            self.assertEqual(catalog['deg1000'], 3)
            self.assertFalse(catalog['eligible'])
            self.assertEqual(len(reader.table('TOX', 'studyA')), 2)
            self.assertEqual(reader.rank('TOX', 'studyA')['rank'].tolist()[:8],
                             [3, 1, 2, 3, 3, 3, 3, 3])
            self.assertEqual(len(tox.ToxCandidate(root).audit('studyA')), 4)
            self.assertEqual(tox.ToxCandidate(root).deg1000().gene_index.tolist(), [2, 1, 3])
            tox.assemble_candidate(root / 'replay', records, axis, {}, 'f' * 64)
            self.assertEqual(tox.ToxCandidate(root).digests(),
                             tox.ToxCandidate(root / 'replay').digests())


if __name__ == '__main__':
    unittest.main()
