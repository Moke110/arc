"""Independent drug admission through the public ledger seam (#215)."""
import sys
import unittest
import tempfile
import json
import hashlib
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import public_api
from bases import tox_admission


class AdmissionTest(unittest.TestCase):
    def test_positive_identity_passes_without_mechanism_or_peer_drugs(self):
        records = [{'source': 'single', 'source_id': 'x', 'drug_name': 'A'}]
        annotation = pd.DataFrame([dict(lookup_name='A', dataset_inchikey='',
            parent_inchikey='KEY', smiles='CC', atc='', moa='', mesh_classes='', notes='')])
        evidence = [{'drug_name': 'A', 'parent_inchikey': 'KEY', 'label': 'positive',
                     'url': 'https://example.org/primary', 'locator': 'Table 1',
                     'retrieved_at': '2026-11-16', 'rationale': 'author positive label',
                     'uncertainty': 'drug property; no exposure-level claim'}]
        result = tox_admission.assess_drugs(records, annotation, evidence)
        self.assertEqual(result.iloc[0].disposition, 'admitted')
        self.assertEqual(result.iloc[0].evidence_label, 'positive')

    def test_duplicate_source_code_is_rejected_instead_of_double_voting(self):
        records = [dict(source='one', source_id='x', drug_name='A')] * 2
        annotation = pd.DataFrame([dict(lookup_name='A', dataset_inchikey='',
            parent_inchikey='KEY', smiles='CC', atc='', moa='')])
        with self.assertRaisesRegex(tox_admission.ReleaseError, 'duplicate'):
            tox_admission.assess_drugs(records, annotation, [])

    def test_conflicting_structural_names_do_not_transfer_positive_evidence(self):
        rows = [dict(lookup_name='A', dataset_inchikey=d, parent_inchikey=k,
                     smiles='CC', atc='', moa='') for d,k in [('', 'CANON'), ('OTHER','OTHER')]]
        fact = dict(parent_inchikey='OTHER', label='positive', url='https://example.org',
                    locator='label', retrieved_at='2026-11-16', rationale='positive', uncertainty='risk')
        result = tox_admission.assess_drugs([dict(source='one',source_id='x',drug_name='A')],
                                          pd.DataFrame(rows), [fact])
        self.assertEqual(result.iloc[0].disposition, 'excluded')
        self.assertIn('identity', result.iloc[0].reason)

    def test_reader_rejects_modified_admission_payload(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = dict(schema='stage2-tox-admission-v1', inputs={}, drugs=[])
            payload = json.dumps(data).encode()
            (root/'admission.json').write_bytes(payload)
            (root/'manifest.json').write_text(json.dumps(dict(schema='stage2-tox-admission-manifest-v1',
                admission_sha256=hashlib.sha256(payload).hexdigest())))
            data['drugs'] = [dict(disposition='admitted')]
            (root/'admission.json').write_text(json.dumps(data))
            with self.assertRaisesRegex(tox_admission.ReleaseError, 'payload'):
                tox_admission.ToxAdmission(root)

    def test_negative_pending_conflict_and_non_small_molecules_are_never_admitted(self):
        names = ['negative', 'unknown', 'conflict', 'Trastuzumab', 'Carfilzomib', 'Cyclosporine', 'unresolved']
        annotation = pd.DataFrame([dict(lookup_name=n, dataset_inchikey='', parent_inchikey=n,
                        smiles='CC', atc='ambiguous', moa='multiple') for n in names[:-1]])
        records = [dict(source='one', source_id=n, drug_name=n) for n in names]
        facts = []
        for n, labels in [('negative',['negative']),('conflict',['positive','negative']),
                          ('Trastuzumab',['positive']),('Carfilzomib',['positive']),('Cyclosporine',['positive'])]:
            facts.extend(dict(parent_inchikey=n,label=l,url='https://example.org/primary',
                              locator='reference',retrieved_at='2026-11-16',rationale='label',
                              uncertainty='drug property') for l in labels)
        result = tox_admission.assess_drugs(records, annotation, facts).set_index('source_id')
        self.assertEqual(result.loc['negative','evidence_label'], 'negative')
        self.assertEqual(result.loc['unknown','disposition'], 'pending')
        self.assertEqual(result.loc['conflict','disposition'], 'pending')
        self.assertEqual(result.loc['unresolved','disposition'], 'excluded')
        self.assertEqual(result.loc['Trastuzumab','disposition'], 'excluded')
        self.assertEqual(result.loc['Carfilzomib','disposition'], 'excluded')
        self.assertEqual(result.loc['Cyclosporine','disposition'], 'excluded')
        self.assertFalse(result.loc['Cyclosporine','small_molecule'])
        self.assertFalse(result.disposition.eq('admitted').any())

if __name__ == '__main__':
    unittest.main()
