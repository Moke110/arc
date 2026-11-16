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

    def test_peptide_identity_remains_excluded_under_an_alias(self):
        annotation = pd.DataFrame([dict(lookup_name=n, parent_inchikey='PEPTIDE', atc='', moa='')
                                   for n in ['Cyclosporine', 'Cyclosporin A']])
        fact = dict(parent_inchikey='PEPTIDE', label='positive', url='https://example.org',
                    locator='label', retrieved_at='2026-11-16', rationale='risk', uncertainty='risk')
        result = tox_admission.assess_drugs([dict(source='panel', source_id='alias', drug_name='Cyclosporin A')],
                                          annotation, [fact])
        self.assertFalse(result.iloc[0].small_molecule)
        self.assertEqual(result.iloc[0].disposition, 'excluded')

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

    def test_reassessed_panels_enter_by_independent_identity_evidence(self):
        import gzip
        import io
        import tarfile
        from bases import cardiotox
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            def put(name, text):
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                if name.endswith('.gz'):
                    with gzip.open(path, 'wt') as handle:
                        handle.write(text)
                else:
                    path.write_text(text)
            put(cardiotox.DTOXS_SERIES, '!Sample_characteristics_ch1\t"treatment condition: X"\n!Sample_characteristics_ch1\t"drug name: known"\n')
            path = root / cardiotox.DTOXS_PARAMS
            path.parent.mkdir(parents=True, exist_ok=True)
            payload = b'Group\tDrug\tType\na\tX\ta\n'
            with tarfile.open(path, 'w:gz') as archive:
                info = tarfile.TarInfo('Conv-RNAseq-Drug-Groups.tsv')
                info.size = len(payload)
                archive.addfile(info, io.BytesIO(payload))
            put(cardiotox.GSE244740_SERIES, '!Sample_description\t"p_w"\n!Sample_characteristics_ch1\t"treatment: DMSO, 0"\n!Sample_characteristics_ch1\t"time: 24h"\n')
            put('panel.txt.gz', '!Sample_characteristics_ch1\t"treatment: known_0.2_uM"\t"treatment: unknown_1_uM"\t"treatment: DMSO"\n')
            annotation = pd.DataFrame([dict(lookup_name=n, parent_inchikey=k, atc='', moa='')
                                       for n,k in [('known','KEY'),('unknown','OTHER')]])
            put(tox_admission.ANNOTATION, annotation.to_csv(index=False))
            fact = dict(parent_inchikey='KEY', label='positive', url='https://example.org/primary',
                        locator='reference', retrieved_at='2026-11-16', rationale='positive label', uncertainty='risk')
            put(tox_admission.EVIDENCE, json.dumps(dict(drugs=[fact], pinned_sources=[], reviewed_at='2026-11-16', review_scope='fixed independent reference labels')))
            source = dict(source='reassessed', disposition='admitted', tier=2,
                          drug_panel=dict(path='panel.txt.gz', characteristic='treatment', dose_suffix=True))
            put(tox_admission.SOURCES, json.dumps(dict(sources=[source], discovery={})))
            output = root / 'candidate'
            tox_admission.build_ledger(root, output)
            reader = tox_admission.ToxAdmission(output, root)
            drugs = reader.drugs('reassessed').set_index('drug_name')
            self.assertEqual(set(drugs.index), {'known', 'unknown'})
            self.assertEqual(drugs.loc['known','disposition'], 'admitted')
            self.assertEqual(drugs.loc['unknown','disposition'], 'pending')
            self.assertIn('panel.txt.gz', reader.value['inputs'])
            self.assertEqual(drugs.loc['unknown','reviewed_at'], '2026-11-16')
            self.assertEqual(drugs.loc['unknown','review_scope'], 'fixed independent reference labels')
            self.assertEqual(drugs.loc['unknown','reviewed_evidence'], tox_admission.EVIDENCE)

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
