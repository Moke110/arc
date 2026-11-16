"""Independent drug admission through the public ledger seam (#215)."""
import sys
import unittest
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

if __name__ == '__main__':
    unittest.main()
