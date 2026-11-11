import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases.identity import GeneIdentity, IdentityError

class IdentityContract(unittest.TestCase):
    def test_index_roundtrip_and_source_specific_native_identity(self):
        axis = [{'symbol': 'A', 'hgnc_id': 'HGNC:1', 'ensembl_id': 'ENSG1'},
                {'symbol': 'FAM153B', 'hgnc_id': 'HGNC:27323', 'ensembl_id': 'ENSG00000182230'}]
        crosswalk = [{'dataset_id':'tahoe', 'source_gene_id':'ENSG00000289731', 'hgnc_id':'HGNC:27323', 'match_type':'direct'},
                     {'dataset_id':'tahoe', 'source_gene_id':'61482', 'hgnc_id':'HGNC:27323', 'match_type':'direct'}]
        identity = GeneIdentity(axis, crosswalk)
        self.assertEqual(identity.indices(['FAM153B', 'A']), [1, 0])
        self.assertEqual(identity.labels([1, 0]), ['FAM153B', 'A'])
        self.assertEqual(identity.indices(['ENSG00000289731', '61482'], source='tahoe'), [1, 1])
        self.assertEqual(identity.native_ids(1, 'tahoe'), ['61482', 'ENSG00000289731'])
        with self.assertRaises(IdentityError):
            identity.indices(['ENSG00000289731'], source='tenv')


class AmbiguityContract(unittest.TestCase):
    def test_duplicate_evidence_is_safe_but_conflicting_aliases_fail(self):
        axis = [{'symbol':'A','hgnc_id':'HGNC:1','ensembl_id':'ENSG1'},
                {'symbol':'B','hgnc_id':'HGNC:2','ensembl_id':'ENSG2'}]
        annotation = [{'hgnc_id':'HGNC:1','alias_symbol':'DUP|OLD'},
                      {'hgnc_id':'HGNC:2','alias_symbol':'DUP'},
                      {'hgnc_id':'HGNC:1','alias_symbol':'OLD'}]
        identity = GeneIdentity(axis, annotation=annotation)
        self.assertEqual(identity.indices(['OLD']), [0])
        with self.assertRaisesRegex(IdentityError, 'ambiguous'):
            identity.indices(['DUP'])
        with self.assertRaises(IdentityError):
            identity.labels([-1])

    def test_conflicting_native_identity_fails(self):
        axis = [{'symbol':'A','hgnc_id':'HGNC:1'}, {'symbol':'B','hgnc_id':'HGNC:2'}]
        rows = [{'dataset_id':'x','source_gene_id':'same','hgnc_id': h}
                for h in ('HGNC:1','HGNC:2')]
        with self.assertRaisesRegex(IdentityError, 'ambiguous'):
            GeneIdentity(axis, rows).indices(['same'], source='x')

    def test_duplicate_axis_symbols_are_rejected(self):
        with self.assertRaises(IdentityError):
            GeneIdentity([{'symbol':'A','hgnc_id':'HGNC:1'},
                          {'symbol':'A','hgnc_id':'HGNC:2'}])

    def test_source_collision_with_off_axis_gene_fails(self):
        identity = GeneIdentity([{'symbol':'A','hgnc_id':'HGNC:1'}], [
            {'dataset_id':'x','source_gene_id':'id','hgnc_id':'HGNC:1'},
            {'dataset_id':'x','source_gene_id':'id','hgnc_id':'HGNC:999'}])
        with self.assertRaisesRegex(IdentityError, 'ambiguous'):
            identity.indices(['id'], source='x')

if __name__ == '__main__': unittest.main()
