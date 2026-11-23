"""Verify the Stage 2 successor fingerprint contract at its public boundary."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'util/output/api'))
from stage2_util.ecfp4 import open_ecfp4, publish_ecfp4


class Ecfp4ContractTest(unittest.TestCase):
    def test_legacy_and_successor_keep_identical_bits_and_parent_order(self):
        for schema, chirality in (
            ('tenv5k-inputs-v2', 'includeChirality'),
            ('tsg-dose-round1-axis-local-inputs-v1', 'include_chirality'),
        ):
            with self.subTest(schema=schema), tempfile.TemporaryDirectory() as temporary:
                legacy = Path(temporary) / 'legacy'
                legacy.mkdir()
                bits = np.zeros((2, 2048), dtype=np.uint8)
                bits[0, [10, 145, 319]] = 1
                bits[1, [10, 1008, 2047]] = 1
                np.save(legacy / 'morgan.npy', bits)
                (legacy / 'parents.csv').write_text(
                    'drug_parent_inchikey,drug_parent_smiles\nFIRST,CCO\nSECOND,CCC\n'
                )
                (legacy / 'config.json').write_text(json.dumps({
                    'morgan': {'radius': 2, 'bits': 2048, chirality: True, 'rdkit': 'test'}
                }))
                files = {path.name: path.read_bytes() for path in legacy.iterdir()}
                (legacy / 'manifest.json').write_text(json.dumps({
                    'schema': schema,
                    'outputs': {name: hashlib.sha256(data).hexdigest()
                                for name, data in files.items()},
                }))
                files['manifest.json'] = (legacy / 'manifest.json').read_bytes()

                candidate = Path(temporary) / 'candidate'
                publish_ecfp4(legacy, candidate)
                np.testing.assert_array_equal(open_ecfp4(legacy), bits)
                np.testing.assert_array_equal(
                    open_ecfp4(candidate, source_manifest=legacy / 'manifest.json'), bits
                )
                unbound_manifest = Path(temporary) / 'unbound.json'
                unbound_manifest.write_text('{}')
                with self.assertRaisesRegex(ValueError, 'source manifest mismatch'):
                    open_ecfp4(candidate, source_manifest=unbound_manifest)
                self.assertEqual((candidate / 'ecfp4.npy').read_bytes(), files['morgan.npy'])
                self.assertEqual((candidate / 'parents.csv').read_bytes(), files['parents.csv'])
                self.assertEqual({path.name: path.read_bytes() for path in legacy.iterdir()}, files)
                contract = json.loads((candidate / 'manifest.json').read_text())
                self.assertEqual(contract['ecfp4'], {
                    'atom_invariants': 'connectivity', 'binary': True, 'bits': 2048,
                    'include_chirality': True, 'radius': 2, 'rdkit': 'test',
                })
                (candidate / 'parents.csv').write_text('corrupted')
                with self.assertRaisesRegex(ValueError, 'digest mismatch'):
                    open_ecfp4(candidate)
