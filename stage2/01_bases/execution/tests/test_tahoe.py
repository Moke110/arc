"""Public raw-release build/replay boundary; counts are hand-worked source truth."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

import h5py
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from anndata.io import read_elem
from scipy.sparse import csr_matrix

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bases.tahoe import build_tahoe, verify_tahoe


class TahoeReleaseTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / 'source'
        (self.source / 'metadata').mkdir(parents=True)
        (self.source / 'data').mkdir()
        pd.DataFrame([
            {'sample': 'A', 'plate': 'plate1', 'drugname_drugconc': "[('drug', 1, 'uM')]"},
            {'sample': 'C', 'plate': 'plate1', 'drugname_drugconc': "[('DMSO_TF', 0, 'uM')]"},
            {'sample': 'B', 'plate': 'plate2', 'drugname_drugconc': "[('drug', 1, 'uM')]"},
        ]).to_parquet(self.source / 'metadata/sample_metadata.parquet', index=False)
        pd.DataFrame([{'Cell_ID_Cellosaur': 'L1', 'cell_name': 'line1'}]).to_parquet(
            self.source / 'metadata/cell_line_metadata.parquet', index=False)
        rows = [
            [{'plate': 'plate1', 'sample': 'A', 'cell_line_id': 'L1', 'genes': [0, 3, 4, 5], 'expressions': [-2, 2, 7, 4]},
             {'plate': 'plate2', 'sample': 'B', 'cell_line_id': 'L1', 'genes': [0, 3, 5], 'expressions': [-2, 8, 1]}],
            [{'plate': 'plate1', 'sample': 'A', 'cell_line_id': 'L1', 'genes': [0, 3, 5], 'expressions': [-2, 5, 2]},
             {'plate': 'plate1', 'sample': 'C', 'cell_line_id': 'L1', 'genes': [0, 3, 5], 'expressions': [-2, 1, 3]}],
        ]
        for index, shard in enumerate(rows):
            pq.write_table(pa.Table.from_pylist(shard), self.source / f'data/train-{index:05d}-of-00002.parquet')
        self.mapping = pd.DataFrame([
            {'gene_index': 0, 'symbol': 'G1', 'hgnc_id': 'HGNC:1', 'ensembl_id': 'global1', 'native_token_id': '3', 'native_ensembl_id': 'source1'},
            {'gene_index': 1, 'symbol': 'G2', 'hgnc_id': 'HGNC:2', 'ensembl_id': 'global2', 'native_token_id': '5', 'native_ensembl_id': 'source2'},
        ])
        self.seal = self.root / 'manifest.sha256'
        self.seal.write_text(''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(self.source)}\n'
                                     for p in sorted(self.source.rglob('*.parquet'))))
        self.output = self.root / 'output'

    def test_mixed_shards_build_raw_counts_and_independently_replay(self):
        result = build_tahoe(self.source, self.mapping, self.seal, self.output,
                            expected_shards=2, expected_genes=2, batch_size=1)
        self.assertEqual(result['profiles'], 3)
        with h5py.File(self.output / 'tahoe.h5ad') as f:
            obs = read_elem(f['obs'])
            matrix = read_elem(f['X']).toarray()
            self.assertEqual(obs.index.tolist(), ['plate1:A:L1', 'plate1:C:L1', 'plate2:B:L1'])
            self.assertEqual(matrix.tolist(), [[7, 6], [1, 3], [8, 1]])
            self.assertEqual(obs.n_cells.tolist(), [2, 1, 1])
            self.assertEqual(obs.is_control.tolist(), [False, True, False])
            self.assertEqual(read_elem(f['var']).token_id.tolist(), ['3', '5'])
        replay = verify_tahoe(self.source, self.mapping, self.seal, self.output,
                              expected_shards=2, expected_genes=2, batch_size=2)
        self.assertEqual(replay['verified_profiles'], 3)
        self.assertEqual(replay['verified_cells'], 4)


    def test_completed_plates_resume_and_changed_source_is_rejected(self):
        first = build_tahoe(self.source, self.mapping, self.seal, self.output,
                            expected_shards=2, expected_genes=2)
        (self.output / 'tahoe.h5ad').unlink()
        resumed = build_tahoe(self.source, self.mapping, self.seal, self.output,
                              expected_shards=2, expected_genes=2)
        self.assertEqual(resumed['resumed_plates'], 2)
        self.assertEqual(resumed['h5ad_sha256'], first['h5ad_sha256'])
        source = next((self.source / 'data').glob('*.parquet'))
        with source.open('ab') as handle:
            handle.write(b'changed')
        with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
            build_tahoe(self.source, self.mapping, self.seal, self.output,
                        expected_shards=2, expected_genes=2)

    def test_competing_native_token_identities_cannot_publish(self):
        self.mapping.loc[1, 'native_token_id'] = '3'
        with self.assertRaisesRegex(ValueError, 'ambiguous'):
            build_tahoe(self.source, self.mapping, self.seal, self.output,
                        expected_shards=2, expected_genes=2)
        self.assertFalse((self.output / 'manifest.json').exists())


if __name__ == '__main__':
    unittest.main()
