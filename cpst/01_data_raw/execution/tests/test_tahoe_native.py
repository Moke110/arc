import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from anndata.io import read_elem

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_corpus import sha256
from native_tahoe import build_tahoe, verify_tahoe


class NativeTahoeTest(unittest.TestCase):
    def test_bulk_replays_off_axis_counts_even_when_symbols_collide(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source"
            (source / "metadata").mkdir(parents=True)
            (source / "data").mkdir()
            genes = pd.DataFrame({"gene_symbol": ["A", "A"], "ensembl_id": ["ens3", "ens4"],
                                  "token_id": [3, 4]})
            genes.to_parquet(source / "metadata/gene_metadata.parquet", index=False)
            pd.DataFrame({
                "sample": ["drug", "ctrl"], "plate": ["plate1", "plate1"],
                "drugname_drugconc": ["[('drug', 1, 'uM')]", "[('DMSO_TF', 0, 'uM')]"],
            }).to_parquet(source / "metadata/sample_metadata.parquet", index=False)
            pd.DataFrame({
                "Cell_ID_Cellosaur": ["L"], "cell_name": ["line"],
            }).to_parquet(source / "metadata/cell_line_metadata.parquet", index=False)
            rows = [
                {"plate": "plate1", "sample": "drug", "cell_line_id": "L",
                 "genes": [0, 3, 4], "expressions": [-2, 5, 9]},
                {"plate": "plate1", "sample": "ctrl", "cell_line_id": "L",
                 "genes": [0, 3], "expressions": [-2, 8]},
                {"plate": "plate1", "sample": "drug", "cell_line_id": "L",
                 "genes": [0, 4], "expressions": [-2, 7]},
            ]
            pq.write_table(pa.Table.from_pylist(rows), source / "data/train-00000.parquet")
            source_manifest = root / "manifest.sha256"
            source_manifest.write_text("".join(
                f"{sha256(path)}  {path.relative_to(source)}\n"
                for path in sorted(source.rglob("*.parquet"))
            ))
            mapping = pd.DataFrame({
                "gene_index": [0, 1], "symbol": ["A", "A"],
                "native_ensembl_id": ["ens3", "ens4"], "native_token_id": [3, 4],
            })
            output = root / "out"
            options = {"expected_genes": 2, "expected_shards": 1}
            built = build_tahoe(source, mapping, source_manifest, output, **options)
            replay = verify_tahoe(source, mapping, source_manifest, output, **options)
            self.assertEqual(replay["h5ad_sha256"], built["h5ad_sha256"])
            with h5py.File(output / "tahoe.h5ad", "r") as handle:
                self.assertEqual(read_elem(handle["var"]).index.tolist(), ["ens3", "ens4"])
                self.assertEqual(read_elem(handle["X"]).toarray().tolist(), [[8, 0], [5, 16]])
