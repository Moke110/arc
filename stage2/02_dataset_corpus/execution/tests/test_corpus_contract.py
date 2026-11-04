"""G1 structure + count contract for the four native-axis CM corpora.

Run from the repository root (the four h5ads are local-only artifacts, so the
suite skips when they are absent):

    python -m unittest discover -s stage2/02_dataset_corpus/execution/tests -v

Every assertion here reads the written files the way a consumer would: h5py
header reads and ``read_elem`` on ``obs``/``var`` only, never a whole matrix.
The layer formula is checked on a few sampled rows; non-negativity is checked in
bounded row chunks.
"""

from __future__ import annotations

import hashlib
import sys
import unittest
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from cmdata import paths  # noqa: E402

CORPORA = list(paths.CORPUS_IDS)
FILES = {corpus_id: paths.corpus_h5ad(corpus_id) for corpus_id in CORPORA}
RECORD = paths.CORPUS_BUILD_RECORD
EXCLUDED = paths.EXCLUDED_COMPOUNDS_RECORD


def _skip_reason() -> str:
    missing = [str(path) for path in FILES.values() if not path.exists()]
    if not RECORD.exists():
        missing.append(str(RECORD))
    if not EXCLUDED.exists():
        missing.append(str(EXCLUDED))
    return "" if not missing else f"corpus artifacts not built: {missing}"


@unittest.skipUnless(not _skip_reason(), _skip_reason())
class CorpusContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = __import__("json").loads(RECORD.read_text())

    def test_header_structure(self):
        for corpus_id, path in FILES.items():
            with self.subTest(corpus=corpus_id), h5py.File(path, "r") as handle:
                expected = paths.EXPECTED_CORPUS_COUNTS[corpus_id]
                self.assertEqual(handle["X"].shape, (expected["rows"],
                                                     expected["axis"]))
                self.assertEqual(handle["X"].dtype, np.int32)
                self.assertEqual(handle["layers/log1p_cp10k"].dtype, np.float32)
                obs = read_elem(handle["obs"])
                var = read_elem(handle["var"])
                self.assertEqual(list(obs.columns), paths.CORPUS_OBS)
                self.assertTrue(obs.index.is_unique)
                self.assertTrue((obs.index.astype(str) != "").all())
                self.assertEqual(list(var.columns),
                                 ["hgnc_id", "entrez_id", "ensembl_id",
                                  "match_type"])
                self.assertTrue(var.index.is_unique)
                self.assertEqual(int(obs["lib_size"].min()) > 0, True)
                # Empty containers are written as empty HDF5 groups.
                for name in ("obsm", "varm", "obsp", "varp", "uns"):
                    if name in handle:
                        self.assertEqual(len(handle[name].keys()), 0,
                                         f"{corpus_id}: {name} not empty")
                for forbidden in ("pert_idx", "cell_idx", "dose_idx"):
                    self.assertNotIn(forbidden, obs.columns)

    def test_counts(self):
        for corpus_id, path in FILES.items():
            with self.subTest(corpus=corpus_id), h5py.File(path, "r") as handle:
                obs = read_elem(handle["obs"])
                expected = paths.EXPECTED_CORPUS_COUNTS[corpus_id]
                control = obs["is_control"].to_numpy()
                blank = obs["is_blank"].to_numpy()
                perturbed = ~control & ~blank
                self.assertEqual(len(obs), expected["rows"])
                self.assertEqual(handle["X"].shape[1], expected["axis"])
                self.assertEqual(int(perturbed.sum()), expected["perturbed"])
                self.assertEqual(int(control.sum()), expected["control"])
                self.assertEqual(int(blank.sum()), expected["blank"])
                self.assertEqual(
                    int(obs.loc[perturbed, "drug_key"].astype(str).nunique()),
                    expected["chemicals"])
                self.assertEqual(set(obs["dataset_id"].astype(str)),
                                 {corpus_id})
                # chemicals never include controls or blanks
                self.assertFalse((obs.loc[~perturbed, "drug_key"]
                                  .astype(str) != "").any())

    def test_counts_non_negative_and_layer_formula(self):
        for corpus_id, path in FILES.items():
            with self.subTest(corpus=corpus_id), h5py.File(path, "r") as handle:
                n_rows = handle["X"].shape[0]
                # Bounded chunks: never materialise the whole matrix.
                for start in range(0, n_rows, 256):
                    stop = min(start + 256, n_rows)
                    block = handle["X"][start:stop]
                    self.assertGreaterEqual(int(block.min()), 0)
                positions = sorted({0, n_rows // 2, n_rows - 1})
                for row in positions:
                    counts = np.asarray(handle["X"][row])
                    layer = np.asarray(handle["layers/log1p_cp10k"][row])
                    expected = np.log1p(counts / counts.sum() * 1e4)
                    np.testing.assert_allclose(layer, expected,
                                               rtol=1e-5, atol=1e-6)

    def test_build_record_matches_files(self):
        for corpus_id, path in FILES.items():
            with self.subTest(corpus=corpus_id):
                entry = self.record["corpora"][corpus_id]
                digest = hashlib.sha256()
                with open(path, "rb") as handle:
                    for block in iter(lambda: handle.read(1 << 20), b""):
                        digest.update(block)
                self.assertEqual(entry["file_sha256"], digest.hexdigest())
                self.assertEqual(entry["bytes"], path.stat().st_size)
                self.assertEqual(entry["rows"],
                                 paths.EXPECTED_CORPUS_COUNTS[corpus_id]["rows"])
                self.assertEqual(entry["axis_length"],
                                 paths.EXPECTED_CORPUS_COUNTS[corpus_id]["axis"])
                self.assertEqual(len(entry["content_sha256"]), 64)
                self.assertTrue(all(gate["passed"] for gate in entry["gates"]))
                # Provenance: the record carries the source accession the
                # coverage README shows, and no CM corpus is left blank.
                self.assertEqual(entry["accession"],
                                 paths.source_accession(entry["source_dataset"]))
                self.assertNotEqual(entry["accession"], "")

    def test_excluded_compounds_accounting(self):
        frame = pd.read_csv(EXCLUDED)
        compounds = frame[frame["entity_type"] != "technical_blank"]
        blanks = frame[frame["entity_type"] == "technical_blank"]
        self.assertEqual(len(compounds), 8)          # 8 non-SMILES keys
        self.assertEqual(int(compounds["n_wells"].sum()), 143)
        self.assertEqual(len(blanks), 1)
        self.assertEqual(int(blanks["n_wells"].sum()), 32)
        expected = {"Bevacizumab", "Cetuximab", "Rituximab", "Trastuzumab",
                    "Endothelin-1", "IGF-1", "TNF-alpha", "ET-1"}
        self.assertEqual(set(compounds["compound_name"]), expected)
        self.assertTrue((frame["reason"].astype(str) != "").all())

    def test_tenv_blanks_are_retained_rows(self):
        with h5py.File(FILES["tenv"], "r") as handle:
            obs = read_elem(handle["obs"])
        blanks = obs[obs["is_blank"]]
        self.assertEqual(len(blanks), 32)
        self.assertFalse(blanks["is_control"].any())
        self.assertTrue((blanks["drug_key"].astype(str) == "").all())
        self.assertTrue((blanks["dose_um"] == 0).all())


if __name__ == "__main__":
    unittest.main()
