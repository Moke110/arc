"""h5ad writer tests: round-trip, dtypes, the frozen obs/var contract.

CM and Tahoe writers are exercised on small synthetic matrices (the 18,164-gene
var is real, the cell counts are not).  The pinned reader check runs in a
subprocess only when the `arc-stage2-bench-prnet` interpreter exists.
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

EXECUTION = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXECUTION))

import anndata as ad
import h5py
import numpy as np
import pandas as pd

import contract
import writers

PINNED_READER = Path("/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python")
HAVE_CM_DATA = contract.CM_VIEW.exists()

SAMPLE_HASH_ROWS = [
    (1, "B2", "VEH", 0.0, 1, "DMSO", "VEH"),
    (1, "B3", "VEH", 0.0, 2, "DMSO", "VEH"),
    (1, "B4", "Media", 0.0, 1, "Media", "Media"),
    (1, "B7", "CHEM1", 0.2, 1, "DMSO", "Acetaminophen"),
    (1, "B8", "CHEM1", 1.0, 2, "DMSO", "Acetaminophen"),
    (2, "B2", "VEH", 0.0, 1, "DMSO", "VEH"),
    (2, "B3", "CHEM2", 10.0, 1, "DMSO", "Citric acid"),
    (2, "B4", "CHEM2", 0.2, 2, "DMSO", "Citric acid"),
]


def cm_well(plate_no: int, raw_well: str) -> str:
    letter, digits = raw_well[0], raw_well[1:]
    return f"temp_env::Plate{plate_no}-{letter}{int(digits):02d}"


def cm_inherited_obs() -> pd.DataFrame:
    rows = []
    for plate_no, raw_well, chem, dose, replicate, solvent, name in SAMPLE_HASH_ROWS:
        is_control = chem in ("VEH", "Media")
        rows.append({
            "well": cm_well(plate_no, raw_well),
            "dataset_id": "temp_env",
            "compound_name": "DMSO" if chem == "VEH" else name,
            "drug_key": "" if is_control else f"KEY-{chem}",
            "inchikey": "" if is_control else f"AAAA-BBBB-{chem}",
            "smiles": "" if is_control else "CC",
            "dose_um": float(dose),
            "dose_source": "" if is_control else "measured",
            "time_h": 24.0,
            "is_control": bool(is_control),
            "vehicle": solvent,
            "plate": f"Plate{plate_no}",
            "control_group": f"temp_env::Plate{plate_no}",
            "cell_type": "ipsc_cm",
            "donor": "NA",
            "replicate": int(replicate),
            "atc": "",
            "moa": "",
            "mesh_classes": "",
            "pharm_atc_moa": False,
            "pharm_mesh": False,
            "annotation_source": "",
            "lib_size": 1000,
            "n_genes_detected": 10,
        })
    obs = pd.DataFrame(rows)
    obs.index = pd.Index(obs["well"].astype(str), name="well")
    return obs


def cm_hash_table() -> pd.DataFrame:
    return pd.DataFrame(SAMPLE_HASH_ROWS, columns=[
        "Plate_ID", "Well_ID", "Chemical_Index", "Chemical_Concentration_uM",
        "Replicate_Number", "Solvent", "Chemical_name"])


def cm_assembled_obs():
    obs = cm_inherited_obs()
    treated_keys = sorted({k for k in obs["drug_key"] if k})
    eligibility = pd.DataFrame({"drug_key": treated_keys,
                                "eligibility_decision": ["include"] * len(treated_keys),
                                "elig_class": ["stub"] * len(treated_keys)})
    qc = pd.DataFrame({"well": obs.index.astype(str),
                       "qc_pass": True, "qc_reason": "ok"})
    return writers.assemble_cm_obs(obs, cm_hash_table(), eligibility, qc), eligibility, qc


TAHOE_LINES = ("CVCL_0001", "CVCL_0002")


def tahoe_inherited_obs() -> pd.DataFrame:
    rows = []
    for plate in ("plate1", "plate13"):
        for line in TAHOE_LINES:
            for well in ("smp_1", "smp_2", "smp_3"):
                is_control = well == "smp_1"
                rows.append({
                    "well": f"tahoe::{plate}:{well}:{line}",
                    "dataset_id": "tahoe",
                    "compound_name": "DMSO_TF" if is_control else "DrugA",
                    "drug_key": "" if is_control else "KEYA",
                    "inchikey": "" if is_control else "AAAA-BBBB-CCCC",
                    "smiles": "" if is_control else "CC",
                    "dose_um": 0.0 if is_control else 0.5,
                    "dose_source": "nominal_atlas",
                    "time_h": np.nan,
                    "is_control": bool(is_control),
                    "vehicle": "DMSO_TF" if is_control else "",
                    "plate": plate,
                    "control_group": f"tahoe::{plate}::{line}",
                    "cell_type": np.nan,
                    "donor": np.nan,
                    "replicate": 0,
                    "atc": "", "moa": "", "mesh_classes": "",
                    "pharm_atc_moa": False, "pharm_mesh": False,
                    "annotation_source": "",
                    "lib_size": 1000, "n_genes_detected": 10,
                    "cell_line_id": line,
                    "cell_name": f"line-{line}",
                    "n_cells": 500,
                    "moa_broad": "", "moa_fine": "", "pubchem_cid": "",
                    "has_structure": True,
                    "has_matched_control": True,
                    "SMILES": "" if is_control else "CC",
                    "dose": 0.0 if is_control else 0.5,
                    "paired_control_index": "",
                })
    obs = pd.DataFrame(rows)
    obs.index = pd.Index(obs["well"].astype(str), name="well")
    return obs


def tahoe_assembled_obs():
    raw = tahoe_inherited_obs()
    qc = pd.DataFrame({"well": raw.index.astype(str), "qc_pass": True,
                       "qc_reason": "ok"})
    return writers.assemble_tahoe_obs(raw, qc)


def small_var(n_rows: int = 6):
    if HAVE_CM_DATA:
        axis = contract.load_gene_axis()
        measured = contract.gene_measured_cm(axis)
        nnz = np.arange(len(axis), dtype=np.int64) % 7
        return writers.var_frame(axis, measured, nnz)
    names = [f"G{i}" for i in range(5)]
    axis = pd.DataFrame({"hgnc_id": [f"HGNC:{i}" for i in range(5)],
                         "entrez_id": [str(i) for i in range(5)],
                         "ensembl_id": [f"ENSG{i:011d}" for i in range(5)],
                         "match_type": ["direct"] * 5}, index=pd.Index(names, name="symbol"))
    return writers.var_frame(axis, pd.Series(True, index=axis.index), np.ones(5, dtype=int))


class CmWriterTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.obs, self.eligibility, self.qc = cm_assembled_obs()
        self.var = small_var()
        self.X = np.arange(len(self.obs) * len(self.var), dtype=np.int32).reshape(
            len(self.obs), len(self.var))
        self.log1p = np.log1p(self.X.astype(np.float32))

    def tearDown(self):
        self.tmp.cleanup()

    def test_roundtrip(self):
        path = self.dir / "cm.h5ad"
        writers.write_cm_h5ad(path, self.obs, self.X, self.log1p, self.var,
                              {"schema": "06.3_input/1"})
        back = ad.read_h5ad(path)
        self.assertEqual(back.shape, (len(self.obs), len(self.var)))
        self.assertEqual(list(back.obs.columns), contract.CM_OBS)
        self.assertEqual(back.obs.index.name, "well")
        self.assertEqual(back.X.dtype, np.int32)
        self.assertEqual(back.layers["log1p_cp10k"].dtype, np.float32)
        np.testing.assert_array_equal(np.asarray(back.X), self.X)
        self.assertEqual(list(back.var.columns), contract.VAR_COLUMNS)
        self.assertEqual(list(back.var.index.astype(str)),
                         list(self.var.index.astype(str)))
        self.assertNotIn("control_mean", back.uns)
        self.assertEqual(back.uns["schema"], "06.3_input/1")

    def test_rejects_wrong_dtype(self):
        with self.assertRaises(contract.ContractError):
            writers.write_cm_h5ad(self.dir / "x.h5ad", self.obs,
                                  self.X.astype(np.int64), self.log1p, self.var, {})
        with self.assertRaises(contract.ContractError):
            writers.write_cm_h5ad(self.dir / "x.h5ad", self.obs, self.X,
                                  self.log1p.astype(np.float64), self.var, {})

    def test_rejects_missing_columns(self):
        broken = self.obs.drop(columns=["main_set"])
        with self.assertRaises(contract.ContractError):
            writers.write_cm_h5ad(self.dir / "x.h5ad", broken, self.X, self.log1p,
                                  self.var, {})

    def test_rejects_negative_counts(self):
        bad = self.X.copy()
        bad[0, 0] = -1
        with self.assertRaises(contract.ContractError):
            writers.write_cm_h5ad(self.dir / "x.h5ad", self.obs, bad, self.log1p,
                                  self.var, {})


class CmAssemblyTest(unittest.TestCase):
    @unittest.skipUnless(HAVE_CM_DATA, f"frozen data not under {contract.DATA_ROOT}")
    def test_raw_design_rejoin(self):
        obs = cm_inherited_obs()
        out = writers.join_raw_design(obs, cm_hash_table())
        self.assertEqual(out["physical_well"].tolist(), ["B2", "B3", "B4", "B7",
                                                        "B8", "B2", "B3", "B4"])
        self.assertEqual(out["replicate_raw"].tolist(), [1, 2, 1, 1, 2, 1, 1, 2])
        self.assertEqual(out["plate_no"].tolist(), [1, 1, 1, 1, 1, 2, 2, 2])
        self.assertEqual(set(out["solvent"].astype(str)), {"DMSO", "Media"})

    @unittest.skipUnless(HAVE_CM_DATA, f"frozen data not under {contract.DATA_ROOT}")
    def test_keep_mask_only_qc_pass(self):
        obs, _, _ = cm_assembled_obs()
        mask = writers.cm_keep_mask(obs)
        self.assertEqual(int(mask.sum()), len(obs))
        obs.loc[obs.index[0], "qc_pass"] = False
        mask = writers.cm_keep_mask(obs)
        self.assertEqual(int(mask.sum()), len(obs) - 1)


class TahoeAssemblyTest(unittest.TestCase):
    def test_identity_columns(self):
        obs = tahoe_assembled_obs()
        self.assertEqual(list(obs.columns), contract.TAHOE_OBS)
        self.assertEqual(obs.index.name, "well")
        self.assertTrue(obs["physical_well"].str.startswith("smp_").all())
        self.assertEqual(set(obs["tahoe_split"].astype(str)), {"train", "valid"})
        self.assertTrue(obs.loc[obs["plate"].astype(str) == "plate13",
                                "tahoe_split"].astype(str).eq("valid").all())
        self.assertTrue(obs["time_h"].isna().all())
        self.assertEqual(set(obs["time_source"].astype(str)),
                         {contract.TIME_SOURCE_TAHOE})
        treated = obs.loc[~obs["is_control"].astype(bool)]
        self.assertTrue(treated["tahoe_valid_group"].str.contains(r"KEYA\|0\.5\|CVCL").all())
        self.assertTrue(obs.loc[obs["is_control"].astype(bool), "tahoe_usable"].eq(False).all())

    def test_raw_time_column_is_not_overwritten(self):
        raw = tahoe_inherited_obs()
        raw.loc[raw.index[0], "time_h"] = 24.0
        qc = pd.DataFrame({"well": raw.index.astype(str), "qc_pass": True,
                           "qc_reason": "ok"})
        with self.assertRaises(contract.ContractError):
            writers.assemble_tahoe_obs(raw, qc)

    def test_missing_qc_coverage_fails_closed(self):
        raw = tahoe_inherited_obs()
        qc = pd.DataFrame({"well": raw.index.astype(str)[:2], "qc_pass": True,
                           "qc_reason": "ok"})
        with self.assertRaises(contract.ContractError):
            writers.assemble_tahoe_obs(raw, qc)

    def test_control_groups_and_pairs(self):
        obs = tahoe_assembled_obs()
        groups = writers.tahoe_control_groups(obs)
        self.assertEqual(len(groups), 4)      # 2 plates x 2 lines
        self.assertTrue((groups["n_controls"] == 1).all())
        self.assertTrue((groups["n_qc_pass"] == 1).all())
        self.assertEqual([len(rows) for rows in groups["control_rows"]], [1, 1, 1, 1])
        pairs = writers.tahoe_pairs(obs)
        treated = int((~obs["is_control"].astype(bool)).sum())
        self.assertEqual(len(pairs["treated_row"]), treated)
        self.assertTrue(set(pairs["control_group"].astype(str))
                        <= set(groups["control_group"].astype(str)))


class TahoeWriterTest(unittest.TestCase):
    def setUp(self):
        if not HAVE_CM_DATA:
            self.skipTest(f"frozen data not under {contract.DATA_ROOT}")
        AnnData = ad.AnnData
        import scipy.sparse as sp
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.obs = tahoe_assembled_obs()
        self.var = small_var()
        n = len(self.obs)
        rng = np.random.default_rng(0)
        dense = rng.integers(0, 3, size=(n, len(self.var))).astype(np.int32)
        X = sp.csr_matrix(dense)
        layer = sp.csr_matrix(np.log1p(dense.astype(np.float32)))
        raw = tahoe_inherited_obs()
        self.src = self.dir / "source.h5ad"
        src = AnnData(X=X, obs=raw, var=self.var,
                      layers={"log1p_cp10k": layer},
                      uns={"control_mean": pd.DataFrame(
                          np.zeros((1, len(self.var))),
                          index=["tahoe::plate1::CVCL_0001"],
                          columns=self.var.index.astype(str))})
        src.write_h5ad(self.src)

    def tearDown(self):
        self.tmp.cleanup()

    def test_full_copy_path(self):
        dst = self.dir / "tahoe_transfer.h5ad"
        writers.write_tahoe_h5ad(self.src, dst, self.obs, self.var,
                                 {"schema": "06.3_input/1"})
        facts = writers.verify_tahoe_h5ad(dst, expected_rows=len(self.obs))
        self.assertEqual(facts["rows"], len(self.obs))
        self.assertEqual(facts["cols"], len(self.var))
        self.assertGreater(facts["nnz"], 0)
        back = ad.read_h5ad(dst)
        self.assertEqual(list(back.obs.columns), contract.TAHOE_OBS)
        self.assertEqual(back.X.dtype, np.int32)
        self.assertEqual(back.layers["log1p_cp10k"].dtype, np.float32)
        self.assertNotIn("control_mean", back.uns)
        # the source matrix survived the streaming copy byte-for-byte
        src = ad.read_h5ad(self.src)
        np.testing.assert_array_equal(np.asarray(back.X.todense()),
                                      np.asarray(src.X.todense()))
        with h5py.File(dst) as fh:
            self.assertEqual(fh["X"].attrs["encoding-type"], "csr_matrix")

    def test_slice_path(self):
        dst = self.dir / "tahoe_slice.h5ad"
        writers.write_tahoe_h5ad(self.src, dst, self.obs, self.var,
                                 {"schema": "06.3_input/1"}, row_limit=3)
        facts = writers.verify_tahoe_h5ad(dst, expected_rows=3)
        self.assertEqual(facts["rows"], 3)

    def test_rejects_wrong_obs_contract(self):
        broken = self.obs.drop(columns=["tahoe_usable"])
        with self.assertRaises(contract.ContractError):
            writers.write_tahoe_h5ad(self.src, self.dir / "x.h5ad", broken,
                                     self.var, {})


class PinnedReaderTest(unittest.TestCase):
    """The trainer env (anndata 0.10.8) must read what the writer produced."""

    def setUp(self):
        if not HAVE_CM_DATA:
            self.skipTest(f"frozen data not under {contract.DATA_ROOT}")
        if not PINNED_READER.exists():
            self.skipTest(f"pinned reader not installed at {PINNED_READER}")

    def test_pinned_reader_reads_cm_artifact(self):
        import tempfile
        AnnData = ad.AnnData
        obs, _, _ = cm_assembled_obs()
        var = small_var()
        X = np.ones((len(obs), len(var)), dtype=np.int32)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cm.h5ad"
            ad.AnnData(X=X, obs=obs, var=var,
                       layers={"log1p_cp10k": np.zeros_like(X, dtype=np.float32)},
                       uns={"schema": "06.3_input/1"}).write_h5ad(path, compression="gzip")
            code = (
                "import anndata as ad, sys;"
                f"a=ad.read_h5ad({str(path)!r});"
                "print(a.shape, a.obs.columns.tolist()[-1], sorted(a.layers))"
            )
            result = subprocess.run([str(PINNED_READER), "-c", code],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("(8,", result.stdout)
            self.assertIn("main_set", result.stdout)
            self.assertIn("log1p_cp10k", result.stdout)


if __name__ == "__main__":
    unittest.main()
