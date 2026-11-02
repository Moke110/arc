"""Builder-level tests: provenance, determinism, and the count invariants of the audit.

The invariant tests read the real frozen inputs and are skipped when they are absent
(they live in the main checkout, not in every worktree); point QC_DATA_ROOT at it.
"""

import gzip
import hashlib
import json
import os
import tempfile
import unittest
from pathlib import Path

import pandas as pd

import build_qc


def data_root() -> Path:
    return Path(os.environ.get("QC_DATA_ROOT", "/home/chang/arc"))


class TestRepoRoot(unittest.TestCase):
    def test_walks_up_to_git(self):
        root = build_qc.repo_root()
        self.assertTrue((root / ".git").exists())
        self.assertTrue((root / "stage2").is_dir())
        # a file deep inside the tree resolves to the same root
        deep = root / "stage2" / "06_prnet_transfer" / "qc" / "execution" / "build_qc.py"
        self.assertEqual(build_qc.repo_root(deep), root)

    def test_raises_outside_a_repo(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "x" / "y.txt"
            p.parent.mkdir(parents=True)
            p.write_text("x")
            with self.assertRaises(RuntimeError):
                build_qc.repo_root(p)


class TestObservationIds(unittest.TestCase):
    def test_temp_env_well_ids_match_the_frozen_view_format(self):
        # literals copied from stage2/02_cm_data/output/datasets-1/temp_env_all-smiles.h5ad
        self.assertEqual(build_qc.tempenv_well_id(1, "B2"), "temp_env::Plate1-B02")
        self.assertEqual(build_qc.tempenv_well_id(16, "N23"), "temp_env::Plate16-N23")
        self.assertEqual(build_qc.tempenv_well_id(4, "J10"), "temp_env::Plate4-J10")

    def test_unparsable_well_is_loud(self):
        with self.assertRaises(ValueError):
            build_qc.tempenv_well_id(1, "02")


class TestDeterministicWriters(unittest.TestCase):
    def test_csv_gz_is_byte_identical_across_runs(self):
        df = pd.DataFrame({"a": [1, 2], "b": [0.5, float("nan")], "c": ["x", "y"]})
        with tempfile.TemporaryDirectory() as td:
            p1, p2 = Path(td) / "a.csv.gz", Path(td) / "b.csv.gz"
            m1 = build_qc.write_csv_gz(df, p1)
            m2 = build_qc.write_csv_gz(df, p2)
            self.assertEqual(m1["sha256_gz"], m2["sha256_gz"])
            self.assertEqual(m1["sha256_uncompressed"], m2["sha256_uncompressed"])
            raw = gzip.open(p1, "rb").read()
            self.assertEqual(hashlib.sha256(raw).hexdigest(), m1["sha256_uncompressed"])
            self.assertEqual(raw.decode().splitlines()[0], "a,b,c")

    def test_different_content_has_a_different_hash(self):
        with tempfile.TemporaryDirectory() as td:
            m1 = build_qc.write_csv_gz(pd.DataFrame({"a": [1]}), Path(td) / "a.csv.gz")
            m2 = build_qc.write_csv_gz(pd.DataFrame({"a": [2]}), Path(td) / "b.csv.gz")
            self.assertNotEqual(m1["sha256_uncompressed"], m2["sha256_uncompressed"])

    def test_json_writer_records_a_resolvable_hash(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "x.json"
            meta = build_qc.write_json({"k": [1, 2]}, p)
            self.assertEqual(build_qc.sha256_file(p), meta["sha256"])
            self.assertEqual(json.loads(p.read_text()), {"k": [1, 2]})


class TestRuleLockIdentity(unittest.TestCase):
    def test_payload_hash_is_deterministic(self):
        refs = {"tempenv_blank_ngenes_threshold": 9242.4742}
        a = build_qc.payload_sha256(build_qc.lock_payload(refs, {}))
        b = build_qc.payload_sha256(build_qc.lock_payload(dict(refs), {}))
        self.assertEqual(a, b)

    def test_payload_hash_changes_when_a_threshold_changes(self):
        base = {"tempenv_blank_ngenes_threshold": 9242.4742}
        moved = {"tempenv_blank_ngenes_threshold": 9242.4743}
        self.assertNotEqual(build_qc.payload_sha256(build_qc.lock_payload(base, {})),
                            build_qc.payload_sha256(build_qc.lock_payload(moved, {})))


def _have_inputs() -> bool:
    return (data_root() / build_qc.TEMPENV_HASH).exists()


@unittest.skipUnless(_have_inputs(), f"frozen inputs not found under {data_root()}")
class TestCountInvariants(unittest.TestCase):
    """The headline counts the ticket and the reconciliation note both quote."""

    @classmethod
    def setUpClass(cls):
        cls.design = build_qc.load_tempenv_design(data_root())

    def test_temp_env_design_table(self):
        d = self.design
        self.assertEqual(len(d), 4558)
        self.assertEqual(d.observation_id.nunique(), 4558)
        kinds = d.Chemical_name.value_counts()
        self.assertEqual(int(kinds["VEH"]), 192)
        self.assertEqual(int(kinds["Media"]), 80)
        self.assertEqual(int(kinds["Blank"]), 32)
        self.assertEqual(int(kinds["ET-1"]), 48)
        treated = d[~d.Chemical_name.isin(["VEH", "Media", "Blank", "ET-1"])]
        self.assertEqual(len(treated), 4206)
        self.assertEqual(treated.groupby("Chemical_Concentration_uM").size().to_dict(),
                         {0.2: 1402, 1.0: 1402, 10.0: 1402})
        self.assertEqual(d.Plate_ID.nunique(), 16)
        self.assertEqual(d.duplicated(["Plate_ID", "Well_ID"]).sum(), 0)

    def test_temp_env_raw_replicate_is_the_source_field_not_the_view_index(self):
        v = self.design.Replicate_Number
        self.assertEqual(int(v.min()), 1)
        self.assertLessEqual(int(v.max()), 12)

    def test_tahoe_physical_wells_are_not_observations(self):
        sm = pd.read_parquet(data_root() / build_qc.TAHOE_META)
        self.assertEqual(len(sm), 1344)                     # physical wells
        self.assertEqual(int((sm.drug == "DMSO_TF").sum()), 28)
        self.assertEqual(sm.plate.nunique(), 14)
        self.assertEqual(int(sm.drug.nunique()), 380)        # 379 drugs + DMSO_TF

    def test_tahoe_view_observation_count(self):
        import anndata as ad
        view = ad.read_h5ad(data_root() / build_qc.TAHOE_VIEW, backed="r")
        obs = view.obs
        self.assertEqual(len(obs), 67018)
        self.assertEqual(int(obs["is_control"].sum()), 1400)
        self.assertEqual(obs["control_group"].nunique(), 700)
        self.assertEqual(obs["cell_line_id"].nunique(), 50)
        # 28 DMSO_TF wells x 50 lines = 1,400 control observations; 2 per group
        per_group = obs.loc[obs["is_control"].astype(bool)].groupby(
            "control_group", observed=True).size()
        self.assertEqual(set(per_group.unique()), {2})


if __name__ == "__main__":
    unittest.main()
