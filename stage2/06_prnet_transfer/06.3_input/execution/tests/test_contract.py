"""Contract tests: the frozen gene axis, handoff schema and locked config.

Offline.  Tests that need the frozen 16 GB sources skip when the data root does
not carry them (the implementation worktree carries tracked files only; run with
`ARC_DATA_ROOT=/home/chang/arc` there).
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

EXECUTION = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXECUTION))

import pandas as pd

import contract

# sha256 over "\n".join(genes.csv symbol order) — the frozen axis identity.
FROZEN_SYMBOL_ORDER_SHA256 = "a9bdf72ac03c714c8f017d19fc63f737a7bb10d7501864af75aefb611fa25826"
FROZEN_ENSEMBL_ORDER_SHA256 = "c9fac320c93c4a7c22260cb79f8c922a27b519e3bdec19d6e28deead6bba0cf2"

HAVE_CM_DATA = contract.CM_VIEW.exists()


class GeneAxisTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not HAVE_CM_DATA:
            raise unittest.SkipTest(f"frozen CM view not under {contract.DATA_ROOT}")
        cls.axis = contract.load_gene_axis()

    def test_axis_size_and_uniqueness(self):
        self.assertEqual(len(self.axis), contract.N_GENES)
        self.assertTrue(self.axis.index.is_unique)
        self.assertEqual(list(self.axis.columns),
                         ["hgnc_id", "entrez_id", "ensembl_id", "match_type"])
        self.assertTrue((self.axis.index.astype(str) != "").all())
        self.assertTrue((self.axis["ensembl_id"].astype(str).str.startswith("ENSG")).all())

    def test_axis_order_is_frozen(self):
        self.assertEqual(contract.sha256_text("\n".join(self.axis.index.astype(str))),
                         FROZEN_SYMBOL_ORDER_SHA256)
        self.assertEqual(contract.sha256_text(
            "\n".join(self.axis["ensembl_id"].astype(str))), FROZEN_ENSEMBL_ORDER_SHA256)

    def test_cm_measurement_mask(self):
        measured = contract.gene_measured_cm(self.axis)
        self.assertEqual(int(measured.sum()), contract.N_GENES,
                         "every axis gene is in the temp_env panel by construction")
        self.assertTrue(bool(measured.loc["FAM153B"]),
                        "FAM153B is measured in CM; it is only missing from Tahoe")

    def test_tahoe_measurement_mask_is_membership_not_counts(self):
        ens = contract.tahoe_gene_space_ensembl()
        measured = contract.gene_measured_tahoe(self.axis, ens)
        self.assertEqual(int((~measured).sum()), 1)
        self.assertFalse(bool(measured.loc["FAM153B"]))
        self.assertTrue(bool(measured.loc["NKX2-1"]),
                        "NKX2-1 is measured but all-zero in Tahoe — never infer "
                        "measurement from non-zero counts")
        self.assertEqual(contract.gene_space_digest(ens),
                         contract.gene_space_digest(sorted(ens)))


class SourceRecordTest(unittest.TestCase):
    def setUp(self):
        if not HAVE_CM_DATA:
            self.skipTest(f"frozen sources not under {contract.DATA_ROOT}")

    def test_records_present_and_sized(self):
        records = contract.source_records()
        roles = [r.role for r in records]
        self.assertEqual(roles, ["cm_view", "tahoe_view", "gene_axis",
                                 "source_gene_sets", "gse262419_hash", "tahoe_corpus"])
        for record in records:
            self.assertTrue(record.path.exists(), record.path)
            self.assertEqual(record.path.stat().st_size, record.bytes, record.role)
        by_role = {r.role: r for r in records}
        self.assertEqual(by_role["cm_view"].content_sha256,
                         "f84066c4dfd2c017785d85a3f6bd54cb2a8eb40bfba28e9e99d650e4b8a613fb")
        self.assertEqual(by_role["tahoe_view"].content_sha256,
                         "476812c09d06c29e56069a3a1e670834384abb524a9468dc596229b5cbcbf407")
        # the corpus's 02 manifest is the authority for the CM view digest
        row = pd.read_csv(contract.CM_BUNDLE_MANIFEST, dtype=str)
        row = row.loc[row["view"] == "temp_env_all-smiles.h5ad"].iloc[0]
        self.assertEqual(by_role["cm_view"].content_sha256, row["content_sha256"])
        self.assertEqual(by_role["cm_view"].bytes, int(row["bytes"]))


class HandoffTest(unittest.TestCase):
    def _write(self, table: pd.DataFrame) -> Path:
        handle = tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False)
        table.to_csv(handle.name, index=False)
        handle.close()
        return Path(handle.name)

    def test_eligibility_valid(self):
        path = self._write(pd.DataFrame({
            "drug_key": ["AAAA-BBBB-CCCCCC", "DDDD-EEEE-FFFFFF"],
            "eligibility_decision": ["INCLUDE", "uncertain"],
            "elig_class": ["approved", "tool"],
        }))
        table, handoff = contract.load_eligibility(path)
        self.assertEqual(len(table), 2)
        self.assertEqual(set(table["eligibility_decision"]), {"include", "uncertain"})
        self.assertEqual(handoff.rows, 2)
        self.assertRegex(handoff.sha256, r"^[0-9a-f]{64}$")
        path.unlink()

    def test_eligibility_rejects_bad_decision_and_duplicates(self):
        bad = self._write(pd.DataFrame({
            "drug_key": ["A"], "eligibility_decision": ["maybe"], "elig_class": ["x"]}))
        with self.assertRaises(contract.ContractError):
            contract.load_eligibility(bad)
        bad.unlink()
        dup = self._write(pd.DataFrame({
            "drug_key": ["A", "A"], "eligibility_decision": ["include", "exclude"],
            "elig_class": ["x", "y"]}))
        with self.assertRaises(contract.ContractError):
            contract.load_eligibility(dup)
        dup.unlink()

    def test_qc_alias_columns_and_booleans(self):
        path = self._write(pd.DataFrame({
            "well": ["w1", "w2"], "pass": ["true", "no"], "reason": ["ok", "low library"]}))
        table, handoff = contract.load_qc(path, "cm")
        self.assertEqual(list(table["qc_pass"]), [True, False])
        self.assertEqual(handoff.kind, "cm_qc")
        path.unlink()

    def test_qc_rejects_unknown_boolean(self):
        path = self._write(pd.DataFrame({
            "well": ["w1"], "qc_pass": ["maybe"], "qc_reason": ["?"]}))
        with self.assertRaises(contract.ContractError):
            contract.load_qc(path, "cm")
        path.unlink()

    def test_split_combined_qc(self):
        table = pd.DataFrame({
            "source": ["cm", "tahoe", "tahoe"],
            "well": ["temp_env::Plate1-B02", "tahoe::plate1:smp_1:CVCL_1", "tahoe::plate1:smp_2:CVCL_2"],
            "qc_pass": ["true", "false", "true"],
            "qc_reason": ["ok", "few cells", "ok"],
        })
        parts = contract.split_combined_qc(table, Path("x.csv"))
        self.assertEqual(set(parts), {"cm", "tahoe"})
        self.assertEqual(len(parts["cm"]), 1)
        self.assertEqual(len(parts["tahoe"]), 2)

    def test_resolve_handoff_missing_fails_closed(self):
        with self.assertRaises(contract.ContractError):
            contract.resolve_handoff("eligibility", override="/nonexistent/elig.csv")


class ColumnContractTest(unittest.TestCase):
    def test_obs_column_lists(self):
        self.assertEqual(len(contract.CM_OBS_INHERITED), 24)
        self.assertEqual(len(contract.TAHOE_OBS_INHERITED), 35)
        self.assertEqual(contract.CM_OBS, contract.CM_OBS_INHERITED + contract.CM_OBS_ADDED)
        self.assertEqual(contract.TAHOE_OBS,
                         contract.TAHOE_OBS_INHERITED + contract.TAHOE_OBS_ADDED)
        for column in ("observation_id", "physical_well", "plate", "cell_line_id"):
            self.assertIn(column, contract.TAHOE_OBS)
        for column in ("eligibility_decision", "elig_class", "qc_pass", "main_set"):
            self.assertIn(column, contract.CM_OBS)
        self.assertEqual(contract.VAR_COLUMNS,
                         ["hgnc_id", "entrez_id", "ensembl_id", "match_type",
                          "gene_measured", "n_nonzero"])

    def test_artifact_names(self):
        self.assertEqual(contract.ARTIFACTS["cm"], "cm_transfer.h5ad")
        self.assertEqual(contract.ARTIFACTS["tahoe"], "tahoe_transfer.h5ad")
        self.assertIn("folds.parquet", contract.ARTIFACTS.values())
        self.assertIn("locked_config.json", contract.ARTIFACTS.values())
        self.assertNotIn(contract.ARTIFACTS["tahoe"], contract.TRACKED_RECORD)
        self.assertNotIn(contract.ARTIFACTS["cm"], contract.TRACKED_RECORD)


class EnvTest(unittest.TestCase):
    def test_repo_root_walk(self):
        self.assertTrue((contract.ROOT / ".git").exists())
        self.assertEqual(contract.find_repo_root(contract.ROOT), contract.ROOT)

    def test_git_state_shape(self):
        state = contract.git_state()
        self.assertRegex(state["sha"], r"^[0-9a-f]{40}$")
        self.assertIsInstance(state["dirty"], bool)

    def test_runtime_env_has_versions(self):
        env = contract.runtime_env()
        for key in ("python", "anndata", "numpy", "pandas"):
            self.assertIn(key, env)


if __name__ == "__main__":
    unittest.main()
